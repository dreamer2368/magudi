"""File-aligned PETSc Vec layout + parallel I/O for optim_ver4.

ParallelIOHandler owns:
  - the concatenated-control-vector schema parsed from <prefix>.layout.txt,
  - the per-rank file-aligned local Vec sizes,
  - the per-actuator .dat sub-communicators,
  - per-callback writes and reads of the full slot set.

Resource policy (anything else raises RuntimeError):

  N_petsc == 1
      Rank 0 owns every slot. Serial I/O for both .dat and .q.

  N_petsc = n_ic + n_actuator * nprocs_per_actuator (with nprocs_per_actuator >= 1)
      Ranks [0, n_actuator * nprocs_per_actuator)
                       form n_actuator sub-comms (one per .dat file).
                       Sub-comm s gets ranks
                         [s*nprocs_per_actuator, (s+1)*nprocs_per_actuator)
                       and shares actuator slot s via collective MPI-IO.
      Remaining ranks  each own exactly one ic slot in schema order,
                       handled serially by plot3dnasa.

The intermediate regime is rejected; round-robin "multiple ic slots per rank"
or "uneven nprocs_per_actuator" can be added later if production needs it.

Each ic slot stores a multi-block PLOT3D solution. Geometry (number of
blocks + per-block sizes + nUnknowns) is cached at init by inspecting the
shared <prefix>.norm_ic.q produced by compute_norm -- the same file
loaded into every ic slot of M_diag. No dependency on .xyz files.

Wire format follows schema kind unambiguously:
   actuator -> .dat (raw real64 with the 5D-subarray layout written by
                     ActuatorPatch / read by msforward/msadjoint)
   ic       -> .q   (multi-block PLOT3D solution; compute_norm writes a
                     single shared <prefix>.norm_ic.q for the metric, and
                     Python passes the same path to every ic slot when
                     loading M_diag)
"""

import glob
import os
import re
import subprocess
from collections import namedtuple
from typing import List, Optional

import h5py
import numpy as np

from mpi4py import MPI  # noqa: F401  (import-order side effect)
from petsc4py import PETSc

from . import plot3dnasa
from .inputs import InputParser


def launch_and_wait(executable, args, n, mode):
    """Collectively launch `n` worker processes via the mode's launcher; block.

    All N_petsc parent ranks Barrier; rank 0 subprocess.runs the launcher
    and blocks until it returns; the return code is broadcast so every
    rank either continues or raises SystemExit together; then Barrier
    again. The worker's stdout/stderr land in ./out/<basename>.out.

    mode='base'  -> launcher = ["mpirun"]
    mode='slurm' -> launcher = ["srun", "--overlap"]

    The Fortran child's disconnectParentIfSpawned (MPIHelperImpl.f90:552-574)
    is a no-op here -- MPI_Comm_get_parent returns MPI_COMM_NULL because
    we launched via mpirun/srun, not MPI_Comm_spawn. Sync is via the
    subprocess.run blocking call, not an inter-comm Barrier.
    """
    comm = MPI.COMM_WORLD
    log_path = os.path.abspath(
        os.path.join("out", os.path.basename(executable) + ".out")
    )
    if comm.Get_rank() == 0:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        open(log_path, "w").close()
    comm.Barrier()

    if comm.Get_rank() == 0:
        launcher = ["mpirun"] if mode == "base" else ["srun", "--overlap"]
        cmd = launcher + ["-n", str(n), executable] + list(args)
        with open(log_path, "a") as logf:
            rc = subprocess.run(
                cmd, stdout=logf, stderr=subprocess.STDOUT
            ).returncode
    else:
        rc = None
    rc = comm.bcast(rc, root=0)
    comm.Barrier()
    if rc != 0:
        raise SystemExit(
            f"{executable} (mode={mode}) exited with code {rc}; "
            f"see {log_path} for details"
        )


FileSpec = namedtuple("FileSpec", ["kind", "identifier", "size"])
# kind ∈ {"actuator", "ic"} — semantic role in the layout, not wire format.


def _flatten_magudi_dict(d, parent="", sep="/"):
    """{'a': {'b': 1}, 'c': 2}  ->  {'a/b': 1, 'c': 2}.

    Magudi.inp stores hierarchical keys flat with '/' separators
    (e.g. time_splitting/number_of_segments). YAML nested dicts become these
    slash-joined keys.
    """
    out = {}
    for k, v in d.items():
        key = f"{parent}{sep}{k}" if parent else str(k)
        if isinstance(v, dict):
            out.update(_flatten_magudi_dict(v, key, sep))
        else:
            out[key] = v
    return out


def _format_magudi_val(v):
    """Render a Python value as a magudi.inp RHS.

    bool -> lowercase 'true'/'false' (matches examples/OneDWave/magudi.inp:14).
    str  -> always quoted; magudi requires quotes when the value contains
            whitespace (examples/magudi.inp:8) and accepts quotes otherwise.
    int / float / other -> str(v).
    """
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, str):
        return f'"{v}"'
    return str(v)


def apply_magudi_inp(cfg) -> None:
    """Rewrite magudi.inp in place from cfg['magudi']['forced_inputs'].

    Standalone (no ParallelIOHandler / MPI required). Used by run_parallel.sh
    before compute_norm so that time_splitting/* and the other forced_inputs
    are visible to compute_norm's getRequiredOption reads, and used internally
    by ParallelIOHandler.set_magudi_inp on rank 0. No-op when forced_inputs
    is empty or absent.
    """
    input_file = cfg.getInput(
        ["magudi", "input_file"], fallback="magudi.inp")
    forced = cfg.getInput(
        ["magudi", "forced_inputs"], fallback={})
    if not forced:
        return
    flat = _flatten_magudi_dict(forced)
    with open(input_file) as fh:
        lines = fh.readlines()
    seen = set()
    key_re = re.compile(r"^\s*([^\s#=][^\s=]*)\s*=")
    for i, line in enumerate(lines):
        m = key_re.match(line)
        if not m:
            continue
        key = m.group(1)
        if key in flat:
            lines[i] = f"{key} = {_format_magudi_val(flat[key])}\n"
            seen.add(key)
    for key, val in flat.items():
        if key not in seen:
            lines.append(f"{key} = {_format_magudi_val(val)}\n")
    with open(input_file, "w") as fh:
        fh.writelines(lines)


def parse_layout(path: str) -> List[FileSpec]:
    """Parse the rank-0 layout.txt written by bin/compute_norm.

    Records: ``actuator <name> <size>`` or ``ic <k> <size>``. Order in the
    file is the canonical concatenation order Python uses for the Vec.
    """
    schema = []
    with open(path) as fh:
        for line in fh:
            payload = line.split("#", 1)[0].strip()
            if not payload:
                continue
            parts = payload.split()
            if len(parts) != 3:
                raise ValueError(
                    f"{path}: expected 'kind identifier size_in_real64'; got {parts!r}"
                )
            kind, identifier, size_str = parts
            if kind not in ("actuator", "ic"):
                raise ValueError(
                    f"{path}: unknown kind {kind!r}; expected 'actuator' or 'ic'"
                )
            try:
                size = int(size_str)
            except ValueError as exc:
                raise ValueError(
                    f"{path}: size must be integer; got {size_str!r}"
                ) from exc
            schema.append(FileSpec(kind=kind, identifier=identifier, size=size))
    return schema


class ParallelIOHandler:
    """File-aligned PETSc Vec layout + parallel I/O.

    Public API:
        create_vec()           -> PETSc.Vec with the file-aligned local size
        global_size, local_size
        report_balance()       -> log max/min(nonzero) local-size ratio
        write(vec, paths)      -> write Vec contents into schema slots
        read(vec, paths)       -> read schema slots into the Vec
        slice_of(slot_index)   -> (lo, hi) global-index range
    """

    def __init__(self, cfg: InputParser,
                 comm: MPI.Comm = MPI.COMM_WORLD,
                 pcomm = PETSc.COMM_WORLD):
        """cfg: parsed YAML. Every config-dependent value (prefix, Nsplit,
        schema, ic_norm_q_path, log path, history-file paths) is derived from
        cfg here; nothing else needs to be passed in.
        comm: MPI comm used for parent-side rank classification and
        per-actuator MPI-IO sub-communicators.
        pcomm: PETSc comm used for binary .petsc Viewer / Vec.createMPI calls
        (read_petsc / write_petsc / save_tao_state / load_tao_state).
        Defaults to PETSc.COMM_WORLD.
        """
        self.cfg = cfg
        self.prefix = cfg.getInput(["global_prefix"], datatype=str)
        self.Nsplit = cfg.getInput(
            ["magudi", "forced_inputs", "time_splitting", "number_of_segments"],
            datatype=int)
        self.log_file_path = cfg.getInput(
            ["optimization", "log_file"],
            fallback=f"{self.prefix}.optim.h5")

        layout_path = f"{self.prefix}.layout.txt"
        self.schema = parse_layout(layout_path)
        if not self.schema:
            raise ValueError(f"ParallelIOHandler: empty schema in {layout_path}")
        self.ic_norm_q_path = f"{self.prefix}.norm_ic.q"

        # msforward / msadjoint mutate solver%outputPrefix = <prefix>-<k> per
        # segment (bin/msforward.f90:240-241, bin/msadjoint.f90:296-297), so
        # the per-timestep history files are per-segment.
        self.fwd_hist_paths = [f"{self.prefix}-{k}.cost_functional.txt"
                               for k in range(self.Nsplit)]
        self.adj_hist_paths = [f"{self.prefix}-{k}.cost_sensitivity.txt"
                               for k in range(self.Nsplit)]
        self.sub_j_path   = f"{self.prefix}.sub_J.txt"            # msforward subOutputFilename
        self.sub_adj_path = f"{self.prefix}.sub_adjoint_run.txt"  # msadjoint subOutputFilename
        self.j_path       = f"{self.prefix}.forward_run.txt"      # msforward total J output
        self.gg_path      = f"{self.prefix}.adjoint_run.txt"      # msadjoint total <g,g>_M output

        self.comm = comm
        self.pcomm = pcomm
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()

        self._classify_slots()
        if self.n_ic > 0 and self.ic_norm_q_path is None:
            raise RuntimeError(
                "ParallelIOHandler: ic_norm_q_path is required when the schema "
                f"has ic slots (n_ic={self.n_ic})."
            )
        self._assign_ownership()
        self._build_dat_subcomm()
        self._cache_q_geometry(self.ic_norm_q_path)
        self._validate_ic_geometry()

        self.x_paths = self._build_paths("x")
        self.grad_paths = self._build_paths("grad")
        self.metric_paths = self._build_paths("metric")

        self._init_log_file()

    # ------------------------------------------------------------------ layout

    def _classify_slots(self):
        self.actuator_slots = [i for i, s in enumerate(self.schema)
                               if s.kind == "actuator"]
        self.ic_slots = [i for i, s in enumerate(self.schema) if s.kind == "ic"]
        self.n_actuator = len(self.actuator_slots)
        self.n_ic = len(self.ic_slots)
        # Canonical order: all actuator slots precede all ic slots.
        for i in self.actuator_slots:
            for j in self.ic_slots:
                if j < i:
                    raise ValueError(
                        "schema order violated: all 'actuator' slots must precede 'ic' slots"
                    )
        offsets = [0]
        for s in self.schema:
            offsets.append(offsets[-1] + s.size)
        self.global_size = offsets[-1]
        self._slot_range = {
            i: (offsets[i], offsets[i + 1]) for i in range(len(self.schema))
        }

    def _assign_ownership(self):
        n_actuator = self.n_actuator
        n_ic = self.n_ic
        N = self.size

        # Serial fallback.
        if N == 1:
            self.is_actuator_rank = (n_actuator > 0)
            self.is_ic_rank = (n_ic > 0)
            self.actuator_slot_owned = None
            self.actuator_subrank = 0
            self.nprocs_per_actuator = max(1, n_actuator)  # for symmetry; unused
            self.ic_owned = list(self.ic_slots)
            self.local_offset = 0
            self.local_size = self.global_size
            return

        # n_actuator > 0 requires n_dat * nprocs_per_actuator ranks; n_ic requires n_ic.
        if n_actuator == 0:
            # Pure-IC optimization: every rank owns one ic slot.
            if N != n_ic:
                raise RuntimeError(
                    f"N_petsc must equal n_ic = {n_ic} when n_actuator = 0; "
                    f"got N_petsc = {N}"
                )
            self.nprocs_per_actuator = 0
        else:
            excess = N - n_ic
            if excess <= 0 or excess % n_actuator != 0:
                raise RuntimeError(
                    "N_petsc = n_ic + n_actuator * nprocs_per_actuator required "
                    "(nprocs_per_actuator >= 1); "
                    f"got N_petsc = {N}, n_ic = {n_ic}, n_actuator = {n_actuator}"
                )
            self.nprocs_per_actuator = excess // n_actuator

        n_actuator_total = n_actuator * self.nprocs_per_actuator

        if self.rank < n_actuator_total:
            self.is_actuator_rank = True
            self.is_ic_rank = False
            s = self.rank // self.nprocs_per_actuator
            sr = self.rank % self.nprocs_per_actuator
            self.actuator_slot_owned = self.actuator_slots[s]
            self.actuator_subrank = sr
            slot_lo, slot_hi = self._slot_range[self.actuator_slot_owned]
            slot_size = slot_hi - slot_lo
            lo_in_slot = (sr * slot_size) // self.nprocs_per_actuator
            hi_in_slot = ((sr + 1) * slot_size) // self.nprocs_per_actuator
            self.local_offset = slot_lo + lo_in_slot
            self.local_size = hi_in_slot - lo_in_slot
            self.ic_owned = []
        else:
            self.is_actuator_rank = False
            self.is_ic_rank = True
            ic_index = self.rank - n_actuator_total
            slot_id = self.ic_slots[ic_index]
            self.ic_owned = [slot_id]
            self.actuator_slot_owned = None
            self.actuator_subrank = 0
            lo, hi = self._slot_range[slot_id]
            self.local_offset = lo
            self.local_size = hi - lo

    def _build_dat_subcomm(self):
        # One sub-comm per actuator. Color = index in self.actuator_slots so each
        # actuator's sub-comm is independent and writes/reads its own file.
        if self.size == 1:
            # Serial: no MPI sub-comm needed; rank 0 uses POSIX I/O for everything.
            self.actuator_subcomm = MPI.COMM_NULL
            return
        if self.is_actuator_rank:
            color = self.actuator_slots.index(self.actuator_slot_owned)
        else:
            color = MPI.UNDEFINED
        # key=self.rank preserves the parent rank ordering inside the sub-comm,
        # so each sub-rank's slab is contiguous in the actuator's global byte range.
        self.actuator_subcomm = self.comm.Split(color=color, key=self.rank)

    def _cache_q_geometry(self, ic_norm_q_path: Optional[str]):
        # Parse <prefix>.norm_ic.q's PLOT3D header once so each ic-rank can
        # construct a Solution of the right shape when reading/writing slabs.
        self.q_nblocks = 0
        self.grid_sizes = []     # list of (Nx, Ny, Nz)
        self.q_nunknowns = None
        if not self.ic_slots or ic_norm_q_path is None:
            return
        if plot3dnasa is None:
            return  # Deferred error: raised at the first .q I/O call.
        fmt = plot3dnasa.FileFormat(ic_norm_q_path)
        if fmt.file_type != "solution":
            raise ValueError(
                f"{ic_norm_q_path}: expected a PLOT3D solution file (.q); got file_type={fmt.file_type!r}"
            )
        self.q_nblocks = int(fmt.nblocks)
        self.grid_sizes = [tuple(int(x) for x in fmt.size[b])
                           for b in range(self.q_nblocks)]
        total_pts = int(sum(np.prod(s) for s in self.grid_sizes))
        if total_pts == 0:
            raise ValueError(f"{ic_norm_q_path}: total grid points is zero")
        sample = self.schema[self.ic_slots[0]].size
        if sample % total_pts != 0:
            raise ValueError(
                f"ic slot size {sample} not divisible by total grid points {total_pts}; "
                f"ic_norm_q_path={ic_norm_q_path!r}"
            )
        self.q_nunknowns = sample // total_pts

    def _validate_ic_geometry(self):
        # All ic slots must agree on size.
        if not self.ic_slots:
            return
        sample = self.schema[self.ic_slots[0]].size
        for i in self.ic_slots[1:]:
            if self.schema[i].size != sample:
                raise ValueError(
                    f"all ic slots must have identical size; "
                    f"slot {self.ic_slots[0]} size {sample}, slot {i} size {self.schema[i].size}"
                )

    # ------------------------------------------------------------------ public

    def slice_of(self, slot_index: int):
        return self._slot_range[slot_index]

    def create_vec(self) -> PETSc.Vec:
        vec = PETSc.Vec().createMPI(
            size=(self.local_size, self.global_size), comm=self.comm
        )
        vec.set(0.0)
        return vec

    def report_balance(self) -> None:
        sizes = self.comm.allgather(self.local_size)
        nonzero = [s for s in sizes if s > 0]
        ratio = (max(sizes) / min(nonzero)) if nonzero else float("nan")
        PETSc.Sys.Print(
            f"ParallelIOHandler: N_petsc={self.size} local_sizes={sizes} "
            f"max/min(nonzero)={ratio:.3g}"
        )

    def write(self, vec: PETSc.Vec, paths: List[str]) -> None:
        self._check_paths(paths)
        local_arr = np.ascontiguousarray(vec.getArray(readonly=True), dtype="<f8")
        if self.size == 1:
            self._serial_write(local_arr, paths)
        else:
            if self.is_actuator_rank:
                self._write_owned_actuator(local_arr, paths)
            if self.is_ic_rank:
                self._write_owned_ic(local_arr, paths)
        self.comm.Barrier()

    def read(self, paths: List[str]) -> PETSc.Vec:
        """Allocate a fresh file-aligned Vec, populate it from `paths`, return it."""
        self._check_paths(paths)
        vec = self.create_vec()
        local_arr = np.empty(self.local_size, dtype="<f8")
        if self.size == 1:
            self._serial_read(local_arr, paths)
        else:
            if self.is_actuator_rank:
                self._read_owned_actuator(local_arr, paths)
            if self.is_ic_rank:
                self._read_owned_ic(local_arr, paths)
        self.comm.Barrier()
        vec.setArray(local_arr)
        vec.assemble()
        return vec

    # ---- convenience wrappers using cached schema/prefix/ic_norm_q_path ----

    def read_x(self) -> PETSc.Vec:
        """Allocate and return a Vec populated from x slabs
        (.control_forcing_<name>.dat + -<k>.ic.q)."""
        return self.read(self.x_paths)

    def write_x(self, vec: PETSc.Vec) -> None:
        """Write vec into x slabs (.control_forcing_<name>.dat + -<k>.ic.q)."""
        self.write(vec, self.x_paths)

    def read_grad(self) -> PETSc.Vec:
        """Allocate and return a Vec populated from raw gradient slabs
        (.gradient_<name>.dat + -<k>.ic.adjoint.q)."""
        return self.read(self.grad_paths)

    def read_metric(self) -> PETSc.Vec:
        """Allocate and return a Vec populated from M_diag slabs
        (.norm_<name>.dat + shared .norm_ic.q)."""
        return self.read(self.metric_paths)

    def cleanup_iteration_artifacts(self) -> None:
        """Remove per-iteration files produced by msforward / msadjoint.

        Every rank independently builds the same `targets` list (the glob
        runs on a shared filesystem and `targets` is sorted so the round-robin
        partition is identical across ranks). Rank r then unlinks
        targets[r::N_petsc]; we end on Barrier so the next spawn sees a clean
        filesystem.

        Inputs that the driver itself writes (.ic.q baseline ICs and
        .control_forcing_<name>.dat) are preserved.

        Targets (driven by self.schema + self.prefix):
          - <prefix>-<ts:08d>.q          per-segment forward snapshots (glob)
          - <prefix>-<ts:08d>.adjoint.q  per-segment adjoint snapshots (glob)
          - <prefix>-<k>.ic.adjoint.q    per-IC adjoint state (grad_paths ic entries)
          - <prefix>.gradient_<name>.dat per-actuator gradient (grad_paths actuator entries)
          - <prefix>.cost_functional.txt
          - <prefix>.cost_sensitivity.txt
          - <prefix>.forward_run.txt
          - <prefix>.adjoint_run.txt
          - <prefix>.sub_J.txt
          - <prefix>.sub_adjoint_run.txt
        """
        targets = []
        # Per-segment .q snapshots: <prefix>-<8-digit-ts>.q and the matching
        # .adjoint.q. The glob is filtered to skip the .ic.q baseline ICs
        # (driver-owned inputs) and the .ic.adjoint.q files (covered via
        # self.grad_paths below).
        stem_lo = len(self.prefix) + 1  # strip "<prefix>-"
        for path in glob.glob(f"{self.prefix}-*.q"):
            base = os.path.basename(path)
            stem = base[stem_lo:-2]  # strip "<prefix>-" prefix and ".q" suffix
            if stem.endswith(".ic") or stem.endswith(".ic.adjoint"):
                continue
            if stem.isdigit() and len(stem) == 8:
                targets.append(path)
            elif stem.endswith(".adjoint") and stem[:-len(".adjoint")].isdigit():
                targets.append(path)
        # Slot-level outputs (one per actuator / ic in schema order).
        targets.extend(self.grad_paths)
        # Fixed text artifacts.
        for name in ("cost_functional.txt", "cost_sensitivity.txt",
                     "forward_run.txt", "adjoint_run.txt",
                     "sub_J.txt", "sub_adjoint_run.txt"):
            targets.append(f"{self.prefix}.{name}")

        # Sort so the round-robin partition agrees across ranks regardless of
        # glob.glob() order (which is OS-dependent and not guaranteed across
        # nodes on a shared FS).
        targets.sort()
        for path in targets[self.rank::self.size]:
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass
        self.comm.Barrier()

    # ---- per-iteration HDF5 history log (rank-0 only) ---------------------

    def _init_log_file(self) -> None:
        """Create empty resizable datasets and labels in self.log_file_path.

        If self.log_file_path already exists, return immediately — the file
        keeps its existing structure and history across resumes. Otherwise
        create the file from scratch with all groups/datasets/labels in one
        pass.
        """
        if self.rank == 0 and not os.path.exists(self.log_file_path):
            Nsplit = self.Nsplit
            str_dt = h5py.string_dtype(encoding="utf-8")
            with h5py.File(self.log_file_path, "w") as f:
                obj_grp = f.create_group("Objectives")
                obj_labels = (
                    [f"J{k}" for k in range(Nsplit)]
                    + [f"L2sq{k}" for k in range(Nsplit)]
                    + ["penalty_weight"]
                )
                obj_grp.create_dataset(
                    "labels", data=np.array(obj_labels, dtype=object), dtype=str_dt)
                obj_grp.create_dataset(
                    "subJ", shape=(0, 2 * Nsplit + 1),
                    maxshape=(None, 2 * Nsplit + 1),
                    chunks=(1, 2 * Nsplit + 1), dtype="f8")

                grad_grp = f.create_group("Gradients")
                grad_labels = (
                    [f"ctrl{k}" for k in range(Nsplit)]
                    + [f"ic{k}" for k in range(Nsplit)]
                )
                grad_grp.create_dataset(
                    "labels", data=np.array(grad_labels, dtype=object), dtype=str_dt)
                grad_grp.create_dataset(
                    "subgrad", shape=(0, 2 * Nsplit),
                    maxshape=(None, 2 * Nsplit),
                    chunks=(1, 2 * Nsplit), dtype="f8")

                th_grp = f.create_group("TimeHistory")
                # time / timesteps width is set on the first log_iteration call
                # (we don't know the total timestep count here without reading
                # the history files).
                th_grp.create_dataset(
                    "time", shape=(0,), maxshape=(None,), dtype="f8")
                th_grp.create_dataset(
                    "timesteps", shape=(0,), maxshape=(None,), dtype="i8")
                th_grp.create_dataset(
                    "forward", shape=(0, 0), maxshape=(None, None),
                    chunks=True, dtype="f8")
                th_grp.create_dataset(
                    "adjoint", shape=(0, 0), maxshape=(None, None),
                    chunks=True, dtype="f8")

                f.create_dataset(
                    "line_steps", shape=(0,), maxshape=(None,),
                    chunks=(1,), dtype="f8")
        self.comm.Barrier()

    def set_magudi_inp(self) -> None:
        """Rank-0 wrapper around apply_magudi_inp; comm barrier at end."""
        if self.rank == 0:
            apply_magudi_inp(self.cfg)
        self.comm.Barrier()

    def _read_segment_table(self, path: str) -> np.ndarray:
        """Rank-0 helper. Read the msforward / msadjoint sub-output table.

        Returns shape (Nsplit, n_cols); column 0 is the segment index.
        """
        table = np.loadtxt(path, comments="#")
        if table.ndim == 1:
            table = table.reshape(1, -1)
        if table.shape[0] != self.Nsplit:
            raise RuntimeError(
                f"{path}: expected {self.Nsplit} segment rows, got {table.shape[0]}"
            )
        return table

    def _read_history_files(self, paths):
        """Rank-0 helper. Read 4-column per-segment history files; concat & sort.

        Each file is a (rows, 4) table written by writeFunctionalToFile /
        writeSensitivityToFile (cols: timestep, time, instant, running_quad).
        Returns 1D (timesteps, time, instant) sorted by time.
        """
        tables = []
        for p in paths:
            tbl = np.loadtxt(p, comments="#")
            if tbl.ndim == 1:
                tbl = tbl.reshape(1, -1)
            if tbl.shape[1] != 4:
                raise RuntimeError(
                    f"{p}: expected 4 columns (timestep, time, instant, running_quad); "
                    f"got {tbl.shape[1]}"
                )
            tables.append(tbl)
        merged = np.concatenate(tables, axis=0)
        order = np.argsort(merged[:, 1], kind="stable")
        merged = merged[order]
        return (merged[:, 0].astype("i8"),
                merged[:, 1].astype("f8"),
                merged[:, 2].astype("f8"))

    def log_iteration(self, penalty_weight: float, line_step: float) -> None:
        """Append one row per dataset for the current accepted TAO iteration.

        Reads:
          - <prefix>.sub_J.txt              (segmentCost, segmentL2sq per segment)
          - <prefix>.sub_adjoint_run.txt    (segmentCtrlIP, segmentICIP per segment)
          - <prefix>-<k>.cost_functional.txt   (per-timestep forward values)
          - <prefix>-<k>.cost_sensitivity.txt  (per-timestep adjoint values)
        Appends to Objectives/subJ, Gradients/subgrad, TimeHistory/forward,
        TimeHistory/adjoint, and line_steps.
        """
        if self.rank == 0:
            sub_j   = self._read_segment_table(self.sub_j_path)    # cols: idx, J, L2sq
            sub_adj = self._read_segment_table(self.sub_adj_path)  # cols: idx, ctrlIP, icIP
            segmentCost   = sub_j[:, 1]
            segmentL2sq   = sub_j[:, 2]
            segmentCtrlIP = sub_adj[:, 1]
            segmentICIP   = sub_adj[:, 2]

            ts_fwd, time_fwd, inst_fwd = self._read_history_files(self.fwd_hist_paths)
            _, time_adj, inst_adj = self._read_history_files(self.adj_hist_paths)
            if len(time_adj) != len(time_fwd):
                raise RuntimeError(
                    f"forward/adjoint history length mismatch: "
                    f"forward={len(time_fwd)} vs adjoint={len(time_adj)}"
                )

            with h5py.File(self.log_file_path, "a") as f:
                subJ = f["Objectives/subJ"]
                row = np.concatenate([segmentCost, segmentL2sq, [penalty_weight]])
                subJ.resize(subJ.shape[0] + 1, axis=0)
                subJ[-1] = row

                subgrad = f["Gradients/subgrad"]
                row = np.concatenate([segmentCtrlIP, segmentICIP])
                subgrad.resize(subgrad.shape[0] + 1, axis=0)
                subgrad[-1] = row

                time_ds = f["TimeHistory/time"]
                ts_ds   = f["TimeHistory/timesteps"]
                fwd_ds  = f["TimeHistory/forward"]
                adj_ds  = f["TimeHistory/adjoint"]
                N = len(time_fwd)
                if time_ds.shape[0] == 0:
                    time_ds.resize((N,))
                    time_ds[:] = time_fwd
                    ts_ds.resize((N,))
                    ts_ds[:] = ts_fwd
                    fwd_ds.resize((0, N))
                    adj_ds.resize((0, N))
                fwd_ds.resize(fwd_ds.shape[0] + 1, axis=0)
                fwd_ds[-1] = inst_fwd
                adj_ds.resize(adj_ds.shape[0] + 1, axis=0)
                adj_ds[-1] = inst_adj

                ls_ds = f["line_steps"]
                ls_ds.resize(ls_ds.shape[0] + 1, axis=0)
                ls_ds[-1] = line_step

        self.comm.Barrier()

    # ---- .petsc binary I/O for PETSc Vecs ---------------------------------

    def write_petsc(self, vec: PETSc.Vec, path: str) -> None:
        """Write a PETSc Vec to a binary .petsc file collectively on self.pcomm."""
        viewer = PETSc.Viewer().createBinary(path, mode="w", comm=self.pcomm)
        vec.view(viewer)
        viewer.destroy()

    def read_petsc(self, path: str) -> PETSc.Vec:
        """Load a Vec from `path` with the file-aligned layout.

        Pre-allocates the target via self.create_vec() so PETSc reads into
        the ParallelIOHandler partition (matches what every other Vec the
        driver uses). Without this, VecLoad would distribute the file
        contents evenly across ranks and downstream VecCopy / pointwiseMult
        would fail with PETSc error code 75 (incompatible local lengths).
        """
        viewer = PETSc.Viewer().createBinary(path, mode="r", comm=self.pcomm)
        vec = self.create_vec()
        vec.load(viewer)
        viewer.destroy()
        return vec

    # ---- scalar .txt I/O (rank 0 reads/writes, broadcasts on read) --------

    def read_scalar(self, path: str) -> float:
        """Rank 0 reads a single float from `path`; broadcast to all ranks."""
        if self.rank == 0:
            with open(path) as fh:
                val = float(fh.read().strip())
        else:
            val = None
        return self.comm.bcast(val, root=0)

    def write_scalar(self, path: str, value: float) -> None:
        """Rank 0 writes `value` to `path` as text. No broadcast (value is
        already known on every rank)."""
        if self.rank == 0:
            with open(path, "w") as fh:
                fh.write(f"{value}\n")

    # ---- driver-side bootstrap and TAO state checkpointing ----------------

    def seed_zero_actuator_files(self) -> None:
        """Create zero-filled .control_forcing_<name>.dat for every actuator
        slot, sized to its layout entry. Idempotent: skipped when the file
        already exists at the expected size. Needed only on the first run,
        before the very first read_x for the baseline x.
        """
        if self.rank == 0:
            for s, path in zip(self.schema, self.x_paths):
                if s.kind != "actuator":
                    continue
                expected = s.size * 8
                actual = os.path.getsize(path) if os.path.exists(path) else -1
                if actual != expected:
                    with open(path, "wb") as fh:
                        fh.write(b"\x00" * expected)
        self.comm.Barrier()

    def save_tao_state(self, tao, checkpoint_path: str, history_path: str,
                       snapshots) -> None:
        """Persist iterate + L-BFGS replay buffer + J0 diagonal.

        `snapshots` is a deque of (vx, vg) tuples of distributed Vecs with the
        file-aligned layout. J0 is persisted explicitly because
        M.updateLMVM replay reconstructs (s_k, y_k) but not the nested
        MATLMVMDIAGBRDN diagonal.
        """
        x = tao.getSolution()
        self.write_petsc(x, checkpoint_path)
        PETSc.Sys.Print(f"Saved {checkpoint_path}: |y| = {x.norm():.6e}")

        viewer = PETSc.Viewer().createBinary(history_path, mode="w",
                                             comm=self.pcomm)
        local_size = 3 if self.rank == 0 else 0
        hdr = PETSc.Vec().createMPI((local_size, 3), comm=self.pcomm)
        if self.rank == 0:
            hdr.setValues([0, 1, 2],
                          [float(tao.getIterationNumber()),
                           float(len(snapshots)),
                           float(self.global_size)])
        hdr.assemble()
        hdr.view(viewer)
        hdr.destroy()
        for vx, vg in snapshots:
            vx.view(viewer)
            vg.view(viewer)
        J0_diag = tao.getLMVMMat().getLMVMJ0().getDiagonal()
        J0_diag.view(viewer)
        J0_diag.destroy()
        viewer.destroy()
        PETSc.Sys.Print(
            f"Saved {history_path}: {len(snapshots)} snapshots + J0 diag, "
            f"iter={tao.getIterationNumber()}"
        )

    def load_tao_state(self, tao, checkpoint_path: str,
                       history_path: str) -> bool:
        """Resume iterate; replay L-BFGS history + J0 if present.

        Caller must have done tao.setSolution(y) AND tao.setUp() so the inner
        LMVM Mat is allocated and updateLMVM / setLMVMJ0 target the right
        object. Returns True if a checkpoint was loaded, False otherwise.
        """
        if not os.path.exists(checkpoint_path):
            PETSc.Sys.Print(
                f"No checkpoint at {checkpoint_path}; starting from current y"
            )
            return False
        x = tao.getSolution()
        x_loaded = self.read_petsc(checkpoint_path)
        if x_loaded.getSize() != self.global_size:
            raise SystemExit(
                f"{checkpoint_path} size {x_loaded.getSize()} "
                f"!= expected {self.global_size}"
            )
        x_loaded.copy(x)
        x_loaded.destroy()
        PETSc.Sys.Print(
            f"Resumed iterate from {checkpoint_path}: |y| = {x.norm():.6e}"
        )

        if not (os.path.exists(history_path)
                and os.path.getsize(history_path) > 0):
            PETSc.Sys.Print(
                f"No L-BFGS history at {history_path}; using fresh L-BFGS"
            )
            return True

        viewer = PETSc.Viewer().createBinary(history_path, mode="r",
                                             comm=self.pcomm)
        hdr = PETSc.Vec().load(viewer)
        local_vals = list(hdr.getArray(readonly=True))
        gathered = self.comm.allgather(local_vals)
        flat = [v for sub in gathered for v in sub]
        iter_no, n_snap, n_dim = int(flat[0]), int(flat[1]), int(flat[2])
        hdr.destroy()
        if n_dim != self.global_size:
            raise SystemExit(
                f"history N={n_dim} != current N={self.global_size}")
        M = tao.getLMVMMat()
        for _ in range(n_snap):
            # Same layout-pinning rationale as read_petsc: snapshots and J0
            # diag share the file-aligned partition.
            vx = self.create_vec(); vx.load(viewer)
            vg = self.create_vec(); vg.load(viewer)
            M.updateLMVM(vx, vg)
            vx.destroy()
            vg.destroy()
        J0_diag = self.create_vec(); J0_diag.load(viewer)
        viewer.destroy()
        M.setLMVMJ0(PETSc.Mat().createDiagonal(J0_diag))
        tao.setIterationNumber(iter_no)
        PETSc.Sys.Print(
            f"Replayed {n_snap} snapshots + J0 diag; iter counter set to {iter_no}"
        )
        return True

    def _build_paths(self, mode: str) -> List[str]:
        """One path per slot for a given access mode.

        mode = 'x'      -> control vector slabs
                           actuator -> <prefix>.control_forcing_<name>.dat
                           ic       -> <prefix>-<k>.ic.q
        mode = 'grad'   -> raw gradient slabs (msadjoint outputs)
                           actuator -> <prefix>.gradient_<name>.dat
                           ic       -> <prefix>-<k>.ic.adjoint.q
        mode = 'metric' -> M_diag slabs (compute_norm outputs)
                           actuator -> <prefix>.norm_<name>.dat
                           ic       -> self.ic_norm_q_path (shared)
        """
        paths = []
        for s in self.schema:
            if s.kind == "actuator":
                if mode == "x":
                    paths.append(f"{self.prefix}.control_forcing_{s.identifier}.dat")
                elif mode == "grad":
                    paths.append(f"{self.prefix}.gradient_{s.identifier}.dat")
                elif mode == "metric":
                    paths.append(f"{self.prefix}.norm_{s.identifier}.dat")
                else:
                    raise ValueError(f"unknown mode {mode!r}")
            else:  # kind == "ic"
                k = int(s.identifier)
                if mode == "x":
                    paths.append(f"{self.prefix}-{k}.ic.q")
                elif mode == "grad":
                    paths.append(f"{self.prefix}-{k}.ic.adjoint.q")
                elif mode == "metric":
                    paths.append(self.ic_norm_q_path)
                else:
                    raise ValueError(f"unknown mode {mode!r}")
        return paths

    # --------------------------------------------------------------- internals

    def _check_paths(self, paths: List[str]) -> None:
        if len(paths) != len(self.schema):
            raise ValueError(
                f"paths length {len(paths)} != schema length {len(self.schema)}"
            )

    # ---- actuator I/O (per-actuator sub-comm collective MPI-IO) ----

    def _write_owned_actuator(self, local_arr, paths):
        slot_id = self.actuator_slot_owned
        slot_lo, _ = self._slot_range[slot_id]
        file_off_doubles = self.local_offset - slot_lo
        fh = MPI.File.Open(
            self.actuator_subcomm, paths[slot_id],
            MPI.MODE_WRONLY | MPI.MODE_CREATE,
        )
        fh.Set_view(0, etype=MPI.DOUBLE, filetype=MPI.DOUBLE, datarep="native")
        fh.Write_at_all(file_off_doubles, local_arr)
        fh.Close()

    def _read_owned_actuator(self, local_arr, paths):
        slot_id = self.actuator_slot_owned
        slot_lo, _ = self._slot_range[slot_id]
        file_off_doubles = self.local_offset - slot_lo
        fh = MPI.File.Open(self.actuator_subcomm, paths[slot_id], MPI.MODE_RDONLY)
        fh.Set_view(0, etype=MPI.DOUBLE, filetype=MPI.DOUBLE, datarep="native")
        fh.Read_at_all(file_off_doubles, local_arr)
        fh.Close()

    # ---- ic I/O (single-rank owner; PLOT3D solution) ----

    def _write_owned_ic(self, local_arr, paths):
        # In the parallel regime, each ic-owner owns exactly one slot whose
        # entire data lives in this rank's local_arr.
        slot_id = self.ic_owned[0]
        self._write_one_q(local_arr, paths[slot_id])

    def _read_owned_ic(self, local_arr, paths):
        slot_id = self.ic_owned[0]
        self._read_one_q(local_arr, paths[slot_id])

    # ---- serial fallback (rank 0 owns everything) ----

    def _serial_write(self, local_arr, paths):
        cursor = 0
        for i, s in enumerate(self.schema):
            slab = local_arr[cursor:cursor + s.size]
            if s.kind == "ic":
                self._write_one_q(slab, paths[i])
            else:
                with open(paths[i], "wb") as fh:
                    fh.write(slab.astype("<f8").tobytes())
            cursor += s.size

    def _serial_read(self, local_arr, paths):
        cursor = 0
        for i, s in enumerate(self.schema):
            slab_view = local_arr[cursor:cursor + s.size]
            if s.kind == "ic":
                self._read_one_q(slab_view, paths[i])
            else:
                with open(paths[i], "rb") as fh:
                    raw = fh.read(s.size * 8)
                if len(raw) != s.size * 8:
                    raise IOError(
                        f"{paths[i]}: expected {s.size * 8} bytes, got {len(raw)}"
                    )
                slab_view[:] = np.frombuffer(raw, dtype="<f8")
            cursor += s.size

    # ---- multi-block PLOT3D solution I/O via plot3dnasa ----

    def _require_plot3d(self):
        if plot3dnasa is None:
            raise RuntimeError(
                f"plot3dnasa import failed; cannot handle .q files: {_plot3d_import_error}"
            )
        if self.q_nblocks == 0 or self.q_nunknowns is None:
            raise RuntimeError(
                "ParallelIOHandler: .q geometry not cached (grid_file=None?)"
            )

    @staticmethod
    def _get_p3d_unknown_indices(nU):
        """File-slot indices that carry magudi data for a given nU.

        Matches plot3dWriteSingleSolution.f90:756-775 and
        plot3dReadSingleSolution.f90:1136-1167: magudi writes vars (1..nU-1)
        to PLOT3D slots (1..nU-1), advances past slots (nU..4) without
        writing, then writes var nU to slot 5. plot3dnasa.Solution.q[b] is
        always shape (Nx, Ny, Nz, 5); the inactive channel (only relevant
        for nU < 5) is padding that magudi ignores on read.
        """
        if nU == 5:
            return [0, 1, 2, 3, 4]
        return list(range(nU - 1)) + [4]

    def _write_one_q(self, flat_slab, path):
        self._require_plot3d()
        nU = self.q_nunknowns
        p3d_unknown_indices = self._get_p3d_unknown_indices(nU)
        sizes_2d = np.array([list(s) for s in self.grid_sizes], dtype=np.int32)
        sol = plot3dnasa.Solution()
        sol.set_size(sizes_2d, allocate=True)
        # ParallelIOHandler only supports OVERWRITING an existing .q file. The
        # PLOT3D aux header (timestep + time) is copied from that file so the
        # rewritten .q preserves it -- magudi reads auxiliaryData[0] as the
        # starting timestep and auxiliaryData[3] as the simulation time
        # (src/RegionImpl.f90:1381-1383), so a fresh zero header would put the
        # controller buffer at the wrong frame offset in subsequent segments.
        # Callers must pre-stage the target .q files (run_parallel.sh copies the
        # baseline snapshots into <prefix>-<k>.ic.q before magudi-optim starts).
        if not os.path.exists(path):
            raise RuntimeError(
                f"ParallelIOHandler._write_one_q requires {path} to exist so the "
                "PLOT3D aux header can be copied from it. Stage the baseline "
                ".q files before invoking the driver."
            )
        existing = plot3dnasa.FileFormat(path)
        if existing.aux_header is None:
            raise RuntimeError(
                f"{path}: expected a PLOT3D aux header to copy; got None."
            )
        sol._format.aux_header = existing.aux_header.copy()
        sol.time = float(existing.aux_header[-1])
        cursor = 0
        for b in range(self.q_nblocks):
            Nx, Ny, Nz = self.grid_sizes[b]
            block_size = Nx * Ny * Nz * nU
            sol.q[b][...] = 0.0  # ensure padding channel(s) carry a defined value
            block_data = flat_slab[cursor:cursor + block_size].reshape(
                (Nx, Ny, Nz, nU), order="F"
            )
            for u, c in enumerate(p3d_unknown_indices):
                sol.q[b][..., c] = block_data[..., u]
            cursor += block_size
        if cursor != len(flat_slab):
            raise RuntimeError(
                f"slab size {len(flat_slab)} != sum of block sizes {cursor}"
            )
        sol.save(path)

    def _read_one_q(self, flat_slab, path):
        self._require_plot3d()
        nU = self.q_nunknowns
        p3d_unknown_indices = self._get_p3d_unknown_indices(nU)
        sol = plot3dnasa.Solution(path, forceread=True)
        if sol.nblocks != self.q_nblocks:
            raise RuntimeError(
                f"{path}: nblocks {sol.nblocks} != expected {self.q_nblocks}"
            )
        cursor = 0
        for b in range(self.q_nblocks):
            Nx, Ny, Nz = self.grid_sizes[b]
            block_size = Nx * Ny * Nz * nU
            block = sol.q[b]
            if block.shape[:3] != (Nx, Ny, Nz):
                raise RuntimeError(
                    f"{path}: block {b} shape {block.shape[:3]} != expected ({Nx},{Ny},{Nz})"
                )
            block_data = np.empty((Nx, Ny, Nz, nU), order="F",
                                  dtype=flat_slab.dtype)
            for u, c in enumerate(p3d_unknown_indices):
                block_data[..., u] = block[..., c]
            flat_slab[cursor:cursor + block_size] = block_data.reshape(
                block_size, order="F"
            )
            cursor += block_size
        if cursor != len(flat_slab):
            raise RuntimeError(
                f"slab size {len(flat_slab)} != sum of block sizes {cursor}"
            )


# ---------------------------------------------------------------- self-test

def _selftest():
    """Round-trip a Vec through write + read for a synthetic schema.

    Stages layout.txt + norm_ic.q + a tiny YAML in a tempdir, then calls the
    cfg-driven ctor. Run with `mpirun -n N python3 parallel_io.py`.
    """
    import argparse
    import shutil
    import tempfile

    import yaml

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-actuator", type=int, default=1)
    parser.add_argument("--actuator-size", type=int, default=64)
    parser.add_argument("--nx", type=int, default=5)
    parser.add_argument("--ny", type=int, default=4)
    parser.add_argument("--nz", type=int, default=1)
    parser.add_argument("--nU", type=int, default=5)
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    workdir = tempfile.mkdtemp() if comm.Get_rank() == 0 else None
    workdir = comm.bcast(workdir, root=0)

    N = comm.Get_size()
    if N == 1:
        n_ic = 2
    else:
        n_ic = max(0, N - args.n_actuator)
        if args.n_actuator > 0 and (N - n_ic) % args.n_actuator != 0:
            raise SystemExit(
                f"selftest: with N={N} and n_actuator={args.n_actuator}, "
                f"choose --n-actuator dividing (N - n_ic)"
            )
    Nx, Ny, Nz, nU = args.nx, args.ny, args.nz, args.nU
    ic_size = Nx * Ny * Nz * nU
    schema = [FileSpec("actuator", f"act{i}", args.actuator_size)
              for i in range(args.n_actuator)]
    schema += [FileSpec("ic", str(k), ic_size) for k in range(1, n_ic + 1)]

    # Prefix the ctor will use. The ctor reads <prefix>.layout.txt and
    # <prefix>.norm_ic.q, both of which we stage below.
    selftest_prefix = os.path.join(workdir, "selftest")
    layout_path = f"{selftest_prefix}.layout.txt"
    ic_norm_q = f"{selftest_prefix}.norm_ic.q"
    yaml_path = os.path.join(workdir, "selftest.yml")
    log_path  = os.path.join(workdir, "selftest.optim.h5")

    if comm.Get_rank() == 0:
        # 1) layout.txt — same format parse_layout expects.
        with open(layout_path, "w") as fh:
            for s in schema:
                fh.write(f"{s.kind} {s.identifier} {s.size}\n")
        # 2) norm_ic.q — single-block PLOT3D solution for geometry caching.
        if n_ic > 0:
            if plot3dnasa is None:
                raise SystemExit(
                    f"selftest: plot3dnasa import failed: {_plot3d_import_error}"
                )
            sol = plot3dnasa.Solution()
            sol.set_size(np.array([[Nx, Ny, Nz]], dtype=np.int32), allocate=True)
            sol.q[0][:] = 1.0
            sol.save(ic_norm_q)
        # 3) YAML — minimal config the ctor will read.
        with open(yaml_path, "w") as fh:
            yaml.safe_dump({
                "global_prefix": selftest_prefix,
                "magudi": {
                    "forced_inputs": {
                        "time_splitting": {"number_of_segments": 2},
                    },
                },
                "optimization": {"log_file": log_path},
            }, fh)
    comm.Barrier()

    cfg = InputParser(yaml_path)
    io = ParallelIOHandler(cfg, comm)
    io.report_balance()

    v = io.create_vec()
    lo, hi = v.getOwnershipRange()
    v.getArray()[:] = np.arange(lo, hi, dtype="<f8")
    v.assemble()

    # Build paths -- .dat for actuator slots, .q for ic slots by default.
    paths = []
    for i, s in enumerate(schema):
        ext = "dat" if s.kind == "actuator" else "q"
        paths.append(os.path.join(workdir, f"{s.kind}_{s.identifier}.{ext}"))

    # _write_one_q requires the target file to exist so it can copy the aux
    # header. Pre-stage each ic slot path with a zero-valued .q at the right
    # geometry. Done on rank 0; barrier so every rank sees the files.
    if comm.Get_rank() == 0:
        for i, s in enumerate(schema):
            if s.kind != "ic":
                continue
            template = plot3dnasa.Solution()
            template.set_size(
                np.array([list(g) for g in io.grid_sizes], dtype=np.int32),
                allocate=True,
            )
            for b in range(io.q_nblocks):
                template.q[b][...] = 0.0
            template.save(paths[i])
    comm.Barrier()

    io.write(v, paths)

    v2 = io.read(paths)

    diff = v.duplicate()
    v.copy(diff)
    diff.axpy(-1.0, v2)
    err = diff.norm()
    PETSc.Sys.Print(f"round-trip |v - v2| = {err:.3e}")
    if comm.Get_rank() == 0:
        shutil.rmtree(workdir, ignore_errors=True)
    if err > 1e-12:
        raise SystemExit(f"FAIL: round-trip error {err:.3e}")
    PETSc.Sys.Print("OK: parallel_io self-test passed")


if __name__ == "__main__":
    _selftest()
