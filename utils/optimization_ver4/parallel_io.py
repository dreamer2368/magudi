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

import os
import sys
from collections import namedtuple
from typing import List, Optional

import numpy as np

from mpi4py import MPI  # noqa: F401  (import-order side effect)
from petsc4py import PETSc

# plot3dnasa is typically copied into the staged run directory (= cwd) by the
# run wrapper. Add cwd to sys.path so the import finds it without requiring
# PYTHONPATH manipulation in shell.
if os.getcwd() not in sys.path:
    sys.path.insert(0, os.getcwd())

try:
    import plot3dnasa  # noqa: F401  (probed for .q handling)
    _plot3d_import_error = None
except Exception as exc:  # pragma: no cover - exercised at runtime via .q path
    plot3dnasa = None
    _plot3d_import_error = exc


FileSpec = namedtuple("FileSpec", ["kind", "identifier", "size"])
# kind ∈ {"actuator", "ic"} — semantic role in the layout, not wire format.


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

    def __init__(self, schema: List[FileSpec], ic_norm_q_path: Optional[str],
                 comm: MPI.Comm = MPI.COMM_WORLD):
        """ic_norm_q_path: path to <prefix>.norm_ic.q (compute_norm's output).
        Inspected once at init to learn (nblocks, per-block sizes, nUnknowns).
        Optional if the schema has no ic slots.
        """
        if not schema:
            raise ValueError("ParallelIOHandler: empty schema")
        self.schema = list(schema)
        self.comm = comm
        self.rank = comm.Get_rank()
        self.size = comm.Get_size()

        self._classify_slots()
        self._assign_ownership()
        self._build_dat_subcomm()
        self._cache_q_geometry(ic_norm_q_path)
        self._validate_ic_geometry()

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

    def read(self, vec: PETSc.Vec, paths: List[str]) -> None:
        self._check_paths(paths)
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

    Run with `mpirun -n N python3 parallel_io.py [--layout ...] [--ic-norm-q ...]`
    """
    import argparse
    import shutil
    import tempfile

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--layout", default=None,
                        help="Optional path to layout.txt; synthetic schema otherwise.")
    parser.add_argument("--ic-norm-q", default=None,
                        help="Optional path to compute_norm's .norm_ic.q for geometry; "
                             "synthetic single-block .q otherwise.")
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

    if args.layout is not None:
        schema = parse_layout(args.layout)
        ic_norm_q = args.ic_norm_q
    else:
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
        # Synthesize a .norm_ic.q-like file (rank 0) so geometry caching works
        # without depending on a .xyz file.
        if n_ic > 0 and args.ic_norm_q is None:
            ic_norm_q = os.path.join(workdir, "synth.norm_ic.q")
            if comm.Get_rank() == 0:
                if plot3dnasa is None:
                    raise SystemExit(
                        f"selftest: plot3dnasa import failed: {_plot3d_import_error}"
                    )
                sol = plot3dnasa.Solution()
                sol.set_size(np.array([[Nx, Ny, Nz]], dtype=np.int32), allocate=True)
                sol.q[0][:] = 1.0
                sol.save(ic_norm_q)
            comm.Barrier()
        else:
            ic_norm_q = args.ic_norm_q

    io = ParallelIOHandler(schema, ic_norm_q, comm)
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

    io.write(v, paths)

    v2 = io.create_vec()
    io.read(v2, paths)

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
