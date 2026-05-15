"""optim_ver4 MPI-parallel TAO driver.

Standalone parallel counterpart to optim.py. Run with:

    mpirun -n N_petsc python3 optim.parallel.py optim.parallel.yml [--max-iter K]

N_petsc is the rank count of this Python driver itself; PETSc Vec/Mat ops live
on its COMM_WORLD. The driver spawns ./forward and ./adjoint via MPI_Comm_spawn
at N_forward and N_adjoint ranks respectively (read from resource_distribution.
jobs.{forward,adjoint,petsc} in the YAML).

Sync mechanism: forward/adjoint call `disconnectParentIfSpawned` (see
include/MPIHelper.f90) immediately before MPI_Finalize. That makes the parent's
inter.Disconnect() return only after the child has done a collective
MPI_Comm_disconnect, which in turn happens only after all child MPI_File_close
calls have flushed. The same binaries work under plain mpirun -- get_parent
returns MPI_COMM_NULL and the disconnect is skipped.

The actuator .dat files (norm, control_forcing, gradient) are flat real64
streams. As long as control_space_norm, forward, and adjoint all run at the
same N_forward, the element ordering is consistent across the three files and
PETSc pointwise multiplies / norms are valid even though Python distributes
the data over a different rank count (N_petsc).

Do not change N_petsc between a save and a resume: history.petsc is laid out
according to the saving run's rank count.
"""
import argparse
import os
import sys
from collections import deque

import numpy as np
import yaml

# mpi4py MUST import before petsc4py so PETSc binds to MPI.COMM_WORLD.
from mpi4py import MPI  # noqa: F401  (import-order side effect)
from petsc4py import PETSc

LMVM_DEPTH = 5  # matches BQNLS default Max. storage = 5


def _spawn_and_wait(executable, args, n):
    """Collectively spawn `n` child MPI processes; block until they finalize.

    All N_petsc parent ranks must call this with identical arguments. Only one
    set of `n` children is launched (not N_petsc * n).

    Sync is via an inter-communicator `Barrier()`, NOT `Disconnect()`.
    Per the MPI standard, MPI_Comm_disconnect only waits for pending traffic on
    the inter-comm to finish; since we never send/recv on it, both sides
    disconnect immediately and no rendezvous occurs. MPI_Barrier on the
    inter-comm IS a true rendezvous over the union of parent+child groups.
    The child calls a matching MPI_Barrier(parent) from disconnectParentIfSpawned
    in MPIHelperImpl.f90, placed immediately before MPI_Comm_disconnect and
    MPI_Finalize -- so the parent's Barrier() returns only after the child has
    finished all its work.
    """
    inter = MPI.COMM_WORLD.Spawn(
        executable,
        args=args,
        maxprocs=n,
        info=MPI.INFO_NULL,
        root=0,
    )
    inter.Barrier()
    inter.Disconnect()
    # Parent-only barrier for filesystem cache visibility before the next
    # collective MPI-IO read of the child's output.
    MPI.COMM_WORLD.Barrier()


def parse_actuators(bc_path="./bc.dat"):
    """Return list of actuator patch names from bc.dat (rows with type ACTUATOR).

    Every rank parses independently; results are assumed identical. bc.dat is
    whitespace-delimited (Name Type ...); '#' marks an inline comment.
    """
    actuators = []
    with open(bc_path) as f:
        for line in f:
            line = line.split("#", 1)[0]
            tokens = line.split()
            if len(tokens) >= 2 and tokens[1] == "ACTUATOR":
                actuators.append(tokens[0])
    return actuators


def _read_flat_dat(path, total_size, comm, vec=None):
    """Collectively read a flat real64 file into a distributed PETSc Vec.

    The file is treated as a contiguous stream of `total_size` little-endian
    real64 values (the layout written by Fortran's MPI_File_write_all). Each
    rank reads its PETSc-default slab via MPI-IO collective Read_at_all.
    """
    if vec is None:
        vec = PETSc.Vec().createMPI(total_size, comm=comm)
    lo, hi = vec.getOwnershipRange()
    local_n = hi - lo
    buf = np.empty(local_n, dtype="<f8")
    fh = MPI.File.Open(comm, path, MPI.MODE_RDONLY)
    fh.Set_view(0, etype=MPI.DOUBLE, filetype=MPI.DOUBLE, datarep="native")
    fh.Read_at_all(lo, buf)
    fh.Close()
    vec.setArray(buf)
    vec.assemble()
    return vec


def _write_flat_dat(path, vec, comm):
    """Collectively write a distributed PETSc Vec to a flat real64 file."""
    lo, hi = vec.getOwnershipRange()
    arr = np.ascontiguousarray(vec.getArray(readonly=True), dtype="<f8")
    fh = MPI.File.Open(
        comm, path, MPI.MODE_WRONLY | MPI.MODE_CREATE
    )
    fh.Set_view(0, etype=MPI.DOUBLE, filetype=MPI.DOUBLE, datarep="native")
    fh.Write_at_all(lo, arr)
    fh.Close()


def save_tao_state(tao, checkpoint_path, history_path, snapshots, N, comm):
    """Persist iterate + replay-buffer + J0 diagonal collectively on `comm`.

    `snapshots` is a deque of (vx, vg) tuples of distributed PETSc Vecs (each
    on `comm`, size N). See optim.py for why J0 must be persisted explicitly --
    M.updateLMVM replay reconstructs (s_k, y_k) but not the nested
    MATLMVMDIAGBRDN diagonal.
    """
    x = tao.getSolution()
    viewer = PETSc.Viewer().createBinary(checkpoint_path, mode="w", comm=comm)
    x.view(viewer)
    viewer.destroy()
    PETSc.Sys.Print(f"Saved {checkpoint_path}: |y| = {x.norm():.6e}")

    viewer = PETSc.Viewer().createBinary(history_path, mode="w", comm=comm)
    # Header: 3-element distributed Vec, all elements pinned to rank 0 so the
    # layout is well-defined regardless of comm size.
    local_size = 3 if comm.Get_rank() == 0 else 0
    hdr = PETSc.Vec().createMPI((local_size, 3), comm=comm)
    if comm.Get_rank() == 0:
        hdr.setValues([0, 1, 2],
                      [float(tao.getIterationNumber()),
                       float(len(snapshots)),
                       float(N)])
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


def load_tao_state(tao, checkpoint_path, history_path, N, comm):
    """Resume iterate; replay L-BFGS history + J0 if present.

    Caller must have done tao.setSolution(y) AND tao.setUp() so the inner
    LMVM Mat is allocated and updateLMVM / setLMVMJ0 target the right object.
    """
    if not os.path.exists(checkpoint_path):
        PETSc.Sys.Print(
            f"No checkpoint at {checkpoint_path}; starting from current y"
        )
        return False
    x = tao.getSolution()
    viewer = PETSc.Viewer().createBinary(checkpoint_path, mode="r", comm=comm)
    x_loaded = PETSc.Vec().load(viewer)
    viewer.destroy()
    if x_loaded.getSize() != x.getSize():
        raise SystemExit(
            f"{checkpoint_path} size {x_loaded.getSize()} "
            f"!= expected {x.getSize()}"
        )
    x_loaded.copy(x)
    x_loaded.destroy()
    PETSc.Sys.Print(
        f"Resumed iterate from {checkpoint_path}: |y| = {x.norm():.6e}"
    )

    if not (os.path.exists(history_path) and os.path.getsize(history_path) > 0):
        PETSc.Sys.Print(
            f"No L-BFGS history at {history_path}; using fresh L-BFGS"
        )
        return True

    viewer = PETSc.Viewer().createBinary(history_path, mode="r", comm=comm)
    hdr = PETSc.Vec().load(viewer)
    # Header is laid out with all 3 elements on rank 0; allgather to make all
    # ranks see (iter_no, n_snap, n_dim) consistently. comm is a PETSc.Comm;
    # tompi4py() returns the underlying mpi4py communicator for allgather().
    local_vals = list(hdr.getArray(readonly=True))
    gathered = comm.tompi4py().allgather(local_vals)
    flat = [v for sub in gathered for v in sub]
    iter_no, n_snap, n_dim = int(flat[0]), int(flat[1]), int(flat[2])
    hdr.destroy()
    if n_dim != N:
        raise SystemExit(f"history N={n_dim} != current N={N}")
    M = tao.getLMVMMat()
    for _ in range(n_snap):
        vx = PETSc.Vec().load(viewer)
        vg = PETSc.Vec().load(viewer)
        M.updateLMVM(vx, vg)
        vx.destroy()
        vg.destroy()
    J0_diag = PETSc.Vec().load(viewer)
    viewer.destroy()
    M.setLMVMJ0(PETSc.Mat().createDiagonal(J0_diag))
    tao.setIterationNumber(iter_no)
    PETSc.Sys.Print(
        f"Replayed {n_snap} snapshots + J0 diag; iter counter set to {iter_no}"
    )
    return True


def _read_jobs_procs(cfg, key):
    """Extract `procs` (index 1) from resource_distribution.jobs[<key>] = [N,P]."""
    try:
        nodes_procs = cfg["resource_distribution"]["jobs"][key]
    except (KeyError, TypeError) as exc:
        raise SystemExit(
            f"resource_distribution.jobs.{key} missing from config"
        ) from exc
    if not (isinstance(nodes_procs, list) and len(nodes_procs) == 2):
        raise SystemExit(
            f"resource_distribution.jobs.{key} must be [nodes, procs]; "
            f"got {nodes_procs!r}"
        )
    return int(nodes_procs[1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", help="Path to optim.parallel.yml")
    parser.add_argument(
        "--max-iter", type=int, default=5,
        help="iterations this run (added to loaded iter if resuming)"
    )
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    pcomm = PETSc.COMM_WORLD

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    N_forward = _read_jobs_procs(cfg, "forward")
    N_adjoint = _read_jobs_procs(cfg, "adjoint")
    N_petsc = _read_jobs_procs(cfg, "petsc")
    if comm.Get_size() != N_petsc:
        raise SystemExit(
            f"mpirun -n mismatch: launched with {comm.Get_size()} ranks but "
            f"resource_distribution.jobs.petsc requests {N_petsc}"
        )

    prefix = cfg["global_prefix"]
    actuators = parse_actuators("./bc.dat")
    if len(actuators) != 1:
        raise SystemExit(
            "optim.parallel.py supports a single actuator; "
            f"got {len(actuators)} from ./bc.dat: {actuators}"
        )
    actuator = actuators[0]

    norm_path = f"{prefix}.norm_{actuator}.dat"
    control_path = f"{prefix}.control_forcing_{actuator}.dat"
    gradient_path = f"{prefix}.gradient_{actuator}.dat"
    j_path = f"{prefix}.forward_run.txt"
    checkpoint_path = "y.petsc"
    history_path = "history.petsc"

    # Total size of M_diag in real64 elements -- rank 0 stats, then bcast.
    if comm.Get_rank() == 0:
        size_bytes = os.path.getsize(norm_path)
        if size_bytes == 0 or size_bytes % 8 != 0:
            N_total = -1
        else:
            N_total = size_bytes // 8
    else:
        N_total = None
    N_total = comm.bcast(N_total, root=0)
    if N_total <= 0:
        raise SystemExit(f"{norm_path} is empty or not a multiple of 8 bytes")

    M_diag = _read_flat_dat(norm_path, N_total, comm)
    # Collective non-positive check.
    local_min = float(M_diag.getArray(readonly=True).min()) if M_diag.getLocalSize() > 0 else float("inf")
    global_min = comm.allreduce(local_min, op=MPI.MIN)
    if global_min <= 0.0:
        raise SystemExit(
            f"{norm_path} has non-positive entries (min={global_min}); "
            "cannot form sqrt(M)."
        )

    # D_inv = 1 / sqrt(M_diag) on COMM_WORLD.
    D_inv = M_diag.duplicate()
    M_diag.copy(D_inv)
    D_inv.sqrtabs()
    D_inv.reciprocal()

    # Optimization iterate and gradient template, both distributed.
    y = PETSc.Vec().createMPI(N_total, comm=pcomm)
    y.set(0.0)
    g_y_template = y.duplicate()
    g_y_template.set(0.0)

    iter_log = []
    snapshots = deque(maxlen=LMVM_DEPTH + 1)
    fg_count = [0]
    n_spawn = [0]
    pending = [None]

    def fg(tao, y_vec, g_vec):
        # x = y * D_inv  (distributed); write to control_forcing.dat
        x_dist = y_vec.duplicate()
        x_dist.pointwiseMult(y_vec, D_inv)
        _write_flat_dat(control_path, x_dist, comm)
        x_dist.destroy()
        comm.Barrier()

        _spawn_and_wait("./forward", ["--input", "magudi.inp"], N_forward)
        n_spawn[0] += 1
        _spawn_and_wait("./adjoint", ["--input", "magudi.inp"], N_adjoint)
        n_spawn[0] += 1

        # J: scalar, read on rank 0, bcast.
        if comm.Get_rank() == 0:
            with open(j_path) as fh:
                J_local = float(fh.read().strip())
        else:
            J_local = None
        J = comm.bcast(J_local, root=0)

        # Gradient: collective read, then g_y = raw_grad * D_inv.
        raw_grad = _read_flat_dat(gradient_path, N_total, comm)
        if raw_grad.getSize() != N_total:
            raise SystemExit(
                f"gradient size {raw_grad.getSize()} != expected {N_total}"
            )
        g_vec.pointwiseMult(raw_grad, D_inv)
        raw_grad.destroy()
        fg_count[0] += 1
        iter_log.append(J)
        return J

    def monitor(tao):
        it = tao.getIterationNumber()
        try:
            its, f_val, gnorm, cnorm, xdiff, reason = tao.getSolutionStatus()
        except AttributeError:
            f_val = iter_log[-1] if iter_log else float("nan")
            gnorm = float("nan")
        PETSc.Sys.Print(
            f"  iter {it:4d}  J = {f_val: .6e}  |g_y| = {gnorm: .6e}"
        )
        # Deferred-commit snapshot capture: monitor fires N+1 times for N TAO
        # iters; drop the trailing fire because its (x_N, g_N) was never fed
        # to MatLMVMUpdate. See optim.py:229-240 for the full rationale.
        if pending[0] is not None:
            snapshots.append(pending[0])
        x_now = tao.getSolution()
        g_now = tao.getGradient()[0]
        vx = x_now.duplicate(); x_now.copy(vx)
        vg = g_now.duplicate(); g_now.copy(vg)
        pending[0] = (vx, vg)

    tao = PETSc.TAO().create(comm=pcomm)
    tao.setType("bqnls")
    tao.setObjectiveGradient(fg, g_y_template)
    tao.setMaximumIterations(args.max_iter)
    tao.setTolerances(grtol=1.0e-6)
    PETSc.Options().setValue("tao_recycle_history", True)
    try:
        tao.setMonitor(monitor)
    except (AttributeError, TypeError):
        PETSc.Options().setValue("tao_monitor", "")
    tao.setFromOptions()

    PETSc.Sys.Print(
        f"optim.parallel: prefix={prefix} actuator={actuator} N={N_total} "
        f"N_petsc={N_petsc} N_forward={N_forward} N_adjoint={N_adjoint}"
    )

    tao.setSolution(y)
    tao.setUp()
    resumed = load_tao_state(tao, checkpoint_path, history_path, N_total, pcomm)
    if not resumed:
        # First run: pre-seed a zero control_forcing.dat so the first TAO
        # callback's ./forward finds the file. Collective write.
        zero_x = y.duplicate(); zero_x.set(0.0)
        _write_flat_dat(control_path, zero_x, comm)
        zero_x.destroy()
        comm.Barrier()

    tao.solve(y)

    save_tao_state(tao, checkpoint_path, history_path, snapshots, N_total, pcomm)

    reason = tao.getConvergedReason()
    final_J = iter_log[-1] if iter_log else float("nan")
    initial_J = iter_log[0] if iter_log else float("nan")
    PETSc.Sys.Print("")
    PETSc.Sys.Print(f"Converged reason: {reason}")
    PETSc.Sys.Print(f"Initial J = {initial_J: .6e}")
    PETSc.Sys.Print(f"Final   J = {final_J: .6e}")
    PETSc.Sys.Print(
        f"fg calls = {fg_count[0]};  spawn calls = {n_spawn[0]}"
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
