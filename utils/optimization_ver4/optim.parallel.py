"""optim_ver4 Phase 3 MPI-parallel TAO driver.

Wraps `./msforward` and `./msadjoint` (bin/msforward.f90, bin/msadjoint.f90)
inside a TAO L-BFGS outer loop. Layout + parallel I/O live in
ParallelIOHandler (parallel_io.py); this file owns the optimizer state,
spawn coordination, and resume/checkpoint logic.

Run:
    mpirun -n N_petsc python3 optim.parallel.py optim.parallel.yml \\
                                                [--max-iter K]

N_petsc satisfies the ParallelIOHandler policy:
  N_petsc == 1  (serial fallback)
or
  N_petsc = n_ic + n_actuator * nprocs_per_actuator     (nprocs_per_actuator >= 1)

The schema is parsed from <prefix>.layout.txt (written by bin/compute_norm).
M_diag is loaded from per-actuator .norm_<name>.dat plus the shared
<prefix>.norm_ic.q file (the same file fed to every ic slot of M_diag).

State persists across allocations the same way as the strawman:
  y.petsc        TAO iterate in y-space
  history.petsc  replay-buffer (snapshot Vecs) + J0 diagonal
Do not change N_petsc between a save and a resume: the layout (and thus
the Vec's local sizes per rank) is layout-bound.
"""
import argparse
import os
import shlex
import sys
from collections import deque

# mpi4py MUST import before petsc4py so PETSc binds to MPI.COMM_WORLD.
from mpi4py import MPI  # noqa: F401  (import-order side effect)
from petsc4py import PETSc

from inputs import InputParser
from parallel_io import ParallelIOHandler, parse_layout

LMVM_DEPTH = 5  # matches BQNLS default Max. storage = 5


def _spawn_and_wait(executable, args, n):
    """Collectively spawn `n` child MPI processes; block until they finalize.

    All N_petsc parent ranks must call this with identical arguments. Only one
    set of `n` children is launched (not N_petsc * n).

    Child stdout/stderr go to ./out/<basename(executable)>.out; rank 0
    truncates the file before each call. Suppression uses
    /bin/bash -c "exec <bin> ... >> log 2>&1": the shell exec replaces
    itself with the real binary, preserving the PID and PMI env vars so
    MPI_Init in the child still completes the spawn handshake.

    Sync is via an inter-communicator `Barrier()`, NOT `Disconnect()`.
    Per the MPI standard, MPI_Comm_disconnect only waits for pending traffic
    on the inter-comm to finish; since we never send/recv on it, both sides
    disconnect immediately and no rendezvous occurs. MPI_Barrier on the
    inter-comm IS a true rendezvous over the union of parent+child groups.
    The child calls a matching MPI_Barrier(parent) from disconnectParentIfSpawned
    in MPIHelperImpl.f90, placed immediately before MPI_Comm_disconnect and
    MPI_Finalize -- so the parent's Barrier() returns only after the child
    has finished all its work.
    """
    comm = MPI.COMM_WORLD
    log_path = os.path.abspath(
        os.path.join("out", os.path.basename(executable) + ".out")
    )
    if comm.Get_rank() == 0:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        open(log_path, "w").close()
    comm.Barrier()

    quoted_args = " ".join(shlex.quote(a) for a in args)
    redir_cmd = (
        f"exec {shlex.quote(executable)} {quoted_args} "
        f">> {shlex.quote(log_path)} 2>&1"
    )

    inter = MPI.COMM_WORLD.Spawn(
        "/bin/bash",
        args=["-c", redir_cmd],
        maxprocs=n,
        info=MPI.INFO_NULL,
        root=0,
    )
    inter.Barrier()
    inter.Disconnect()
    MPI.COMM_WORLD.Barrier()


def save_tao_state(tao, checkpoint_path, history_path, snapshots, N, comm):
    """Persist iterate + replay-buffer + J0 diagonal collectively on `comm`.

    `snapshots` is a deque of (vx, vg) tuples of distributed PETSc Vecs
    with the file-aligned layout. See optim.py for why J0 must be persisted
    explicitly -- M.updateLMVM replay reconstructs (s_k, y_k) but not the
    nested MATLMVMDIAGBRDN diagonal.
    """
    x = tao.getSolution()
    viewer = PETSc.Viewer().createBinary(checkpoint_path, mode="w", comm=comm)
    x.view(viewer)
    viewer.destroy()
    PETSc.Sys.Print(f"Saved {checkpoint_path}: |y| = {x.norm():.6e}")

    viewer = PETSc.Viewer().createBinary(history_path, mode="w", comm=comm)
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


def _seed_zero_actuator_files(schema, prefix, comm):
    """Rank 0 creates zero-filled .control_forcing_<name>.dat for each actuator
    slot, sized to its layout entry. Idempotent: skipped when the file already
    exists at the expected size. Needed only on the first run, before the very
    first io.read for the baseline x -- subsequent iterations overwrite the
    file via io.write.
    """
    if comm.Get_rank() == 0:
        for s in schema:
            if s.kind != "actuator":
                continue
            path = f"{prefix}.control_forcing_{s.identifier}.dat"
            expected = s.size * 8
            actual = os.path.getsize(path) if os.path.exists(path) else -1
            if actual != expected:
                with open(path, "wb") as fh:
                    fh.write(b"\x00" * expected)
    comm.Barrier()


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

    cfg = InputParser(args.config)
    N_forward = int(cfg.getInput(
        ["resource_distribution", "jobs", "forward"], datatype=list)[1])
    N_adjoint = int(cfg.getInput(
        ["resource_distribution", "jobs", "adjoint"], datatype=list)[1])
    N_petsc = int(cfg.getInput(
        ["resource_distribution", "jobs", "petsc"], datatype=list)[1])
    if comm.Get_size() != N_petsc:
        raise SystemExit(
            f"mpirun -n mismatch: launched with {comm.Get_size()} ranks but "
            f"resource_distribution.jobs.petsc requests {N_petsc}"
        )

    prefix = cfg.getInput(["global_prefix"], datatype=str)
    layout_path = f"{prefix}.layout.txt"
    ic_norm_q_path = f"{prefix}.norm_ic.q"
    checkpoint_path = "y.petsc"
    history_path = "history.petsc"
    j_path = f"{prefix}.forward_run.txt"
    gg_path = f"{prefix}.adjoint_run.txt"

    schema = parse_layout(layout_path)
    io = ParallelIOHandler(schema, prefix, ic_norm_q_path, comm=comm)
    io.report_balance()

    N_total = io.global_size

    # Load M_diag (per-actuator norm files + shared .norm_ic.q) and form D_inv.
    M_diag = io.read_metric()
    local_min = (
        float(M_diag.getArray(readonly=True).min())
        if M_diag.getLocalSize() > 0 else float("inf")
    )
    global_min = comm.allreduce(local_min, op=MPI.MIN)
    if global_min <= 0.0:
        raise SystemExit(
            f"M_diag has non-positive entries (min={global_min}); cannot form sqrt(M). "
            f"Check controller mollifier coverage and state_controllability."
        )
    D_inv = M_diag.duplicate()
    M_diag.copy(D_inv)
    D_inv.sqrtabs()
    D_inv.reciprocal()

    # Optimization iterate and gradient template.
    y = io.create_vec()
    g_y_template = y.duplicate(); g_y_template.set(0.0)

    iter_log = []
    snapshots = deque(maxlen=LMVM_DEPTH + 1)
    fg_count = [0]
    n_spawn = [0]
    pending = [None]

    def fg(tao, y_vec, g_vec):
        # Wipe per-iteration artifacts from the previous msforward/msadjoint
        # run so this evaluation starts on a clean filesystem (barrier inside).
        io.cleanup_iteration_artifacts()

        x_dist = y_vec.duplicate()
        x_dist.pointwiseMult(y_vec, D_inv)
        io.write_x(x_dist)
        x_dist.destroy()
        comm.Barrier()

        _spawn_and_wait("./msforward", ["--input", "magudi.inp"], N_forward)
        n_spawn[0] += 1
        _spawn_and_wait("./msadjoint", ["--input", "magudi.inp"], N_adjoint)
        n_spawn[0] += 1

        # J: scalar, read on rank 0, bcast. msforward folds the matching-condition
        # penalty into the total before writing (bin/msforward.f90:237).
        # gg: <g,g>_M scalar from msadjoint (.sub_adjoint_run.txt sum).
        if comm.Get_rank() == 0:
            with open(j_path) as fh:
                J_local = float(fh.read().strip())
            with open(gg_path) as fh:
                gg_local = float(fh.read().strip())
            print("fwd eval: J = %.5e   <g,g>_M = %.5e" % (J_local, gg_local),
                  flush=True)
        else:
            J_local = None
            gg_local = None
        J = comm.bcast(J_local, root=0)
        comm.bcast(gg_local, root=0)

        raw_grad = io.read_grad()
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
        f"optim.parallel: prefix={prefix} N_total={N_total} "
        f"n_actuator={io.n_actuator} n_ic={io.n_ic} "
        f"N_petsc={N_petsc} N_forward={N_forward} N_adjoint={N_adjoint}"
    )

    tao.setSolution(y)
    tao.setUp()
    resumed = load_tao_state(tao, checkpoint_path, history_path, N_total, pcomm)
    if not resumed:
        # Initial x = (zero control, baseline IC). The baseline IC slabs were staged
        # from a controller-off forward run by run_parallel.sh as <prefix>-<k>.ic.q.
        # .control_forcing_<name>.dat doesn't exist yet -- rank 0 seeds zero files
        # so io.read can stitch them into baseline_x. Then y_init = D * x_init.
        _seed_zero_actuator_files(schema, prefix, comm)
        baseline_x = io.read_x()
        y.pointwiseDivide(baseline_x, D_inv)   # y = baseline_x / D_inv = D * baseline_x
        baseline_x.destroy()
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
