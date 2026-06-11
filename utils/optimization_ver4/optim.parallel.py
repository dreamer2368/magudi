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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", help="Path to optim.parallel.yml")
    parser.add_argument(
        "--max-iter", type=int, default=5,
        help="iterations this run (added to loaded iter if resuming)"
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="print extra progress information"
    )
    parser.add_argument(
        "--debug", action="store_true",
        help="print debug-level diagnostics"
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

    # Optimization tolerance: TAO's relative gradient-norm test (grtol).
    grtol = cfg.getInput(["optimization", "tolerance"], fallback=1.0e-8)
    # Line-search overrides. Each fallback matches the TAO documented default
    # so that an absent YAML entry leaves the line search at its built-in
    # configuration (setting a PETSc option to its default is a no-op).
    ls_type       = cfg.getInput(
        ["optimization", "line_search", "type"], fallback="more-thuente")
    ls_init_step  = cfg.getInput(
        ["optimization", "line_search", "initial_step_size"], fallback=1.0)
    ls_step_max   = cfg.getInput(
        ["optimization", "line_search", "max_step"], fallback=1.0e+15)
    ls_step_min   = cfg.getInput(
        ["optimization", "line_search", "min_step"], fallback=1.0e-20)
    ls_stol       = cfg.getInput(
        ["optimization", "line_search", "bracket_tol"], fallback=1.0e-4)
    ls_max_funcs  = cfg.getInput(
        ["optimization", "line_search", "max_funcs"], fallback=30)
    ls_log_file   = cfg.getInput(
        ["optimization", "line_search", "log_file"],
        fallback=f"{prefix}.line_search.h5")

    schema = parse_layout(layout_path)
    io = ParallelIOHandler(schema, prefix, ic_norm_q_path, comm=comm, pcomm=pcomm)
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
    # Preconditioning: y = D * x with D = M^(1/2), so that the y-space
    # Hessian D_inv * H_x * D_inv ≈ I (assuming H_x ≈ M). msadjoint returns
    # the M-Riesz gradient g (defined by J(x+ah) ≈ J(x) + a * g^T M h), so
    # the chain-rule gives g_y = D * g. Mind the asymmetry between the
    # iterate (x = D_inv * y) and the gradient (g_y = D * g).
    D = M_diag.duplicate()
    M_diag.copy(D)
    D.sqrtabs()                 # D = sqrt(M_diag) = M^(1/2)
    D_inv = D.duplicate()
    D.copy(D_inv)
    D_inv.reciprocal()          # D_inv = D^(-1) = M^(-1/2)

    # Optimization iterate and gradient template.
    y = io.create_vec()
    g_y_template = y.duplicate(); g_y_template.set(0.0)

    iter_log = []
    snapshots = deque(maxlen=LMVM_DEPTH + 1)
    fg_count = [0]
    n_spawn = [0]
    pending = [None]
    # Diagnostic state for the line-search bracket. ls_iter_seen tracks the
    # outer-iteration counter at which we last captured the bracket start;
    # ls_bracket holds (y_init, g_init) -- copies of y_vec / g_vec at the
    # bracket start, used to reconstruct alpha*d := y_vec - y_init at trials.
    ls_iter_seen = [-1]
    ls_bracket = [None]

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

        # J: msforward folds the matching-condition penalty into the total
        # before writing (bin/msforward.f90:237).
        # gg: <g,g>_M scalar from msadjoint (.sub_adjoint_run.txt sum).
        J = io.read_scalar(j_path)
        gg_local = io.read_scalar(gg_path)

        raw_grad = io.read_grad()
        # raw_grad is the M-Riesz gradient g (msadjoint convention); the
        # y-space gradient is g_y = D * g, NOT D_inv * g. See the comment
        # at the D / D_inv definition above for the chain-rule derivation.
        g_vec.pointwiseMult(raw_grad, D)
        raw_grad.destroy()

        if args.debug:
            # DIAGNOSTIC (1): cross-check g_y · g_y == g^T M g == gg_local.
            # The y-space Euclidean norm of g_y must equal the M-norm of the
            # M-Riesz gradient that msadjoint reported. A mismatch flags an
            # error in the D-scaling, the layout, or the M_diag file.
            gg_check = g_vec.dot(g_vec)
            rel = abs(gg_check - gg_local) / max(abs(gg_local), 1e-30)
            PETSc.Sys.Print(
                f"fwd eval #{fg_count[0]+1:4d}  J = {J: .5e}  "
                f"<g_y,g_y> = {gg_check: .5e}  vs gg_local = {gg_local: .5e}  "
                f"rel = {rel: .2e}"
            )

            # DIAGNOSTIC (2): reconstruct (Dg)^T d at each line-search trial.
            # petsc4py does not expose the LS direction or step length, so we
            # detect the bracket start via the outer-iteration counter (which
            # increments only after a LS accepts) and store (y_init, g_init).
            # At every subsequent trial in the same LS, y_vec - y_init = alpha*d
            # so g_vec . (y_vec - y_init) = alpha * (g_now . d) and similarly
            # g_init . (y_vec - y_init) = alpha * dginit. The ratio cancels
            # alpha, giving dg_now / dginit -- the strong-Wolfe gate (must
            # shrink below gtol = 0.9 in magnitude for the LS to accept).
            it = tao.getIterationNumber()
            if it != ls_iter_seen[0]:
                ls_iter_seen[0] = it
                if ls_bracket[0] is not None:
                    ls_bracket[0][0].destroy()
                    ls_bracket[0][1].destroy()
                y_init = y_vec.duplicate(); y_vec.copy(y_init)
                g_init = g_vec.duplicate(); g_vec.copy(g_init)
                ls_bracket[0] = (y_init, g_init)
                PETSc.Sys.Print(f"    [LS#{it}] bracket start captured")
            else:
                y_init, g_init = ls_bracket[0]
                dy = y_vec.duplicate(); y_vec.copy(dy); dy.axpy(-1.0, y_init)
                ginit_dot_dy = g_init.dot(dy)          # = alpha * dginit
                gnow_dot_dy  = g_vec.dot(dy)           # = alpha * dg(alpha)
                dy_norm = dy.norm()
                g_init_norm = g_init.norm()
                dy.destroy()
                # cos(dy, g_init) should be -1 at iter 0 (pure steepest descent
                # along the initial scaled-identity Hessian). Deviation flags
                # a non-steepest direction (e.g. stale L-BFGS history active).
                denom = dy_norm * g_init_norm
                if denom > 0.0:
                    cos_dy_g0 = ginit_dot_dy / denom
                    cos_str = f"cos(dy,g0) = {cos_dy_g0: .6f}"
                else:
                    cos_str = "cos(dy,g0) = NaN"
                if abs(ginit_dot_dy) > 0.0:
                    ratio = gnow_dot_dy / ginit_dot_dy
                    ratio_str = f"|dg/dg0| = {abs(ratio):.3f}"
                else:
                    ratio_str = "|dg/dg0| = NaN"
                PETSc.Sys.Print(
                    f"    [LS#{it} trial] |alpha*d| = {dy_norm: .3e}  "
                    f"alpha*g0.d = {ginit_dot_dy: .3e}  "
                    f"alpha*g.d  = {gnow_dot_dy: .3e}  {ratio_str}  {cos_str}"
                )

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
    tao.setTolerances(grtol=grtol)
    PETSc.Options().setValue("tao_recycle_history", True)
    # Canonical option names per PETSc/src/tao/linesearch/interface/taolinesearch.c:497-503.
    # NOTE: stepinit / stepmax / stepmin are one word (no underscore); the
    # bracket relative tolerance is tao_ls_rtol (tao_ls_stol does not exist in
    # this PETSc release).
    PETSc.Options().setValue("tao_ls_type", ls_type)
    PETSc.Options().setValue("tao_ls_stepinit", ls_init_step)
    PETSc.Options().setValue("tao_ls_stepmax",  ls_step_max)
    PETSc.Options().setValue("tao_ls_stepmin",  ls_step_min)
    PETSc.Options().setValue("tao_ls_rtol",     ls_stol)
    PETSc.Options().setValue("tao_ls_max_funcs", ls_max_funcs)
    # Per-trial line-search log -- writes "LS step <alpha> f <f(alpha)>" for
    # every Wolfe trial inside an outer iteration. Useful for diagnosing -6
    # (TAO_DIVERGED_LS_FAILURE) without re-running. `:append` keeps records
    # across the multiple optim.parallel.py invocations that run_parallel.sh
    # chains, instead of truncating on each invocation.
    PETSc.Options().setValue("tao_ls_monitor", ls_log_file)
    # Surface any options the user set that PETSc never consumed (typically
    # a spelling mistake against the PETSc option registry).
    PETSc.Options().setValue("options_left", True)
    try:
        tao.setMonitor(monitor)
    except (AttributeError, TypeError):
        PETSc.Options().setValue("tao_monitor", "")
    tao.setFromOptions()

    # Dump the resolved TAO config (type, tolerances, line search subtype +
    # its parameters) so the run log records the actual setting picked up from
    # the YAML / command-line overlay. Equivalent to passing -tao_view.
    if args.verbose:
        PETSc.Sys.Print("--- TAO configuration ---")
        tao.view()
        PETSc.Sys.Print("-------------------------")

    PETSc.Sys.Print(
        f"optim.parallel: prefix={prefix} N_total={N_total} "
        f"n_actuator={io.n_actuator} n_ic={io.n_ic} "
        f"N_petsc={N_petsc} N_forward={N_forward} N_adjoint={N_adjoint}"
    )

    tao.setSolution(y)
    tao.setUp()
    resumed = io.load_tao_state(tao, checkpoint_path, history_path)
    if not resumed:
        # Initial x = (zero control, baseline IC). The baseline IC slabs were staged
        # from a controller-off forward run by run_parallel.sh as <prefix>-<k>.ic.q.
        # .control_forcing_<name>.dat doesn't exist yet -- rank 0 seeds zero files
        # so io.read can stitch them into baseline_x. Then y_init = D * x_init.
        io.seed_zero_actuator_files()
        baseline_x = io.read_x()
        y.pointwiseMult(baseline_x, D)         # y_init = D * baseline_x
        baseline_x.destroy()
        comm.Barrier()

    tao.solve(y)

    io.save_tao_state(tao, checkpoint_path, history_path, snapshots)

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
