"""Toy driver for investigating TAO BQNLS restart behavior.

Replaces magudi's forward/adjoint subprocess pair with an in-process numpy
Rosenbrock function so the save_tao_state / load_tao_state path can be
exercised in isolation, with no PLOT3D files, no bc.dat, and no MPI.

Usage:
    python3 toy_optim.py --max-iter 5                  # first run, makes checkpoint
    python3 toy_optim.py --max-iter 5                  # resumes from checkpoint

    python3 toy_optim.py --max-iter 10 --no-resume     # baseline: 10 iters in one shot

Compare the (iter, J, |g|) trajectories: with working L-BFGS history
persistence, the resumed run's iter-k state should match the baseline's
iter-(prev+k) state.
"""
import argparse
import os
import sys

import numpy as np
from petsc4py import PETSc


def rosenbrock(x):
    """N-D Rosenbrock f(x) = sum_i [100 (x_{i+1} - x_i^2)^2 + (1 - x_i)^2].

    Returns (J, grad). Minimum is 0 at x = [1, ..., 1].
    """
    x = np.asarray(x, dtype=np.float64)
    a = x[:-1]
    b = x[1:]
    J = np.sum(100.0 * (b - a * a) ** 2 + (1.0 - a) ** 2)

    g = np.zeros_like(x)
    g[:-1] += -400.0 * a * (b - a * a) - 2.0 * (1.0 - a)
    g[1:]  +=  200.0 * (b - a * a)
    return float(J), g


def _report_lmvm_storage(M, N):
    """Empirically probe the LMVM Mat's storage footprint.

    Question we're answering: is MATLMVMBFGS a dense N x N matrix (infeasible
    at magudi-scale N), or a shell matrix backed by ~2m vectors of length N
    (cheap)?
    """
    PETSc.Sys.Print("-" * 70)
    PETSc.Sys.Print(f"LMVM Mat: type = {M.getType()!r}")
    PETSc.Sys.Print(f"  getSize()       = {M.getSize()}     (logical N x N)")
    PETSc.Sys.Print(f"  getLocalSize()  = {M.getLocalSize()}")
    try:
        info = M.getInfo()
        PETSc.Sys.Print(f"  getInfo()       = {info}")
        nz = info.get("nz_allocated", float("nan"))
        mem = info.get("memory", float("nan"))
        PETSc.Sys.Print(
            f"  -> nz_allocated = {nz:g}   memory = {mem:g} bytes"
        )
        PETSc.Sys.Print(
            f"  Dense N x N would need {N*N*8} bytes "
            f"({N*N*8/1024/1024:.2f} MiB) for N={N}."
        )
    except Exception as exc:
        PETSc.Sys.Print(f"  getInfo() failed: {exc!r}")
    for name in ("getLMVMUpdateCount", "getUpdateCount"):
        if hasattr(M, name):
            try:
                PETSc.Sys.Print(f"  M.{name}() = {getattr(M, name)()}")
            except Exception as exc:
                PETSc.Sys.Print(f"  M.{name}() raised {exc!r}")
    PETSc.Sys.Print("-" * 70)


def save_tao_state(tao, checkpoint_path, history_path):
    """Persist TAO iterate and (attempted) L-BFGS curvature pairs.

    Mirrors optim.py: the iterate write works; the MatLMVMBFGS view is a
    silent no-op in PETSc 3.25, kept here so behavior matches optim.py.
    """
    x = tao.getSolution()
    viewer = PETSc.Viewer().createBinary(checkpoint_path, mode="w", comm=x.getComm())
    x.view(viewer)
    viewer.destroy()
    PETSc.Sys.Print(f"Saved {checkpoint_path}: |x| = {x.norm():.6e}")

    M = tao.getLMVMMat()
    _report_lmvm_storage(M, x.getSize())

    viewer = PETSc.Viewer().createBinary(history_path, mode="w", comm=x.getComm())
    M.view(viewer)
    viewer.destroy()
    on_disk = os.path.getsize(history_path) if os.path.exists(history_path) else 0
    PETSc.Sys.Print(
        f"Attempted L-BFGS save to {history_path} (size on disk: {on_disk} bytes)"
    )


def load_tao_state(tao, checkpoint_path, history_path):
    """Resume iterate, and L-BFGS curvature if a non-empty history file exists.

    Returns True iff iterate was loaded. Caller must have done
    tao.setSolution(x) and tao.setUp() so the inner LMVM Mat is allocated.
    """
    if not os.path.exists(checkpoint_path):
        PETSc.Sys.Print(f"No checkpoint at {checkpoint_path}; starting fresh")
        return False
    x = tao.getSolution()
    viewer = PETSc.Viewer().createBinary(checkpoint_path, mode="r", comm=x.getComm())
    x_loaded = PETSc.Vec().load(viewer)
    viewer.destroy()
    if x_loaded.getSize() != x.getSize():
        raise SystemExit(
            f"{checkpoint_path} size {x_loaded.getSize()} != expected {x.getSize()}"
        )
    x_loaded.copy(x)
    x_loaded.destroy()
    PETSc.Sys.Print(f"Resumed from {checkpoint_path}: |x| = {x.norm():.6e}")

    if os.path.exists(history_path) and os.path.getsize(history_path) > 0:
        M = tao.getLMVMMat()
        viewer = PETSc.Viewer().createBinary(history_path, mode="r", comm=x.getComm())
        M.load(viewer)
        viewer.destroy()
        PETSc.Sys.Print(f"Resumed L-BFGS history from {history_path}")
    else:
        PETSc.Sys.Print(f"No usable L-BFGS history at {history_path}; using fresh L-BFGS")
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dim", type=int, default=10, help="problem dimension N")
    parser.add_argument("--max-iter", type=int, default=5, help="TAO max iterations this run")
    parser.add_argument("--grtol", type=float, default=1.0e-8, help="TAO gradient rtol")
    parser.add_argument("--checkpoint", default="toy_x.petsc", help="iterate checkpoint path")
    parser.add_argument("--history", default="toy_history.petsc", help="L-BFGS history path")
    parser.add_argument("--no-resume", action="store_true",
                        help="ignore (and remove) any existing checkpoint/history")
    parser.add_argument("--x0", type=float, default=-1.2,
                        help="initial guess: x[0]=value, x[1:]=1 (classic Rosenbrock start)")
    args = parser.parse_args()

    if args.no_resume:
        for p in (args.checkpoint, args.history):
            if os.path.exists(p):
                os.remove(p)
                PETSc.Sys.Print(f"Removed {p}")

    N = args.dim
    # Backing arrays kept alive in main()'s scope: PETSc Vecs alias them via
    # createWithArray, so they must outlive tao.solve().
    x_arr = np.ones(N, dtype=np.float64)
    x_arr[0] = args.x0
    g_arr = np.zeros(N, dtype=np.float64)
    x = PETSc.Vec().createWithArray(x_arr, comm=PETSc.COMM_SELF)
    g_template = PETSc.Vec().createWithArray(g_arr, comm=PETSc.COMM_SELF)

    iter_log = []

    def fg(tao, x_vec, g_vec):
        x_np = x_vec.getArray(readonly=True)
        J, grad = rosenbrock(x_np)
        g_vec.setArray(grad)
        iter_log.append(J)
        return J

    def monitor(tao):
        it = tao.getIterationNumber()
        try:
            _, f_val, gnorm, _, _, _ = tao.getSolutionStatus()
        except AttributeError:
            f_val = iter_log[-1] if iter_log else float("nan")
            gnorm = float("nan")
        PETSc.Sys.Print(f"  iter {it:4d}  J = {f_val: .6e}  |g| = {gnorm: .6e}")

    tao = PETSc.TAO().create(comm=PETSc.COMM_SELF)
    tao.setType("bqnls")
    tao.setObjectiveGradient(fg, g_template)
    tao.setMaximumIterations(args.max_iter)
    tao.setTolerances(grtol=args.grtol)
    PETSc.Options().setValue("tao_recycle_history", True)
    try:
        tao.setMonitor(monitor)
    except (AttributeError, TypeError):
        PETSc.Options().setValue("tao_monitor", "")
    tao.setFromOptions()

    PETSc.Sys.Print(
        f"Toy optim: N={N} max_iter={args.max_iter} "
        f"checkpoint={args.checkpoint} history={args.history}"
    )

    tao.setSolution(x)
    tao.setUp()

    load_tao_state(tao, args.checkpoint, args.history)

    tao.solve(x)

    save_tao_state(tao, args.checkpoint, args.history)

    reason = tao.getConvergedReason()
    final_J = iter_log[-1] if iter_log else float("nan")
    initial_J = iter_log[0] if iter_log else float("nan")
    PETSc.Sys.Print("")
    PETSc.Sys.Print(f"Converged reason: {reason}")
    PETSc.Sys.Print(f"Initial J = {initial_J: .6e}")
    PETSc.Sys.Print(f"Final   J = {final_J: .6e}")
    PETSc.Sys.Print(f"f-evals this run: {len(iter_log)}")
    PETSc.Sys.Print(f"|x - 1| = {np.linalg.norm(x.getArray(readonly=True) - 1.0):.6e}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
