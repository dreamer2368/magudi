"""Toy driver for TAO BQNLS restart via replay-buffer L-BFGS history.

PETSc 3.25's `Mat.view`/`Mat.load` is a silent no-op for MATLMVMBFGS, so
the curvature pairs cannot be persisted via the natural Mat binary path.
Two things need to be reconstructed for bitwise-equivalent restart:

  1. The (s_k, y_k) BFGS curvature pairs.
     Capture (x_k, g_k) at iterate boundaries via `tao.setUpdate` (which
     fires once per accepted iter, matching TAO's internal MatLMVMUpdate
     count). On resume, call `M.updateLMVM(x_k, g_k)` for each saved pair
     in order — this reproduces the same (s_k, y_k) by construction since
     LMVM forms s = x - x_prev, y = g - g_prev from successive inputs.

  2. The J0 (initial-Hessian) diagonal rescale.
     MATLMVMBFGS maintains J0 via a nested MATLMVMDIAGBRDN matrix whose
     update path doesn't reproduce externally even with identical (x, g)
     inputs. Without restoring J0, M_b's H * g diverges from M_ref's by
     ~90% — even when the pair count matches. Fix: persist the J0
     diagonal Vec alongside the snapshots and reseat it via
     `M.setLMVMJ0(PETSc.Mat().createDiagonal(j0_vec))` on resume.

With both pieces restored, `--verify-restart` confirms rel diff = 0 on a
probe-vector H * g comparison and max |dJ| = |dg| = 0 across the full
post-restart trajectory.

Usage:
    # Normal use: run, checkpoint, resume
    python3 toy_optim.py --max-iter 5                 # first run
    python3 toy_optim.py --max-iter 5                 # resumes
    python3 toy_optim.py --max-iter 10 --no-resume    # baseline (one shot)

    # Self-contained restart-equivalence test (no leftover files)
    python3 toy_optim.py --verify-restart --dim 50 --total-iters 12 --split-at 8
"""
import argparse
import os
import sys
from collections import deque

import numpy as np
from petsc4py import PETSc

LMVM_DEPTH = 5  # matches BQNLS default Max. storage = 5


def rosenbrock(x):
    """N-D Rosenbrock f(x) = sum_i [100 (x_{i+1} - x_i^2)^2 + (1 - x_i)^2]."""
    x = np.asarray(x, dtype=np.float64)
    a = x[:-1]
    b = x[1:]
    J = np.sum(100.0 * (b - a * a) ** 2 + (1.0 - a) ** 2)
    g = np.zeros_like(x)
    g[:-1] += -400.0 * a * (b - a * a) - 2.0 * (1.0 - a)
    g[1:]  +=  200.0 * (b - a * a)
    return float(J), g


def save_tao_state(tao, checkpoint_path, history_path, snapshots, N):
    """Persist iterate, replay-buffer snapshots, and J0 diagonal.

    Replay of (x_k, g_k) reconstructs the BFGS (s_k, y_k) pairs by
    construction, but NOT the diagonal J0 (initial-Hessian rescale)
    that PETSc's MATLMVMBFGS maintains via a nested MATLMVMDIAGBRDN
    matrix. M.updateLMVM called externally on a fresh Mat enters that
    nested update path from a different state than during TAO's own
    solve, so J0 diverges even with identical (x, g) inputs. Persisting
    J0's diagonal Vec and reseating it via setLMVMJ0 on resume restores
    bitwise equivalence (rel diff = 0 vs a non-restarted run).
    """
    x = tao.getSolution()
    viewer = PETSc.Viewer().createBinary(checkpoint_path, mode="w", comm=x.getComm())
    x.view(viewer)
    viewer.destroy()
    PETSc.Sys.Print(f"Saved {checkpoint_path}: |x| = {x.norm():.6e}")

    viewer = PETSc.Viewer().createBinary(history_path, mode="w", comm=PETSc.COMM_SELF)
    hdr = PETSc.Vec().createWithArray(
        np.array([tao.getIterationNumber(), len(snapshots), N], dtype="<f8"),
        comm=PETSc.COMM_SELF,
    )
    hdr.view(viewer)
    hdr.destroy()
    for x_arr, g_arr in snapshots:
        vx = PETSc.Vec().createWithArray(x_arr, comm=PETSc.COMM_SELF)
        vg = PETSc.Vec().createWithArray(g_arr, comm=PETSc.COMM_SELF)
        vx.view(viewer)
        vg.view(viewer)
        vx.destroy()
        vg.destroy()
    J0_diag = tao.getLMVMMat().getLMVMJ0().getDiagonal()
    J0_diag.view(viewer)
    J0_diag.destroy()
    viewer.destroy()
    PETSc.Sys.Print(
        f"Saved {history_path}: {len(snapshots)} snapshots + J0 diag, "
        f"iter={tao.getIterationNumber()}"
    )


def load_tao_state(tao, checkpoint_path, history_path, N):
    """Resume iterate; replay L-BFGS history if present. Returns True iff iterate loaded."""
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
    PETSc.Sys.Print(f"Resumed iterate from {checkpoint_path}: |x| = {x.norm():.6e}")

    if not (os.path.exists(history_path) and os.path.getsize(history_path) > 0):
        PETSc.Sys.Print(f"No L-BFGS history at {history_path}; using fresh L-BFGS")
        return True

    viewer = PETSc.Viewer().createBinary(history_path, mode="r", comm=PETSc.COMM_SELF)
    hdr = PETSc.Vec().load(viewer)
    iter_no, n_snap, n_dim = (int(v) for v in hdr.getArray(readonly=True))
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


def _build_tao(N, max_iter, x0_value, snapshots=None, log=None,
               fg_count=None, grtol=1.0e-8):
    """Construct a BQNLS TAO with capture/monitor wired in.

    Returns (tao, x_vec). The numpy arrays backing the createWithArray
    Vecs need no explicit keepalive — petsc4py's createWithArray holds a
    Python reference to the underlying numpy buffer for the Vec's lifetime.

    `fg_count`, if provided, is a 1-element list that fg increments on
    every call. Lets the caller distinguish outer iters from total
    objective+gradient evaluations (the setUpdate-triggered
    TaoComputeObjective after each iter adds one extra).
    """
    x_arr = np.ones(N, dtype=np.float64)
    x_arr[0] = x0_value
    g_arr = np.zeros(N, dtype=np.float64)
    x = PETSc.Vec().createWithArray(x_arr, comm=PETSc.COMM_SELF)
    g_template = PETSc.Vec().createWithArray(g_arr, comm=PETSc.COMM_SELF)

    def fg(_tao, x_vec, g_vec):
        J, grad = rosenbrock(x_vec.getArray(readonly=True))
        g_vec.setArray(grad)
        if fg_count is not None:
            fg_count[0] += 1
        return J

    # Deferred-commit slot for snapshot capture. Monitor fires N+1 times for
    # N TAO iterations (once at niter=0 with the initial iterate via
    # bnk.c:72, then once per accepted iter after ++tao->niter at
    # bnls.c:167). We want the first N captures (x_0..x_{N-1}) which match
    # TAO's internal MatLMVMUpdate inputs; the final fire's (x_N, g_N) was
    # never fed to any MatLMVMUpdate. Defer committing each fire's data
    # until the NEXT fire confirms it wasn't the last; the unconfirmed
    # trailing entry is discarded when solve() returns.
    pending = [None]

    def monitor(tao):
        it = tao.getIterationNumber()
        try:
            _, f_val, gnorm, *_ = tao.getSolutionStatus()
        except AttributeError:
            f_val = float("nan")
            gnorm = float("nan")
        if log is not None:
            log.append((it, f_val, gnorm))
        PETSc.Sys.Print(f"  iter {it:4d}  J = {f_val: .6e}  |g| = {gnorm: .6e}")
        if snapshots is not None:
            if pending[0] is not None:
                snapshots.append(pending[0])
            pending[0] = (
                tao.getSolution().getArray(readonly=True).copy(),
                tao.getGradient()[0].getArray(readonly=True).copy(),
            )

    tao = PETSc.TAO().create(comm=PETSc.COMM_SELF)
    tao.setType("bqnls")
    tao.setObjectiveGradient(fg, g_template)
    tao.setMaximumIterations(max_iter)
    tao.setTolerances(grtol=grtol)
    PETSc.Options().setValue("tao_recycle_history", True)
    try:
        tao.setMonitor(monitor)
    except (AttributeError, TypeError):
        PETSc.Options().setValue("tao_monitor", "")
    tao.setFromOptions()
    return tao, x


def verify_restart(N, total_iters, split_at, x0):
    """Two-layer equivalence test for replay-based restart.

    Layer A: H_k probe — show that the LMVM operator's action on a fixed
    probe vector matches between a baseline TAO paused at iter `split_at`
    and a fresh TAO whose LMVM was reconstructed via replay.

    Layer B: trajectory match — run full save -> file -> load wiring and
    assert that the resumed run's (J, |g|) at iter (split_at + k) matches
    the baseline's at the same iter, bitwise.
    """
    PETSc.Sys.Print(
        f"verify_restart: N={N} total_iters={total_iters} split_at={split_at}"
    )

    # ------------------------------------------------------------------
    # Layer A: H * g_probe equivalence
    # ------------------------------------------------------------------
    PETSc.Sys.Print("--- Layer A: H * g_probe ---")
    snap_ref = deque(maxlen=LMVM_DEPTH + 1)
    tao_ref, x_ref = _build_tao(N, split_at, x0, snapshots=snap_ref)
    tao_ref.setSolution(x_ref)
    tao_ref.setUp()
    tao_ref.solve(x_ref)
    M_ref = tao_ref.getLMVMMat()
    x_at_split = x_ref.getArray(readonly=True).copy()

    J0_ref_diag = M_ref.getLMVMJ0().getDiagonal()

    tao_b, x_b = _build_tao(N, 1, x0)
    x_b.setArray(x_at_split)
    tao_b.setSolution(x_b)
    tao_b.setUp()
    M_b = tao_b.getLMVMMat()
    for x_arr, g_arr in snap_ref:
        vx = PETSc.Vec().createWithArray(x_arr, comm=PETSc.COMM_SELF)
        vg = PETSc.Vec().createWithArray(g_arr, comm=PETSc.COMM_SELF)
        M_b.updateLMVM(vx, vg)
        vx.destroy()
        vg.destroy()
    M_b.setLMVMJ0(PETSc.Mat().createDiagonal(J0_ref_diag))
    tao_b.setIterationNumber(split_at)

    probe = np.arange(N, dtype=np.float64) + 1.0
    g_ref = PETSc.Vec().createWithArray(probe.copy(), comm=PETSc.COMM_SELF)
    g_b = PETSc.Vec().createWithArray(probe.copy(), comm=PETSc.COMM_SELF)
    d_ref = M_ref.createVecLeft()
    d_b = M_b.createVecLeft()
    M_ref.mult(g_ref, d_ref)
    M_b.mult(g_b, d_b)
    diff = d_ref.duplicate()
    d_ref.copy(diff)
    diff.axpy(-1.0, d_b)
    rel = diff.norm() / max(d_ref.norm(), 1e-30)
    layer_a_ok = rel < 1e-12
    PETSc.Sys.Print(
        f"  |H_ref g| = {d_ref.norm():.6e}  |H_b g| = {d_b.norm():.6e}  "
        f"rel diff = {rel:.2e}  {'OK' if layer_a_ok else 'FAIL'}"
    )

    # ------------------------------------------------------------------
    # Layer B: full save/load trajectory match
    # ------------------------------------------------------------------
    PETSc.Sys.Print("--- Layer B: post-restart trajectory ---")
    baseline_log = []
    tao_base, x_base = _build_tao(N, total_iters, x0, log=baseline_log)
    tao_base.setSolution(x_base)
    tao_base.setUp()
    tao_base.solve(x_base)

    snap_split = deque(maxlen=LMVM_DEPTH + 1)
    tao_a, x_a = _build_tao(N, split_at, x0, snapshots=snap_split)
    tao_a.setSolution(x_a)
    tao_a.setUp()
    tao_a.solve(x_a)
    ckpt = ".verify_x.petsc"
    hist = ".verify_history.petsc"
    save_tao_state(tao_a, ckpt, hist, snap_split, N)

    split_log = []
    tao_c, x_c = _build_tao(N, total_iters, x0, log=split_log)
    tao_c.setSolution(x_c)
    tao_c.setUp()
    load_tao_state(tao_c, ckpt, hist, N)
    tao_c.solve(x_c)

    n_match = min(len(baseline_log) - split_at, len(split_log))
    dJ_max = 0.0
    dg_max = 0.0
    for k in range(n_match):
        _, J_base, gnorm_base = baseline_log[split_at + k]
        _, J_split, gnorm_split = split_log[k]
        dJ_max = max(dJ_max, abs(J_base - J_split))
        dg_max = max(dg_max, abs(gnorm_base - gnorm_split))
    layer_b_ok = (n_match > 0) and dJ_max < 1e-10 and dg_max < 1e-10
    PETSc.Sys.Print(
        f"  matched {n_match} post-restart iters  "
        f"max |dJ| = {dJ_max:.2e}  max |dg| = {dg_max:.2e}  "
        f"{'OK' if layer_b_ok else 'FAIL'}"
    )

    for p in (ckpt, hist):
        if os.path.exists(p):
            os.remove(p)

    return 0 if (layer_a_ok and layer_b_ok) else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dim", type=int, default=10, help="problem dimension N")
    parser.add_argument("--max-iter", type=int, default=5,
                        help="iterations this run (added to loaded iter if resuming)")
    parser.add_argument("--grtol", type=float, default=1.0e-8)
    parser.add_argument("--checkpoint", default="toy_x.petsc")
    parser.add_argument("--history", default="toy_history.petsc")
    parser.add_argument("--no-resume", action="store_true",
                        help="remove any existing checkpoint/history first")
    parser.add_argument("--x0", type=float, default=-1.2,
                        help="initial guess: x[0]=value, x[1:]=1")
    parser.add_argument("--verify-restart", action="store_true",
                        help="run the two-layer equivalence test and exit")
    parser.add_argument("--total-iters", type=int, default=12,
                        help="(verify only) baseline total iterations")
    parser.add_argument("--split-at", type=int, default=8,
                        help="(verify only) iteration at which to restart")
    args = parser.parse_args()

    if args.verify_restart:
        return verify_restart(
            N=args.dim, total_iters=args.total_iters,
            split_at=args.split_at, x0=args.x0,
        )

    if args.no_resume:
        for p in (args.checkpoint, args.history):
            if os.path.exists(p):
                os.remove(p)
                PETSc.Sys.Print(f"Removed {p}")

    snapshots = deque(maxlen=LMVM_DEPTH + 1)
    fg_count = [0]
    tao, x = _build_tao(
        args.dim, args.max_iter, args.x0,
        snapshots=snapshots, fg_count=fg_count, grtol=args.grtol,
    )
    PETSc.Sys.Print(
        f"Toy optim: N={args.dim} max_iter={args.max_iter} "
        f"checkpoint={args.checkpoint} history={args.history}"
    )
    tao.setSolution(x)
    tao.setUp()

    load_tao_state(tao, args.checkpoint, args.history, args.dim)
    # TAO BQNLS resets niter to 0 inside tao.solve() regardless of
    # setIterationNumber, so args.max_iter is exactly this run's budget.

    tao.solve(x)
    save_tao_state(tao, args.checkpoint, args.history, snapshots, args.dim)

    PETSc.Sys.Print("")
    PETSc.Sys.Print(f"Converged reason: {tao.getConvergedReason()}")
    PETSc.Sys.Print(f"fg calls = {fg_count[0]}")
    PETSc.Sys.Print(
        f"|x - 1| = {np.linalg.norm(x.getArray(readonly=True) - 1.0):.6e}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
