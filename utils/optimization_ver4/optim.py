"""optim_ver4 Phase 1 strawman driver.

Drives PETSc/TAO BQNLS (Bound-constrained Quasi-Newton with Line Search; the
L-BFGS-style replacement for the deprecated LMVM solver) on the OneDWave
example with Nsplit=1 (single time window). The objective and gradient come
from the existing `forward` and `adjoint` binaries via subprocess; the
variable transformation `y = D x` (with D = diag(sqrt(M_diag))) lives in
this callback.

Pin to a single MPI rank: the actuator `.dat` files written by the Fortran
exes use MPI_File_write_all with per-rank slabs concatenated. Phase 1 keeps
np=1 so numpy.fromfile reads the file as a single contiguous real64 buffer.
"""
import argparse
import os
import subprocess
import sys
from collections import deque

import numpy as np
import yaml
from petsc4py import PETSc

LMVM_DEPTH = 5  # matches BQNLS default Max. storage = 5

def _run_silent(cmd):
    """Run a subprocess silently. On failure, print its captured output and raise."""
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        sys.stderr.write(f"FAILED ({exc.returncode}): {' '.join(exc.cmd)}\n")
        if exc.stdout:
            sys.stderr.write(exc.stdout)
        if exc.stderr:
            sys.stderr.write(exc.stderr)
        raise

def parse_actuators(bc_path="./bc.dat"):
    """Return list of actuator patch names from bc.dat (rows with type ACTUATOR).

    bc.dat is whitespace-delimited with columns:
        Name Type Grid normDir iMin iMax jMin jMax kMin kMax
    Lines starting with '#' (or with '#' anywhere — comment marker) are skipped.
    """
    actuators = []
    with open(bc_path) as f:
        for line in f:
            line = line.split("#", 1)[0]
            tokens = line.split()
            if len(tokens) >= 2 and tokens[1] == "ACTUATOR":
                actuators.append(tokens[0])
    return actuators


def _make_capture(snapshots):
    """Return a `tao.setUpdate` callback that appends (x_k, g_k) snapshots.

    setUpdate fires once per accepted TAO iter, matching the count of
    TAO's internal MatLMVMUpdate calls — capturing (x_0, g_0) through
    (x_{N-1}, g_{N-1}). Monitor-based capture would over-collect by one
    (it also fires at iter N AFTER the last MatLMVMUpdate); the extra
    replay call forms a spurious BFGS pair on resume.

    Arrays are copied because the Vec storage aliases TAO's internals
    and mutates on the next iteration.
    """
    def capture(tao, _it):
        x_vec = tao.getSolution()
        g_vec, _ = tao.getGradient()
        snapshots.append((
            x_vec.getArray(readonly=True).copy(),
            g_vec.getArray(readonly=True).copy(),
        ))
    return capture


def save_tao_state(tao, checkpoint_path, history_path, snapshots, N):
    """Persist iterate, replay-buffer snapshots, and J0 diagonal.

    Replay of (x_k, g_k) reconstructs the BFGS (s_k, y_k) pairs by
    construction, but NOT the diagonal J0 (initial-Hessian rescale)
    that MATLMVMBFGS maintains via a nested MATLMVMDIAGBRDN matrix.
    M.updateLMVM called externally on a fresh Mat enters that nested
    update path from a different state than during TAO's own solve,
    so J0 diverges even with identical (x, g) inputs. Persisting J0's
    diagonal Vec and reseating it via setLMVMJ0 on resume restores
    bitwise equivalence (verified in toy_optim.py: rel diff = 0).
    """
    x = tao.getSolution()
    viewer = PETSc.Viewer().createBinary(checkpoint_path, mode="w", comm=x.getComm())
    x.view(viewer)
    viewer.destroy()
    PETSc.Sys.Print(f"Saved {checkpoint_path}: |y| = {x.norm():.6e}")

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
    """Resume iterate; replay L-BFGS history + J0 if present. Returns True iff iterate loaded.

    Caller must have done `tao.setSolution(y)` AND `tao.setUp()` so the
    inner LMVM Mat is allocated and `M.updateLMVM` / `M.setLMVMJ0`
    target the right object.
    """
    if not os.path.exists(checkpoint_path):
        PETSc.Sys.Print(f"No checkpoint at {checkpoint_path}; starting from current y")
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
    PETSc.Sys.Print(f"Resumed iterate from {checkpoint_path}: |y| = {x.norm():.6e}")

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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", help="Path to optim.single.yml")
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    prefix = cfg["global_prefix"]
    actuators = parse_actuators("./bc.dat")
    if len(actuators) != 1:
        raise SystemExit(
            "Phase 1 strawman supports a single actuator; "
            f"got {len(actuators)} from ./bc.dat: {actuators}"
        )
    actuator = actuators[0]

    # Filenames must match the defaults built in:
    #   src/ActuatorPatchImpl.f90:53-54  -> <prefix>.gradient_<patch>.dat
    #   src/ActuatorPatchImpl.f90:69-70  -> <prefix>.control_forcing_<patch>.dat
    #   bin/control_space_norm.f90:214   -> <prefix>.norm_<patch>.dat
    #   bin/Forward.f90:96               -> <prefix>.forward_run.txt
    norm_path = f"{prefix}.norm_{actuator}.dat"
    control_path = f"{prefix}.control_forcing_{actuator}.dat"
    gradient_path = f"{prefix}.gradient_{actuator}.dat"
    j_path = f"{prefix}.forward_run.txt"
    checkpoint_path = "y.petsc"      # TAO iterate, in y-coordinates (DESIGN.md §4)
    history_path = "history.petsc"   # replay buffer + J0 diagonal

    M_diag = np.fromfile(norm_path, dtype="<f8")
    if M_diag.size == 0:
        raise SystemExit(f"{norm_path} is empty; run control_space_norm first.")
    if np.any(M_diag <= 0.0):
        raise SystemExit(f"{norm_path} has non-positive entries; cannot form sqrt(M).")
    D_inv = 1.0 / np.sqrt(M_diag)

    # Backing arrays: PETSc Vecs alias these via createWithArray, so they must
    # outlive tao.solve(). Declare at module scope of main() to keep them alive.
    y_arr = np.zeros_like(M_diag)
    g_arr = np.zeros_like(M_diag)
    y = PETSc.Vec().createWithArray(y_arr, comm=PETSc.COMM_SELF)
    g_y_template = PETSc.Vec().createWithArray(g_arr, comm=PETSc.COMM_SELF)

    iter_log = []
    snapshots = deque(maxlen=LMVM_DEPTH + 1)

    def fg(tao, y_vec, g_vec):
        # TAO marks the input vec read-only inside the callback; getArray() with
        # the default readonly=False would fail with VecSetErrorIfLocked.
        x = y_vec.getArray(readonly=True) * D_inv
        x.tofile(control_path)

        _run_silent(["./forward", "--input", "magudi.inp"])
        _run_silent(["./adjoint", "--input", "magudi.inp"])

        with open(j_path) as f:
            J = float(f.read().strip())

        raw_grad = np.fromfile(gradient_path, dtype="<f8")
        if raw_grad.size != M_diag.size:
            raise SystemExit(
                f"gradient size {raw_grad.size} != M_diag size {M_diag.size}"
            )
        g_vec.setArray(raw_grad * D_inv)
        iter_log.append(J)
        return J

    def monitor(tao):
        it = tao.getIterationNumber()
        try:
            its, f_val, gnorm, cnorm, xdiff, reason = tao.getSolutionStatus()
        except AttributeError:
            f_val = iter_log[-1] if iter_log else float("nan")
            gnorm = float("nan")
        PETSc.Sys.Print(f"  iter {it:4d}  J = {f_val: .6e}  |g_y| = {gnorm: .6e}")

    tao = PETSc.TAO().create(comm=PETSc.COMM_SELF)
    tao.setType("bqnls")
    tao.setObjectiveGradient(fg, g_y_template)
    tao.setMaximumIterations(5)
    tao.setTolerances(grtol=1.0e-6)
    PETSc.Options().setValue("tao_recycle_history", True)
    try:
        tao.setMonitor(monitor)
    except (AttributeError, TypeError):
        # petsc4py monitor signature varies across versions; fall back to the
        # built-in PETSc monitor via options.
        PETSc.Options().setValue("tao_monitor", "")
    tao.setUpdate(_make_capture(snapshots))
    tao.setFromOptions()

    PETSc.Sys.Print(
        f"Phase 1 strawman: prefix={prefix} actuator={actuator} N={M_diag.size}"
    )

    # Bind y so tao.getSolution() returns it inside the save/load helpers,
    # then setUp so the inner LMVM Mat is allocated before load_tao_state
    # replays into it.
    tao.setSolution(y)
    tao.setUp()
    resumed = load_tao_state(tao, checkpoint_path, history_path, M_diag.size)
    if not resumed:
        # First run: pre-seed a zero controlForcing.dat so the first TAO
        # callback's ./forward call finds the file. Subsequent callbacks
        # (and resumed runs) overwrite it from x = y * D_inv.
        np.zeros_like(M_diag).tofile(control_path)

    tao.solve(y)

    save_tao_state(tao, checkpoint_path, history_path, snapshots, M_diag.size)

    reason = tao.getConvergedReason()
    final_J = iter_log[-1] if iter_log else float("nan")
    initial_J = iter_log[0] if iter_log else float("nan")
    PETSc.Sys.Print("")
    PETSc.Sys.Print(f"Converged reason: {reason}")
    PETSc.Sys.Print(f"Initial J = {initial_J: .6e}")
    PETSc.Sys.Print(f"Final   J = {final_J: .6e}")
    PETSc.Sys.Print(f"Iterations evaluated: {len(iter_log)}")
    print(iter_log)

    return 0


if __name__ == "__main__":
    sys.exit(main())
