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
import subprocess
import sys

import numpy as np
import yaml
from petsc4py import PETSc

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

    # Write zero controlForcing.dat so any forward call (including the first TAO
    # callback's, before x.tofile() runs) finds the file. TAO callbacks overwrite it.
    # This is strawman-example-specific command.
    np.zeros_like(M_diag).tofile(control_path)

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
    try:
        tao.setMonitor(monitor)
    except (AttributeError, TypeError):
        # petsc4py monitor signature varies across versions; fall back to the
        # built-in PETSc monitor via options.
        PETSc.Options().setValue("tao_monitor", "")
    tao.setFromOptions()

    PETSc.Sys.Print(
        f"Phase 1 strawman: prefix={prefix} actuator={actuator} N={M_diag.size}"
    )
    tao.solve(y)

    reason = tao.getConvergedReason()
    final_J = iter_log[-1] if iter_log else float("nan")
    initial_J = iter_log[0] if iter_log else float("nan")
    PETSc.Sys.Print("")
    PETSc.Sys.Print(f"Converged reason: {reason}")
    PETSc.Sys.Print(f"Initial J = {initial_J: .6e}")
    PETSc.Sys.Print(f"Final   J = {final_J: .6e}")
    PETSc.Sys.Print(f"Iterations evaluated: {len(iter_log)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
