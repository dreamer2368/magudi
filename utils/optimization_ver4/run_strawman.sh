#!/usr/bin/env bash
# optim_ver4 Phase 1 strawman: stage OneDWave with Nsplit=1, run TAO L-BFGS.
#
# Pinned to np=1: the actuator .dat files (controlForcing, gradient, norm)
# are MPI_File_write_all output. With np=1 they are a single contiguous
# real64 buffer that numpy.fromfile reads directly. Changing the rank count
# requires reworking optim.py's I/O.
#
# Run from the build/ directory after `make -j` has produced bin/forward,
# bin/adjoint, bin/control_space_norm.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$(pwd)}"
WORK_DIR="${BUILD_DIR}/OneDWave_strawman"

if [ ! -x "${BUILD_DIR}/bin/forward" ]; then
    echo "error: ${BUILD_DIR}/bin/forward not found; run cmake + make first" >&2
    exit 1
fi

rm -rf "${WORK_DIR}"
mkdir -p "${WORK_DIR}"
cp -r "${REPO_ROOT}/examples/OneDWave/." "${WORK_DIR}/"
cp "${REPO_ROOT}/utils/python/plot3dnasa.py" "${WORK_DIR}/"

cd "${WORK_DIR}"

python3 config.py

ln -sf "${BUILD_DIR}/bin/forward" .
ln -sf "${BUILD_DIR}/bin/adjoint" .
ln -sf "${BUILD_DIR}/bin/control_space_norm" .

# Patch magudi.inp for the Nsplit=1 strawman:
#   save_interval=480                  -> matches adjoint_save_timesteps for Nsplit=1
#   baseline_prediction_available=true -> ./adjoint reuses the just-written forward snapshots
# controller_switch is left at its initial value (false) for control_space_norm,
# which doesn't read controlForcing.dat. We flip it to true just before optim.py.
sed -i 's/^save_interval = .*/save_interval = 480/' magudi.inp
sed -i 's/^baseline_prediction_available = .*/baseline_prediction_available = true/' magudi.inp

# control_space_norm walks the grid to compute the diagonal of the SBP norm
# (already weighted by timeIntegrator%norm * timeStepSize per
# bin/control_space_norm.f90:237). Output: OneDWave.norm_<actuator>.dat.
mpirun -n 1 ./control_space_norm

# Now enable the controller. The first TAO callback is the y=0 baseline:
# it writes a zero controlForcing.dat, runs ./forward (populating snapshots
# and J.txt), then ./adjoint (which consumes the just-written snapshots).
sed -i 's/^controller_switch = .*/controller_switch = true/' magudi.inp

# Run python3 directly (NOT under mpirun). If we wrap it in `mpirun -n 1`,
# the MPICH PMI env vars (PMI_FD, PMI_RANK, ...) leak into the subprocess
# environment of ./forward and ./adjoint, and their MPI_Init fails trying to
# register with the parent PMI manager. Without mpirun, both this process and
# its children use MPICH's singleton-init path independently. Phase 1 pins
# np=1; multi-rank will need an env-stripping or MPI_Comm_spawn approach.

BASELINE=2
SPLIT=1

# Baseline: one BASELINE-iter run. Kept small to keep CI fast; the restart
# wiring is shallow at SPLIT=1 (one TAO step per leg), so this catches an
# I/O break but not bugs that only surface after several L-BFGS updates.
python3 "${REPO_ROOT}/utils/optimization_ver4/optim.py" optim.single.yml --max-iter ${BASELINE}
J_baseline=$(tr -d '[:space:]' < OneDWave.forward_run.txt)

rm -f *.petsc

# Split: SPLIT iters, checkpoint, resume for SPLIT more. Final J should match
# the baseline to near-bitwise precision (toy_optim.py --verify-restart proved
# rel diff = 0 on the algorithm; this checks the magudi I/O wiring).
python3 "${REPO_ROOT}/utils/optimization_ver4/optim.py" optim.single.yml --max-iter ${SPLIT}
python3 "${REPO_ROOT}/utils/optimization_ver4/optim.py" optim.single.yml --max-iter ${SPLIT}
J_split=$(tr -d '[:space:]' < OneDWave.forward_run.txt)

python3 - "${J_baseline}" "${J_split}" "${BASELINE}" "${SPLIT}" <<'PY'
import sys
Jb, Js = float(sys.argv[1]), float(sys.argv[2])
baseline_iters, split_iters = sys.argv[3], sys.argv[4]
rel = abs(Jb - Js) / max(abs(Jb), 1e-30)
print(f"baseline J ({baseline_iters}-iter)                = {Jb:.12e}")
print(f"split    J ({split_iters} + resumed {split_iters}) = {Js:.12e}")
print(f"relative difference         = {rel:.3e}")
TOL = 1e-12
if rel > TOL:
    print(f"FAIL: rel diff {rel:.3e} exceeds tolerance {TOL:.0e}",
          file=sys.stderr)
    sys.exit(1)
print(f"OK: restart equivalence within {TOL:.0e}")
PY
