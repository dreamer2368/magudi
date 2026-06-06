#!/usr/bin/env bash
# Frozen Phase 1/2 strawman wrapper — kept for regression vs the Phase 3 rewrite.
#
# optim_ver4 MPI-parallel driver: stage OneDWave, run TAO L-BFGS with
# N_petsc Python ranks spawning N_forward / N_adjoint child workers.
#
# All three rank counts come from optim.parallel.yml's resource_distribution.
# jobs.{forward,adjoint,petsc} block. control_space_norm runs at N_forward
# so the .dat layout matches what forward/adjoint will read/write.
#
# Child binaries are launched via MPI_Comm_spawn from inside Python -- no
# env-stripping, no PMI leak, and no modification to forward/adjoint. If the
# pre-flight Spawn check below fails, the host's MPI runtime does not support
# dynamic process management (most likely SLURM srun) and the only paths
# forward are env-stripping or running mpirun inside the SLURM allocation.
#
# Run from the build/ directory after `make -j` has produced bin/forward,
# bin/adjoint, bin/control_space_norm.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$(pwd)}"
WORK_DIR="${BUILD_DIR}/OneDWave_parallel"
CFG="optim.parallel.yml"

if [ ! -x "${BUILD_DIR}/bin/forward" ]; then
    echo "error: ${BUILD_DIR}/bin/forward not found; run cmake + make first" >&2
    exit 1
fi

# # Pre-flight: does this MPI runtime support Comm_spawn with a child that
# # only calls MPI_Init/MPI_Finalize (no MPI_Comm_disconnect on the parent)?
# # That's exactly forward/adjoint's lifecycle. Catch hangs/errors here rather
# # than mid-optimization. Note: child must call MPI_Init, so /bin/echo would
# # hang -- use a python child that calls MPI_Init via mpi4py import.
# echo "Pre-flight: MPI_Comm_spawn smoke test"
# SPAWN_CHILD=$(sudo mktemp /tmp/spawn_child.XXXX.py)
# cat > "${SPAWN_CHILD}" <<'PY'
# from mpi4py import MPI
# # import triggers MPI_Init; module finalizer calls MPI_Finalize.
# PY
# mpirun -n 2 python3 - "${SPAWN_CHILD}" <<'PY'
# import sys
# from mpi4py import MPI
# inter = MPI.COMM_WORLD.Spawn("python3", args=[sys.argv[1]], maxprocs=1)
# inter.Disconnect()
# PY
# rm -f "${SPAWN_CHILD}"

rm -rf "${WORK_DIR}"
mkdir -p "${WORK_DIR}"
cp -r "${REPO_ROOT}/examples/OneDWave/." "${WORK_DIR}/"
cp "${REPO_ROOT}/utils/python/plot3dnasa.py" "${WORK_DIR}/"

cd "${WORK_DIR}"

python3 config.py

ln -sf "${BUILD_DIR}/bin/forward" .
ln -sf "${BUILD_DIR}/bin/adjoint" .
ln -sf "${BUILD_DIR}/bin/control_space_norm" .

# Parse N_forward and N_petsc from the YAML once -- bash invocations below
# need them as integers.
N_FORWARD=$(python3 -c "
import yaml
with open('${CFG}') as f: c = yaml.safe_load(f)
print(c['resource_distribution']['jobs']['forward'][1])
")
N_PETSC=$(python3 -c "
import yaml
with open('${CFG}') as f: c = yaml.safe_load(f)
print(c['resource_distribution']['jobs']['petsc'][1])
")
echo "N_forward = ${N_FORWARD}, N_petsc = ${N_PETSC}"

# Patch magudi.inp for the Nsplit=1 layout, matching run_strawman.sh:
sed -i 's/^save_interval = .*/save_interval = 480/' magudi.inp
sed -i 's/^baseline_prediction_available = .*/baseline_prediction_available = true/' magudi.inp

# control_space_norm at N_forward so the .dat ordering matches forward/adjoint.
mpirun -n "${N_FORWARD}" ./control_space_norm

sed -i 's/^controller_switch = .*/controller_switch = true/' magudi.inp

BASELINE=2
SPLIT=1

# Baseline: BASELINE iters in one go.
mpirun -n "${N_PETSC}" python3 \
    "${REPO_ROOT}/utils/optimization_ver4/optim.parallel.strawman.py" \
    "${CFG}" --max-iter "${BASELINE}"
J_baseline=$(tr -d '[:space:]' < OneDWave.forward_run.txt)

rm -f *.petsc

# Split: SPLIT iters, checkpoint, resume for SPLIT more. Final J should match
# the baseline to near-bitwise precision (toy_optim.py --verify-restart proved
# rel diff = 0 on the algorithm; this checks the parallel-I/O wiring).
mpirun -n "${N_PETSC}" python3 \
    "${REPO_ROOT}/utils/optimization_ver4/optim.parallel.strawman.py" \
    "${CFG}" --max-iter "${SPLIT}"
mpirun -n "${N_PETSC}" python3 \
    "${REPO_ROOT}/utils/optimization_ver4/optim.parallel.strawman.py" \
    "${CFG}" --max-iter "${SPLIT}"
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
