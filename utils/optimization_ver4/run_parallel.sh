#!/usr/bin/env bash
# optim_ver4 Phase 3 driver: stage OneDWave, run baseline forward to seed the
# per-segment ICs, then drive TAO L-BFGS via msforward / msadjoint with
# N_petsc Python ranks spawning N_forward / N_adjoint child workers.
#
# All three rank counts come from optim.parallel.yml's resource_distribution.
# jobs.{forward,adjoint,petsc} block. compute_norm runs at N_forward so the
# .norm_<actuator>.dat layout matches what msforward/msadjoint will produce.
#
# Run from the build/ directory after `make -j` has produced
#   bin/forward bin/msforward bin/msadjoint bin/compute_norm

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$(pwd)}"
WORK_DIR="${BUILD_DIR}/OneDWave_msparallel"
CFG="optim.parallel.yml"

for bin in forward msforward msadjoint compute_norm; do
    if [ ! -x "${BUILD_DIR}/bin/${bin}" ]; then
        echo "error: ${BUILD_DIR}/bin/${bin} not found; run cmake + make first" >&2
        exit 1
    fi
done

rm -rf "${WORK_DIR}"
mkdir -p "${WORK_DIR}"
cp -r "${REPO_ROOT}/examples/OneDWave/." "${WORK_DIR}/"
cp "${REPO_ROOT}/utils/python/plot3dnasa.py" "${WORK_DIR}/"

cd "${WORK_DIR}"

python3 config.py

ln -sf "${BUILD_DIR}/bin/forward" .
ln -sf "${BUILD_DIR}/bin/msforward" .
ln -sf "${BUILD_DIR}/bin/msadjoint" .
ln -sf "${BUILD_DIR}/bin/compute_norm" .

# Parse rank counts and the time-splitting layout from optim.parallel.yml once.
read N_FORWARD N_PETSC NSPLIT NTS START_TS PENALTY_WEIGHT STATE_CTRL PREFIX <<EOF
$(python3 - <<PY
import yaml
with open("${CFG}") as f: c = yaml.safe_load(f)
print(
    c["resource_distribution"]["jobs"]["forward"][1],
    c["resource_distribution"]["jobs"]["petsc"][1],
    c["time_splitting"]["number_of_segments"],
    c["time_splitting"]["segment_length"],
    c["time_splitting"]["start_timestep"],
    c["time_splitting"]["penalty_weight"],
    c["time_splitting"]["state_controllability"],
    c["global_prefix"],
)
PY
)
EOF
echo "N_forward=${N_FORWARD} N_petsc=${N_PETSC} Nsplit=${NSPLIT} Nts=${NTS} start_ts=${START_TS} prefix=${PREFIX}"

# Patch magudi.inp for the multi-segment layout. The baseline forward runs with
# controller_switch=false; controller_switch=true is set later for compute_norm
# and the optim run.
sed -i "s/^number_of_timesteps = .*/number_of_timesteps = $((NSPLIT * NTS))/" magudi.inp
sed -i "s/^save_interval = .*/save_interval = ${NTS}/" magudi.inp
sed -i 's/^baseline_prediction_available = .*/baseline_prediction_available = true/' magudi.inp
sed -i 's/^controller_switch = .*/controller_switch = false/' magudi.inp

# msadjoint requires adjoint_nonzero_initial_condition in the dict (it mutates the
# value per segment). Add a placeholder line if it isn't already there.
if ! grep -q '^adjoint_nonzero_initial_condition' magudi.inp; then
    echo 'adjoint_nonzero_initial_condition = false' >> magudi.inp
fi

# msforward / msadjoint / compute_norm read time_splitting/* via getRequiredOption.
# magudi.inp uses literal flat keys, so write the slash-style names directly.
# Replace any pre-existing entries with the YAML-driven values.
sed -i '/^time_splitting\//d' magudi.inp
cat >> magudi.inp <<EOF
time_splitting/number_of_segments = ${NSPLIT}
time_splitting/segment_length = ${NTS}
time_splitting/start_timestep = ${START_TS}
time_splitting/penalty_weight = ${PENALTY_WEIGHT}
time_splitting/state_controllability = ${STATE_CTRL}
EOF

# 1. Baseline forward run with controller_switch=false: produces snapshots
#    at every segment boundary (ts = NTS, 2*NTS, ..., NSPLIT*NTS).
mpirun -n "${N_FORWARD}" ./forward

# 2. Stage per-segment IC files. Segment 0 IC is the canonical initial state
#    (read from the .ic.q produced by config.py). Segments 1..NSPLIT-1 take the
#    baseline snapshot at the corresponding timestep.
cp "${PREFIX}.ic.q" "${PREFIX}-0.ic.q"
for k in $(seq 1 $((NSPLIT - 1))); do
    ts=$(printf "%08d" $((START_TS + k * NTS)))
    cp "${PREFIX}-${ts}.q" "${PREFIX}-${k}.ic.q"
done

# 3. Flip controller_switch=true so compute_norm (and the optim run) see the
#    actuator patch with its controlForcing buffer allocated.
sed -i 's/^controller_switch = .*/controller_switch = true/' magudi.inp

# 4. compute_norm writes .norm_<actuator>.dat (5D-subarray), .norm_ic.q (PLOT3D
#    solution, shared across all ic slabs of M_diag), and .layout.txt for Python.
mpirun -n "${N_FORWARD}" ./compute_norm --input magudi.inp

BASELINE=2
SPLIT=1

# Baseline: BASELINE iters in one go.
mpirun -n "${N_PETSC}" python3 \
    "${REPO_ROOT}/utils/optimization_ver4/optim.parallel.py" \
    "${CFG}" --max-iter "${BASELINE}"
J_baseline=$(tr -d '[:space:]' < "${PREFIX}.forward_run.txt")

rm -f *.petsc

# Split: SPLIT iters, checkpoint, resume for SPLIT more. Final J should match
# the baseline to near-bitwise precision; this checks the parallel-I/O wiring
# under Phase 3's file-aligned layout.
mpirun -n "${N_PETSC}" python3 \
    "${REPO_ROOT}/utils/optimization_ver4/optim.parallel.py" \
    "${CFG}" --max-iter "${SPLIT}"
mpirun -n "${N_PETSC}" python3 \
    "${REPO_ROOT}/utils/optimization_ver4/optim.parallel.py" \
    "${CFG}" --max-iter "${SPLIT}"
J_split=$(tr -d '[:space:]' < "${PREFIX}.forward_run.txt")

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
