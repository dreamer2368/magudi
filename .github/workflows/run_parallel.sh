#!/usr/bin/env bash
# optim_ver4 driver: stage OneDWave, run baseline forward to seed the
# per-segment ICs, then drive TAO L-BFGS via msforward / msadjoint with
# N_petsc Python ranks spawning N_forward / N_adjoint child workers.
#
# All three rank counts come from optim.yml's resource_distribution.
# jobs.{forward,adjoint,petsc} block. compute_norm runs at N_forward so the
# .norm_<actuator>.dat layout matches what msforward/msadjoint will produce.
#
# Run from the build/ directory after `make -j` has produced
#   bin/forward bin/msforward bin/msadjoint bin/compute_norm

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$(pwd)}"
WORK_DIR="${BUILD_DIR}/OneDWave_msparallel"
CFG="optim.yml"

# Launch mode. base = mpirun (default, dev container / interactive). slurm =
# srun (production allocation). Used by both the standalone forward / compute_norm
# steps and by magudi-optim (which is passed through as --mode).
MODE="${MODE:-base}"
case "${MODE}" in
    base)  LAUNCHER="mpirun" ;;
    slurm) LAUNCHER="srun"   ;;
    *) echo "error: MODE must be 'base' or 'slurm', got '${MODE}'" >&2; exit 1 ;;
esac

for bin in forward msforward msadjoint compute_norm; do
    if [ ! -x "${BUILD_DIR}/bin/${bin}" ]; then
        echo "error: ${BUILD_DIR}/bin/${bin} not found; run cmake + make first" >&2
        exit 1
    fi
done

rm -rf "${WORK_DIR}"
mkdir -p "${WORK_DIR}"
cp -r "${REPO_ROOT}/examples/OneDWave/." "${WORK_DIR}/"

cd "${WORK_DIR}"

python3 config.py

ln -sf "${BUILD_DIR}/bin/forward" .
ln -sf "${BUILD_DIR}/bin/msforward" .
ln -sf "${BUILD_DIR}/bin/msadjoint" .
ln -sf "${BUILD_DIR}/bin/compute_norm" .

# Parse rank counts and the time-splitting layout from optim.yml once.
read N_FORWARD N_PETSC NSPLIT NTS START_TS PENALTY_WEIGHT STATE_CTRL PREFIX <<EOF
$(python3 - <<PY
import yaml
with open("${CFG}") as f: c = yaml.safe_load(f)
ts = c["magudi"]["forced_inputs"]["time_splitting"]
print(
    c["resource_distribution"]["jobs"]["forward"][1],
    c["resource_distribution"]["jobs"]["petsc"][1],
    ts["number_of_segments"],
    ts["segment_length"],
    ts["start_timestep"],
    ts["penalty_weight"],
    ts["state_controllability"],
    c["global_prefix"],
)
PY
)
EOF
echo "N_forward=${N_FORWARD} N_petsc=${N_PETSC} Nsplit=${NSPLIT} Nts=${NTS} start_ts=${START_TS} prefix=${PREFIX}"

# Magudi.inp mutations needed BEFORE the optim driver runs (which only fires
# set_magudi_inp() at startup). The baseline forward + compute_norm need:
#   - number_of_timesteps = NSPLIT*NTS so the baseline spans every segment
#     (solver%setup reads it as the per-run step count; src/SolverImpl.f90:988
#     uses region%timestep + nTimesteps to locate the terminal snapshot).
#   - controller_switch = false for the unforced baseline.
# Everything else (save_interval, baseline_prediction_available, controller_switch
# for the optim run, time_splitting/*, adjoint_nonzero_initial_condition) is
# declared in optim.yml's magudi.forced_inputs and applied by
# ParallelIOHandler.set_magudi_inp() at the start of magudi-optim.
sed -i "s/^number_of_timesteps = .*/number_of_timesteps = $((NSPLIT * NTS))/" magudi.inp
sed -i 's/^controller_switch = .*/controller_switch = false/' magudi.inp

# 1. Baseline forward run with controller_switch=false: produces snapshots
#    at every segment boundary (ts = NTS, 2*NTS, ..., NSPLIT*NTS).
${LAUNCHER} -n "${N_FORWARD}" ./forward

# 2. Stage per-segment IC files. Segment 0 IC is the canonical initial state
#    (read from the .ic.q produced by config.py). Segments 1..NSPLIT-1 take the
#    baseline snapshot at the corresponding timestep.
cp "${PREFIX}.ic.q" "${PREFIX}-0.ic.q"
for k in $(seq 1 $((NSPLIT - 1))); do
    ts=$(printf "%08d" $((START_TS + k * NTS)))
    cp "${PREFIX}-${ts}.q" "${PREFIX}-${k}.ic.q"
done

# 3. Apply optim.yml's magudi.forced_inputs to magudi.inp once, so compute_norm
#    sees time_splitting/* (it reads them via getRequiredOption and would
#    gracefulExit otherwise) and the actuator patch is enabled
#    (controller_switch=true is in forced_inputs). Then bump number_of_timesteps
#    back to NSPLIT*NTS just for compute_norm -- the YAML value (= NTS,
#    per-segment) is wrong here because compute_norm needs the full-trajectory
#    step count to allocate the norm buffer. The optim driver re-applies
#    set_magudi_inp at startup, restoring number_of_timesteps to NTS.
python3 - <<PY
from magudi_utils.inputs import InputParser
from magudi_utils.parallel_io import apply_magudi_inp
apply_magudi_inp(InputParser("${CFG}"))
PY
sed -i "s/^number_of_timesteps = .*/number_of_timesteps = $((NSPLIT * NTS))/" magudi.inp

# 4. compute_norm writes .norm_<actuator>.dat (5D-subarray, NSPLIT*NTS frames),
#    .norm_ic.q (PLOT3D solution, shared across all ic slabs of M_diag), and
#    .layout.txt for Python.
${LAUNCHER} -n "${N_FORWARD}" ./compute_norm --input magudi.inp

BASELINE=2
SPLIT=1

# Snapshot the pristine baseline ICs (and the .control_forcing_*.dat seed if
# Python created one) before any optim run mutates them. Both the baseline-run
# and the split-run start by reading these to initialize y, so we must reset
# them between runs to compare apples to apples.
for k in $(seq 0 $((NSPLIT - 1))); do
    cp "${PREFIX}-${k}.ic.q" "${PREFIX}-${k}.ic.q.bak"
done

restore_baseline_inputs() {
    rm -f *.petsc "${PREFIX}.control_forcing_"*.dat
    for k in $(seq 0 $((NSPLIT - 1))); do
        cp "${PREFIX}-${k}.ic.q.bak" "${PREFIX}-${k}.ic.q"
    done
}

# Baseline: BASELINE iters in one go.
${LAUNCHER} -n "${N_PETSC}" magudi-optim \
    "${CFG}" --max-iter "${BASELINE}" --mode "${MODE}"
J_baseline=$(tr -d '[:space:]' < "${PREFIX}.forward_run.txt")

restore_baseline_inputs

# Split: SPLIT iters, checkpoint, resume for SPLIT more. Final J should match
# the baseline to near-bitwise precision; this checks the parallel-I/O wiring
# under Phase 3's file-aligned layout.
${LAUNCHER} -n "${N_PETSC}" magudi-optim \
    "${CFG}" --max-iter "${SPLIT}" --mode "${MODE}"
${LAUNCHER} -n "${N_PETSC}" magudi-optim \
    "${CFG}" --max-iter "${SPLIT}" --mode "${MODE}"
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
