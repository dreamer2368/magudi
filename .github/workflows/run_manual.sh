#!/usr/bin/env bash
# ManualBQNLS driver smoke test: stage OneDWave, run the TAO baseline for a
# reference J, then run the manual driver via a bash resubmit loop (one
# msforward and at most one msadjoint per invocation) and compare the
# converged J.
#
# The manual driver targets PLAN_case_a.md "Case A" (one fwd+adj per SLURM
# allocation). In CI we don't have SLURM, so the bash loop simulates the
# allocation boundary by exiting Python and restarting it between every
# state-machine step.
#
# Run from the build/ directory after `make -j` has produced
#   bin/forward bin/msforward bin/msadjoint bin/compute_norm

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$(pwd)}"
WORK_DIR="${BUILD_DIR}/OneDWave_manual"
CFG="optim.yml"

# Launch mode: base = mpirun (default, CI / dev). slurm = srun (production).
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

# Same magudi.inp preparation as run_parallel.sh: full-trajectory step count
# for the baseline forward + compute_norm, controller off.
sed -i "s/^number_of_timesteps = .*/number_of_timesteps = $((NSPLIT * NTS))/" magudi.inp
sed -i 's/^controller_switch = .*/controller_switch = false/' magudi.inp

# 1. Baseline forward (controller off) -> snapshots at every segment boundary.
${LAUNCHER} -n "${N_FORWARD}" ./forward

# 2. Stage per-segment IC files from the baseline.
cp "${PREFIX}.ic.q" "${PREFIX}-0.ic.q"
for k in $(seq 1 $((NSPLIT - 1))); do
    ts=$(printf "%08d" $((START_TS + k * NTS)))
    cp "${PREFIX}-${ts}.q" "${PREFIX}-${k}.ic.q"
done

# 3. Apply forced_inputs and run compute_norm (writes .norm_* and .layout.txt).
python3 - <<PY
from magudi_utils.inputs import InputParser
from magudi_utils.parallel_io import apply_magudi_inp
apply_magudi_inp(InputParser("${CFG}"))
PY
sed -i "s/^number_of_timesteps = .*/number_of_timesteps = $((NSPLIT * NTS))/" magudi.inp

${LAUNCHER} -n "${N_FORWARD}" ./compute_norm --input magudi.inp

# 4. Snapshot the pristine baseline ICs so both runs start from the same y_0.
for k in $(seq 0 $((NSPLIT - 1))); do
    cp "${PREFIX}-${k}.ic.q" "${PREFIX}-${k}.ic.q.bak"
done

restore_baseline_inputs() {
    rm -f *.petsc "${PREFIX}.control_forcing_"*.dat lbfgs_state.json
    for k in $(seq 0 $((NSPLIT - 1))); do
        cp "${PREFIX}-${k}.ic.q.bak" "${PREFIX}-${k}.ic.q"
    done
}

# K outer iterations for both drivers; chosen small enough to fit comfortably
# in CI but large enough that L-BFGS history actually starts to matter.
K=2

# Sanity cap on the manual driver's bash loop. With golden-section LS the
# bracket+linmin pass takes ~6-10 trial forwards plus one accept allocation
# per outer iter, so K=2 fits well under 30 allocations on OneDWave.
MAX_ALLOC_LOOPS=60

# ===========================================================================
# 5. Baseline: TAO driver for K iterations in one process. Captures J_tao.
# ===========================================================================
${LAUNCHER} -n "${N_PETSC}" magudi-optim \
    "${CFG}" --max-iter "${K}" --mode "${MODE}"
J_tao=$(tr -d '[:space:]' < "${PREFIX}.forward_run.txt")

restore_baseline_inputs

# ===========================================================================
# 6. Manual driver: bash loop with one msforward (or one msforward+msadjoint)
#    per invocation. Resubmits while exit == 42 (EXIT_RESUME); stops on
#    exit == 0 (EXIT_CONVERGED). Patches optim.yml to select the manual path.
# ===========================================================================
python3 - <<PY
import yaml
with open("${CFG}") as f: c = yaml.safe_load(f)
opt = c.setdefault("optimization", {})
opt["type"] = "manual"
opt.setdefault("line_search", {})["type"] = "golden_section"
with open("${CFG}", "w") as f: yaml.safe_dump(c, f)
PY

n_alloc=0
converged=0
for i in $(seq 1 $MAX_ALLOC_LOOPS); do
    n_alloc=$i
    set +e
    ${LAUNCHER} -n "${N_PETSC}" magudi-optim \
        "${CFG}" --max-iter "${K}" --mode "${MODE}"
    status=$?
    set -e
    case "$status" in
        0)
            echo "manual driver converged after ${n_alloc} allocations"
            converged=1
            break
            ;;
        42)
            continue
            ;;
        *)
            echo "FATAL: manual driver exited with code $status" >&2
            exit 1
            ;;
    esac
done

if [ $converged -ne 1 ]; then
    echo "FAIL: manual driver did not converge in ${MAX_ALLOC_LOOPS} allocations" >&2
    exit 1
fi

J_manual=$(tr -d '[:space:]' < "${PREFIX}.forward_run.txt")

# ===========================================================================
# 7. Compare. The TAO (BQNLS+More-Thuente) and manual (L-BFGS+golden-section)
#    drivers trace different paths through the same J landscape, so an exact
#    match is not expected. Test gates:
#       (a) manual J is finite (no crash sentinel),
#       (b) within a wide relative band of the TAO J (sanity, not equivalence).
# ===========================================================================
python3 - "${J_tao}" "${J_manual}" "${K}" "${n_alloc}" <<'PY'
import math
import sys
J_tao, J_manual = float(sys.argv[1]), float(sys.argv[2])
K, n_alloc = sys.argv[3], sys.argv[4]
print(f"baseline TAO    J ({K}-iter, 1 process)            = {J_tao:.12e}")
print(f"manual          J ({K}-iter, {n_alloc} allocations)        = {J_manual:.12e}")
if not math.isfinite(J_manual):
    print("FAIL: manual J is not finite", file=sys.stderr)
    sys.exit(1)
rel = abs(J_tao - J_manual) / max(abs(J_tao), 1e-30)
print(f"relative difference                                 = {rel:.3e}")
# Wide band: golden-section + manual L-BFGS converges along a different path
# from BQNLS+MT; the gate is "in the same neighborhood", not "equivalent".
TOL = 1.0e-2
if rel > TOL:
    print(f"FAIL: rel diff {rel:.3e} exceeds tolerance {TOL:.0e}",
          file=sys.stderr)
    sys.exit(1)
print(f"OK: manual J within {TOL:.0e} of TAO J after {K} outer iterations")
PY
