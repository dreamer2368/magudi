#!/usr/bin/env bash
# msforward/msadjoint gradient accuracy test on OneDWave.
#
# Stage build/OneDWave_msgrad/, generate per-segment IC files via a warmup
# forward, then run the Taylor test driver (msgrad_test.py).
#
# Pinned to np=1 for the same MPI_File_write_all reasons as run_strawman.sh:
# actuator .dat files are per-rank slabs concatenated, so np=1 gives the
# single contiguous real64 buffer the Python perturbation code assumes.
#
# Run from the build/ directory after `make -j` has produced bin/forward,
# bin/msforward, bin/msadjoint, bin/zaxpy, bin/qfile_zaxpy, bin/control_space_norm.
#
# Tunables (env vars):
#   NSPLIT          (default 4)        number of msforward/msadjoint segments
#   NTS             (default 120)      timesteps per segment
#   PENALTY_WEIGHT  (default 1.0e-2)   matching penalty weight (0 turns it off)
#   MODE            (default full)     full | ctrl | ic
#                                      ctrl = perturb only control_forcing,
#                                             reference = sum control_forcing IPs.
#                                      ic   = perturb only ICs (k>=1),
#                                             reference = sum ic IPs.
#                                      full = perturb both, reference = total <g,g>.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$(pwd)}"
WORK_DIR="${BUILD_DIR}/OneDWave_msgrad"

NSPLIT="${NSPLIT:-4}"
NTS="${NTS:-120}"
PENALTY_WEIGHT="${PENALTY_WEIGHT:-1.0e-2}"
MODE="${MODE:-full}"
TOTAL_TS=$(( NSPLIT * NTS ))

for b in forward msforward msadjoint zaxpy qfile_zaxpy control_space_norm; do
    if [ ! -x "${BUILD_DIR}/bin/${b}" ]; then
        echo "error: ${BUILD_DIR}/bin/${b} not found; run cmake + make first" >&2
        exit 1
    fi
done

rm -rf "${WORK_DIR}"
mkdir -p "${WORK_DIR}"
cp -r "${REPO_ROOT}/examples/OneDWave/." "${WORK_DIR}/"
cp "${REPO_ROOT}/utils/python/plot3dnasa.py" "${WORK_DIR}/"

cd "${WORK_DIR}"

python3 config.py

for b in forward msforward msadjoint zaxpy qfile_zaxpy control_space_norm; do
    ln -sf "${BUILD_DIR}/bin/${b}" .
done

# Mutate magudi.inp for ms test: enable controller, save snapshots every NTS,
# set number_of_timesteps = TOTAL_TS for the warmup forward (which is one big
# pass), append the time_splitting block. We will narrow number_of_timesteps
# back to NTS after the warmup, since msforward calls runForward per segment.
sed -i "s/^number_of_timesteps = .*/number_of_timesteps = ${TOTAL_TS}/" magudi.inp
sed -i "s/^save_interval = .*/save_interval = ${NTS}/" magudi.inp
sed -i 's/^controller_switch = .*/controller_switch = true/' magudi.inp
cat >> magudi.inp <<EOF

# msadjoint requires this key to exist (it mutates the value per segment).
adjoint_nonzero_initial_condition = false

# Time splitting (multi-segment) parameters for msforward/msadjoint.
time_splitting/number_of_segments = ${NSPLIT}
time_splitting/segment_length = ${NTS}
time_splitting/start_timestep = 0
time_splitting/penalty_weight = ${PENALTY_WEIGHT}
EOF

# control_space_norm writes OneDWave.norm_controlRegion.dat for the full
# trajectory (it reads number_of_timesteps = TOTAL_TS from magudi.inp).
mpirun -n 1 ./control_space_norm

# Zero controlForcing.dat at the matching size (so the warmup forward applies
# no control and produces the natural trajectory).
python3 - <<'PY'
import numpy as np
M = np.fromfile("OneDWave.norm_controlRegion.dat", dtype="<f8")
np.zeros_like(M).tofile("OneDWave.control_forcing_controlRegion.dat")
PY

# Warmup: one forward over TOTAL_TS timesteps with zero control. save_interval=NTS
# produces OneDWave-{NTS:08d}.q, ..., OneDWave-{TOTAL_TS:08d}.q.
mpirun -n 1 ./forward --input magudi.inp

# Build per-segment IC files. ic_0 = original IC (timestep 0); ic_k for k>=1
# comes from the k-th segment-boundary snapshot of the warmup trajectory.
cp OneDWave.ic.q OneDWave-0.ic.q
for k in $(seq 1 $((NSPLIT - 1))); do
    ts=$(printf "%08d" $((k * NTS)))
    cp "OneDWave-${ts}.q" "OneDWave-${k}.ic.q"
done

# Narrow number_of_timesteps to NTS so each msforward segment's runForward
# advances by Nts timesteps. The controlForcing.dat / norm file sizes still
# cover the full TOTAL_TS trajectory; msforward's controlTimestepOffset = k*Nts
# seeks the actuator into the right window for segment k.
sed -i "s/^number_of_timesteps = .*/number_of_timesteps = ${NTS}/" magudi.inp

# Drop the test config the Python driver reads.
cat > msgrad.yml <<EOF
output_prefix: OneDWave
patches: [controlRegion]
time_splitting:
  number_of_segments: ${NSPLIT}
  segment_length: ${NTS}
  penalty_weight: ${PENALTY_WEIGHT}
finite_difference:
  step_sizes: [1.0e-1, 3.0e-2, 1.0e-2, 3.0e-3, 1.0e-3, 3.0e-4, 1.0e-4, 3.0e-5, 1.0e-5, 3.0e-6, 1.0e-6, 3.0e-7, 1.0e-7]
  order_threshold: 0.5
EOF

# Direct python3 (no mpirun); see run_strawman.sh for PMI env-leak rationale.
python3 "${REPO_ROOT}/utils/optimization_ver4/msgrad_test.py" msgrad.yml --mode "${MODE}"