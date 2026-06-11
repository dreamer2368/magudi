#!/usr/bin/env bash
# msforward/msadjoint Phase 3 gradient accuracy test on OneDWave.
#
# Stage build/OneDWave_msgrad/, run a warmup forward to generate per-segment
# IC files, run compute_norm (parallel) to produce M_diag / layout.txt, then
# launch msgrad_test.py under mpirun -n N_PETSC. msgrad_test.py spawns
# msforward / msadjoint / zaxpy / qfile_zaxpy via MPI_Comm_spawn using the
# resource_distribution block we inject into msgrad.yml.
#
# Run from the build/ directory after `make -j` has produced bin/forward,
# bin/msforward, bin/msadjoint, bin/zaxpy, bin/qfile_zaxpy, bin/compute_norm.
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
#   ADD_IC_NOISE    (default unset)    if non-empty, pass --add-ic-noise to
#                                      msgrad_test.py: pointwise uniform
#                                      [-1e-4, 1e-4] noise on intermediate
#                                      ic slabs before the baseline pass.
#   N_FORWARD       (default 2)        msforward child rank count
#   N_ADJOINT       (default 2)        msadjoint child rank count
#   N_PETSC         (default NSPLIT)   Python driver rank count; must satisfy
#                                      the ParallelIOHandler policy:
#                                        n_actuator * nprocs_per_actuator + n_ic.
#                                      For OneDWave (1 actuator), N_PETSC = NSPLIT
#                                      gives nprocs_per_actuator = 1.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$(pwd)}"
WORK_DIR="${BUILD_DIR}/OneDWave_msgrad"

NSPLIT="${NSPLIT:-4}"
NTS="${NTS:-120}"
PENALTY_WEIGHT="${PENALTY_WEIGHT:-1.0e-2}"
MODE="${MODE:-full}"
ADD_IC_NOISE="${ADD_IC_NOISE:-1}"
TOTAL_TS=$(( NSPLIT * NTS ))

NOISE_FLAG=""
if [ -n "${ADD_IC_NOISE}" ]; then
    NOISE_FLAG="--add-ic-noise"
fi

N_FORWARD="${N_FORWARD:-7}"
N_ADJOINT="${N_ADJOINT:-7}"
N_PETSC="${N_PETSC:-1}"

for b in forward msforward msadjoint zaxpy qfile_zaxpy compute_norm; do
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

for b in forward msforward msadjoint zaxpy qfile_zaxpy compute_norm; do
    ln -sf "${BUILD_DIR}/bin/${b}" .
done

# Mutate magudi.inp for the ms test: enable controller, save snapshots every
# NTS, set number_of_timesteps = TOTAL_TS for the warmup forward and for
# compute_norm (which gracefulExits unless solver%nTimesteps == Nts*Nsplit).
# Append the time_splitting block. number_of_timesteps narrows to NTS only
# right before launching msgrad_test.py.
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

# Warmup: one forward over TOTAL_TS timesteps with zero control. We have not
# yet seeded OneDWave.control_forcing_controlRegion.dat, but msgrad_test.py
# writes a zero one of the right size later anyway. For the warmup the
# controller path needs a file to load -- create a zero one sized from the
# *intended* layout (nPatchPoints * nUnknowns * TOTAL_TS * nStages * 8); the
# simplest path: run compute_norm FIRST so we know the size, then truncate the
# norm file's bytes to zero in a copy.
mpirun -n ${N_FORWARD} ./compute_norm --input magudi.inp
python3 - <<'PY'
import os
size = os.path.getsize("OneDWave.norm_controlRegion.dat")
with open("OneDWave.control_forcing_controlRegion.dat", "wb") as fh:
    fh.write(b"\x00" * size)
PY

# save_interval = NTS produces OneDWave-{NTS:08d}.q, ..., OneDWave-{TOTAL_TS:08d}.q.
mpirun -n ${N_FORWARD} ./forward --input magudi.inp

# Build per-segment IC files. ic_0 = original IC (timestep 0); ic_k for k>=1
# comes from the k-th segment-boundary snapshot of the warmup trajectory.
cp OneDWave.ic.q OneDWave-0.ic.q
for k in $(seq 1 $((NSPLIT - 1))); do
    ts=$(printf "%08d" $((k * NTS)))
    cp "OneDWave-${ts}.q" "OneDWave-${k}.ic.q"
done

# Narrow number_of_timesteps to NTS so each msforward segment's runForward
# advances by Nts timesteps. The .control_forcing.dat and .norm files still
# cover the full TOTAL_TS trajectory; msforward's controlTimestepOffset = k*Nts
# seeks the actuator into the right window for segment k.
sed -i "s/^number_of_timesteps = .*/number_of_timesteps = ${NTS}/" magudi.inp

# Drop the test config the Python driver reads. resource_distribution mirrors
# optim.parallel.yml's convention; msgrad_test.py consumes the procs (index 1).
cat > msgrad.yml <<EOF
output_prefix: OneDWave
patches: [controlRegion]
time_splitting:
  number_of_segments: ${NSPLIT}
  segment_length: ${NTS}
  penalty_weight: ${PENALTY_WEIGHT}
  state_controllability: 1.0
finite_difference:
  step_sizes: [1.0e-1, 3.0e-2, 1.0e-2, 3.0e-3, 1.0e-3, 3.0e-4, 1.0e-4, 3.0e-5, 1.0e-5, 3.0e-6, 1.0e-6, 3.0e-7, 1.0e-7]
  order_threshold: 0.5
resource_distribution:
  jobs:
    forward: [1, ${N_FORWARD}]
    adjoint: [1, ${N_ADJOINT}]
    petsc:   [1, ${N_PETSC}]
EOF

# msgrad_test.py runs under mpirun -n N_PETSC. All child binaries
# (msforward, msadjoint, zaxpy, qfile_zaxpy) go through MPI_Comm_spawn so PMI
# env vars don't leak. zaxpy / qfile_zaxpy are spawned at n=1 (.dat / .q are
# rank-portable; perturbation is O(N) memory work).
mpirun -n "${N_PETSC}" python3 \
    "${REPO_ROOT}/utils/optimization_ver4/msgrad_test.py" \
    msgrad.yml --mode "${MODE}" ${NOISE_FLAG}