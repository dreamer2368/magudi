#!/usr/bin/env bash
# msforward/msadjoint gradient accuracy test on OneDWave.
#
# Stage build/OneDWave_msgrad/, run a warmup forward to generate per-segment
# IC files, run compute_norm (parallel) to produce M_diag / layout.txt, then
# launch magudi-msgrad under ${LAUNCHER} -n N_PETSC. magudi-msgrad launches
# msforward / msadjoint via subprocess + ${LAUNCHER} (selected by --exec-mode).
# Control / IC perturbation is done in Python via PETSc.Vec.axpy, so zaxpy /
# qfile_zaxpy binaries are no longer needed.
#
# Run from the build/ directory after `make -j` has produced bin/forward,
# bin/msforward, bin/msadjoint, bin/compute_norm.
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
#                                      magudi-msgrad: pointwise uniform
#                                      [-1e-4, 1e-4] noise on intermediate
#                                      ic slabs before the baseline pass.
#   N_FORWARD       (default 2)        msforward child rank count
#   N_ADJOINT       (default 2)        msadjoint child rank count
#   N_PETSC         (default NSPLIT)   Python driver rank count; must satisfy
#                                      the ParallelIOHandler policy:
#                                        n_actuator * nprocs_per_actuator + n_ic.
#                                      For OneDWave (1 actuator), N_PETSC = NSPLIT
#                                      gives nprocs_per_actuator = 1.
#   EXEC_MODE       (default base)     base | slurm
#                                      base  = mpirun (dev container / interactive)
#                                      slurm = srun (production SLURM allocation);
#                                              passed through to magudi-msgrad
#                                              as --exec-mode.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$(pwd)}"
WORK_DIR="${BUILD_DIR}/OneDWave_msgrad"

NSPLIT="${NSPLIT:-2}"
NTS="${NTS:-240}"
PENALTY_WEIGHT="${PENALTY_WEIGHT:-1.0e-3}"
MODE="${MODE:-full}"
ADD_IC_NOISE="${ADD_IC_NOISE:-1}"
TOTAL_TS=$(( NSPLIT * NTS ))

NOISE_FLAG=""
if [ -n "${ADD_IC_NOISE}" ]; then
    NOISE_FLAG="--add-ic-noise"
fi

N_FORWARD="${N_FORWARD:-1}"
N_ADJOINT="${N_ADJOINT:-1}"
N_PETSC="${N_PETSC:-3}"

EXEC_MODE="${EXEC_MODE:-base}"
case "${EXEC_MODE}" in
    base)  LAUNCHER="mpirun" ;;
    slurm) LAUNCHER="srun"   ;;
    *) echo "error: EXEC_MODE must be 'base' or 'slurm', got '${EXEC_MODE}'" >&2; exit 1 ;;
esac

for b in forward msforward msadjoint compute_norm; do
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

for b in forward msforward msadjoint compute_norm; do
    ln -sf "${BUILD_DIR}/bin/${b}" .
done

# Mutate magudi.inp for the ms test: enable controller, save snapshots every
# NTS, set number_of_timesteps = TOTAL_TS for the warmup forward and for
# compute_norm (which gracefulExits unless solver%nTimesteps == Nts*Nsplit).
# Append the time_splitting block. number_of_timesteps narrows to NTS only
# right before launching magudi-msgrad.
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
time_splitting/penalty_type = huber
time_splitting/state_mollifier/enabled = true
time_splitting/state_mollifier/uniform_in_time = true
EOF

# Warmup: one forward over TOTAL_TS timesteps with zero control. We have not
# yet seeded OneDWave.control_forcing_controlRegion.dat, but magudi-msgrad
# writes a zero one of the right size later anyway. For the warmup the
# controller path needs a file to load -- create a zero one sized from the
# *intended* layout (nPatchPoints * nUnknowns * TOTAL_TS * nStages * 8); the
# simplest path: run compute_norm FIRST so we know the size, then truncate the
# norm file's bytes to zero in a copy.
${LAUNCHER} -n ${N_FORWARD} ./compute_norm --input magudi.inp
python3 - <<'PY'
import os
size = os.path.getsize("OneDWave.norm_controlRegion.dat")
with open("OneDWave.control_forcing_controlRegion.dat", "wb") as fh:
    fh.write(b"\x00" * size)
PY

# save_interval = NTS produces OneDWave-{NTS:08d}.q, ..., OneDWave-{TOTAL_TS:08d}.q.
${LAUNCHER} -n ${N_FORWARD} ./forward --input magudi.inp

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

# Drop the test config the Python driver reads. Schema matches optim.yml so
# ParallelIOHandler(cfg, ...) can derive prefix / Nsplit / log path directly.
# resource_distribution mirrors optim.yml; magudi-msgrad consumes
# the procs (index 1).
cat > msgrad.yml <<EOF
global_prefix: OneDWave
magudi:
  forced_inputs:
    time_splitting:
      number_of_segments: ${NSPLIT}
      segment_length: ${NTS}
      penalty_weight: ${PENALTY_WEIGHT}
      penalty_type: huber
      state_controllability: 1.0
      state_mollifier:
        enabled: true
        uniform_in_time: true
finite_difference:
  step_sizes: [1.0e-1, 3.0e-2, 1.0e-2, 3.0e-3, 1.0e-3, 3.0e-4, 1.0e-4, 3.0e-5, 1.0e-5, 3.0e-6, 1.0e-6, 3.0e-7, 1.0e-7]
  order_threshold: 0.5
resource_distribution:
  jobs:
    forward: [1, ${N_FORWARD}]
    adjoint: [1, ${N_ADJOINT}]
    petsc:   [1, ${N_PETSC}]
EOF

# magudi-msgrad runs under ${LAUNCHER} -n N_PETSC. Child binaries
# (msforward, msadjoint) are launched by magudi-msgrad via subprocess +
# ${LAUNCHER} (mpirun or srun --overlap, selected by --exec-mode).
${LAUNCHER} -n "${N_PETSC}" magudi-msgrad \
    msgrad.yml --mode "${MODE}" --exec-mode "${EXEC_MODE}" ${NOISE_FLAG}