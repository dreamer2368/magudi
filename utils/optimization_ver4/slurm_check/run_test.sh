#!/usr/bin/env bash
# Standalone PMIx spawn smoke test.
#
# Must be invoked from inside this directory (./spawn_child is path-relative).
#
# Interactive use:
#   cd utils/optimization_ver4/slurm_check
#   make
#   salloc -N 1 -n 8 bash run_test.sh
#
# Batch use: wrap with an sbatch script that cd's here, then calls run_test.sh.
#
# Tunables (env):
#   N        number of children to spawn (default 4)
#   TIMEOUT  seconds before the run is considered hung (default 60)
set -uo pipefail

N=${N:-4}
TIMEOUT=${TIMEOUT:-60}
LOG=$(mktemp)
trap 'rm -f "${LOG}"' EXIT

timeout "${TIMEOUT}" srun --mpi=pmix -n 1 python3 test_spawn.py "${N}" | tee "${LOG}"
rc=${PIPESTATUS[0]}

if [ "${rc}" -eq 124 ]; then
    echo "FAIL: srun timed out after ${TIMEOUT}s (PMIx accepted spawn but Barrier never completed)" >&2
    exit 1
fi
if [ "${rc}" -ne 0 ]; then
    echo "FAIL: test_spawn.py exited with code ${rc}" >&2
    exit 1
fi

# Confirm all N child ranks reported in (one line each).
expected=${N}
got=$(grep -c '^hello from rank ' "${LOG}")
if [ "${got}" -ne "${expected}" ]; then
    echo "FAIL: expected ${expected} child ranks, saw ${got}" >&2
    exit 1
fi
echo "PASS: PMIx spawn handshake completed for N=${N}"
