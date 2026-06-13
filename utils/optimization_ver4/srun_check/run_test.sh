#!/usr/bin/env bash
# Standalone nested-srun smoke test.
#
# Must be invoked from inside this directory (./worker is path-relative).
#
# Interactive use:
#   cd utils/optimization_ver4/srun_check
#   make
#   salloc -N 1 -n 8 bash run_test.sh
#
# Batch use: wrap with an sbatch script that cd's here, then calls run_test.sh.
#
# Tunables (env):
#   N_PETSC     number of parent ranks (default 1)
#   N_WORKERS   number of worker ranks (default 4)
#   TIMEOUT     seconds before the run is considered hung (default 60)
#   PYLAUNCHER  command to launch the parent (default "srun"; set to
#               "mpirun" for non-SLURM environments, e.g. dev container)
#   LAUNCHER    command the parent uses for the worker (default
#               "srun --overlap"; set to "mpirun" for non-SLURM)
set -uo pipefail

N_PETSC=${N_PETSC:-1}
N_WORKERS=${N_WORKERS:-4}
TIMEOUT=${TIMEOUT:-60}
PYLAUNCHER=${PYLAUNCHER:-srun}
export LAUNCHER=${LAUNCHER:-srun --overlap}

mkdir -p out
rm -f out/worker.out

# shellcheck disable=SC2086
timeout "${TIMEOUT}" ${PYLAUNCHER} -n "${N_PETSC}" python3 test_subprocess.py "${N_WORKERS}"
rc=$?

if [ "${rc}" -eq 124 ]; then
    echo "FAIL: timed out after ${TIMEOUT}s (nested srun hung)" >&2
    [ -f out/worker.out ] && { echo "--- captured worker output ---" >&2; cat out/worker.out >&2; }
    exit 1
fi
if [ "${rc}" -ne 0 ]; then
    echo "FAIL: test_subprocess.py exited with code ${rc}" >&2
    [ -f out/worker.out ] && { echo "--- captured worker output ---" >&2; cat out/worker.out >&2; }
    exit 1
fi

log="out/worker.out"
if [ ! -f "${log}" ]; then
    echo "FAIL: ${log} not found" >&2
    exit 1
fi
got=$(grep -c '^hello from rank ' "${log}")
if [ "${got}" -ne "${N_WORKERS}" ]; then
    echo "FAIL: expected ${N_WORKERS} worker ranks in ${log}, saw ${got}" >&2
    cat "${log}" >&2
    exit 1
fi
echo "PASS: nested-srun launch completed for N_PETSC=${N_PETSC}, N_WORKERS=${N_WORKERS}"
