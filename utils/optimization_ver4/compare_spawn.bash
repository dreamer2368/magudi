#!/usr/bin/env bash
# Compare two ways of launching the OneDWave forward run on 8 ranks:
#   1) direct mpirun
#   2) MPI_Comm_spawn via _spawn_and_wait
# Outputs are not redirected.
#
# Run from the build/ directory after `make -j` has produced bin/forward.

set -euo pipefail

NPROCS=7
REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$(pwd)}"
WORK_DIR="${BUILD_DIR}/OneDWave_compare_spawn"

if [ ! -x "${BUILD_DIR}/bin/forward" ]; then
    echo "error: ${BUILD_DIR}/bin/forward not found; run cmake + make first" >&2
    exit 1
fi

rm -rf "${WORK_DIR}"
mkdir -p "${WORK_DIR}"
cp -r "${REPO_ROOT}/examples/OneDWave/." "${WORK_DIR}/"
cp "${REPO_ROOT}/utils/python/plot3dnasa.py" "${WORK_DIR}/"
cp "${REPO_ROOT}/utils/optimization_ver4/spawn_one_forward.py" "${WORK_DIR}/"

cd "${WORK_DIR}"
python3 config.py
ln -sf "${BUILD_DIR}/bin/forward" .

echo "=== Run 1: direct mpirun -n ${NPROCS} ./forward ==="
mpirun -n "${NPROCS}" ./forward --input magudi.inp

echo "=== Run 2: 1 parent rank, spawn ${NPROCS} children via _spawn_and_wait ==="
mpirun -n 1 python3 spawn_one_forward.py "${NPROCS}"
