#!/usr/bin/env bash
# Read-only cluster probe for the nested-srun approach; no allocation needed.
#
# Confirms that:
#   1. srun is on PATH and supports the --overlap flag (SLURM >= 20.11)
#   2. mpi4py imports cleanly
#   3. The cluster MPI is identified (informational)
#
# Exits 0 on PASS, 1 on FAIL.
set -uo pipefail

echo "--- srun --version ---"
srun --version 2>&1 || { echo "FAIL: srun not on PATH"; exit 1; }
echo
echo "--- srun supports --overlap? ---"
if srun --help 2>&1 | grep -q -- '--overlap'; then
    echo "yes"
else
    echo "no"
    echo "FAIL: this SLURM does not recognize --overlap (need SLURM >= 20.11)" >&2
    exit 1
fi
echo
echo "--- MPI runtime ---"
if command -v mpichversion >/dev/null 2>&1; then
    mpichversion | head -2
elif command -v ompi_info >/dev/null 2>&1; then
    ompi_info | grep -iE "open mpi:|ident string" | head -2
elif command -v mpiname >/dev/null 2>&1; then
    mpiname -a 2>&1 | head -2
else
    echo "no known MPI version probe in PATH (informational only)"
fi
echo
echo "--- scontrol MpiDefault (informational) ---"
scontrol show config 2>/dev/null | grep -i MpiDefault || echo "scontrol not available"
echo
echo "--- mpi4py ---"
python3 -c "from mpi4py import MPI; print('mpi4py OK, library_version=', MPI.Get_library_version().split('\n')[0])" \
    || { echo "FAIL: mpi4py import failed; install with: pip install --no-binary=mpi4py mpi4py"; exit 1; }

echo
echo "PASS: srun + --overlap + mpi4py are all available"
exit 0
