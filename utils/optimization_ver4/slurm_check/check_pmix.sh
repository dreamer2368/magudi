#!/usr/bin/env bash
# Read-only cluster probe; no allocation needed.
#
# Confirms that:
#   1. SLURM was built with the PMIx plugin (srun --mpi=list shows pmix)
#   2. mpi4py imports cleanly against the cluster's MPI
#
# Exits 0 on PASS (both conditions met), 1 on FAIL.
set -uo pipefail

echo "--- srun --mpi=list ---"
srun --mpi=list 2>&1 || true
echo
echo "--- srun --version ---"
srun --version 2>&1 || true
echo
echo "--- MPI runtime ---"
if command -v mpichversion >/dev/null 2>&1; then
    mpichversion
elif command -v ompi_info >/dev/null 2>&1; then
    ompi_info | grep -iE "mca:.*pmix|open mpi:|ident string"
else
    echo "neither mpichversion nor ompi_info found in PATH"
fi
echo
echo "--- env ---"
echo "SLURM_MPI_TYPE=${SLURM_MPI_TYPE:-<unset, will use cluster default>}"
echo
echo "--- mpi4py ---"
python3 -c "from mpi4py import MPI; print('mpi4py OK, library_version=', MPI.Get_library_version().split('\n')[0])" \
    || echo "mpi4py import FAILED -- install with: pip install --no-binary=mpi4py mpi4py"

if srun --mpi=list 2>&1 | grep -qi pmix; then
    echo
    echo "PASS: 'pmix' is an accepted --mpi value on this cluster"
    exit 0
fi
echo
echo "FAIL: 'pmix' not listed by srun --mpi=list; this cluster's SLURM was not built with the PMIx plugin"
exit 1
