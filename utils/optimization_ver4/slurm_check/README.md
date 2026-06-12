# SLURM / PMIx spawn smoke tests

Standalone validation that the target cluster's SLURM PMIx plugin
implements MPI-2 dynamic process management (`MPI_Comm_spawn`). A PASS here
green-lights running `utils/optimization_ver4/optim.parallel.py` under
`srun --mpi=pmix`; a FAIL means the production code needs a different
worker-launch design (nested-srun via subprocess).

The directory is **self-contained**: no magudi build, no PETSc, no Fortran
required. Only `mpicc` for the child binary and Python+mpi4py for the
parent.

## Files

- `check_pmix.sh` — read-only cluster probe; no allocation needed.
- `Makefile` — builds `spawn_child` via `MPICC` (default `mpicc`).
- `spawn_child.c` — MPI child binary; mirrors
  `src/MPIHelperImpl.f90:disconnectParentIfSpawned` (parent Barrier +
  Disconnect, gated on `MPI_Comm_get_parent != MPI_COMM_NULL`).
- `test_spawn.py` — Python+mpi4py parent that calls `MPI.COMM_WORLD.Spawn`
  and runs the exact handshake used by `optim.parallel.py:_spawn_and_wait`.
- `run_test.sh` — orchestrates `timeout 60 srun --mpi=pmix -n 1 python3
  test_spawn.py 4` and PASS/FAIL-grades the output.

## Run

```
# 1. Cluster has PMIx and mpi4py imports?  (no allocation)
bash check_pmix.sh

# 2. Build the C child.
make

# 3. Spawn handshake works under srun?  (needs an allocation)
salloc -N 1 -n 8 bash run_test.sh
```

If `mpi4py` is not yet on the cluster Python, install it against the
cluster MPI with:

```
pip install --no-binary=mpi4py mpi4py
```

(matches the recipe in the dev container's `docker/Dockerfile`).

## Interpreting results

| `check_pmix.sh` | `run_test.sh` | Meaning |
|---|---|---|
| PASS | PASS | Production OK to run `optim.parallel.py` under `srun --mpi=pmix`. |
| PASS | FAIL (timeout) | PMIx accepts the spawn call but the inter-comm Barrier never returns. The cluster's PMIx is incomplete; switch the production launcher to nested-srun. |
| PASS | FAIL (non-zero) | mpi4py or the spawn call errored before the handshake. Inspect the captured stdout. |
| FAIL | — | Cluster's SLURM was not built with PMIx. Switch the production launcher to nested-srun (or use a different cluster). |
