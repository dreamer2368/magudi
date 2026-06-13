# Nested-srun (subprocess) spawn smoke test

Standalone validation of the **nested-srun** launch pattern: a parent
Python process launched by `srun` invokes `subprocess.run(["srun",
"--overlap", "-n", N, "./worker"])` from rank 0 to spawn worker tasks
as a new SLURM job step inside the same allocation.

Use this when the PMIx-spawn smoke test in [`../slurm_check/`](../slurm_check/)
fails — typically because the cluster MPI (e.g. MVAPICH2 on PSM) does
not implement MPI-2 `MPI_Comm_spawn`, even though SLURM's PMIx plugin
is present. A PASS here green-lights the deferred `--mode slurm` patch
to `optim.parallel.py` / `msgrad_test.py`, which would swap their
`_spawn_and_wait` for a subprocess-based variant whose function
signature already matches (see `launch_and_wait` in `test_subprocess.py`).

The directory is **self-contained**: no magudi build, no PETSc, no
Fortran. Only `mpicc` for the worker binary and Python+mpi4py for the
parent.

## Files

- `check_srun.sh` — read-only probe: `srun --version`, `--overlap`
  support, MPI runtime, `MpiDefault`, `mpi4py` import.
- `Makefile` — builds `worker` via `MPICC` (default `mpicc`).
- `worker.c` — MPI binary: Init / `printf("hello from rank %d / %d")` /
  Finalize. No parent-comm handshake (no `MPI_Comm_spawn` involved).
- `test_subprocess.py` — Python+mpi4py parent. All N_PETSC ranks
  Barrier; rank 0 `subprocess.run`s `srun --overlap -n N_WORKERS
  ./worker` with stdout/stderr captured to `out/worker.out`; the
  return code is broadcast; all ranks Barrier again.
- `run_test.sh` — orchestrator. Launches the parent under `srun -n
  N_PETSC` and PASS/FAIL-grades the captured worker stdout.

## Run

```
# 1. (no allocation) Cluster has srun + --overlap + mpi4py?
bash check_srun.sh

# 2. Build the worker.
make

# 3. (needs allocation) Nested srun job step works?
salloc -N 1 -n 8 bash run_test.sh
```

The allocation must have at least `max(N_PETSC, N_WORKERS)` tasks so the
worker job step has resources to land on under `--overlap`.

## Tunables (env vars for `run_test.sh`)

| Var | Default | Notes |
|---|---|---|
| `N_PETSC` | 1 | Parent rank count. |
| `N_WORKERS` | 4 | Worker rank count. |
| `TIMEOUT` | 60 | Seconds before the test is considered hung. |
| `PYLAUNCHER` | `srun` | Launcher for the parent. Set to `mpirun` for non-SLURM (e.g. dev container). |
| `LAUNCHER` | `srun --overlap` | Launcher the parent uses for the worker. Set to `mpirun` for non-SLURM. |

## Interpreting results

| `check_srun.sh` | `run_test.sh` | Meaning |
|---|---|---|
| PASS | PASS | Production OK to use nested-srun under `--mode slurm`. |
| PASS | FAIL (timeout) | The nested `srun` step hung — usually the allocation doesn't have spare tasks to land the worker step on under `--overlap`. Re-allocate with more tasks, or check cluster policy on overlapping job steps. |
| PASS | FAIL (non-zero) | `srun` rejected the command or the worker errored. Inspect `out/worker.out` and the script's stderr. |
| FAIL | — | `srun` or `--overlap` not available, or mpi4py broken. Address before running step 3. |

## Local dev-container sanity check

The dev container has no SLURM, so the test must use `mpirun` as both
launchers:

```
PYLAUNCHER=mpirun LAUNCHER=mpirun bash run_test.sh
```

This validates the Python/C wiring (parent Barriers, subprocess
launch, stdout capture, return-code broadcast) before deploying to a
real cluster.
