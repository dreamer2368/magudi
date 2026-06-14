# magudi_optimizer

PETSc/TAO-based PDE-constrained optimization driver for [`magudi`](../../README.md), the
MPI-parallel compressible Navier-Stokes solver. This is the `ver4` production framework;
it replaces the per-segment bash orchestration in `utils/optimization_ver3/` with a single
TAO L-BFGS outer loop running on top of the multi-segment Fortran binaries `msforward`
and `msadjoint` (`bin/msforward.f90`, `bin/msadjoint.f90`).

See [`DESIGN.md`](DESIGN.md) for the architecture (problem-scale constraints, variable
transformation, the Python ↔ Fortran interface, async-evaluation cases).

## Install

From the repository root:

```
pip install ./utils/optimization_ver4
```

## Console scripts

| Command | Purpose |
|---|---|
| `magudi-optim`  | TAO L-BFGS outer loop driving `./msforward` + `./msadjoint`. |
| `magudi-msgrad` | Finite-difference gradient-accuracy test for `msforward` / `msadjoint`. |

Both are launched under MPI:

```
mpirun -n N_petsc magudi-optim optim.yml --max-iter 5
mpirun -n N_petsc magudi-msgrad msgrad.yml --mode full
```

`N_petsc` must satisfy the `ParallelIOHandler` resource policy (see
`src/magudi_optimizer/parallel_io.py`).

## CI driver scripts

End-to-end test wrappers live next to the GitHub Actions workflow:

- [`.github/workflows/run_parallel.sh`](../../.github/workflows/run_parallel.sh) —
  restart-equivalence check via the `OneDWave` example.
- [`.github/workflows/run_msgrad.sh`](../../.github/workflows/run_msgrad.sh) —
  multi-segment gradient-accuracy check via the `OneDWave` example.

Both run from `build/` after `make -j` has produced
`bin/{forward,msforward,msadjoint,compute_norm}`.