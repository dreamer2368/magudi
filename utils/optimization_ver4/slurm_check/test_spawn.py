"""Standalone PMIx spawn smoke test, parent side.

Run as:
    srun --mpi=pmix -n 1 python3 test_spawn.py N

Mirrors utils/optimization_ver4/optim.parallel.py:_spawn_and_wait one-for-one
so a PASS here directly green-lights the production code path. Uses only
mpi4py -- no petsc4py, no magudi imports, no build dependency on magudi.

Build the sibling spawn_child binary with `make` before running.
"""
import sys

from mpi4py import MPI


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2

    inter = MPI.COMM_WORLD.Spawn(
        "./spawn_child",
        args=[],
        maxprocs=n,
        info=MPI.INFO_NULL,
        root=0,
    )
    inter.Barrier()
    inter.Disconnect()
    MPI.COMM_WORLD.Barrier()

    if MPI.COMM_WORLD.Get_rank() == 0:
        print(
            f"parent: spawn-and-wait completed for n={n} children",
            flush=True,
        )


if __name__ == "__main__":
    main()
