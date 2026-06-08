"""Single-shot spawner: launch ./forward with N children via _spawn_and_wait.

Run as:
    mpirun -n 1 python3 spawn_one_forward.py N

The body of _spawn_and_wait is copied from utils/optimization_ver4/optim.parallel.py
to keep this comparison script self-contained (no petsc4py / parallel_io import).
"""
import os
import shlex
import sys

from mpi4py import MPI


def _spawn_and_wait(executable, args, n):
    comm = MPI.COMM_WORLD
    log_path = os.path.abspath(
        os.path.join("out", os.path.basename(executable) + ".out")
    )
    if comm.Get_rank() == 0:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        open(log_path, "w").close()
    comm.Barrier()

    quoted_args = " ".join(shlex.quote(a) for a in args)
    redir_cmd = f"exec {shlex.quote(executable)} {quoted_args}"

    inter = MPI.COMM_WORLD.Spawn(
        "/bin/bash",
        args=["-c", redir_cmd],
        maxprocs=n,
        info=MPI.INFO_NULL,
        root=0,
    )
    inter.Barrier()
    inter.Disconnect()
    MPI.COMM_WORLD.Barrier()


if __name__ == "__main__":
    n = int(sys.argv[1])
    _spawn_and_wait("./forward", ["--input", "magudi.inp"], n)
