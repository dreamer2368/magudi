"""Standalone nested-srun smoke test, parent side.

Run as:
    srun -n N_PETSC python3 test_subprocess.py N_WORKERS

All N_PETSC parent ranks Barrier; rank 0 calls subprocess.run on
`srun --overlap -n N_WORKERS ./worker`, blocks until it returns; the
return code is broadcast so every rank either continues or raises
together; then all ranks Barrier again.

This is the candidate drop-in replacement for `_spawn_and_wait` under a
future `--mode slurm` patch to optim.parallel.py / msgrad_test.py.
Function signature is intentionally identical to that helper so the
swap is mechanical.

Env vars:
  LAUNCHER  child launcher (default "srun --overlap"). Set to "mpirun"
            for non-SLURM environments (e.g. the dev container).
"""
import os
import shlex
import subprocess
import sys

from mpi4py import MPI


def launch_and_wait(executable, args, n):
    """Collectively launch `n` worker processes via the LAUNCHER; block."""
    comm = MPI.COMM_WORLD
    log_path = os.path.abspath(
        os.path.join("out", os.path.basename(executable) + ".out")
    )
    if comm.Get_rank() == 0:
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        open(log_path, "w").close()
    comm.Barrier()

    if comm.Get_rank() == 0:
        launcher = shlex.split(os.environ.get("LAUNCHER", "srun --overlap"))
        cmd = launcher + ["-n", str(n), executable] + list(args)
        with open(log_path, "a") as logf:
            rc = subprocess.run(
                cmd, stdout=logf, stderr=subprocess.STDOUT
            ).returncode
    else:
        rc = None
    rc = comm.bcast(rc, root=0)
    comm.Barrier()
    if rc != 0:
        raise SystemExit(f"worker exited with code {rc}")


def main():
    n_workers = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    launch_and_wait("./worker", [], n_workers)

    if MPI.COMM_WORLD.Get_rank() == 0:
        print(
            f"parent: launch-and-wait completed for n={n_workers} workers",
            flush=True,
        )


if __name__ == "__main__":
    main()
