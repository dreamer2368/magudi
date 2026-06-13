/* Standalone nested-srun smoke test, worker side.
 *
 * Built by `make` (sibling Makefile). Run by test_subprocess.py via
 * `srun --overlap -n N ./worker` (or `mpirun -n N ./worker` in non-SLURM
 * environments via the LAUNCHER env var).
 *
 * Unlike spawn_child in ../slurm_check, this binary has no parent
 * inter-communicator -- it's a fresh srun job step, not an MPI_Comm_spawn
 * target. So there's no disconnectParentIfSpawned-style handshake; the
 * parent Python waits on subprocess.run() returning, not on an MPI Barrier.
 */
#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    printf("hello from rank %d / %d\n", rank, size);
    fflush(stdout);

    MPI_Finalize();
    return 0;
}
