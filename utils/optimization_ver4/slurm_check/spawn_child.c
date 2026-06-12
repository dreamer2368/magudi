/* Standalone PMIx spawn smoke test, child side.
 *
 * Built by `make` (sibling Makefile). Run only via MPI_Comm_spawn from
 * test_spawn.py -- not directly.
 *
 * Mirrors src/MPIHelperImpl.f90:disconnectParentIfSpawned (lines 552-574):
 * if a parent inter-communicator exists, synchronize with the parent via a
 * Barrier (true rendezvous) then Disconnect (returns immediately when no
 * traffic has flowed). This makes the parent's inter.Barrier() the actual
 * wait point for child completion.
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

    MPI_Comm parent;
    MPI_Comm_get_parent(&parent);
    if (parent != MPI_COMM_NULL) {
        MPI_Barrier(parent);
        MPI_Comm_disconnect(&parent);
    }

    MPI_Finalize();
    return 0;
}
