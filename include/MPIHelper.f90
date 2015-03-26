#include "config.h"

module MPIHelper

  implicit none
  public

  interface

     subroutine pigeonhole(nPigeons, nHoles, holeIndex, offset, nPigeonsInHole)

       !> Puts `nPigeons` pigeons in `nHoles` ordered holes, such that each hole gets an equal
       !> number of pigeons except the first `mod(nPigeons, nHoles)`, which get an extra
       !> pigeon. The number of pigeons in the first `holeIndex` holes (zero if
       !> `holeIndex=0`), and the number of pigeons in the hole with zero-based index
       !> `holeIndex` are returned as `offset` and `nPigeonsInHole`, respectively.

       integer, intent(in) :: nPigeons, nHoles, holeIndex
       integer, intent(out) :: offset, nPigeonsInHole

     end subroutine pigeonhole

  end interface

  interface

     subroutine splitCommunicatorMultigrid(comm, gridSizes,                                  &
          gridCommunicators, numProcsInGrid)

       !> Splits the MPI communicator `comm` into multiple communicators `gridCommunicators`
       !> according to the process distribution map `numProcsInGrid`, where
       !> `gridCommunicators(i)` is a grid-level communicator for the `i`-th grid. If
       !> `numProcsInGrid` is not specified, it is chosen such that the number of processes
       !> assigned to a grid is proportional to the number of grid points.

       integer, intent(in) :: comm, gridSizes(:,:)
       integer, intent(out) :: gridCommunicators(:)

       integer, intent(in), optional :: numProcsInGrid(:)

     end subroutine splitCommunicatorMultigrid

  end interface

  interface

     subroutine fillGhostPoints(comm, arrayWithGhostPoints,                              &
          direction, numGhostPoints, periodicOffset)

       !> The first `numGhostPoints(1)` and last `numGhostPoints(2)` points along dimension
       !> `direction` of `arrayWithGhostPoints` are assumed to be ``ghost'' points, and are
       !> filled with values received from the _previous_ and _next_ processes along direction
       !> `direction`, respectively. Among the remaining ``physical`` points along dimension
       !> `direction` of `arrayWithGhostPoints`, the first `numGhostPoints(1)` and last
       !> `numGhostPoints(2)` points are sent to the ``previous'' and ``next'' process along
       !> direction `direction`, respectively. Periodicity information is obtained from the
       !> MPI communicator `comm`, which is assumed to contain Cartesian topology
       !> information. If the direction `direction` is periodic and the optional argument
       !> `periodicOffset` is present, the first `periodicOffset(2)` and last
       !> `periodicOffset(1)` ``physical'' points are skipped when sending data (this is
       !> useful for a ``ring'' topology along direction `direction` with overlapping points).

       integer, intent(in) :: comm
       SCALAR_TYPE, intent(inout) :: arrayWithGhostPoints(:,:,:,:)
       integer, intent(in) :: direction, numGhostPoints(2)

       integer, intent(in), optional :: periodicOffset(2)

     end subroutine fillGhostPoints

  end interface

  interface

     subroutine gatherAlongDirection(cartesianCommunicator, localArray,                      &
          localSize, direction, offsetAlongDirection, gatheredArray)

       !> Gathers a grid vector `localArray` from all processes along direction `direction`
       !> and places them sequentially in the grid vector `gatheredArray`. The first dimension
       !> of `localArray` and `gatheredArray` are flattened from three-dimensional shapes
       !> using Fortran-ordering. `localSize` is the ``original'' shape of `localArray` prior
       !> to flattening. The ``original'' shape of `gatheredArray` is identical to
       !> `localArray` except along dimension `direction`, where it is the sum of
       !> `localSize(direction)` taken across all processes along direction
       !> `direction`. `offsetAlongDirection` is the zero-based index of the first element of
       !> `gatheredArray` along dimension `direction` received from the ``current'' process.

       integer, intent(in) :: cartesianCommunicator
       SCALAR_TYPE, intent(in) :: localArray(:,:)
       integer, intent(in) :: localSize(3), direction, offsetAlongDirection
       SCALAR_TYPE, intent(out) :: gatheredArray(:,:)

     end subroutine gatherAlongDirection

  end interface

end module MPIHelper
