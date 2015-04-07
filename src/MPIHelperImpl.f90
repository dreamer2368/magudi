#include "config.h"

subroutine pigeonhole(nPigeons, nHoles, holeIndex, offset, nPigeonsInHole)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: nPigeons, nHoles, holeIndex
  integer, intent(out) :: offset, nPigeonsInHole

  assert(nHoles > 0)
  assert(nPigeons >= 0)
  assert(holeIndex >= 0 .and. holeIndex < nHoles)

  offset = holeIndex * (nPigeons / nHoles) + min(holeIndex, mod(nPigeons, nHoles))
  nPigeonsInHole = nPigeons / nHoles
  if (holeIndex < mod(nPigeons, nHoles)) nPigeonsInHole = nPigeonsInHole + 1

end subroutine pigeonhole

subroutine splitCommunicatorMultigrid(comm, gridSizes, gridCommunicators, numProcsInGrid)

  ! <<< External modules >>>
  use MPI

  ! <<< Public members >>>
  use MPIHelper, only : pigeonhole

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm, gridSizes(:,:)
  integer, intent(out) :: gridCommunicators(:)
  integer, intent(in), optional :: numProcsInGrid(:)

  ! <<< Local variables >>>
  integer :: i, procRank, nProcs, group, newgroup, offset, count, color, ierror
  integer, allocatable :: ranks(:), numProcsInGrid_(:)

  assert(size(gridSizes, 1) >= 1 .and. size(gridSizes, 1) <= 3)
  assert(size(gridSizes, 2) > 0)
  assert(size(gridSizes, 2) == size(gridCommunicators))

  ! Get the size, rank and group of `comm`.
  call MPI_Comm_size(comm, nProcs, ierror)
  assert(nProcs > 0)
  call MPI_Comm_rank(comm, procRank, ierror)
  call MPI_Comm_group(comm, group, ierror)

  ! Set all elements of `gridCommunicators` to `MPI_COMM_NULL` by default.
  gridCommunicators = MPI_COMM_NULL

  if (size(gridSizes, 2) > nProcs) then

     ! Partition the grids into `nProcs` processes; more than one grid may be distributed
     ! to the same process.
     call pigeonhole(size(gridSizes, 2), nProcs, procRank, offset, count)

     do i = 1, size(gridSizes, 2)
        allocate(ranks(1), source = -1)
        if (i >= offset + 1 .and. i <= offset + count) ranks = procRank
        call MPI_Allreduce(MPI_IN_PLACE, ranks, 1, MPI_INTEGER, MPI_MAX, comm, ierror)
        call MPI_Group_incl(group, 1, ranks, newgroup, ierror)
        call MPI_Comm_create(comm, newgroup, gridCommunicators(i), ierror)
        call MPI_Group_free(newgroup, ierror)
        SAFE_DEALLOCATE(ranks)
     end do

  else

     if (present(numProcsInGrid)) then

        ! Use manual distribution, if specified.
        assert(size(numProcsInGrid) == size(gridSizes, 2))
        assert(sum(numProcsInGrid) == nProcs)
        allocate(numProcsInGrid_(size(gridSizes, 2)), source = numProcsInGrid)

     else

        ! Assign processes to a grid proportional to number of grid points.
        allocate(numProcsInGrid_(size(gridSizes, 2)))
        numProcsInGrid_ = ceiling(real(nProcs) * real(product(gridSizes, dim = 1)) /         &
             real(sum(product(gridSizes, dim = 1))))
        i = 1
        do while (sum(numProcsInGrid_) /= nProcs) !... total should match
           numProcsInGrid_(i) = max(1, numProcsInGrid_(i) - 1)
           i = mod(i, size(gridSizes, 2)) + 1
        end do

     end if

     ! Split using MPI color/key mechanism:

     color = MPI_UNDEFINED
     if (procRank < numProcsInGrid_(1)) color = 1
     call MPI_Comm_split(comm, color, procRank, gridCommunicators(1), ierror)

     do i = 2, size(gridSizes, 2)
        color = MPI_UNDEFINED
        if (procRank >= sum(numProcsInGrid_(:i-1)) .and.                                     &
             procRank < sum(numProcsInGrid_(:i))) color = i
        call MPI_Comm_split(comm, color, procRank, gridCommunicators(i), ierror)
     end do

  end if

  call MPI_Group_free(group, ierror)

  SAFE_DEALLOCATE(numProcsInGrid_)

end subroutine splitCommunicatorMultigrid

subroutine fillGhostPoints(cartesianCommunicator, arrayWithGhostPoints,                      &
     direction, numGhostPoints, periodicOffset)

  ! <<< External modules >>>
  use MPI

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: cartesianCommunicator
  SCALAR_TYPE, intent(inout) :: arrayWithGhostPoints(:,:,:,:)
  integer, intent(in) :: direction, numGhostPoints(2)
  integer, intent(in), optional :: periodicOffset(2)

  ! <<< Local variables >>>
  integer :: i, j, k, l, nScalars, nDimensions, gridSizeWithGhostPoints(3), gridSize(3),     &
       normalPlaneSize, periodicOffset_(2), request(4), tagPrev, tagNext,                    &
       processDistribution(3), processCoordinates(3), rankOfPreviousProcess,                 &
       rankOfNextProcess, ierror
  logical :: isPeriodic(3)
  SCALAR_TYPE, allocatable :: receivedFromPreviousProcess(:,:,:),                            &
       receivedFromNextProcess(:,:,:),                                                       &
       sentToPreviousProcess(:,:,:),                                                         &
       sentToNextProcess(:,:,:)

  assert(all(shape(arrayWithGhostPoints) > 0))

  call startTiming("fillGhostPoints")

  ! Get number of dimensions from communicator.
  call MPI_Cartdim_get(cartesianCommunicator, nDimensions, ierror)
  assert_key(nDimensions, (1, 2, 3))
  assert(direction >= 1 .and. direction <= nDimensions)

  if (all(numGhostPoints <= 0)) return

  ! Get process distribution, coordinates and periodicity.
  processDistribution = 1
  isPeriodic = .false.
  processCoordinates = 0
  call MPI_Cart_get(cartesianCommunicator, nDimensions, processDistribution, isPeriodic,     &
       processCoordinates, ierror)
  if (processDistribution(direction) == 1 .and. numGhostPoints(1) /= numGhostPoints(2)) return

  ! Get array dimensions.
  nScalars = size(arrayWithGhostPoints, 4)
  gridSizeWithGhostPoints = shape(arrayWithGhostPoints(:,:,:,1))
  normalPlaneSize = product(gridSizeWithGhostPoints) / size(arrayWithGhostPoints, direction)
  gridSize = gridSizeWithGhostPoints - sum(numGhostPoints)
  periodicOffset_ = 0
  if (present(periodicOffset)) periodicOffset_ = periodicOffset

  !. Initialize MPI variables.
  request = MPI_REQUEST_NULL
  tagPrev = 1
  tagNext = 2
  call MPI_Cart_shift(cartesianCommunicator, direction - 1, 1,                               &
       rankOfPreviousProcess, rankOfNextProcess, ierror)

  if (numGhostPoints(1) > 0) then

     ! Allocate temporary buffers which hold send/receive data.
     allocate(receivedFromPreviousProcess(numGhostPoints(1), normalPlaneSize, nScalars))
     allocate(sentToPreviousProcess(numGhostPoints(1), normalPlaneSize, nScalars))

     ! Non-blocking receive.
     call MPI_Irecv(receivedFromPreviousProcess, size(receivedFromPreviousProcess),          &
          SCALAR_TYPE_MPI, rankOfPreviousProcess, tagNext, cartesianCommunicator,            &
          request(1), ierror)

  end if

  if (numGhostPoints(2) > 0) then

     ! Allocate temporary buffers which hold send/receive data.
     allocate(receivedFromNextProcess(numGhostPoints(2), normalPlaneSize, nScalars))
     allocate(sentToNextProcess(numGhostPoints(2), normalPlaneSize, nScalars))

     ! Non-blocking receive.
     call MPI_Irecv(receivedFromNextProcess, size(receivedFromNextProcess), SCALAR_TYPE_MPI, &
          rankOfNextProcess, tagPrev, cartesianCommunicator, request(2), ierror)

  end if

  ! Copy ghost point data to send buffers:

  select case (direction)

  case (1)

     do l = 1, nScalars
        do k = 1, gridSizeWithGhostPoints(3)
           do j = 1, gridSizeWithGhostPoints(2)
              do i = 1, numGhostPoints(1)
                 sentToPreviousProcess(i, j + gridSizeWithGhostPoints(2) * (k - 1), l) =     &
                      arrayWithGhostPoints(i + numGhostPoints(1) +                           &
                      periodicOffset_(2), j, k, l)
              end do
           end do
        end do
     end do

     do l = 1, nScalars
        do k = 1, gridSizeWithGhostPoints(3)
           do j = 1, gridSizeWithGhostPoints(2)
              do i = 1, numGhostPoints(2)
                 sentToNextProcess(numGhostPoints(2) + 1 - i, j +                            &
                      gridSizeWithGhostPoints(2) * (k - 1), l) =                             &
                      arrayWithGhostPoints(gridSize(1) + 1 - i + numGhostPoints(1) -         &
                      periodicOffset_(1), j, k, l)
              end do
           end do
        end do
     end do

  case (2)

     do l = 1, nScalars
        do k = 1, gridSizeWithGhostPoints(3)
           do i = 1, gridSizeWithGhostPoints(1)
              do j = 1, numGhostPoints(1)
                 sentToPreviousProcess(j, i + gridSizeWithGhostPoints(1) * (k - 1), l) =     &
                      arrayWithGhostPoints(i, j + numGhostPoints(1) +                        &
                      periodicOffset_(2), k, l)
              end do
           end do
        end do
     end do

     do l = 1, nScalars
        do k = 1, gridSizeWithGhostPoints(3)
           do i = 1, gridSizeWithGhostPoints(1)
              do j = 1, numGhostPoints(2)
                 sentToNextProcess(numGhostPoints(2) + 1 - j, i +                            &
                 gridSizeWithGhostPoints(1) * (k - 1), l) =                                  &
                      arrayWithGhostPoints(i, gridSize(2) + 1 - j + numGhostPoints(1) -      &
                      periodicOffset_(1), k, l)
              end do
           end do
        end do
     end do

  case (3)

     do l = 1, nScalars
        do j = 1, gridSizeWithGhostPoints(2)
           do i = 1, gridSizeWithGhostPoints(1)
              do k = 1, numGhostPoints(1)
                 sentToPreviousProcess(k, i + gridSizeWithGhostPoints(1) * (j - 1), l) =     &
                      arrayWithGhostPoints(i, j, k + numGhostPoints(1) +                     &
                      periodicOffset_(2), l)
              end do
           end do
        end do
     end do

     do l = 1, nScalars
        do j = 1, gridSizeWithGhostPoints(2)
           do i = 1, gridSizeWithGhostPoints(1)
              do k = 1, numGhostPoints(2)
                 sentToNextProcess(numGhostPoints(2) + 1 - k, i +                            &
                      gridSizeWithGhostPoints(1) * (j - 1), l) =                             &
                      arrayWithGhostPoints(i, j, gridSize(3) + 1 - k + numGhostPoints(1)     &
                      - periodicOffset_(1), l)
              end do
           end do
        end do
     end do

  end select

  ! Non-blocking send followed by `MPI_Waitall`.
  if (numGhostPoints(2) > 0) then
     call MPI_Isend(sentToNextProcess, size(sentToNextProcess), SCALAR_TYPE_MPI,             &
          rankOfNextProcess, tagNext, cartesianCommunicator, request(3), ierror)
  end if
  if (numGhostPoints(1) > 0) then
     call MPI_Isend(sentToPreviousProcess, size(sentToPreviousProcess), SCALAR_TYPE_MPI,     &
       rankOfPreviousProcess, tagPrev, cartesianCommunicator, request(4), ierror)
  end if
  call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierror)

  ! Copy receive buffer data to ghost points:

  select case (direction)

  case (1)

     do l = 1, nScalars
        do k = 1, gridSizeWithGhostPoints(3)
           do j = 1, gridSizeWithGhostPoints(2)
              do i = 1, numGhostPoints(1)
                 arrayWithGhostPoints(i,j,k,l) =                                             &
                      receivedFromPreviousProcess(i, j +                                     &
                      gridSizeWithGhostPoints(2) * (k - 1), l)
              end do
           end do
        end do
     end do

     do l = 1, nScalars
        do k = 1, gridSizeWithGhostPoints(3)
           do j = 1, gridSizeWithGhostPoints(2)
              do i = 1, numGhostPoints(2)
                 arrayWithGhostPoints(gridSizeWithGhostPoints(1) + 1 - i, j, k, l) =         &
                      receivedFromNextProcess(numGhostPoints(2) + 1 - i,                     &
                      j + gridSizeWithGhostPoints(2) * (k - 1), l)
              end do
           end do
        end do
     end do

  case (2)

     do l = 1, nScalars
        do k = 1, gridSizeWithGhostPoints(3)
           do i = 1, gridSizeWithGhostPoints(1)
              do j = 1, numGhostPoints(1)
                 arrayWithGhostPoints(i,j,k,l) =                                             &
                      receivedFromPreviousProcess(j, i +                                     &
                      gridSizeWithGhostPoints(1) * (k - 1), l)
              end do
           end do
        end do
     end do

     do l = 1, nScalars
        do k = 1, gridSizeWithGhostPoints(3)
           do i = 1, gridSizeWithGhostPoints(1)
              do j = 1, numGhostPoints(2)
                 arrayWithGhostPoints(i, gridSizeWithGhostPoints(2) + 1 - j, k, l) =         &
                      receivedFromNextProcess(numGhostPoints(2) + 1 - j,                     &
                      i + gridSizeWithGhostPoints(1) * (k - 1), l)
              end do
           end do
        end do
     end do

  case (3)

     do l = 1, nScalars
        do j = 1, gridSizeWithGhostPoints(2)
           do i = 1, gridSizeWithGhostPoints(1)
              do k = 1, numGhostPoints(1)
                 arrayWithGhostPoints(i,j,k,l) =                                             &
                      receivedFromPreviousProcess(k, i +                                     &
                      gridSizeWithGhostPoints(1) * (j - 1), l)
              end do
           end do
        end do
     end do

     do l = 1, nScalars
        do j = 1, gridSizeWithGhostPoints(2)
           do i = 1, gridSizeWithGhostPoints(1)
              do k = 1, numGhostPoints(2)
                 arrayWithGhostPoints(i, j, gridSizeWithGhostPoints(3) + 1 - k, l) =         &
                      receivedFromNextProcess(numGhostPoints(2) + 1 - k,                     &
                      i + gridSizeWithGhostPoints(1) * (j - 1), l)
              end do
           end do
        end do
     end do

  end select

  SAFE_DEALLOCATE(sentToPreviousProcess)
  SAFE_DEALLOCATE(sentToNextProcess)
  SAFE_DEALLOCATE(receivedFromPreviousProcess)
  SAFE_DEALLOCATE(receivedFromNextProcess)

  call endTiming("fillGhostPoints")

end subroutine fillGhostPoints

subroutine gatherAlongDirection(cartesianCommunicator, localArray,                           &
     localSize, direction, offsetAlongDirection, gatheredArray)

  ! <<< External modules >>>
  use MPI

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: cartesianCommunicator
  SCALAR_TYPE, intent(in) :: localArray(:,:)
  integer, intent(in) :: localSize(3), direction, offsetAlongDirection
  SCALAR_TYPE, intent(out) :: gatheredArray(:,:)

  ! <<< Local variables >>>
  integer :: i, j, k, l, nScalars, globalSizeAlongDirection, pencilCommunicator,             &
       numProcsInPencil, nDimensions, ierror
  logical :: keepThisDirection(3)
  integer, allocatable :: localSizes(:), offsets(:)
  SCALAR_TYPE, allocatable :: sendBuffer(:,:), receiveBuffer(:,:)

  assert(all(localSize > 0))
  assert(size(localArray, 1) == product(localSize))
  assert(size(gatheredArray, 1) > 0)
  assert(size(localArray, 2) == size(gatheredArray, 2))
  assert(mod(size(gatheredArray, 1), size(localArray, 1) / localSize(direction)) == 0)

  ! Get the dimensions of the local and gathered arrays.
  globalSizeAlongDirection = size(gatheredArray, 1) /                                        &
       (size(localArray, 1) / localSize(direction))
  nScalars = size(localArray, 2)

  ! Get number of dimensions from the grid communicator.
  call MPI_Cartdim_get(cartesianCommunicator, nDimensions, ierror)
  assert_key(nDimensions, (1, 2, 3))
  assert(direction >= 1 .and. direction <= nDimensions)

  ! Create a pencil communicator including the ``current'' process and all neighboring
  ! processes along direction `direction`.
  keepThisDirection = .false.; keepThisDirection(direction) = .true.
  call MPI_Cart_sub(cartesianCommunicator, keepThisDirection(1:nDimensions),                 &
       pencilCommunicator, ierror)

  ! Gather the offsets and sizes of `localArray` along dimension `direction` of all processes
  ! in the pencil communicator.
  call MPI_Comm_size(pencilCommunicator, numProcsInPencil, ierror)
  allocate(localSizes(numProcsInPencil), offsets(numProcsInPencil))
  call MPI_Allgather(localSize(direction), 1, MPI_INTEGER,                                   &
       localSizes, 1, MPI_INTEGER, pencilCommunicator, ierror)
  call MPI_Allgather(offsetAlongDirection, 1, MPI_INTEGER, offsets,                          &
       1, MPI_INTEGER, pencilCommunicator, ierror)

  ! Find the sizes and displacements of data to be gathered.
  localSizes = localSizes * size(localArray) / localSize(direction)
  offsets = offsets * size(localArray) / localSize(direction)

  ! Allocate send and receive buffers with shapes that can be gathered easily using
  ! `MPI_Allgatherv`.
  allocate(sendBuffer(nScalars * size(localArray) / localSize(direction),                    &
       localSize(direction)))
  allocate(receiveBuffer(nScalars * size(localArray) / localSize(direction),                 &
       globalSizeAlongDirection))

  ! Pack the data to be gathered:

  select case (direction)

  case (1)
     do l = 1, nScalars
        do k = 1, localSize(3)
           do j = 1, localSize(2)
              do i = 1, localSize(1)
                 sendBuffer(j + localSize(2) * (k - 1 + localSize(3) * (l - 1)), i) =        &
                      localArray(i + localSize(1) * (j - 1 + localSize(2) * (k - 1)), l)
              end do
           end do
        end do
     end do

  case (2)
     do l = 1, nScalars
        do k = 1, localSize(3)
           do j = 1, localSize(2)
              do i = 1, localSize(1)
                 sendBuffer(i + localSize(1) * (k - 1 + localSize(3) * (l - 1)), j) =        &
                      localArray(i + localSize(1) * (j - 1 + localSize(2) * (k - 1)), l)
              end do
           end do
        end do
     end do

  case (3)
     do l = 1, nScalars
        do k = 1, localSize(3)
           do j = 1, localSize(2)
              do i = 1, localSize(1)
                 sendBuffer(i + localSize(1) * (j - 1 + localSize(2) * (l - 1)), k) =        &
                      localArray(i + localSize(1) * (j - 1 + localSize(2) * (k - 1)), l)
              end do
           end do
        end do
     end do

  end select

  call MPI_Allgatherv(sendBuffer, size(sendBuffer), SCALAR_TYPE_MPI, receiveBuffer,          &
       localSizes, offsets, SCALAR_TYPE_MPI, pencilCommunicator, ierror)

  ! Unpack the data that was gathered:

  select case (direction)

  case (1)
     do l = 1, nScalars
        do k = 1, localSize(3)
           do j = 1, localSize(2)
              do i = 1, globalSizeAlongDirection
                 gatheredArray(i + globalSizeAlongDirection * (j - 1 +                       &
                      localSize(2) * (k - 1)), l) =                                          &
                      receiveBuffer(j + localSize(2) * (k - 1 + localSize(3) * (l - 1)), i)
              end do
           end do
        end do
     end do

  case (2)
     do l = 1, nScalars
        do k = 1, localSize(3)
           do j = 1, globalSizeAlongDirection
              do i = 1, localSize(1)
                 gatheredArray(i + localSize(1) * (j - 1 +                                   &
                      globalSizeAlongDirection * (k - 1)), l) =                              &
                      receiveBuffer(i + localSize(1) * (k - 1 + localSize(3) * (l - 1)), j)
              end do
           end do
        end do
     end do

  case (3)
     do l = 1, nScalars
        do k = 1, globalSizeAlongDirection
           do j = 1, localSize(2)
              do i = 1, localSize(1)
                 gatheredArray(i + localSize(1) * (j - 1 + localSize(2) * (k - 1)), l) =     &
                      receiveBuffer(i + localSize(1) * (j - 1 + localSize(2) * (l - 1)), k)
              end do
           end do
        end do
     end do

  end select

  call MPI_Comm_free(pencilCommunicator, ierror)

  SAFE_DEALLOCATE(receiveBuffer)
  SAFE_DEALLOCATE(sendBuffer)
  SAFE_DEALLOCATE(offsets)
  SAFE_DEALLOCATE(localSizes)

end subroutine gatherAlongDirection
