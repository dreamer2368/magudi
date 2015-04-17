#include "config.h"

program patch_collectives

  use MPI
  use, intrinsic :: iso_fortran_env, only : real64

  use Grid_mod, only : t_Grid
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  real(real64), parameter :: testDuration = 2.0_real64
  integer :: i, direction, nDimensions, gridSize(3,1), extent(6), processDistribution(3),    &
       errorCode, color, patchCommunicator = MPI_COMM_NULL, procRank, numProcs, ierror
  logical :: success, success_
  real(real64) :: startTime
  character(len = STRING_LENGTH) :: str

  type(t_SimulationFlags) :: simulationFlags
  type(t_SolverOptions) :: solverOptions
  type(t_Grid) :: grid
  type(t_PatchDescriptor) :: patchDescriptor
  type(t_BlockInterfacePatch) :: patch

  interface

     subroutine testScatterAndGatherBack(patch, success)

       use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

       type(t_BlockInterfacePatch) :: patch
       logical, intent(out) :: success

     end subroutine testScatterAndGatherBack

  end interface

  call MPI_Init(ierror)

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  startTime = MPI_Wtime()

  do

     do nDimensions = 1, 3

        direction = random(1, nDimensions)
        if (random(0, 1) == 0) direction = -direction
        call MPI_Bcast(direction, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

        call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)

        gridSize(:,:) = 1
        do i = 1, nDimensions
           gridSize(i,1) = random(2 * numProcs, max(2 * numProcs, 60))
        end do
        call MPI_Bcast(gridSize, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

        processDistribution = 1
        where (gridSize(:,1) > 1)
           processDistribution = 0
        end where

        call solverOptions%initialize(nDimensions, simulationFlags)
        call grid%setup(1, gridSize(1:nDimensions,1), MPI_COMM_WORLD,                        &
             processDistribution = processDistribution(1:nDimensions),                       &
             simulationFlags = simulationFlags)
        call grid%setupSpatialDiscretization(simulationFlags, solverOptions)

        call MPI_Comm_rank(grid%comm, procRank, ierror)

        extent(:) = 1
        do i = 1, nDimensions
           extent(1+2*(i-1)) = random(1, gridSize(i,1) - 1)
           extent(2+2*(i-1)) = random(extent(1+2*(i-1)) + 1, gridSize(i,1))
        end do
        extent(2+2*(abs(direction)-1)) = extent(1+2*(abs(direction)-1))
        call MPI_Bcast(extent, 6, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

        patchDescriptor = t_PatchDescriptor("testPatch", "SAT_BLOCK_INTERFACE", 1,           &
             direction, extent(1), extent(2), extent(3), extent(4), extent(5), extent(6))

        call patchDescriptor%validate(gridSize, simulationFlags,                             &
             solverOptions, errorCode, str)
        success = success .and. (errorCode == 0)
        call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL,                            &
             MPI_LAND, MPI_COMM_WORLD, ierror)
        if (.not. success) exit

        color = 1
        if (extent(2) < grid%offset(1) + 1 .or.                                              &
             extent(1) > grid%offset(1) + grid%localSize(1) .or.                             &
             extent(4) < grid%offset(2) + 1 .or.                                             &
             extent(3) > grid%offset(2) + grid%localSize(2) .or.                             &
             extent(6) < grid%offset(3) + 1 .or.                                             &
             extent(5) > grid%offset(3) + grid%localSize(3)) then
           color = MPI_UNDEFINED
        end if
        call MPI_Comm_split(grid%comm, color, procRank, patchCommunicator, ierror)
        call patch%setup(1, patchCommunicator, patchDescriptor,                              &
             grid, simulationFlags, solverOptions)

        call testScatterAndGatherBack(patch, success_)
        success = success .and. success_
        call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL,                            &
             MPI_LAND, MPI_COMM_WORLD, ierror)
        if (.not. success) exit

        if (patchCommunicator /= MPI_COMM_NULL) call MPI_Comm_free(patchCommunicator, ierror)
        patchCommunicator = MPI_COMM_NULL

     end do !... nDimensions = 1, 3

     if (.not. success) exit
     call MPI_Barrier(MPI_COMM_WORLD, ierror)
     if (MPI_Wtime() - startTime > testDuration) exit !... run for `testDuration` seconds

  end do

  call cleanupErrorHandler()

  call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL,                                  &
       MPI_LAND, MPI_COMM_WORLD, ierror)
  call MPI_Finalize(ierror)
  if (.not. success) stop -1
  stop 0

end program patch_collectives

subroutine testScatterAndGatherBack(patch, success)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  ! <<< Internal modules >>>
  use RandomNumber, only : random

  implicit none

  ! <<< Arguments >>>
  type(t_BlockInterfacePatch) :: patch
  logical, intent(out) :: success

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(SCALAR_KIND), allocatable :: f(:,:,:)
  SCALAR_TYPE, allocatable :: fLocal1(:), fGlobal1(:), fLocal2(:,:), fGlobal2(:,:), &
       fLocal3(:,:,:), fGlobal3(:,:,:)
  integer :: n, m, procRank, ierror

  success = .true.

  if (patch%comm == MPI_COMM_NULL) return

  n = random(1, 5)
  m = random(1, 5)
  call MPI_Bcast(n, 1, MPI_INTEGER, 0, patch%comm, ierror)
  call MPI_Bcast(m, 1, MPI_INTEGER, 0, patch%comm, ierror)

  call MPI_Comm_rank(patch%comm, procRank, ierror)

  if (procRank == 0) then
     allocate(fGlobal1(product(patch%globalSize)))
     allocate(fGlobal2(product(patch%globalSize), n))
     allocate(fGlobal3(product(patch%globalSize), n, m))
     allocate(f(product(patch%globalSize), n, m))
     call random_number(f)
     fGlobal1 = f(:,1,1)
     fGlobal2 = f(:,:,1)
     fGlobal3 = f
  end if

  if (patch%nPatchPoints > 0) then
     allocate(fLocal1(patch%nPatchPoints))
     allocate(fLocal2(patch%nPatchPoints, n))
     allocate(fLocal3(patch%nPatchPoints, n, m))
  end if

  call patch%scatterData(fGlobal1, fLocal1)
  if (procRank == 0) fGlobal1 = 0.0_wp
  call patch%gatherData(fLocal1, fGlobal1)
  if (procRank == 0) success = success .and. (all(abs(fGlobal1 - f(:,1,1)) < epsilon(0.0_wp)))

  call patch%scatterData(fGlobal2, fLocal2)
  if (procRank == 0) fGlobal2 = 0.0_wp
  call patch%gatherData(fLocal2, fGlobal2)
  if (procRank == 0) success = success .and. (all(abs(fGlobal2 - f(:,:,1)) < epsilon(0.0_wp)))

  call patch%scatterData(fGlobal3, fLocal3)
  if (procRank == 0) fGlobal3 = 0.0_wp
  call patch%gatherData(fLocal3, fGlobal3)
  if (procRank == 0) success = success .and. (all(abs(fGlobal3 - f) < epsilon(0.0_wp)))

  SAFE_DEALLOCATE(fLocal3)
  SAFE_DEALLOCATE(fLocal2)
  SAFE_DEALLOCATE(fLocal1)
  SAFE_DEALLOCATE(f)
  SAFE_DEALLOCATE(fGlobal3)
  SAFE_DEALLOCATE(fGlobal2)
  SAFE_DEALLOCATE(fGlobal1)

  call MPI_Bcast(success, 1, MPI_LOGICAL, 0, patch%comm, ierror)

end subroutine testScatterAndGatherBack
