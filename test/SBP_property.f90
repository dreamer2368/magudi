#include "config.h"

program SBP_property

  use MPI
  use, intrinsic :: iso_fortran_env, only : real64

  use MPIHelper, only : pigeonhole
  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random
  use StencilOperator_mod, only : t_StencilOperator

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  real(real64) :: testDuration = 2.0_real64
  character(len = STRING_LENGTH), parameter :: discretizationTypes(*) =                      &
       (/ "SBP 1-2", "SBP 2-4", "SBP 3-6", "SBP 4-8" /)
  integer :: i, n, direction, nDimensions, cartesianCommunicator, numProcs, procRank, ierror
  integer, dimension(:), allocatable :: localSize, offset, globalSize,                       &
       processDistribution, processCoordinates
  logical :: success
  logical, allocatable :: isPeriodic(:)
  real(real64) :: startTime
  character(len = STRING_LENGTH) :: str
  type(t_StencilOperator) :: D, adjointOfD

  interface

     function verifyFirstDerivativeSBP(D, adjointOfD,                                        &
          localSize, tolerance) result(satisfiesSBP)

       use StencilOperator_mod, only : t_StencilOperator

       type(t_StencilOperator), intent(in) :: D, adjointOfD
       integer, intent(in) :: localSize(:)

       real(SCALAR_KIND), intent(in), optional :: tolerance

       logical :: satisfiesSBP

     end function verifyFirstDerivativeSBP

  end interface

  call MPI_Init(ierror)

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)

  startTime = MPI_Wtime()

  do

     str = trim(discretizationTypes(random(1, size(discretizationTypes))))
     call MPI_Bcast(str, len(str), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)

     call D%setup(trim(str) // " first derivative")
     call D%getAdjoint(adjointOfD)

     do nDimensions = 1, 3

        allocate(localSize(nDimensions))
        allocate(offset(nDimensions))
        allocate(globalSize(nDimensions))
        allocate(processDistribution(nDimensions))
        allocate(processCoordinates(nDimensions))
        allocate(isPeriodic(nDimensions))

        direction = random(1, nDimensions)
        call MPI_Bcast(direction, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

        do i = 1, nDimensions
           isPeriodic(i) = (random(0, 1) == 0)
        end do
        call MPI_Bcast(isPeriodic, nDimensions, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)

        processDistribution = 0
        call MPI_Dims_create(numProcs, nDimensions, processDistribution, ierror)
        call MPI_Cart_create(MPI_COMM_WORLD, nDimensions, processDistribution,               &
             isPeriodic, .true., cartesianCommunicator, ierror)
        call D%update(cartesianCommunicator, direction)
        call adjointOfD%update(cartesianCommunicator, direction)

        do i = 1, nDimensions
           if (i == direction) then
              n = numProcs * (2 * adjointOfD%boundaryWidth + 1)
              globalSize(i) = random(n, max(2 * n, 2 ** 6))
           else
              globalSize(i) = random(n, max(2 * n, 2 ** 6))
           end if
        end do
        call MPI_Bcast(globalSize, nDimensions, MPI_INTEGER, 0, cartesianCommunicator, ierror)

        call MPI_Comm_rank(cartesianCommunicator, procRank, ierror)
        call MPI_Cart_coords(cartesianCommunicator, procRank,                                &
             nDimensions, processCoordinates, ierror)

        do i = 1, nDimensions
           call pigeonhole(globalSize(i), processDistribution(i),                            &
                processCoordinates(i), offset(i), localSize(i))
        end do

        success = success .and. verifyFirstDerivativeSBP(D, adjointOfD,                      &
             localSize, tolerance = 10.0_wp * epsilon(0.0_wp))

        call MPI_Barrier(MPI_COMM_WORLD, ierror)
        call MPI_Comm_free(cartesianCommunicator, ierror)

        SAFE_DEALLOCATE(isPeriodic)
        SAFE_DEALLOCATE(processCoordinates)
        SAFE_DEALLOCATE(processDistribution)
        SAFE_DEALLOCATE(globalSize)
        SAFE_DEALLOCATE(offset)
        SAFE_DEALLOCATE(localSize)

        call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL,                            &
             MPI_LAND, MPI_COMM_WORLD, ierror)
        if (.not. success) exit

     end do !... nDimensions = 1, 3

     if (.not. success) exit
     call MPI_Barrier(MPI_COMM_WORLD, ierror)
     if (MPI_Wtime() - startTime > testDuration) exit !... run for `testDuration` seconds

  end do

  call D%cleanup()
  call adjointOfD%cleanup()
  call cleanupErrorHandler()

  call MPI_Finalize(ierror)
  if (.not. success) stop -1
  stop 0

end program SBP_property

function verifyFirstDerivativeSBP(D, adjointOfD, localSize, tolerance) result(satisfiesSBP)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use StencilOperator_mod, only : t_StencilOperator

  ! <<< Internal modules >>>
  use RandomNumber, only : random

  implicit none

  ! <<< Arguments >>>
  type(t_StencilOperator), intent(in) :: D, adjointOfD
  integer, intent(in) :: localSize(:)
  real(SCALAR_KIND), intent(in), optional :: tolerance

  ! <<< Result >>>
  logical :: satisfiesSBP

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp) :: tolerance_
  integer :: i, j, k, gridIndex, is, ie, js, je, ks, ke, nDimensions, nGridPoints, gridSize(3)
  real(wp), allocatable :: f(:,:)
  SCALAR_TYPE, allocatable :: u(:,:), v(:,:)

  nDimensions = size(localSize)
  assert_key(nDimensions, (1, 2, 3))

  assert(D%direction == adjointOfD%direction)
  assert(D%direction >= 1 .and. D%direction <= nDimensions)

  tolerance_ = epsilon(0.0_wp)
  if (present(tolerance)) tolerance_ = tolerance

  gridSize(:) = 1
  gridSize(1:nDimensions) = localSize
  nGridPoints = product(localSize)

  allocate(f(nGridPoints, 1))
  allocate(u(nGridPoints, 1))
  allocate(v(nGridPoints, 1))

  call random_number(f)

  u = f
  v = f
  call D%apply(u, gridSize)
  call adjointOfD%apply(v, gridSize)
  u = D%normBoundary(1) * (u + v)

  satisfiesSBP = .true.

  is = 1
  js = 1
  ks = 1

  ie = gridSize(1)
  je = gridSize(2)
  ke = gridSize(3)

  if (D%hasDomainBoundary(1)) then

     select case (D%direction)
     case (1)
        ie = 1
     case (2)
        je = 1
     case (3)
        ke = 1
     end select

     do k = ks, ke
        do j = js, je
           do i = is, ie
              gridIndex = i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1))
              u(gridIndex, 1) = u(gridIndex, 1) + f(gridIndex, 1)
           end do
        end do
     end do

  end if

  is = 1
  js = 1
  ks = 1

  ie = gridSize(1)
  je = gridSize(2)
  ke = gridSize(3)

  if (D%hasDomainBoundary(2)) then

     select case (D%direction)
     case (1)
        is = gridSize(1)
     case (2)
        js = gridSize(2)
     case (3)
        ks = gridSize(3)
     end select

     do k = ks, ke
        do j = js, je
           do i = is, ie
              gridIndex = i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1))
              u(gridIndex, 1) = u(gridIndex, 1) - f(gridIndex, 1)
           end do
        end do
     end do

  end if

  satisfiesSBP = (maxval(abs(u)) < tolerance_)

  SAFE_DEALLOCATE(v)
  SAFE_DEALLOCATE(u)
  SAFE_DEALLOCATE(f)

end function verifyFirstDerivativeSBP
