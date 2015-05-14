#include "config.h"

program composite_dissipation_sanity

  use MPI
  use, intrinsic :: iso_fortran_env, only : real64

  use MPIHelper, only : pigeonhole
  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random
  use StencilOperator_mod, only : t_StencilOperator

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  real(real64), parameter :: testDuration = 2.0_real64
  character(len = STRING_LENGTH), parameter :: discretizationTypes(2) =                      &
       (/ "SBP 2-4", "SBP 4-8" /)
  real(wp), parameter :: compositeAmountFactors(2) = (/ 16.0_wp, 256.0_wp /)
  integer :: i, n, direction, nDimensions, cartesianCommunicator, numProcs, procRank, ierror
  integer, dimension(:), allocatable :: localSize, offset, globalSize,                       &
       processDistribution, processCoordinates
  logical :: success
  logical, allocatable :: isPeriodic(:)
  real(real64) :: startTime
  character(len = STRING_LENGTH) :: str
  type(t_StencilOperator) :: compositeDissipation, dissipation, dissipationTranspose
  real(wp) :: compositeAmountFactor

  interface

     function checkCompositeDissipationOperator(compositeDissipation, dissipation,           &
          dissipationTranspose, localSize, compositeAmountFactor, tolerance) result(success)

       use StencilOperator_mod, only : t_StencilOperator

       type(t_StencilOperator), intent(in) :: compositeDissipation,                          &
            dissipation, dissipationTranspose
       integer, intent(in) :: localSize(:)
       real(SCALAR_KIND), intent(in) :: compositeAmountFactor

       real(SCALAR_KIND), intent(in), optional :: tolerance

       logical :: success

     end function checkCompositeDissipationOperator

  end interface

  call MPI_Init(ierror)

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)

  startTime = MPI_Wtime()

  do

     i = random(1, size(discretizationTypes))

     str = trim(discretizationTypes(i))
     call MPI_Bcast(str, len(str), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)

     compositeAmountFactor = compositeAmountFactors(i)
     call MPI_Bcast(compositeAmountFactor, 1, REAL_TYPE_MPI, 0, MPI_COMM_WORLD, ierror)

     call dissipation%setup(trim(str) // " dissipation")
     call dissipationTranspose%setup(trim(str) // " dissipation transpose")
     call compositeDissipation%setup(trim(str) // " composite dissipation")

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
        call dissipation%update(cartesianCommunicator, direction)
        call dissipationTranspose%update(cartesianCommunicator, direction)
        call compositeDissipation%update(cartesianCommunicator, direction)

        do i = 1, nDimensions
           if (i == direction) then
              n = numProcs * (2 * compositeDissipation%boundaryWidth + 1)
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

        success = success .and. checkCompositeDissipationOperator(compositeDissipation,      &
             dissipation, dissipationTranspose, localSize, compositeAmountFactor,            &
             tolerance = compositeAmountFactor * epsilon(0.0_wp))

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

  call dissipation%cleanup()
  call dissipationTranspose%cleanup()
  call compositeDissipation%cleanup()
  call cleanupErrorHandler()

  call MPI_Finalize(ierror)
  if (.not. success) stop -1
  stop 0

end program composite_dissipation_sanity

function checkCompositeDissipationOperator(compositeDissipation, dissipation,                &
     dissipationTranspose, localSize, compositeAmountFactor, tolerance) result(success)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use StencilOperator_mod, only : t_StencilOperator

  ! <<< Internal modules >>>
  use RandomNumber, only : random

  implicit none

  ! <<< Arguments >>>
  type(t_StencilOperator), intent(in) :: compositeDissipation,                               &
       dissipation, dissipationTranspose
  integer, intent(in) :: localSize(:)
  real(SCALAR_KIND), intent(in) :: compositeAmountFactor
  real(SCALAR_KIND), intent(in), optional :: tolerance

  ! <<< Result >>>
  logical :: success

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp) :: tolerance_
  integer :: nDimensions, nGridPoints, gridSize(3)
  real(wp), allocatable :: f(:,:)
  SCALAR_TYPE, allocatable :: u(:,:), v(:,:)

  nDimensions = size(localSize)
  assert_key(nDimensions, (1, 2, 3))

  assert(dissipation%direction == dissipationTranspose%direction)
  assert(dissipation%direction == compositeDissipation%direction)
  assert(dissipation%direction >= 1 .and. dissipation%direction <= nDimensions)

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
  call compositeDissipation%apply(u, gridSize)
  u = u * compositeAmountFactor

  call dissipation%apply(v, gridSize)
  call dissipationTranspose%apply(v, gridSize)
  call compositeDissipation%applyNormInverse(v, gridSize)

  success = (maxval(abs(u + v)) / real(compositeDissipation%boundaryWidth, wp) < tolerance_)

  SAFE_DEALLOCATE(v)
  SAFE_DEALLOCATE(u)
  SAFE_DEALLOCATE(f)

end function checkCompositeDissipationOperator
