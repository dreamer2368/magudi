#include "config.h"

program SBP_operators_grid

  use MPI

  use CNSHelper
  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random

  implicit none

  logical :: success, success_, isPeriodic
  integer :: i, j, nDimensions, ierror
  integer :: procRank
  character(len = STRING_LENGTH), parameter :: discretizationTypes(4) =                      &
       (/ "SBP 1-2", "SBP 2-4", "SBP 3-6", "SBP 4-8" /)

  interface

     real(SCALAR_KIND) function meanTrimmed(a)
       real(SCALAR_KIND), intent(inout) :: a(:)
     end function meanTrimmed

     subroutine sort(a)
       real(SCALAR_KIND), intent(inout) :: a(:)
     end subroutine sort

     subroutine testAdjointRelation(identifier, nDimensions, success, isPeriodic, tolerance)

       character(len = *), intent(in) :: identifier
       integer, intent(in) :: nDimensions
       logical, intent(out) :: success

       logical, intent(in), optional :: isPeriodic
       real(SCALAR_KIND), intent(in), optional :: tolerance

     end subroutine testAdjointRelation

  end interface

  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  do nDimensions = 1, 3
    do j = 1, 4 !... for each discretizationTypes
      success = .true.
      do i = 1, 10 !... test multiple times
        ! Didn't test periodic grid yet!!
        ! isPeriodic = .true.
        ! call testAdjointRelation(discretizationTypes(j), nDimensions,           &
        !                          success_, isPeriodic)
        ! success = success .and. success_
        ! if( .not. success) then
        !   if( procRank == 0 ) then
        !     print *, 'Failed, ', trim(discretizationTypes(j))
        !     print *, 'dimension: ', nDimensions
        !     print *, 'periodicity: ', isPeriodic
        !   end if
        !   exit
        ! end if

        isPeriodic = .false.
        call testAdjointRelation(discretizationTypes(j), nDimensions,           &
                                 success_, isPeriodic)
        success = success .and. success_
        if( .not. success_) then
          if( procRank == 0 ) then
            print *, 'Failed, ', trim(discretizationTypes(j))
            print *, 'dimension: ', nDimensions
            print *, 'periodicity: ', isPeriodic
          end if
          exit
        end if
      end do
      if( procRank == 0 .and. success ) then
        print *, 'Success, ', trim(discretizationTypes(j))
        print *, 'dimension: ', nDimensions
      end if
    end do
  end do

  call cleanupErrorHandler()

  call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierror)
  call MPI_Finalize(ierror)
  if (.not. success) stop -1
  stop 0

end program SBP_operators_grid

subroutine sort(a)

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(inout) :: a(:)

  ! <<< Local variables >>>
  integer :: i, j
  real(SCALAR_KIND) :: temp

  do i = 2, size(a)
     j = i - 1
     temp = a(i)
     do while (a(j) > temp)
        a(j+1) = a(j)
        j = j - 1
        if (j < 1) exit
     end do
     a(j+1) = temp
  end do

end subroutine sort

real(SCALAR_KIND) function median(a)

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(in) :: a(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: n

  n = size(a)
  if (mod(n, 2) == 1) then
     median = a((n + 1) / 2)
  else
     median = 0.5_wp * (a(n / 2) + a(n / 2 + 1))
  end if

end function median

real(SCALAR_KIND) function meanTrimmed(a)

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(inout) :: a(:)

  ! <<< Scalar variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, n
  real(wp) :: firstQuartile, thirdQuartile

  interface
     real(SCALAR_KIND) function median(a)
       real(SCALAR_KIND), intent(in) :: a(:)
     end function median
  end interface

  n = size(a)

  if (mod(n, 2) == 0) then
     firstQuartile = median(a(1:n/2))
     thirdQuartile = median(a(n/2+1:n))
  else
     firstQuartile = median(a(1:(n-1)/2))
     thirdQuartile = median(a((n+1)/2+1:n))
  end if

  meanTrimmed = 0.0_wp
  n = 0

  do i = 1, size(a)
     if (a(i) >= firstQuartile .and. a(i) <= thirdQuartile) then
        meanTrimmed = meanTrimmed + a(i)
        n = n + 1
     end if
  end do

  if (n == 0) return
  meanTrimmed = meanTrimmed / real(n, wp)

end function meanTrimmed

subroutine testAdjointRelation(identifier, nDimensions, success, isPeriodic, tolerance)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use StencilOperator_mod, only : t_StencilOperator
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  use Region_enum, only : FORWARD, ADJOINT

  use CNSHelper

  ! <<< Internal modules >>>
  use MPIHelper, only : pigeonhole
  use RandomNumber, only : random

  ! <<< Arguments >>>
  character(len = *), intent(in) :: identifier
  integer, intent(in) :: nDimensions
  logical, intent(out) :: success
  logical, intent(in), optional :: isPeriodic
  real(SCALAR_KIND), intent(in), optional :: tolerance

  ! <<< interface >>>
  interface
     real(SCALAR_KIND) function meanTrimmed(a)
       real(SCALAR_KIND), intent(inout) :: a(:)
     end function meanTrimmed

     subroutine sort(a)
       real(SCALAR_KIND), intent(inout) :: a(:)
     end subroutine sort
  end interface

  ! <<< Local derived type variables >>>
  type(t_SimulationFlags) :: simulationFlags
  type(t_SolverOptions) :: solverOptions
  type(t_Grid) :: grid

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp) :: scalar1, scalar2, scalar3, scalar4,                                  &
              stepSizes(32), errorHistory(32), convergenceHistory(31)
  integer :: i, j, k, gridSize(nDimensions, 1), nUnknowns
  real(SCALAR_KIND), allocatable :: u(:,:), v(:,:), Du(:,:), DTv(:,:), norm(:,:)
  SCALAR_TYPE, dimension(nDimensions) :: h, gridPerturbation

  success = .true.

  ! set up simulation flags
  call simulationFlags%initialize()
  simulationFlags%enableController = .true.
  simulationFlags%enableFunctional = .true.
  simulationFlags%enableAdjoint = .true.
  simulationFlags%isDomainCurvilinear = .true.
  simulationFlags%viscosityOn = .true.
  simulationFlags%repeatFirstDerivative = .true. ! this is default value.

  ! randomize grid size
  ! Note that too small grid size will yield null matrix for stencil operators.
  gridSize(:,:) = 1
  do i = 1, nDimensions
     gridSize(i,1) = random(20, 40)
  end do

  ! initialize solver option, grid, and states
  call solverOptions%initialize(nDimensions, simulationFlags)
  solverOptions%discretizationType = trim(identifier)
  solverOptions%reynoldsNumberInverse = 1.0_wp / 5.0_wp
  solverOptions%prandtlNumberInverse = 1.0_wp / 0.72_wp
  solverOptions%powerLawExponent = 0.666_wp
  solverOptions%bulkViscosityRatio = 0.6_wp
  call grid%setup(1, gridSize(1:nDimensions,1), MPI_COMM_WORLD,                        &
       simulationFlags = simulationFlags)
  call grid%setupSpatialDiscretization(simulationFlags, solverOptions)

  ! randomize grid coordinates (didn't reflect periodicity)
  h = 1.0_wp / real(grid%globalSize(1:nDimensions)-1,wp)
  do i = 1, grid%nGridPoints
    call random_number(gridPerturbation)
    gridPerturbation = (2.0_wp * gridPerturbation - 1.0_wp) * 0.05_wp * h
    grid%coordinates(i,:) = grid%coordinates(i,:) + gridPerturbation
  end do

  ! grid is not randomized: do not put any argument in updateGrid!!
  call grid%update()
  ! print *, 'Jacobian range: (', minval(grid%jacobian), ', ', maxval(grid%jacobian), ')'
  ! print *, 'Metrics range: (', minval(grid%metrics), ', ', maxval(grid%metrics), ')'

  allocate(norm(grid%nGridPoints,1))

  allocate(u(grid%nGridPoints,1))
  allocate(v(grid%nGridPoints,1))
  allocate(Du(grid%nGridPoints,1))
  allocate(DTv(grid%nGridPoints,1))
  call random_number(u)
  call random_number(v)
  do i = 1, nDimensions
    Du = u
    DTv = v
    call grid%firstDerivative(i)%apply(Du, grid%localSize)
    call grid%adjointFirstDerivative(i)%apply(DTv, grid%localSize)

    norm = 1.0_wp
    call grid%firstDerivative(i)%applyNorm(norm, grid%localSize)

    print *, 'grid%norm: ', grid%norm
    print *, 'norm: ', norm

    scalar3 = sum(Du * norm * v)
    scalar4 = sum(u * norm * DTv)
    if (abs(scalar3 - scalar4) > 1.0e-10) then
      print *, 'scalar3 =/= scalar4, ', scalar3, scalar4
    end if
    assert(abs(scalar3 - scalar4) < 1.0e-10)
    
    scalar1 = grid%computeInnerProduct(Du, v)
    scalar2 = grid%computeInnerProduct(u, DTv)

    if (abs(scalar1 - scalar2) > 1.0e-10) then
      print *, 'scalar1 =/= scalar2, ', scalar1, scalar2
    end if
    assert(abs(scalar1 - scalar2) < 1.0e-10)
  end do

  call grid%cleanup()

end subroutine testAdjointRelation
