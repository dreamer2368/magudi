#include "config.h"

program time_integrator

  use MPI

  use CNSHelper
  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random

  implicit none

  logical :: success, success_, isPeriodic
  integer :: i, j, k, nDimensions, ierror
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

       logical, intent(in) :: isPeriodic
       real(SCALAR_KIND), intent(in), optional :: tolerance

     end subroutine testAdjointRelation

  end interface

  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  do nDimensions = 1, 1
    do j = 1, 1 !... for each discretizationTypes
      success = .true.
      do i = 1, 1 !... test multiple times
        isPeriodic = .false.
        call testAdjointRelation(discretizationTypes(j), nDimensions,           &
                                 success_, isPeriodic)
        success = success .and. success_
      end do
      if( procRank == 0 .and. success ) then
        print *, 'Success'
        print *, 'dimension: ', nDimensions
      end if
    end do
  end do

  call cleanupErrorHandler()

  call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierror)
  call MPI_Finalize(ierror)
  if (.not. success) stop -1
  stop 0

end program time_integrator

subroutine testAdjointRelation(identifier, nDimensions, success, isPeriodic, tolerance)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use StencilOperator_mod, only : t_StencilOperator
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use testRegion_mod, only : t_testRegion
  use TimeIntegrator_mod, only : t_TimeIntegrator
  use TimeIntegrator_factory, only : t_TimeIntegratorFactory

  use Region_enum, only : FORWARD, ADJOINT

  use CNSHelper

  ! <<< Internal modules >>>
  use MPIHelper, only : pigeonhole
  use RandomNumber, only : random

  ! <<< Arguments >>>
  character(len = *), intent(in) :: identifier
  integer, intent(in) :: nDimensions
  logical, intent(out) :: success
  logical, intent(in) :: isPeriodic
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
  type(t_testRegion) :: region
  type(t_TimeIntegratorFactory) :: timeIntegratorFactory
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  logical :: isPeriodic_(3), hasNegativeJacobian
  real(wp) :: scalar1, scalar2, tolerance_,&
              stepSizes(32), errorHistory(32), convergenceHistory(31)
  integer, allocatable :: gridSize(:,:)
  integer :: i, j, k, nUnknowns, nTimesteps, timestep, startTimestep, timemarchDirection
  real(wp) :: timeStepSize, time, error
  real(SCALAR_KIND), allocatable :: forwardState(:), adjointState(:),           &
                                    forwardRhs(:), adjointRhs(:)
  character(len = STRING_LENGTH) :: errorMessage

  tolerance_ = 100*epsilon(0.0_wp)
  if( present(tolerance) ) tolerance_ = tolerance

  success = .true.

  ! set up simulation flags
  call simulationFlags%initialize()
  simulationFlags%enableController = .true.
  simulationFlags%enableFunctional = .true.
  simulationFlags%enableAdjoint = .true.
  ! randomize curvilinear domain
  ! simulationFlags%isDomainCurvilinear = (random(0, 2) == 0)
  simulationFlags%isDomainCurvilinear = .false.
  simulationFlags%viscosityOn = .false.
  simulationFlags%useConstantCfl = .false.

  ! randomize grid size
  ! Note that too small grid size will yield null matrix for stencil operators.
  allocate(gridSize(nDimensions,1))
  gridSize = 1
  do i = 1, nDimensions
     gridSize(i,1) = 1
  end do
  ! randomize timesteps
  nTimesteps = random(1,4800)

  ! initialize solver option and region
  simulationFlags%useConstantCfl = .true.
  call solverOptions%initialize(nDimensions, simulationFlags)
  simulationFlags%useConstantCfl = .false.
  solverOptions%discretizationType = trim(identifier)
  solverOptions%timeStepSize = random(0.01_wp, 1.0_wp)
  call region%setup(MPI_COMM_WORLD,gridSize,simulationFlags,solverOptions)

  ! initialize time integrator
  call timeIntegratorFactory%connect(timeIntegrator,                                    &
       trim(region%solverOptions%timeIntegratorType))
  assert(associated(timeIntegrator))
  call timeIntegrator%setup(region)

  ! allocate time-dependeant variables
  allocate(forwardState(nTimesteps*timeIntegrator%nStages+1))
  allocate(adjointState(nTimesteps*timeIntegrator%nStages+1))
  allocate(forwardRhs(nTimesteps*timeIntegrator%nStages))
  allocate(adjointRhs(nTimesteps*timeIntegrator%nStages))
  allocate(region%tempRhs(nTimesteps*timeIntegrator%nStages))

  ! randomize RHS
  call random_number(forwardRhs)
  call random_number(adjointRhs)

  ! Randomize initial condition: uniform over space, unknowns
  region%states(1)%conservedVariables = random(-10.0_wp, 10.0_wp)
  forwardState(1) = region%states(1)%conservedVariables(1,1)

  ! Forward run: Dx = f
  region%tempRhs = forwardRhs
  region%tempRhsIndex = 1
  ! scalar1 = (D^{\dagger}x^{\dagger}) dot x
  scalar1 = 0.0_wp

  time = 0.0_wp
  region%states(:)%time = time
  do timestep = 1, nTimesteps

    region%timestep = timestep
    timeStepSize = region%getTimeStepSize()

    do i = 1, timeIntegrator%nStages
      ! Take a single sub-step using the time integrator.
      call timeIntegrator%substepForward(region, time, timeStepSize, timestep, i)

      ! Save forward states
      forwardState((timestep-1)*timeIntegrator%nStages+i+1) = region%states(1)%conservedVariables(1,1)

      ! Update the cost functional.
       scalar1 = scalar1 + adjointRhs((timestep-1)*timeIntegrator%nStages+i)                  &
            * timeIntegrator%norm(i) * timeStepSize * region%states(1)%conservedVariables(1,1)
    end do

  end do

  ! Adjoint run = D^{\dagger}x^{\dagger} = g
  region%tempRhs = adjointRhs
  region%tempRhsIndex = size(adjointRhs)
  ! scalar2 = x^{\dagger} dot (D x)
  scalar2 = 0.0_wp

  region%states(1)%adjointVariables = -region%solverOptions%timeStepSize/6.0_wp                &
                                          * adjointRhs(size(adjointRhs))
  adjointState(size(adjointState)) = region%states(1)%adjointVariables(1,1)

  time = nTimesteps * region%solverOptions%timeStepSize
  region%states(:)%time = time
  startTimestep = nTimesteps
  timemarchDirection = -1
  do timestep = startTimestep + sign(1, timemarchDirection),&
      startTimestep + sign(nTimesteps, timemarchDirection), timemarchDirection

    region%timestep = timestep
    timeStepSize = region%getTimeStepSize()

    do i = timeIntegrator%nStages, 1, -1
      ! Update the cost functional.
       scalar2 = scalar2 + forwardRhs(timestep * timeIntegrator%nStages + i)                  &
            * timeIntegrator%norm(i) * timeStepSize * region%states(1)%adjointVariables(1,1)

      ! Take a single sub-step using the time integrator.
      call timeIntegrator%substepAdjoint(region, time, timeStepSize, timestep, i)

      ! Save forward states
      adjointState(timestep*timeIntegrator%nStages+i) = region%states(1)%adjointVariables(1,1)
    end do

  end do

  error = abs( (scalar1 + (scalar2 + adjointState(1)*forwardState(1)))/scalar1 )
  success = (error < tolerance_)

  SAFE_DEALLOCATE(forwardState)
  SAFE_DEALLOCATE(forwardRhs)
  SAFE_DEALLOCATE(adjointState)
  SAFE_DEALLOCATE(adjointRhs)
  SAFE_DEALLOCATE(region%tempRhs)

  call timeIntegratorFactory%cleanup()
  call region%cleanup()

  SAFE_DEALLOCATE(gridSize)

end subroutine testAdjointRelation

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
