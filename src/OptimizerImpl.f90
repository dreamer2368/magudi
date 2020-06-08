#include "config.h"

subroutine setupOptimizer(this, region)

  ! <<< Derived types >>>
  use Optimizer_mod, only : t_Optimizer
  use Region_mod, only : t_Region
  use Functional_mod, only : t_Functional

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, bufferSize
  class(t_Functional), pointer :: functional => null()

  call this%cleanup()

  this%goldenRatio = 0.5_wp * ( 1.0_wp + sqrt(5.0_wp) )
  this%initialStep = getOption("optimization/initial_",0.1_wp)
  this%linminTol = getOption("optimization/linmin_tolerance",0.1_wp)
  this%cgTol = getOption("optimization/cg_tolerance",(1.0_wp)**(-5))

  this%numGrids = size(region%grids)
  assert(region%simulationFlags%enableFunctional)
  this%numParams = 0
  if (allocated(region%params%buffer)) this%numParams = size(region%params%buffer,1)

  allocate(this%base%states(this%numGrids))
  allocate(this%grad%states(this%numGrids))
  allocate(this%conjGrad%states(this%numGrids))
  allocate(this%prevGrad%states(this%numGrids))
  allocate(this%prevCG%states(this%numGrids))
  do i = 1, this%numGrids
    allocate(this%base%states(i)%conservedVariables,MOLD=region%states(i)%conservedVariables)
    allocate(this%grad%states(i)%conservedVariables,MOLD=region%states(i)%conservedVariables)
    allocate(this%conjGrad%states(i)%conservedVariables,MOLD=region%states(i)%conservedVariables)
    allocate(this%prevGrad%states(i)%conservedVariables,MOLD=region%states(i)%conservedVariables)
    allocate(this%prevCG%states(i)%conservedVariables,MOLD=region%states(i)%conservedVariables)
  end do
  allocate(this%base%params(this%numParams))
  allocate(this%grad%params(this%numParams))
  allocate(this%conjGrad%params(this%numParams))
  allocate(this%prevGrad%params(this%numParams))
  allocate(this%prevCG%params(this%numParams))

end subroutine setupOptimizer

subroutine cleanupOptimizer(this)

  ! <<< Derived types >>>
  use Optimizer_mod, only : t_Optimizer
  use RegionVector_mod, only : t_RegionVector

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this

  call this%base%cleanup()
  call this%grad%cleanup()
  call this%conjGrad%cleanup()
  call this%prevGrad%cleanup()
  call this%prevCG%cleanup()

end subroutine cleanupOptimizer

subroutine verifyNLCG(this, region, solver)

  ! <<< External modules >>>
  use iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver
  use Optimizer_mod, only : t_Optimizer
  use Functional_mod, only : t_Functional
  use ActuatorPatch_mod, only : t_ActuatorPatch

  ! <<< Enumerations >>>
  use State_enum, only : QOI_FORWARD_STATE, QOI_ADJOINT_STATE, QOI_TIME_AVERAGED_STATE
  use Region_enum, only : FORWARD, LINEARIZED

  ! <<< Private members >>>
  use SolverImpl, only : showProgress, checkSolutionLimits, loadInitialCondition
  use RegionImpl, only : computeRegionIntegral

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming
  use ErrorHandler, only : writeAndFlush
  use InputHelper, only : getOption, getRequiredOption

  ! <<< SeungWhan: debug >>>
  use InputHelper, only : getOption
  use, intrinsic :: iso_fortran_env, only : output_unit
  use ErrorHandler, only : writeAndFlush

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this
  class(t_Solver) :: solver
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: filename, message
  class(t_Functional), pointer :: functional => null()
  integer :: i, j
  SCALAR_TYPE :: instantaneousCostFunctional, costFunctional
  logical :: solutionCrashes = .false.

  call startTiming("verifyNLCG")

  ! Connect to the previously allocated functional.
  assert(region%simulationFlags%enableFunctional)
  call solver%functionalFactory%connect(functional)
  assert(associated(functional))
  functional%runningTimeQuadrature = 0.0_wp

  ! Load the initial condition.
  call loadInitialCondition(solver, region, FORWARD)                              ! for baseline simulation.

  if (region%simulationFlags%enableBodyForce .or. region%simulationFlags%checkConservation) then
    region%oneOverVolume = computeRegionIntegral(region)
    region%oneOverVolume = 1.0_wp / region%oneOverVolume

    region%momentumLossPerVolume = 0.0_wp
  end if

  ! ! Save the initial condition if it was not specified as a restart file.
  ! write(filename, '(2A,I8.8,A)') trim(solver%outputPrefix), "-", startTimestep, ".q"
  ! call region%saveData(QOI_FORWARD_STATE, filename)

  ! Reset probes.
  if (solver%probeInterval > 0) call region%resetProbes()

  do i = 1, size(region%states) !... update state
     call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
          region%solverOptions)
  end do

  call region%computeRhs(FORWARD,1,1)

  !TODO: Update the delta functional.
  instantaneousCostFunctional = functional%compute(region)
  functional%runningTimeQuadrature = instantaneousCostFunctional

  ! ! Report simulation progess.
  ! call showProgress(solver, region, LINEARIZED, startTimestep, timestep,                       &
  !     time, instantaneousCostFunctional)

  ! Save solution on probe patches.
  if (solver%probeInterval > 0) call region%saveProbeData(FORWARD)
  !
  ! ! Filter solution if required.
  ! if (region%simulationFlags%filterOn) then
  !   do j = 1, size(region%grids)
  !      call region%grids(j)%applyFilter(region%states(j)%conservedVariables, timestep)
  !   end do
  ! end if

  ! Finish writing remaining data gathered on probes.
  if (solver%probeInterval > 0) call region%saveProbeData(FORWARD, finish = .true.)

  costFunctional = functional%runningTimeQuadrature
  write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))') 'Forward run: cost functional = ', &
                                                           costFunctional
  call writeAndFlush(region%comm, output_unit, message)

  call endTiming("verifyNLCG")

end subroutine verifyNLCG

subroutine mnbrakConjugateGradient(this, solver, region)

  ! <<< Derived types >>>
  use Optimizer_mod, only : t_Optimizer
  use Solver_mod, only : t_Solver
  use Region_mod, only : t_Region

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this
  class(t_Solver) :: solver
  class(t_Region) :: region

  ! <<< Local variables >>>

end subroutine mnbrakConjugateGradient
