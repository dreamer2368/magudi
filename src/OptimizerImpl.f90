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
  use State_enum, only : QOI_FORWARD_STATE, QOI_ADJOINT_STATE, QOI_RIGHT_HAND_SIDE
  use Region_enum, only : FORWARD, ADJOINT, LINEARIZED

  ! <<< Private members >>>
  use SolverImpl, only : showProgress, checkSolutionLimits, loadInitialCondition
  use RegionImpl, only : computeRegionIntegral
  use RegionVector_mod
  use RegionVectorImpl

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
  integer :: i, j, k, nDimensions, nUnknowns
  SCALAR_TYPE :: costFunctional0, costFunctional1, costSensitivity
  SCALAR_TYPE :: stepSizes(32), errorHistory(32), convergenceHistory(31)
  logical :: solutionCrashes = .false., success_

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

  do i = 1, size(region%states)
    call random_number(region%states(i)%adjointVariables)
  end do

  call region%computeRhs(FORWARD,1,1)
  costFunctional0 = functional%compute(region)
  write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))')                                      &
                              'Forward run: cost functional = ', costFunctional0
  call writeAndFlush(region%comm, output_unit, message)

  call saveRegionVector(this%base,region,QOI_FORWARD_STATE)
  ! call saveRegionVector(this%prevGrad,region,QOI_RIGHT_HAND_SIDE)
  ! call loadRegionVector(region,this%prevGrad,QOI_ADJOINT_STATE)

  call functional%updateAdjointForcing(region,.true.)
  call region%computeRhs(ADJOINT,1,1)
  call saveRegionVector(this%grad,region,QOI_RIGHT_HAND_SIDE)
  costSensitivity = regionInnerProduct(this%grad,this%grad,region)
  write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))')                                      &
                             'Adjoint run: cost sensitivity = ', costSensitivity
  call writeAndFlush(region%comm, output_unit, message)

  stepSizes(1) = 0.001_wp
  do k = 2, size(stepSizes)
     stepSizes(k) = stepSizes(k-1) * 10.0_wp**(-0.25_wp)
  end do
  errorHistory = 0.0_wp
  do k = 1, size(stepSizes)
    call loadRegionVector(region,this%base + this%grad*stepSizes(k),QOI_FORWARD_STATE)
    do i = 1, size(region%states) !... update state
       call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
            region%solverOptions)
    end do
    call region%computeRhs(FORWARD,1,1)
    costFunctional1 = functional%compute(region)

    errorHistory(k) = abs( ((costFunctional1-costFunctional0)/stepSizes(k) + costSensitivity)/costSensitivity )
    write(message, '(4(' // SCALAR_FORMAT // ',1X))') stepSizes(k), -costSensitivity,                   &
                  (costFunctional1-costFunctional0)/stepSizes(k), errorHistory(k)
    call writeAndFlush(region%comm, output_unit, message)

    if (k > 1) then
       convergenceHistory(k-1) = log(errorHistory(k) / errorHistory(k-1)) /              &
            log(stepSizes(k) / stepSizes(k-1))
       if (k > 5) then
           if (sum(convergenceHistory(k-3:k-1))/3.0_wp < 0.0_wp) exit
       end if
    end if
  end do
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
