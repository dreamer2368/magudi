#include "config.h"

module SolverImpl

  implicit none
  public

contains

  subroutine showProgress(this, region, mode, startTimestep,                                 &
       timestep, time, instantaneousFunctional)

    ! <<< External modules >>>
    use iso_fortran_env, only : output_unit

    ! <<< Derived types >>>
    use Region_mod, only : t_Region
    use Solver_mod, only : t_Solver
    use Controller_mod, only : t_Controller
    use Functional_mod, only : t_Functional

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE, QOI_ADJOINT_STATE
    use Region_enum, only : FORWARD, ADJOINT

    ! <<< Internal modules >>>
    use ErrorHandler, only : writeAndFlush

    ! <<< Arguments >>>
    class(t_Solver) :: this
    class(t_Region) :: region
    integer, intent(in) :: mode, timestep, startTimestep
    real(SCALAR_KIND), intent(in) :: time, instantaneousFunctional

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: nDimensions
    real(SCALAR_KIND) :: timeStepSize, cfl, residuals(3)
    character(len = STRING_LENGTH) :: str, str_, filename
    class(t_Controller), pointer :: controller => null()
    class(t_Functional), pointer :: functional => null()

    assert_key(mode, (FORWARD, ADJOINT))
    assert(startTimestep >= 0)

    nDimensions = size(region%globalGridSizes, 1)
    assert_key(nDimensions, (1, 2, 3))

    ! Report time step size in constant CFL mode
    if (region%simulationFlags%useConstantCfl) then
       timeStepSize = region%getTimeStepSize()
    else
       cfl = region%getCfl()
    end if

    if (this%reportInterval > 0 .and. mod(timestep, max(1, this%reportInterval)) == 0) then

       if (region%simulationFlags%useConstantCfl) then
          write(str, '(2A,I8,2(A,E13.6))') PROJECT_NAME, ": timestep = ", timestep,          &
               ", dt = ", timeStepSize, ", time = ", abs(time)
       else
          write(str, '(2A,I8,2(A,E13.6))') PROJECT_NAME, ": timestep = ", timestep,          &
               ", CFL = ", cfl, ", time = ", abs(time)
       end if

       if (.not. region%simulationFlags%predictionOnly) then

          select case (mode)
          case (FORWARD)
             write(str_, '(A,E13.6)') ", cost = ", instantaneousFunctional
          case (ADJOINT)
             write(str_, '(A,E13.6)') ", gradient = ", instantaneousFunctional
          end select

          str = trim(str) // trim(str_)

       end if

       call writeAndFlush(region%comm, output_unit, str)

       if (.not. region%simulationFlags%predictionOnly .and. region%outputOn) then

          select case (mode)

          case (FORWARD)
             call this%functionalFactory%connect(functional)
             assert(associated(functional))
             call functional%writeToFile(region%comm, trim(this%outputPrefix) //             &
                  ".cost_functional.txt", timestep, time,                                    &
                  timestep - startTimestep > this%reportInterval)

          case (ADJOINT)
             call this%controllerFactory%connect(controller)
             assert(associated(controller))
             call controller%writeSensitivityToFile(region%comm, trim(this%outputPrefix) //  &
                  ".cost_sensitivity.txt", timestep, time,                                   &
                  startTimestep - timestep > this%reportInterval)

          end select

       end if

    end if

    if (this%saveInterval > 0 .and. mod(timestep, max(1, this%saveInterval)) == 0) then

       select case (mode)
       case (FORWARD)
          write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", timestep, ".q"
          call region%saveData(QOI_FORWARD_STATE, filename)
       case (ADJOINT)
          write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", timestep, ".adjoint.q"
          call region%saveData(QOI_ADJOINT_STATE, filename)
       end select

    end if

    if (region%simulationFlags%steadyStateSimulation .and.                                   &
         this%residualManager%reportInterval > 0 .and.                                       &
         mod(timestep, max(1, this%residualManager%reportInterval)) == 0) then

       call this%residualManager%compute(region)

       select case (mode)

       case (FORWARD)
          call this%residualManager%writeToFile(region%comm, trim(this%outputPrefix) //      &
               ".residuals.txt", timestep, time,                                             &
               timestep - startTimestep > this%residualManager%reportInterval)
       case (ADJOINT)
          call this%residualManager%writeToFile(region%comm, trim(this%outputPrefix) //      &
               ".adjoint_residuals.txt", timestep, time,                                     &
               timestep - startTimestep > this%residualManager%reportInterval)
       end select

       residuals(1) = this%residualManager%residuals(1)
       residuals(2) = maxval(this%residualManager%residuals(1:nDimensions))
       residuals(3) = this%residualManager%residuals(nDimensions+2)

       write(str, '(2X,3(A,(ES11.4E2)))') "residuals: density = ", residuals(1),             &
            ", momentum = ", residuals(2), ", energy = ", residuals(3)
       call writeAndFlush(region%comm, output_unit, str)

       if (this%residualManager%hasSimulationConverged) then

          call writeAndFlush(region%comm, output_unit, "Solution has converged!")

          select case (mode)
          case (FORWARD)
             call region%saveData(QOI_FORWARD_STATE,                                         &
                  trim(this%outputPrefix) // ".steady_state.q")
          case (ADJOINT)
             call region%saveData(QOI_ADJOINT_STATE,                                         &
                  trim(this%outputPrefix) // ".steady_state.adjoint.q")
          end select

       end if

    end if

  end subroutine showProgress

  subroutine checkSolutionLimits(region, mode, outputPrefix)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use State_mod, only : t_State
    use Region_mod, only : t_Region

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE
    use Region_enum, only : FORWARD

    ! <<< Internal modules >>>
    use ErrorHandler, only : gracefulExit

    ! <<< Arguments >>>
    class(t_Region) :: region
    integer, intent(in) :: mode
    character(len = *), intent(in) :: outputPrefix

    ! <<< Local variables >>>
    integer :: i, iGlobal, jGlobal, kGlobal, rankReportingError, procRank, ierror
    character(len = STRING_LENGTH) :: message
    SCALAR_TYPE :: fOutsideRange

    rankReportingError = -1
    call MPI_Comm_rank(region%comm, procRank, ierror)

    do i = 1, size(region%states)

       if (.not. region%grids(i)%isVariableWithinRange(                                      &
            region%states(i)%conservedVariables(:,1),                                        &
            fOutsideRange, iGlobal, jGlobal, kGlobal,                                        &
            minValue = region%solverOptions%densityRange(1),                                 &
            maxValue = region%solverOptions%densityRange(2))) then
          write(message, '(4(A,I0.0),3(A,(SS,ES9.2E2)),A)') "Density on grid ",              &
               region%grids(i)%index, " at (", iGlobal, ", ", jGlobal, ", ", kGlobal, "): ", &
               fOutsideRange, " out of range (",                                             &
               region%solverOptions%densityRange(1), ", ",                                   &
               region%solverOptions%densityRange(2), ")!"
          rankReportingError = procRank
          exit
       end if

       if (.not. region%grids(i)%isVariableWithinRange(region%states(i)%temperature(:,1),    &
            fOutsideRange, iGlobal, jGlobal, kGlobal,                                        &
            minValue = region%solverOptions%temperatureRange(1),                             &
            maxValue = region%solverOptions%temperatureRange(2))) then
          write(message, '(4(A,I0.0),3(A,(SS,ES9.2E2)),A)') "Temperature on grid ",          &
               region%grids(i)%index, " at (", iGlobal, ", ", jGlobal, ", ", kGlobal, "): ", &
               fOutsideRange, " out of range (",                                             &
               region%solverOptions%temperatureRange(1), ", ",                               &
               region%solverOptions%temperatureRange(2), ")!"
          rankReportingError = procRank
          exit
       end if

    end do

    call MPI_Allreduce(MPI_IN_PLACE, rankReportingError, 1,                                  &
         MPI_INTEGER, MPI_MAX, region%comm, ierror)

    if (rankReportingError /= -1) then

       if (procRank == 0 .and. rankReportingError /= 0)                                      &
            call MPI_Recv(message, STRING_LENGTH, MPI_CHARACTER, rankReportingError,         &
            rankReportingError, region%comm, MPI_STATUS_IGNORE, ierror)
       if (procRank == rankReportingError .and. rankReportingError /= 0)                     &
            call MPI_Send(message, STRING_LENGTH, MPI_CHARACTER, 0, procRank,                &
            region%comm, ierror)

       select case (mode)
       case (FORWARD)
          call region%saveData(QOI_FORWARD_STATE, trim(outputPrefix) // "-crashed.q")
       end select

       call gracefulExit(region%comm, message)

    end if

  end subroutine checkSolutionLimits

  subroutine loadInitialCondition(this, region, mode, restartFilename)

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use Region_mod, only : t_Region
    use Solver_mod, only : t_Solver
    use Functional_mod, only : t_Functional
    use CostTargetPatch_mod, only : t_CostTargetPatch

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE, QOI_ADJOINT_STATE
    use Region_enum, only : FORWARD, ADJOINT

    ! <<< Internal modules >>>
    use InputHelper, only : getOption, getRequiredOption

    implicit none

    ! <<< Arguments >>>
    class(t_Solver) :: this
    class(t_Region) :: region
    integer, intent(in) :: mode
    character(len = *), intent(in), optional :: restartFilename

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    character(len = STRING_LENGTH) :: filename
    integer :: i, j
    real(wp) :: timeStepSize
    class(t_Functional), pointer :: functional => null()
    class(t_Patch), pointer :: patch => null()

    select case (mode)

    case (FORWARD) !... initialize conserved variables.

       if (present(restartFilename) .and. region%simulationFlags%predictionOnly) then
          call region%loadData(QOI_FORWARD_STATE, restartFilename)
       else if (region%simulationFlags%useTargetState) then
          filename = getOption("initial_condition_file", "")
          if (len_trim(filename) == 0) then
             region%timestep = 0
             do i = 1, size(region%states) !... initialize from target state.
                region%states(i)%conservedVariables = region%states(i)%targetState
                region%states(i)%time = 0.0_wp
             end do
          else
             call region%loadData(QOI_FORWARD_STATE, filename) !... initialize from file.
          end if
       else
          call getRequiredOption("initial_condition_file", filename, region%comm)
          call region%loadData(QOI_FORWARD_STATE, filename) !... initialize from file.
       end if

    case (ADJOINT) !... initialize adjoint variables.

       if (.not. region%simulationFlags%predictionOnly) then
          filename = getOption("adjoint_initial_condition_file", "")
          if (len_trim(filename) == 0) then
             do i = 1, size(region%states)
                region%states(i)%adjointVariables = 0.0_wp
             end do
          else
             call region%loadData(QOI_ADJOINT_STATE, filename)
          end if
       end if

       if (allocated(region%patchFactories) .and.                                            &
            .not. region%simulationFlags%steadyStateSimulation .and.                         &
            .not. region%simulationFlags%useContinuousAdjoint) then

          timeStepSize = region%getTimeStepSize()
          region%states(:)%adjointForcingFactor = - timeStepSize / 6.0_wp !... RK4 only.

          ! Connect to the previously allocated functional.
          call this%functionalFactory%connect(functional)
          assert(associated(functional))
          call functional%updateAdjointForcing(region)

          do i = 1, size(region%states)
             region%states(i)%rightHandSide = 0.0_wp
          end do

          do i = 1, size(region%patchFactories)
             call region%patchFactories(i)%connect(patch)
             if (.not. associated(patch)) cycle
             select type (patch)
                class is (t_CostTargetPatch)
                do j = 1, size(region%states)
                   if (patch%gridIndex /= region%grids(j)%index) cycle
                   call patch%updateRhs(ADJOINT, region%simulationFlags,                     &
                        region%solverOptions, region%grids(j), region%states(j))
                end do
             end select
          end do

          do i = 1, size(region%states)
             region%states(i)%adjointVariables = region%states(i)%rightHandSide
          end do
          region%states(:)%adjointForcingFactor = 1.0_wp !... restore

       end if

    end select

  end subroutine loadInitialCondition

end module SolverImpl

subroutine setupSolver(this, region, restartFilename, outputPrefix)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver
  use Controller_mod, only : t_Controller
  use Functional_mod, only : t_Functional
  use TimeIntegrator_mod, only : t_TimeIntegrator

  ! <<< Enumerations >>>
  use Grid_enum, only : QOI_CONTROL_MOLLIFIER, QOI_TARGET_MOLLIFIER
  use State_enum, only : QOI_FORWARD_STATE, QOI_TARGET_STATE, QOI_ADJOINT_STATE

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use Patch_factory, only : computeSpongeStrengths, updatePatchFactories
  use InterfaceHelper, only : checkFunctionContinuityAtInterfaces

  implicit none

  ! <<< Arguments >>>
  class(t_Solver) :: this
  class(t_Region) :: region
  character(len = *), intent(in), optional :: restartFilename, outputPrefix

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: filename
  integer :: i
  class(t_Controller), pointer :: controller => null()
  class(t_Functional), pointer :: functional => null()
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()

  if (present(outputPrefix)) then
     this%outputPrefix = outputPrefix
  else
     this%outputPrefix = getOption("output_prefix", PROJECT_NAME)
  end if

  this%nTimesteps = getOption("number_of_timesteps", 1000)
  this%nTimesteps = max(0, this%nTimesteps)

  this%saveInterval = getOption("save_interval", -1)
  if (this%saveInterval == 0) this%saveInterval = -1

  this%reportInterval = getOption("report_interval", 1)
  if (this%reportInterval == 0) this%reportInterval = -1

  call this%timeIntegratorFactory%connect(timeIntegrator,                                    &
       trim(region%solverOptions%timeIntegratorType))
  assert(associated(timeIntegrator))
  call timeIntegrator%setup(region)

  ! If a target state file was specified, read the target state. Otherwise, initialize the
  ! target state to a quiescent state by default.
  if (region%simulationFlags%useTargetState) then
     filename = getOption("target_state_file", "")
     if (len_trim(filename) == 0) then
        do i = 1, size(region%states)
           call region%states(i)%makeQuiescent(size(region%globalGridSizes, 1),              &
                region%solverOptions%ratioOfSpecificHeats, region%states(i)%targetState)
        end do
     else
        call region%loadData(QOI_TARGET_STATE, filename)
     end if
  end if

  if (.not. region%simulationFlags%predictionOnly) then

     ! Initialize control mollifier.
     filename = getOption("control_mollifier_file", "")
     if (len_trim(filename) == 0) then
        do i = 1, size(region%grids)
           region%grids(i)%controlMollifier = 1.0_wp
        end do
     else
        call region%loadData(QOI_CONTROL_MOLLIFIER, filename)
     end if

     ! Target mollifier.
     filename = getOption("target_mollifier_file", "")
     if (len_trim(filename) == 0) then
        do i = 1, size(region%grids)
           region%grids(i)%targetMollifier = 1.0_wp
        end do
     else
        call region%loadData(QOI_TARGET_MOLLIFIER, filename)
     end if

  end if

  ! Setup boundary conditions.
  call getRequiredOption("boundary_condition_file", filename)
  call region%setupBoundaryConditions(filename)

  ! Compute damping strength on sponge patches.
  do i = 1, size(region%grids)
     call computeSpongeStrengths(region%patchFactories, region%grids(i))
  end do

  ! Check continuity at block interfaces.
  if (getOption("check_interface_continuity", .false.))                                      &
       call checkFunctionContinuityAtInterfaces(region, epsilon(0.0_wp))

  ! Update patches.
  do i = 1, size(region%grids)
     call updatePatchFactories(region%patchFactories, region%simulationFlags,                &
          region%solverOptions, region%grids(i), region%states(i))
  end do

  if (.not. region%simulationFlags%predictionOnly) then

     call this%controllerFactory%connect(controller,                                         &
          trim(region%solverOptions%controllerType))
     assert(associated(controller))
     call controller%setup(region)

     call this%functionalFactory%connect(functional,                                         &
          trim(region%solverOptions%costFunctionalType))
     assert(associated(functional))
     call functional%setup(region)

  end if

end subroutine setupSolver

subroutine cleanupSolver(this)

  ! <<< Derived types >>>
  use Solver_mod, only : t_Solver

  implicit none

  ! <<< Arguments >>>
  class(t_Solver) :: this

  call this%timeIntegratorFactory%cleanup()
  call this%controllerFactory%cleanup()
  call this%functionalFactory%cleanup()

end subroutine cleanupSolver

function runForward(this, region, actuationAmount, restartFilename,desiredPrecision) result(costFunctional)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver
  use Controller_mod, only : t_Controller
  use Functional_mod, only : t_Functional
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use TimeIntegrator_mod, only : t_TimeIntegrator

  ! <<< Enumerations >>>
  use State_enum, only : QOI_FORWARD_STATE, QOI_TIME_AVERAGED_STATE
  use Region_enum, only : FORWARD

  ! <<< Private members >>>
  use SolverImpl, only : showProgress, checkSolutionLimits, loadInitialCondition

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_Solver) :: this
  class(t_Region) :: region
  real(SCALAR_KIND), intent(in), optional :: actuationAmount
  real(SCALAR_KIND), intent(in), optional :: desiredPrecision
  character(len = *), intent(in), optional :: restartFilename

  ! <<< Result >>>
  SCALAR_TYPE :: costFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: filename
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()
  class(t_Controller), pointer :: controller => null()
  class(t_Functional), pointer :: functional => null()
  integer :: i, j, timestep, startTimestep
  real(wp) :: time, startTime, timeStepSize
  SCALAR_TYPE :: instantaneousCostFunctional

  call startTiming("runForward")

  costFunctional = 0.0_wp

  region%states(:)%actuationAmount = 0.0_wp
  if (present(actuationAmount) .and. .not. region%simulationFlags%predictionOnly)            &
       region%states(:)%actuationAmount = actuationAmount

  ! Connect to the previously allocated time integrator.
  call this%timeIntegratorFactory%connect(timeIntegrator)
  assert(associated(timeIntegrator))

  ! Connect to the previously allocated controller.
  if (.not. region%simulationFlags%predictionOnly) then
     call this%controllerFactory%connect(controller)
     assert(associated(controller))
  end if

  ! Connect to the previously allocated functional.
  if (.not. region%simulationFlags%predictionOnly) then
     call this%functionalFactory%connect(functional)
     assert(associated(functional))
  end if

  ! Setup residual manager if this is a steady-state simulation.
  if (region%simulationFlags%steadyStateSimulation)                                          &
       call this%residualManager%setup("", region)

  ! Load the initial condition.
  if (present(restartFilename)) then
     call loadInitialCondition(this, region, FORWARD, restartFilename)
  else
     call loadInitialCondition(this, region, FORWARD)
  end if

  startTimestep = region%timestep
  startTime = region%states(1)%time

  ! Save the initial condition if it was not specified as a restart file.
  if (.not. present(restartFilename)) then
     write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", startTimestep, ".q"
     call region%saveData(QOI_FORWARD_STATE, filename)
  end if

  ! Call controller hooks before time marching starts.
  if (.not. region%simulationFlags%predictionOnly .and.                                      &
       abs(region%states(1)%actuationAmount) > 0.0_wp)                                       &
       call controller%hookBeforeTimemarch(region, FORWARD)

  time = startTime
  region%states(:)%time = time
  do i = 1, size(region%states) !... update state
     call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
          region%solverOptions)
  end do

  do timestep = startTimestep + 1, startTimestep + this%nTimesteps

     region%timestep = timestep
     timeStepSize = region%getTimeStepSize()

     do i = 1, timeIntegrator%nStages

        ! Check if physical quantities are within allowed limits.
        if (region%simulationFlags%enableSolutionLimits)                                     &
             call checkSolutionLimits(region, FORWARD, this%outputPrefix)

        ! Update control forcing.
        if (.not. region%simulationFlags%predictionOnly .and.                                &
             abs(region%states(1)%actuationAmount) > 0.0_wp)                                 &
             call controller%updateForcing(region)

        ! Take a single sub-step using the time integrator.
        call timeIntegrator%substepForward(region, time, timeStepSize, timestep, i)

        do j = 1, size(region%states) !... update state
           call region%states(j)%update(region%grids(j), region%simulationFlags,             &
                region%solverOptions)
        end do

        ! Update the cost functional.
        if (.not. region%simulationFlags%predictionOnly) then
           instantaneousCostFunctional = functional%compute(region)
           costFunctional = costFunctional +                                                 &
                timeIntegrator%norm(i) * timeStepSize * instantaneousCostFunctional
        end if

        ! Update the time average.
        if (region%simulationFlags%computeTimeAverage) then
           do j = 1, size(region%states)
              region%states(j)%timeAverage = region%states(j)%timeAverage +                  &
                   timeIntegrator%norm(i) * timeStepSize *                                   &
                   region%states(j)%conservedVariables
           end do
        end if

     end do !... i = 1, timeIntegrator%nStages

     ! Report simulation progess.
     call showProgress(this, region, FORWARD, startTimestep, timestep,                       &
          time, instantaneousCostFunctional)

     ! Stop if this is a steady-state simulation and solution has converged.
     if (this%residualManager%hasSimulationConverged) exit

     ! Filter solution if required.
     if (region%simulationFlags%filterOn) then
        do j = 1, size(region%grids)
           call region%grids(j)%applyFilter(region%states(j)%conservedVariables, timestep)
        end do
     end if

     ! Add random precision nuber if set
     if (present(desiredPrecision)) then
        do j = 1, size(region%grids)
           call region%grids(j)%applyRandFluctuation(region%states(j)%conservedVariables,&
               desiredPrecision)
        end do
     end if 

  end do !... timestep = startTimestep + 1, startTimestep + this%nTimesteps

  ! Call controller hooks after time marching ends.
  if (.not. region%simulationFlags%predictionOnly .and.                                      &
       abs(region%states(1)%actuationAmount) > 0.0_wp)                                       &
       call controller%hookAfterTimemarch(region, FORWARD)

  call this%residualManager%cleanup()

  if (region%simulationFlags%computeTimeAverage) then
     do i = 1, size(region%states)
        region%states(i)%timeAverage = region%states(i)%timeAverage / (time - startTime)
     end do
     call region%saveData(QOI_TIME_AVERAGED_STATE, trim(this%outputPrefix) // ".mean.q")
  end if

  call endTiming("runForward")

end function runForward

function runAdjoint(this, region) result(costSensitivity)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver
  use Controller_mod, only : t_Controller
  use Functional_mod, only : t_Functional
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use TimeIntegrator_mod, only : t_TimeIntegrator
  use ReverseMigrator_mod, only : t_ReverseMigrator
  use ReverseMigrator_factory, only : t_ReverseMigratorFactory

  ! <<< Enumerations >>>
  use State_enum, only : QOI_ADJOINT_STATE, QOI_FORWARD_STATE
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Private members >>>
  use SolverImpl, only : showProgress, checkSolutionLimits, loadInitialCondition

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_Solver) :: this
  class(t_Region) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: costSensitivity

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: filename
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()
  class(t_Controller), pointer :: controller => null()
  class(t_Functional), pointer :: functional => null()
  type(t_ReverseMigratorFactory) :: reverseMigratorFactory
  class(t_ReverseMigrator), pointer :: reverseMigrator => null()
  integer :: i, j, timestep, startTimestep, timemarchDirection
  real(SCALAR_KIND) :: time, startTime, timeStepSize
  SCALAR_TYPE :: instantaneousCostSensitivity

  assert(.not. region%simulationFlags%predictionOnly)

  call startTiming("runAdjoint")

  costSensitivity = 0.0_wp

  ! Connect to the previously allocated time integrator.
  call this%timeIntegratorFactory%connect(timeIntegrator)
  assert(associated(timeIntegrator))

  ! Connect to the previously allocated controller.
  call this%controllerFactory%connect(controller)
  assert(associated(controller))

  ! Connect to the previously allocated functional.
  call this%functionalFactory%connect(functional)
  assert(associated(functional))

  ! Setup residual manager if this is a steady-state simulation
  if (region%simulationFlags%steadyStateSimulation)                                          &
       call this%residualManager%setup("adjoint_residuals", region)

  ! Load the initial condition.
  call loadInitialCondition(this, region, FORWARD) !... for control horizon end timestep.

  ! Load the adjoint coefficients corresponding to the end of the control time horizon.
  if (region%simulationFlags%steadyStateSimulation) then
     write(filename, '(2A)') trim(this%outputPrefix), ".steady_state.q"
  else
     write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-",                            &
          region%timestep + this%nTimesteps, ".q"
  end if
  call region%loadData(QOI_FORWARD_STATE, filename)
  do i = 1, size(region%states) !... update state
     call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
          region%solverOptions)
  end do

  startTimestep = region%timestep
  startTime = region%states(1)%time

  ! Connect to the previously allocated reverse migrator.
  call reverseMigratorFactory%connect(reverseMigrator,                                       &
       region%solverOptions%checkpointingScheme)
  assert(associated(reverseMigrator))

  ! Setup the revese-time migrator if this is not a steady-state simulation.
  if (.not. region%simulationFlags%steadyStateSimulation)                                    &
       call reverseMigrator%setup(region, timeIntegrator, this%outputPrefix,                 &
       startTimestep - this%nTimesteps, startTimestep, this%saveInterval,                    &
       this%saveInterval * timeIntegrator%nStages)

  ! March forward for adjoint steady-state simulation.
  timemarchDirection = -1
  if (region%simulationFlags%steadyStateSimulation) timemarchDirection = 1

  ! Adjoint initial condition (if specified).
  call loadInitialCondition(this, region, ADJOINT)

  write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", region%timestep, ".adjoint.q"
  call region%saveData(QOI_ADJOINT_STATE, filename)

  ! Call controller hooks before time marching starts.
  call controller%hookBeforeTimemarch(region, ADJOINT)

  time = startTime

  do timestep = startTimestep + sign(1, timemarchDirection),                                 &
       startTimestep + sign(this%nTimesteps, timemarchDirection), timemarchDirection

     region%timestep = timestep
     timeStepSize = region%getTimeStepSize()

     do i = timeIntegrator%nStages, 1, -1

        ! Load adjoint coefficients.
        if (.not. region%simulationFlags%steadyStateSimulation) then !... unsteady simulation.
           if (i == 1) then
              call reverseMigrator%migrateTo(region, timeIntegrator,                         &
                   timestep, timeIntegrator%nStages)
           else
              call reverseMigrator%migrateTo(region, timeIntegrator, timestep + 1, i - 1)
           end if
        end if

        ! Update gradient.
        call controller%updateGradient(region)

        ! Update cost sensitivity.
        instantaneousCostSensitivity = controller%computeSensitivity(region)
        costSensitivity = costSensitivity +                                                  &
             timeIntegrator%norm(i) * timeStepSize * instantaneousCostSensitivity

        ! Update adjoint forcing on cost target patches.
        call functional%updateAdjointForcing(region)

        ! Take a single adjoint sub-step using the time integrator.
        call timeIntegrator%substepAdjoint(region, time, timeStepSize, timestep, i)

        ! TODO: how to enforce limits on adjoint variables... check for NaN?
        if (region%simulationFlags%enableSolutionLimits)                                     &
             call checkSolutionLimits(region, ADJOINT, this%outputPrefix)

     end do

     ! Report simulation progress.
     call showProgress(this, region, ADJOINT, startTimestep, timestep,                       &
          time, instantaneousCostSensitivity)

     ! Stop if this is a steady-state simulation and solution has converged.
     if (this%residualManager%hasSimulationConverged) exit

     ! Filter solution if required.
     if (region%simulationFlags%filterOn) then
        do j = 1, size(region%grids)
           call region%grids(j)%applyFilter(region%states(j)%adjointVariables, timestep)
        end do
     end if

  end do !... timestep = startTimestep + sign(1, timemarchDirection), ...

  ! Call controller hooks after time marching ends.
  call controller%hookAfterTimemarch(region, ADJOINT)

  call this%residualManager%cleanup()
  call reverseMigratorFactory%cleanup()

  call endTiming("runAdjoint")

end function runAdjoint

subroutine checkGradientAccuracy(this, region)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver

  ! <<< Enumerations >>>
  use State_enum, only : QOI_FORWARD_STATE

  ! <<< Internal modules >>>
  use InputHelper, only : getFreeUnit, getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_Solver) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nIterations, restartIteration, fileUnit, iostat, procRank, ierror
  character(len = STRING_LENGTH) :: filename, message
  real(wp) :: actuationAmount, baselineCostFunctional, costFunctional, costSensitivity,      &
       initialActuationAmount, geometricGrowthFactor, gradientError,dummyValue,&
       gradientAccuracyPrecision

  call getRequiredOption("number_of_control_iterations", nIterations)
  if (nIterations < 0) then
     write(message, '(A)') "Number of control iterations must be a non-negative number!"
     call gracefulExit(region%comm, message)
  end if

  restartIteration = getOption("restart_control_iteration", 0)
  restartIteration = max(restartIteration, 0)

  if (restartIteration > 0 .and. .not. region%simulationFlags%isBaselineAvailable) then
     write(message, '(A)') "Can't restart with controlled prediction without baseline!"
     call gracefulExit(region%comm, message)
  end if

  call MPI_Comm_rank(region%comm, procRank, ierror)

  write(filename, '(2A)') trim(this%outputPrefix), ".gradient_error.txt"
  if (procRank == 0) then
     if (restartIteration == 0 .and. .not. region%simulationFlags%isBaselineAvailable) then
        open(unit = getFreeUnit(fileUnit), file = trim(filename), action = 'write',          &
             status = 'unknown', iostat = iostat)
     else
        open(unit = getFreeUnit(fileUnit), file = trim(filename), action = 'readwrite',      &
             status = 'old', position = 'rewind', iostat = iostat)
     end if
  end if

  call MPI_Bcast(iostat, 1, MPI_INTEGER, 0, region%comm, ierror)
  if (iostat /= 0) then
     write(message, "(2A)") trim(filename), ": Failed to open file for writing!"
     call gracefulExit(region%comm, message)
  end if

  if (nIterations > 0) then
     call getRequiredOption("initial_actuation_amount", initialActuationAmount)
     call getRequiredOption("precision_for_gradient_accuracy",gradientAccuracyPrecision)
     if (nIterations > 1) then
        call getRequiredOption("actuation_amount_geometric_growth", geometricGrowthFactor)
     else
        geometricGrowthFactor = getOption("actuation_amount_geometric_growth", 1.0_wp)
     end if
  end if

  ! Find (or load from file) the cost functional for the baseline prediction.
  if (region%simulationFlags%isBaselineAvailable) then
     if (procRank == 0)                                                                      &
          read(fileUnit, *, iostat = iostat) i, actuationAmount,                             &
          baselineCostFunctional, costSensitivity, gradientError
     call MPI_Bcast(iostat, 1, MPI_INTEGER, 0, region%comm, ierror)
     if (iostat /= 0) then
        write(message, "(2A)") trim(filename),                                               &
             ": Failed to read baseline cost functional from file!"
        call gracefulExit(region%comm, message)
     end if
     call MPI_Bcast(baselineCostFunctional, 1, REAL_TYPE_MPI, 0, region%comm, ierror)
  else
     baselineCostFunctional = this%runForward(region)
  end if

  ! Find the sensitivity gradient (this is the only time the adjoint simulation will be run).
  if (restartIteration == 0) then
     costSensitivity = this%runAdjoint(region)
  else
     call MPI_Bcast(costSensitivity, 1, REAL_TYPE_MPI, 0, region%comm, ierror)
  end if

  if (procRank == 0 .and. .not. region%simulationFlags%isBaselineAvailable)                  &
       write(fileUnit, '(I4,4(1X,SP,' // SCALAR_FORMAT // '))') 0, 0.0_wp,                   &
       baselineCostFunctional, costSensitivity, 0.0_wp

  if (nIterations == 0) return

  ! Turn off output for controlled predictions.
  region%outputOn = .false.

  if (restartIteration == 0) restartIteration = restartIteration + 1

  do i = 1, restartIteration - 1
     if (procRank == 0)                                                                      &
          read(fileUnit, *, iostat = iostat) j, actuationAmount, costFunctional,             &
          dummyValue, gradientError
     call MPI_Bcast(iostat, 1, MPI_INTEGER, 0, region%comm, ierror)
     if (iostat /= 0) then
        write(message, "(2A)") trim(filename),                                               &
             ": Cost functional history is too short for the specified restart iteration!"
        call gracefulExit(region%comm, message)
     end if
  end do

  do i = restartIteration, restartIteration + nIterations - 1
     actuationAmount = initialActuationAmount * geometricGrowthFactor ** real(i - 1, wp)
     costFunctional = this%runForward(region, actuationAmount =&
          actuationAmount,desiredPrecision=gradientAccuracyPrecision)
     gradientError = (costFunctional - baselineCostFunctional) / actuationAmount +           &
          costSensitivity
     if (procRank == 0)                                                                      &
          write(fileUnit, '(I4,4(1X,SP,' // SCALAR_FORMAT // '))') i, actuationAmount,       &
          costFunctional, -(costFunctional - baselineCostFunctional) / actuationAmount,      &
          abs(gradientError)
  end do

  if (procRank == 0) close(fileUnit)

end subroutine checkGradientAccuracy
