#include "config.h"

module SolverImpl

  implicit none
  public

contains

  subroutine showProgress(this, region, mode, startTimestep,                                 &
       timestep, time, instantaneousFunctional, controlIteration)

    ! <<< External modules >>>
    use iso_fortran_env, only : output_unit
    use MPI

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
    integer, intent(in), optional :: controlIteration
    real(SCALAR_KIND), intent(in) :: time, instantaneousFunctional

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, nDimensions, ierror, controlIteration_
    real(SCALAR_KIND) :: timeStepSize, cfl, maxTemperature
    character(len = STRING_LENGTH) :: str, str_, filename
    class(t_Controller), pointer :: controller => null()
    class(t_Functional), pointer :: functional => null()

    assert_key(mode, (FORWARD, ADJOINT))
    assert(startTimestep >= 0)

    nDimensions = size(region%globalGridSizes, 1)
    assert_key(nDimensions, (1, 2, 3))

    controlIteration_ = 0
    if (present(controlIteration)) controlIteration_ = controlIteration
    assert(controlIteration_ >= 0)

    if (this%reportInterval > 0 .and. mod(timestep, max(1, this%reportInterval)) == 0) then

       timeStepSize = region%getTimeStepSize()
       cfl = region%getCfl()

       do i = 1, size(region%states)
          maxTemperature = max(maxTemperature,maxval(region%states(i)%temperature(:,1)))
       end do
       call MPI_Allreduce(MPI_IN_PLACE, maxTemperature, 1, REAL_TYPE_MPI, MPI_MAX,           &
            region%comm, ierror)
       maxTemperature = maxTemperature *                                                     &
            (region%solverOptions%ratioOfSpecificHeats - 1.0_wp) * 293.15_wp

       if (timestep <= this%reportInterval) then
          write(str,'(A12,A2,A12,A2,3A12)') 'Step','  ','Time','   ','dt','CFL','Tmax [K]'
          call writeAndFlush(region%comm, output_unit, str)
       end if

       write(str, '(I12,A2,1ES12.5,A2,1ES12.5,1F12.4,1F12.4)') timestep, '   ', abs(time),   &
            '   ',timeStepSize, cfl, maxTemperature

       if (.not. region%simulationFlags%predictionOnly) then

          select case (mode)
          case (FORWARD)
             if (controlIteration_ > 0) then
                write(str_, '(A,1ES12.5,A,I2)') ", cost = ", instantaneousFunctional,        &
                     " iteration = ",controlIteration_
             else
                write(str_, '(A,1ES12.5)') ", cost = ", instantaneousFunctional
             end if
          case (ADJOINT)
             write(str_, '(A,1ES12.5)') ", gradient = ", instantaneousFunctional
          end select

          str = trim(str) // trim(str_)

       end if

       call writeAndFlush(region%comm, output_unit, str)

       if (.not. region%simulationFlags%predictionOnly .and. region%outputOn) then

          select case (mode)

          case (FORWARD)
             call this%functionalFactory%connect(functional)
             assert(associated(functional))
             if (controlIteration_ > 0) then
                write(filename, '(2A,I2.2,A)') trim(this%outputPrefix), ".cost_functional_", &
                     controlIteration_, ".txt"
                call functional%writeToFile(region%comm, filename, timestep, time,           &
                     timestep - startTimestep > this%reportInterval)
             else
                call functional%writeToFile(region%comm, trim(this%outputPrefix) //          &
                     ".cost_functional.txt", timestep, time,                                 &
                     timestep - startTimestep > this%reportInterval)
             end if

          case (ADJOINT)
             call this%controllerFactory%connect(controller)
             assert(associated(controller))
             if (controlIteration_ > 0) then
                write(filename, '(2A,I2.2,A)') trim(this%outputPrefix), ".cost_sensitivity_",&
                     controlIteration_, ".txt"
                call controller%writeSensitivityToFile(region%comm, filename, timestep, time,&
                     timestep - startTimestep > this%reportInterval)
             else
                call controller%writeSensitivityToFile(region%comm, trim(this%outputPrefix)  &
                     // ".cost_sensitivity.txt", timestep, time,                             &
                     startTimestep - timestep > this%reportInterval)
             end if

          end select

       end if

    end if

    if (this%saveInterval > 0 .and. mod(timestep, max(1, this%saveInterval)) == 0) then

       select case (mode)
       case (FORWARD)
          if (controlIteration_ > 0) then
             write(filename, '(2A,I2.2,A,I8.8,A)') trim(this%outputPrefix), "_",             &
                  controlIteration_, "-", timestep, ".q"
             call region%saveData(QOI_FORWARD_STATE, filename)
          else
             write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", timestep, ".q"
             call region%saveData(QOI_FORWARD_STATE, filename)
          end if
       case (ADJOINT)
          if (controlIteration_ > 0) then
             write(filename, '(2A,I2.2,A,I8.8,A)') trim(this%outputPrefix), "_",             &
                  controlIteration_, "-", timestep, ".adjoint.q"
             call region%saveData(QOI_ADJOINT_STATE, filename)
          else
             write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", timestep,          &
                  ".adjoint.q"
             call region%saveData(QOI_ADJOINT_STATE, filename)
          end if
       end select

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
    integer :: i, k, iGlobal, jGlobal, kGlobal, rankReportingError, procRank, ierror
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

       do k = 1, region%solverOptions%nSpecies
          if (.not. region%grids(i)%isVariableWithinRange(region%states(i)%massFraction(:,k),&
               fOutsideRange, iGlobal, jGlobal, kGlobal,                                     &
               minValue = region%solverOptions%massFractionRange(1),                         &
               maxValue = region%solverOptions%massFractionRange(2))) then
             write(message, '(4(A,I0.0),3(A,(SS,ES9.2E2)),A)') "Mass fraction on grid ",     &
                  region%grids(i)%index,                                                     &
                  " at (", iGlobal, ", ", jGlobal, ", ", kGlobal, "): ",                     &
                  fOutsideRange, " out of range (",                                          &
                  region%solverOptions%massFractionRange(1), ", ",                            &
                  region%solverOptions%massFractionRange(2), ")!"
             rankReportingError = procRank
             exit
          end if
       end do

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

  if (region%simulationFlags%outputToEnsight) then
     call getRequiredOption("ensight_frequency", this%ensightFrequency)
     this%ensightSave = int(region%states(1)%time / this%ensightFrequency)
     allocate(this%ensight(size(region%grids)))
  end if

  this%adjointIterations = getOption("adjoint_iterations", this%nTimesteps)
  this%adjointIterations = max(0, this%adjointIterations)
  this%adjointIterations = min(this%nTimesteps, this%adjointIterations)

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
                region%solverOptions%nSpecies, region%solverOptions%ratioOfSpecificHeats,    &
                region%states(i)%targetState)
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
  SAFE_DEALLOCATE(this%ensight)

end subroutine cleanupSolver

function runForward(this, region, actuationAmount, controlIteration, restartFilename)        &
     result(costFunctional)

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
  integer, intent(in), optional :: controlIteration
  real(SCALAR_KIND), intent(in), optional :: actuationAmount
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
     functional%runningTimeQuadrature = 0.0_wp
  end if

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

  ! Setup EnSight output.
  if (region%simulationFlags%outputToEnsight) then
     do i = 1, size(region%grids)
        call this%ensight(i)%setup(region%grids(i), i, region%states(i)%time)
     end do
     ! Output the initial condition to EnSight.
     if (startTimestep == 0) then
        do i = 1, size(region%states)
           call this%ensight(i)%output(region%states(i), region%grids(i)%comm, i,            &
                FORWARD, time, region%solverOptions%nSpecies)
        end do
     end if
  end if

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
           functional%runningTimeQuadrature = functional%runningTimeQuadrature +             &
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

     ! Filter solution if required.
     if (region%simulationFlags%filterOn) then
        do j = 1, size(region%grids)
           call region%grids(j)%applyFilter(region%states(j)%conservedVariables, timestep)
        end do
     end if

     ! Report simulation progess.
     if (present(controlIteration)) then
        call showProgress(this, region, FORWARD, startTimestep, timestep,                    &
             time, instantaneousCostFunctional, controlIteration)
     else
        call showProgress(this, region, FORWARD, startTimestep, timestep,                    &
             time, instantaneousCostFunctional)
     end if

     ! Ensight output.
     if (region%simulationFlags%outputToEnsight) then
        if (int(time / this%ensightFrequency) .ne. this%ensightSave) then
           this%ensightSave = int(time / this%ensightFrequency)
           do i = 1, size(region%states)
              call this%ensight(i)%output(region%states(i), region%grids(i)%comm, i,         &
                   FORWARD, time, region%solverOptions%nSpecies)
           end do
        end if
     end if

  end do !... timestep = startTimestep + 1, startTimestep + this%nTimesteps

  ! Call controller hooks after time marching ends.
  if (.not. region%simulationFlags%predictionOnly .and.                                      &
       abs(region%states(1)%actuationAmount) > 0.0_wp)                                       &
       call controller%hookAfterTimemarch(region, FORWARD)

  if (region%simulationFlags%computeTimeAverage) then
     do i = 1, size(region%states)
        region%states(i)%timeAverage = region%states(i)%timeAverage / (time - startTime)
     end do
     call region%saveData(QOI_TIME_AVERAGED_STATE, trim(this%outputPrefix) // ".mean.q")
  end if

  if (.not. region%simulationFlags%predictionOnly)                                           &
       costFunctional = functional%runningTimeQuadrature

  call endTiming("runForward")

end function runForward

function runAdjoint(this, region, controlIteration) result(costSensitivity)

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
  integer, intent(in), optional :: controlIteration

  ! <<< Result >>>
  SCALAR_TYPE, dimension(:), allocatable :: costSensitivity

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: filename
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()
  class(t_Controller), pointer :: controller => null()
  class(t_Functional), pointer :: functional => null()
  type(t_ReverseMigratorFactory) :: reverseMigratorFactory
  class(t_ReverseMigrator), pointer :: reverseMigrator => null()
  integer :: i, j, timestep, startTimestep, timemarchDirection
  SCALAR_TYPE, dimension(:), allocatable :: instantaneousCostSensitivity
  real(SCALAR_KIND) :: time, startTime, timeStepSize

  assert(.not. region%simulationFlags%predictionOnly)

  call startTiming("runAdjoint")

  ! Connect to the previously allocated time integrator.
  call this%timeIntegratorFactory%connect(timeIntegrator)
  assert(associated(timeIntegrator))

  ! Connect to the previously allocated controller.
  call this%controllerFactory%connect(controller)
  assert(associated(controller))

  ! Connect to the previously allocated functional.
  call this%functionalFactory%connect(functional)
  assert(associated(functional))

  ! Load the initial condition.
  call loadInitialCondition(this, region, FORWARD) !... for control horizon end timestep.

  ! Load the adjoint coefficients corresponding to the end of the control time horizon.
  write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-",                               &
       region%timestep + this%adjointIterations, ".q"
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

  ! Setup the revese-time migrator.
  call reverseMigrator%setup(region, timeIntegrator, this%outputPrefix,                      &
       startTimestep - this%adjointIterations, startTimestep, this%saveInterval,             &
       this%saveInterval * timeIntegrator%nStages)
  timemarchDirection = -1

  ! Adjoint initial condition (if specified).
  call loadInitialCondition(this, region, ADJOINT)

  write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", region%timestep, ".adjoint.q"
  call region%saveData(QOI_ADJOINT_STATE, filename)

  ! Call controller hooks before time marching starts.
  call controller%hookBeforeTimemarch(region, ADJOINT)

  allocate(costSensitivity(controller%nParameters))
  allocate(instantaneousCostSensitivity(controller%nParameters))
  costSensitivity = 0.0_wp
  instantaneousCostSensitivity = 0.0_wp

  time = startTime

  do timestep = startTimestep + sign(1, timemarchDirection),                                 &
       startTimestep + sign(this%adjointIterations, timemarchDirection), timemarchDirection

     region%timestep = timestep
     timeStepSize = region%getTimeStepSize()

     do i = timeIntegrator%nStages, 1, -1

        ! Load adjoint coefficients.
        if (i == 1) then
           call reverseMigrator%migrateTo(region, timeIntegrator,                            &
                timestep, timeIntegrator%nStages)
        else
           call reverseMigrator%migrateTo(region, timeIntegrator, timestep + 1, i - 1)
        end if

        ! Update gradient.
        call controller%updateGradient(region)

        ! Update cost sensitivity.
        call controller%computeSensitivity(region)
        instantaneousCostSensitivity = controller%cachedValue
        controller%runningTimeQuadrature = controller%runningTimeQuadrature +                &
             timeIntegrator%norm(i) * timeStepSize * instantaneousCostSensitivity

        ! Update adjoint forcing on cost target patches.
        call functional%updateAdjointForcing(region)

        ! Take a single adjoint sub-step using the time integrator.
        call timeIntegrator%substepAdjoint(region, time, timeStepSize, timestep, i)

        ! TODO: how to enforce limits on adjoint variables... check for NaN?
        if (region%simulationFlags%enableSolutionLimits)                                     &
             call checkSolutionLimits(region, ADJOINT, this%outputPrefix)

     end do

     ! Report simulation progess.
     if (present(controlIteration)) then
        call showProgress(this, region, ADJOINT, startTimestep, timestep,                    &
             time, sum(instantaneousCostSensitivity**2), controlIteration)
     else
        call showProgress(this, region, ADJOINT, startTimestep, timestep,                    &
             time, sum(instantaneousCostSensitivity**2))
     end if

     ! Filter solution if required.
     if (region%simulationFlags%filterOn) then
        do j = 1, size(region%grids)
           call region%grids(j)%applyFilter(region%states(j)%adjointVariables, timestep)
        end do
     end if

  end do !... timestep = startTimestep + sign(1, timemarchDirection), ...

  ! Call controller hooks after time marching ends.
  call controller%hookAfterTimemarch(region, ADJOINT)

  call reverseMigratorFactory%cleanup()

  SAFE_DEALLOCATE(instantaneousCostSensitivity)

  costSensitivity = controller%runningTimeQuadrature

  call endTiming("runAdjoint")

end function runAdjoint

subroutine checkGradientAccuracy(this, region)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver
  use Controller_mod, only : t_Controller

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
  class(t_Controller), pointer :: controller => null()
  real(wp) :: actuationAmount, baselineCostFunctional, costFunctional, costSensitivity,      &
       initialActuationAmount, geometricGrowthFactor, gradientError, dummyValue
  real(WP), dimension(:), allocatable :: individualSensitivities
  logical :: minimizeCost

  ! Connect to the previously allocated controller.
  call this%controllerFactory%connect(controller)
  assert(associated(controller))

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
     region%states(:)%gradientDirection = -1
     minimizeCost = getOption("minimize_cost_functional", .true.)
     if (.not. minimizeCost) region%states(:)%gradientDirection = 1
     if (nIterations > 1) then
        call getRequiredOption("actuation_amount_geometric_growth", geometricGrowthFactor)
     else
        geometricGrowthFactor = getOption("actuation_amount_geometric_growth", 1.0_wp)
     end if
  end if

  ! Find (or load from file) the cost functional & sensitivity for the baseline prediction.
  allocate(individualSensitivities(controller%nParameters))
  if (region%simulationFlags%isBaselineAvailable) then
     if (procRank == 0) then
        read(fileUnit, *, iostat = iostat)
        read(fileUnit, *, iostat = iostat) i, actuationAmount,                               &
             baselineCostFunctional, costSensitivity, gradientError, individualSensitivities
     end if
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
     individualSensitivities = this%runAdjoint(region)
  else
     call MPI_Bcast(individualSensitivities, size(individualSensitivities), REAL_TYPE_MPI, 0,&
          region%comm, ierror)
  end if
  assert(size(individualSensitivities) == controller%nParameters)
  costSensitivity = sum(individualSensitivities**2)

  ! Store gradient to be used for control forcing.
  do i = 1, size(region%grids)
     allocate(region%states(i)%controlGradient(controller%nParameters))
     region%states(i)%controlGradient = individualSensitivities
  end do

   if (procRank == 0 .and. .not. region%simulationFlags%isBaselineAvailable) then
     write(fileUnit, '(A4,100A24)') 'i', 'Actuation amount', 'Cost functional',              &
          'Cost sensitivity', 'Gradient error',                                              &
          ('dJ/d '//trim(controller%sensitivityParameter(i)), i = 1, controller%nParameters)
     write(fileUnit, '(I4,100(1X,SP,' // SCALAR_FORMAT // '))') 0, 0.0_wp,                   &
          baselineCostFunctional, costSensitivity, 0.0_wp,                                   &
          (individualSensitivities(i), i = 1, controller%nParameters)
  end if

  if (nIterations == 0) return

  ! Turn off output for controlled predictions.
  region%outputOn = getOption("output_control_iterations", .false.)
  region%simulationFlags%outputToEnsight = .false.

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
     costFunctional = this%runForward(region, actuationAmount = actuationAmount,             &
          controlIteration = i)
     gradientError = abs( abs(costFunctional - baselineCostFunctional) / actuationAmount -   &
          costSensitivity )
     if (procRank == 0) then
        write(fileUnit, '(I4,4(1X,SP,' // SCALAR_FORMAT // '))') i, actuationAmount,         &
          costFunctional, (costFunctional - baselineCostFunctional) / actuationAmount,       &
          gradientError
        flush(fileUnit)
     end if
  end do

  if (procRank == 0) close(fileUnit)

end subroutine checkGradientAccuracy

subroutine findOptimalForcing(this, region)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver
  use Controller_mod, only : t_Controller
  use Functional_mod, only : t_Functional

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
  integer :: i, j, nIterations, restartIteration, controlIteration, nForward, nAdjoint,      &
       p1, p2, parameterDirection, fileUnit, iostat, procRank, ierror
  character(len = STRING_LENGTH) :: optimizationType, filename, message
  class(t_Controller), pointer :: controller => null()
  class(t_Functional), pointer :: functional => null()
  real(wp) :: baselineCostFunctional, costFunctional, previousCostFunctional,                &
       indicatorFunction, costSensitivity, actuationAmount, previousActuationAmount,         &
       baselineActuationAmount, tangentActuationAmount, burnValue, minimumTolerance
  real(WP), dimension(:), allocatable :: individualSensitivities, parameters, tangentVector
  logical:: done, foundNewMinimum, burning, minimizeParameter

  ! Connect to the previously allocated controller.
  call this%controllerFactory%connect(controller)
  assert(associated(controller))

  ! Connect to the previously allocated functional.
  call this%functionalFactory%connect(functional)
  assert(associated(functional))

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

  write(filename, '(2A)') trim(this%outputPrefix), ".cost_optimization.txt"
  if (procRank == 0) then
     if (restartIteration == 0 .and. .not.region%simulationFlags%isBaselineAvailable) then
        open(unit = getFreeUnit(fileUnit), file = trim(filename), action ='write',           &
             status = 'unknown', iostat = iostat)
     else
        open(unit = getFreeUnit(fileUnit), file = trim(filename), action ='readwrite',       &
             status = 'old', position = 'rewind', iostat = iostat)
     end if
  end if

  call MPI_Bcast(iostat, 1, MPI_INTEGER, 0, region%comm, ierror)
  if (iostat /= 0) then
     write(message, "(2A)") trim(filename), ": Failed to open file for writing!"
     call gracefulExit(region%comm, message)
  end if

  ! Find (or load from file) useful data from the baseline prediction.
  allocate(parameters(controller%nParameters))
  allocate(individualSensitivities(controller%nParameters))
  if (region%simulationFlags%isBaselineAvailable) then
     if (procRank == 0) then
        read(fileUnit, *, iostat = iostat)
        do i = 1, restartIteration - 1
           read(fileUnit, *, iostat = iostat) j, actuationAmount,                            &
                baselineCostFunctional, indicatorFunction, parameters, individualSensitivities
        end do
     end if
     call MPI_Bcast(iostat, 1, MPI_INTEGER, 0, region%comm, ierror)
     if (iostat /= 0) then
        write(message, "(2A)") trim(filename),                                               &
             ": Failed to read baseline cost functional from file!"
        call gracefulExit(region%comm, message)
     end if
     call MPI_Bcast(baselineCostFunctional, 1, REAL_TYPE_MPI, 0, region%comm, ierror)
     call MPI_Bcast(indicatorFunction, 1, REAL_TYPE_MPI, 0, region%comm, ierror)
     call MPI_Bcast(parameters, 1, REAL_TYPE_MPI, size(parameters), region%comm, ierror)
     controller%baselineValue = parameters
  else
     baselineCostFunctional = this%runForward(region)
     indicatorFunction = functional%auxilaryFunctional
  end if

  ! Find the initial sensitivity gradient.
  if (restartIteration == 0) then
     individualSensitivities = this%runAdjoint(region)
  else
     call MPI_Bcast(individualSensitivities, size(individualSensitivities), REAL_TYPE_MPI,   &
          0, region%comm, ierror)
  end if
  assert(size(individualSensitivities) == controller%nParameters)
  costSensitivity = sum(individualSensitivities**2)

  if (restartIteration == 0) restartIteration = restartIteration + 1

  ! Find the previous actuation amount.
  if (restartIteration <= 1) then
     call getRequiredOption("initial_actuation_amount", actuationAmount)
  else
     call MPI_Bcast(actuationAmount, 1, REAL_TYPE_MPI, 0, region%comm, ierror)
  end if

  ! Store gradient to be used for control forcing.
  do i = 1, size(region%grids)
     allocate(region%states(i)%controlGradient(controller%nParameters))
     region%states(i)%controlGradient = individualSensitivities
  end do

  if (procRank == 0 .and. .not. region%simulationFlags%isBaselineAvailable) then
     write(fileUnit, '(A4,1000A24)') 'i', 'Actuation amount', 'Cost functional',             &
          'Inidicator functional',                                                           &
          (trim(controller%sensitivityParameter(i)), i = 1, controller%nParameters),         &
          ('dJ/d ' // trim(controller%sensitivityParameter(i)), i = 1,                       &
          controller%nParameters)
     write(fileUnit, '(I4,1000(1X,SP,' // SCALAR_FORMAT // '))') 0, 0.0_wp,                  &
          baselineCostFunctional, indicatorFunction,                                         &
          (controller%baselineValue(i), i = 1, controller%nParameters),                      &
          (individualSensitivities(i), i = 1, controller%nParameters)
  end if

  if (nIterations == 0) return

  ! Turn off output for controlled predictions.
  region%outputOn = getOption("output_control_iterations", .false.)
  region%simulationFlags%outputToEnsight = .false.

  call getRequiredOption("optimization_type", optimizationType)

  select case(trim(optimizationType))

  case ('IGNITION_BOUNDARY')

     ! Determine if the initial run was burning based on the last value of the instantaneous
     ! cost functional.
     call getRequiredOption("burn_value", burnValue)
     burning = indicatorFunction > burnValue

     ! Normalize the gradient.
     do i = 1, size(region%grids)
        region%states(i)%controlGradient = region%states(i)%controlGradient /                &
             sqrt(costSensitivity)
     end do

     ! We have at this point a baseline cost functional and the first gradient with respect
     ! to each parameter.
     nForward = 1
     nAdjoint = 1
     done = .false.
     controlIteration = restartIteration
     minimumTolerance = getOption("minimum_actuation_tolerance", 1.0E-9_wp)
     do while (controlIteration < nIterations .and. .not.done)

        ! Perform line search.
        do i = controlIteration, restartIteration + nIterations - 1

           ! Choose a direction to march.
           region%states(:)%gradientDirection = -1
           if (.not. burning) region%states(:)%gradientDirection = 1

           ! Compute a new cost functional.
           costFunctional = this%runForward(region, actuationAmount = actuationAmount,       &
                controlIteration = nForward)
           nForward = nForward + 1
           indicatorFunction = functional%auxilaryFunctional
           controlIteration = controlIteration + 1

          ! Output progress.
           if (procRank == 0) then
              write(fileUnit, '(I4,1000(1X,SP,' // SCALAR_FORMAT // '))') i,                 &
                   actuationAmount, costFunctional, indicatorFunction,                       &
                   (controller%baselineValue(j) +                                            &
                   real(region%states(1)%gradientDirection, wp) * actuationAmount *          &
                   region%states(1)%controlGradient(j), j = 1, controller%nParameters),      &
                   (individualSensitivities(j), j = 1, controller%nParameters)
              flush(fileUnit)
           end if

           ! Exit loop and recompute the gradient if we didn't passed the threshold.
           if (burning .and. indicatorFunction > burnvalue) then
              exit
           else if (.not.burning .and. indicatorFunction < burnValue) then
              exit
           end if

           ! If we made it this far, reduce the actuation amount and try again.
           actuationAmount = 0.5_wp * actuationAmount

           ! Check actuation tolerance.
           if (actuationAmount < minimumTolerance) done = .true.

        end do

        ! Update the baseline values and compute a new sensitivity gradient.
        if (.not.done .and. controlIteration < nIterations) then
           do i = 1, controller%nParameters
              controller%baselineValue(i) = controller%baselineValue(i) +                    &
                   real(region%states(1)%gradientDirection, wp) * actuationAmount *          &
                   region%states(1)%controlGradient(i)
           end do
           individualSensitivities = this%runAdjoint(region, controlIteration = nAdjoint)
           nAdjoint = nAdjoint + 1
           costSensitivity = sum(individualSensitivities**2)
           do i = 1, size(region%grids)
              region%states(i)%controlGradient = individualSensitivities /                   &
                   sqrt(costSensitivity)
           end do
        end if

     end do

     ! Dump optimizatiom summary and close.
     if (procRank == 0) then
        write(fileUnit, *) ''
        write(fileUnit, '(A28,I4)') 'Number of forward runs:', nForward
        write(fileUnit, '(A28,I4)') 'Number of adjoint runs:', nAdjoint
        write(fileUnit, '(A28,L4)') 'Ignition threshold found:', done
        close(fileUnit)
     end if

  case ('MAP_IGNITION')

     ! Assume we are at the burn boundary.
     ! We have at this point a baseline cost functional and a gradient with respect
     ! to each parameter.
     nForward = 1
     nAdjoint = 1
     controlIteration = restartIteration
     minimumTolerance = getOption("minimum_actuation_tolerance", 1.0E-9_wp)
     call getRequiredOption("burn_value", burnValue)
     tangentActuationAmount = actuationAmount
     allocate(tangentVector(controller%nParameters))

     ! Determine which parameters to control/adjust while keeping all others constant.
     call getRequiredOption("control_parameter", p1, region%comm)
     call getRequiredOption("adjust_parameter", p2, region%comm)
     assert(p1 <= controller%nParameters)
     assert(p2 <= controller%nParameters)
     minimizeParameter = getOption("minimize_control_parameter", .true.)
     parameterDirection = 1
     if (minimizeParameter) parameterDirection = -1

     ! Perform line search.
     do i = controlIteration, restartIteration + nIterations - 1

        ! Compute the tangent vector.
        tangentVector = 0.0_wp
        tangentVector(p1) = real(parameterDirection, wp)
        tangentVector(p2) = - tangentVector(p1) * individualSensitivities(p1) /              &
             individualSensitivities(p2)
        do j = 1, size(region%grids)
           region%states(j)%controlGradient = tangentVector / sqrt(sum(tangentVector**2))
        end do

        ! Compute a new cost functional in the tangent space.
        region%states(:)%gradientDirection = 1
        costFunctional = this%runForward(region, actuationAmount = tangentActuationAmount,   &
             controlIteration = nForward)
        nForward = nForward + 1
        indicatorFunction = functional%auxilaryFunctional
        burning = indicatorFunction > burnValue

        ! Update the baseline values.
        do j = 1, controller%nParameters
           controller%baselineValue(j) = controller%baselineValue(j) +                       &
                real(region%states(1)%gradientDirection, wp) * tangentActuationAmount *      &
                region%states(1)%controlGradient(j)
        end do

        ! Check if we are close to the ignition boundary and correct if needed.
        actuationAmount = minimumTolerance
        previousActuationAmount = 0.0_wp
        foundNewMinimum = .false.
        do j = 1, size(region%grids)
           region%states(j)%controlGradient = individualSensitivities / sqrt(costSensitivity)
        end do
        if (burning) region%states(:)%gradientDirection = -1
        do
           costFunctional = this%runForward(region, actuationAmount = actuationAmount,       &
                controlIteration = nForward)
           nForward = nForward + 1
           indicatorFunction = functional%auxilaryFunctional
           if ((burning .and. indicatorFunction < burnvalue) .or. &
                (.not.burning .and. indicatorFunction > burnvalue)) then
              exit
           end if
           foundNewMinimum = .true.
           previousActuationAmount = actuationAmount
           actuationAmount = actuationAmount + minimumTolerance
        end do

        if (foundNewMinimum) then
           ! Update the baseline values and compute a new sensitivity gradient.
           do j = 1, controller%nParameters
              controller%baselineValue(j) = controller%baselineValue(j) +                    &
                   real(region%states(1)%gradientDirection, wp) * previousActuationAmount *  &
                   region%states(1)%controlGradient(j)
           end do
           individualSensitivities = this%runAdjoint(region, controlIteration = nAdjoint)
           nAdjoint = nAdjoint + 1
           costSensitivity = sum(individualSensitivities**2)
        end if

        ! Output progress.
        if (procRank == 0) then
           write(fileUnit, '(I4,1000(1X,SP,' // SCALAR_FORMAT // '))') i,                    &
                actuationAmount, costFunctional, indicatorFunction,                          &
                (controller%baselineValue(j), j = 1, controller%nParameters),                &
                (individualSensitivities(j), j = 1, controller%nParameters)
           flush(fileUnit)
        end if

        controlIteration = controlIteration + 1

     end do

     ! Dump optimizatiom summary and close.
     if (procRank == 0) then
        write(fileUnit, *) ''
        write(fileUnit, '(A28,I4)') 'Number of forward runs:', nForward
        write(fileUnit, '(A28,I4)') 'Number of adjoint runs:', nAdjoint
        close(fileUnit)
     end if

  case ('MINIMIZATION')

     ! We have at this point a baseline cost functional and the first gradient with respect
     ! to each parameter.
     nForward = 1
     nAdjoint = 1
     done = .false.
     foundNewMinimum = .false.
     controlIteration = restartIteration
     actuationAmount = actuationAmount / baselineCostFunctional
     baselineActuationAmount = actuationAmount
     minimumTolerance = getOption("minimum_actuation_tolerance", 1.0E-9_wp)
     previousCostFunctional = baselineCostFunctional
     do while (controlIteration < nIterations .and. .not.done)

        ! Perform line search.
        do i = controlIteration + 1, restartIteration + nIterations

           ! Compute a new cost functional.
           costFunctional = this%runForward(region, actuationAmount = actuationAmount,       &
                controlIteration = nForward)
           indicatorFunction = functional%auxilaryFunctional
           nForward = nForward + 1
           controlIteration = controlIteration + 1

          ! Output progress.
           if (procRank == 0) then
              write(fileUnit, '(I4,1000(1X,SP,' // SCALAR_FORMAT // '))') i,                 &
                   actuationAmount, costFunctional, indicatorFunction,                       &
                   (controller%baselineValue(j) +                                            &
                   real(region%states(1)%gradientDirection, wp) * actuationAmount *          &
                   individualSensitivities(j), j = 1, controller%nParameters),               &
                   (individualSensitivities(j), j = 1, controller%nParameters)
              flush(fileUnit)
           end if

           ! Reset the actuation amount and exit if we start climbing back up the gradient.
           if (foundNewMinimum .and. costFunctional > previousCostFunctional) then
              actuationAmount = previousActuationAmount
              exit
           end if

           ! Check if a new minimum was found.
           if (costFunctional < previousCostFunctional) then
              foundNewMinimum = .true.
              previousCostFunctional = costFunctional
              previousActuationAmount = actuationAmount
           end if

           ! If we made it this far, reduce the actuation amount and try again.
           actuationAmount = 0.5_wp * actuationAmount

           ! Check actuation tolerance.
           if (actuationAmount < minimumTolerance) done = .true.

        end do

        ! Update the baseline values and compute a new sensitivity gradient.
        if (.not.done .and. controlIteration < nIterations) then
           do i = 1, controller%nParameters
              controller%baselineValue(i) = controller%baselineValue(i) +                    &
                   real(region%states(1)%gradientDirection, wp) * actuationAmount *          &
                   individualSensitivities(i)
           end do
           individualSensitivities = this%runAdjoint(region, controlIteration = nAdjoint)
           nAdjoint = nAdjoint + 1
           costSensitivity = sum(individualSensitivities**2)
           do i = 1, size(region%grids)
              region%states(i)%controlGradient = individualSensitivities
           end do
           actuationAmount = baselineActuationAmount
           foundNewMinimum = .false.
        end if

     end do

     ! Dump optimizatiom summary and close.
     if (procRank == 0) then
        write(fileUnit, *) ''
        write(fileUnit, '(A28,I4)') 'Number of forward runs:', nForward
        write(fileUnit, '(A28,I4)') 'Number of adjoint runs:', nAdjoint
        write(fileUnit, '(A28,L4)') 'Minimization found:', done
        close(fileUnit)
     end if

  case default
     write(message, '(A)') "Unknown optimization type!"
     call gracefulExit(region%comm, message)
  end select

  ! Clean up.
  SAFE_DEALLOCATE(parameters)
  SAFE_DEALLOCATE(individualSensitivities)

end subroutine findOptimalForcing
