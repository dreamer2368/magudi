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
    use State_enum, only : QOI_FORWARD_STATE, QOI_ADJOINT_STATE, QOI_RIGHT_HAND_SIDE    !SeungWhan: added rhs
    use Region_enum, only : FORWARD, ADJOINT, LINEARIZED

    ! <<< Internal modules >>>
    use RegionImpl, only : computeBulkQuantities
    use ErrorHandler, only : writeAndFlush
    use InputHelper, only : getOption, getRequiredOption

    ! <<< Arguments >>>
    class(t_Solver) :: this
    class(t_Region) :: region
    integer, intent(in) :: mode, timestep, startTimestep
    real(SCALAR_KIND), intent(in) :: time, instantaneousFunctional

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: nDimensions
    real(SCALAR_KIND) :: timeStepSize, cfl, residuals(3)
    real(SCALAR_KIND), allocatable :: bulkConservedVariables(:)
    character(len = STRING_LENGTH) :: str, str_, filename, strFormat
    class(t_Controller), pointer :: controller => null()
    class(t_Functional), pointer :: functional => null()
    logical :: adjointRestart, append
    integer :: accumulatedTimesteps = -1

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

       select case (mode)
       case (FORWARD)
          if (region%simulationFlags%enableFunctional) then
             write(str_, '(A,E13.6)') ", cost = ", instantaneousFunctional
             str = trim(str) // trim(str_)
          end if
       case (ADJOINT)
          write(str_, '(A,E13.6)') ", gradient = ", instantaneousFunctional
          str = trim(str) // trim(str_)
       end select

       call writeAndFlush(region%comm, output_unit, str)

       if (region%simulationFlags%checkConservation) then
         allocate(bulkConservedVariables(nDimensions+2))
         call computeBulkQuantities(region,bulkConservedVariables)
         write(strFormat,'(A,I1,A)') "(A,", nDimensions+2, "(X,E13.6),A)"
         write(str, strFormat) "Conservation status: (", bulkConservedVariables, ")"
         call writeAndFlush(region%comm, output_unit, str)
         SAFE_DEALLOCATE(bulkConservedVariables)
       end if

       if (region%outputOn) then

          select case (mode)

          case (FORWARD)
             if (region%simulationFlags%enableFunctional) then
                call this%functionalFactory%connect(functional)
                assert(associated(functional))
                call functional%writeToFile(region%comm, trim(this%outputPrefix) //             &
                     ".cost_functional.txt", timestep, time,                                    &
                     timestep - startTimestep > this%reportInterval)
             end if

          case (ADJOINT)
             call this%controllerFactory%connect(controller)
             assert(associated(controller))

             append = (startTimestep - timestep > this%reportInterval)
             if (.not. append ) then
                adjointRestart = getOption("enable_adjoint_restart", .false.)
                if (adjointRestart) then ! This is relative timesteps: timestep at initial condition must be added.
                    call getRequiredOption("adjoint_restart/accumulated_timesteps",accumulatedTimesteps)
                    assert( accumulatedTimesteps.ge.0 )
                end if
                append = adjointRestart .and. (accumulatedTimesteps>0)
             end if
             call controller%writeSensitivityToFile(region%comm, trim(this%outputPrefix) //  &
                  ".cost_sensitivity.txt", timestep, time, append)

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
       case (LINEARIZED)
          write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", timestep, ".linearized.q"
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

  function checkSolutionLimits(region, mode, outputPrefix) result(solutionCrashes)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env, only : output_unit

    ! <<< Derived types >>>
    use State_mod, only : t_State
    use Region_mod, only : t_Region

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE
    use Region_enum, only : FORWARD

    ! <<< Internal modules >>>
    use ErrorHandler, only : gracefulExit, writeAndFlush

    ! <<< Arguments >>>
    class(t_Region) :: region
    integer, intent(in) :: mode
    character(len = *), intent(in) :: outputPrefix

    logical :: solutionCrashes

    ! <<< Local variables >>>
    integer :: i, iGlobal, jGlobal, kGlobal, rankReportingError, procRank, ierror
    character(len = STRING_LENGTH) :: message
    SCALAR_TYPE :: fOutsideRange

    solutionCrashes = .false.

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
       call writeAndFlush(region%comm, output_unit, message)
       solutionCrashes = .true.

    end if

  end function checkSolutionLimits

  subroutine loadInitialCondition(this, region, mode, restartFilename)
    ! <<< External modules >>>
    use iso_fortran_env, only : output_unit

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use Region_mod, only : t_Region
    use Solver_mod, only : t_Solver
    use Functional_mod, only : t_Functional
    use CostTargetPatch_mod, only : t_CostTargetPatch

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE, QOI_ADJOINT_STATE
    use Region_enum, only : FORWARD, ADJOINT, LINEARIZED

    ! <<< Internal modules >>>
    use InputHelper, only : getOption, getRequiredOption
    use ErrorHandler, only : writeAndFlush, gracefulExit

    implicit none

    ! <<< Arguments >>>
    class(t_Solver) :: this
    class(t_Region) :: region
    integer, intent(in) :: mode
    character(len = *), intent(in), optional :: restartFilename

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    character(len = STRING_LENGTH) :: filename, message
    integer :: i, j
    real(wp) :: timeStepSize
    class(t_Functional), pointer :: functional => null()
    class(t_Patch), pointer :: patch => null()
    logical :: noAdjointForcing

    select case (mode)

    case (FORWARD) !... initialize conserved variables.

       !SeungWhan: unclear use of predictionOnly. revisit later
       if (present(restartFilename)) then ! .and. region%simulationFlags%predictionOnly) then
          call region%loadData(QOI_FORWARD_STATE, restartFilename)
       else if (region%simulationFlags%useTargetState) then
          filename = getOption("initial_condition_file", "")
          if (len_trim(filename) == 0) then
             write(message, '(A)') "Initial condition not specified. Using target state.."
             call writeAndFlush(region%comm, output_unit, message)
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

      if (present(restartFilename)) then
        call region%loadData(QOI_ADJOINT_STATE, restartFilename)
      else
        do i = 1, size(region%states)
           region%states(i)%adjointVariables = 0.0_wp
        end do
      end if

       if (allocated(region%patchFactories) .and.                                            &
            .not. region%simulationFlags%steadyStateSimulation .and.                         &
            .not. region%simulationFlags%useContinuousAdjoint) then

          timeStepSize = region%getTimeStepSize()
          region%states(:)%adjointForcingFactor = - timeStepSize / 6.0_wp !... RK4 only.

          ! Connect to the previously allocated functional.
          call this%functionalFactory%connect(functional)
          assert(associated(functional))
          noAdjointForcing = .not. region%simulationFlags%adjointForcingSwitch
          call functional%updateAdjointForcing(region,noAdjointForcing) !...SeungWhan:obviously not final step

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

          ! The last step of adjoint run will not include adjoint forcing term.
          ! If the adjoint continues from the previous adjoint run,
          ! then the corresponding adjoint forcing term is included here.
          ! If this is the initial adjoint step,
          ! then the corresponding adjoint forcing term is the initial condition.
          do i = 1, size(region%states)
             region%states(i)%adjointVariables = region%states(i)%adjointVariables            &
                                                + region%states(i)%rightHandSide
          end do
          region%states(:)%adjointForcingFactor = 1.0_wp !... restore

       end if

     case (LINEARIZED) !... initialize linearized variables. use adjoint variable.

       if (present(restartFilename)) then
         call region%loadData(QOI_ADJOINT_STATE, restartFilename) !... initialize from file.
       else
         write(message,'(A)') "Initial condition file for linearized run is not given!"
         call gracefulExit(region%comm,message)
       end if

    end select

  end subroutine loadInitialCondition

end module SolverImpl

subroutine setupSolver(this, region, restartFilename, outputPrefix)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver
  use Patch_mod, only : t_Patch
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch
  use Controller_mod, only : t_Controller
  use Functional_mod, only : t_Functional
  use TimeIntegrator_mod, only : t_TimeIntegrator

  ! <<< Enumerations >>>
  use Grid_enum, only : QOI_CONTROL_MOLLIFIER, QOI_TARGET_MOLLIFIER
  use State_enum, only : QOI_FORWARD_STATE, QOI_TARGET_STATE, QOI_ADJOINT_STATE
  use BlockInterfacePatch_enum, only : METRICS

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use Patch_factory, only : computeSpongeStrengths, updatePatchFactories
  use InterfaceHelper, only : checkFunctionContinuityAtInterfaces, exchangeInterfaceData

  implicit none

  ! <<< Arguments >>>
  class(t_Solver) :: this
  class(t_Region) :: region
  character(len = *), intent(in), optional :: restartFilename, outputPrefix

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: filename
  integer :: i, j
  class(t_Patch), pointer :: patch => null()
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

  this%probeInterval = getOption("probe_interval", -1)
  if (this%probeInterval == 0) this%probeInterval = -1

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

  if (region%simulationFlags%enableController) then

     ! Initialize control mollifier.
     filename = getOption("control_mollifier_file", "")
     if (len_trim(filename) == 0) then
        do i = 1, size(region%grids)
           region%grids(i)%controlMollifier = 1.0_wp
        end do
     else
        call region%loadData(QOI_CONTROL_MOLLIFIER, filename)
     end if

  end if

  if (region%simulationFlags%enableFunctional) then

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

  ! Exchange metrics data at block interfaces.
  if (allocated(region%patchFactories)) then
     do i = 1, size(region%patchFactories)
        call region%patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        do j = 1, size(region%states)
           if (patch%gridIndex /= region%grids(j)%index) cycle
           select type (patch)
           class is (t_BlockInterfacePatch)
              call patch%collectInterfaceData(METRICS, region%simulationFlags,                    &
                   region%solverOptions, region%grids(j), region%states(j))
           end select
        end do
     end do
  end if

  call exchangeInterfaceData(region)

  ! Disperse received metrics data at block interfaces.
  if (allocated(region%patchFactories)) then
     do i = 1, size(region%patchFactories)
        call region%patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        do j = 1, size(region%states)
           if (patch%gridIndex /= region%grids(j)%index) cycle
           select type (patch)
           class is (t_BlockInterfacePatch)
              call patch%disperseInterfaceData(METRICS, region%simulationFlags,                   &
                   region%solverOptions)
           end select
        end do
     end do
  end if

  if (region%simulationFlags%enableController) then

     call this%controllerFactory%connect(controller,                                         &
          trim(region%solverOptions%controllerType))
     assert(associated(controller))
     call controller%setup(region)

  end if

  if (region%simulationFlags%enableFunctional) then

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

function runForward(this, region, restartFilename) result(costFunctional)

  ! <<< External modules >>>
  use iso_fortran_env, only : output_unit

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
  class(t_Solver) :: this
  class(t_Region) :: region
  character(len = *), intent(in), optional :: restartFilename

  ! <<< Result >>>
  SCALAR_TYPE :: costFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: filename, message
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()
  class(t_Controller), pointer :: controller => null()
  class(t_Functional), pointer :: functional => null()
  integer :: i, j, timestep, startTimestep
  real(wp) :: time, startTime, timeStepSize
  SCALAR_TYPE :: instantaneousCostFunctional
  logical :: controllerSwitch = .false., solutionCrashes = .false.

  call startTiming("runForward")

  costFunctional = 0.0_wp

  ! Connect to the previously allocated time integrator.
  call this%timeIntegratorFactory%connect(timeIntegrator)
  assert(associated(timeIntegrator))

  ! Connect to the previously allocated controller.
  if (region%simulationFlags%enableController) then
     call this%controllerFactory%connect(controller)
     assert(associated(controller))
  end if

  ! Connect to the previously allocated functional.
  if (region%simulationFlags%enableFunctional) then
     call this%functionalFactory%connect(functional)
     assert(associated(functional))
     functional%runningTimeQuadrature = 0.0_wp
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

  if (region%simulationFlags%enableBodyForce) then
    call getRequiredOption("body_force/initial_momentum", region%initialXmomentum)
    region%oneOverVolume = computeRegionIntegral(region)
    region%initialXmomentum = region%initialXmomentum * region%oneOverVolume
    region%oneOverVolume = 1.0_wp / region%oneOverVolume
    region%momentumLossPerVolume = 0.0_wp
  end if

  startTimestep = region%timestep
  startTime = region%states(1)%time

  ! Save the initial condition if it was not specified as a restart file.
  if (.not. present(restartFilename)) then
     write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", startTimestep, ".q"
     call region%saveData(QOI_FORWARD_STATE, filename)
  end if

  ! Call controller hooks before time marching starts. SeungWhan: changed
  ! duration, onsetTime.
  if (region%simulationFlags%enableController) then
    controllerSwitch = controller%controllerSwitch
    write(message,'(A,L4,A,F16.8,A,F16.8)') 'Control Forcing Flag: ',                          &
                                            controller%controllerSwitch
    call writeAndFlush(region%comm, output_unit, message)
    if (controller%controllerSwitch) then
       controller%onsetTime = startTime
       controller%duration = this%nTimesteps * region%solverOptions%timeStepSize
       call controller%hookBeforeTimemarch(region, FORWARD)
    end if
  end if

  ! Reset probes.
  if (this%probeInterval > 0) call region%resetProbes()

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
        if (region%simulationFlags%enableSolutionLimits) then
          if ( checkSolutionLimits(region, FORWARD, this%outputPrefix) ) then
            costFunctional = HUGE(0.0_wp)
            return
          end if
        end if

        ! Update control forcing.
        !SeungWhan: flush out previous control forcing
        if (controllerSwitch) then
          call controller%cleanupForcing(region)
          call controller%updateForcing(region)
        end if

        ! Take a single sub-step using the time integrator.
        call timeIntegrator%substepForward(region, time, timeStepSize, timestep, i)

        do j = 1, size(region%states) !... update state
           call region%states(j)%update(region%grids(j), region%simulationFlags,             &
                region%solverOptions)
        end do

        ! Update the cost functional.
        if (region%simulationFlags%enableFunctional) then
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

     ! Report simulation progess.
     call showProgress(this, region, FORWARD, startTimestep, timestep,                       &
          time, instantaneousCostFunctional)

     ! Save solution on probe patches.
     if (this%probeInterval > 0 .and. mod(timestep, max(1, this%probeInterval)) == 0)        &
          call region%saveProbeData(FORWARD)

     ! Stop if this is a steady-state simulation and solution has converged.
     if (this%residualManager%hasSimulationConverged) exit

     ! Filter solution if required.
     if (region%simulationFlags%filterOn) then
        do j = 1, size(region%grids)
           call region%grids(j)%applyFilter(region%states(j)%conservedVariables, timestep)
        end do
     end if

  end do !... timestep = startTimestep + 1, startTimestep + this%nTimesteps

  ! Finish writing remaining data gathered on probes.
  if (this%probeInterval > 0) call region%saveProbeData(FORWARD, finish = .true.)

  ! Call controller hooks after time marching ends.
  if (controllerSwitch) call controller%hookAfterTimemarch(region, FORWARD)

  call this%residualManager%cleanup()

  if (region%simulationFlags%computeTimeAverage) then
     do i = 1, size(region%states)
        region%states(i)%timeAverage = region%states(i)%timeAverage / (time - startTime)
     end do
     call region%saveData(QOI_TIME_AVERAGED_STATE, trim(this%outputPrefix) // ".mean.q")
  end if

  if (region%simulationFlags%enableFunctional) then
       costFunctional = functional%runningTimeQuadrature
       write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))') 'Forward run: cost functional = ', &
                                                               costFunctional
       call writeAndFlush(region%comm, output_unit, message)
  end if

  call endTiming("runForward")

end function runForward

function runAdjoint(this, region) result(costSensitivity)

  ! <<< External modules >>>
  use iso_fortran_env, only : output_unit

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
  use RegionImpl, only : computeRegionIntegral

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming
  use ErrorHandler, only : writeAndFlush
  use InputHelper, only : getOption, getRequiredOption

  ! <<< SeungWhan: debug >>>
  use, intrinsic :: iso_fortran_env, only : output_unit

  implicit none

  ! <<< Arguments >>>
  class(t_Solver) :: this
  class(t_Region) :: region

  ! <<< Result >>> - This value is currently not meaningful. left for future purpose.
  SCALAR_TYPE :: costSensitivity

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: filename, message
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()
  class(t_Controller), pointer :: controller => null()
  class(t_Functional), pointer :: functional => null()
  type(t_ReverseMigratorFactory) :: reverseMigratorFactory
  class(t_ReverseMigrator), pointer :: reverseMigrator => null()
  integer :: i, j, timestep, startTimestep, timemarchDirection
  real(SCALAR_KIND) :: time, startTime, timeStepSize
  logical :: adjointRestart, nonzeroAdjointInitialCondition,                    &
             IS_FINAL_STEP, noAdjointForcing
  integer :: accumulatedNTimesteps, intermediateEndTimestep
  integer :: procRank, ierror
  SCALAR_TYPE :: instantaneousCostSensitivity

  assert(region%simulationFlags%enableController)
  assert(region%simulationFlags%enableFunctional)
  assert(region%simulationFlags%enableAdjoint)

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

  ! Adjoint initial condition option.
  nonzeroAdjointInitialCondition =                                                           &
                  getOption("adjoint_nonzero_initial_condition",.false.)

  ! Load adjoint restart options.
  adjointRestart = getOption("enable_adjoint_restart", .false.)
  accumulatedNTimesteps = -1
  intermediateEndTimestep = -1
  if (adjointRestart) then ! This is relative timesteps: timestep at initial condition must be added.
    call getRequiredOption("adjoint_restart/accumulated_timesteps",accumulatedNTimesteps)
    call getRequiredOption("adjoint_restart/intermediate_end_timestep",intermediateEndTimestep)
    assert( (accumulatedNTimesteps.ge.0).and.(intermediateEndTimestep.ge.0) )

    nonzeroAdjointInitialCondition = (accumulatedNTimesteps>0) .or.                          &
                                                            nonzeroAdjointInitialCondition
  end if

  ! Load the initial condition: this does not require restart filename even for adjoint restart.
  call loadInitialCondition(this, region, FORWARD) !... for control horizon end timestep.
  controller%onsetTime = region%states(1)%time
  controller%duration = this%nTimesteps * region%solverOptions%timeStepSize

  if (region%simulationFlags%enableBodyForce) then
    call getRequiredOption("body_force/initial_momentum", region%initialXmomentum)
    region%oneOverVolume = computeRegionIntegral(region)
    region%initialXmomentum = region%initialXmomentum * region%oneOverVolume
    region%oneOverVolume = 1.0_wp / region%oneOverVolume

    region%momentumLossPerVolume = 0.0_wp
    region%adjointMomentumLossPerVolume = 0.0_wp
  end if

  ! Load the adjoint coefficients corresponding to the end of the control time horizon.
  if (region%simulationFlags%steadyStateSimulation) then
     write(filename, '(2A)') trim(this%outputPrefix), ".steady_state.q"
  elseif (adjointRestart) then
    write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-",                            &
         region%timestep + intermediateEndTimestep + this%nTimesteps, ".q"
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

  region%states(:)%time = startTime
  region%states(:)%timeProgressive = startTime

  ! Adjoint initial condition (if specified).
  write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", region%timestep, ".adjoint.q"
  if (nonzeroAdjointInitialCondition) then
    call loadInitialCondition(this, region, ADJOINT, filename)
  else
    call loadInitialCondition(this, region, ADJOINT)
    call region%saveData(QOI_ADJOINT_STATE, filename)
  end if

  ! Call controller hooks before time marching starts.
  if (adjointRestart) then
    call controller%hookBeforeTimemarch(region, ADJOINT, accumulatedNTimesteps)
  else
    call controller%hookBeforeTimemarch(region, ADJOINT)
  end if
  !!!SeungWhan: need additional execution with FORWARD, in case of non-zero control forcing.
  if (controller%controllerSwitch) then
    controller%duration = this%nTimesteps * region%solverOptions%timeStepSize
    controller%onsetTime = startTime - controller%duration
    call controller%hookBeforeTimemarch(region, FORWARD)
  end if

  ! Reset probes.
  if (this%probeInterval > 0) call region%resetProbes()

  ! SeungWhan: time must be synchronized between functional in forward and adjoint forcing.
  timeStepSize = region%getTimeStepSize()
  time = startTime - 0.5_wp * timeStepSize
  region%states(:)%time = time

  do timestep = startTimestep + sign(1, timemarchDirection),                                 &
       startTimestep + sign(this%nTimesteps, timemarchDirection), timemarchDirection

     region%timestep = timestep
     timeStepSize = region%getTimeStepSize()

     do i = timeIntegrator%nStages, 1, -1

        ! Load adjoint coefficients.
        if (.not. region%simulationFlags%steadyStateSimulation) then !... unsteady simulation.
           if (i == 1) then
              call reverseMigrator%migrateTo(region, controller, timeIntegrator,             &
                   timestep, timeIntegrator%nStages)
           else
              call reverseMigrator%migrateTo(region, controller, timeIntegrator,             &
                   timestep + 1, i - 1)
           end if
        end if

        call startTiming("controllerUpdateGradient")

        ! Update gradient.
        call controller%updateGradient(region)

        call endTiming("controllerUpdateGradient")

        ! Update cost sensitivity.
        instantaneousCostSensitivity = controller%computeSensitivity(region)
        controller%runningTimeQuadrature = controller%runningTimeQuadrature +                &
             timeIntegrator%norm(i) * timeStepSize * instantaneousCostSensitivity

        ! Check the switch for adjoint forcing.
        if( (timestep.eq.startTimestep+sign(this%nTimesteps,timemarchDirection))             &
            .and. (i.eq.1) ) then
            IS_FINAL_STEP = .true.
        else
            IS_FINAL_STEP = .false.
        end if
        noAdjointForcing = (.not. region%simulationFlags%adjointForcingSwitch)               &
                            .or. (IS_FINAL_STEP)
        ! Update adjoint forcing on cost target patches.
        call functional%updateAdjointForcing(region,noAdjointForcing)

        ! Take a single adjoint sub-step using the time integrator.
        call timeIntegrator%substepAdjoint(region, time, timeStepSize, timestep, i)

        ! TODO: how to enforce limits on adjoint variables... check for NaN?
        if (region%simulationFlags%enableSolutionLimits) then
          if ( checkSolutionLimits(region, ADJOINT, this%outputPrefix) ) then
            costSensitivity = HUGE(0.0_wp)
            return
          end if
        end if

     end do

     ! Report simulation progress.
     call showProgress(this, region, ADJOINT, startTimestep, timestep,                       &
          time, instantaneousCostSensitivity)

     ! Save solution on probe patches.
     if (this%probeInterval > 0 .and. mod(timestep, max(1, this%probeInterval)) == 0)        &
          call region%saveProbeData(ADJOINT)

     ! Stop if this is a steady-state simulation and solution has converged.
     if (this%residualManager%hasSimulationConverged) exit

     ! Filter solution if required.
     if (region%simulationFlags%filterOn) then
        do j = 1, size(region%grids)
           call region%grids(j)%applyFilter(region%states(j)%adjointVariables, timestep)
        end do
     end if

  end do !... timestep = startTimestep + sign(1, timemarchDirection), ...

  call MPI_Comm_rank(region%comm, procRank, ierror)

  ! Finish writing remaining data gathered on probes.
  if (this%probeInterval > 0) call region%saveProbeData(ADJOINT, finish = .true.)

  ! Call controller hooks after time marching ends.
  if (controller%controllerSwitch) call controller%hookAfterTimemarch(region, FORWARD)
  call controller%hookAfterTimemarch(region, ADJOINT)

  call MPI_Barrier(region%comm,ierror)

  call this%residualManager%cleanup()
  call reverseMigratorFactory%cleanup()

  costSensitivity = controller%runningTimeQuadrature

  write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))') 'Adjoint run: cost sensitivity = ', &
                                                          costSensitivity
  call writeAndFlush(region%comm, output_unit, message)

  call endTiming("runAdjoint")

end function runAdjoint

function runLinearized(this, region) result(costFunctional)

  ! <<< External modules >>>
  use iso_fortran_env, only : output_unit

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
  class(t_Solver) :: this
  class(t_Region) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: costFunctional !TODO: compute delta functional.

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: filename, message
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()
  class(t_Controller), pointer :: controller => null()
  class(t_Functional), pointer :: functional => null()
  type(t_ReverseMigratorFactory) :: reverseMigratorFactory
  class(t_ReverseMigrator), pointer :: reverseMigrator => null()
  integer :: i, j, timestep, startTimestep
  real(wp) :: time, startTime, timeStepSize
  SCALAR_TYPE :: instantaneousCostFunctional
  logical :: controllerSwitch = .false., solutionCrashes = .false.

  call startTiming("runLinearized")

  costFunctional = 0.0_wp

  ! Connect to the previously allocated time integrator.
  call this%timeIntegratorFactory%connect(timeIntegrator)
  assert(associated(timeIntegrator))

  ! Connect to the previously allocated controller.
  if (region%simulationFlags%enableController) then
     call this%controllerFactory%connect(controller)
     assert(associated(controller))
  end if

  ! Connect to the previously allocated functional.
  if (region%simulationFlags%enableFunctional) then
     call this%functionalFactory%connect(functional)
     assert(associated(functional))
     functional%runningTimeQuadrature = 0.0_wp
  end if

  ! Load the initial condition.
  call loadInitialCondition(this, region, FORWARD)                              ! for baseline simulation.
  write(filename, '(2A)') trim(this%outputPrefix), ".ic.linearized.q"
  call loadInitialCondition(this, region, LINEARIZED, trim(filename))           ! for linearized simulation.

  if (region%simulationFlags%enableBodyForce) then
    call getRequiredOption("body_force/initial_momentum", region%initialXmomentum)
    region%oneOverVolume = computeRegionIntegral(region)
    region%initialXmomentum = region%initialXmomentum * region%oneOverVolume
    region%oneOverVolume = 1.0_wp / region%oneOverVolume

    region%momentumLossPerVolume = 0.0_wp
    region%adjointMomentumLossPerVolume = 0.0_wp
  end if

  startTimestep = region%timestep
  startTime = region%states(1)%time

  ! Connect to the previously allocated reverse migrator.
  call reverseMigratorFactory%connect(reverseMigrator,                                       &
       region%solverOptions%checkpointingScheme)
  assert(associated(reverseMigrator))

  ! Setup the revese-time migrator.
  call reverseMigrator%setup(region, timeIntegrator, this%outputPrefix,                      &
  startTimestep, startTimestep + this%nTimesteps, this%saveInterval,                         &
  this%saveInterval * timeIntegrator%nStages)

  ! Save the initial condition if it was not specified as a restart file.
  write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", startTimestep, ".q"
  call region%saveData(QOI_FORWARD_STATE, filename)
  write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", startTimestep, ".linearized.q"
  call region%saveData(QOI_ADJOINT_STATE, filename)

  ! Call controller hooks before time marching starts. SeungWhan: changed
  ! duration, onsetTime.
  if (region%simulationFlags%enableController) then
    controllerSwitch = controller%controllerSwitch
    write(message,'(A,L4,A,F16.8,A,F16.8)') 'Control Forcing Flag: ',                          &
                                            controller%controllerSwitch
    call writeAndFlush(region%comm, output_unit, message)
    controller%onsetTime = startTime
    controller%duration = this%nTimesteps * region%solverOptions%timeStepSize
    if (controller%controllerSwitch) then
       call controller%hookBeforeTimemarch(region, FORWARD)
    end if
    call controller%hookBeforeTimemarch(region, LINEARIZED) ! TODO: test linearized mode
  end if

  ! Reset probes.
  if (this%probeInterval > 0) call region%resetProbes()

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

       ! Load forward baseline state.
       if (i == 1) then
          call reverseMigrator%migrateTo(region, controller, timeIntegrator,             &
               timestep - 1, timeIntegrator%nStages)
       else
          call reverseMigrator%migrateTo(region, controller, timeIntegrator,             &
               timestep, i - 1)
       end if

        ! linearized step first
        call controller%cleanupDeltaForcing(region)
        if(controllerSwitch) call controller%updateDeltaForcing(region)

        call timeIntegrator%substepLinearized(region, time, timeStepSize, timestep, i)

        !TODO: Update the delta functional.
        if (region%simulationFlags%enableFunctional) then
           instantaneousCostFunctional = functional%compute(region)
           functional%runningTimeQuadrature = functional%runningTimeQuadrature +             &
                timeIntegrator%norm(i) * timeStepSize * instantaneousCostFunctional
        end if

        ! ! Update the time average.
        ! if (region%simulationFlags%computeTimeAverage) then
        !    do j = 1, size(region%states)
        !       region%states(j)%timeAverage = region%states(j)%timeAverage +                  &
        !            timeIntegrator%norm(i) * timeStepSize *                                   &
        !            region%states(j)%conservedVariables
        !    end do
        ! end if

     end do !... i = 1, timeIntegrator%nStages

     ! Report simulation progess.
     call showProgress(this, region, LINEARIZED, startTimestep, timestep,                       &
          time, instantaneousCostFunctional)

     ! Save solution on probe patches.
     if (this%probeInterval > 0 .and. mod(timestep, max(1, this%probeInterval)) == 0)        &
          call region%saveProbeData(FORWARD)

     ! Filter solution if required.
     if (region%simulationFlags%filterOn) then
        do j = 1, size(region%grids)
           call region%grids(j)%applyFilter(region%states(j)%conservedVariables, timestep)
        end do
     end if

  end do !... timestep = startTimestep + 1, startTimestep + this%nTimesteps

  ! Finish writing remaining data gathered on probes.
  if (this%probeInterval > 0) call region%saveProbeData(FORWARD, finish = .true.)

  ! Call controller hooks after time marching ends.
  if (controllerSwitch) then
    call controller%hookAfterTimemarch(region, FORWARD)
    call controller%hookAfterTimemarch(region, LINEARIZED)
  end if

  ! if (region%simulationFlags%computeTimeAverage) then
  !    do i = 1, size(region%states)
  !       region%states(i)%timeAverage = region%states(i)%timeAverage / (time - startTime)
  !    end do
  !    call region%saveData(QOI_TIME_AVERAGED_STATE, trim(this%outputPrefix) // ".mean.q")
  ! end if

  call reverseMigratorFactory%cleanup()

  ! TODO: compute delta functional.
  if (region%simulationFlags%enableFunctional) then
       costFunctional = functional%runningTimeQuadrature
       write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))') 'Forward run: cost functional = ', &
                                                               costFunctional
       call writeAndFlush(region%comm, output_unit, message)
  end if

  call endTiming("runLinearized")

end function runLinearized

! checkGradientAccuracy is obsolete!!
subroutine checkGradientAccuracy(this, region)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver
  use ActuatorPatch_mod, only : t_ActuatorPatch

  ! <<< Enumerations >>>
  use State_enum, only : QOI_FORWARD_STATE

  ! <<< Internal modules >>>
  use ControlSpaceAdvancer, only : ZAXPY
  use InputHelper, only : getFreeUnit, getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_Solver) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nIterations, restartIteration, fileUnit, iostat, procRank, ierror
  integer :: numberOfActuatorPatches
  class(t_Patch), pointer :: patch => null()
  character(len = STRING_LENGTH) :: filename, message
  character(len = STRING_LENGTH) :: gradientFilename, controlForcingFilename
  real(wp) :: actuationAmount, baselineCostFunctional, costFunctional, costSensitivity,      &
       initialActuationAmount, geometricGrowthFactor, gradientError, dummyValue

  write(message, '(A)') "Subroutine checkGradientAccuracy is obsolete. Use utils/python/checkGradientAccuracy.py."
  call gracefulExit(region%comm, message)

  ! call getRequiredOption("number_of_control_iterations", nIterations)
  ! if (nIterations < 0) then
  !    write(message, '(A)') "Number of control iterations must be a non-negative number!"
  !    call gracefulExit(region%comm, message)
  ! end if
  !
  ! restartIteration = getOption("restart_control_iteration", 0)
  ! restartIteration = max(restartIteration, 0)
  !
  ! if (restartIteration > 0 .and. .not. region%simulationFlags%isBaselineAvailable) then
  !    write(message, '(A)') "Can't restart with controlled prediction without baseline!"
  !    call gracefulExit(region%comm, message)
  ! end if
  !
  ! call MPI_Comm_rank(region%comm, procRank, ierror)
  !
  ! !SeungWhan: add option for gradient_error_filename
  ! filename = getOption("gradient_error_file", "")
  ! if (len_trim(filename) == 0) then
  !    write(filename, '(2A)') trim(this%outputPrefix), ".gradient_error.txt"
  ! end if
  !
  ! if (procRank == 0) then
  !    if (restartIteration == 0 .and. .not. region%simulationFlags%isBaselineAvailable) then
  !       open(unit = getFreeUnit(fileUnit), file = trim(filename), action = 'write',          &
  !            status = 'unknown', iostat = iostat)
  !    else
  !       open(unit = getFreeUnit(fileUnit), file = trim(filename), action = 'readwrite',      &
  !            status = 'old', position = 'append', iostat = iostat)
  !    end if
  ! end if
  !
  ! call MPI_Bcast(iostat, 1, MPI_INTEGER, 0, region%comm, ierror)
  ! if (iostat /= 0) then
  !    write(message, "(2A)") trim(filename), ": Failed to open file for writing!"
  !    call gracefulExit(region%comm, message)
  ! end if
  !
  ! if (nIterations > 0) then
  !    call getRequiredOption("initial_actuation_amount", initialActuationAmount)
  !    if (nIterations > 1) then
  !       call getRequiredOption("actuation_amount_geometric_growth", geometricGrowthFactor)
  !    else
  !       geometricGrowthFactor = getOption("actuation_amount_geometric_growth", 1.0_wp)
  !    end if
  ! end if
  !
  ! ! Find (or load from file) the cost functional for the baseline prediction.
  ! if (region%simulationFlags%isBaselineAvailable) then
  !    if (procRank == 0)                                                                      &
  !         read(fileUnit, *, iostat = iostat) i, actuationAmount,                             &
  !         baselineCostFunctional, costSensitivity, gradientError
  !    call MPI_Bcast(iostat, 1, MPI_INTEGER, 0, region%comm, ierror)
  !    if (iostat /= 0) then
  !       write(message, "(2A)") trim(filename),                                               &
  !            ": Failed to read baseline cost functional from file!"
  !       call gracefulExit(region%comm, message)
  !    end if
  !    call MPI_Bcast(baselineCostFunctional, 1, REAL_TYPE_MPI, 0, region%comm, ierror)
  ! else
  !    baselineCostFunctional = this%runForward(region)
  ! end if
  !
  ! ! Find the sensitivity gradient (this is the only time the adjoint simulation will be run).
  ! if (restartIteration == 0) then
  !    costSensitivity = this%runAdjoint(region)
  ! else
  !    call MPI_Bcast(costSensitivity, 1, REAL_TYPE_MPI, 0, region%comm, ierror)
  ! end if
  !
  ! if (procRank == 0 .and. .not. region%simulationFlags%isBaselineAvailable)                  &
  !      write(fileUnit, '(I4,4(1X,SP,' // SCALAR_FORMAT // '))') 0, 0.0_wp,                   &
  !      baselineCostFunctional, costSensitivity, 0.0_wp
  !
  ! if (nIterations == 0) return
  !
  ! ! Find filenames for gradient and control forcing
  ! numberOfActuatorPatches = 0
  ! do i = 1, size(region%patchFactories)
  !   call region%patchFactories(i)%connect(patch)
  !   print *, procRank, i, trim(region%patchData(i)%patchType)
  !   if (.not. associated(patch)) cycle
  !   select type (patch)
  !   class is (t_ActuatorPatch)
  !     if (patch%comm == MPI_COMM_NULL) cycle
  !
  !     numberOfActuatorPatches = numberOfActuatorPatches + 1
  !     gradientFilename = trim(patch%gradientFilename)
  !     controlForcingFilename = trim(patch%controlForcingFilename)
  !   end select
  ! end do
  !
  ! ! Turn off output for controlled predictions.
  ! region%outputOn = .false.
  !
  ! if (restartIteration == 0) restartIteration = restartIteration + 1
  !
  ! do i = 1, restartIteration - 1
  !    if (procRank == 0)                                                                      &
  !         read(fileUnit, *, iostat = iostat) j, actuationAmount, costFunctional,             &
  !         dummyValue, gradientError
  !    call MPI_Bcast(iostat, 1, MPI_INTEGER, 0, region%comm, ierror)
  !    if (iostat /= 0) then
  !       write(message, "(2A)") trim(filename),                                               &
  !            ": Cost functional history is too short for the specified restart iteration!"
  !       call gracefulExit(region%comm, message)
  !    end if
  ! end do
  !
  ! do i = restartIteration, restartIteration + nIterations - 1
  !    actuationAmount = initialActuationAmount * geometricGrowthFactor ** real(i - 1, wp)
  !    call ZAXPY(region%comm,trim(controlForcingFilename),-actuationAmount,trim(gradientFilename))
  !    costFunctional = this%runForward(region, actuationAmount = actuationAmount)
  !    gradientError = (costFunctional - baselineCostFunctional) / actuationAmount +           &
  !         costSensitivity
  !    if (procRank == 0) then
  !       write(fileUnit, '(I4,4(1X,SP,' // SCALAR_FORMAT // '))') i, actuationAmount,         &
  !            costFunctional, -(costFunctional - baselineCostFunctional) / actuationAmount,   &
  !            gradientError
  !       flush(fileUnit)
  !    end if
  ! end do
  !
  ! if (procRank == 0) close(fileUnit)

end subroutine checkGradientAccuracy
