#include "config.h"

module Solver_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  use Controller_mod, only : t_Controller
  use Functional_mod, only : t_Functional
  use TimeIntegrator_mod, only : t_TimeIntegrator
  use ResidualManager_mod, only : t_ResidualManager
  use ReverseMigrator_mod, only : t_ReverseMigrator

  implicit none
  private

  type, public :: t_Solver

     type(t_ResidualManager) :: residualManager
     class(t_TimeIntegrator), pointer :: timeIntegrator => null()
     class(t_ReverseMigrator), pointer :: reverseMigrator => null()
     class(t_Functional), pointer :: functional => null()
     class(t_Controller), pointer :: controller => null()
     integer :: nTimesteps, saveInterval, reportInterval, probeInterval, filterInterval
     character(len = STRING_LENGTH) :: outputPrefix

   contains

     procedure, pass :: setup
     procedure, pass :: cleanup
     procedure, pass :: runForward
     procedure, pass :: runAdjoint
     procedure, pass :: checkGradientAccuracy

  end type t_Solver

contains

  subroutine setup(this, timeIntegrator, functional, controller,                             &
       reverseMigrator, outputPrefix)

    ! <<< Internal modules >>>
    use InputHelper, only : getOption

    implicit none

    ! <<< Arguments >>>
    class(t_Solver) :: this
    class(t_TimeIntegrator), pointer, intent(in) :: timeIntegrator
    class(t_Functional), pointer, intent(in), optional :: functional
    class(t_Controller), pointer, intent(in), optional :: controller
    class(t_ReverseMigrator), pointer, intent(in), optional :: reverseMigrator
    character(len = *), intent(in), optional :: outputPrefix

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
    this%filterInterval = getOption("filter_interval", 1)

    this%timeIntegrator => timeIntegrator
    assert(associated(this%timeIntegrator))

    if (present(functional)) then
       this%functional => functional
       assert(associated(this%functional))
    end if

    if (present(controller)) then
       this%controller => controller
       assert(associated(this%controller))
    end if

    if (present(reverseMigrator)) then
       this%reverseMigrator => reverseMigrator
       assert(associated(this%reverseMigrator))
    end if

  end subroutine setup

  subroutine cleanup(this)

    implicit none

    ! <<< Arguments >>>
    class(t_Solver) :: this

    nullify(this%timeIntegrator)
    nullify(this%controller)
    nullify(this%functional)
    nullify(this%reverseMigrator)

  end subroutine cleanup

  function runForward(this, region, actuationAmount, restartFilename) result(functional)

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE, QOI_TIME_AVERAGED_STATE
    use SolverOptions_enum, only : FORWARD

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    implicit none

    ! <<< Arguments >>>
    class(t_Solver) :: this
    class(t_Region) :: region
    real(SCALAR_KIND), intent(in), optional :: actuationAmount
    character(len = *), intent(in), optional :: restartFilename

    ! <<< Result >>>
    real(SCALAR_KIND) :: functional

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, timestep, startTimestep
    real(wp) :: time, startTime, timeStepSize
    character(len = STRING_LENGTH) :: filename
    real(SCALAR_KIND) :: instantaneousFunctional

    call startTiming("Compute functional")

    functional = 0.0_wp

    region%states(:)%actuationAmount = 0.0_wp
    if (present(actuationAmount) .and. .not. region%simulationFlags%predictionOnly)          &
         region%states(:)%actuationAmount = actuationAmount

    if (.not. region%simulationFlags%predictionOnly)                                         &
         this%functional%runningTimeQuadrature = 0.0_wp

    ! Setup residual manager if this is a steady-state simulation.
    if (region%simulationFlags%steadyStateSimulation)                                        &
         call this%residualManager%setup("", region)

    ! Load the initial condition.
    if (present(restartFilename)) then
       call initializeState(this, region, FORWARD, restartFilename)
    else
       call initializeState(this, region, FORWARD)
    end if

    startTimestep = region%timestep
    startTime = region%states(1)%time

    ! Save the initial condition if it was not specified as a restart file.
    if (.not. present(restartFilename)) then
       write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", startTimestep, ".q"
       call region%saveData(QOI_FORWARD_STATE, filename)
    end if

    ! Call controller hooks before time marching starts.
    if (.not. region%simulationFlags%predictionOnly .and.                                    &
         abs(region%states(1)%actuationAmount) > 0.0_wp) then
       this%controller%onsetTime = startTime
       this%controller%duration = this%nTimesteps * region%solverOptions%timeStepSize
       call this%controller%hookBeforeTimemarch(region, FORWARD)
    end if

    ! Reset probes.
    if (this%probeInterval > 0) call region%resetProbes()

    time = startTime
    region%states(:)%time = time
    do i = 1, size(region%states) !... update state
       call region%states(i)%update(region%grids(i), region%simulationFlags,                 &
            region%solverOptions)
    end do

    do timestep = startTimestep + 1, startTimestep + this%nTimesteps

       region%timestep = timestep
       timeStepSize = region%computeTimeStepSize()

       do i = 1, this%timeIntegrator%nStages

          ! Check if physical quantities are within allowed limits.
          if (region%simulationFlags%enableSolutionLimits)                                   &
               call checkSolutionLimits(region, FORWARD, this%outputPrefix)

          ! Update control forcing.
          if (.not. region%simulationFlags%predictionOnly .and.                              &
               abs(region%states(1)%actuationAmount) > 0.0_wp)                               &
               call this%controller%update(region)

          ! Take a single sub-step using the time integrator.
          call this%timeIntegrator%substepForward(region, time, timeStepSize, timestep, i)

          do j = 1, size(region%states) !... update state
             call region%states(j)%update(region%grids(j), region%simulationFlags,           &
                  region%solverOptions)
          end do

          ! Update the cost functional.
          if (.not. region%simulationFlags%predictionOnly) then
             instantaneousFunctional = this%functional%compute(region)
             this%functional%runningTimeQuadrature = this%functional%runningTimeQuadrature + &
                  this%timeIntegrator%norm(i) * timeStepSize * instantaneousFunctional
          end if

          ! Update the time average.
          if (region%simulationFlags%computeTimeAverage) then
             do j = 1, size(region%states)
                region%states(j)%timeAverage = region%states(j)%timeAverage +                &
                     this%timeIntegrator%norm(i) * timeStepSize *                            &
                     region%states(j)%conservedVariables
             end do
          end if

       end do !... i = 1, timeIntegrator%nStages

       ! Book-keeping, including simulation progress, residual monitoring, I/O, etc.
       call bookKeeping(this, region, FORWARD, startTimestep,                                &
            timestep, time, instantaneousFunctional)

       ! Save solution on probe patches.
       if (this%probeInterval > 0 .and. mod(timestep, max(1, this%probeInterval)) == 0)      &
            call region%saveProbeData(FORWARD)

       ! Stop if this is a steady-state simulation and solution has converged.
       if (this%residualManager%hasSimulationConverged) exit

       ! Filter solution if required.
       if (region%simulationFlags%filterOn .and.                                             &
            mod(timestep, max(1, this%filterInterval)) == 0) then
          do j = 1, size(region%grids)
             call region%grids(j)%applyFilter(region%states(j)%conservedVariables, timestep)
          end do
       end if

    end do !... timestep = startTimestep + 1, startTimestep + this%nTimesteps

    ! Finish writing remaining data gathered on probes.
    if (this%probeInterval > 0) call region%saveProbeData(FORWARD, finish = .true.)

    ! Call controller hooks after time marching ends.
    if (.not. region%simulationFlags%predictionOnly .and.                                    &
         abs(region%states(1)%actuationAmount) > 0.0_wp)                                     &
         call this%controller%hookAfterTimemarch(region, FORWARD)

    call this%residualManager%cleanup()

    if (region%simulationFlags%computeTimeAverage) then
       do i = 1, size(region%states)
          region%states(i)%timeAverage = region%states(i)%timeAverage / (time - startTime)
       end do
       call region%saveData(QOI_TIME_AVERAGED_STATE, trim(this%outputPrefix) // ".mean.q")
    end if

    if (.not. region%simulationFlags%predictionOnly)                                         &
         functional = this%functional%runningTimeQuadrature

    call endTiming("Compute functional")

  end function runForward

  function runAdjoint(this, region) result(sensitivity)

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE, QOI_ADJOINT_STATE
    use SolverOptions_enum, only : FORWARD, ADJOINT

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    implicit none

    ! <<< Arguments >>>
    class(t_Solver) :: this
    class(t_Region) :: region

    ! <<< Result >>>
    real(SCALAR_KIND) :: sensitivity

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    character(len = STRING_LENGTH) :: filename
    integer :: i, j, timestep, startTimestep, timemarchDirection
    real(wp) :: time, startTime, timeStepSize, instantaneousSensitivity

    assert(.not. region%simulationFlags%predictionOnly)

    call startTiming("Compute sensitivity")

    sensitivity = 0.0_wp
    this%controller%runningTimeQuadrature = 0.0_wp

    ! Setup residual manager if this is a steady-state simulation
    if (region%simulationFlags%steadyStateSimulation)                                        &
         call this%residualManager%setup("adjoint_residuals", region)

    ! Load the initial condition.
    call initializeState(this, region, FORWARD) !... for control horizon end timestep.

    this%controller%onsetTime = region%states(1)%time
    this%controller%duration = this%nTimesteps * region%solverOptions%timeStepSize

    ! Load the adjoint coefficients corresponding to the end of the control time horizon.
    if (region%simulationFlags%steadyStateSimulation) then
       write(filename, '(2A)') trim(this%outputPrefix), ".steady_state.q"
    else
       write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-",                          &
            region%timestep + this%nTimesteps, ".q"
    end if
    call region%loadData(QOI_FORWARD_STATE, filename)
    do i = 1, size(region%states) !... update state
       call region%states(i)%update(region%grids(i), region%simulationFlags,                 &
            region%solverOptions)
    end do

    startTimestep = region%timestep
    startTime = region%states(1)%time

    ! Setup the revese-time migrator if this is not a steady-state simulation.
    if (.not. region%simulationFlags%steadyStateSimulation)                                  &
         call this%reverseMigrator%setup(region, this%timeIntegrator, this%outputPrefix,     &
         startTimestep - this%nTimesteps, startTimestep, this%saveInterval,                  &
         this%saveInterval * this%timeIntegrator%nStages)

    ! March forward for adjoint steady-state simulation.
    timemarchDirection = -1
    if (region%simulationFlags%steadyStateSimulation) timemarchDirection = +1

    ! Adjoint initial condition (if specified).
    call initializeState(this, region, ADJOINT)

    write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", region%timestep, ".adjoint.q"
    call region%saveData(QOI_ADJOINT_STATE, filename)

    ! Call controller hooks before time marching starts.
    call this%controller%hookBeforeTimemarch(region, ADJOINT)

    ! Reset probes.
    if (this%probeInterval > 0) call region%resetProbes()

    time = startTime

    do timestep = startTimestep + sign(1, timemarchDirection),                               &
         startTimestep + sign(this%nTimesteps, timemarchDirection), timemarchDirection

       region%timestep = timestep
       timeStepSize = region%computeTimeStepSize()

       do i = this%timeIntegrator%nStages, 1, -1

          ! Load adjoint coefficients.
          if (.not. region%simulationFlags%steadyStateSimulation) then
             if (i == 1) then
                call this%reverseMigrator%migrateTo(region, this%timeIntegrator,             &
                     timestep, this%timeIntegrator%nStages)
             else
                call this%reverseMigrator%migrateTo(region, this%timeIntegrator,             &
                     timestep + 1, i - 1)
             end if
          end if

          ! Update gradient.
          call this%controller%updateGradient(region)

          ! Update cost sensitivity.
          instantaneousSensitivity = this%controller%computeSensitivity(region)
          this%controller%runningTimeQuadrature = this%controller%runningTimeQuadrature +    &
               this%timeIntegrator%norm(i) * timeStepSize * instantaneousSensitivity

          ! Update adjoint forcing on cost target patches.
          call this%functional%updateAdjointForcing(region)

          ! Take a single adjoint sub-step using the time integrator.
          call this%timeIntegrator%substepAdjoint(region, time, timeStepSize, timestep, i)

          ! TODO: how to enforce limits on adjoint variables... check for NaN?
          if (region%simulationFlags%enableSolutionLimits)                                   &
               call checkSolutionLimits(region, ADJOINT, this%outputPrefix)

       end do

       ! Book-keeping, including simulation progress, residual monitoring, I/O, etc.
       call bookKeeping(this, region, ADJOINT, startTimestep, timestep,                      &
            time, instantaneousSensitivity)

       ! Save solution on probe patches.
       if (this%probeInterval > 0 .and. mod(timestep, max(1, this%probeInterval)) == 0)      &
            call region%saveProbeData(ADJOINT)

       ! Stop if this is a steady-state simulation and solution has converged.
       if (this%residualManager%hasSimulationConverged) exit

       ! Filter solution if required.
       if (region%simulationFlags%filterOn .and.                                             &
            mod(timestep, max(1, this%filterInterval)) == 0) then
          do j = 1, size(region%grids)
             call region%grids(j)%applyFilter(region%states(j)%adjointVariables, timestep)
          end do
       end if

    end do !... timestep = startTimestep + sign(1, timemarchDirection), ...

    ! Finish writing remaining data gathered on probes.
    if (this%probeInterval > 0) call region%saveProbeData(ADJOINT, finish = .true.)

    ! Call controller hooks after time marching ends.
    call this%controller%hookAfterTimemarch(region, ADJOINT)

    call this%residualManager%cleanup()
    call this%reverseMigrator%cleanup()

    sensitivity = this%controller%runningTimeQuadrature

    call endTiming("Compute sensitivity")

  end function runAdjoint

  subroutine checkGradientAccuracy(this, region)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

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
    real(wp) :: actuationAmount, functionalBaseline, functional, sensitivity,                &
         initialActuationAmount, geometricGrowthFactor, gradientError, dummyValue

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
          open(unit = getFreeUnit(fileUnit), file = trim(filename), action = 'write',        &
               status = 'unknown', iostat = iostat)
       else
          open(unit = getFreeUnit(fileUnit), file = trim(filename), action = 'readwrite',    &
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
       if (nIterations > 1) then
          call getRequiredOption("actuation_amount_geometric_growth", geometricGrowthFactor)
       else
          geometricGrowthFactor = getOption("actuation_amount_geometric_growth", 1.0_wp)
       end if
    end if

    ! Find (or load from file) the cost functional for the baseline prediction.
    if (region%simulationFlags%isBaselineAvailable) then
       if (procRank == 0)                                                                    &
            read(fileUnit, *, iostat = iostat) i, actuationAmount,                           &
            functionalBaseline, sensitivity, gradientError
       call MPI_Bcast(iostat, 1, MPI_INTEGER, 0, region%comm, ierror)
       if (iostat /= 0) then
          write(message, "(2A)") trim(filename),                                             &
               ": Failed to read baseline cost functional from file!"
          call gracefulExit(region%comm, message)
       end if
       call MPI_Bcast(functionalBaseline, 1, SCALAR_TYPE_MPI, 0, region%comm, ierror)
    else
       functionalBaseline = this%runForward(region)
    end if

    ! Find the sensitivity gradient (this is the only time the adjoint simulation will be run).
    if (restartIteration == 0) then
       sensitivity = this%runAdjoint(region)
    else
       call MPI_Bcast(sensitivity, 1, SCALAR_TYPE_MPI, 0, region%comm, ierror)
    end if

    if (procRank == 0 .and. .not. region%simulationFlags%isBaselineAvailable)                &
         write(fileUnit, '(I4,4(1X,SP,' // SCALAR_FORMAT // '))') 0, 0.0_wp,                 &
         functionalBaseline, sensitivity, 0.0_wp

    if (nIterations == 0) return

    ! Turn off output for controlled predictions.
    region%outputOn = .false.

    if (restartIteration == 0) restartIteration = restartIteration + 1

    do i = 1, restartIteration - 1
       if (procRank == 0)                                                                    &
            read(fileUnit, *, iostat = iostat) j, actuationAmount, functional,               &
            dummyValue, gradientError
       call MPI_Bcast(iostat, 1, MPI_INTEGER, 0, region%comm, ierror)
       if (iostat /= 0) then
          write(message, "(2A)") trim(filename),                                             &
               ": Cost functional history is too short for the specified restart iteration!"
          call gracefulExit(region%comm, message)
       end if
    end do

    do i = restartIteration, restartIteration + nIterations - 1
       actuationAmount = initialActuationAmount * geometricGrowthFactor ** real(i - 1, wp)
       functional = this%runForward(region, actuationAmount = actuationAmount)
       gradientError = (functional - functionalBaseline) / actuationAmount + sensitivity
       if (procRank == 0) then
          write(fileUnit, '(I4,4(1X,SP,' // SCALAR_FORMAT // '))') i, actuationAmount,       &
               functional, -(functional - functionalBaseline) / actuationAmount, gradientError
          flush(fileUnit)
       end if
    end do

    if (procRank == 0) close(fileUnit)

  end subroutine checkGradientAccuracy

  subroutine bookKeeping(this, region, mode, startTimestep,                                  &
       timestep, time, instantaneousQoI)

    ! <<< External modules >>>
    use iso_fortran_env, only : output_unit

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE, QOI_ADJOINT_STATE
    use SolverOptions_enum, only : FORWARD, ADJOINT

    ! <<< Internal modules >>>
    use ErrorHandler, only : writeAndFlush

    ! <<< Arguments >>>
    class(t_Solver) :: this
    class(t_Region) :: region
    integer, intent(in) :: mode, timestep, startTimestep
    real(SCALAR_KIND), intent(in) :: time, instantaneousQoI

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: nDimensions
    real(wp) :: timeStepSize, cfl, residuals(3)
    character(len = STRING_LENGTH) :: str, message, filename

    assert_key(mode, (FORWARD, ADJOINT))
    assert(startTimestep >= 0)

    nDimensions = size(region%globalGridSizes, 1)
    assert_key(nDimensions, (1, 2, 3))

    ! Report time step size in constant CFL mode
    if (region%simulationFlags%useConstantCfl) then
       timeStepSize = region%computeTimeStepSize()
    else
       cfl = region%computeCfl()
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
             write(message, '(2A,E13.6)') trim(str), ", cost = ", instantaneousQoI
          case (ADJOINT)
             write(message, '(2A,E13.6)') trim(str), ", gradient = ", instantaneousQoI
          end select

       else

         write(message, '(A)') trim(str)

       end if

       call writeAndFlush(region%comm, output_unit, message)

       if (.not. region%simulationFlags%predictionOnly .and. region%outputOn) then

          select case (mode)

          case (FORWARD)
             call this%functional%writeToFile(region%comm, trim(this%outputPrefix) //        &
                  ".cost_functional.txt", timestep, time,                                    &
                  timestep - startTimestep > this%reportInterval)

          case (ADJOINT)
             call this%controller%writeToFile(region%comm, trim(this%outputPrefix) //        &
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

  end subroutine bookKeeping

  subroutine checkSolutionLimits(region, mode, outputPrefix)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use State_mod, only : t_State
    use Region_mod, only : t_Region

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE
    use SolverOptions_enum, only : FORWARD

    ! <<< Internal modules >>>
    use ErrorHandler, only : gracefulExit

    ! <<< Arguments >>>
    class(t_Region) :: region
    integer, intent(in) :: mode
    character(len = *), intent(in) :: outputPrefix

    ! <<< Local variables >>>
    integer :: i, iGlobal, jGlobal, kGlobal, rankReportingError, procRank, ierror
    character(len = STRING_LENGTH) :: message
    real(SCALAR_KIND) :: fOutsideRange

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

  subroutine initializeState(this, region, mode, restartFilename)

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use Region_mod, only : t_Region
    use CostTargetPatch_mod, only : t_CostTargetPatch

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE, QOI_ADJOINT_STATE
    use SolverOptions_enum, only : FORWARD, ADJOINT

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

       filename = getOption("adjoint_initial_condition_file", "")
       if (len_trim(filename) == 0) then
          do i = 1, size(region%states)
             region%states(i)%adjointVariables = 0.0_wp
          end do
       else
          call region%loadData(QOI_ADJOINT_STATE, filename)
       end if

       if (allocated(region%patchFactories) .and.                                            &
            .not. region%simulationFlags%steadyStateSimulation .and.                         &
            .not. region%simulationFlags%useContinuousAdjoint) then

          timeStepSize = region%computeTimeStepSize()
          region%states(:)%adjointForcingFactor = - timeStepSize / 6.0_wp !... RK4 only.

          ! Connect to the previously allocated functional.
          call this%functional%updateAdjointForcing(region)

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

  end subroutine initializeState

end module Solver_mod
