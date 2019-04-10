#include "config.h"

program control_space_norm

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver

  use Grid_enum
  use State_enum

  use InputHelper, only : parseInputFile, getFreeUnit, getOption, getRequiredOption
  use InputHelperImpl, only: dict, find
  use ErrorHandler
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, stat, fileUnit, dictIndex, procRank, numProcs, ierror, STATUS
  character(len = STRING_LENGTH) :: filename, resultFilename, outputPrefix, message
  logical :: success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_Solver) :: solver

  ! << output variables >>
  integer :: inputNumber, simulationNumber
  SCALAR_TYPE :: dummyValue = 0.0_wp

  ! Initialize MPI.
  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)

  call initializeErrorHandler()

  if (command_argument_count() > 0) then
     write(message, '(A)') "Usage: magudi [FILE]"
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)') "High-performance Fortran-based adjoint optimization tool."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)')                                                                   &
          "FILE is an optional restart file used if running in prediction-only mode."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     call cleanupErrorHandler()
     call MPI_Finalize(ierror)
     stop -1
  end if

  call startTiming("total")

  ! Parse options from the input file.
  filename = PROJECT_NAME // ".inp"
  call parseInputFile(filename)

  ! adjoint run options
  call find("enable_adjoint_solver",dictIndex)
  dict(dictIndex)%val = "true"

  outputPrefix = getOption("output_prefix", PROJECT_NAME)

  ! Verify that the grid file is in valid PLOT3D format and fetch the grid dimensions:
  ! `globalGridSizes(i,j)` is the number of grid points on grid `j` along dimension `i`.
  call getRequiredOption("grid_file", filename)
  call plot3dDetectFormat(MPI_COMM_WORLD, filename, success,                                 &
       globalGridSizes = globalGridSizes)
  if (.not. success) call gracefulExit(MPI_COMM_WORLD, plot3dErrorMessage)

  ! Setup the region.
  call region%setup(MPI_COMM_WORLD, globalGridSizes)

  ! Read the grid file.
  call getRequiredOption("grid_file", filename)
  call region%loadData(QOI_GRID, filename)

  ! Update the grids by computing the Jacobian, metrics, and norm.
  do i = 1, size(region%grids)
     call region%grids(i)%update()
  end do
  call MPI_Barrier(region%comm, ierror)

  ! Write out some useful information.
  call region%reportGridDiagnostics()

  ! Save the Jacobian and normalized metrics.
  write(filename, '(2A)') trim(outputPrefix), ".Jacobian.f"
  call region%saveData(QOI_JACOBIAN, filename)
  write(filename, '(2A)') trim(outputPrefix), ".metrics.f"
  call region%saveData(QOI_METRICS, filename)

  ! Initialize the solver.
  call solver%setup(region, outputPrefix = outputPrefix)

  ! Save the control and target mollifier if using code-generated values.
  if (region%simulationFlags%enableController) then
     filename = getOption("control_mollifier_file", "")
     if (len_trim(filename) == 0) call region%saveData(QOI_CONTROL_MOLLIFIER,                &
          trim(outputPrefix) // ".control_mollifier.f")
  end if
  if (region%simulationFlags%enableFunctional) then
     filename = getOption("target_mollifier_file", "")
     if (len_trim(filename) == 0) call region%saveData(QOI_TARGET_MOLLIFIER,                 &
          trim(outputPrefix) // ".target_mollifier.f")
  end if

  ! Main code logic.
  call collectControlSpaceNorm(solver, region)
  ! dummyValue = solver%runAdjoint(region)
  ! call get_command_argument(1, resultFilename, STATUS)
  ! if( STATUS .eq. 0 )                                                                      &
  !   resultFilename = trim(outputPrefix) // ".adjoint_run.txt"
  ! if (procRank == 0) then
  !   open(unit = getFreeUnit(fileUnit), file = trim(resultFilename), action='write',          &
  !     iostat = stat, status = 'replace')
  !   write(fileUnit, '(1X,SP,' // SCALAR_FORMAT // ')') dummyValue
  !   close(fileUnit)
  ! end if

  call solver%cleanup()
  call region%cleanup()

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  ! Finalize MPI.
  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

contains

  subroutine collectControlSpaceNorm(this, region)

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

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    character(len = STRING_LENGTH) :: filename, message
    class(t_Patch), pointer :: patch => null()
    class(t_TimeIntegrator), pointer :: timeIntegrator => null()
    class(t_Controller), pointer :: controller => null()
    class(t_Functional), pointer :: functional => null()
    type(t_ReverseMigratorFactory) :: reverseMigratorFactory
    class(t_ReverseMigrator), pointer :: reverseMigrator => null()
    integer :: i, j, timestep, startTimestep, timemarchDirection
    real(SCALAR_KIND) :: time, startTime, timeStepSize

    assert(region%simulationFlags%enableController)
    assert(region%simulationFlags%enableFunctional)
    assert(region%simulationFlags%enableAdjoint)

    call startTiming("collectControlSpaceNorm")

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
    controller%onsetTime = region%states(1)%time
    controller%duration = this%nTimesteps * region%solverOptions%timeStepSize

    ! Load the adjoint coefficients corresponding to the end of the control time horizon.
    if (region%simulationFlags%steadyStateSimulation) then
       write(filename, '(2A)') trim(this%outputPrefix), ".steady_state.q"
    else
       write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-",                            &
            region%timestep + this%nTimesteps, ".q"
    end if
    call region%loadData(QOI_FORWARD_STATE, filename)
    ! do i = 1, size(region%states) !... update state
    !    call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
    !         region%solverOptions)
    ! end do

    startTimestep = region%timestep
    startTime = region%states(1)%time

    ! ! Connect to the previously allocated reverse migrator.
    ! call reverseMigratorFactory%connect(reverseMigrator,                                       &
    !      region%solverOptions%checkpointingScheme)
    ! assert(associated(reverseMigrator))

    ! ! Setup the revese-time migrator if this is not a steady-state simulation.
    ! if (.not. region%simulationFlags%steadyStateSimulation)                                    &
    !      call reverseMigrator%setup(region, timeIntegrator, this%outputPrefix,                 &
    !      startTimestep - this%nTimesteps, startTimestep, this%saveInterval,                    &
    !      this%saveInterval * timeIntegrator%nStages)

    ! March forward for adjoint steady-state simulation.
    timemarchDirection = -1
    if (region%simulationFlags%steadyStateSimulation) timemarchDirection = 1

    ! ! Adjoint initial condition (if specified).
    ! call loadInitialCondition(this, region, ADJOINT)

    ! write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", region%timestep, ".adjoint.q"
    ! call region%saveData(QOI_ADJOINT_STATE, filename)

    ! Put norm filename into gradient filename and use gradient filename.
    do i = 1, size(region%patchFactories)
      call region%patchFactories(i)%connect(patch)
      if (.not. associated(patch)) cycle
      do j = 1, size(region%states)
        if (patch%gridIndex /= region%grids(j)%index .or. patch%nPatchPoints <= 0) cycle
        select type (patch)
        class is (t_ActuatorPatch)
          write(filename, '(4A)') trim(this%outputPrefix), ".norm_", trim(patch%name), ".dat"
          patch%gradientFilename = trim(filename)
        end select
      end do
    end do
    ! Call controller hooks before time marching starts.
    call controller%hookBeforeTimemarch(region, ADJOINT)
    ! !!!SeungWhan: need additional execution with FORWARD, in case of non-zero control forcing.
    ! if (controller%controllerSwitch) then
    !   controller%duration = this%nTimesteps * region%solverOptions%timeStepSize
    !   controller%onsetTime = startTime - controller%duration
    !   call controller%hookBeforeTimemarch(region, FORWARD)
    ! end if

    ! Reset probes.
    if (this%probeInterval > 0) call region%resetProbes()

    time = startTime

    do timestep = startTimestep + sign(1, timemarchDirection),                                 &
         startTimestep + sign(this%nTimesteps, timemarchDirection), timemarchDirection

       region%timestep = timestep
       timeStepSize = region%getTimeStepSize() !This currently only works for constant time step size.

       do i = timeIntegrator%nStages, 1, -1

          ! ! Load adjoint coefficients.
          ! if (.not. region%simulationFlags%steadyStateSimulation) then !... unsteady simulation.
          !    if (i == 1) then
          !       call reverseMigrator%migrateTo(region, controller, timeIntegrator,             &
          !            timestep, timeIntegrator%nStages)
          !    else
          !       call reverseMigrator%migrateTo(region, controller, timeIntegrator,             &
          !            timestep + 1, i - 1)
          !    end if
          ! end if

          ! Collect norm
          call controller%collectNorm(region, timeIntegrator%norm(i)*timeStepSize)

          ! ! Update cost sensitivity.
          ! instantaneousCostSensitivity = controller%computeSensitivity(region)
          ! controller%runningTimeQuadrature = controller%runningTimeQuadrature +                &
          !      timeIntegrator%norm(i) * timeStepSize * instantaneousCostSensitivity

          ! ! Update adjoint forcing on cost target patches.
          ! ! SeungWhan: Bug fix for final step
          ! if( (timestep.eq.startTimestep+sign(this%nTimesteps,timemarchDirection)) .and.       &
          !     (i.eq.1) ) then
          !     IS_FINAL_STEP = .true.
          ! else
          !     IS_FINAL_STEP = .false.
          ! end if
          ! call functional%updateAdjointForcing(region,IS_FINAL_STEP)
          !
          ! ! Take a single adjoint sub-step using the time integrator.
          ! call timeIntegrator%substepAdjoint(region, time, timeStepSize, timestep, i)
          !
          ! ! TODO: how to enforce limits on adjoint variables... check for NaN?
          ! if (region%simulationFlags%enableSolutionLimits)                                     &
          !      call checkSolutionLimits(region, ADJOINT, this%outputPrefix)

       end do

       ! Report simulation progress.
       write(message, '(A,I8)') 'Collected control space norm at the time step = ', timestep
       call writeAndFlush(region%comm, output_unit, message)

       ! ! Save solution on probe patches.
       ! if (this%probeInterval > 0 .and. mod(timestep, max(1, this%probeInterval)) == 0)        &
       !      call region%saveProbeData(ADJOINT)

       ! ! Stop if this is a steady-state simulation and solution has converged.
       ! if (this%residualManager%hasSimulationConverged) exit

       ! ! Filter solution if required.
       ! if (region%simulationFlags%filterOn) then
       !    do j = 1, size(region%grids)
       !       call region%grids(j)%applyFilter(region%states(j)%adjointVariables, timestep)
       !    end do
       ! end if

    end do !... timestep = startTimestep + sign(1, timemarchDirection), ...

    ! ! Finish writing remaining data gathered on probes.
    ! if (this%probeInterval > 0) call region%saveProbeData(ADJOINT, finish = .true.)

    ! Call controller hooks after time marching ends.
    if (controller%controllerSwitch) call controller%hookAfterTimemarch(region, FORWARD)
    call controller%hookAfterTimemarch(region, ADJOINT)

    call this%residualManager%cleanup()
    ! call reverseMigratorFactory%cleanup()

    ! costSensitivity = controller%runningTimeQuadrature

    write(message, '(A)') 'Control space norm collection is finished.'
    call writeAndFlush(region%comm, output_unit, message)

    call endTiming("collectControlSpaceNorm")

  end subroutine collectControlSpaceNorm

end program control_space_norm
