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
  integer :: i, dictIndex, procRank, numProcs, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix, message
  logical :: success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_Solver) :: solver

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

    startTimestep = region%timestep + this%nTimesteps
    startTime = region%states(1)%time + this%nTimesteps * region%solverOptions%timeStepSize

    ! March forward for adjoint steady-state simulation.
    timemarchDirection = -1
    if (region%simulationFlags%steadyStateSimulation) timemarchDirection = 1

    ! Put norm filename into gradient filename and use gradient filename.
    if (allocated(region%patchFactories)) then
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
    end if
    ! Call controller hooks before time marching starts.
    call controller%hookBeforeTimemarch(region, ADJOINT)

    ! Reset probes.
    if (this%probeInterval > 0) call region%resetProbes()

    time = startTime

    do timestep = startTimestep + sign(1, timemarchDirection),                                 &
         startTimestep + sign(this%nTimesteps, timemarchDirection), timemarchDirection

       region%timestep = timestep
       timeStepSize = region%getTimeStepSize() !This currently only works for constant time step size.

       do i = timeIntegrator%nStages, 1, -1

          ! Collect norm
          call controller%collectNorm(region, timeIntegrator%norm(i) * timeStepSize)

       end do

       ! Report simulation progress.
       write(message, '(A,I8)') 'Collected control space norm at the time step = ', timestep
       call writeAndFlush(region%comm, output_unit, message)

    end do !... timestep = startTimestep + sign(1, timemarchDirection), ...

    ! Call controller hooks after time marching ends.
    if (controller%controllerSwitch) call controller%hookAfterTimemarch(region, FORWARD)
    call controller%hookAfterTimemarch(region, ADJOINT)

    call this%residualManager%cleanup()

    write(message, '(A)') 'Control space norm collection is finished.'
    call writeAndFlush(region%comm, output_unit, message)

    call endTiming("collectControlSpaceNorm")

  end subroutine collectControlSpaceNorm

end program control_space_norm
