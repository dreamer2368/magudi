#include "config.h"

program main

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver
  use Functional_mod, only : t_Functional
  use Controller_mod, only : t_Controller
  use TimeIntegrator_mod, only : t_TimeIntegrator
  use Functional_factory, only : t_FunctionalFactory
  use Controller_factory, only : t_ControllerFactory
  use ReverseMigrator_mod, only : t_ReverseMigrator
  use TimeIntegrator_factory, only : t_TimeIntegratorFactory
  use ReverseMigrator_factory, only : t_ReverseMigratorFactory

  use Grid_enum
  use State_enum

  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, procRank, numProcs, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix, message
  logical :: success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_TimeIntegratorFactory) :: timeIntegratorFactory
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()
  type(t_FunctionalFactory) :: functionalFactory
  class(t_Functional), pointer :: functional => null()
  type(t_ControllerFactory) :: controllerFactory
  class(t_Controller), pointer :: controller => null()
  type(t_ReverseMigratorFactory) :: reverseMigratorFactory
  class(t_ReverseMigrator), pointer :: reverseMigrator => null()
  type(t_Solver) :: solver
  real(SCALAR_KIND) :: dummyValue

  ! Initialize MPI.
  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)

  call initializeErrorHandler()

  if (command_argument_count() > 1) then
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

  ! Read the target state if specified, otherwise use a quiescent state.
  if (region%simulationFlags%useTargetState) then
     filename = getOption("target_state_file", "")
     if (len_trim(filename) == 0) then
        do i = 1, size(region%states)
           call region%states(i)%makeQuiescent(region%grids(i)%nDimensions,                  &
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

  ! Setup time integrator.
  call timeIntegratorFactory%connect(timeIntegrator,                                         &
       trim(region%solverOptions%timeIntegratorType))
  if (.not. associated(timeIntegrator)) then
     write(message, '(3A)') "Invalid time integration scheme '",                             &
          trim(region%solverOptions%timeIntegratorType), "'!"
     call gracefulExit(MPI_COMM_WORLD, message)
  end if
  call timeIntegrator%setup(region)

  if (.not. region%simulationFlags%predictionOnly) then

     ! Connect to and setup functional.
     call functionalFactory%connect(functional, trim(region%solverOptions%costFunctionalType))
     if (.not. associated(functional)) then
        write(message, '(3A)') "Invalid cost functional '",                                  &
             trim(region%solverOptions%costFunctionalType), "'!"
        call gracefulExit(MPI_COMM_WORLD, message)
     end if
     call functional%setup(region)

     ! Connect to and setup controller.
     call controllerFactory%connect(controller, trim(region%solverOptions%controllerType))
     if (.not. associated(controller)) then
        write(message, '(3A)') "Invalid controller '",                                       &
             trim(region%solverOptions%controllerType), "'!"
        call gracefulExit(MPI_COMM_WORLD, message)
     end if
     call controller%setup(region)

     ! Connect to reverse migrator (setup will be handled by solver).
     call reverseMigratorFactory%connect(reverseMigrator,                                    &
          trim(region%solverOptions%checkpointingScheme))
     if (.not. associated(reverseMigrator)) then
        write(message, '(3A)') "Invalid checkpointing scheme '",                             &
             trim(region%solverOptions%checkpointingScheme), "'!"
        call gracefulExit(MPI_COMM_WORLD, message)
     end if

     ! Initialize the solver with specified functional, controller and reverse migrator.
     call solver%setup(timeIntegrator, functional, controller, reverseMigrator, outputPrefix)

  else

     ! Initialize the solver.
     call solver%setup(timeIntegrator, outputPrefix = outputPrefix)

  end if

  ! Save the control and target mollifier if using code-generated values.
  if (.not. region%simulationFlags%predictionOnly) then
     filename = getOption("control_mollifier_file", "")
     if (len_trim(filename) == 0) call region%saveData(QOI_CONTROL_MOLLIFIER,                &
          trim(outputPrefix) // ".control_mollifier.f")
     filename = getOption("target_mollifier_file", "")
     if (len_trim(filename) == 0) call region%saveData(QOI_TARGET_MOLLIFIER,                 &
          trim(outputPrefix) // ".target_mollifier.f")
  end if

  ! Main code logic.
  if (region%simulationFlags%predictionOnly) then !... just a predictive simulation.
     if (command_argument_count() == 1) then
        call get_command_argument(1, filename)
        dummyValue = solver%runForward(region, restartFilename = filename)
     else
        dummyValue = solver%runForward(region)
     end if
  else if (getOption("single_controlled_prediction", .false.)) then
     dummyValue = solver%runForward(region, actuationAmount =                                &
          getOption("actuation_amount", 1.0_wp))
  else
     if (.not. region%simulationFlags%isBaselineAvailable)                                   &
          dummyValue = solver%runForward(region)
     if (.not. getOption("gradient_available", .false.))                                     &
          dummyValue = solver%runAdjoint(region)
  end if

  call solver%cleanup()

  if (associated(reverseMigrator)) then
     call reverseMigrator%cleanup()
     nullify(reverseMigrator)
  end if

  if (associated(controller)) then
     call controller%cleanup()
     nullify(controller)
  end if

  if (associated(functional)) then
     call functional%cleanup()
     nullify(functional)
  end if

  if (associated(timeIntegrator)) then
     call timeIntegrator%cleanup()
     nullify(timeIntegrator)
  end if

  call region%cleanup()

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  ! Finalize MPI.
  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

end program main
