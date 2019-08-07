#include "config.h"

program adjoint

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
  integer :: i, stat, fileUnit, dictIndex, procRank, numProcs, ierror
  character(len = STRING_LENGTH) :: filename, resultFilename, outputPrefix, message
  logical :: adjointRestart, success
  integer :: accumulatedNTimesteps
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
  if (.not. region%simulationFlags%isBaselineAvailable) then
    dummyValue = solver%runForward(region)
    if (procRank == 0) then
      resultFilename = trim(outputPrefix) // ".forward_run.txt"
      open(unit = getFreeUnit(fileUnit), file = trim(resultFilename), action='write',          &
        iostat = stat, status = 'replace')
      write(fileUnit, '(1X,SP,' // SCALAR_FORMAT // ')') dummyValue
      close(fileUnit)
    end if
  end if

  call get_command_argument(1, resultFilename, stat)
  if( (stat.ne.0) .or. (len_trim(resultFilename).le.0) )                                     &
    resultFilename = trim(outputPrefix) // ".adjoint_run.txt"

  dummyValue = 0.0_wp

  adjointRestart = getOption("adjoint_restart", .false.)
  accumulatedNTimesteps = -1
  if (adjointRestart) then ! This is relative timesteps: timestep at initial condition must be added.
    call getRequiredOption("adjoint_restart/accumulated_timesteps",accumulatedNTimesteps)
    assert(accumulatedNTimesteps.ge.0)
  end if
  if (adjointRestart .and. (accumulatedNTimesteps>0)) then
    if (procRank==0) then
      open(unit = getFreeUnit(fileUnit), file = trim(resultFilename), action='read',         &
        iostat = stat, status = 'old')
      read(fileUnit, '(1X,SP,' // SCALAR_FORMAT // ')') dummyValue
      close(fileUnit)
    end if
    call MPI_Bcast(dummyValue, 1, SCALAR_TYPE_MPI, 0, MPI_COMM_WORLD, ierror)
  end if
  dummyValue = dummyValue + solver%runAdjoint(region)
print *, procRank, ': adjoint run is finished.'
  if (procRank == 0) then
    open(unit = getFreeUnit(fileUnit), file = trim(resultFilename), action='write',          &
      iostat = stat, status = 'replace')
    write(fileUnit, '(1X,SP,' // SCALAR_FORMAT // ')') dummyValue
    close(fileUnit)
print *, procRank, ': adjoint sensitivity is saved.'
  end if

  call solver%cleanup()
  call region%cleanup()

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  ! Finalize MPI.
  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

end program adjoint
