#include "config.h"

program traveling_wave

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use Region_mod, only : t_Region
  use Optimizer_mod, only : t_Optimizer

  use Grid_enum
  use State_enum
  use Optimier_enum

  use InputHelper, only : parseInputFile, getFreeUnit, getOption, getRequiredOption
  use ErrorHandler
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers
  use TravelingWaveImpl, only : setupTravelingWave

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, stat, fileUnit, procRank, numProcs, ierror
  integer :: kthArgument, numberOfArguments, mode = -1
  logical :: lookForInput = .false., lookForOutput = .false., lookForRestart = .false.,      &
            inputFlag = .false., outputFlag = .false., restartFlag = .false., saveMetricsFlag = .false.
  character(len = STRING_LENGTH) :: argument, inputFilename, outputFilename, restartFilename
  character(len = STRING_LENGTH) :: filename, outputPrefix, message
  logical :: fileExists, success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_Optimizer) :: optimizer

  ! << output variables >>
  integer :: inputNumber, simulationNumber
  SCALAR_TYPE :: dummyValue = 0.0_wp

  ! Initialize MPI.
  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)

  call initializeErrorHandler()

  numberOfArguments = command_argument_count()
  do kthArgument = 1, numberOfArguments
    call get_command_argument(kthArgument,argument)
    select case(trim(adjustl(argument)))
    case("--input")
      lookForInput = .true.
    case("--output")
      lookForOutput = .true.
    case("--restart")
      lookForRestart = .true.
    case("--verify")
      mode = VERIFY
    case("--nlcg")
      mode = NLCG
    case("--save_metrics")
      saveMetricsFlag = .true.
    case default
      if (lookForInput) then
        inputFilename = trim(adjustl(argument))
        if (procRank==0) inquire(file=inputFilename,exist=fileExists)
        call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
        if (.not.fileExists) then
           write(message, '(3A)') "input file ", trim(inputFilename), " does not exists!"
           call gracefulExit(MPI_COMM_WORLD, message)
        end if
        lookForInput = .false.
        inputFlag = .true.
      elseif (lookForOutput) then
        outputFilename = trim(adjustl(argument))
        lookForOutput = .false.
        outputFlag = .true.
      elseif (lookForRestart) then
        restartFilename = trim(adjustl(argument))
        if (procRank==0) inquire(file=restartFilename,exist=fileExists)
        call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
        if (.not.fileExists) then
           write(message, '(3A)') "restart file ", trim(restartFilename), " does not exists!"
           call gracefulExit(MPI_COMM_WORLD, message)
        end if
        lookForRestart = .false.
        restartFlag = .true.
      else
        write(message, '(3A)') "option ", trim(argument), " unknown!"
        call gracefulExit(MPI_COMM_WORLD, message)
      end if
    end select
  end do

  call startTiming("total")

  if ( .not. inputFlag ) inputFilename = PROJECT_NAME // ".inp"
  ! Parse options from the input file.
  call parseInputFile(inputFilename)

  outputPrefix = getOption("output_prefix", PROJECT_NAME)

  if ( .not. outputFlag ) outputFilename = trim(outputPrefix) // ".forward_run.txt"
  write(message, '(2A)') "Input file: ", trim(inputFilename)
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  write(message, '(2A)') "Output file: ", trim(outputFilename)
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

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
  if (saveMetricsFlag) then
    write(filename, '(2A)') trim(outputPrefix), ".Jacobian.f"
    call region%saveData(QOI_JACOBIAN, filename)
    write(filename, '(2A)') trim(outputPrefix), ".metrics.f"
    call region%saveData(QOI_METRICS, filename)
  end if

  call setupTravelingWave(region)

  ! Initialize the optimizer.
  call optimizer%setup(region)

  ! Save the control and target mollifier if using code-generated values.
  assert(region%simulationFlags%enableFunctional)
  filename = getOption("target_mollifier_file", "")
  if (len_trim(filename) == 0) call region%saveData(QOI_TARGET_MOLLIFIER,                 &
       trim(outputPrefix) // ".target_mollifier.f")

  ! Main code logic.
  select case (mode)
  case(VERIFY)
    write(message, '(A)') "Running on the verification mode."
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
    call optimizer%verifyAdjoint(region)
  case(NLCG)
    if (restartFlag) then
      call optimizer%runNLCG(region, trim(restartFilename))
    else
      call optimizer%runNLCG(region)
    end if
  case default
    write(message, '(A)') 'Mode is not specified correctly!: --verify, --nlcg, --newton'
    call gracefulExit(MPI_COMM_WORLD, message)
  end select

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  ! if (procRank == 0) then
  !   open(unit = getFreeUnit(fileUnit), file = trim(outputFilename), action='write',          &
  !        iostat = stat, status = 'replace')
  !   write(fileUnit, '(1X,SP,' // SCALAR_FORMAT // ')') dummyValue
  !   close(fileUnit)
  ! end if

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  call region%cleanup()
  call optimizer%cleanup()

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  ! Finalize MPI.
  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

end program traveling_wave
