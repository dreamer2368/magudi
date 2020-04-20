#include "config.h"

program control_mollifier_factor

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver

  use Grid_enum
  use State_enum

  use InputHelper, only : parseInputFile, getFreeUnit, getOption, getRequiredOption
  use ErrorHandler
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  use RegionImpl, only : normalizeControlMollifier, computeRegionIntegral

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, stat, fileUnit, procRank, numProcs, ierror
  integer :: kthArgument, numberOfArguments
  logical :: lookForInput = .false., lookForOutput = .false.,                                             &
              inputFlag = .false., outputFlag = .false.
  character(len = STRING_LENGTH) :: argument, inputFilename, outputFilename, restartFilename
  character(len = STRING_LENGTH) :: filename, outputPrefix, message
  logical :: fileExists, success
  integer, dimension(:,:), allocatable :: globalGridSizes
  SCALAR_TYPE, dimension(:,:), allocatable :: F
  logical, dimension(:), allocatable :: mask
  type(t_Region) :: region
  type(t_Solver) :: solver

  ! << output variables >>
  integer :: inputNumber, simulationNumber, intValue = 0
  SCALAR_TYPE :: dummyValue = 0.0_wp, maxValue = 0.0_wp

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

  if ( .not. outputFlag ) outputFilename = trim(outputPrefix) // ".control_mollifier_factor.txt"
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

  ! Read control mollifier file.
  call getRequiredOption("control_mollifier_file", filename)
  call region%loadData(QOI_CONTROL_MOLLIFIER, filename)

  ! Get maximum value of control mollifier before normalization
  maxValue = 0.0_wp
  do i = 1, size(region%grids)
    maxValue = MAX(maxValue, MAXVAL(region%grids(i)%controlMollifier))
  end do

  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE, maxValue, 1, REAL_TYPE_MPI,                        &
       MPI_MAX, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(maxValue, 1, REAL_TYPE_MPI, 0, region%grids(i)%comm, ierror)
  end do

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  write(message,'(A,1X,SP,' // SCALAR_FORMAT // ')') 'Before normalization: Maximum control mollifier value=', maxValue
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

  dummyValue = computeRegionIntegral(region)

  write(message,'(A,1X,SP,' // SCALAR_FORMAT // ')') 'Before normalization: Volume=', dummyValue
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

  ! Setup boundary conditions.
  call getRequiredOption("boundary_condition_file", filename)
  call region%setupBoundaryConditions(filename)

  ! Normalize control mollifier.
  call normalizeControlMollifier(region)

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  ! Get maximum value of control mollifier.
  maxValue = 0.0_wp
  do i = 1, size(region%grids)
    maxValue = MAX(maxValue, MAXVAL(region%grids(i)%controlMollifier))
  end do

  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE, maxValue, 1, REAL_TYPE_MPI,                        &
       MPI_MAX, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(maxValue, 1, REAL_TYPE_MPI, 0, region%grids(i)%comm, ierror)
  end do

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  write(message,'(A,1X,SP,' // SCALAR_FORMAT // ')') 'Maximum control mollifier value=', maxValue
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

  if ( outputFlag ) then
    if (procRank == 0) then
      open(unit = getFreeUnit(fileUnit), file = trim(outputFilename), action='write',          &
           iostat = stat, status = 'replace')
      write(fileUnit, '(1X,SP,' // SCALAR_FORMAT // ')') maxValue
      close(fileUnit)
    end if
  end if

  dummyValue = 0.0_wp
  do i = 1, size(region%grids)
    assert(region%grids(i)%nGridPoints > 0)
    allocate(F(region%grids(i)%nGridPoints, 2))
    F(:,2) = 1.0_wp
    F(:,1) = region%grids(i)%controlMollifier(:,1)

    dummyValue = dummyValue + region%grids(i)%computeInnerProduct(F(:,1),F(:,2))
    SAFE_DEALLOCATE(F)
  end do

  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE, dummyValue, 1, SCALAR_TYPE_MPI,                      &
       MPI_SUM, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(dummyValue, 1, SCALAR_TYPE_MPI, 0, region%grids(i)%comm, ierror)
  end do

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  write(message,'(A,1X,SP,' // SCALAR_FORMAT // ')') 'Control mollifier volume=', dummyValue/maxValue
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

  intValue = 0
  do i = 1, size(region%grids)
    allocate(mask(region%grids(i)%nGridPoints))
    mask = region%grids(i)%controlMollifier(:,1) > 0.0_wp
    intValue = intValue + count(mask)
    SAFE_DEALLOCATE(mask)
  end do

  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE, intValue, 1, MPI_INTEGER,                      &
       MPI_SUM, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(intValue, 1, MPI_INTEGER, 0, region%grids(i)%comm, ierror)
  end do

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  write(message,'(A,1X,SP,I10.5)') 'Control mollifier degree of freedom=', intValue
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

  call solver%cleanup()
  call region%cleanup()

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  ! Finalize MPI.
  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

end program control_mollifier_factor
