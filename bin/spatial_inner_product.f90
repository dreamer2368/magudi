#include "config.h"

program spatial_inner_product

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use Region_mod, only : t_Region
  use RegionImpl, only : normalizeControlMollifier, readBoundaryConditions

  use Grid_enum
  use State_enum

  use InputHelper, only : parseInputFile, getFreeUnit, getOption, getRequiredOption
  use ErrorHandler
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  type :: t_VectorInternal
     SCALAR_TYPE, allocatable :: F(:,:)
  end type t_VectorInternal

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, stat, fileUnit, procRank, numProcs, ierror
  integer :: kthArgument, numberOfArguments
  logical :: lookForInput = .false., lookForOutput = .false., lookForMollifier = .false.,                   &
              inputFlag = .false., outputFlag = .false., mollifierFlag = .false.
  character(len = STRING_LENGTH) :: argument, inputFilename, outputFilename, mollifierFilename
  character(len = STRING_LENGTH) :: filename, outputPrefix, message, Qfilenames(2)
  logical :: fileExists, success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_VectorInternal), allocatable :: temp(:)

  ! << output variables >>
  integer :: inputNumber, simulationNumber
  SCALAR_TYPE :: dummyValue = 0.0_wp

  ! Initialize MPI.
  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)

  call initializeErrorHandler()

  numberOfArguments = command_argument_count()
  if (numberOfArguments<2) then
    write(message, '(2A)') "spatial_inner_product requires the first 2 arguments ",                         &
                           "to be the filenames of two solution vectors."
    call gracefulExit(MPI_COMM_WORLD, message)
  end if
  do i=1,2
    call get_command_argument(i,argument)
    Qfilenames(i) = trim(adjustl(argument))
    if (procRank==0) inquire(file=Qfilenames(i),exist=fileExists)
    call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
    if (.not.fileExists) then
       write(message, '(3A)') "Q file ", trim(Qfilenames(i)), " does not exists!"
       call gracefulExit(MPI_COMM_WORLD, message)
    else
      write(message, '(2A)') "Q file: ", trim(Qfilenames(i))
      call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
    end if
  end do
  if (numberOfArguments>2) then
    do kthArgument = 3, numberOfArguments
      call get_command_argument(kthArgument,argument)
      select case(trim(adjustl(argument)))
      case("--input")
        lookForInput = .true.
      case("--output")
        lookForOutput = .true.
      case("--mollifier")
        lookForMollifier = .true.
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
        elseif (lookForMollifier) then
          mollifierFilename = trim(adjustl(argument))
          if (procRank==0) inquire(file=mollifierFilename,exist=fileExists)
          call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
          if (.not.fileExists) then
             write(message, '(3A)') "mollifier file ", trim(mollifierFilename), " does not exists!"
             call gracefulExit(MPI_COMM_WORLD, message)
          end if
          lookForMollifier = .false.
          mollifierFlag = .true.
        else
          write(message, '(3A)') "option ", trim(argument), " unknown!"
          call gracefulExit(MPI_COMM_WORLD, message)
        end if
      end select
    end do
  end if

  call startTiming("total")

  if ( .not. inputFlag ) inputFilename = PROJECT_NAME // ".inp"
  ! Parse options from the input file.
  call parseInputFile(inputFilename)

  outputPrefix = getOption("output_prefix", PROJECT_NAME)

  if ( .not. outputFlag ) outputFilename = trim(outputPrefix) // ".spatial_inner_product.txt"
  write(message, '(2A)') "Input file: ", trim(inputFilename)
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  write(message, '(2A)') "Output file: ", trim(outputFilename)
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  if (mollifierFlag) then
    write(message, '(2A)') "Mollifier file: ", trim(mollifierFilename)
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  end if

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

  call startTiming("load solution variables")

  ! Load Q solution vectors.
  do i=1,2
    call region%loadData(QOI_FORWARD_STATE, Qfilenames(i))
    if(i==2) exit
    allocate(temp(size(region%states)))
    do j=1,size(region%states)
      allocate(temp(j)%F(region%grids(j)%nGridPoints,region%solverOptions%nUnknowns))
      temp(j)%F = region%states(j)%conservedVariables
      ! temp(j)%F(:,1:region%solverOptions%nUnknowns-1) = 0.0_wp
    end do
  end do

  call endTiming("load solution variables")

  call startTiming("load mollifier")

  if (mollifierFlag) then
    ! Setup boundary conditions.
    call getRequiredOption("boundary_condition_file", filename)
    call region%setupBoundaryConditions(filename)

    call region%loadData(QOI_CONTROL_MOLLIFIER, mollifierFilename)
    ! call normalizeControlMollifier(region)
  else
    do i = 1, size(region%grids)
       region%grids(i)%controlMollifier = 1.0_wp
    end do
  end if

  call endTiming("load mollifier")

  call startTiming("compute inner products")

  dummyValue = 0.0_wp

  do i=1,size(region%states)
    dummyValue = dummyValue +                                                                &
                  region%grids(i)%computeInnerProduct(region%states(i)%conservedVariables,   &
                                        temp(i)%F, region%grids(i)%controlMollifier(:,1) )
  end do

  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE, dummyValue, 1,                                       &
       SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(dummyValue, 1, SCALAR_TYPE_MPI,                                          &
          0, region%grids(i)%comm, ierror)
  end do

  call endTiming("compute inner products")

  write(message, '(A,1X,SP,' // SCALAR_FORMAT // ')') 'Inner product = ', dummyValue
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

  if (procRank == 0) then
    fileUnit = getFreeUnit()
    open(unit = fileUnit, file = trim(outputFilename), action='write',          &
         iostat = stat, status = 'replace')
    write(fileUnit, '(1X,SP,' // SCALAR_FORMAT // ')') dummyValue
    close(fileUnit)
  end if

  do i=1,size(region%states)
    SAFE_DEALLOCATE(temp(i)%F)
  end do
  SAFE_DEALLOCATE(temp)
  call region%cleanup()

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  ! Finalize MPI.
  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

end program spatial_inner_product
