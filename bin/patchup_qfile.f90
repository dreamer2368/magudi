#include "config.h"

program patchup_qfile

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use Region_mod, only : t_Region

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
  integer :: i, j, procRank, numProcs, ierror
  integer :: kthArgument, numberOfArguments
  logical :: lookForInput = .false., inputFlag = .false.
  character(len = STRING_LENGTH) :: argument, inputFilename, mollifierFilename
  character(len = STRING_LENGTH) :: filename, outputPrefix, message,                                          &
                                    Zfilename, Xfilename, Yfilename
  logical :: fileExists, zeroYfile = .false., success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_VectorInternal), allocatable :: X(:), Y(:), antiW(:)

  ! Initialize MPI.
  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)

  call initializeErrorHandler()

  numberOfArguments = command_argument_count()
  if (numberOfArguments<4) then
    write(message, '(2A)') "patch-up_qfile requires the first 4 arguments: ",                         &
                           "z = w * x + (1-w) * y."
    call gracefulExit(MPI_COMM_WORLD, message)
  end if
  call get_command_argument(1,argument)
  Zfilename = trim(adjustl(argument))
  write(message, '(2A)') "Z file: ", trim(Zfilename)
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

  call get_command_argument(2, argument)
  mollifierFilename = trim(adjustl(argument))
  if (procRank==0) inquire(file=mollifierFilename,exist=fileExists)
  call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
  if (.not.fileExists) then
     write(message, '(3A)') "Mollifier file ", trim(mollifierFilename), " does not exists!"
     call gracefulExit(MPI_COMM_WORLD, message)
  else
    write(message, '(2A)') "Mollifier file: ", trim(mollifierFilename)
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  end if

  call get_command_argument(3,argument)
  Xfilename = trim(adjustl(argument))
  if (procRank==0) inquire(file=Xfilename,exist=fileExists)
  call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
  if (.not.fileExists) then
     write(message, '(3A)') "X file ", trim(Xfilename), " does not exists!"
     call gracefulExit(MPI_COMM_WORLD, message)
  else
    write(message, '(2A)') "X file: ", trim(Xfilename)
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  end if

  call get_command_argument(4,argument)
  if (trim(adjustl(argument))=="--zero") then
    zeroYfile = .true.
  else
    Yfilename = trim(adjustl(argument))
    if (procRank==0) inquire(file=Yfilename,exist=fileExists)
    call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
    if (.not.fileExists) then
       write(message, '(3A)') "Y file ", trim(Yfilename), " does not exists!"
       call gracefulExit(MPI_COMM_WORLD, message)
    else
      write(message, '(2A)') "Y file: ", trim(Yfilename)
      call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
    end if
    zeroYfile = .false.
  end if

  if (numberOfArguments>4) then
    do kthArgument = 5, numberOfArguments
      call get_command_argument(kthArgument,argument)
      select case(trim(adjustl(argument)))
      case("--input")
        lookForInput = .true.
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

  write(message, '(2A)') "Input file: ", trim(inputFilename)
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

  call startTiming("load solutions")

  ! Load Q solution vectors.
  call region%loadData(QOI_FORWARD_STATE, Xfilename)
  allocate(X(size(region%states)))
  do j=1,size(region%states)
    allocate(X(j)%F(region%grids(j)%nGridPoints,region%solverOptions%nUnknowns))
    X(j)%F = region%states(j)%conservedVariables
  end do

  if (.not. zeroYfile) then
    call region%loadData(QOI_FORWARD_STATE, Yfilename)
    allocate(Y(size(region%states)))
    do j=1,size(region%states)
      allocate(Y(j)%F(region%grids(j)%nGridPoints,region%solverOptions%nUnknowns))
      Y(j)%F = region%states(j)%conservedVariables
    end do
  end if

  call endTiming("load solutions")

  call startTiming("load mollifier")

  call region%loadData(QOI_CONTROL_MOLLIFIER, mollifierFilename)
  allocate(antiW(size(region%states)))
  do j=1,size(region%states)
    allocate(antiW(j)%F(region%grids(j)%nGridPoints,1))
    antiW(j)%F = 1.0_wp - region%grids(j)%controlMollifier
  end do

  call endTiming("load mollifier")

  call startTiming("compute Z = W * X + (1-W) * Y")

  do i = 1, size(region%grids)
    do j = 1, region%solverOptions%nUnknowns
      region%states(i)%conservedVariables(:,j) = region%grids(i)%controlMollifier(:,1) *                           &
                                                    X(i)%F(:,j)
      if (.not. zeroYfile) then
        region%states(i)%conservedVariables(:,j) = region%states(i)%conservedVariables(:,j) +                      &
                                                      antiW(i)%F(:,1) * Y(i)%F(:,j)
      end if
    end do
  end do
  call region%saveData(QOI_FORWARD_STATE, Zfilename)

  call endTiming("compute Z = W * X + (1-W) * Y")

  do i=1,size(region%states)
    SAFE_DEALLOCATE(X(i)%F)
    SAFE_DEALLOCATE(antiW(i)%F)
  end do
  SAFE_DEALLOCATE(X)
  SAFE_DEALLOCATE(antiW)
  if (allocated(Y)) then
    do i = 1, size(Y)
      SAFE_DEALLOCATE(Y(i)%F)
    end do
  end if
  SAFE_DEALLOCATE(Y)
  call region%cleanup()

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  ! Finalize MPI.
  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

end program patchup_qfile
