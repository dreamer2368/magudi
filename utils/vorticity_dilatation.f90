#include "config.h"

program vorticity_dilatation

  use MPI
  use, intrinsic :: iso_fortran_env

  use Grid_enum
  use State_enum

  use Region_mod, only : t_Region

  use CNSHelper, only : computeDependentVariables
  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush, gracefulExit
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage

  implicit none

  integer :: i, j, startTimestep, endTimestep, saveInterval, procRank, ierror
  character(len = STRING_LENGTH) :: filename, message, programName, outputPrefix
  logical :: success
  type(t_Region) :: region
  integer, allocatable :: globalGridSizes(:,:)

  ! Initialize MPI.
  call MPI_Init(ierror)

  ! Parse options from the input file.
  filename = "vorticity_dilatation.inp"
  call parseInputFile(filename)

  ! Verify that the grid file is in valid PLOT3D format and fetch the grid dimensions:
  ! `globalGridSizes(i,j)` is the number of grid points on grid `j` along dimension `i`.
  call getRequiredOption("grid_file", filename)
  call plot3dDetectFormat(MPI_COMM_WORLD, filename,                                          &
       success, globalGridSizes = globalGridSizes)
  if (.not. success) call gracefulExit(MPI_COMM_WORLD, plot3dErrorMessage)

  ! Setup the region and load the grid file.
  call region%setup(MPI_COMM_WORLD, globalGridSizes)
  call region%loadData(QOI_GRID, filename)

  ! Compute normalized metrics, norm matrix and Jacobian.
  do i = 1, size(region%grids)
     call region%grids(i)%update()
  end do
  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  if (command_argument_count() >= 1) then !... only one solution file to process.

     ! Load the solution file.
     call get_command_argument(1, filename)
     call region%loadData(QOI_FORWARD_STATE, filename)

     ! Compute velocity
     do i = 1, size(region%grids)
        call computeDependentVariables(size(globalGridSizes, 1),                                &
             region%states(i)%conservedVariables, velocity = region%states(i)%velocity)
     end do
     
     ! Save vorticity and dilatation
     i = len_trim(filename)
     if (filename(i-1:i) == ".q") then
        filename = filename(:i-2) // ".vorticity_dilatation.f"
     else
        filename = PROJECT_NAME // ".vorticity_dilatation.f"
     end if
     call region%saveData(QOI_VORTICITY_DILATATION, filename)

  else

     call getRequiredOption("save_interval", saveInterval)
     call getRequiredOption("start_timestep", startTimestep)
     call getRequiredOption("end_timestep", endTimestep)

     outputPrefix = getOption("output_prefix", PROJECT_NAME)
     
     do i = startTimestep, endTimestep, saveInterval

        write(filename, '(2A,I8.8,A)') trim(outputPrefix), "-", i, ".q"
        call region%loadData(QOI_FORWARD_STATE, filename)

        ! Compute velocity
        do j = 1, size(region%grids)
           call computeDependentVariables(size(globalGridSizes, 1),                          &
                region%states(j)%conservedVariables, velocity = region%states(j)%velocity)
        end do
     
        ! Save vorticity and dilatation
        write(filename, '(2A,I8.8,A)') trim(outputPrefix), "-", i, "..vorticity_dilatation.f"
        call region%saveData(QOI_VORTICITY_DILATATION, filename)        
        
     end do

  end if

  ! Cleanup.
  call region%cleanup()

  ! Finalize MPI.
  call MPI_Finalize(ierror)

end program
