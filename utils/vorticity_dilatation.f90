#include "config.h"

program vorticity_dilatation

  use MPI
  use, intrinsic :: iso_fortran_env

  use Grid_enum
  use State_enum
  use Region_type

  use CNSHelper, only : computeDependentVariables
  use Region_mod, only : setupRegion, cleanupRegion, loadRegionData, saveRegionData
  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush, gracefulExit
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage

  implicit none

  integer :: i, procRank, ierror
  character(len = STRING_LENGTH) :: filename, message, programName
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
  call setupRegion(region, MPI_COMM_WORLD, globalGridSizes)
  call loadRegionData(region, QOI_GRID, filename)

  ! Load the solution file.
  if (command_argument_count() >= 1) then
     call get_command_argument(1, filename)
  else
     call getRequiredOption("solution_file", filename)
  end if
  call loadRegionData(region, QOI_FORWARD_STATE, filename)

  ! Setup spatial discretization.
  do i = 1, size(region%grids)
     call region%grids(i)%setupSpatialDiscretization()
  end do
  call MPI_Barrier(region%comm, ierror)

  ! Compute normalized metrics, norm matrix and Jacobian.
  do i = 1, size(region%grids)
     call region%grids(i)%update()
  end do
  call MPI_Barrier(MPI_COMM_WORLD, ierror)

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
  call saveRegionData(region, QOI_VORTICITY_DILATATION, filename)

  ! Cleanup.
  call cleanupRegion(region)

  ! Finalize MPI.
  call MPI_Finalize(ierror)

end program
