#include "config.h"

program vorticity_dilatation

  use MPI
  use, intrinsic :: iso_fortran_env

  use Grid_type
  use State_type
  use Region_type

  use Grid_mod, only : setupSpatialDiscretization, updateGrid
  use CNSHelper, only : computeDependentVariables
  use MPIHelper, only : writeAndFlush
  use Region_mod, only : setupRegion, loadRegionData, saveRegionData
  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use PLOT3DHelper, only : plot3dDetectFormat, errorMessage

  implicit none

  integer :: i, procRank, ierror
  character(len = STRING_LENGTH) :: filename, message, programName
  logical :: success
  type(t_Region) :: region
  integer, allocatable :: globalGridSizes(:,:)

  ! Initialize MPI.
  call MPI_Init(ierror)

  if (command_argument_count() == 2) then

     ! Setup the region and load the grid file.
     call get_command_argument(1, filename)
     call plot3dDetectFormat(MPI_COMM_WORLD, filename, success,                              &
          globalGridSizes = globalGridSizes)
     if (.not. success) call gracefulExit(MPI_COMM_WORLD, errorMessage)
     call setupRegion(region, MPI_COMM_WORLD, globalGridSizes)
     call loadRegionData(region, QOI_GRID, filename)

     ! Load the solution file.
     call get_command_argument(2, filename)
     call loadRegionData(region, QOI_FORWARD_STATE, filename)

     ! Setup spatial discretization.
     do i = 1, size(region%grids)
        call setupSpatialDiscretization(region%grids(i))
     end do
     call MPI_Barrier(region%comm, ierror)

     ! Compute normalized metrics, norm matrix and Jacobian.
     do i = 1, size(region%grids)
        call updateGrid(region%grids(i))
     end do
     call MPI_Barrier(MPI_COMM_WORLD, ierror)

     ! Compute velocity
     do i = 1, size(region%grids)
        call computeDependentVariables(size(globalGridSizes, 1),                             &
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

  else

     call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
     call get_command_argument(0, programName)
     write(message, '(3A)') "Usage: ", trim(programName), " <grid file> <solution file>"
     call writeAndFlush(MPI_COMM_WORLD, error_unit, message)

  end if

  ! Finalize MPI.
  call MPI_Finalize(ierror)

end program
