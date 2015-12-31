#include "config.h"

program init_flow

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env

  ! <<< Enumerations >>>
  use Grid_enum
  use State_enum

  ! <<< Derived types >>>
  use Region_mod, only : t_Region

  ! <<< Internal modules >>>
  use InputHelper, only : getInputName, parseInputFile, getOption, getRequiredOption

  implicit none

  ! <<< Local variables >>>
  type(t_Region) :: region
  character(len = STRING_LENGTH) :: filename, key
  integer, allocatable :: globalGridSizes(:,:)
  integer :: i, j, procRank, numProcs, numGrids, nDimensions, ierror

  interface

     subroutine getSimulation(region)

       use Region_mod, only : t_Region

       class(t_Region) :: region

     end subroutine getSimulation

  end interface

  ! Initialize MPI.
  call MPI_Init(ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  ! Parse options from the input file.
  call getInputName(filename)
  call parseInputFile(filename)

  ! Get the grid dimensions.
  ! `globalGridSizes(i,j)` is the number of grid points on grid `j` along dimension `i`.
  numGrids = getOption("number_of_grids", 1)
  nDimensions = 0
  do j = 1, numGrids
     do i = 1, 3
        write(key, '(A,I3.3,A,I1.1,A)') "grid", j, "/dir", i, "/size"
        if (getOption(trim(key), 1) > 1) nDimensions = nDimensions + 1
     end do
  end do
  allocate(globalGridSizes(nDimensions, numGrids))
  do j = 1, numGrids
     do i = 1, nDimensions
        write(key, '(A,I3.3,A,I1.1,A)') "grid", j, "/dir", i, "/size"
        call getRequiredOption(trim(key), globalGridSizes(i, j))
     end do
  end do

  ! Setup the region.
  call region%setup(MPI_COMM_WORLD, globalGridSizes)
  deallocate(globalGridSizes)

  ! Initialize the simulation.
  call getSimulation(region)

  ! Save the grid.
  call getRequiredOption("grid_file", filename)
  call region%saveData(QOI_GRID, filename)

  ! Save initial condition.
  call getRequiredOption("grid_file", filename)
  i = len_trim(filename)
  if (filename(i-3:i) == ".xyz") then
     filename = filename(:i-4) // ".ic.q"
  else
     filename = PROJECT_NAME // ".ic.q"
  end if
  call region%saveData(QOI_FORWARD_STATE, filename)

  ! Cleanup.
  call region%cleanup()

  ! Finalize MPI.
  call MPI_Finalize(ierror)

end program init_flow

subroutine getSimulation(region)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region

  ! <<< Local variables >>>
  character(len = STRING_LENGTH) :: simulation, message

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(size(region%grids) == size(region%states))

  assert(allocated(region%globalGridSizes))

  ! Detect the simulation type.
  call getRequiredOption("simulation_name", simulation, MPI_COMM_WORLD)

  select case (trim(simulation))

  case('COUNTERFLOW_DIFFUSION')
     call counterflowDiffusionGrid(region)
     call counterflowDiffusionQ(region)
     call counterflowDiffusionTargetMollifier(region)
     call counterflowDiffusionTargetMollifier(region)
     call counterflowDiffusionBC

  case('JET')
     call jetGrid(region)
     call jetQ(region)
     call jetTargetMollifier(region)
     call jetControlMollifier(region)
     call jetBC

  case('JET_IN_CROSSFLOW')
     call jetCrossFlowGrid(region)
     call jetCrossFlowQ(region)
     call jetCrossFlowTargetMollifier(region)
     call jetCrossFlowControlMollifier(region)
     call jetCrossFlowBC

  case('MIXING_LAYER')
     call mixingLayerGrid(region)
     call mixingLayerQ(region)
     call mixingLayerTargetMollifier(region)
     call mixingLayerControlMollifier(region)
     call mixingLayerBC

  case('PREMIXED')
     call premixedGrid(region)
     call premixedQ(region)
     call premixedTargetMollifier(region)
     call premixedControlMollifier(region)
     call premixedBC

  case default
     write(message, '(3A)') "Unknown simulation '", trim(simulation), "'"
     call gracefulExit(MPI_COMM_WORLD, message)
  end select

end subroutine getSimulation
