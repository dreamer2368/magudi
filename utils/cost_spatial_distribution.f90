#include "config.h"

program cost_spatial_distribution

  use MPI

  use Grid_enum, only : QOI_GRID
  use State_enum, only : QOI_FORWARD_STATE

  use Region_mod, only : t_Region
  use Functional_factory, only : t_FunctionalFactory
  use Functional_mod, only : t_Functional

  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush, gracefulExit
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage

  implicit none

  integer :: i, j, startTimestep, endTimestep, saveInterval, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix
  logical :: success
  type(t_Region) :: region
  type(t_FunctionalFactory) :: functionalFactory
  class(t_Functional), pointer :: functional => null()
  integer, allocatable :: globalGridSizes(:,:)

  interface

     subroutine saveSpatialDistribution(region, functionalFactory, filename)

       use Region_mod, only : t_Region
       use Functional_factory, only : t_FunctionalFactory

       class(t_Region) :: region
       type(t_FunctionalFactory) :: functionalFactory
       character(len = *), intent(in) :: filename

     end subroutine saveSpatialDistribution

  end interface

  ! Initialize MPI.
  call MPI_Init(ierror)

  ! Parse options from the input file.
  filename = PROJECT_NAME // ".inp"
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

  ! Setup the functional.
  call functionalFactory%connect(functional,                                                    &
                                 trim(region%solverOptions%costFunctionalType))
  assert(associated(functional))
  call functional%setup(region)

  if (command_argument_count() >= 1) then !... only one solution file to process.

     ! Load the solution file.
     call get_command_argument(1, filename)
     call region%loadData(QOI_FORWARD_STATE, filename)

     do i = 1, size(region%states) !... update state
        call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
             region%solverOptions)
     end do

     ! Save vorticity and dilatation
     i = len_trim(filename)
     if (filename(i-1:i) == ".q") then
        filename = filename(:i-2) // ".cost_spatial_distribution.f"
     else
        filename = PROJECT_NAME // ".cost_spatial_distribution.f"
     end if
     call saveSpatialDistribution(region, functionalFactory, filename)

  else

     call getRequiredOption("cost_spatial_distribution/save_interval", saveInterval)
     call getRequiredOption("cost_spatial_distribution/start_timestep", startTimestep)
     call getRequiredOption("cost_spatial_distribution/end_timestep", endTimestep)

     outputPrefix = getOption("output_prefix", PROJECT_NAME)

     do i = startTimestep, endTimestep, saveInterval

        write(filename, '(2A,I8.8,A)') trim(outputPrefix), "-", i, ".q"
        call region%loadData(QOI_FORWARD_STATE, filename)

        do j = 1, size(region%states) !... update state
           call region%states(j)%update(region%grids(j), region%simulationFlags,                   &
                region%solverOptions)
        end do

        write(filename, '(2A,I8.8,A)') trim(outputPrefix), "-", i, ".cost_spatial_distribution.f"
        call saveSpatialDistribution(region, functionalFactory, filename)

     end do

  end if

  ! Cleanup.
  call region%cleanup()

  ! Finalize MPI.
  call MPI_Finalize(ierror)

end program cost_spatial_distribution

subroutine saveSpatialDistribution(region, functionalFactory, filename)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Functional_factory, only : t_FunctionalFactory
  use Functional_mod, only : t_Functional

  ! <<< Enumerations >>>
  use State_enum, only : QOI_DUMMY_FUNCTION

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables, computeVorticityMagnitudeAndDilatation

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region
  type(t_FunctionalFactory) :: functionalFactory
  character(len = *), intent(in) :: filename

  ! <<< Local variables >>>
  type :: t_CostSpatialDistributionInternal
     SCALAR_TYPE, pointer :: buffer(:,:) => null()
  end type t_CostSpatialDistributionInternal
  class(t_Functional), pointer :: functional => null()
  type(t_CostSpatialDistributionInternal), allocatable, save :: data_(:)
  integer, save :: nDimensions = 0
  integer :: i

  call functionalFactory%connect(functional)
  assert(associated(functional))

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(size(region%grids) == size(region%states))

  assert(allocated(region%globalGridSizes))

  if (allocated(data_)) then
     if (size(data_) /= size(region%grids) .or.                                              &
          size(region%globalGridSizes, 1) /= nDimensions) then
        do i = 1, size(data_)
           if (associated(data_(i)%buffer)) deallocate(data_(i)%buffer)
           nullify(data_(i)%buffer)
        end do
        SAFE_DEALLOCATE(data_)
     end if
  end if

  nDimensions = size(region%globalGridSizes, 1)
  assert_key(nDimensions, (1, 2, 3))

  if (.not. allocated(data_)) then
     allocate(data_(size(region%grids)))
     do i = 1, size(data_)
        allocate(data_(i)%buffer(region%grids(i)%nGridPoints, 1))
     end do
  end if

  do i = 1, size(region%states)

     call functional%computeSpatialDistribution(region%grids(i),                              &
                                 region%states(i),data_(i)%buffer)

     region%states(i)%dummyFunction => data_(i)%buffer

  end do

  call region%saveData(QOI_DUMMY_FUNCTION, filename)

end subroutine saveSpatialDistribution
