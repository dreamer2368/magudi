#include "config.h"

program spectral_radius

  use MPI

  use Grid_enum, only : QOI_GRID
  use State_enum, only : QOI_FORWARD_STATE

  use Region_mod, only : t_Region

  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush, gracefulExit
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage

  implicit none

  integer :: i, startTimestep, endTimestep, saveInterval, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix
  logical :: success
  type(t_Region) :: region
  integer, allocatable :: globalGridSizes(:,:)

  interface

     subroutine saveSpectralRadius(region, filename)

       use Region_mod, only : t_Region

       class(t_Region) :: region
       character(len = *), intent(in) :: filename

     end subroutine saveSpectralRadius

  end interface

  ! Initialize MPI.
  call MPI_Init(ierror)

  ! Parse options from the input file.
  filename = "spectral_radius.inp"
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

     ! Save vorticity and dilatation
     i = len_trim(filename)
     if (filename(i-1:i) == ".q") then
        filename = filename(:i-2) // ".spectral_radius.f"
     else
        filename = PROJECT_NAME // ".spectral_radius.f"
     end if
     call saveSpectralRadius(region, filename)

  else

     call getRequiredOption("save_interval", saveInterval)
     call getRequiredOption("start_timestep", startTimestep)
     call getRequiredOption("end_timestep", endTimestep)

     outputPrefix = getOption("output_prefix", PROJECT_NAME)

     do i = startTimestep, endTimestep, saveInterval

        write(filename, '(2A,I8.8,A)') trim(outputPrefix), "-", i, ".q"
        call region%loadData(QOI_FORWARD_STATE, filename)

        write(filename, '(2A,I8.8,A)') trim(outputPrefix), "-", i, ".spectral_radius.f"
        call saveSpectralRadius(region, filename)

     end do

  end if

  ! Cleanup.
  call region%cleanup()

  ! Finalize MPI.
  call MPI_Finalize(ierror)

end program spectral_radius

subroutine saveSpectralRadius(region, filename)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region

  ! <<< Enumerations >>>
  use State_enum, only : QOI_DUMMY_FUNCTION

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables, computeSpectralRadius

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region
  character(len = *), intent(in) :: filename

  ! <<< Local variables >>>
  type :: t_SpectralRadiusInternal
     SCALAR_TYPE, pointer :: buffer(:,:) => null()
  end type t_SpectralRadiusInternal
  type(t_SpectralRadiusInternal), allocatable, save :: data_(:)
  integer, save :: nDimensions = 0
  integer :: i

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
        allocate(data_(i)%buffer(region%grids(i)%nGridPoints, nDimensions))
     end do
  end if

  do i = 1, size(region%states)

     if (.not. allocated(region%states(i)%velocity))                                         &
          allocate(region%states(i)%velocity(region%grids(i)%nGridPoints, nDimensions))
     if (.not. allocated(region%states(i)%temperature))                                      &
          allocate(region%states(i)%temperature(region%grids(i)%nGridPoints, 1))


     call computeDependentVariables(nDimensions, 0, region%states(i)%conservedVariables,     &
          velocity = region%states(i)%velocity,                                              &
          temperature = region%states(i)%temperature(:,1))
     call computeSpectralRadius(nDimensions, region%solverOptions%ratioOfSpecificHeats,      &
          region%states(i)%specificVolume(:,1), region%states(i)%velocity,                   &
          region%states(i)%pressure(:,1), region%grids(i)%metrics,                           &
          spectralRadius = data_(i)%buffer)

     region%states(i)%dummyFunction => data_(i)%buffer

  end do

  call region%saveData(QOI_DUMMY_FUNCTION, filename)

end subroutine saveSpectralRadius
