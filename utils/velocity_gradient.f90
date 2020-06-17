#include "config.h"

program velocity_gradient

  use MPI

  use Grid_enum, only : QOI_GRID
  use State_enum, only : QOI_FORWARD_STATE

  use Region_mod, only : t_Region

  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush, gracefulExit
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use Patch_factory, only : computeSpongeStrengths, updatePatchFactories
  use InterfaceHelper, only : checkFunctionContinuityAtInterfaces

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix, message
  logical :: success
  type(t_Region) :: region
  integer, allocatable :: globalGridSizes(:,:)

  interface

     subroutine saveVelocityGradient(region,filename)

       use Region_mod, only : t_Region

       class(t_Region) :: region
       character(len=*), intent(in) :: filename

     end subroutine saveVelocityGradient

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

  ! Setup boundary conditions.
  call getRequiredOption("boundary_condition_file", filename)
  call region%setupBoundaryConditions(filename)

  ! Compute damping strength on sponge patches.
  do i = 1, size(region%grids)
     call computeSpongeStrengths(region%patchFactories, region%grids(i))
  end do

  ! Check continuity at block interfaces.
  if (getOption("check_interface_continuity", .false.))                                     &
       call checkFunctionContinuityAtInterfaces(region, epsilon(0.0_wp))

  ! Update patches.
  do i = 1, size(region%grids)
     call updatePatchFactories(region%patchFactories, region%simulationFlags,               &
          region%solverOptions, region%grids(i), region%states(i))
  end do

  ! Reset probes.
  call region%resetProbes()

  if (command_argument_count() < 1) then !... only one solution file to process.
    write(message,'(A)') 'Requires a filename!'
    call gracefulExit(MPI_COMM_WORLD, message)
  end if

  ! Load the solution file.
  call get_command_argument(1, filename)
  call region%loadData(QOI_FORWARD_STATE, filename)

  i = len_trim(filename)
  if (filename(i-1:i) == ".q") then
     filename = filename(:i-2) // ".velocity_gradient.f"
  else
     filename = PROJECT_NAME // ".velocity_gradient.f"
  end if

  call saveVelocityGradient(region, filename)

  ! Cleanup.
  call region%cleanup()

  ! Finalize MPI.
  call MPI_Finalize(ierror)

end program velocity_gradient

subroutine saveVelocityGradient(region, filename)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use ProbePatch_mod, only : t_ProbePatch
  use Region_mod, only : t_Region

  ! <<< Enumerations >>>
  use State_enum, only : QOI_DUMMY_FUNCTION

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables, computeStressTensor
  use ErrorHandler, only : writeAndFlush

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region
  character(len = *), intent(in) :: filename

  ! <<< Local variables >>>
  type :: t_VelocityGradientInternal
     SCALAR_TYPE, pointer :: buffer(:,:) => null()
  end type t_VelocityGradientInternal
  type(t_VelocityGradientInternal), allocatable :: data_(:)
  integer, parameter :: wp = SCALAR_KIND
  class(t_Patch), pointer :: patch => null()
  integer, save :: nDimensions = 0
  integer :: i, j, k, l, m, ierror, index_, direction
  SCALAR_TYPE, allocatable :: mask(:), F(:,:), localStressTensor(:)
  SCALAR_TYPE :: volume, density, shearStress
  character(len=STRING_LENGTH) :: message

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(size(region%grids) == size(region%states))

  assert(allocated(region%globalGridSizes))

  nDimensions = size(region%globalGridSizes, 1)
  assert_key(nDimensions, (1, 2, 3))
  if (.not. allocated(data_)) then
     allocate(data_(size(region%grids)))
     do i = 1, size(data_)
        allocate(data_(i)%buffer(region%grids(i)%nGridPoints, nDimensions**2))
     end do
  end if

  do i = 1, size(region%states)
    call region%states(i)%update(region%grids(i), region%simulationFlags,       &
                                 region%solverOptions)
    call region%grids(i)%computeGradient(region%states(i)%velocity,             &
         region%states(i)%stressTensor)
    data_(i)%buffer = region%states(i)%stressTensor
    region%states(i)%dummyFunction => data_(i)%buffer
  end do

  call region%saveData(QOI_DUMMY_FUNCTION, filename)

end subroutine saveVelocityGradient
