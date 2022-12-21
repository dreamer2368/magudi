#include "config.h"

program rhs

  use MPI

  use Grid_enum, only : QOI_GRID
  use State_enum, only : QOI_FORWARD_STATE

  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver

  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush, gracefulExit
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use Patch_factory, only : computeSpongeStrengths, updatePatchFactories
  use InterfaceHelper, only : checkFunctionContinuityAtInterfaces

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, startTimestep, endTimestep, saveInterval, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix
  logical :: success
  type(t_Region) :: region
  type(t_Solver) :: solver
  integer, allocatable :: globalGridSizes(:,:)

  interface

     subroutine saveRhs(region, filename)

       use Region_mod, only : t_Region

       class(t_Region) :: region
       character(len = *), intent(in) :: filename

     end subroutine saveRhs

  end interface

  ! Initialize MPI.
  call MPI_Init(ierror)

  call initializeErrorHandler()

  ! Parse options from the input file.
  filename = PROJECT_NAME // ".inp"
  call parseInputFile(filename)

  outputPrefix = getOption("output_prefix", PROJECT_NAME)

  ! Verify that the grid file is in valid PLOT3D format and fetch the grid dimensions:
  ! `globalGridSizes(i,j)` is the number of grid points on grid `j` along dimension `i`.
  call getRequiredOption("grid_file", filename)
  call plot3dDetectFormat(MPI_COMM_WORLD, filename,                                          &
       success, globalGridSizes = globalGridSizes)
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

  ! Initialize the solver.
  call solver%setup(region, outputPrefix = outputPrefix)

  if (command_argument_count() >= 1) then !... only one solution file to process.

     ! Load the solution file.
     call get_command_argument(1, filename)
     call region%loadData(QOI_FORWARD_STATE, filename)

     ! Save RHS
     i = len_trim(filename)
     if (filename(i-1:i) == ".q") then
        filename = filename(:i-2)
     else
        filename = PROJECT_NAME
     end if
     call saveRhs(region, filename)

  else

     call getRequiredOption("rhs/save_interval", saveInterval)
     call getRequiredOption("rhs/start_timestep", startTimestep)
     call getRequiredOption("rhs/end_timestep", endTimestep)

     outputPrefix = getOption("output_prefix", PROJECT_NAME)

     do i = startTimestep, endTimestep, saveInterval

        write(filename, '(2A,I8.8,A)') trim(outputPrefix), "-", i, ".q"
        call region%loadData(QOI_FORWARD_STATE, filename)

        write(filename, '(2A,I8.8)') trim(outputPrefix), "-", i
        call saveRhs(region, filename)

     end do

  end if

  ! Cleanup.
  call region%cleanup()

  ! Finalize MPI.
  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

end program rhs

subroutine saveRhs(region, filename)

  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Patch_mod, only : t_Patch

  ! <<< Enumerations >>>
  use State_enum, only : QOI_DUMMY_FUNCTION
  use Region_enum, only : FORWARD

  ! <<< External modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region
  character(len = *), intent(in) :: filename

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer, save :: nDimensions = 0
  integer :: i, j, maxPatchSize, ierror
  character(len = STRING_LENGTH) :: patchFile
  type :: t_RhsInternal
     SCALAR_TYPE, pointer :: buffer(:,:) => null()
  end type t_RhsInternal
  type(t_RhsInternal), allocatable, save :: data_(:)
  class(t_Patch), pointer :: patch => null()

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
        allocate(data_(i)%buffer(region%grids(i)%nGridPoints, region%solverOptions%nUnknowns))
     end do
  end if

  call region%computeRhs(FORWARD, region%timestep, 1)

  do i = 1, size(region%states)
     data_(i)%buffer = region%states(i)%rightHandSide
     region%states(i)%dummyFunction => data_(i)%buffer
  end do

  call region%saveData(QOI_DUMMY_FUNCTION, trim(filename) // ".rhs.f")

  if (getOption("rhs/save_patch_rhs", .false.) .and. (allocated(region%patchFactories))) then
    maxPatchSize = size(region%patchFactories)
    call MPI_Allreduce(MPI_IN_PLACE, maxPatchSize, 1,                                           &
                       MPI_INTEGER, MPI_MAX, region%comm, ierror)

    ! Add patch penalties.
    do i = 1, maxPatchSize
      ! Zero out the right-hand side.
      do j = 1, size(region%states)
        region%states(j)%rightHandSide = 0.0_wp
      end do

      if (i .le. size(region%patchFactories)) then
        call region%patchFactories(i)%connect(patch)
        if (associated(patch)) then
          do j = 1, size(region%states)
            if (patch%gridIndex /= region%grids(j)%index) cycle
            call patch%updateRhs(FORWARD, region%simulationFlags, region%solverOptions,              &
                                 region%grids(j), region%states(j))
          end do
        end if
      end if

      do j = 1, size(region%states)
         data_(j)%buffer = region%states(j)%rightHandSide
         region%states(j)%dummyFunction => data_(j)%buffer
      end do

      write(patchFile, '(3A)') trim(filename) // ".", trim(patch%name), ".f"
      call region%saveData(QOI_DUMMY_FUNCTION, patchFile)
    end do
  end if

  if (getOption("enable_immersed_boundary", .false.)) then
    ! resize the data buffer.
    do i = 1, size(data_)
      deallocate(data_(i)%buffer)
      allocate(data_(i)%buffer(region%grids(i)%nGridPoints, 1 + 2 * region%grids(i)%nDimensions))
    end do

    do j = 1, size(region%states)
       data_(j)%buffer(:, 1) = region%states(j)%levelset(:, 1)
       data_(j)%buffer(:, 2:region%grids(j)%nDimensions + 1) = region%states(j)%levelsetNormal
       data_(j)%buffer(:, region%grids(j)%nDimensions + 2 : 2 * region%grids(j)%nDimensions + 1) = region%states(j)%objectVelocity
       region%states(j)%dummyFunction => data_(j)%buffer
    end do

    write(patchFile, '(A)') trim(filename) // ".levelset.f"
    call region%saveData(QOI_DUMMY_FUNCTION, patchFile)
  end if

end subroutine saveRhs
