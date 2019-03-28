#include "config.h"

program grid_generator

  use MPI
  use, intrinsic :: iso_fortran_env

  use Grid_enum
  use State_enum

  use Region_mod, only : t_Region
  use MappingFunction_mod, only : t_MappingFunction
  use MappingFunction_factory, only : t_MappingFunctionFactory

  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush, gracefulExit

  !> Generates a grid and writes it to a PLOT3D multi-block whole-format grid file.

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, direction, nBlocks, nDimensions, gridSize(3), offset(3), ierror
  character(len = STRING_LENGTH) :: filename, message, key, mappingFunctionType, outputPrefix
  type(t_Region) :: region
  integer, allocatable :: globalGridSizes(:,:)
  real(wp), allocatable :: coordinate(:)
  type(t_MappingFunctionFactory) :: mappingFunctionFactory
  class(t_MappingFunction), pointer :: mappingFunction => null()
  real(wp) :: minMaxRange(2)

  ! Initialize MPI.
  call MPI_Init(ierror)

  ! Parse options from the input file.
  filename = PROJECT_NAME // ".inp"
  call parseInputFile(filename)

  nBlocks = getOption("number_of_blocks", 1)
  if (nBlocks < 1 .or. nBlocks > 999) then
     write(message, '(2A)') trim(filename), ": number of blocks is invalid (range: 1-999)!"
     call gracefulExit(MPI_COMM_WORLD, message)
  end if

  call getRequiredOption("number_of_dimensions", nDimensions)
  if (nDimensions < 1 .or. nDimensions > 3) then
     write(message, '(2A)') trim(filename), ": number of dimensions is invalid (range: 1-3)!"
     call gracefulExit(MPI_COMM_WORLD, message)
  end if

  allocate(globalGridSizes(nDimensions, nBlocks), source = 0)

  do i = 1, nBlocks
     do j = 1, nDimensions
        write(key, '(A,I3.3,A,I1.1,A)') "grid", i, "/dir", j, "/number_of_points"
        call getRequiredOption(trim(key), globalGridSizes(j,i))
     end do
  end do

  ! Setup the region and load the grid file.
  call region%setup(MPI_COMM_WORLD, globalGridSizes)

  do l = 1, size(region%grids)
     do direction = 1, nDimensions

        write(key, '(A,I3.3,A,I1.1,A)') "grid", l, "/dir", direction, "/"

        mappingFunctionType = getOption(trim(key) // "mapping_function", "UNIFORM")
        call mappingFunctionFactory%connect(mappingFunction,                                 &
             trim(mappingFunctionType), createNew = .true.)
        if (.not. associated(mappingFunction)) then
           write(message, '(4A)') trim(filename), ": invalid mapping function '",            &
                trim(mappingFunctionType), "'!"
           call gracefulExit(region%grids(l)%comm, message)
        end if

        minMaxRange(1) = getOption(trim(key) //                                              &
             achar(iachar('x') + direction - 1) // "_min", 0.0_wp)
        if (region%grids(1)%periodicityType(direction) == PLANE) then
           minMaxRange(2) = minMaxRange(1) + region%grids(1)%periodicLength(direction)
        else
           minMaxRange(2) = getOption(trim(key) //                                           &
                achar(iachar('x') + direction - 1) // "_max", 1.0_wp)
        end if

        allocate(coordinate(region%grids(l)%globalSize(direction)))
        call mappingFunction%compute(coordinate, region%grids(l)%index, direction,           &
             minMaxRange, region%grids(1)%periodicityType(direction) == PLANE)

        gridSize = region%grids(l)%localSize
        offset = region%grids(l)%offset

        select case (direction)

        case (1)
           do k = 1, gridSize(3)
              do j = 1, gridSize(2)
                 do i = 1, gridSize(1)
                    region%grids(l)%coordinates(i + gridSize(1) * (j - 1 +                   &
                         gridSize(2) * (k - 1)), direction) = coordinate(i + offset(direction))
                 end do
              end do
           end do

        case (2)
           do k = 1, gridSize(3)
              do j = 1, gridSize(2)
                 do i = 1, gridSize(1)
                    region%grids(l)%coordinates(i + gridSize(1) * (j - 1 +                   &
                         gridSize(2) * (k - 1)), direction) = coordinate(j + offset(direction))
                 end do
              end do
           end do

        case (3)
           do k = 1, gridSize(3)
              do j = 1, gridSize(2)
                 do i = 1, gridSize(1)
                    region%grids(l)%coordinates(i + gridSize(1) * (j - 1 +                   &
                         gridSize(2) * (k - 1)), direction) = coordinate(k + offset(direction))
                 end do
              end do
           end do

        end select

        SAFE_DEALLOCATE(coordinate)
        call mappingFunctionFactory%cleanup()

     end do
  end do

  ! Compute normalized metrics, norm matrix and Jacobian.
  do i = 1, size(region%grids)
     call region%grids(i)%update()
  end do
  call MPI_Barrier(region%comm, ierror)

  ! Write out some useful information.
  call region%reportGridDiagnostics()

  ! Save the grid.
  outputPrefix = getOption("output_prefix", PROJECT_NAME)
  call region%saveData(QOI_GRID, trim(outputPrefix) // ".xyz")

  ! Cleanup.
  call region%cleanup()

  ! Finalize MPI.
  call MPI_Finalize(ierror)

end program grid_generator
