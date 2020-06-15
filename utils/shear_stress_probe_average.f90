#include "config.h"

program shear_stress_probe_average

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

     subroutine saveShearStressOnProbe(region)

       use Region_mod, only : t_Region

       class(t_Region) :: region

     end subroutine saveShearStressOnProbe

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

  call saveShearStressOnProbe(region)

  ! Cleanup.
  call region%cleanup()

  ! Finalize MPI.
  call MPI_Finalize(ierror)

end program shear_stress_probe_average

subroutine saveShearStressOnProbe(region)

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

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  class(t_Patch), pointer :: patch => null()
  integer, save :: nDimensions = 0
  integer :: i, j, k, l, m, ierror
  SCALAR_TYPE, allocatable :: mask(:), F(:,:)
  SCALAR_TYPE :: volume, density, shearStress
  character(len=STRING_LENGTH) :: message

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(size(region%grids) == size(region%states))

  assert(allocated(region%globalGridSizes))

  nDimensions = size(region%globalGridSizes, 1)
  assert_key(nDimensions, (1, 2, 3))
  assert((component(1)>=1).and.(component(1)<=nDimensions))
  assert((component(2)>=1).and.(component(2)<=nDimensions))

  do i = 1, size(region%states)
    call region%states(i)%update(region%grids(i), region%simulationFlags,       &
                                 region%solverOptions)
  end do

  volume = 0.0_wp
  shearStress = 0.0_wp
  density = 0.0_wp
  do l = 1, size(region%patchFactories)
     call region%patchFactories(l)%connect(patch)
     if (.not. associated(patch)) cycle
     do m = 1, size(region%states)
        if (patch%gridIndex /= region%grids(m)%index .or. patch%nPatchPoints <= 0) cycle
        select type (patch)
        class is (t_ProbePatch)
          allocate(mask(region%grids(m)%nGridPoints))
          allocate(F(region%grids(m)%nGridPoints,nDimensions))
          mask = 0.0_wp
          F = 0.0_wp

          assert(all(patch%gridLocalSize == region%grids(m)%localSize))
          assert(all(patch%gridOffset == region%grids(m)%offset))

          do k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)
             do j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
                do i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
                   mask(i - region%grids(m)%offset(1) + region%grids(m)%localSize(1) * (j - 1 - region%grids(m)%offset(2) +     &
                        region%grids(m)%localSize(2) * (k - 1 - region%grids(m)%offset(3)))) = 1.0_wp
                end do
             end do
          end do

          where (region%grids(m)%iblank == 0)
             mask = 0.0_wp
          end where

          F = region%states(m)%stressTensor(:,                                        &
          1+(patch%normalDirection-1)*nDimensions:patch%normalDirection*nDimensions)
          F(:,patch%normalDirection) = 0.0_wp

          volume = volume + region%grids(m)%computeInnerProduct(mask,mask)
          shearStress = shearStress + region%grids(m)%computeInnerProduct(mask,sqrt(sum(F**2,dim=2)))
          density = density + region%grids(m)%computeInnerProduct(mask,region%states(m)%conservedVariables(:,1))

          SAFE_DEALLOCATE(mask)
          SAFE_DEALLOCATE(F)
        end select
     end do
  end do

  if (region%commGridMasters /= MPI_COMM_NULL) then
    call MPI_Allreduce(MPI_IN_PLACE, volume, 1,                                 &
    SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)
    call MPI_Allreduce(MPI_IN_PLACE, shearStress, 1,                            &
    SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)
    call MPI_Allreduce(MPI_IN_PLACE, density, 1,                                &
    SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)
  end if

  do i = 1, size(region%grids)
    call MPI_Bcast(volume, 1, SCALAR_TYPE_MPI, 0, region%grids(i)%comm, ierror)
    call MPI_Bcast(shearStress, 1, SCALAR_TYPE_MPI, 0, region%grids(i)%comm, ierror)
    call MPI_Bcast(density, 1, SCALAR_TYPE_MPI, 0, region%grids(i)%comm, ierror)
  end do

  write(message,'(A,1X,'// SCALAR_FORMAT //')') 'volume: ', volume
  call writeAndFlush(region%comm, output_unit, message)
  write(message,'(A,1X,'// SCALAR_FORMAT //')') 'shear stress: ', shearStress
  call writeAndFlush(region%comm, output_unit, message)
  write(message,'(A,1X,'// SCALAR_FORMAT //')') 'density: ', density
  call writeAndFlush(region%comm, output_unit, message)

end subroutine saveShearStressOnProbe
