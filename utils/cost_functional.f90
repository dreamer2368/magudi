#include "config.h"

program cost_functional

  use MPI

  use Grid_enum, only : QOI_GRID, QOI_TARGET_MOLLIFIER
  use State_enum, only : QOI_FORWARD_STATE

  use Region_mod, only : t_Region
  use Functional_mod, only : t_Functional
  use Functional_factory, only : t_FunctionalFactory

  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush, gracefulExit
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use Patch_factory, only : computeSpongeStrengths, updatePatchFactories
  use InterfaceHelper, only : checkFunctionContinuityAtInterfaces 

  use, intrinsic :: iso_fortran_env, only : output_unit
  use CNSHelper, only : computeDependentVariables

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, startTimestep, endTimestep, saveInterval, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix
  logical :: success
  type(t_Region) :: region
  integer, allocatable :: globalGridSizes(:,:)
  type(t_FunctionalFactory) :: functionalFactory
  class(t_Functional), pointer :: functional => null()
  SCALAR_TYPE :: instantaneousFunctional
  character(len = STRING_LENGTH) :: message  !SeungWhan: for debugging

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
  call plot3dDetectFormat(MPI_COMM_WORLD, filename,                                         &
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

  ! SeungWhan: allocate Target Mollifier
  do i = 1, size(region%grids)
     if ( .not. allocated(region%grids(i)%targetMollifier) ) then
        allocate(region%grids(i)%targetMollifier(region%grids(i)%nGridPoints, 1))
     end if
  end do

  ! SeungWhan: Load Target mollifier.
  filename = getOption("target_mollifier_file", "")
  if (len_trim(filename) == 0) then
     do i = 1, size(region%grids)
        assert(allocated(region%grids(i)%targetMollifier))
        region%grids(i)%targetMollifier = 1.0_wp
     end do
  else
     call region%loadData(QOI_TARGET_MOLLIFIER, filename)
  end if

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

  ! Compute normalized metrics, norm matrix and Jacobian.
  do i = 1, size(region%grids)
     call region%grids(i)%update()
  end do
  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  if (.not. region%simulationFlags%predictionOnly) then
     call functionalFactory%connect(functional,                                             &
          trim(region%solverOptions%costFunctionalType))
     assert(associated(functional))
     call functional%setup(region)
  end if

  if (command_argument_count() >= 1) then !... only one solution file to process.

     ! Load the solution file.
     call get_command_argument(1, filename)
     call region%loadData(QOI_FORWARD_STATE, filename)
     do j = 1, size(region%states)
        if (.not. allocated(region%states(j)%pressure))                                         &
           allocate(region%states(j)%pressure(region%grids(j)%nGridPoints, 1))
        call computeDependentVariables(region%grids(j)%nDimensions, region%states(j)%conservedVariables,        &
                                        pressure = region%states(j)%pressure(:,1))
     end do

     instantaneousFunctional = 0.0_wp
     if (.not. region%simulationFlags%predictionOnly) then
        instantaneousFunctional = functional%compute(region)
!SeungWhan
write(message,*) instantaneousFunctional
call writeAndFlush(region%comm, output_unit, message, advance = 'no')
!=========
        call functional%writeToFile(region%comm, trim(outputPrefix) //                      &
             ".cost_functional.txt", region%timestep, region%states(1)%time,                &
             .false.)
     end if

  else

     call getRequiredOption("cost_functional/save_interval", saveInterval)
     call getRequiredOption("cost_functional/start_timestep", startTimestep)
     call getRequiredOption("cost_functional/end_timestep", endTimestep)

     outputPrefix = getOption("output_prefix", PROJECT_NAME)

     do i = startTimestep, endTimestep, saveInterval

        write(filename, '(2A,I8.8,A)') trim(outputPrefix), "-", i, ".q"
        call region%loadData(QOI_FORWARD_STATE, filename)
        do j = 1, size(region%states)
           if (.not. allocated(region%states(j)%pressure))                                         &
              allocate(region%states(j)%pressure(region%grids(j)%nGridPoints, 1))
           call computeDependentVariables(region%grids(j)%nDimensions, region%states(j)%conservedVariables,        &
                                           pressure = region%states(j)%pressure(:,1))
        end do

        instantaneousFunctional = 0.0_wp
        if (.not. region%simulationFlags%predictionOnly) then
           instantaneousFunctional = functional%compute(region)
!SeungWhan
write(message,*) i,'-th step: ',instantaneousFunctional
call writeAndFlush(region%comm, output_unit, message, advance = 'no')
!=========
           call functional%writeToFile(region%comm, trim(outputPrefix) //                   &
                ".cost_functional.txt", i, region%states(1)%time,                           &
                i > startTimestep)
        end if

     end do

  end if

  ! Cleanup.
  call region%cleanup()

  ! Finalize MPI.
  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

end program cost_functional
