#include "config.h"

program probe_from_qfile

  use MPI

  use Grid_enum, only : QOI_GRID, QOI_TARGET_MOLLIFIER
  use State_enum, only : QOI_FORWARD_STATE
  use Region_enum, only : FORWARD, ADJOINT

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

  ! Reset probes.
  call region%resetProbes()

  if (command_argument_count() >= 1) then !... only one solution file to process.

     ! Load the solution file.
     call get_command_argument(1, filename)
     call region%loadData(QOI_FORWARD_STATE, filename)

     ! Save solution on probe patches.
     call region%saveProbeData(FORWARD)

     write(message,'(A)') 'probe data collected.'
     call writeAndFlush(region%comm, output_unit, message)

  else

     call getRequiredOption("probe_from_qfile/qfile_interval", saveInterval)
     call getRequiredOption("probe_from_qfile/start_timestep", startTimestep)
     call getRequiredOption("probe_from_qfile/end_timestep", endTimestep)

     outputPrefix = getOption("output_prefix", PROJECT_NAME)

     do i = startTimestep, endTimestep, saveInterval

        write(filename, '(2A,I8.8,A)') trim(outputPrefix), "-", i, ".q"
        call region%loadData(QOI_FORWARD_STATE, filename)

        ! Save solution on probe patches.
        call region%saveProbeData(FORWARD)

        write(message,'(I8.8,A)') i,'-th step probe data collected.'
        call writeAndFlush(region%comm, output_unit, message)

     end do

  end if

  ! Finish writing remaining data gathered on probes.
  call region%saveProbeData(FORWARD, finish = .true.)

  write(message,'(A)') 'Data collection is finished.'
  call writeAndFlush(region%comm, output_unit, message)

  ! Cleanup.
  call region%cleanup()

  ! Finalize MPI.
  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

end program probe_from_qfile
