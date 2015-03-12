#include "config.h"

program main

  use MPI

  use Grid_type
  use State_type
  use Region_type
  use RK4Integrator_type

  use Grid_mod, only : setupSpatialDiscretization, updateGrid
  use State_mod, only : updatePatches
  use Region_mod
  use MPIHelper, only : gracefulExit
  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use PLOT3DHelper, only : plot3dDetectFormat, errorMessage

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix
  logical :: success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_RK4Integrator) :: integrator

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
  if (.not. success) call gracefulExit(MPI_COMM_WORLD, errorMessage)

  ! Setup the region.
  call getRequiredOption("boundary_condition_file", filename)
  call setupRegion(region, MPI_COMM_WORLD, globalGridSizes, filename)

  ! Read the grid file.
  call getRequiredOption("grid_file", filename)
  call loadRegionData(region, QOI_GRID, filename)

  ! Read the target state.
  if (region%simulationFlags%useTargetState) then
     call getRequiredOption("target_state_file", filename)
     call loadRegionData(region, QOI_TARGET_STATE, filename)
  end if

  ! Setup spatial discretization.
  do i = 1, size(region%grids)
     call setupSpatialDiscretization(region%grids(i))
  end do
  call MPI_Barrier(region%comm, ierror)

  ! Update the grids by computing the Jacobian, metrics, and norm.
  do i = 1, size(region%grids)
     call updateGrid(region%grids(i))
     call updatePatches(region%states(i), region%grids(i),                                   &
          region%patches, region%simulationFlags, region%solverOptions)
  end do
  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  ! Write out some useful information.
  call reportGridDiagnostics(region)

  ! Save the Jacobian and normalized metrics.
  outputPrefix = getOption("output_prefix", PROJECT_NAME)
  write(filename, '(2A)') trim(outputPrefix), ".Jacobian.f"
  call saveRegionData(region, QOI_JACOBIAN, filename)
  write(filename, '(2A)') trim(outputPrefix), ".metrics.f"
  call saveRegionData(region, QOI_METRICS, filename)

  ! Setup the RK4 integrator.
  call setupRK4Integrator(integrator, region)

  ! Initial condition.
  if (region%simulationFlags%predictionOnly .and. command_argument_count() == 1) then
     call get_command_argument(1, filename)
  else
     call getRequiredOption("initial_condition_file", filename)
  end if

  if (region%simulationFlags%predictionOnly) then
     call loadRegionData(region, QOI_FORWARD_STATE, filename)
     call solveForward(region, integrator,                                                   &
          real(region%states(1)%plot3dAuxiliaryData(4), wp),                                 &
          nint(real(region%states(1)%plot3dAuxiliaryData(1), wp)),                           &
          getOption("number_of_timesteps", 1000), getOption("report_interval", 1),           &
          getOption("save_interval", 20), outputPrefix)
  end if

  call cleanupRK4Integrator(integrator)
  call cleanupRegion(region)

  ! Finalize MPI.
  call MPI_Finalize(ierror)

contains

  subroutine solveForward(region, integrator, startTime, startTimestep,                      &
       nTimesteps, reportInterval, saveInterval, outputPrefix)

    ! <<< Derived types >>>
    use Region_type
    use RK4Integrator_type

    ! <<< Internal modules >>>
    use RK4Integrator_mod, only : stepForward

    implicit none

    ! <<< Arguments >>>
    type(t_Region) :: region
    type(t_RK4Integrator) :: integrator
    real(SCALAR_KIND), intent(in) :: startTime
    integer, intent(in) :: startTimestep, nTimesteps
    integer, intent(in), optional :: reportInterval, saveInterval
    character(len = STRING_LENGTH), intent(in), optional :: outputPrefix

    ! <<< Local variables >>>
    character(len = STRING_LENGTH) :: outputPrefix_, filename
    real(SCALAR_KIND) :: time
    integer :: timestep
    logical :: verbose

    outputPrefix_ = PROJECT_NAME
    if (present(outputPrefix)) outputPrefix_ = outputPrefix

    time = startTime
    verbose = .false.

    do timestep = startTimestep + 1, startTimestep + nTimesteps

       if (present(reportInterval)) verbose = (reportInterval > 0 .and.                      &
            mod(timestep, reportInterval) == 0)
       call stepForward(integrator, region, time, timestep, verbose)

       if (present(saveInterval)) then
          if (saveInterval > 0) then
             if (mod(timestep, saveInterval) == 0) then
                write(filename, '(2A,I8.8,A)') trim(outputPrefix_), "-", timestep, ".q"
                call saveRegionData(region, QOI_FORWARD_STATE, filename)
             end if
          end if
       end if

    end do

  end subroutine solveForward

end program main
