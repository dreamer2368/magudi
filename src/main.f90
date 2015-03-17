#include "config.h"

program main

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use Grid_type
  use State_type
  use Region_type
  use RK4Integrator_type

  use Grid_mod, only : setupSpatialDiscretization, updateGrid, computeSpongeStrengths
  use State_mod, only : updatePatches, makeQuiescent
  use Region_mod
  use MPIHelper, only : gracefulExit, writeAndFlush
  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, startTimestep, nTimesteps, reportInterval, saveInterval, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix, message
  logical :: success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_RK4Integrator) :: integrator
  real(wp) :: time

  ! Initialize MPI.
  call MPI_Init(ierror)

  if (command_argument_count() > 1) then
     write(message, '(A)') "Usage: magudi [FILE]"
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)') "High-performance Fortran-based adjoint optimization tool."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)') &
          "FILE is an optional restart file used if running in prediction-only mode."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     call MPI_Finalize(ierror)
     stop -1
  end if

  call startTiming("total")

  ! Parse options from the input file.
  filename = PROJECT_NAME // ".inp"
  call parseInputFile(filename)

  ! Verify that the grid file is in valid PLOT3D format and fetch the grid dimensions:
  ! `globalGridSizes(i,j)` is the number of grid points on grid `j` along dimension `i`.
  call getRequiredOption("grid_file", filename)
  call plot3dDetectFormat(MPI_COMM_WORLD, filename, success,                                 &
       globalGridSizes = globalGridSizes)
  if (.not. success) call gracefulExit(MPI_COMM_WORLD, plot3dErrorMessage)

  ! Setup the region.
  call getRequiredOption("boundary_condition_file", filename)
  call setupRegion(region, MPI_COMM_WORLD, globalGridSizes, filename)

  ! Read the grid file.
  call getRequiredOption("grid_file", filename)
  call loadRegionData(region, QOI_GRID, filename)

  ! If a target state file was specified, read the target state. Otherwise, initialize the
  ! target state to a quiescent state by default.
  if (region%simulationFlags%useTargetState) then
     filename = getOption("target_state_file", "")
     if (len_trim(filename) == 0) then
        do i = 1, size(region%states)
           call makeQuiescent(region%states(i), size(globalGridSizes, 1),                    &
                region%solverOptions%ratioOfSpecificHeats, region%states(i)%targetState)
        end do
     else
        call loadRegionData(region, QOI_TARGET_STATE, filename)
     end if
  end if

  ! Setup spatial discretization.
  do i = 1, size(region%grids)
     call setupSpatialDiscretization(region%grids(i))
  end do
  call MPI_Barrier(region%comm, ierror)

  ! Update the grids by computing the Jacobian, metrics, and norm.
  do i = 1, size(region%grids)
     call updateGrid(region%grids(i))
     call computeSpongeStrengths(region%grids(i), region%patches)
     call updatePatches(region%states(i), region%grids(i),                                   &
          region%patches, region%simulationFlags, region%solverOptions)
  end do
  call MPI_Barrier(region%comm, ierror)

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

  ! Initialize conserved variables.
  if (command_argument_count() == 1) then
     if (region%simulationFlags%predictionOnly) then
        call get_command_argument(1, filename) ! initialize from restart file.
        call loadRegionData(region, QOI_FORWARD_STATE, filename)
     end if
  else if (region%simulationFlags%predictionOnly .or. .not. &
       region%simulationFlags%isBaselineAvailable) then
     if (region%simulationFlags%useTargetState) then
        filename = getOption("initial_condition_file", "")
        if (len_trim(filename) == 0) then
           do i = 1, size(region%states) !... initialize from target state.
              region%states(i)%conservedVariables = region%states(i)%targetState
           end do
        else
           call loadRegionData(region, QOI_FORWARD_STATE, filename) !... initialize from file.
        end if
     else
        call getRequiredOption("initial_condition_file", filename)
        call loadRegionData(region, QOI_FORWARD_STATE, filename) !... initialize from file.
     end if
  end if

  ! Time advancement options.
  nTimesteps = getOption("number_of_timesteps", 1000)
  reportInterval = getOption("report_interval", 1)
  saveInterval = getOption("save_interval", 1000)

  if (region%simulationFlags%predictionOnly) then !... just a predictive simulation.
     time = real(region%states(1)%plot3dAuxiliaryData(4), wp)
     startTimestep = nint(real(region%states(1)%plot3dAuxiliaryData(1), wp))
     call solveForward(region, integrator, time, startTimestep, nTimesteps,                  &
          reportInterval, saveInterval, outputPrefix)
  else

     ! Baseline forward.
     time = 0.0_wp
     if (.not. region%simulationFlags%isBaselineAvailable) then
        call solveForward(region, integrator, time, 0, nTimesteps,                           &
             reportInterval, saveInterval, outputPrefix)
     else
        write(filename, '(2A,I8.8,A)')                                                       &
             trim(outputPrefix), "-", nTimesteps, ".q"
        call loadRegionData(region, QOI_FORWARD_STATE, filename)
        time = real(region%states(1)%plot3dAuxiliaryData(4), wp)
        startTimestep = nint(real(region%states(1)%plot3dAuxiliaryData(1), wp))
     end if

     ! Baseline adjoint.
     call solveAdjoint(region, integrator, time, startTimestep, nTimesteps,                  &
          reportInterval, saveInterval, outputPrefix)

  end if

  call cleanupRK4Integrator(integrator)
  call cleanupRegion(region)

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  ! Finalize MPI.
  call MPI_Finalize(ierror)

contains

  subroutine solveForward(region, integrator, time, startTimestep,                           &
       nTimesteps, reportInterval, saveInterval, outputPrefix)

    ! <<< Derived types >>>
    use State_type
    use Region_type
    use RK4Integrator_type

    ! <<< Internal modules >>>
    use Region_mod, only : saveRegionData
    use RK4Integrator_mod, only : stepForward

    implicit none

    ! <<< Arguments >>>
    type(t_Region) :: region
    type(t_RK4Integrator) :: integrator
    real(SCALAR_KIND), intent(inout) :: time
    integer, intent(in) :: startTimestep, nTimesteps
    integer, intent(in), optional :: reportInterval, saveInterval
    character(len = STRING_LENGTH), intent(in), optional :: outputPrefix

    ! <<< Local variables >>>
    character(len = STRING_LENGTH) :: outputPrefix_, filename
    integer :: timestep
    logical :: verbose

    outputPrefix_ = PROJECT_NAME
    if (present(outputPrefix)) outputPrefix_ = outputPrefix

    verbose = .false.

    timestep = startTimestep
    write(filename, '(2A,I8.8,A)') trim(outputPrefix_), "-", timestep, ".q"
    call saveRegionData(region, QOI_FORWARD_STATE, filename)

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

  subroutine solveAdjoint(region, integrator, time, startTimestep,                           &
       nTimesteps, reportInterval, saveInterval, outputPrefix)

    ! <<< Derived types >>>
    use State_type
    use Region_type
    use RK4Integrator_type

    ! <<< Internal modules >>>
    use Region_mod, only : saveRegionData
    use RK4Integrator_mod, only : stepAdjoint

    implicit none

    ! <<< Arguments >>>
    type(t_Region) :: region
    type(t_RK4Integrator) :: integrator
    real(SCALAR_KIND), intent(inout) :: time
    integer, intent(in) :: startTimestep, nTimesteps
    integer, intent(in), optional :: reportInterval, saveInterval
    character(len = STRING_LENGTH), intent(in), optional :: outputPrefix

    ! <<< Local variables >>>
    character(len = STRING_LENGTH) :: outputPrefix_, filename
    integer :: timestep
    logical :: verbose

    outputPrefix_ = PROJECT_NAME
    if (present(outputPrefix)) outputPrefix_ = outputPrefix

    verbose = .false.

    timestep = startTimestep
    write(filename, '(2A,I8.8,A)') trim(outputPrefix_), "-", timestep, ".adjoint.q"
    call saveRegionData(region, QOI_ADJOINT_STATE, filename)

    do timestep = startTimestep - 1, startTimestep - nTimesteps, -1

       if (present(reportInterval)) verbose = (reportInterval > 0 .and.                      &
            mod(timestep, reportInterval) == 0)
       call stepAdjoint(integrator, region, time, timestep, verbose)

       if (present(saveInterval)) then
          if (saveInterval > 0) then
             if (mod(timestep, saveInterval) == 0) then
                write(filename, '(2A,I8.8,A)')                                               &
                     trim(outputPrefix_), "-", timestep, ".adjoint.q"
                call saveRegionData(region, QOI_ADJOINT_STATE, filename)
             end if
          end if
       end if

    end do

  end subroutine solveAdjoint

end program main
