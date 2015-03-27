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
  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, timestep, nTimesteps, reportInterval, saveInterval, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix, message
  logical :: success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_RK4Integrator) :: integrator
  real(wp) :: time

  ! Initialize MPI.
  call MPI_Init(ierror)

  call initializeErrorHandler()

  if (command_argument_count() > 1) then
     write(message, '(A)') "Usage: magudi [FILE]"
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)') "High-performance Fortran-based adjoint optimization tool."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)') &
          "FILE is an optional restart file used if running in prediction-only mode."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     call cleanupErrorHandler()
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
     timestep = nint(real(region%states(1)%plot3dAuxiliaryData(1), wp))
     call solveForward(region, integrator, time, timestep, nTimesteps,                       &
          saveInterval, reportInterval, outputPrefix)
  else

     ! Control mollifier.
     filename = getOption("control_mollifier_file", "")
     if (len_trim(filename) == 0) then
        do i = 1, size(region%grids)
           region%grids(i)%controlMollifier = 1.0_wp
        end do
     else
        call loadRegionData(region, QOI_CONTROL_MOLLIFIER, filename)
     end if

     ! Target mollifier.
     filename = getOption("target_mollifier_file", "")
     if (len_trim(filename) == 0) then
        do i = 1, size(region%grids)
           region%grids(i)%targetMollifier = 1.0_wp
        end do
     else
        call loadRegionData(region, QOI_TARGET_MOLLIFIER, filename)
     end if

     ! Baseline forward.
     if (.not. region%simulationFlags%isBaselineAvailable) then
        time = 0.0_wp
        timestep = 0
        call solveForward(region, integrator, time, timestep, nTimesteps,                    &
             saveInterval, reportInterval, outputPrefix)
     else
        write(filename, '(2A,I8.8,A)')                                                       &
             trim(outputPrefix), "-", nTimesteps, ".q"
        call loadRegionData(region, QOI_FORWARD_STATE, filename)
        time = real(region%states(1)%plot3dAuxiliaryData(4), wp)
        timestep = nint(real(region%states(1)%plot3dAuxiliaryData(1), wp))
     end if

     ! Baseline adjoint.
     call solveAdjoint(region, integrator, time, timestep, nTimesteps,                       &
          saveInterval, reportInterval, outputPrefix)

  end if

  call cleanupRK4Integrator(integrator)
  call cleanupRegion(region)

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  ! Finalize MPI.
  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

contains

  subroutine solveForward(region, integrator, time, timestep, nTimesteps,                    &
       saveInterval, reportInterval, outputPrefix)

    ! <<< Derived types >>>
    use State_type, only : t_State
    use Region_type, only : t_Region
    use RK4Integrator_type, only : t_RK4Integrator

    ! <<< Internal modules >>>
    use Region_mod, only : saveRegionData, reportResiduals
    use ErrorHandler, only : writeAndFlush
    use RK4Integrator_mod, only : substepForward

    implicit none

    ! <<< Arguments >>>
    type(t_Region) :: region
    type(t_RK4Integrator) :: integrator
    real(SCALAR_KIND), intent(inout) :: time
    integer, intent(inout) :: timestep
    integer, intent(in) :: nTimesteps
    integer, intent(in), optional :: saveInterval, reportInterval
    character(len = STRING_LENGTH), intent(in), optional :: outputPrefix

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    character(len = STRING_LENGTH) :: outputPrefix_, filename, str
    integer :: i, timestep_
    logical :: verbose

    assert(timestep >= 0)

    outputPrefix_ = PROJECT_NAME
    if (present(outputPrefix)) outputPrefix_ = outputPrefix

    write(filename, '(2A,I8.8,A)') trim(outputPrefix_), "-", timestep, ".q"
    call saveRegionData(region, QOI_FORWARD_STATE, filename)

    verbose = .false.

    do timestep_ = timestep + 1, timestep + nTimesteps

       if (present(reportInterval)) verbose = (reportInterval > 0 .and.                      &
            mod(timestep_, reportInterval) == 0)

       do i = 1, 4
          call substepForward(integrator, region, time, timestep_, i)
       end do

       if (verbose) then
          if (region%simulationFlags%useConstantCfl) then
             write(str, '(2A,I8,2(A,D13.6))') PROJECT_NAME, ": timestep = ", timestep_,      &
                  ", dt = ", region%states(1)%timeStepSize, ", time = ", time
          else
             write(str, '(2A,I8,2(A,D13.6))') PROJECT_NAME, ": timestep = ", timestep_,      &
                  ", CFL = ", region%states(1)%cfl, ", time = ", time
          end if
          call writeAndFlush(region%comm, output_unit, str)
          if (region%simulationFlags%steadyStateSimulation) call reportResiduals(region)
       end if

       if (present(saveInterval)) then
          if (saveInterval > 0) then
             if (mod(timestep_, saveInterval) == 0) then
                do i = 1, size(region%states)
                   region%states(i)%plot3dAuxiliaryData(1) = real(timestep_, wp)
                   region%states(i)%plot3dAuxiliaryData(4) = time
                end do
                write(filename, '(2A,I8.8,A)') trim(outputPrefix_), "-", timestep_, ".q"
                call saveRegionData(region, QOI_FORWARD_STATE, filename)
             end if
          end if
       end if

    end do

    timestep = timestep + nTimesteps

  end subroutine solveForward

  subroutine solveAdjoint(region, integrator, time, timestep, nTimesteps,                    &
       saveInterval, reportInterval, outputPrefix)

    ! <<< Derived types >>>
    use State_type, only : t_State
    use Region_type, only : t_Region
    use RK4Integrator_type, only : t_RK4Integrator
    use ReverseMigrator_type, only : t_ReverseMigrator

    ! <<< Internal modules >>>
    use Region_mod, only : saveRegionData, reportResiduals
    use InputHelper, only : getOption
    use ErrorHandler, only : writeAndFlush
    use RK4Integrator_mod, only : substepAdjoint
    use ReverseMigrator_mod

    implicit none

    ! <<< Arguments >>>
    type(t_Region) :: region
    type(t_RK4Integrator) :: integrator
    real(SCALAR_KIND), intent(inout) :: time
    integer, intent(inout) :: timestep
    integer, intent(in) :: nTimesteps, saveInterval
    integer, intent(in), optional :: reportInterval
    character(len = STRING_LENGTH), intent(in), optional :: outputPrefix

    ! <<< Local variables >>>
    character(len = STRING_LENGTH) :: outputPrefix_, filename, str
    type(t_ReverseMigrator) :: reverseMigrator
    integer :: i, timestep_
    logical :: verbose

    assert(timestep >= nTimesteps)

    outputPrefix_ = PROJECT_NAME
    if (present(outputPrefix)) outputPrefix_ = outputPrefix

    call setupReverseMigrator(reverseMigrator, region, outputPrefix_,                        &
         getOption("checkpointing_scheme", "uniform checkpointing"),                         &
         timestep - nTimesteps, timestep,                                                    &
         saveInterval, saveInterval * 4)

    write(filename, '(2A,I8.8,A)') trim(outputPrefix_), "-", timestep, ".adjoint.q"
    call saveRegionData(region, QOI_ADJOINT_STATE, filename)

    verbose = .false.

    do timestep_ = timestep - 1, timestep - nTimesteps, -1

       if (present(reportInterval)) verbose = (reportInterval > 0 .and.                      &
            mod(timestep_, reportInterval) == 0)

       do i = 4, 1, -1
          if (.not. region%simulationFlags%steadyStateSimulation) then
             if (i == 1) then
                call migrateToSubstep(reverseMigrator, region,                               &
                     integrator, timestep_, 4)
             else
                call migrateToSubstep(reverseMigrator, region,                               &
                     integrator, timestep_ + 1, i - 1)
             end if
          end if
          call substepAdjoint(integrator, region, time, timestep_, i)
       end do

       if (verbose) then
          if (region%simulationFlags%useConstantCfl) then
             write(str, '(2A,I8,2(A,D13.6))') PROJECT_NAME, ": timestep = ", timestep_,      &
                  ", dt = ", region%states(1)%timeStepSize, ", time = ", time
          else
             write(str, '(2A,I8,2(A,D13.6))') PROJECT_NAME, ": timestep = ", timestep_,      &
                  ", CFL = ", region%states(1)%cfl, ", time = ", time
          end if
          call writeAndFlush(region%comm, output_unit, str)
          if (region%simulationFlags%steadyStateSimulation) call reportResiduals(region)
       end if

       if (saveInterval > 0) then
          if (mod(timestep_, saveInterval) == 0) then
             do i = 1, size(region%states)
                region%states(i)%plot3dAuxiliaryData(1) = real(timestep_, wp)
                region%states(i)%plot3dAuxiliaryData(4) = time
             end do
             write(filename, '(2A,I8.8,A)')                                                  &
                  trim(outputPrefix_), "-", timestep_, ".adjoint.q"
             call saveRegionData(region, QOI_ADJOINT_STATE, filename)
          end if
       end if

    end do

    call cleanupReverseMigrator(reverseMigrator)

    timestep = timestep - nTimesteps

  end subroutine solveAdjoint

end program main
