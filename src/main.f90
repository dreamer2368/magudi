#include "config.h"

program main

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver

  use Grid_enum
  use State_enum

  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use Patch_factory, only : computeSpongeStrengths, updatePatchFactories
  use InterfaceHelper, only : checkInterfaceContinuity
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, timestep, nTimesteps, saveInterval, procRank, numProcs, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix, message
  logical :: success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_Solver) :: solver
  real(wp) :: time
  SCALAR_TYPE :: costFunctional, costSensitivity

  ! Initialize MPI.
  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)

  call initializeErrorHandler()

  if (command_argument_count() > 1) then
     write(message, '(A)') "Usage: magudi [FILE]"
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)') "High-performance Fortran-based adjoint optimization tool."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)')                                                                   &
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

  outputPrefix = getOption("output_prefix", PROJECT_NAME)

  ! Verify that the grid file is in valid PLOT3D format and fetch the grid dimensions:
  ! `globalGridSizes(i,j)` is the number of grid points on grid `j` along dimension `i`.
  call getRequiredOption("grid_file", filename)
  call plot3dDetectFormat(MPI_COMM_WORLD, filename, success,                                 &
       globalGridSizes = globalGridSizes)
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

  ! Save the Jacobian and normalized metrics.
  write(filename, '(2A)') trim(outputPrefix), ".Jacobian.f"
  call region%saveData(QOI_JACOBIAN, filename)
  write(filename, '(2A)') trim(outputPrefix), ".metrics.f"
  call region%saveData(QOI_METRICS, filename)

  ! Initialize the solver.
  if (command_argument_count() == 1) then
     call get_command_argument(1, filename)
     call solver%setup(region, outputPrefix = outputPrefix, restartFilename = filename)
  else
     call solver%setup(region, outputPrefix = outputPrefix)
  end if

  ! Setup boundary conditions.
  call getRequiredOption("boundary_condition_file", filename)
  call region%setupBoundaryConditions(filename)

  ! Compute damping strength on sponge patches.
  do i = 1, size(region%grids)
     call computeSpongeStrengths(region%patchFactories, region%grids(i))
  end do
  call MPI_Barrier(region%comm, ierror)

  ! Check continuity of grid coordinates at block interfaces.
  call checkInterfaceContinuity(region, sqrt(epsilon(0.0_wp)), success)
  if (.not. success) then
     write(message, '(A)') "Grid coordinates are discontinuous across block interfaces!"
     call issueWarning(region%comm, message)
  end if

  ! Update patches.
  do i = 1, size(region%grids)
     call updatePatchFactories(region%patchFactories, region%simulationFlags,                &
          region%solverOptions, region%grids(i), region%states(i))
  end do
  call MPI_Barrier(region%comm, ierror)

  ! Time advancement options.
  time = real(region%states(1)%plot3dAuxiliaryData(4), wp)
  timestep = nint(real(region%states(1)%plot3dAuxiliaryData(1), wp))
  nTimesteps = getOption("number_of_timesteps", 1000)
  saveInterval = getOption("save_interval", 1000)

  if (region%simulationFlags%predictionOnly) then !... just a predictive simulation.
     costFunctional = solver%runForward(region, time, timestep, nTimesteps)
  else

     ! Baseline forward.
     if (.not. region%simulationFlags%isBaselineAvailable) then
        costFunctional = solver%runForward(region, time, timestep, nTimesteps)
     else

        if (region%simulationFlags%steadyStateSimulation) then
           write(filename, '(2A)') trim(outputPrefix), ".steady_state.q"
        else
           timestep = timestep + nTimesteps
           write(filename, '(2A,I8.8,A)')                                                    &
                trim(outputPrefix), "-", timestep, ".q"
        end if

        call region%loadData(QOI_FORWARD_STATE, filename)
        if (.not. region%simulationFlags%steadyStateSimulation) then
           time = region%states(1)%time
           timestep = nint(real(region%states(1)%plot3dAuxiliaryData(1), wp))
        end if

     end if

     ! Baseline adjoint.
     costSensitivity = solver%runAdjoint(region, time, timestep, nTimesteps)

  end if

  call solver%cleanup()
  call region%cleanup()

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  ! Finalize MPI.
  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

end program main
