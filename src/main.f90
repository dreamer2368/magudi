#include "config.h"

program main

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use Region_type
  use TimeIntegrator_mod, only : t_TimeIntegrator

  use Grid_enum
  use State_enum

  use Solver, only : initializeSolver, solveForward, solveAdjoint
  use Region_mod
  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers
  use TimeIntegrator_factory, only : createTimeIntegrator

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, timestep, nTimesteps, reportInterval, saveInterval, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix, message
  logical :: success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()
  real(wp) :: time
  SCALAR_TYPE :: costFunctional

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

  ! Setup spatial discretization.
  do i = 1, size(region%grids)
     call region%grids(i)%setupSpatialDiscretization()
  end do
  call MPI_Barrier(region%comm, ierror)

  ! Update the grids by computing the Jacobian, metrics, and norm.
  do i = 1, size(region%grids)
     call region%grids(i)%update()
     call region%grids(i)%computeSpongeStrengths(region%patches)
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
  call createTimeIntegrator(timeIntegrator, getOption("time_integration_scheme", "RK4"))
  call timeIntegrator%setup(region)

  ! Initialize the solver.
  call initializeSolver(region)

  ! Time advancement options.
  time = real(region%states(1)%plot3dAuxiliaryData(4), wp)
  timestep = nint(real(region%states(1)%plot3dAuxiliaryData(1), wp))
  nTimesteps = getOption("number_of_timesteps", 1000)
  reportInterval = getOption("report_interval", 1)
  saveInterval = getOption("save_interval", 1000)

  ! Update patches.
  do i = 1, size(region%grids)
     call region%states(i)%updatePatches(region%grids(i), region%patches,                    &
          region%simulationFlags, region%solverOptions)
  end do
  call MPI_Barrier(region%comm, ierror)

  if (region%simulationFlags%predictionOnly) then !... just a predictive simulation.
     call solveForward(region, timeIntegrator, time, timestep, nTimesteps,                   &    
          saveInterval, reportInterval, outputPrefix)
  else

     ! Baseline forward.
     if (.not. region%simulationFlags%isBaselineAvailable) then
        call solveForward(region, timeIntegrator, time, timestep, nTimesteps,                &    
             saveInterval, reportInterval, outputPrefix, costFunctional)
     else
        timestep = timestep + nTimesteps
        write(filename, '(2A,I8.8,A)')                                                       &
             trim(outputPrefix), "-", timestep, ".q"
        call loadRegionData(region, QOI_FORWARD_STATE, filename)
        time = real(region%states(1)%plot3dAuxiliaryData(4), wp)
     end if

     ! Baseline adjoint.
     call solveAdjoint(region, timeIntegrator, time, timestep, nTimesteps,                   &    
          saveInterval, reportInterval, outputPrefix)

  end if

  call timeIntegrator%cleanup()
  call cleanupRegion(region)

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  ! Finalize MPI.
  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

end program main
