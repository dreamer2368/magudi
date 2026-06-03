#include "config.h"

program msadjoint

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver
  use Controller_mod, only : t_Controller

  use Grid_enum
  use State_enum
  use Region_enum, only : ADJOINT

  use InputHelper, only : parseInputFile, getFreeUnit, getOption, getRequiredOption
  use InputHelperImpl, only: dict, find
  use ErrorHandler
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers
  use MPIHelper, only : disconnectParentIfSpawned

  implicit none

  type :: t_StateBuffer
     SCALAR_TYPE, allocatable :: F(:,:)
  end type t_StateBuffer

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, k, dictIndex, icDictIndex, nonzeroDictIndex, stat, fileUnit
  integer :: procRank, numProcs, ierror
  integer :: kthArgument, numberOfArguments
  logical :: lookForInput = .false., inputFlag = .false., saveMetricsFlag = .false.
  character(len = STRING_LENGTH) :: argument, inputFilename
  character(len = STRING_LENGTH) :: filename, outputPrefix, message
  character(len = STRING_LENGTH) :: icFilename, endFilename, terminalFilename, icAdjointFilename
  logical :: fileExists, success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_Solver) :: solver
  class(t_Controller), pointer :: controller => null()
  ! Scratch for state arithmetic (matching adjoint pre-computation, IC gradient combination).
  type(t_StateBuffer), allocatable :: scratch(:)

  ! << time-splitting parameters >>
  integer :: Nsplit, Nts, startTimestep
  real(wp) :: penaltyWeight

  SCALAR_TYPE :: dummyValue = 0.0_wp

  ! Initialize MPI.
  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)

  call initializeErrorHandler()

  numberOfArguments = command_argument_count()
  do kthArgument = 1, numberOfArguments
    call get_command_argument(kthArgument,argument)
    select case(trim(adjustl(argument)))
    case("--input")
      lookForInput = .true.
    case("--save_metrics")
      saveMetricsFlag = .true.
    case default
      if (lookForInput) then
        inputFilename = trim(adjustl(argument))
        if (procRank==0) inquire(file=inputFilename,exist=fileExists)
        call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
        if (.not.fileExists) then
           write(message, '(3A)') "input file ", trim(inputFilename), " does not exists!"
           call gracefulExit(MPI_COMM_WORLD, message)
        end if
        lookForInput = .false.
        inputFlag = .true.
      else
        write(message, '(3A)') "option ", trim(argument), " unknown!"
        call gracefulExit(MPI_COMM_WORLD, message)
      end if
    end select
  end do

  call startTiming("total")

  if ( .not. inputFlag ) inputFilename = PROJECT_NAME // ".inp"
  call parseInputFile(inputFilename)

  ! Force the adjoint code path on inside the solver setup, regardless of input file.
  call find("enable_adjoint_solver", dictIndex)
  if (dictIndex < 0) then
    call gracefulExit(MPI_COMM_WORLD,                                                       &
         "magudi.inp must contain 'enable_adjoint_solver' to use msadjoint.")
  end if
  dict(dictIndex)%val = "true"

  ! Per-segment dict-mutation indices (must exist in magudi.inp).
  call find("initial_condition_file", icDictIndex)
  if (icDictIndex < 0) then
    call gracefulExit(MPI_COMM_WORLD,                                                       &
         "magudi.inp must contain 'initial_condition_file' (placeholder ok) for msadjoint.")
  end if
  call find("adjoint_nonzero_initial_condition", nonzeroDictIndex)
  if (nonzeroDictIndex < 0) then
    call gracefulExit(MPI_COMM_WORLD,                                                       &
         "magudi.inp must contain 'adjoint_nonzero_initial_condition' for msadjoint.")
  end if

  outputPrefix = getOption("output_prefix", PROJECT_NAME)

  write(message, '(2A)') "Input file: ", trim(inputFilename)
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

  ! Time-splitting parameters.
  call getRequiredOption("time_splitting/number_of_segments", Nsplit)
  call getRequiredOption("time_splitting/segment_length", Nts)
  call getRequiredOption("time_splitting/start_timestep", startTimestep)
  call getRequiredOption("time_splitting/penalty_weight", penaltyWeight)
  assert(Nsplit >= 1)
  assert(Nts >= 1)
  assert(startTimestep >= 0)

  ! Verify that the grid file is in valid PLOT3D format and fetch the grid dimensions.
  call getRequiredOption("grid_file", filename)
  call plot3dDetectFormat(MPI_COMM_WORLD, filename, success,                                 &
       globalGridSizes = globalGridSizes)
  if (.not. success) call gracefulExit(MPI_COMM_WORLD, plot3dErrorMessage)

  ! Setup the region and load the grid.
  call region%setup(MPI_COMM_WORLD, globalGridSizes)
  call getRequiredOption("grid_file", filename)
  call region%loadData(QOI_GRID, filename)

  ! Update the grids by computing the Jacobian, metrics, and norm.
  do i = 1, size(region%grids)
     call region%grids(i)%update()
  end do
  call MPI_Barrier(region%comm, ierror)

  call region%reportGridDiagnostics()

  if (saveMetricsFlag) then
    write(filename, '(2A)') trim(outputPrefix), ".Jacobian.f"
    call region%saveData(QOI_JACOBIAN, filename)
    write(filename, '(2A)') trim(outputPrefix), ".metrics.f"
    call region%saveData(QOI_METRICS, filename)
  end if

  ! Initialize the solver.
  call solver%setup(region, outputPrefix = outputPrefix)

  ! Save the control and target mollifier if using code-generated values.
  if (region%simulationFlags%enableController) then
     filename = getOption("control_mollifier_file", "")
     if (len_trim(filename) == 0) call region%saveData(QOI_CONTROL_MOLLIFIER,                &
          trim(outputPrefix) // ".control_mollifier.f")
  end if
  if (region%simulationFlags%enableFunctional) then
     filename = getOption("target_mollifier_file", "")
     if (len_trim(filename) == 0) call region%saveData(QOI_TARGET_MOLLIFIER,                 &
          trim(outputPrefix) // ".target_mollifier.f")
  end if

  ! Pre-pass: clear the global gradient .dat files. The segment loop below passes
  ! deleteGradientFile=.false. so accumulation works across segments; this one-shot call
  ! starts each msadjoint invocation from an empty gradient file (default deleteGradientFile=.true.).
  if (region%simulationFlags%enableController) then
    call solver%controllerFactory%connect(controller)
    assert(associated(controller))
    call controller%hookBeforeTimemarch(region, ADJOINT)
  end if

  ! Allocate scratch buffer used for matching-adjoint construction and IC-gradient assembly.
  allocate(scratch(size(region%grids)))
  do i = 1, size(region%grids)
    allocate(scratch(i)%F(region%grids(i)%nGridPoints, region%solverOptions%nUnknowns))
  end do

  ! Pre-compute matching adjoint terminals: for k = 0..Nsplit-2, save
  !   matching_k = penaltyWeight * (end_k - ic_{k+1})
  ! as an adjoint q-file at <prefix>-<startTs + (k+1)*Nts:08d>.adjoint.q. runAdjoint
  ! will load this as the adjoint terminal when adjoint_nonzero_initial_condition = "true".
  do k = 0, Nsplit-2
    ! Load end_k.
    write(endFilename, '(2A,I8.8,A)') trim(outputPrefix), "-",                              &
         startTimestep + (k+1) * Nts, ".q"
    if (procRank == 0) inquire(file = endFilename, exist = fileExists)
    call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
    if (.not. fileExists) then
      write(message, '(3A)') "End-of-segment snapshot ", trim(endFilename),                 &
           " missing (run msforward first)."
      call gracefulExit(MPI_COMM_WORLD, message)
    end if
    call region%loadData(QOI_FORWARD_STATE, endFilename)
    do i = 1, size(region%states)
      scratch(i)%F = region%states(i)%conservedVariables
    end do

    ! Load ic_{k+1}, form matching_k = w * (end_k - ic_{k+1}) into adjointVariables.
    write(icFilename, '(2A,I0,A)') trim(outputPrefix), "-", k+1, ".ic.q"
    if (procRank == 0) inquire(file = icFilename, exist = fileExists)
    call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
    if (.not. fileExists) then
      write(message, '(3A)') "Per-segment IC file ", trim(icFilename), " does not exist!"
      call gracefulExit(MPI_COMM_WORLD, message)
    end if
    call region%loadData(QOI_FORWARD_STATE, icFilename)
    do i = 1, size(region%states)
      region%states(i)%adjointVariables = penaltyWeight *                                   &
           (scratch(i)%F - region%states(i)%conservedVariables)
    end do

    ! Save as adjoint state at the segment-boundary timestep. region%timestep is now
    ! startTs + (k+1)*Nts (from the ic_{k+1} load? actually ic_{k+1}'s embedded timestep
    ! should match startTs + (k+1)*Nts since it's the IC at start of segment k+1).
    ! Build the filename explicitly to match runAdjoint's loadInitialCondition convention.
    write(terminalFilename, '(2A,I8.8,A)') trim(outputPrefix), "-",                         &
         startTimestep + (k+1) * Nts, ".adjoint.q"
    call region%saveData(QOI_ADJOINT_STATE, terminalFilename)
  end do

  ! Adjoint segment loop: process segments in reverse order.
  do k = Nsplit-1, 0, -1
    write(message, '(A,I0,A,I0,A)') "=== msadjoint: segment ", k, " of ", Nsplit, " ==="
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

    write(icFilename, '(2A,I0,A)') trim(outputPrefix), "-", k, ".ic.q"
    if (procRank == 0) inquire(file = icFilename, exist = fileExists)
    call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
    if (.not. fileExists) then
      write(message, '(3A)') "Per-segment IC file ", trim(icFilename), " does not exist!"
      call gracefulExit(MPI_COMM_WORLD, message)
    end if

    dict(icDictIndex)%val = trim(icFilename)
    if (k == Nsplit-1) then
      ! No segment to the right -> zero adjoint terminal (runAdjoint synthesizes it).
      dict(nonzeroDictIndex)%val = "false"
    else
      ! Load the matching adjoint terminal pre-computed above.
      dict(nonzeroDictIndex)%val = "true"
    end if

    dummyValue = solver%runAdjoint(region, controlTimestepOffset = k * Nts,                     &
                                   deleteGradientFile = .false.)

    ! Save the IC-side adjoint (= dJ_time_integral_via_segment_k / d(ic_k); the matching
    ! penalty's direct contribution to ic_k is added in the post-pass below).
    write(icAdjointFilename, '(2A,I0,A)') trim(outputPrefix), "-", k, ".ic.adjoint.q"
    call region%saveData(QOI_ADJOINT_STATE, icAdjointFilename)

    call MPI_Barrier(MPI_COMM_WORLD, ierror)
  end do

  ! Post-pass: add the matching penalty's direct contribution to each IC gradient for k >= 1.
  !   dJ_p/d(ic_k) = penaltyWeight * (ic_k - end_{k-1})
  do k = 1, Nsplit-1
    write(icAdjointFilename, '(2A,I0,A)') trim(outputPrefix), "-", k, ".ic.adjoint.q"
    call region%loadData(QOI_ADJOINT_STATE, icAdjointFilename)
    ! adjointVariables now holds the time-integral adjoint at start of segment k.

    write(icFilename, '(2A,I0,A)') trim(outputPrefix), "-", k, ".ic.q"
    call region%loadData(QOI_FORWARD_STATE, icFilename)
    do i = 1, size(region%states)
      scratch(i)%F = region%states(i)%conservedVariables
    end do

    write(endFilename, '(2A,I8.8,A)') trim(outputPrefix), "-",                              &
         startTimestep + k * Nts, ".q"
    call region%loadData(QOI_FORWARD_STATE, endFilename)
    do i = 1, size(region%states)
      region%states(i)%adjointVariables = region%states(i)%adjointVariables +               &
           penaltyWeight * (scratch(i)%F - region%states(i)%conservedVariables)
    end do

    call region%saveData(QOI_ADJOINT_STATE, icAdjointFilename)
  end do

  do i = 1, size(scratch)
    SAFE_DEALLOCATE(scratch(i)%F)
  end do
  SAFE_DEALLOCATE(scratch)

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  call solver%cleanup()
  call region%cleanup()

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  call cleanupErrorHandler()
  call disconnectParentIfSpawned()
  call MPI_Finalize(ierror)

end program msadjoint
