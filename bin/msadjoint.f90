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
  integer :: i, k, m, dictIndex, icDictIndex, nonzeroDictIndex, stat, fileUnit
  integer :: procRank, numProcs, ierror
  integer :: kthArgument, numberOfArguments
  logical :: lookForInput = .false., lookForOutput = .false., lookForSubOutput = .false.,   &
              inputFlag = .false., outputFlag = .false., subOutputFlag = .false.,           &
              saveMetricsFlag = .false.
  character(len = STRING_LENGTH) :: argument, inputFilename, outputFilename, subOutputFilename
  character(len = STRING_LENGTH) :: filename, outputPrefix, segPrefix, message
  character(len = STRING_LENGTH) :: icFilename, endFilename, terminalFilename, icAdjointFilename
  character(len = STRING_LENGTH) :: stateMollifierFilename
  logical :: fileExists, success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_Solver) :: solver
  class(t_Controller), pointer :: controller => null()
  ! Scratch for state arithmetic (matching adjoint pre-computation, IC gradient combination).
  type(t_StateBuffer), allocatable :: scratch(:)
  ! Buffer that carries the in-memory adjoint_start_of_{k+1} from the previous
  ! iteration's runAdjoint into the current iteration's on-the-fly IC-gradient
  ! assembly and (with state mollifier) matching-terminal blend.
  type(t_StateBuffer), allocatable :: icAdjBuf(:)

  ! << time-splitting parameters >>
  integer :: Nsplit, Nts, startTimestep
  real(wp) :: penaltyWeight

  ! << state-mollifier flags >>
  logical :: useStateMollifier, stateMollifierUniformInTime

  ! << per-segment gradient inner-product accumulators >>
  ! segmentCtrlIP(k) = <g_ctrl, g_ctrl>_M restricted to segment k's time window (from runAdjoint).
  ! segmentICIP(k)   = <g_ic_k, g_ic_k>_SBP; left at 0 for k=0 since ic_0 is not optimized.
  SCALAR_TYPE, allocatable :: segmentCtrlIP(:), segmentICIP(:)
  SCALAR_TYPE :: L2sq

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
    case("--output")
      lookForOutput = .true.
    case("--segments-output")
      lookForSubOutput = .true.
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
      elseif (lookForOutput) then
        outputFilename = trim(adjustl(argument))
        lookForOutput = .false.
        outputFlag = .true.
      elseif (lookForSubOutput) then
        subOutputFilename = trim(adjustl(argument))
        lookForSubOutput = .false.
        subOutputFlag = .true.
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

  if ( .not. outputFlag )    outputFilename    = trim(outputPrefix) // ".adjoint_run.txt"
  if ( .not. subOutputFlag ) subOutputFilename = trim(outputPrefix) // ".sub_adjoint_run.txt"

  write(message, '(2A)') "Input file: ", trim(inputFilename)
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  write(message, '(2A)') "Output file (total inner product): ", trim(outputFilename)
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  write(message, '(2A)') "Per-segment output file: ", trim(subOutputFilename)
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

  ! Time-splitting parameters.
  call getRequiredOption("time_splitting/number_of_segments", Nsplit)
  call getRequiredOption("time_splitting/segment_length", Nts)
  call getRequiredOption("time_splitting/start_timestep", startTimestep)
  call getRequiredOption("time_splitting/penalty_weight", penaltyWeight)
  assert(Nsplit >= 1)
  assert(Nts >= 1)
  assert(startTimestep >= 0)

  ! State-mollifier flags (must agree with msforward's reading of the same options).
  useStateMollifier = getOption(                                                            &
       "time_splitting/state_mollifier/enabled", .false.)
  stateMollifierUniformInTime = getOption(                                                  &
       "time_splitting/state_mollifier/uniform_in_time", .true.)

  ! Per-segment gradient inner-product accumulators (segmentICIP(0) stays 0; ic_0 is fixed).
  allocate(segmentCtrlIP(0:Nsplit-1)); segmentCtrlIP = 0.0_wp
  allocate(segmentICIP  (0:Nsplit-1)); segmentICIP   = 0.0_wp

  ! Verify that the grid file is in valid PLOT3D format and fetch the grid dimensions.
  call getRequiredOption("grid_file", filename)
  call plot3dDetectFormat(MPI_COMM_WORLD, filename, success,                                 &
       globalGridSizes = globalGridSizes)
  if (.not. success) call gracefulExit(MPI_COMM_WORLD, plot3dErrorMessage)

  ! Setup the region and load the grid.
  call region%setup(MPI_COMM_WORLD, globalGridSizes)
  call getRequiredOption("grid_file", filename)
  call region%loadData(QOI_GRID, filename)

  ! Allocate the state-mollifier field per grid and, for uniform-in-time mode,
  ! load it once. (Per-segment mode loads inside the segment loop below.)
  if (useStateMollifier) then
    do i = 1, size(region%grids)
      allocate(region%grids(i)%stateMollifier(region%grids(i)%nGridPoints, 1))
    end do
    if (stateMollifierUniformInTime) then
      write(stateMollifierFilename, '(2A)') trim(outputPrefix), ".ic_mollifier.f"
      call region%loadData(QOI_STATE_MOLLIFIER, stateMollifierFilename)
    end if
  end if

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

  ! Buffer that carries adjoint_start_of_{k+1} from iteration k+1 into iteration k
  ! (used by both step A1 IC-grad assembly and step A2 matching-terminal blend).
  allocate(icAdjBuf(size(region%grids)))
  do i = 1, size(region%grids)
    allocate(icAdjBuf(i)%F(region%grids(i)%nGridPoints, region%solverOptions%nUnknowns))
  end do

  ! Unified adjoint segment loop. Two roles per iteration k (going backward):
  !
  !   Step A (k < Nsplit-1): reuse the in-memory adjoint_start_of_{k+1} produced
  !     by the previous iteration's runAdjoint to (A1) assemble segment k+1's
  !     final IC gradient and (A2) build segment k's matching adjoint terminal,
  !     each in one pass with the file I/O folded together. This collapses the
  !     old segment-loop + post-pass into one loop and skips the write-then-
  !     reread of the bare adjoint that ver3 / earlier ver4 paid for.
  !   Step B: runAdjoint over segment k. Leaves adjoint_start_of_k in
  !     region%states(:)%adjointVariables for the next iteration to consume.
  !
  ! Order matches ver3's adjointRunCommand (base_extension.py:113-145): the
  ! matching-terminal patchup uses the BARE adjoint_start_of_{k+1} (not the
  ! penalty-augmented / mollifier-weighted gradient); the penalty add + w
  ! multiplication only enter the final IC gradient.
  do k = Nsplit-1, 0, -1
    write(message, '(A,I0,A,I0,A)') "=== msadjoint: segment ", k, " of ", Nsplit, " ==="
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

    !==========================================================================
    ! STEP A. Build segment k+1's final IC gradient and segment k's matching
    ! adjoint terminal on the fly. Skip on the first iteration (k=Nsplit-1):
    ! no segment to the right, runAdjoint synthesizes the zero terminal.
    !==========================================================================
    if (k < Nsplit - 1) then

      ! Preserve adjoint_start_of_{k+1} from the previous iteration's runAdjoint.
      ! Needed by both A1 (penalty add) and A2 (mollifier blend).
      do i = 1, size(region%states)
        icAdjBuf(i)%F = region%states(i)%adjointVariables
      end do

      ! Load end_k -> conservedVariables, copy to scratch. (end_k lives under
      ! segment k's prefix, written by msforward's runForward via showProgress.)
      write(endFilename, '(2A,I0,A,I8.8,A)') trim(outputPrefix), "-", k, "-",               &
           startTimestep + (k+1) * Nts, ".q"
      if (procRank == 0) inquire(file = endFilename, exist = fileExists)
      call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
      if (.not. fileExists) then
        write(message, '(3A)') "End-of-segment snapshot ", trim(endFilename),               &
             " missing (run msforward first)."
        call gracefulExit(MPI_COMM_WORLD, message)
      end if
      call region%loadData(QOI_FORWARD_STATE, endFilename)
      do i = 1, size(region%states)
        scratch(i)%F = region%states(i)%conservedVariables    ! = end_k
      end do

      ! Load ic_{k+1} -> conservedVariables (the un-blended optimization variable).
      write(icFilename, '(2A,I0,A)') trim(outputPrefix), "-", k+1, ".ic.q"
      if (procRank == 0) inquire(file = icFilename, exist = fileExists)
      call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
      if (.not. fileExists) then
        write(message, '(3A)') "Per-segment IC file ", trim(icFilename), " does not exist!"
        call gracefulExit(MPI_COMM_WORLD, message)
      end if
      call region%loadData(QOI_FORWARD_STATE, icFilename)

      ! Per-segment mollifier for index (k+1): drives both A1 and A2.
      if (useStateMollifier .and. .not. stateMollifierUniformInTime) then
        write(stateMollifierFilename, '(2A,I0,A)')                                          &
             trim(outputPrefix), "-", k+1, ".ic_mollifier.f"
        call region%loadData(QOI_STATE_MOLLIFIER, stateMollifierFilename)
      end if

      ! --- A1. Final IC gradient of segment k+1 ---------------------------------
      ! dJ/d(ic_{k+1}) = w · (adjoint_start_of_{k+1} + pw · (ic_{k+1} - end_k))
      ! (w = stateMollifier when enabled, else identity.)
      do i = 1, size(region%states)
        region%states(i)%adjointVariables = icAdjBuf(i)%F +                                 &
             penaltyWeight * (region%states(i)%conservedVariables - scratch(i)%F)
        if (useStateMollifier) then
          do m = 1, region%solverOptions%nUnknowns
            region%states(i)%adjointVariables(:,m) =                                        &
                 region%grids(i)%stateMollifier(:,1) *                                       &
                 region%states(i)%adjointVariables(:,m)
          end do
        end if
      end do

      ! segmentICIP(k+1) = <grad, grad>_SBP (plain norm of the finalized gradient,
      ! matching ver3's gradient inner-product convention).
      L2sq = 0.0_wp
      do i = 1, size(region%states)
        L2sq = L2sq + region%grids(i)%computeInnerProduct(                                  &
             region%states(i)%adjointVariables, region%states(i)%adjointVariables)
      end do
      if (region%commGridMasters /= MPI_COMM_NULL)                                          &
           call MPI_Allreduce(MPI_IN_PLACE, L2sq, 1, SCALAR_TYPE_MPI, MPI_SUM,              &
                              region%commGridMasters, ierror)
      do i = 1, size(region%grids)
        call MPI_Bcast(L2sq, 1, SCALAR_TYPE_MPI, 0, region%grids(i)%comm, ierror)
      end do
      segmentICIP(k+1) = L2sq

      write(icAdjointFilename, '(2A,I0,A)') trim(outputPrefix), "-", k+1, ".ic.adjoint.q"
      call region%saveData(QOI_ADJOINT_STATE, icAdjointFilename)

      ! --- A2. Segment k's matching adjoint terminal ----------------------------
      ! Without mollifier: terminal_k = pw · (end_k - ic_{k+1})
      ! With mollifier:    terminal_k = w · pw · (end_k - ic_{k+1})
      !                               + (1-w) · adjoint_start_of_{k+1}
      ! The terminal file lives under segment k's prefix so runAdjoint
      ! (with solver%outputPrefix = <prefix>-<k> below) picks it up via its
      ! "<this%outputPrefix>-<region%timestep>.adjoint.q" path.
      do i = 1, size(region%states)
        region%states(i)%adjointVariables =                                                 &
             penaltyWeight * (scratch(i)%F - region%states(i)%conservedVariables)
        if (useStateMollifier) then
          do m = 1, region%solverOptions%nUnknowns
            region%states(i)%adjointVariables(:,m) =                                        &
                 region%grids(i)%stateMollifier(:,1) *                                       &
                 region%states(i)%adjointVariables(:,m) +                                   &
                 (1.0_wp - region%grids(i)%stateMollifier(:,1)) * icAdjBuf(i)%F(:,m)
          end do
        end if
      end do

      write(terminalFilename, '(2A,I0,A,I8.8,A)') trim(outputPrefix), "-", k, "-",          &
           startTimestep + (k+1) * Nts, ".adjoint.q"
      call region%saveData(QOI_ADJOINT_STATE, terminalFilename)

      dict(nonzeroDictIndex)%val = "true"
    else
      ! No segment to the right -> zero adjoint terminal (runAdjoint synthesizes it).
      dict(nonzeroDictIndex)%val = "false"
    end if

    !==========================================================================
    ! STEP B. runAdjoint for segment k.
    !==========================================================================
    write(icFilename, '(2A,I0,A)') trim(outputPrefix), "-", k, ".ic.q"
    if (procRank == 0) inquire(file = icFilename, exist = fileExists)
    call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
    if (.not. fileExists) then
      write(message, '(3A)') "Per-segment IC file ", trim(icFilename), " does not exist!"
      call gracefulExit(MPI_COMM_WORLD, message)
    end if

    dict(icDictIndex)%val = trim(icFilename)

    ! Mutate solver%outputPrefix to <prefix>-<k> so runAdjoint, showProgress, and
    ! the reverse-time migrator all key their I/O off segment k's namespace.
    ! That makes the migrator's "<prefix>-<startTimestep>.q" read pick up the
    ! IC snapshot written by msforward (which is the blended IC when state-
    ! mollifier blending is on, or the raw ic_k otherwise).
    !
    ! Mutating solver%outputPrefix here does NOT change the actuator's
    ! .control_forcing_<name>.dat / .gradient_<name>.dat paths -- those were
    ! resolved once during setupActuatorPatch (src/ActuatorPatchImpl.f90:52,68)
    ! via getOption("output_prefix", PROJECT_NAME) and cached on the patch as
    ! controlForcingFilename / gradientFilename. They stay under the global
    ! prefix, so segments keep reading/writing one shared .dat. The driver
    ! supplies (controlTimestepOffset, controlTotalTimesteps) = (k*Nts, Nsplit*Nts);
    ! runAdjoint derives the FORWARD referenceTimestep (= controlTimestepOffset,
    ! counted from start) and the ADJOINT referenceTimestep
    ! (= controlTotalTimesteps - controlTimestepOffset - Nts, counted from end
    ! in reverse) per the time-reversed file convention.
    write(segPrefix, '(2A,I0)') trim(outputPrefix), "-", k
    solver%outputPrefix = trim(segPrefix)

    segmentCtrlIP(k) = solver%runAdjoint(region,                                            &
                                         controlTimestepOffset = k * Nts,                   &
                                         controlTotalTimesteps = Nsplit * Nts,              &
                                         deleteGradientFile = .false.)

    call MPI_Barrier(MPI_COMM_WORLD, ierror)
  end do

  ! Post-loop: After k=0's runAdjoint, adjointVariables = adjoint_start_of_0
  ! = dJ/d(ic_blended_0). ic_0 is not optimized (no penalty contribution and
  ! segmentICIP(0) stays 0), but we still write the file because parallel_io's
  ! _build_paths("grad") lists <prefix>-0.ic.adjoint.q whenever the schema
  ! includes ic slot 0.
  write(icAdjointFilename, '(2A,I0,A)') trim(outputPrefix), "-", 0, ".ic.adjoint.q"
  call region%saveData(QOI_ADJOINT_STATE, icAdjointFilename)

  do i = 1, size(scratch)
    SAFE_DEALLOCATE(scratch(i)%F)
  end do
  SAFE_DEALLOCATE(scratch)

  do i = 1, size(icAdjBuf)
    SAFE_DEALLOCATE(icAdjBuf(i)%F)
  end do
  SAFE_DEALLOCATE(icAdjBuf)

  ! Aggregate gradient inner product across segments:
  !   total = sum_k (control_forcing_IP_k + ic_IP_k); ic_IP_0 = 0 since ic_0 is not optimized.
  dummyValue = 0.0_wp
  do k = 0, Nsplit-1
    dummyValue = dummyValue + segmentCtrlIP(k) + segmentICIP(k)
  end do

  write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))')                                     &
       'msadjoint: total gradient inner product = ', dummyValue
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  if (procRank == 0) then
    open(unit = getFreeUnit(fileUnit), file = trim(outputFilename), action='write',         &
         iostat = stat, status = 'replace')
    write(fileUnit, '(1X,SP,' // SCALAR_FORMAT // ')') dummyValue
    close(fileUnit)

    open(unit = getFreeUnit(fileUnit), file = trim(subOutputFilename), action='write',      &
         iostat = stat, status = 'replace')
    write(fileUnit, '(A)') "# segment  control_forcing_inner_product   ic_inner_product"
    do k = 0, Nsplit-1
      write(fileUnit, '(I8,2(1X,SP,' // SCALAR_FORMAT // '))')                              &
           k, segmentCtrlIP(k), segmentICIP(k)
    end do
    close(fileUnit)
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  SAFE_DEALLOCATE(segmentCtrlIP)
  SAFE_DEALLOCATE(segmentICIP)

  call solver%cleanup()
  call region%cleanup()

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  call cleanupErrorHandler()
  call disconnectParentIfSpawned()
  call MPI_Finalize(ierror)

end program msadjoint
