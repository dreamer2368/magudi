#include "config.h"

program msforward

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver
  use MSPenalty_mod, only : t_MSPenalty, connectMSPenalty

  use Grid_enum
  use State_enum

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
  integer :: i, k, m, icDictIndex, stat, fileUnit, procRank, numProcs, ierror
  integer :: kthArgument, numberOfArguments
  logical :: lookForInput = .false., lookForOutput = .false., lookForSubOutput = .false.,   &
              inputFlag = .false., outputFlag = .false., subOutputFlag = .false.,           &
              saveMetricsFlag = .false.
  character(len = STRING_LENGTH) :: argument, inputFilename, outputFilename, subOutputFilename
  character(len = STRING_LENGTH) :: filename, outputPrefix, segPrefix, message
  character(len = STRING_LENGTH) :: icFilename
  character(len = STRING_LENGTH) :: stateMollifierFilename, blendedIcFilename
  logical :: fileExists, success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_Solver) :: solver
  ! After segment k completes, scratch(:)%F holds end_k for use in segment k+1's penalty.
  type(t_StateBuffer), allocatable :: scratch(:)
  ! Auxiliary buffer for state-mollifier IC blending: holds ic_k while we
  ! overwrite conservedVariables with the blended IC for the solver.
  type(t_StateBuffer), allocatable :: blendBuf(:)

  ! << time-splitting parameters >>
  integer :: Nsplit, Nts, startTimestep
  real(wp) :: penaltyWeight

  ! << penalty norm (quadratic by default, Huber optional) >>
  character(len = STRING_LENGTH) :: penaltyType
  real(wp) :: huberThreshold
  class(t_MSPenalty), pointer :: penalty => null()

  ! << state-mollifier flags >>
  logical :: useStateMollifier, stateMollifierUniformInTime

  ! << per-segment accumulators >>
  SCALAR_TYPE, allocatable :: segmentCost(:), segmentL2sq(:)
  SCALAR_TYPE :: J, JtimeIntegral, Jpenalty, L2sq, L2sqSum

  ! Set when any segment's runForward trips the hard solution-limit check
  ! and returns HUGE(0.0_wp). Causes the segment loop to bail and the
  ! output writer to emit HUGE to outputFilename for the Python driver.
  logical :: solutionCrashed

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

  ! Per-segment dict mutation: ic file is rewritten each iteration so runForward's
  ! built-in loadInitialCondition path picks up ic_k. The segment-prefixed
  ! solver%outputPrefix makes runForward save ic_k as <prefix>-<k>-<k*Nts>.q,
  ! which the reverse-time migrator inside msadjoint later loads as the
  ! segment start state.
  call find("initial_condition_file", icDictIndex)
  if (icDictIndex < 0) then
    call gracefulExit(MPI_COMM_WORLD,                                                       &
         "magudi.inp must contain 'initial_condition_file' (placeholder ok) for msforward.")
  end if

  outputPrefix = getOption("output_prefix", PROJECT_NAME)

  if ( .not. outputFlag )    outputFilename    = trim(outputPrefix) // ".forward_run.txt"
  if ( .not. subOutputFlag ) subOutputFilename = trim(outputPrefix) // ".sub_J.txt"

  write(message, '(2A)') "Input file: ", trim(inputFilename)
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  write(message, '(2A)') "Output file (total J): ", trim(outputFilename)
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

  ! Penalty norm: quadratic recovers the legacy 0.5*pw*sum_k(L2sq_k); huber applies
  ! a smooth-absolute-value norm to the summed mismatch, with a kink at threshold.
  penaltyType    = getOption("time_splitting/penalty_type", "quadratic")
  huberThreshold = getOption("time_splitting/huber/threshold", real(1.0e-12, wp))
  call connectMSPenalty(penalty, trim(penaltyType), huberThreshold)

  ! State-mollifier flags. When enabled, the matching penalty is
  ! mollifier-weighted and segment k's IC is blended with end_{k-1} so the
  ! solver sees a continuous field in the passive region.
  useStateMollifier = getOption(                                                            &
       "time_splitting/state_mollifier/enabled", .false.)
  stateMollifierUniformInTime = getOption(                                                  &
       "time_splitting/state_mollifier/uniform_in_time", .true.)

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

  ! Allocate per-segment accumulators.
  allocate(segmentCost(0:Nsplit-1))
  allocate(segmentL2sq(0:Nsplit-1))
  segmentCost     = 0.0_wp
  segmentL2sq     = 0.0_wp
  L2sqSum         = 0.0_wp
  solutionCrashed = .false.

  ! Scratch buffer for the end state of segment k-1 / the diff (ic_k - end_{k-1}).
  ! With state mollifier, the inner product is w-weighted; without, it is the SBP norm.
  allocate(scratch(size(region%grids)))
  do i = 1, size(region%grids)
    allocate(scratch(i)%F(region%grids(i)%nGridPoints, region%solverOptions%nUnknowns))
  end do

  ! Auxiliary buffer for state-mollifier IC blending. Only allocated when needed.
  if (useStateMollifier) then
    allocate(blendBuf(size(region%grids)))
    do i = 1, size(region%grids)
      allocate(blendBuf(i)%F(region%grids(i)%nGridPoints, region%solverOptions%nUnknowns))
    end do
  end if

  ! Segment loop: forward solve + on-the-fly matching-condition penalty.
  do k = 0, Nsplit-1
    ! Per-segment IC file: <prefix>-<k>.ic.q. ic_0.ic.q must be staged
    ! (run_msgrad.sh copies the original .ic.q into it); ic_k for k>=1 comes
    ! from the warmup trajectory or msgrad_test's perturbation.
    write(icFilename, '(2A,I0,A)') trim(outputPrefix), "-", k, ".ic.q"
    if (procRank == 0) inquire(file = icFilename, exist = fileExists)
    call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
    if (.not. fileExists) then
      write(message, '(3A)') "Per-segment IC file ", trim(icFilename), " does not exist!"
      call gracefulExit(MPI_COMM_WORLD, message)
    end if

    write(message, '(A,I0,A,I0,A)') "=== msforward: segment ", k, " of ", Nsplit, " ==="
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

    ! For k > 0: load ic_k now, form diff against scratch (= end_{k-1}), compute L2sq,
    ! and (with state mollifier) build the blended IC the solver will actually run from.
    ! runForward will re-load the IC immediately after, so this load is purely diagnostic
    ! in the no-mollifier path; with mollifier we point the dict to a shadow file.
    if (k > 0) then
      call region%loadData(QOI_FORWARD_STATE, icFilename)

      if (useStateMollifier) then
        if (.not. stateMollifierUniformInTime) then
          write(stateMollifierFilename, '(2A,I0,A)')                                        &
               trim(outputPrefix), "-", k, ".ic_mollifier.f"
          call region%loadData(QOI_STATE_MOLLIFIER, stateMollifierFilename)
        end if

        ! Stash ic_k (we are about to overwrite conservedVariables with the blend).
        do i = 1, size(region%states)
          blendBuf(i)%F = region%states(i)%conservedVariables
        end do

        ! Mollifier-weighted L2sq of the un-blended raw mismatch:
        !   L2sq = <w · (ic_k - end_{k-1}), (ic_k - end_{k-1})>_SBP
        L2sq = 0.0_wp
        do i = 1, size(region%states)
          scratch(i)%F = blendBuf(i)%F - scratch(i)%F            ! = ic_k - end_{k-1}
          L2sq = L2sq +                                                                     &
               region%grids(i)%computeInnerProduct(scratch(i)%F, scratch(i)%F,               &
                                                   region%grids(i)%stateMollifier(:,1))
        end do

        ! Build the blended IC the solver will see:
        !   ic_blended_k = w · ic_k + (1-w) · end_{k-1}
        !                = blendBuf - (1-w) · (blendBuf - end_{k-1})
        !                = blendBuf - (1-w) · scratch          (scratch now = ic_k - end_{k-1})
        do i = 1, size(region%states)
          do m = 1, region%solverOptions%nUnknowns
            region%states(i)%conservedVariables(:,m) =                                      &
                 blendBuf(i)%F(:,m) -                                                        &
                 (1.0_wp - region%grids(i)%stateMollifier(:,1)) * scratch(i)%F(:,m)
          end do
        end do

        ! Persist the blended IC; route runForward's loadInitialCondition there.
        write(blendedIcFilename, '(2A,I0,A)')                                               &
             trim(outputPrefix), "-", k, ".ic.blended.q"
        call region%saveData(QOI_FORWARD_STATE, blendedIcFilename)
      else
        ! Plain SBP-norm L2sq (no blending).
        L2sq = 0.0_wp
        do i = 1, size(region%states)
          scratch(i)%F = region%states(i)%conservedVariables - scratch(i)%F
          L2sq = L2sq +                                                                     &
               region%grids(i)%computeInnerProduct(scratch(i)%F, scratch(i)%F)
        end do
      end if

      if (region%commGridMasters /= MPI_COMM_NULL)                                          &
           call MPI_Allreduce(MPI_IN_PLACE, L2sq, 1, SCALAR_TYPE_MPI, MPI_SUM,              &
                              region%commGridMasters, ierror)
      do i = 1, size(region%grids)
        call MPI_Bcast(L2sq, 1, SCALAR_TYPE_MPI, 0, region%grids(i)%comm, ierror)
      end do
      segmentL2sq(k) = L2sq
      L2sqSum        = L2sqSum + L2sq
    end if

    ! Route runForward through its own loadInitialCondition path (dict mutation
    ! instead of restartFilename) and mutate solver%outputPrefix so all .q
    ! snapshots are written under <prefix>-<k>-<ts:08d>.q. This isolates each
    ! segment's start and end-of-segment files so the reverse-time migrator
    ! inside msadjoint loads ic_k -- not end_{k-1} from the previous segment --
    ! when it re-runs forward over segment k.
    !
    ! With state-mollifier blending we point at the shadow blended IC; otherwise
    ! at the raw <prefix>-<k>.ic.q (which Python wrote and msadjoint will reread).
    !
    ! Mutating solver%outputPrefix here does NOT change the actuator's
    ! .control_forcing_<name>.dat / .gradient_<name>.dat paths -- those were
    ! resolved once during setupActuatorPatch (src/ActuatorPatchImpl.f90:52,68)
    ! via getOption("output_prefix", PROJECT_NAME) and cached on the patch as
    ! controlForcingFilename / gradientFilename. They stay under the global
    ! prefix, so msforward/msadjoint can keep accumulating into a single .dat
    ! across segments via controlTimestepOffset = k*Nts.
    if (useStateMollifier .and. k > 0) then
      dict(icDictIndex)%val = trim(blendedIcFilename)
    else
      dict(icDictIndex)%val = trim(icFilename)
    end if
    write(segPrefix, '(2A,I0)') trim(outputPrefix), "-", k
    solver%outputPrefix = trim(segPrefix)

    segmentCost(k) = solver%runForward(region, controlTimestepOffset = k * Nts)

    ! runForward returns HUGE(0.0_wp) when checkSolutionLimits aborts a
    ! mid-segment step. Bail out: subsequent segments would propagate the
    ! crashed state and pollute JtimeIntegral / Jpenalty.
    if (segmentCost(k) > 0.5_wp * HUGE(0.0_wp)) then
      solutionCrashed = .true.
      write(message, '(A,I0,A)') "msforward: segment ", k,                                  &
           " returned HUGE (runForward solution-limit crash); skipping remaining segments."
      call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
      exit
    end if

    ! Cache the end-state in scratch for the next segment's penalty computation.
    ! The end-of-segment snapshot has already been written by showProgress under
    ! the segment-scoped prefix at <prefix>-<k>-<(k+1)*Nts>.q.
    do i = 1, size(region%states)
      scratch(i)%F = region%states(i)%conservedVariables
    end do

    call MPI_Barrier(MPI_COMM_WORLD, ierror)
  end do

  do i = 1, size(scratch)
    SAFE_DEALLOCATE(scratch(i)%F)
  end do
  SAFE_DEALLOCATE(scratch)

  if (allocated(blendBuf)) then
    do i = 1, size(blendBuf)
      SAFE_DEALLOCATE(blendBuf(i)%F)
    end do
    SAFE_DEALLOCATE(blendBuf)
  end if

  ! Aggregate. Penalty is now a single norm applied to the summed mismatch L2sqSum
  ! rather than a per-segment sum; matches optimization_ver3/penalty_norm.py.
  if (solutionCrashed) then
    J             = HUGE(0.0_wp)
    JtimeIntegral = HUGE(0.0_wp)
    Jpenalty      = 0.0_wp
    write(message, '(A)')                                                                  &
         'msforward: hard crash detected; writing HUGE to outputFilename.'
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  else
    JtimeIntegral = 0.0_wp
    do k = 0, Nsplit-1
      JtimeIntegral = JtimeIntegral + segmentCost(k)
    end do
    Jpenalty = penaltyWeight * penalty%norm(L2sqSum)
    J        = JtimeIntegral + Jpenalty

    write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))') 'msforward: time-integral J  = ', &
                                                            JtimeIntegral
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
    write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))') 'msforward: matching penalty = ', &
                                                            Jpenalty
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
    write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))') 'msforward: total J          = ', &
                                                            J
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  ! Write outputs (rank 0).
  if (procRank == 0) then
    open(unit = getFreeUnit(fileUnit), file = trim(outputFilename), action='write',         &
         iostat = stat, status = 'replace')
    write(fileUnit, '(1X,SP,' // SCALAR_FORMAT // ')') J
    close(fileUnit)

    open(unit = getFreeUnit(fileUnit), file = trim(subOutputFilename), action='write',      &
         iostat = stat, status = 'replace')
    write(fileUnit, '(A)') "# segment  time_integral_J            L2sq_mismatch"
    do k = 0, Nsplit-1
      write(fileUnit, '(I8,2(1X,SP,' // SCALAR_FORMAT // '))')                              &
           k, segmentCost(k), segmentL2sq(k)
    end do
    close(fileUnit)
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  SAFE_DEALLOCATE(segmentCost)
  SAFE_DEALLOCATE(segmentL2sq)

  if (associated(penalty)) deallocate(penalty)
  nullify(penalty)

  call solver%cleanup()
  call region%cleanup()

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  call cleanupErrorHandler()
  call disconnectParentIfSpawned()
  call MPI_Finalize(ierror)

end program msforward
