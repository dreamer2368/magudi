#include "config.h"

program msforward

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver

  use Grid_enum
  use State_enum

  use InputHelper, only : parseInputFile, getFreeUnit, getOption, getRequiredOption
  use ErrorHandler
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers
  use MPIHelper, only : disconnectParentIfSpawned

  implicit none

  type :: t_StateBuffer
     SCALAR_TYPE, allocatable :: F(:,:)
  end type t_StateBuffer

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, k, stat, fileUnit, procRank, numProcs, ierror
  integer :: kthArgument, numberOfArguments
  logical :: lookForInput = .false., lookForOutput = .false., lookForSubOutput = .false.,   &
              inputFlag = .false., outputFlag = .false., subOutputFlag = .false.,           &
              saveMetricsFlag = .false.
  character(len = STRING_LENGTH) :: argument, inputFilename, outputFilename, subOutputFilename
  character(len = STRING_LENGTH) :: filename, outputPrefix, message
  character(len = STRING_LENGTH) :: icFilename, endFilename
  logical :: fileExists, success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_Solver) :: solver
  ! After segment k completes, scratch(:)%F holds end_k for use in segment k+1's penalty.
  type(t_StateBuffer), allocatable :: scratch(:)

  ! << time-splitting parameters >>
  integer :: Nsplit, Nts, startTimestep
  real(wp) :: penaltyWeight

  ! << per-segment accumulators >>
  SCALAR_TYPE, allocatable :: segmentCost(:), segmentL2sq(:), segmentPenalty(:)
  SCALAR_TYPE :: J, JtimeIntegral, Jpenalty, L2sq

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

  ! Allocate per-segment accumulators.
  allocate(segmentCost(0:Nsplit-1))
  allocate(segmentL2sq(0:Nsplit-1))
  allocate(segmentPenalty(0:Nsplit-1))
  segmentCost     = 0.0_wp
  segmentL2sq     = 0.0_wp
  segmentPenalty  = 0.0_wp

  ! Scratch buffer for the end state of segment k-1 / the diff (ic_k - end_{k-1}).
  ! No state mollifier in this first step; the M-weighted inner product is the SBP norm.
  allocate(scratch(size(region%grids)))
  do i = 1, size(region%grids)
    allocate(scratch(i)%F(region%grids(i)%nGridPoints, region%solverOptions%nUnknowns))
  end do

  ! Segment loop: forward solve + on-the-fly matching-condition penalty.
  do k = 0, Nsplit-1
    write(icFilename, '(2A,I0,A)') trim(outputPrefix), "-", k, ".ic.q"
    if (procRank == 0) inquire(file = icFilename, exist = fileExists)
    call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
    if (.not. fileExists) then
      write(message, '(3A)') "Per-segment IC file ", trim(icFilename), " does not exist!"
      call gracefulExit(MPI_COMM_WORLD, message)
    end if

    write(message, '(A,I0,A,I0,A)') "=== msforward: segment ", k, " of ", Nsplit, " ==="
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

    ! For k > 0: load ic_k now, form diff against scratch (= end_{k-1}), compute L2sq.
    ! runForward will re-load the IC immediately after, so this load is purely diagnostic.
    if (k > 0) then
      call region%loadData(QOI_FORWARD_STATE, icFilename)
      L2sq = 0.0_wp
      do i = 1, size(region%states)
        scratch(i)%F = region%states(i)%conservedVariables - scratch(i)%F
        L2sq = L2sq +                                                                       &
             region%grids(i)%computeInnerProduct(scratch(i)%F, scratch(i)%F)
      end do
      if (region%commGridMasters /= MPI_COMM_NULL)                                          &
           call MPI_Allreduce(MPI_IN_PLACE, L2sq, 1, SCALAR_TYPE_MPI, MPI_SUM,              &
                              region%commGridMasters, ierror)
      do i = 1, size(region%grids)
        call MPI_Bcast(L2sq, 1, SCALAR_TYPE_MPI, 0, region%grids(i)%comm, ierror)
      end do
      segmentL2sq(k)    = L2sq
      segmentPenalty(k) = 0.5_wp * penaltyWeight * L2sq
    end if

    segmentCost(k) = solver%runForward(region, restartFilename = icFilename,                &
                                       referenceTimestep = k * Nts)

    ! Save the end-of-segment snapshot for msadjoint to consume.
    write(endFilename, '(2A,I8.8,A)') trim(outputPrefix), "-", region%timestep, ".q"
    call region%saveData(QOI_FORWARD_STATE, endFilename)

    ! Cache the end-state in scratch for the next segment's penalty computation.
    do i = 1, size(region%states)
      scratch(i)%F = region%states(i)%conservedVariables
    end do

    call MPI_Barrier(MPI_COMM_WORLD, ierror)
  end do

  do i = 1, size(scratch)
    SAFE_DEALLOCATE(scratch(i)%F)
  end do
  SAFE_DEALLOCATE(scratch)

  ! Aggregate.
  JtimeIntegral = 0.0_wp
  Jpenalty      = 0.0_wp
  do k = 0, Nsplit-1
    JtimeIntegral = JtimeIntegral + segmentCost(k)
    Jpenalty      = Jpenalty      + segmentPenalty(k)
  end do
  J = JtimeIntegral + Jpenalty

  write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))') 'msforward: time-integral J  = ',   &
                                                          JtimeIntegral
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))') 'msforward: matching penalty = ',   &
                                                          Jpenalty
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))') 'msforward: total J          = ',   &
                                                          J
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  ! Write outputs (rank 0).
  if (procRank == 0) then
    open(unit = getFreeUnit(fileUnit), file = trim(outputFilename), action='write',         &
         iostat = stat, status = 'replace')
    write(fileUnit, '(1X,SP,' // SCALAR_FORMAT // ')') J
    close(fileUnit)

    open(unit = getFreeUnit(fileUnit), file = trim(subOutputFilename), action='write',      &
         iostat = stat, status = 'replace')
    write(fileUnit, '(A)') "# segment  time_integral_J            L2sq_mismatch              penalty"
    do k = 0, Nsplit-1
      write(fileUnit, '(I8,3(1X,SP,' // SCALAR_FORMAT // '))')                              &
           k, segmentCost(k), segmentL2sq(k), segmentPenalty(k)
    end do
    close(fileUnit)
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  SAFE_DEALLOCATE(segmentCost)
  SAFE_DEALLOCATE(segmentL2sq)
  SAFE_DEALLOCATE(segmentPenalty)

  call solver%cleanup()
  call region%cleanup()

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  call cleanupErrorHandler()
  call disconnectParentIfSpawned()
  call MPI_Finalize(ierror)

end program msforward
