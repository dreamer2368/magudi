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
  logical :: lookForInput = .false., lookForOutput = .false., lookForSubOutput = .false.,   &
              inputFlag = .false., outputFlag = .false., subOutputFlag = .false.,           &
              saveMetricsFlag = .false.
  character(len = STRING_LENGTH) :: argument, inputFilename, outputFilename, subOutputFilename
  character(len = STRING_LENGTH) :: filename, outputPrefix, segPrefix, message
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

  ! Adjoint segment loop: process segments in reverse order. The matching adjoint
  ! terminal for each segment's right boundary is (re)written immediately before
  ! that segment's runAdjoint call — it CANNOT be done as a one-shot pre-pass
  ! because runAdjoint's internal showProgress save (src/SolverImpl.f90:128-140)
  ! overwrites <prefix>-<ts:08d>.adjoint.q at every save_interval boundary,
  ! including the segment boundaries pre-pass writes target.
  do k = Nsplit-1, 0, -1
    write(message, '(A,I0,A,I0,A)') "=== msadjoint: segment ", k, " of ", Nsplit, " ==="
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

    ! Re-write the matching adjoint terminal at startTs + (k+1)*Nts:
    !   matching_k = penaltyWeight * (end_k - ic_{k+1})
    ! runAdjoint will load this when adjoint_nonzero_initial_condition = "true".
    ! Skip for k=Nsplit-1 (no right neighbor; runAdjoint synthesizes zero terminal).
    if (k < Nsplit - 1) then
      ! end_k lives under segment k's prefix (msforward wrote it via showProgress
      ! after running with solver%outputPrefix = <prefix>-<k>).
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
        scratch(i)%F = region%states(i)%conservedVariables
      end do

      write(icFilename, '(2A,I0,A)') trim(outputPrefix), "-", k+1, ".ic.q"
      if (procRank == 0) inquire(file = icFilename, exist = fileExists)
      call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
      if (.not. fileExists) then
        write(message, '(3A)') "Per-segment IC file ", trim(icFilename), " does not exist!"
        call gracefulExit(MPI_COMM_WORLD, message)
      end if
      call region%loadData(QOI_FORWARD_STATE, icFilename)
      do i = 1, size(region%states)
        region%states(i)%adjointVariables = penaltyWeight *                                 &
             (scratch(i)%F - region%states(i)%conservedVariables)
      end do

      ! Terminal file must live under segment k's prefix so runAdjoint
      ! (with solver%outputPrefix = <prefix>-<k> below) picks it up via its
      ! "<this%outputPrefix>-<region%timestep>.adjoint.q" path.
      write(terminalFilename, '(2A,I0,A,I8.8,A)') trim(outputPrefix), "-", k, "-",          &
           startTimestep + (k+1) * Nts, ".adjoint.q"
      call region%saveData(QOI_ADJOINT_STATE, terminalFilename)
    end if

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
      ! Load the matching adjoint terminal just written above.
      dict(nonzeroDictIndex)%val = "true"
    end if

    ! Mutate solver%outputPrefix to <prefix>-<k> so runAdjoint, showProgress, and
    ! the reverse-time migrator all key their I/O off segment k's namespace.
    ! That makes the migrator's "<prefix>-<startTimestep>.q" read pick up ic_k
    ! (written by msforward's segment k via runForward's IC-save) instead of
    ! end_{k-1} (which used to live at <prefix>-<k*Nts>.q under the global prefix).
    !
    ! Mutating solver%outputPrefix here does NOT change the actuator's
    ! .control_forcing_<name>.dat / .gradient_<name>.dat paths -- those were
    ! resolved once during setupActuatorPatch (src/ActuatorPatchImpl.f90:52,68)
    ! via getOption("output_prefix", PROJECT_NAME) and cached on the patch as
    ! controlForcingFilename / gradientFilename. They stay under the global
    ! prefix, so segments keep reading/writing one shared .dat with
    ! controlTimestepOffset = k*Nts selecting the segment's window.
    write(segPrefix, '(2A,I0)') trim(outputPrefix), "-", k
    solver%outputPrefix = trim(segPrefix)

    segmentCtrlIP(k) = solver%runAdjoint(region, controlTimestepOffset = k * Nts,           &
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

    ! end_{k-1} lives under segment (k-1)'s prefix.
    write(endFilename, '(2A,I0,A,I8.8,A)') trim(outputPrefix), "-", k-1, "-",               &
         startTimestep + k * Nts, ".q"
    call region%loadData(QOI_FORWARD_STATE, endFilename)
    do i = 1, size(region%states)
      region%states(i)%adjointVariables = region%states(i)%adjointVariables +               &
           penaltyWeight * (scratch(i)%F - region%states(i)%conservedVariables)
    end do

    ! Spatial inner product <g_ic_k, g_ic_k>_SBP of the finalized IC gradient, before save.
    L2sq = 0.0_wp
    do i = 1, size(region%states)
      L2sq = L2sq + region%grids(i)%computeInnerProduct(                                    &
           region%states(i)%adjointVariables, region%states(i)%adjointVariables)
    end do
    if (region%commGridMasters /= MPI_COMM_NULL)                                            &
         call MPI_Allreduce(MPI_IN_PLACE, L2sq, 1, SCALAR_TYPE_MPI, MPI_SUM,                &
                            region%commGridMasters, ierror)
    do i = 1, size(region%grids)
      call MPI_Bcast(L2sq, 1, SCALAR_TYPE_MPI, 0, region%grids(i)%comm, ierror)
    end do
    segmentICIP(k) = L2sq

    call region%saveData(QOI_ADJOINT_STATE, icAdjointFilename)
  end do

  do i = 1, size(scratch)
    SAFE_DEALLOCATE(scratch(i)%F)
  end do
  SAFE_DEALLOCATE(scratch)

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
