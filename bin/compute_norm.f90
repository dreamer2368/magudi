#include "config.h"

! compute_norm: builds the diagonal of the L-BFGS metric M_diag for the
! optim_ver4 concatenated control vector.
!
! Outputs:
!   <prefix>.norm_<actuator>.dat   one file per ACTUATOR row in bc.dat.
!                                   5D-subarray real64 stream matching the
!                                   .control_forcing_<actuator>.dat layout.
!                                   Diagonal entries = controllerNorm*dt -- the
!                                   natural L^2_t (x) L^2_x metric for forcing.
!   <prefix>.norm_ic.q             single multi-block PLOT3D solution file
!                                   shared by all Nsplit-1 IC slabs (the metric
!                                   is identical across segments by construction
!                                   of computeInnerProduct). Each conserved-
!                                   variable channel holds state_controllability
!                                   * grid%norm. Written via the standard
!                                   region%saveData(QOI_FORWARD_STATE, ...)
!                                   path so msforward/msadjoint's .ic.q
!                                   readers can ingest it identically.
!   <prefix>.layout.txt            rank-0 text descriptor (canonical slot order):
!                                       actuator <name> <size_in_real64>
!                                       ic <k> <size_in_real64>
!                                   Sizes obtained by inquire on the just-written
!                                   actuator files and computed for ic slabs.
!                                   Python loads this once at startup to build
!                                   its schema.

program compute_norm

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver

  use Grid_enum
  use State_enum

  use InputHelper, only : parseInputFile, getFreeUnit, getOption, getRequiredOption
  use InputHelperImpl, only: dict, find
  use ErrorHandler
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, dictIndex, procRank, numProcs, ierror
  integer :: kthArgument, numberOfArguments
  logical :: lookForInput = .false., lookForLayout = .false.,                                &
             inputFlag = .false., layoutFlag = .false., saveMetricsFlag = .false.
  character(len = STRING_LENGTH) :: argument, inputFilename, layoutFilename
  character(len = STRING_LENGTH) :: filename, outputPrefix, message
  logical :: fileExists, success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_Solver) :: solver

  ! Time-splitting parameters.
  integer :: Nsplit, Nts, startTimestep
  real(wp) :: stateControllability

  ! Initialize MPI.
  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)

  call initializeErrorHandler()

  numberOfArguments = command_argument_count()
  do kthArgument = 1, numberOfArguments
    call get_command_argument(kthArgument, argument)
    select case(trim(adjustl(argument)))
    case("--input")
      lookForInput = .true.
    case("--layout-output")
      lookForLayout = .true.
    case("--save_metrics")
      saveMetricsFlag = .true.
    case default
      if (lookForInput) then
        inputFilename = trim(adjustl(argument))
        if (procRank == 0) inquire(file = inputFilename, exist = fileExists)
        call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
        if (.not. fileExists) then
          write(message, '(3A)') "input file ", trim(inputFilename), " does not exist!"
          call gracefulExit(MPI_COMM_WORLD, message)
        end if
        lookForInput = .false.
        inputFlag = .true.
      elseif (lookForLayout) then
        layoutFilename = trim(adjustl(argument))
        lookForLayout = .false.
        layoutFlag = .true.
      else
        write(message, '(3A)') "option ", trim(argument), " unknown!"
        call gracefulExit(MPI_COMM_WORLD, message)
      end if
    end select
  end do

  call startTiming("total")

  if (.not. inputFlag) inputFilename = PROJECT_NAME // ".inp"
  call parseInputFile(inputFilename)

  ! adjoint run options (collectNorm requires enableAdjoint).
  call find("enable_adjoint_solver", dictIndex)
  if (dictIndex < 0) then
    call gracefulExit(MPI_COMM_WORLD,                                                        &
         "magudi.inp must contain 'enable_adjoint_solver' for compute_norm.")
  end if
  dict(dictIndex)%val = "true"

  outputPrefix = getOption("output_prefix", PROJECT_NAME)

  if (.not. layoutFlag) layoutFilename = trim(outputPrefix) // ".layout.txt"

  write(message, '(2A)') "Input file: ", trim(inputFilename)
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  write(message, '(2A)') "Layout file: ", trim(layoutFilename)
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

  ! Time-splitting parameters.
  call getRequiredOption("time_splitting/number_of_segments", Nsplit)
  call getRequiredOption("time_splitting/segment_length", Nts)
  call getRequiredOption("time_splitting/start_timestep", startTimestep)
  stateControllability = getOption("time_splitting/state_controllability", 1.0_wp)
  assert(Nsplit >= 1)
  assert(Nts >= 1)
  assert(startTimestep >= 0)
  assert(stateControllability > 0.0_wp)

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

  ! Pass 1: actuator norms (writes <prefix>.norm_<actuator>.dat per actuator patch).
  call collectControlSpaceNorm(solver, region)

  ! Pass 2: single shared IC norm file <prefix>.norm_ic.q. The metric is identical
  ! across Nsplit-1 IC slabs by construction (computeInnerProduct weights only the
  ! spatial index), so one .q file is loaded by every ic-owner rank in Python.
  call writeIcNormQ(region, outputPrefix, stateControllability)

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  ! Pass 3: rank-0 text layout descriptor, sizes via inquire on real files.
  call writeLayout(region, outputPrefix, layoutFilename, Nsplit)

  call solver%cleanup()
  call region%cleanup()

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

contains

  subroutine collectControlSpaceNorm(this, region)

    ! Lifted from bin/control_space_norm.f90. The only output side effect is
    ! that each ACTUATOR patch's gradientFilename is redirected to
    ! '<prefix>.norm_<name>.dat' before the time march, so controller%collectNorm
    ! writes the metric there.

    use iso_fortran_env, only : output_unit

    use Patch_mod, only : t_Patch
    use Region_mod, only : t_Region
    use Solver_mod, only : t_Solver
    use Controller_mod, only : t_Controller
    use Functional_mod, only : t_Functional
    use ActuatorPatch_mod, only : t_ActuatorPatch
    use TimeIntegrator_mod, only : t_TimeIntegrator

    use Region_enum, only : FORWARD, ADJOINT

    use SolverImpl, only : loadInitialCondition

    use MPITimingsHelper, only : startTiming, endTiming
    use ErrorHandler, only : writeAndFlush

    implicit none

    class(t_Solver) :: this
    class(t_Region) :: region

    integer, parameter :: wp = SCALAR_KIND
    character(len = STRING_LENGTH) :: filename, message
    class(t_Patch), pointer :: patch => null()
    class(t_TimeIntegrator), pointer :: timeIntegrator => null()
    class(t_Controller), pointer :: controller => null()
    class(t_Functional), pointer :: functional => null()
    integer :: i, j, timestep, startTimestepLocal, timemarchDirection
    real(wp) :: time, startTime, timeStepSize

    assert(region%simulationFlags%enableController)
    assert(region%simulationFlags%enableFunctional)
    assert(region%simulationFlags%enableAdjoint)

    call startTiming("collectControlSpaceNorm")

    call this%timeIntegratorFactory%connect(timeIntegrator)
    assert(associated(timeIntegrator))
    call this%controllerFactory%connect(controller)
    assert(associated(controller))
    call this%functionalFactory%connect(functional)
    assert(associated(functional))

    if (region%simulationFlags%steadyStateSimulation)                                        &
         call this%residualManager%setup("adjoint_residuals", region)

    call loadInitialCondition(this, region, FORWARD)
    controller%onsetTime = region%states(1)%time
    controller%duration = this%nTimesteps * region%solverOptions%timeStepSize

    startTimestepLocal = region%timestep + this%nTimesteps
    startTime = region%states(1)%time + this%nTimesteps * region%solverOptions%timeStepSize

    timemarchDirection = -1
    if (region%simulationFlags%steadyStateSimulation) timemarchDirection = 1

    ! Redirect the gradient filename per actuator so the accumulated norm lands in
    ! <prefix>.norm_<name>.dat instead of <prefix>.gradient_<name>.dat.
    if (allocated(region%patchFactories)) then
      do i = 1, size(region%patchFactories)
        call region%patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        do j = 1, size(region%states)
          if (patch%gridIndex /= region%grids(j)%index .or. patch%nPatchPoints <= 0) cycle
          select type (patch)
          class is (t_ActuatorPatch)
            write(filename, '(4A)') trim(this%outputPrefix), ".norm_",                       &
                 trim(patch%name), ".dat"
            patch%gradientFilename = trim(filename)
          end select
        end do
      end do
    end if

    call controller%hookBeforeTimemarch(region, ADJOINT)

    if (this%probeInterval > 0) call region%resetProbes()

    time = startTime

    do timestep = startTimestepLocal + sign(1, timemarchDirection),                          &
         startTimestepLocal + sign(this%nTimesteps, timemarchDirection), timemarchDirection

       region%timestep = timestep
       timeStepSize = region%getTimeStepSize()

       do i = timeIntegrator%nStages, 1, -1
          call controller%collectNorm(region, timeIntegrator%norm(i) * timeStepSize)
       end do

       write(message, '(A,I8)') 'Collected control space norm at the time step = ', timestep
       call writeAndFlush(region%comm, output_unit, message)
    end do

    if (controller%controllerSwitch) call controller%hookAfterTimemarch(region, FORWARD)
    call controller%hookAfterTimemarch(region, ADJOINT)

    call this%residualManager%cleanup()

    write(message, '(A)') 'Control space norm collection is finished.'
    call writeAndFlush(region%comm, output_unit, message)

    call endTiming("collectControlSpaceNorm")

  end subroutine collectControlSpaceNorm

  subroutine writeIcNormQ(region, outputPrefix, controllabilityWeight)

    ! Writes <prefix>.norm_ic.q -- a single multi-block PLOT3D solution file
    ! shared by all Nsplit-1 IC slabs in Python's M_diag schema. The metric
    ! is identical across segments by construction of computeInnerProduct
    ! (src/GridImpl.f90:1142-1148): the SBP norm weights only the spatial
    ! index, so the per-segment diagonal entries are
    !   controllabilityWeight * grid%norm(p, 1)
    ! for every conserved-variable channel u.
    !
    ! We reuse region%saveData(QOI_FORWARD_STATE, ...) -- the same code path
    ! msforward calls when staging .q snapshots -- so the on-disk layout is
    ! identical to what ParallelIOHandler will read via plot3dnasa.

    use Region_mod, only : t_Region
    use State_enum, only : QOI_FORWARD_STATE

    implicit none

    class(t_Region) :: region
    character(len = *), intent(in) :: outputPrefix
    real(wp), intent(in) :: controllabilityWeight

    integer :: i, u, nUnknowns
    character(len = STRING_LENGTH) :: filename, message

    nUnknowns = region%solverOptions%nUnknowns

    ! Overwrite conservedVariables on every grid block with the broadcast metric.
    ! collectControlSpaceNorm has finished by now, so the state arrays are no
    ! longer needed for time stepping.
    do i = 1, size(region%states)
      if (region%grids(i)%nGridPoints <= 0) cycle
      do u = 1, nUnknowns
        region%states(i)%conservedVariables(:, u) =                                          &
             controllabilityWeight * region%grids(i)%norm(:, 1)
      end do
    end do

    write(filename, '(2A)') trim(outputPrefix), ".norm_ic.q"
    call region%saveData(QOI_FORWARD_STATE, filename)

    write(message, '(2A)') 'Wrote shared IC norm slab to ', trim(filename)
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

  end subroutine writeIcNormQ

  subroutine writeLayout(region, outputPrefix, layoutPath, Nsplit)

    ! Rank-0 text descriptor. One line per slot, in canonical order Python uses:
    !   actuator <name> <size_in_real64>     (one row per ACTUATOR patch in bc.dat;
    !                                          size from inquire on .norm_<name>.dat)
    !   ic <k> <size_in_real64>              (k = 1..Nsplit-1; size computed from
    !                                          total grid points * nUnknowns)
    !
    ! Sizes are in real64 elements. Rank 0 has the canonical view of patch metadata.

    use Region_mod, only : t_Region
    use Patch_mod, only : t_Patch
    use ActuatorPatch_mod, only : t_ActuatorPatch

    implicit none

    class(t_Region), intent(in) :: region
    character(len = *), intent(in) :: outputPrefix, layoutPath
    integer, intent(in) :: Nsplit

    integer :: i, k, fileUnit, stat, procRank, ierror
    integer(kind = MPI_OFFSET_KIND) :: fileSizeBytes, icSize
    class(t_Patch), pointer :: patch => null()
    character(len = STRING_LENGTH) :: dataPath, name

    call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
    if (procRank /= 0) return

    ! Per-segment ic slab size = sum_b(product(globalSize_b)) * nUnknowns.
    icSize = 0_MPI_OFFSET_KIND
    do i = 1, size(region%grids)
      icSize = icSize + product(int(region%grids(i)%globalSize, MPI_OFFSET_KIND))
    end do
    icSize = icSize * int(region%solverOptions%nUnknowns, MPI_OFFSET_KIND)

    open(unit = getFreeUnit(fileUnit), file = trim(layoutPath), action = 'write',            &
         iostat = stat, status = 'replace')

    write(fileUnit, '(A)') "# kind identifier size_in_real64"

    if (allocated(region%patchFactories)) then
      do i = 1, size(region%patchFactories)
        call region%patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        select type (patch)
        class is (t_ActuatorPatch)
          name = trim(patch%name)
          write(dataPath, '(4A)') trim(outputPrefix), ".norm_", trim(name), ".dat"
          call inquireFileSize(dataPath, fileSizeBytes)
          write(fileUnit, '(A,1X,A,1X,I0)') "actuator", trim(name),                          &
               fileSizeBytes / int(SIZEOF_SCALAR, MPI_OFFSET_KIND)
        end select
      end do
    end if

    do k = 1, Nsplit - 1
      write(fileUnit, '(A,1X,I0,1X,I0)') "ic", k, icSize
    end do

    close(fileUnit)

  end subroutine writeLayout

  subroutine inquireFileSize(path, sizeBytes)
    character(len = *), intent(in) :: path
    integer(kind = MPI_OFFSET_KIND), intent(out) :: sizeBytes

    integer(kind = MPI_OFFSET_KIND) :: localSize
    logical :: exists

    inquire(file = trim(path), exist = exists, size = localSize)
    if (.not. exists) then
       call gracefulExit(MPI_COMM_WORLD,                                                     &
            "writeLayout: expected output file " // trim(path) // " missing.")
    end if
    sizeBytes = localSize
  end subroutine inquireFileSize

end program compute_norm
