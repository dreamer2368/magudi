#include "config.h"

program terminal_objective

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver

  use Grid_enum
  use State_enum
  use Region_enum, only : FORWARD, ADJOINT

  use InputHelper, only : parseInputFile, getFreeUnit, getOption, getRequiredOption
  use InputHelperImpl, only: dict, find
  use ErrorHandler
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, stat, fileUnit, dictIndex, procRank, numProcs, ierror
  integer :: kthArgument, numberOfArguments
  logical :: lookForInput = .false., lookForOutput = .false., lookForMode = .false.,          &
              inputFlag = .false., outputFlag = .false., saveMetricsFlag = .false.,           &
              modeFlag = .false.
  integer :: mode
  character(len = STRING_LENGTH) :: argument, inputFilename, outputFilename
  character(len = STRING_LENGTH) :: filename, outputPrefix, message
  logical :: adjointRestart, fileExists, success
  integer :: accumulatedNTimesteps
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_Solver) :: solver

  ! << output variables >>
  integer :: inputNumber, simulationNumber
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
    case("--save_metrics")
      saveMetricsFlag = .true.
    case("--mode")
      lookForMode = .true.
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
      elseif (lookForMode) then
        select case(trim(adjustl(argument)))
        case("forward")
          mode = FORWARD
        case("adjoint")
          mode = ADJOINT
        case default
          write(message, '(3A)') "Mode ", trim(argument), " unknown!"
          call gracefulExit(MPI_COMM_WORLD, message)
        end select
        modeFlag = .true.
      else
        write(message, '(3A)') "option ", trim(argument), " unknown!"
        call gracefulExit(MPI_COMM_WORLD, message)
      end if
    end select
  end do

  if (.not.modeFlag) then
     write(message, '(3A)') "Mode is not specified! Choose 'forward' or 'adjoint'."
     call gracefulExit(MPI_COMM_WORLD, message)
  end if

  call startTiming("total")

  if ( .not. inputFlag ) inputFilename = PROJECT_NAME // ".inp"
  ! Parse options from the input file.
  call parseInputFile(inputFilename)

  ! adjoint run options
  call find("enable_adjoint_solver",dictIndex)
  dict(dictIndex)%val = "true"

  outputPrefix = getOption("output_prefix", PROJECT_NAME)

  if ( .not. outputFlag ) outputFilename = trim(outputPrefix) // ".terminal_objective.txt"
  write(message, '(2A)') "Input file: ", trim(inputFilename)
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  write(message, '(2A)') "Output file: ", trim(outputFilename)
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

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
  if (saveMetricsFlag) then
    write(filename, '(2A)') trim(outputPrefix), ".Jacobian.f"
    call region%saveData(QOI_JACOBIAN, filename)
    write(filename, '(2A)') trim(outputPrefix), ".metrics.f"
    call region%saveData(QOI_METRICS, filename)
  end if

  ! Initialize the solver.
  call solver%setup(region, outputPrefix = outputPrefix)

  ! Save the target mollifier if using code-generated values.
  if (region%simulationFlags%enableFunctional) then
     filename = getOption("target_mollifier_file", "")
     if (len_trim(filename) == 0) call region%saveData(QOI_TARGET_MOLLIFIER,                 &
          trim(outputPrefix) // ".target_mollifier.f")
  end if

  dummyValue = runTerminalObjective(solver,region,mode)

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  if (mode == FORWARD .and. procRank == 0) then
    open(unit = getFreeUnit(fileUnit), file = trim(outputFilename), action='write',          &
      iostat = stat, status = 'replace')
    write(fileUnit, '(1X,SP,' // SCALAR_FORMAT // ')') dummyValue
    close(fileUnit)
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  call solver%cleanup()
  call region%cleanup()

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  ! Finalize MPI.
  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

contains

  function runTerminalObjective(this, region, mode) result(instantaneousCostFunctional)

    ! <<< External modules >>>
    use iso_fortran_env, only : output_unit

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use Region_mod, only : t_Region
    use Solver_mod, only : t_Solver
    use Functional_mod, only : t_Functional
    use CostTargetPatch_mod, only : t_CostTargetPatch

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE, QOI_TIME_AVERAGED_STATE
    use Region_enum, only : FORWARD, ADJOINT

    ! <<< Private members >>>
    use SolverImpl, only : loadInitialCondition
    use RegionImpl, only : computeRegionIntegral

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming
    use ErrorHandler, only : writeAndFlush, gracefulExit
    use InputHelper, only : getOption, getRequiredOption

    implicit none

    ! <<< Arguments >>>
    class(t_Solver) :: this
    class(t_Region) :: region
    integer, intent(in) :: mode

    ! <<< Result >>>
    SCALAR_TYPE :: instantaneousCostFunctional

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    character(len = STRING_LENGTH) :: filename, message
    class(t_Functional), pointer :: functional => null()
    class(t_Patch), pointer :: patch => null()
    integer :: i, j, timestep, startTimestep
    real(wp) :: time, startTime, timeStepSize
    logical :: controllerSwitch = .false., solutionCrashes = .false.

    call startTiming("runTerminalObjective")

    ! Connect to the previously allocated functional.
    assert(region%simulationFlags%enableFunctional)
    call this%functionalFactory%connect(functional)
    assert(associated(functional))
    functional%runningTimeQuadrature = 0.0_wp
    instantaneousCostFunctional = 0.0_wp

    ! Load the initial condition.
    call loadInitialCondition(this, region, FORWARD)

    if (region%simulationFlags%enableBodyForce) then
      call getRequiredOption("body_force/initial_momentum", region%initialXmomentum)
      region%oneOverVolume = computeRegionIntegral(region)
      region%initialXmomentum = region%initialXmomentum * region%oneOverVolume
      region%oneOverVolume = 1.0_wp / region%oneOverVolume
      region%momentumLossPerVolume = 0.0_wp
    end if

    startTimestep = region%timestep
    startTime = region%states(1)%time

    write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-",                            &
         startTimestep + this%nTimesteps, ".q"
    call region%loadData(QOI_FORWARD_STATE, filename)
    do i = 1, size(region%states) !... update state
       call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
            region%solverOptions)
    end do

    select case(mode)
    case(FORWARD)

      instantaneousCostFunctional = functional%compute(region)

      write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))') 'Forward run: terminal objective = ', &
                                                             instantaneousCostFunctional
      call writeAndFlush(region%comm, output_unit, message)

    case(ADJOINT)

      region%states(:)%adjointForcingFactor = - 1.0_wp !... adjoint forcing is defined with negative sign for backward time-integration.
      call functional%updateAdjointForcing(region,.false.) !...always non-zero adjoint forcing.

      do i = 1, size(region%states)
         region%states(i)%rightHandSide = 0.0_wp
      end do

      do i = 1, size(region%patchFactories)
         call region%patchFactories(i)%connect(patch)
         if (.not. associated(patch)) cycle
         select type (patch)
            class is (t_CostTargetPatch)
            do j = 1, size(region%states)
               if (patch%gridIndex /= region%grids(j)%index) cycle
               call patch%updateRhs(ADJOINT, region%simulationFlags,                     &
                    region%solverOptions, region%grids(j), region%states(j))
            end do
         end select
      end do

      do i = 1, size(region%states)
         region%states(i)%adjointVariables = region%states(i)%rightHandSide
      end do

      write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", region%timestep, ".adjoint.q"
      call region%saveData(QOI_ADJOINT_STATE, filename)

    end select

    call endTiming("runTerminalObjective")

  end function runTerminalObjective

end program terminal_objective
