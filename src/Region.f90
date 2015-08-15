#include "config.h"

module Region_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  use MPI, only : MPI_COMM_NULL

  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use Patch_factory, only : t_PatchFactory
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none
  private

  type, public :: t_Region

     type(t_Grid), allocatable :: grids(:)
     type(t_State), allocatable :: states(:)
     type(t_PatchFactory), allocatable :: patchFactories(:)
     type(t_SolverOptions) :: solverOptions
     type(t_SimulationFlags) :: simulationFlags
     integer :: comm = MPI_COMM_NULL, commGridMasters = MPI_COMM_NULL, timestep = 0
     integer, allocatable :: globalGridSizes(:,:), processDistributions(:,:),                &
          gridCommunicators(:)
     character(len = STRING_LENGTH), allocatable :: patchNames(:), patchTypes(:)
     integer, allocatable :: patchData(:,:), patchCommunicators(:), patchMasterRanks(:),     &
          patchInterfaces(:), interfaceIndexReorderings(:,:)
     logical :: outputOn = .true.

   contains

     procedure, pass :: setup
     procedure, pass :: cleanup
     procedure, pass :: setupBoundaryConditions
     procedure, pass :: loadData
     procedure, pass :: saveData
     procedure, pass :: computeCfl
     procedure, pass :: computeTimeStepSize
     procedure, pass :: reportGridDiagnostics
     procedure, pass :: computeRhs
     procedure, pass :: saveSpongeStrength
     procedure, pass :: resetProbes
     procedure, pass :: saveProbeData

  end type t_Region

contains

  subroutine setup(this, comm, globalGridSizes, simulationFlags, solverOptions, verbose)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Internal modules >>>
    use InputHelper, only : getRequiredOption
    use MPITimingsHelper, only : startTiming, endTiming

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: this
    integer, intent(in) :: comm, globalGridSizes(:,:)
    type(t_SimulationFlags), intent(in), optional :: simulationFlags
    type(t_SolverOptions), intent(in), optional :: solverOptions
    logical, intent(in), optional :: verbose

    ! <<< Local variables >>>
    integer :: i, j, color, procRank, numProcs, ierror
    logical :: verbose_
    character(len = STRING_LENGTH) :: decompositionMapFilename

    call startTiming("Setup region")

    ! Clean slate.
    call this%cleanup()
    this%comm = comm
    call MPI_Comm_size(this%comm, numProcs, ierror)

    verbose_ = .true.
    if (present(verbose)) verbose_ = verbose

    SAFE_DEALLOCATE(this%globalGridSizes)
    SAFE_DEALLOCATE(this%processDistributions)

    allocate(this%globalGridSizes(size(globalGridSizes, 1), size(globalGridSizes, 2)),       &
         source = globalGridSizes)
    allocate(this%processDistributions(size(this%globalGridSizes, 1),                        &
         size(this%globalGridSizes, 2)), source = 0)

    ! Initialize simulation flags.
    if (present(simulationFlags)) then
       this%simulationFlags = simulationFlags
    else
       call this%simulationFlags%initialize()
    end if

    ! Initialize solver options.
    if (present(solverOptions)) then
       this%solverOptions = solverOptions
    else
       call this%solverOptions%initialize(size(this%globalGridSizes, 1),                     &
            this%simulationFlags, this%comm)
    end if

    ! Distribute the grids between available MPI processes.
    if (this%simulationFlags%manualDomainDecomp .and.                                        &
         numProcs > size(this%globalGridSizes, 2)) then
       call getRequiredOption("decomposition_map_file",                                      &
            decompositionMapFilename, this%comm)
       call readDecompositionMap(this, decompositionMapFilename)
    end if
    call distributeGrids(this, verbose_)
    call MPI_Barrier(this%comm, ierror)

    ! Setup grids.
    allocate(this%grids(count(this%gridCommunicators /= MPI_COMM_NULL)))
    this%grids(:)%index = 0
    do i = 1, size(this%globalGridSizes, 2)
       do j = 1, size(this%grids)
          if (this%grids(j)%index /= 0) cycle
          if (this%gridCommunicators(i) /= MPI_COMM_NULL) then
             call this%grids(j)%setup(i, this%globalGridSizes(:,i),                          &
                  this%gridCommunicators(i), this%processDistributions(:,i),                 &
                  simulationFlags = this%simulationFlags)
             exit
          end if
       end do
       call MPI_Barrier(this%comm, ierror)
    end do

    ! Setup spatial discretization.
    do i = 1, size(this%grids)
       call this%grids(i)%setupSpatialDiscretization(this%simulationFlags, this%solverOptions)
    end do
    call MPI_Barrier(this%comm, ierror)

    ! Setup states.
    allocate(this%states(size(this%grids)))
    do i = 1, size(this%states)
       call this%states(i)%setup(this%grids(i), this%simulationFlags, this%solverOptions)
    end do
    call MPI_Barrier(this%comm, ierror)

    ! Create a communicator for master processes of grid-level communicators.
    assert(MPI_UNDEFINED /= 1)
    color = MPI_UNDEFINED
    do i = 1, size(this%grids)
       call MPI_Comm_rank(this%grids(i)%comm, procRank, ierror)
       if (procRank == 0) color = 1
    end do
    call MPI_Comm_rank(this%comm, procRank, ierror)
    call MPI_Comm_split(this%comm, color, procRank, this%commGridMasters, ierror)

    call endTiming("Setup region")

  end subroutine setup

  subroutine cleanup(this)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: this

    ! <<< Local variables >>>
    integer :: i, ierror
    class(t_Patch), pointer :: patch => null()

    if (allocated(this%grids)) then
       do i = 1, size(this%grids)
          call this%grids(i)%cleanup()
       end do
    end if
    SAFE_DEALLOCATE(this%grids)

    if (allocated(this%states)) then
       do i = 1, size(this%states)
          call this%states(i)%cleanup()
       end do
    end if
    SAFE_DEALLOCATE(this%states)

    if (allocated(this%patchFactories)) then
       do i = 1, size(this%patchFactories)
          call this%patchFactories(i)%connect(patch)
          if (associated(patch)) call patch%cleanup()
          call this%patchFactories(i)%cleanup()
       end do
    end if
    SAFE_DEALLOCATE(this%patchFactories)

    SAFE_DEALLOCATE(this%globalGridSizes)
    SAFE_DEALLOCATE(this%processDistributions)
    SAFE_DEALLOCATE(this%gridCommunicators)

    SAFE_DEALLOCATE(this%patchData)
    SAFE_DEALLOCATE(this%patchNames)
    SAFE_DEALLOCATE(this%patchTypes)
    SAFE_DEALLOCATE(this%patchCommunicators)
    SAFE_DEALLOCATE(this%patchMasterRanks)
    SAFE_DEALLOCATE(this%patchInterfaces)
    SAFE_DEALLOCATE(this%interfaceIndexReorderings)

    if (this%comm /= MPI_COMM_NULL .and. this%comm /= MPI_COMM_WORLD)                        &
         call MPI_Comm_free(this%comm, ierror)
    this%comm = MPI_COMM_NULL

    if (this%commGridMasters /= MPI_COMM_NULL)                                               &
         call MPI_Comm_free(this%commGridMasters, ierror)
    this%commGridMasters = MPI_COMM_NULL

    this%timestep = 0
    this%outputOn = .true.

  end subroutine cleanup

  subroutine setupBoundaryConditions(this, boundaryConditionFilename)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: this
    character(len = *), intent(in) :: boundaryConditionFilename

    ! <<< Local variables >>>
    integer :: i, j, k, nPatches, ierror
    integer, allocatable :: patchIndices(:)
    class(t_Patch), pointer :: patch => null()
    integer :: extent(6)

    call startTiming("Setup BCs")

    ! Cleanup previously allocated patch information.
    if (allocated(this%patchFactories)) then
       do i = 1, size(this%patchFactories)
          call this%patchFactories(i)%connect(patch)
          if (associated(patch)) call patch%cleanup()
          call this%patchFactories(i)%cleanup()
       end do
    end if
    SAFE_DEALLOCATE(this%patchFactories)
    SAFE_DEALLOCATE(this%patchData)
    SAFE_DEALLOCATE(this%patchNames)
    SAFE_DEALLOCATE(this%patchTypes)
    SAFE_DEALLOCATE(this%patchCommunicators)
    SAFE_DEALLOCATE(this%patchMasterRanks)
    SAFE_DEALLOCATE(this%patchInterfaces)
    SAFE_DEALLOCATE(this%interfaceIndexReorderings)

    call readBoundaryConditions(this, boundaryConditionFilename)
    call readPatchInterfaceInformation(this)
    call distributePatches(this)
    call MPI_Barrier(this%comm, ierror)

    if (allocated(this%patchData)) then

       nPatches = 0
       do i = 1, size(this%grids)
          nPatches = nPatches + count(this%patchData(:,1) == this%grids(i)%index)
       end do

       allocate(this%patchFactories(nPatches))
       allocate(patchIndices(nPatches), source = 0)

       do k = 1, size(this%grids)
          do i = 1, size(this%patchData, 1)
             if (this%patchData(i,1) /= this%grids(k)%index) cycle
             do j = 1, size(this%patchFactories)
                if (patchIndices(j) /= 0) cycle
                call this%patchFactories(j)%connect(patch, trim(this%patchTypes(i)))
                assert(associated(patch))
                extent = this%patchData(i,3:8)
                call patch%setup(this%patchNames(i), this%patchCommunicators(i),             &
                     this%grids(k), this%states(k), extent, this%patchData(i,2),             &
                     this%simulationFlags, this%solverOptions)
                patch%index = i
                patchIndices(j) = i
                exit
             end do
          end do
          call MPI_Barrier(this%grids(k)%comm, ierror)
       end do

       call computeSpongeStrengths(this)

       if (.not. this%simulationFlags%predictionOnly) then
          call normalizeControlMollifier(this)
          call normalizeTargetMollifier(this)
       end if

    end if

    SAFE_DEALLOCATE(patchIndices)

    call endTiming("Setup BCs")

  end subroutine setupBoundaryConditions

  subroutine loadData(this, quantityOfInterest, filename, speciesFilename)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env, only : output_unit

    ! <<< Derived types >>>
    use PLOT3DDescriptor_type, only : t_PLOT3DDescriptor

    ! <<< Enumerations >>>
    use Grid_enum
    use State_enum, only : QOI_DUMMY_FUNCTION, QOI_FORWARD_STATE

    ! <<< Internal modules >>>
    use ErrorHandler, only : gracefulExit, writeAndFlush
    use PLOT3DHelper
    use MPITimingsHelper, only : startTiming, endTiming

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: this
    integer, intent(in) :: quantityOfInterest
    character(len = *), intent(in) :: filename
    character(len = *), intent(in), optional :: speciesFilename

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    character(len = STRING_LENGTH) :: message, speciesFilename_
    logical :: success
    integer :: i, j, errorRank, procRank, ierror
    integer(kind = MPI_OFFSET_KIND) :: offset, speciesFileOffset
    logical :: isSolutionFile
    real(wp) :: solutionHeader(4)

    call startTiming("I/O (load)")

    isSolutionFile = .false.
    speciesFilename_ = ""
    if (present(speciesFilename)) speciesFilename_ = speciesFilename
    speciesFileOffset = int(0, MPI_OFFSET_KIND)

    select case(quantityOfInterest)
    case (QOI_GRID, QOI_JACOBIAN, QOI_TARGET_MOLLIFIER, QOI_CONTROL_MOLLIFIER,               &
         QOI_METRICS, QOI_DUMMY_FUNCTION)
    case default
       isSolutionFile = .true.
       if (this%solverOptions%nSpecies > 0 .and. len_trim(speciesFilename_) == 0 .and.       &
            filename(len_trim(filename)-1:len_trim(filename)) /= ".q") then
          write(message, '(3A)')                                                             &
               "Auto-detection of species filename failed: Solution file '", trim(filename), &
               "' does not have a '.q' extension!"
          call gracefulExit(this%comm, message)
       end if
       if (this%solverOptions%nSpecies > 0)                                                  &
            speciesFilename_ = filename(:len_trim(filename)-2) // ".f"
    end select

    if (present(speciesFilename)) then
       write(message, '(5A)') "Reading '", trim(filename), "', '",                           &
            trim(speciesFilename), "'..."
    else if (len_trim(speciesFilename_) > 0) then
       write(message, '(3A)') "Reading '", filename(:len_trim(filename)-2), ".{q,f}'..."
    else
       write(message, '(3A)') "Reading '", trim(filename), "'..."
    end if
    call writeAndFlush(this%comm, output_unit, message, advance = 'no')

    do i = 1, size(this%gridCommunicators)

       success = .true.

       do j = 1, size(this%grids)
          if (this%grids(j)%index == i) then !... read one grid at a time

             offset = plot3dGetOffset(this%gridCommunicators(i), filename, i, success)
             if (.not. success) exit

             if (len_trim(speciesFilename_) > 0) then
                speciesFileOffset = plot3dGetOffset(this%gridCommunicators(i),               &
                     speciesFilename_, i, success)
                if (.not. success) exit
             end if

             select case(quantityOfInterest)
             case (QOI_GRID, QOI_JACOBIAN, QOI_METRICS, QOI_TARGET_MOLLIFIER,                &
                  QOI_CONTROL_MOLLIFIER)
                call this%grids(j)%loadData(quantityOfInterest,                              &
                     trim(filename), offset, success)
             case default
                isSolutionFile = (quantityOfInterest /= QOI_DUMMY_FUNCTION)
                call this%states(j)%loadData(this%grids(j), quantityOfInterest,              &
                     trim(filename), offset, success,                                        &
                     speciesFilename = trim(speciesFilename_),                               &
                     speciesFileOffset = speciesFileOffset)
             end select

             exit
          end if
       end do

       call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL, MPI_LAND, this%comm, ierror)
       if (.not. success) exit
       call MPI_Barrier(this%comm, ierror)

    end do

    if (isSolutionFile) then
       if (.not. this%simulationFlags%steadyStateSimulation) then
          solutionHeader = real(this%states(1)%plot3dSolutionHeader, wp)
          call MPI_Bcast(solutionHeader, 4, SCALAR_TYPE_MPI, 0, this%comm, ierror)
          do i = 1, 4
             this%states(:)%plot3dSolutionHeader(i) = solutionHeader(i)
          end do
          if (quantityOfInterest == QOI_FORWARD_STATE) then
             this%states(:)%time = real(solutionHeader(4), wp)
             this%timestep = nint(real(solutionHeader(1), wp))
          end if
       else
          this%states(:)%time = 0.0_wp
          this%timestep = 0
       end if
    end if

    if (success) then
       write(message, '(A)') " done!"
    else
       write(message, '(A)') " failed!"
    end if
    call writeAndFlush(this%comm, output_unit, message)

    if (.not. success) then
       call MPI_Comm_rank(this%comm, procRank, ierror)
       errorRank = 0
       if (len_trim(plot3dErrorMessage) > 0) errorRank = procRank
       call MPI_Allreduce(MPI_IN_PLACE, errorRank, 1, MPI_INTEGER, MPI_MAX, this%comm, ierror)
       call MPI_Bcast(plot3dErrorMessage, STRING_LENGTH, MPI_CHARACTER,                      &
            errorRank, this%comm, ierror)
       call gracefulExit(this%comm, plot3dErrorMessage)
    end if

    call endTiming("I/O (load)")

  end subroutine loadData

  subroutine saveData(this, quantityOfInterest, filename)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env, only : output_unit

    ! <<< Derived types >>>
    use PLOT3DDescriptor_type

    ! <<< Enumerations >>>
    use Grid_enum
    use State_enum, only : QOI_DUMMY_FUNCTION

    ! <<< Internal modules >>>
    use ErrorHandler, only : writeAndFlush, gracefulExit
    use PLOT3DHelper
    use MPITimingsHelper, only : startTiming, endTiming

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: this
    integer, intent(in) :: quantityOfInterest
    character(len = *), intent(in) :: filename

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    character(len = STRING_LENGTH) :: message, speciesFilename
    logical :: success
    integer :: i, j, nScalars, errorRank, procRank, ierror
    integer(kind = MPI_OFFSET_KIND) :: offset, speciesFileOffset
    logical :: isSolutionFile

    if (.not. this%OutputOn) return

    call startTiming("I/O (save)")

    isSolutionFile = .false.
    speciesFilename = ""
    speciesFileOffset = int(0, MPI_OFFSET_KIND)

    select case(quantityOfInterest)
    case (QOI_GRID, QOI_JACOBIAN, QOI_TARGET_MOLLIFIER, QOI_CONTROL_MOLLIFIER,               &
         QOI_METRICS, QOI_DUMMY_FUNCTION)
    case default
       isSolutionFile = .true.
       if (filename(len_trim(filename)-1:len_trim(filename)) /= ".q") then
          write(message, '(A)') "Solution files must have extension '.q'!"
          call gracefulExit(this%comm, message)
       end if
       if (this%solverOptions%nSpecies > 0)                                                  &
            speciesFilename = filename(:len_trim(filename)-2) // ".f"
    end select

    if (len_trim(speciesFilename) > 0) then
       write(message, '(3A)') "Writing '", filename(:len_trim(filename)-2), ".{q,f}'..."
    else
       write(message, '(3A)') "Writing '", trim(filename), "'..."
    end if
    call writeAndFlush(this%comm, output_unit, message, advance = 'no')

    select case(quantityOfInterest)

    case (QOI_GRID)
       call plot3dWriteSkeleton(this%comm, trim(filename),                                   &
            PLOT3D_GRID_FILE, this%globalGridSizes, success)
    case (QOI_JACOBIAN, QOI_TARGET_MOLLIFIER, QOI_CONTROL_MOLLIFIER)
       call plot3dWriteSkeleton(this%comm, trim(filename),                                   &
            PLOT3D_FUNCTION_FILE, this%globalGridSizes, success, 1)
    case (QOI_METRICS)
       call plot3dWriteSkeleton(this%comm, trim(filename),                                   &
            PLOT3D_FUNCTION_FILE, this%globalGridSizes, success,                             &
            size(this%globalGridSizes, 1) ** 2)

    case (QOI_DUMMY_FUNCTION)

       nScalars = huge(1)
       do i = 1, size(this%states)
          assert(associated(this%states(i)%dummyFunction))
          nScalars = min(nScalars, size(this%states(i)%dummyFunction, 2))
       end do
       call MPI_Allreduce(MPI_IN_PLACE, nScalars, 1,                                         &
            MPI_INTEGER, MPI_MIN, this%comm, ierror)
#ifndef NDEBUG
       do i = 1, size(this%states)
          assert(size(this%states(i)%dummyFunction, 2) == nScalars)
       end do
#endif
       call plot3dWriteSkeleton(this%comm, trim(filename), PLOT3D_FUNCTION_FILE,             &
            this%globalGridSizes, success, nScalars)

    case default

       call plot3dWriteSkeleton(this%comm, trim(filename), PLOT3D_SOLUTION_FILE,             &
            this%globalGridSizes, success)
       if (len_trim(speciesFilename) > 0)                                                    &
            call plot3dWriteSkeleton(this%comm, trim(speciesFilename), PLOT3D_FUNCTION_FILE, &
            this%globalGridSizes, success, this%solverOptions%nSpecies)

    end select

    if (isSolutionFile) then
       if (.not. this%simulationFlags%steadyStateSimulation) then
          do i = 1, size(this%states)
             this%states(:)%plot3dSolutionHeader(1) =                                        &
                  real(this%timestep, wp)
          end do
       else
          do i = 1, size(this%states)
             this%states(:)%plot3dSolutionHeader(1) = 0.0_wp
          end do
       end if
    end if

    do i = 1, size(this%gridCommunicators)

       success = .true.

       do j = 1, size(this%grids)
          if (this%grids(j)%index == i) then !... read one grid at a time.

             offset = plot3dGetOffset(this%gridCommunicators(i), filename, i, success)
             if (.not. success) exit

             if (len_trim(speciesFilename) > 0) then
                speciesFileOffset = plot3dGetOffset(this%gridCommunicators(i),               &
                     speciesFilename, i, success)
                if (.not. success) exit
             end if

             select case(quantityOfInterest)
             case (QOI_GRID, QOI_JACOBIAN, QOI_METRICS, QOI_TARGET_MOLLIFIER,                &
                  QOI_CONTROL_MOLLIFIER)
                call this%grids(j)%saveData(quantityOfInterest,                              &
                     trim(filename), offset, success)
             case default
                call this%states(j)%saveData(this%grids(j), quantityOfInterest,              &
                     trim(filename), offset, success,                                        &
                     speciesFilename = trim(speciesFilename),                                &
                     speciesFileOffset = speciesFileOffset)
             end select

             exit
          end if
       end do

       call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL,                             &
            MPI_LAND, this%comm, ierror)
       if (.not. success) exit
       call MPI_Barrier(this%comm, ierror)

    end do

    if (success) then
       write(message, '(A)') " done!"
    else
       write(message, '(A)') " failed!"
    end if
    call writeAndFlush(this%comm, output_unit, message)

    if (.not. success) then
       call MPI_Comm_rank(this%comm, procRank, ierror)
       errorRank = 0
       if (len_trim(plot3dErrorMessage) > 0) errorRank = procRank
       call MPI_Allreduce(MPI_IN_PLACE, errorRank, 1, MPI_INTEGER, MPI_MAX, this%comm, ierror)
       call MPI_Bcast(plot3dErrorMessage, STRING_LENGTH, MPI_CHARACTER,                      &
            errorRank, this%comm, ierror)
       call gracefulExit(this%comm, plot3dErrorMessage)
    end if

    call endTiming("I/O (save)")

  end subroutine saveData

  function computeCfl(this) result(cfl)

    ! <<< External modules >>>
    use MPI

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: this

    ! <<< Result >>>
    real(SCALAR_KIND) :: cfl

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, ierror

    if (this%simulationFlags%useConstantCfl) then
       cfl = this%solverOptions%cfl
    else
       cfl = 0.0_wp
       do i = 1, size(this%states)
          cfl = max(cfl, this%states(i)%computeCfl(this%grids(i),                            &
               this%simulationFlags, this%solverOptions))
       end do
       call MPI_Allreduce(MPI_IN_PLACE, cfl, 1, SCALAR_TYPE_MPI, MPI_MAX, this%comm, ierror)
    end if

  end function computeCfl

  function computeTimeStepSize(this) result(timeStepSize)

    ! <<< External modules >>>
    use MPI

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: this

    ! <<< Result >>>
    real(SCALAR_KIND) :: timeStepSize

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, ierror

    if (this%simulationFlags%useConstantCfl) then
       timeStepSize = huge(0.0_wp)
       do i = 1, size(this%states)
          timeStepSize = min(timeStepSize,                                                   &
               this%states(i)%computeTimeStepSize(this%grids(i),                             &
               this%simulationFlags, this%solverOptions))
       end do
       call MPI_Allreduce(MPI_IN_PLACE, timeStepSize, 1, SCALAR_TYPE_MPI,                    &
            MPI_MIN, this%comm, ierror)
    else
       timeStepSize = this%solverOptions%timeStepSize
    end if

  end function computeTimeStepSize

  subroutine reportGridDiagnostics(this)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env, only : output_unit

    ! <<< Internal modules >>>
    use ErrorHandler, only : writeAndFlush

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: this

    ! <<< Local variables >>>
    integer :: i, j, iGlobal, jGlobal, kGlobal, procRank, ierror
    character(len = STRING_LENGTH) :: str
    real(SCALAR_KIND) :: minimumJacobian, maximumJacobian

    call MPI_Comm_rank(this%comm, procRank, ierror)

    call writeAndFlush(this%comm, output_unit, "")
    write(str, '(A,I0.0,A,I2,A)') "Region is ", size(this%globalGridSizes, 1),               &
         "D and has ", size(this%globalGridSizes, 2), " block(s):"
    call writeAndFlush(this%comm, output_unit, str)
    write(str, '(A)') repeat('-', 32)
    call writeAndFlush(this%comm, output_unit, str)

    do i = 1, size(this%globalGridSizes, 2)
       do j = 1, size(this%grids)
          if (this%grids(j)%index == i) then

             write(str, '(2X,A,I2,3(A,I4),A)') "Block ", i, ": ",                            &
                  this%grids(j)%globalSize(1), " x ",                                        &
                  this%grids(j)%globalSize(2), " x ",                                        &
                  this%grids(j)%globalSize(3), " points"
             call writeAndFlush(this%grids(j)%comm, output_unit, str)

             call this%grids(j)%findMinimum(this%grids(j)%jacobian(:,1),                     &
                  minimumJacobian, iGlobal, jGlobal, kGlobal)
             write(str, '(4X,A,(SS,ES9.2E2),3(A,I4),A)') "min. Jacobian = ",                 &
                  real(minimumJacobian, SCALAR_KIND), " at (",                               &
                  iGlobal, ", ", jGlobal, ", ", kGlobal, ")"
             call writeAndFlush(this%grids(j)%comm, output_unit, str)

             call this%grids(j)%findMaximum(this%grids(j)%jacobian(:,1),                     &
                  maximumJacobian, iGlobal, jGlobal, kGlobal)
             write(str, '(4X,A,(SS,ES9.2E2),3(A,I4),A)') "max. Jacobian = ",                 &
                  real(maximumJacobian, SCALAR_KIND), " at (",                               &
                  iGlobal, ", ", jGlobal, ", ", kGlobal, ")"
             call writeAndFlush(this%grids(j)%comm, output_unit, str)

          end if
          call MPI_Barrier(this%grids(j)%comm, ierror)
       end do
       call MPI_Barrier(this%comm, ierror)
    end do

    call writeAndFlush(this%comm, output_unit, "")

  end subroutine reportGridDiagnostics

  subroutine computeRhs(this, mode)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : FORWARD, ADJOINT

    ! <<< Internal modules >>>
    use RhsHelper, only : computeRhsForward, computeRhsAdjoint
    use MPITimingsHelper, only : startTiming, endTiming

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: this
    integer, intent(in) :: mode

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j
    class(t_Patch), pointer :: patch => null()

    call startTiming("Compute RHS")

    ! Semi-discrete right-hand-side operator.
    do i = 1, size(this%states)
       select case (mode)
       case (FORWARD)
          call computeRhsForward(this%simulationFlags, this%solverOptions,                   &
               this%grids(i), this%states(i), this%patchFactories)
       case (ADJOINT)
          call computeRhsAdjoint(this%simulationFlags, this%solverOptions,                   &
               this%grids(i), this%states(i))
       end select
    end do

    call exchangeInterfaceData(this)

    ! Multiply by Jacobian.
    do i = 1, size(this%states)
       do j = 1, this%solverOptions%nUnknowns
          this%states(i)%rightHandSide(:,j) = this%states(i)%rightHandSide(:,j) *            &
               this%grids(i)%jacobian(:,1)
       end do
    end do

    ! Add patch penalties.
    if (allocated(this%patchFactories)) then
       do i = 1, size(this%patchFactories)
          call this%patchFactories(i)%connect(patch)
          if (.not. associated(patch)) cycle
          do j = 1, size(this%states)
             if (patch%gridIndex /= this%grids(j)%index) cycle
             call patch%updateRhs(mode, this%simulationFlags, this%solverOptions,            &
                  this%grids(j), this%states(j))
          end do
       end do
    end if

    ! Source terms.
    if (mode == FORWARD) then
       do i = 1, size(this%states)
          call this%states(i)%addSources(this%grids(i))
       end do
    end if

    ! Zero out right-hand-side in holes.
    do i = 1, size(this%states)
       do j = 1, this%solverOptions%nUnknowns
          where (this%grids(i)%iblank == 0)
             this%states(i)%rightHandSide(:,j) = 0.0_wp
          end where
       end do
    end do

    call endTiming("Compute RHS")

  end subroutine computeRhs

  subroutine saveSpongeStrength(this, filename)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use SpongePatch_mod, only : t_SpongePatch

    ! <<< Enumerations >>>
    use State_enum, only : QOI_DUMMY_FUNCTION

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: this
    character(len = *), intent(in) :: filename

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    type :: t_IntermediateStorage
       real(wp), pointer :: buffer(:,:) => null()
    end type t_IntermediateStorage
    type(t_IntermediateStorage), allocatable :: data(:)
    integer :: i, j, ierror
    class(t_Patch), pointer :: patch => null()

    allocate(data(size(this%grids)))

    do i = 1, size(this%grids)

       allocate(data(i)%buffer(this%grids(i)%nGridPoints, 1), source = 0.0_wp)

       if (allocated(this%patchFactories)) then
          do j = 1, size(this%patchFactories)
             call this%patchFactories(j)%connect(patch)

             if (.not. associated(patch)) cycle
             if (patch%gridIndex /= this%grids(i)%index) cycle

             select type (patch)
                class is (t_SpongePatch)
                call patch%disperseAdd(patch%spongeStrength, data(i)%buffer(:,1))
             end select

          end do
       end if

       call MPI_Allreduce(MPI_IN_PLACE, data(i)%buffer, size(data(i)%buffer),                &
            SCALAR_TYPE_MPI, MPI_SUM, this%grids(i)%comm, ierror)
       this%states(i)%dummyFunction => data(i)%buffer

    end do

    call this%saveData(QOI_DUMMY_FUNCTION, trim(filename))

    do i = 1, size(data)
       if (associated(data(i)%buffer)) deallocate(data(i)%buffer)
       nullify(data(i)%buffer)
    end do
    SAFE_DEALLOCATE(data)

  end subroutine saveSpongeStrength

  subroutine resetProbes(this)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use ProbePatch_mod, only : t_ProbePatch
    use ActuatorPatch_mod, only : t_ActuatorPatch

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: this

    ! <<< Local variables >>>
    integer :: i
    class(t_Patch), pointer :: patch => null()

    if (allocated(this%patchFactories)) then
       do i = 1, size(this%patchFactories)
          call this%patchFactories(i)%connect(patch)
          if (.not. associated(patch)) cycle
          if (patch%comm == MPI_COMM_NULL) cycle
          select type (patch)
          class is (t_ActuatorPatch)
          class is (t_ProbePatch)
             call patch%reset()
          end select
       end do
    end if

  end subroutine resetProbes

  subroutine saveProbeData(this, mode, finish)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use ProbePatch_mod, only : t_ProbePatch
    use ActuatorPatch_mod, only : t_ActuatorPatch

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: this
    integer, intent(in) :: mode
    logical, intent(in), optional :: finish

    ! <<< Local variables >>>
    logical :: finish_
    integer :: i, j
    class(t_Patch), pointer :: patch => null()
    class(t_ProbePatch), pointer :: probePatch => null()

    finish_ = .false.
    if (present(finish)) finish_ = finish

    do i = 1, size(this%patchFactories)
       call this%patchFactories(i)%connect(patch)
       if (.not. associated(patch)) cycle
       do j = 1, size(this%states)
          if (patch%gridIndex /= this%grids(j)%index .or. patch%nPatchPoints <= 0) cycle

          nullify(probePatch)
          select type (patch)
          class is (t_ActuatorPatch)
          class is (t_ProbePatch)
             probePatch => patch
          end select

          if (.not. associated(probePatch)) cycle

          if (finish_) then
             call probePatch%saveData()
             probePatch%iProbeBuffer = 0
             cycle
          end if

          probePatch%iProbeBuffer = probePatch%iProbeBuffer + 1
          assert(probePatch%iProbeBuffer >= 1)
          assert(probePatch%iProbeBuffer <= size(probePatch%probeBuffer, 3))

          call probePatch%update(mode, this%simulationFlags, this%solverOptions,             &
               this%grids(j), this%states(j))

          if (probePatch%iProbeBuffer == size(probePatch%probeBuffer, 3)) then
             call probePatch%saveData()
             probePatch%iProbeBuffer = 0
          end if

       end do
    end do

  end subroutine saveProbeData

  subroutine readDecompositionMap(this, filename)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env

    ! <<< Internal modules >>>
    use InputHelper, only : getFreeUnit, stripComments
    use ErrorHandler, only : gracefulExit, writeAndFlush

    ! <<< Arguments >>>
    class(t_Region) :: this
    character(len = *), intent(in) :: filename

    ! <<< Local variables >>>
    integer :: i, fileUnit, proc, numProcs, lineNo, gridIndex,                               &
         numProcsInGrid(3), istat, ierror
    character(len = STRING_LENGTH) :: line, message
    character(len = 1), parameter :: commentMarker = '#'

    call MPI_Comm_rank(this%comm, proc, ierror)
    call MPI_Comm_size(this%comm, numProcs, ierror)

    write(message, "(3A)") "Reading MPI decomposition map from '", trim(filename), "'..."
    call writeAndFlush(this%comm, output_unit, message)

    ! Check if file exists.
    if (proc == 0) then
       open(unit = getFreeUnit(fileUnit), file = trim(filename), action = 'read',            &
            status = 'old', iostat = istat)
    end if
    call MPI_Bcast(istat, 1, MPI_INTEGER, 0, this%comm, ierror)
    if (istat /= 0) then
       write(message, "(2A)") trim(filename), ": File not found or permission denied!"
       call gracefulExit(this%comm, message)
    end if

    ! Only the root process reads the file.
    if (proc == 0) then
       i = 0; lineNo = 0; istat = 0
       do

          read(fileUnit, '(A)', iostat = istat) line
          if (istat < 0) then
             istat = 0
             exit
          end if
          lineNo = lineNo + 1 !... need line number for reporting errors.
          call stripComments(line, commentMarker)
          if (len_trim(line) == 0) cycle
          i = i + 1

          read(line, *, iostat = istat) gridIndex,                                           &
               numProcsInGrid(1), numProcsInGrid(2), numProcsInGrid(3)

          if (istat /= 0) then
             write(message, "(2A,I0.0,A)") trim(filename),                                   &
                  ": Failed to parse input on line ", lineNo, "!"
             exit
          end if

          if (gridIndex /= i .or. gridIndex > size(this%processDistributions, 2) .or.        &
               any(numProcsInGrid < 0) .or. product(numProcsInGrid) > numProcs) then
             istat = -1
             write(message, "(2A,I0.0,A)") trim(filename),                                   &
                  ": Invalid process distribution on line ", lineNo, "!"
             exit
          end if
          if (size(this%processDistributions, 1) < 3) then
             if (any(numProcsInGrid(size(this%processDistributions,1)+1:) /= 1)) then
                istat = -1
                write(message, "(2A,I0.0,A)") trim(filename),                                &
                     ": Invalid process distribution on line ", lineNo, "!"
                exit
             end if
          end if

          this%processDistributions(:,i) = numProcsInGrid(1:size(this%processDistributions,1))

       end do
       close(fileUnit)
    end if

    ! Check if an error occurred.
    call MPI_Bcast(istat, 1, MPI_INTEGER, 0,                                                 &
         this%comm, ierror) !... broadcast istat to collectively return on error.
    if (istat /= 0) then
       call MPI_Bcast(message, len(message), MPI_CHARACTER, 0, this%comm, ierror)
       call gracefulExit(this%comm, message)
    end if

    ! Broadcast process distribution.
    call MPI_Bcast(this%processDistributions, size(this%processDistributions), MPI_INTEGER,  &
         0, this%comm, ierror)

    ! Validate process distribution.
    if (sum(product(this%processDistributions, dim = 1)) /= numProcs) then
       write(message, '(A,2(A,I0.0),A)') trim(filename),                                     &
            ": Invalid process distribution: expected a total of ", numProcs,                &
            " processes, got ", sum(product(this%processDistributions, dim = 1)),            &
            " processes!"
       call gracefulExit(this%comm, message)
    end if

  end subroutine readDecompositionMap

  subroutine distributeGrids(this, verbose)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env

    ! <<< Internal modules >>>
    use MPIHelper, only : splitCommunicatorMultigrid
    use ErrorHandler, only : writeAndFlush

    ! <<< Arguments >>>
    class(t_Region) :: this
    logical, intent(in) :: verbose

    ! <<< Local variables >>>
    integer :: numProcs, ierror
    character(len = STRING_LENGTH) :: message
    integer, allocatable :: numProcsInGrid(:)

    ! Find the size of the communicator.
    call MPI_Comm_size(this%comm, numProcs, ierror)

    if (verbose) then
       write(message, "(2(A,I0.0),A)") "Distributing ", size(this%globalGridSizes, 2),       &
            " grid(s) across ", numProcs, " process(es)..."
       call writeAndFlush(this%comm, output_unit, message)
    end if

    ! Find the number of processes to be assigned to each grid: `numProcsInGrid(i)` is the
    ! number of processes assigned to grid `i`.
    allocate(numProcsInGrid(size(this%globalGridSizes, 2)))
    numProcsInGrid = product(this%processDistributions, dim = 1)
    if (any(numProcsInGrid <= 0)) numProcsInGrid = 0

    ! `gridCommunicators(i)`: Communicator for grid `i`. If the current process does not
    ! contain a part of grid `i`, then `gridCommunicators(i)` is `MPI_COMM_NULL`.
    SAFE_DEALLOCATE(this%gridCommunicators)
    allocate(this%gridCommunicators(size(this%globalGridSizes, 2)), source = MPI_COMM_NULL)

    ! Split the region communicator into grid-level communicators:
    if (numProcs > size(this%globalGridSizes, 2) .and.                                       &
         this%simulationFlags%manualDomainDecomp) then
       call splitCommunicatorMultigrid(this%comm, this%globalGridSizes,                      &
            this%gridCommunicators, numProcsInGrid) !... manual process distribution.
    else
       call splitCommunicatorMultigrid(this%comm, this%globalGridSizes,                      &
            this%gridCommunicators) !... automatic process distribution.
    end if

    SAFE_DEALLOCATE(numProcsInGrid)

  end subroutine distributeGrids

  subroutine readBoundaryConditions(this, filename)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use Patch_factory, only : t_PatchFactory

    ! <<< Internal modules >>>
    use InputHelper, only : getFreeUnit, stripComments
    use ErrorHandler, only : gracefulExit, writeAndFlush

    ! <<< Arguments >>>
    class(t_Region) :: this
    character(len = *), intent(in) :: filename

    ! <<< Local variables >>>
    integer :: i, fileUnit, proc, nPatches, lineNo, istat, ierror
    character(len = STRING_LENGTH) :: line, message
    character(len = 1), parameter :: commentMarker = '#'
    type(t_PatchFactory) :: patchFactory
    class(t_Patch), pointer :: patch => null()
    logical :: isExtentValid

    call MPI_Comm_rank(this%comm, proc, ierror)

    ! Check if file exists.
    if (proc == 0) then
       open(unit = getFreeUnit(fileUnit), file = trim(filename), action = 'read',            &
            status = 'old', iostat = istat)
    end if
    call MPI_Bcast(istat, 1, MPI_INTEGER, 0, this%comm, ierror)
    if (istat /= 0) then
       write(message, "(2A)") trim(filename), ": File not found or permission denied!"
       call gracefulExit(this%comm, message)
    end if

    write(message, "(3A)") "Reading boundary conditions from '", trim(filename), "'..."
    call writeAndFlush(this%comm, output_unit, message)

    ! Only the root process reads the file.
    if (proc == 0) then
       nPatches = 0
       do !... read once to find the number of patches.
          read(fileUnit, '(A)', iostat = istat) line
          if (istat < 0) exit
          call stripComments(line, commentMarker) !... skip comments.
          if (len_trim(line) == 0) cycle !... skip empty lines.
          nPatches = nPatches + 1
       end do
       close(fileUnit)
    end if

    ! Broadcast the number of patches to all processes.
    call MPI_Bcast(nPatches, 1, MPI_INTEGER, 0, this%comm, ierror)

    if (nPatches == 0) return

    write(message, "(A,I0.0,A)") "Found ", nPatches, " boundary conditions!"
    call writeAndFlush(this%comm, output_unit, message)

    allocate(this%patchData(nPatches, 8))
    allocate(this%patchNames(nPatches))
    allocate(this%patchTypes(nPatches))

    ! Again, only the root process reads the file.
    if (proc == 0) then
       i = 0; lineNo = 0; istat = 0
       open(unit = getFreeUnit(fileUnit), file = trim(filename),                             &
            action = 'read', status = 'old')
       do !... read again to fill patch information.

          read(fileUnit, '(A)', iostat = istat) line
          if (istat < 0) then
             istat = 0
             exit
          end if
          lineNo = lineNo + 1 !... need line number for reporting errors.
          call stripComments(line, commentMarker)
          if (len_trim(line) == 0) cycle
          i = i + 1

          ! Parse patch data.
          read(line, *, iostat = istat) this%patchNames(i), this%patchTypes(i),              &
               this%patchData(i,1), this%patchData(i,2), this%patchData(i,3),                &
               this%patchData(i,4), this%patchData(i,5), this%patchData(i,6),                &
               this%patchData(i,7), this%patchData(i,8)

          if (istat /= 0) then
             write(message, "(2A,I0.0,A)") trim(filename), ":", lineNo,                      &
                  ": Failed to parse input on this line!"
             exit
          end if

          if (any(this%patchNames(:i-1) == this%patchNames(i))) then
             istat = -1
             write(message, "(2A,I0.0,3A)") trim(filename), ":", lineNo,                     &
                  ": A patch with name '", trim(this%patchNames(i)), "' already exists!"
             exit
          end if

          call patchFactory%connect(patch, this%patchTypes(i), .true.)
          if (.not. associated(patch)) then
             istat = -1
             write(message, "(2A,I0.0,5A)") trim(filename), ":", lineNo,                     &
                  ": Invalid type '", trim(this%patchTypes(i)), "' for patch '",             &
                  trim(this%patchNames(i)), "'!"
             exit
          end if
          call patchFactory%cleanup()

          if (this%patchData(i,1) <= 0 .or.                                                  &
               this%patchData(i,1) > size(this%globalGridSizes, 2)) then
             istat = -1
             write(message, '(3A,I0.0,A)') "Patch '", trim(this%patchNames(i)),              &
                  "' has an invalid grid index: ", this%patchData(i,1), "!"
             exit
          end if

          if (abs(this%patchData(i,2)) > size(this%globalGridSizes, 1)) then
             istat = -1
             write(message, '(3A)') "Patch '", trim(this%patchNames(i)),                     &
                  "' has an invalid normal direction!"
             exit
          end if

          call checkPatchExtent(this%patchData(i,3:8),                                       &
               this%globalGridSizes(:,this%patchData(i,1)), isExtentValid)
          if (.not. isExtentValid) then
             istat = -1
             write(message, '(3A)') "Patch '", trim(this%patchNames(i)),                     &
                  "' has an invalid extent!"
             exit
          end if

       end do
       close(fileUnit)
    end if !... proc == 0

    ! Check if an error occurred.
    call MPI_Bcast(istat, 1, MPI_INTEGER, 0,                                                 &
         this%comm, ierror) !... broadcast istat to collectively return on error.
    if (istat /= 0) then
       call MPI_Bcast(message, len(message),                                                 &
            MPI_CHARACTER, 0, this%comm, ierror)
       call gracefulExit(this%comm, message)
    end if

    ! Broadcast and unpack data.
    call MPI_Bcast(this%patchData, size(this%patchData), MPI_INTEGER, 0, this%comm, ierror)
    do i = 1, nPatches
       call MPI_Bcast(this%patchNames(i), len(this%patchNames(i)),                           &
            MPI_CHARACTER, 0, this%comm, ierror)
       call MPI_Bcast(this%patchTypes(i), len(this%patchTypes(i)),                           &
            MPI_CHARACTER, 0, this%comm, ierror)
    end do

  end subroutine readBoundaryConditions

  subroutine checkPatchExtent(extent, globalGridSize, isExtentValid)

    implicit none

    ! <<< Arguments >>>
    integer, intent(inout) :: extent(:)
    integer, intent(in) :: globalGridSize(:)
    logical, intent(out) :: isExtentValid

    ! <<< Local variabels >>>
    integer :: i, nDimensions

    nDimensions = size(globalGridSize)
    assert_key(nDimensions, (1, 2, 3))

    isExtentValid = .true.

    if (any(extent == 0)) then
       isExtentValid = .false.
       return
    end if

    ! Negative values indicate counting backwards from the end.
    do i = 1, nDimensions
       if (extent(1+2*(i-1)) < 0)                                                            &
            extent(1+2*(i-1)) = extent(1+2*(i-1)) + globalGridSize(i) + 1
       if (extent(2+2*(i-1)) < 0)                                                            &
            extent(2+2*(i-1)) = extent(2+2*(i-1)) + globalGridSize(i) + 1
    end do

    ! Check that extent describes a part of the grid.
    do i = 1, nDimensions
       isExtentValid = isExtentValid .and. (extent(2+2*(i-1)) >= extent(1+2*(i-1)))
       isExtentValid = isExtentValid .and.                                                   &
            (extent(2+2*(i-1)) - extent(1+2*(i-1)) + 1 <= globalGridSize(i))
    end do
    do while (i <= 3) !... reset for direction > number of dimensions.
       extent(1+2*(i-1)) = 1
       extent(2+2*(i-1)) = 1
       i = i + 1
    end do

  end subroutine checkPatchExtent

  subroutine distributePatches(this)

    ! <<< External modules >>>
    use MPI

    ! <<< Arguments >>>
    class(t_Region) :: this

    ! <<< Local variables >>>
    integer :: i, j, nPatches, gridOffset(3), gridLocalSize(3), gridIndex, color, comm,      &
         rankInGridCommunicator, rankInRegionCommunicator, rankInPatchCommunicator, ierror

    if (.not. allocated(this%patchData)) return

    assert(allocated(this%patchNames))
    assert(allocated(this%patchTypes))
    assert(size(this%patchNames) == size(this%patchData, 1))
    assert(size(this%patchTypes) == size(this%patchData, 1))

    nPatches = size(this%patchData, 1)

    allocate(this%patchCommunicators(nPatches), source = MPI_COMM_NULL)
    allocate(this%patchMasterRanks(nPatches), source = -1)

    call MPI_Comm_rank(this%comm, rankInRegionCommunicator, ierror)

    do i = 1, size(this%grids)

       gridIndex = this%grids(i)%index
       gridOffset = this%grids(i)%offset
       gridLocalSize = this%grids(i)%localSize

       call MPI_Comm_rank(this%grids(i)%comm, rankInGridCommunicator, ierror)

       do j = 1, nPatches

          if (this%patchData(j,1) /= gridIndex) cycle

          color = 1
          if (this%patchData(j,4) < gridOffset(1) + 1 .or.                                   &
               this%patchData(j,3) > gridOffset(1) + gridLocalSize(1) .or.                   &
               this%patchData(j,6) < gridOffset(2) + 1 .or.                                  &
               this%patchData(j,5) > gridOffset(2) + gridLocalSize(2) .or.                   &
               this%patchData(j,8) < gridOffset(3) + 1 .or.                                  &
               this%patchData(j,7) > gridOffset(3) + gridLocalSize(3)) then
             color = MPI_UNDEFINED
          end if
          call MPI_Comm_split(this%grids(i)%comm, color, rankInGridCommunicator, comm, ierror)
          this%patchCommunicators(j) = comm

          if (comm /= MPI_COMM_NULL) then
             call MPI_Comm_rank(comm, rankInPatchCommunicator, ierror)
             if (rankInPatchCommunicator == 0) then
                assert(this%patchMasterRanks(j) == -1)
                this%patchMasterRanks(j) = rankInRegionCommunicator
             end if
          end if

          call MPI_Barrier(this%grids(i)%comm, ierror)

       end do

    end do

    call MPI_Allreduce(MPI_IN_PLACE, this%patchMasterRanks, size(this%patchMasterRanks),     &
         MPI_INTEGER, MPI_MAX, this%comm, ierror)

  end subroutine distributePatches

  subroutine normalizeControlMollifier(this)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use ActuatorPatch_mod, only : t_ActuatorPatch

    ! <<< Internal modules >>>
    use ErrorHandler, only : gracefulExit, issueWarning

    ! <<< Arguments >>>
    class(t_Region) :: this

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: mollifierNorm
    integer :: i, j, ierror
    logical :: hasNegativeMollifier
    character(len = STRING_LENGTH) :: message
    class(t_Patch), pointer :: patch => null()
    real(wp), allocatable :: ones(:)

    assert(allocated(this%grids))

    mollifierNorm = 0.0_wp

    do i = 1, size(this%grids)

       assert(allocated(this%grids(i)%controlMollifier))

       hasNegativeMollifier = any(this%grids(i)%controlMollifier(:,1) < 0.0_wp)
       call MPI_Allreduce(MPI_IN_PLACE, hasNegativeMollifier, 1,                             &
            MPI_LOGICAL, MPI_LOR, this%grids(i)%comm, ierror)
       if (hasNegativeMollifier) then
          write(message, '(A,I0.0,A)') "Control mollifying support function on grid ",       &
               this%grids(i)%index, " is not non-negative everywhere!"
          call gracefulExit(this%grids(i)%comm, message)
       end if

       allocate(ones(this%grids(i)%nGridPoints), source = 1.0_wp)

       if (allocated(this%patchFactories)) then
          do j = 1, size(this%patchFactories)
             call this%patchFactories(j)%connect(patch)
             if (.not. associated(patch)) cycle
             if (patch%gridIndex /= this%grids(i)%index) cycle
             select type (patch)
             class is (t_ActuatorPatch)
                mollifierNorm = mollifierNorm + patch%computeInnerProduct(this%grids(i),     &
                     this%grids(i)%controlMollifier(:,1), ones)
             end select
          end do
       end if

       SAFE_DEALLOCATE(ones)

    end do

    if (this%commGridMasters /= MPI_COMM_NULL)                                               &
         call MPI_Allreduce(MPI_IN_PLACE, mollifierNorm, 1, SCALAR_TYPE_MPI,                 &
         MPI_SUM, this%commGridMasters, ierror)

    do i = 1, size(this%grids)
       call MPI_Bcast(mollifierNorm, 1, SCALAR_TYPE_MPI, 0, this%grids(i)%comm, ierror)
    end do

    if (mollifierNorm <= 0.0_wp) then
       call issueWarning(this%comm,                                                          &
            "Control mollifying support is trivial! Is an actuator patch present?")
       do i = 1, size(this%grids)
          this%grids(i)%controlMollifier = 0.0_wp
       end do
    else
       do i = 1, size(this%grids)
          this%grids(i)%controlMollifier = this%grids(i)%controlMollifier / mollifierNorm
       end do
    end if

  end subroutine normalizeControlMollifier

  subroutine normalizeTargetMollifier(this)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use CostTargetPatch_mod, only : t_CostTargetPatch

    ! <<< Internal modules >>>
    use ErrorHandler, only : gracefulExit, issueWarning

    ! <<< Arguments >>>
    class(t_Region) :: this

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: mollifierNorm
    integer :: i, j, ierror
    logical :: hasNegativeMollifier
    character(len = STRING_LENGTH) :: message
    class(t_Patch), pointer :: patch => null()
    real(wp), allocatable :: ones(:)

    assert(allocated(this%grids))

    mollifierNorm = 0.0_wp

    do i = 1, size(this%grids)

       assert(allocated(this%grids(i)%targetMollifier))

       hasNegativeMollifier = any(this%grids(i)%targetMollifier(:,1) < 0.0_wp)
       call MPI_Allreduce(MPI_IN_PLACE, hasNegativeMollifier, 1,                             &
            MPI_LOGICAL, MPI_LOR, this%grids(i)%comm, ierror)
       if (hasNegativeMollifier) then
          write(message, '(A,I0.0,A)') "Target mollifying support function on grid ",       &
               this%grids(i)%index, " is not non-negative everywhere!"
          call gracefulExit(this%grids(i)%comm, message)
       end if

       allocate(ones(this%grids(i)%nGridPoints), source = 1.0_wp)

       if (allocated(this%patchFactories)) then
          do j = 1, size(this%patchFactories)
             call this%patchFactories(j)%connect(patch)
             if (.not. associated(patch)) cycle
             if (patch%gridIndex /= this%grids(i)%index) cycle
             select type (patch)
             class is (t_CostTargetPatch)
                mollifierNorm = mollifierNorm + patch%computeInnerProduct(this%grids(i),     &
                     this%grids(i)%targetMollifier(:,1), ones)
             end select
          end do
       end if

       SAFE_DEALLOCATE(ones)

    end do

    if (this%commGridMasters /= MPI_COMM_NULL)                                               &
         call MPI_Allreduce(MPI_IN_PLACE, mollifierNorm, 1, SCALAR_TYPE_MPI,                 &
         MPI_SUM, this%commGridMasters, ierror)

    do i = 1, size(this%grids)
       call MPI_Bcast(mollifierNorm, 1, SCALAR_TYPE_MPI, 0, this%grids(i)%comm, ierror)
    end do

    if (mollifierNorm <= 0.0_wp) then
       call issueWarning(this%comm,                                                          &
            "Target mollifying support is trivial! Is a cost target patch present?")
       do i = 1, size(this%grids)
          this%grids(i)%targetMollifier = 0.0_wp
       end do
    else
       do i = 1, size(this%grids)
          this%grids(i)%targetMollifier = this%grids(i)%targetMollifier / mollifierNorm
       end do
    end if

  end subroutine normalizeTargetMollifier

  subroutine computeSpongeStrengths(this)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use SpongePatch_mod, only : t_SpongePatch

    ! <<< Internal modules >>>
    use MPIHelper, only : gatherAlongDirection

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: this

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, direction, nDimensions
    class(t_Patch), pointer :: patch => null()
    class(t_SpongePatch), pointer :: spongePatch => null()
    logical :: spongesExistAlongDirection
    real(wp), allocatable :: coordinateDerivatives(:,:), arcLength(:,:),                     &
         globalArcLengthsAlongDirection(:,:)

    nDimensions = size(this%globalGridSizes, 1)
    assert_key(nDimensions, (1, 2, 3))

    if (.not. allocated(this%patchFactories)) return

    do direction =  1, nDimensions

       do i = 1, size(this%grids)

          spongesExistAlongDirection = .false.

          ! Check if there are sponge patches along direction `direction`.
          do j = 1, size(this%patchFactories)
             call this%patchFactories(j)%connect(patch)
             if (.not. associated(patch)) cycle
             if (patch%gridIndex /= this%grids(i)%index .or. &
                  abs(patch%normalDirection) /= direction) cycle
             select type (patch)
             class is (t_SpongePatch)
                spongesExistAlongDirection = .true.
             end select
          end do

          if (.not. spongesExistAlongDirection) cycle

          allocate(arcLength(this%grids(i)%nGridPoints, 1))
          allocate(coordinateDerivatives(this%grids(i)%nGridPoints, nDimensions))
          allocate(globalArcLengthsAlongDirection(this%grids(i)%nGridPoints /                &
               this%grids(i)%localSize(direction) * this%grids(i)%globalSize(direction), 1))

          ! Compute local arc length.
          call this%grids(i)%computeCoordinateDerivatives(direction, coordinateDerivatives)
          arcLength(:,1) = sqrt(sum(coordinateDerivatives ** 2, dim = 2))

          ! Gather arc length along direction `direction`.
          call gatherAlongDirection(this%grids(i)%comm, arcLength, this%grids(i)%localSize,  &
               direction, this%grids(i)%offset(direction), globalArcLengthsAlongDirection)

          do j = 1, size(this%patchFactories)
             call this%patchFactories(j)%connect(patch)
             if (.not. associated(patch)) cycle
             if (patch%gridIndex /= this%grids(i)%index .or. &
                  abs(patch%normalDirection) /= direction) cycle

             nullify(spongePatch)
             select type (patch)
             class is (t_SpongePatch)
                spongePatch => patch
             end select

             if (associated(spongePatch)) call spongePatch%computeStrength(this%grids(i),    &
                  globalArcLengthsAlongDirection(:,1))

          end do !... j = 1, size(this%patchFactories)

          SAFE_DEALLOCATE(arcLength)
          SAFE_DEALLOCATE(coordinateDerivatives)
          SAFE_DEALLOCATE(globalArcLengthsAlongDirection)

       end do !... i = 1, size(this%grids)

    end do !... direction = 1, nDimensions

  end subroutine computeSpongeStrengths

  subroutine readPatchInterfaceInformation(this)

    ! <<< Internal modules >>>
    use InputHelper, only : getOption
    use ErrorHandler, only : gracefulExit

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: this

    ! <<< Local variables >>>
    integer :: i, j, k, l, nDimensions
    character(len = STRING_LENGTH) :: key, str, message

    if (.not. allocated(this%patchData)) return

    SAFE_DEALLOCATE(this%patchInterfaces)
    SAFE_DEALLOCATE(this%interfaceIndexReorderings)

    allocate(this%patchInterfaces(size(this%patchData, 1)), source = 0)
    allocate(this%interfaceIndexReorderings(3, size(this%patchData, 1)), source = 0)

    nDimensions = size(this%globalGridSizes, 1)

    ! Read interface information:

    do i = 1, size(this%patchNames)

       write(key, '(A)') "patches/" // trim(this%patchNames(i)) // "/conforms_with"
       str = getOption(key, "")

       if (len_trim(str) > 0) then

          do j = 1, size(this%patchNames)
             if (trim(str) == trim(this%patchNames(j))) this%patchInterfaces(i) = j
          end do
          if (this%patchInterfaces(i) == 0) then
             write(message, '(5A)') "Invalid interface specification for patch '",           &
                  trim(this%patchNames(i)), "': no patch found matching the name '",       &
                  trim(str), "'!"
             call gracefulExit(this%comm, message)
          end if

          do j = 1, 3
             this%interfaceIndexReorderings(j,i) = j
             if (j <= nDimensions) then
                write(key, '(A,I1)') "patches/" // trim(this%patchNames(i)) //             &
                     "/interface_index", j
                this%interfaceIndexReorderings(j,i) = getOption(trim(key), j)
                if (j == 3 .and. this%interfaceIndexReorderings(j,i) /= 3) then
                   write(message, '(3A)') "Interface index reordering for patch '",          &
                        trim(this%patchNames(i)), "' is currently not supported!"
                   call gracefulExit(this%comm, message)
                end if
             end if
          end do

          if (.not. all(this%interfaceIndexReorderings(1:nDimensions,i) /= 0 .and.         &
               abs(this%interfaceIndexReorderings(1:nDimensions,i)) <= nDimensions)) then
             write(message, '(3A)') "Invalid interface index reordering for patch '",        &
                  trim(this%patchNames(i)), "'!"
             call gracefulExit(this%comm, message)
          end if

       end if

    end do

    ! Commutativity of interfaces:

    do i = 1, size(this%patchData, 1)
       if (this%patchInterfaces(i) == 0) cycle
       j = this%patchInterfaces(i)

       if (this%patchInterfaces(j) == 0) then

          this%patchInterfaces(j) = i

          do l = 1, 3
             do k = 1, 3
                if (abs(this%interfaceIndexReorderings(k,i)) == l) then
                   this%interfaceIndexReorderings(l,j) =                                   &
                        sign(k, this%interfaceIndexReorderings(k,i))
                   exit
                end if
             end do
          end do

       else if (this%patchInterfaces(j) /= i) then

          write(message, '(3A)') "Invalid interface specification for patch '",              &
               trim(this%patchNames(j)), "': violates commutativity!"
          call gracefulExit(this%comm, message)

       end if

    end do

  end subroutine readPatchInterfaceInformation

  subroutine exchangeInterfaceData(this)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: this

    ! <<< Local variables >>>
    integer :: i, j, iRequest, mpiTag, procRank, ierror
    integer, allocatable :: mpiSendRequests(:)
    class(t_Patch), pointer :: patch => null()
    class(t_BlockInterfacePatch), pointer :: interfacePatch => null()

    if (.not. allocated(this%patchInterfaces)) return

    call MPI_Comm_rank(this%comm, procRank, ierror)

    if (any(this%patchMasterRanks == procRank))                                              &
         allocate(mpiSendRequests(count(this%patchInterfaces > 0 .and.                       &
         this%patchMasterRanks == procRank)), source = MPI_REQUEST_NULL)

    iRequest = 0

    do i = 1, size(this%patchInterfaces)
       if (procRank == this%patchMasterRanks(i) .and. this%patchInterfaces(i) > 0) then

          assert(allocated(this%patchFactories))

          mpiTag = i + size(this%patchInterfaces) * (this%patchInterfaces(i) - 1)

          nullify(interfacePatch)

          do j = 1, size(this%patchFactories)
             call this%patchFactories(j)%connect(patch)
             if (.not. associated(patch)) cycle
             if (patch%index /= i) cycle
             select type (patch)
             class is (t_BlockInterfacePatch)
                interfacePatch => patch
                exit
             end select
          end do

          assert(associated(interfacePatch))
          assert(allocated(interfacePatch%sendBuffer))

          iRequest = iRequest + 1
          call MPI_Isend(interfacePatch%sendBuffer,                                          &
               size(interfacePatch%sendBuffer), SCALAR_TYPE_MPI,                             &
               this%patchMasterRanks(this%patchInterfaces(i)), mpiTag, this%comm,            &
               mpiSendRequests(iRequest), ierror)

       end if
    end do

    do i = 1, size(this%patchInterfaces)
       if (procRank == this%patchMasterRanks(i) .and. this%patchInterfaces(i) > 0) then

          assert(allocated(this%patchFactories))

          mpiTag = this%patchInterfaces(i) + size(this%patchInterfaces) * (i - 1)

          nullify(interfacePatch)

          do j = 1, size(this%patchFactories)
             call this%patchFactories(j)%connect(patch)
             if (.not. associated(patch)) cycle
             if (patch%index /= i) cycle
             select type (patch)
             class is (t_BlockInterfacePatch)
                interfacePatch => patch
                exit
             end select
          end do

          assert(associated(interfacePatch))
          assert(allocated(interfacePatch%receiveBuffer))

          call MPI_Recv(interfacePatch%receiveBuffer,                                        &
               size(interfacePatch%receiveBuffer), SCALAR_TYPE_MPI,                          &
               this%patchMasterRanks(this%patchInterfaces(i)),                               &
               mpiTag, this%comm, MPI_STATUS_IGNORE, ierror)

       end if
    end do

    if (allocated(mpiSendRequests))                                                          &
         call MPI_Waitall(size(mpiSendRequests), mpiSendRequests, MPI_STATUSES_IGNORE, ierror)

    do i = 1, size(this%patchInterfaces)
       if (procRank == this%patchMasterRanks(i) .and. this%patchInterfaces(i) > 0) then

          assert(allocated(this%patchFactories))

          nullify(interfacePatch)

          do j = 1, size(this%patchFactories)
             call this%patchFactories(j)%connect(patch)
             if (.not. associated(patch)) cycle
             if (patch%index /= i) cycle
             select type (patch)
             class is (t_BlockInterfacePatch)
                interfacePatch => patch
                exit
             end select
          end do

          assert(associated(interfacePatch))

          call interfacePatch%reshapeReceiveBuffer(this%interfaceIndexReorderings(:,i))

       end if
       call MPI_Barrier(this%comm, ierror)
    end do

    SAFE_DEALLOCATE(mpiSendRequests)

  end subroutine exchangeInterfaceData

end module Region_mod
