#include "config.h"

module RegionImpl

  implicit none
  public

contains

  subroutine readDecompositionMap(this, filename)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Internal modules >>>
    use InputHelper, only : stripComments
    use ErrorHandler, only : gracefulExit, writeAndFlush

    ! <<< Arguments >>>
    class(t_Region) :: this
    character(len = *), intent(in) :: filename

    ! <<< Local variables >>>
    integer :: i, fileUnit, proc, nProcs, lineNo, gridIndex, numProcsInGrid(3), istat, ierror
    character(len = STRING_LENGTH) :: line, message
    character(len = 1), parameter :: commentMarker = '#'

    call MPI_Comm_rank(this%comm, proc, ierror)
    call MPI_Comm_size(this%comm, nProcs, ierror)

    write(message, "(3A)") "Reading MPI decomposition map from '", trim(filename), "'..."
    call writeAndFlush(this%comm, output_unit, message)

    ! Check if file exists.
    if (proc == 0) then
       open(newunit = fileUnit, file = trim(filename), action = 'read', status = 'old',      &
            iostat = istat)
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
               any(numProcsInGrid < 0) .or. product(numProcsInGrid) > nProcs) then
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
    if (sum(product(this%processDistributions, dim = 1)) /= nProcs) then
       write(message, '(A,2(A,I0.0),A)') trim(filename),                                     &
            ": Invalid process distribution: expected a total of ", nProcs,                  &
            " processes, got ", sum(product(this%processDistributions, dim = 1)),            &
            " processes!"
       call gracefulExit(this%comm, message)
    end if

  end subroutine readDecompositionMap

  subroutine distributeGrids(this, verbose)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Internal modules >>>
    use MPIHelper, only : splitCommunicatorMultigrid
    use ErrorHandler, only : writeAndFlush

    ! <<< Arguments >>>
    class(t_Region) :: this
    logical, intent(in) :: verbose

    ! <<< Local variables >>>
    integer :: nProcs, ierror
    character(len = STRING_LENGTH) :: message
    integer, allocatable :: numProcsInGrid(:)

    ! Find the size of the communicator.
    call MPI_Comm_size(this%comm, nProcs, ierror)

    if (verbose) then
       write(message, "(2(A,I0.0),A)") "Distributing ", size(this%globalGridSizes, 2),       &
            " grid(s) across ", nProcs, " process(es)..."
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
    if (nProcs > size(this%globalGridSizes, 2) .and.                                         &
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
    use Region_mod, only : t_Region
    use Patch_factory, only : t_PatchFactory

    ! <<< Internal modules >>>
    use InputHelper, only : stripComments
    use ErrorHandler, only : gracefulExit, writeAndFlush

    ! <<< Arguments >>>
    class(t_Region) :: this
    character(len = *), intent(in) :: filename

    ! <<< Local variables >>>
    integer :: i, fileUnit, proc, nPatches, lineNo, istat, ierror
    character(len = STRING_LENGTH) :: line, message
    character(len = 1), parameter :: commentMarker = '#'
    type(t_PatchFactory) :: patchFactory
    class(t_Patch), pointer :: dummyPatch => null()
    integer, allocatable :: tempBuffer(:,:)

    call MPI_Comm_rank(this%comm, proc, ierror)

    ! Check if file exists.
    if (proc == 0) then
       open(newunit = fileUnit, file = trim(filename), action = 'read', status = 'old',      &
            iostat = istat)
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

    if (nPatches > 0) then
       write(message, "(A,I0.0,A)") "Found ", nPatches, " boundary conditions!"
       call writeAndFlush(this%comm, output_unit, message)
    end if

    ! Allocate memory to hold patch information.
    SAFE_DEALLOCATE(this%patchData)
    if (nPatches == 0) return
    allocate(this%patchData(nPatches), stat = istat)
    if (istat /= 0) then
       write(message, "(A)") "Insufficient memory: &
            &Could not allocate storage for patch information!"
       call gracefulExit(this%comm, message)
    end if

    ! Allocate tempBuffer for use in MPI communication.
    allocate(tempBuffer(nPatches, 8))

    ! Again, only the root process reads the file.
    if (proc == 0) then
       i = 0; lineNo = 0; istat = 0
       open(newunit = fileUnit, file = trim(filename), action = 'read', status = 'old')
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
          read(line, *, iostat = istat) this%patchData(i)%name, this%patchData(i)%patchType, &
               this%patchData(i)%gridIndex, this%patchData(i)%normalDirection,               &
               this%patchData(i)%iMin, this%patchData(i)%iMax,                               &
               this%patchData(i)%jMin, this%patchData(i)%jMax,                               &
               this%patchData(i)%kMin, this%patchData(i)%kMax

          if (istat /= 0) then
             write(message, "(2A,I0.0,A)") trim(filename), ":", lineNo,                      &
                  ": Failed to parse input on this line!"
             exit
          end if

          if (any(this%patchData(:i-1)%name == this%patchData(i)%name)) then
             istat = -1
             write(message, "(2A,I0.0,3A)") trim(filename), ":", lineNo,                     &
                  ": A patch with name '", trim(this%patchData(i)%name), "' already exists!"
             exit
          end if

          call patchFactory%connect(dummyPatch, this%patchData(i)%patchType, .true.)
          if (.not. associated(dummyPatch)) then
             istat = -1
             write(message, "(2A,I0.0,3A)") trim(filename), ":", lineNo,                     &
                  ": Invalid type for patch '", trim(this%patchData(i)%name), "'!"
             exit
          end if
          call patchFactory%cleanup()

          ! Pack patch data into tempBuffer for broadcasting.
          tempBuffer(i,:) = (/ this%patchData(i)%gridIndex,                                  &
               this%patchData(i)%normalDirection,                                            &
               this%patchData(i)%iMin, this%patchData(i)%iMax,                               &
               this%patchData(i)%jMin, this%patchData(i)%jMax,                               &
               this%patchData(i)%kMin, this%patchData(i)%kMax /)

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
    call MPI_Bcast(tempBuffer, size(tempBuffer), MPI_INTEGER, 0, this%comm, ierror)
    do i = 1, nPatches
       this%patchData(i)%gridIndex = tempBuffer(i,1)
       this%patchData(i)%normalDirection = tempBuffer(i,2)
       this%patchData(i)%iMin = tempBuffer(i,3)
       this%patchData(i)%iMax = tempBuffer(i,4)
       this%patchData(i)%jMin = tempBuffer(i,5)
       this%patchData(i)%jMax = tempBuffer(i,6)
       this%patchData(i)%kMin = tempBuffer(i,7)
       this%patchData(i)%kMax = tempBuffer(i,8)
       call MPI_Bcast(this%patchData(i)%name, len(this%patchData(i)%name),                   &
            MPI_CHARACTER, 0, this%comm, ierror)
       call MPI_Bcast(this%patchData(i)%patchType, len(this%patchData(i)%patchType),         &
            MPI_CHARACTER, 0, this%comm, ierror)
    end do

    SAFE_DEALLOCATE(tempBuffer)

  end subroutine readBoundaryConditions

  subroutine validatePatches(this)

    ! <<< External modules >>>
    use, intrinsic :: iso_fortran_env

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Internal modules >>>
    use ErrorHandler, only : issueWarning, gracefulExit, writeAndFlush

    ! <<< Arguments >>>
    class(t_Region) :: this

    ! <<< Local variables >>>
    integer :: i, j, errorCode
    character(len = STRING_LENGTH) :: message

    if (.not. allocated(this%patchData)) return

    write(message, "(A,I0.0,A)") "Validating boundary conditions..."
    call writeAndFlush(this%comm, output_unit, message)

    do i = 1, size(this%patchData)
       call this%patchData(i)%validate(this%globalGridSizes,                                 &
            this%simulationFlags, this%solverOptions, errorCode, message)
       if (errorCode == 1) then
          call issueWarning(this%comm, message)
       else if (errorCode == 2) then
          call gracefulExit(this%comm, message)
       end if
    end do

    ! Validate interface connections.
    do i = 1, size(this%patchData)
       if (this%patchInterfaces(i) == 0) cycle
       j = this%patchInterfaces(i)
       call this%patchData(i)%validateInterface(this%globalGridSizes, this%simulationFlags,  &
            this%solverOptions, this%interfaceIndexReorderings(:,i),                         &
            this%patchData(j), errorCode, message)
       if (errorCode == 1) then
          call issueWarning(this%comm, message)
       else if (errorCode == 2) then
          call gracefulExit(this%comm, message)
       end if
    end do

  end subroutine validatePatches

  subroutine distributePatches(this)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_mod, only : t_Region
    use PatchDescriptor_mod, only : t_PatchDescriptor

    ! <<< Arguments >>>
    class(t_Region) :: this

    ! <<< Local variables >>>
    integer :: i, j, gridOffset(3), gridLocalSize(3), gridIndex, color,                      &
         comm, procRank, rankInPatchCommunicator, ierror
    type(t_PatchDescriptor) :: p

    if (.not. allocated(this%patchData)) return

    SAFE_DEALLOCATE(this%patchCommunicators)
    SAFE_DEALLOCATE(this%patchMasterRanks)

    allocate(this%patchCommunicators(size(this%patchData)), source = MPI_COMM_NULL)
    allocate(this%patchMasterRanks(size(this%patchData)), source = -1)

    do i = 1, size(this%grids)

       gridIndex = this%grids(i)%index
       gridOffset = this%grids(i)%offset
       gridLocalSize = this%grids(i)%localSize

       call MPI_Comm_rank(this%grids(i)%comm, procRank, ierror)

       do j = 1, size(this%patchData)

          p = this%patchData(j)
          color = 1
          if (p%gridIndex /= gridIndex .or.                                                  &
               p%iMax < gridOffset(1) + 1 .or.                                               &
               p%iMin > gridOffset(1) + gridLocalSize(1) .or.                                &
               p%jMax < gridOffset(2) + 1 .or.                                               &
               p%jMin > gridOffset(2) + gridLocalSize(2) .or.                                &
               p%kMax < gridOffset(3) + 1 .or.                                               &
               p%kMin > gridOffset(3) + gridLocalSize(3)) then
             color = MPI_UNDEFINED
          end if
          call MPI_Comm_split(this%grids(i)%comm, color, procRank, comm, ierror)
          this%patchCommunicators(j) = comm

          if (comm /= MPI_COMM_NULL) then
             call MPI_Comm_rank(comm, rankInPatchCommunicator, ierror)
             if (rankInPatchCommunicator == 0) then
                assert(this%patchMasterRanks(j) == -1)
                this%patchMasterRanks(j) = procRank
             end if
          end if

          call MPI_Barrier(this%grids(i)%comm, ierror)

       end do

    end do

    call MPI_Allreduce(MPI_IN_PLACE, this%patchMasterRanks, size(this%patchMasterRanks),     &
         MPI_INTEGER, MPI_MAX, this%comm, ierror)

  end subroutine distributePatches

end module RegionImpl

subroutine setupRegion(this, comm, globalGridSizes, boundaryConditionFilename,               &
     simulationFlags, solverOptions, verbose)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Private members >>>
  use RegionImpl

  ! <<< Public members >>>
  use MPITimingsHelper, only : startTiming, endTiming

  ! <<< Internal modules >>>
  use InputHelper, only : getRequiredOption
  use InterfaceHelper, only : readPatchInterfaceInformation

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: this
  integer, intent(in) :: comm, globalGridSizes(:,:)
  character(len = *), intent(in), optional :: boundaryConditionFilename
  type(t_SimulationFlags), intent(in), optional :: simulationFlags
  type(t_SolverOptions), intent(in), optional :: solverOptions
  logical, intent(in), optional :: verbose

  ! <<< Local variables >>>
  integer :: i, j, k, nPatches, color, procRank, nProcs, ierror
  logical :: verbose_
  character(len = STRING_LENGTH) :: decompositionMapFilename
  type(t_PatchDescriptor) :: p
  integer, allocatable :: patchIndices(:)
  class(t_Patch), pointer :: patch => null()

  call startTiming("setupRegion")

  ! Clean slate.
  call this%cleanup()
  this%comm = comm
  call MPI_Comm_size(this%comm, nProcs, ierror)

  verbose_ = .true.
  if (present(verbose)) verbose_ = verbose

  SAFE_DEALLOCATE(this%globalGridSizes)
  SAFE_DEALLOCATE(this%processDistributions)

  allocate(this%globalGridSizes(size(globalGridSizes, 1), size(globalGridSizes, 2)),         &
       source = globalGridSizes)
  allocate(this%processDistributions(size(this%globalGridSizes, 1),                          &
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
     call this%solverOptions%initialize(size(this%globalGridSizes, 1),                       &
          this%simulationFlags, this%comm)
  end if

  ! Distribute the grids between available MPI processes.
  if (this%simulationFlags%manualDomainDecomp .and.                                          &
       nProcs > size(this%globalGridSizes, 2)) then
     call getRequiredOption("decomposition_map_file",                                        &
          decompositionMapFilename, this%comm)
     call readDecompositionMap(this, decompositionMapFilename)
  end if
  call distributeGrids(this, verbose_)
  call MPI_Barrier(this%comm, ierror)

  ! Setup grids:

  allocate(this%grids(count(this%gridCommunicators /= MPI_COMM_NULL)))
  this%grids(:)%index = 0

  do i = 1, size(this%globalGridSizes, 2)
     do j = 1, size(this%grids)
        if (this%grids(j)%index /= 0) cycle
        if (this%gridCommunicators(i) /= MPI_COMM_NULL) then
           call this%grids(j)%setup(i, this%globalGridSizes(:,i), this%gridCommunicators(i), &
                this%processDistributions(:,i), simulationFlags = this%simulationFlags)
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

  ! Setup states:

  allocate(this%states(size(this%grids)))
  do i = 1, size(this%states)
     call this%states(i)%setup(this%grids(i), this%simulationFlags, this%solverOptions)
  end do
  call MPI_Barrier(this%comm, ierror)

  ! Setup patches:

  if (present(boundaryConditionFilename)) then

     call readBoundaryConditions(this, boundaryConditionFilename)
     call readPatchInterfaceInformation(this)
     call validatePatches(this)
     call distributePatches(this)
     call MPI_Barrier(this%comm, ierror)

     nPatches = 0
     if (allocated(this%patchCommunicators))                                                 &
          nPatches = count(this%patchCommunicators /= MPI_COMM_NULL)
     if (nPatches > 0) then
        allocate(this%patchFactories(nPatches))
        allocate(patchIndices(nPatches))
        patchIndices = 0
     end if

     if (allocated(this%patchData)) then

        do k = 1, size(this%grids)
           do i = 1, size(this%patchData)
              p = this%patchData(i)
              if (p%gridIndex /= this%grids(k)%index) cycle
              if (allocated(this%patchFactories)) then
                 do j = 1, size(this%patchFactories)
                    if (patchIndices(j) /= 0 .or.                                            &
                         this%patchCommunicators(i) == MPI_COMM_NULL) cycle
                    call this%patchFactories(j)%connect(patch, trim(p%patchType))
                    assert(associated(patch))
                    call patch%setup(i, this%patchCommunicators(i), p, this%grids(k),        &
                         this%simulationFlags, this%solverOptions)
                    patchIndices(j) = i
                    exit
                 end do
              end if
              call MPI_Barrier(this%grids(k)%comm, ierror)
           end do
        end do

     end if

  end if

  ! Create a communicator for master processes of grid-level communicators.
  assert(MPI_UNDEFINED /= 1)
  color = MPI_UNDEFINED
  do i = 1, size(this%grids)
     call MPI_Comm_rank(this%grids(i)%comm, procRank, ierror)
     if (procRank == 0) color = 1
  end do
  call MPI_Comm_rank(this%comm, procRank, ierror)
  call MPI_Comm_split(this%comm, color, procRank, this%commGridMasters, ierror)

  call endTiming("setupRegion")

end subroutine setupRegion

subroutine cleanupRegion(this)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region

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

  SAFE_DEALLOCATE(this%patchData)
  SAFE_DEALLOCATE(this%globalGridSizes)
  SAFE_DEALLOCATE(this%processDistributions)
  SAFE_DEALLOCATE(this%gridCommunicators)
  SAFE_DEALLOCATE(this%patchCommunicators)
  SAFE_DEALLOCATE(this%patchInterfaces)
  SAFE_DEALLOCATE(this%interfaceIndexReorderings)
  SAFE_DEALLOCATE(this%patchMasterRanks)

  if (this%comm /= MPI_COMM_NULL .and. this%comm /= MPI_COMM_WORLD)                          &
       call MPI_Comm_free(this%comm, ierror)
  this%comm = MPI_COMM_NULL

  if (this%commGridMasters /= MPI_COMM_NULL) call MPI_Comm_free(this%commGridMasters, ierror)
  this%commGridMasters = MPI_COMM_NULL

end subroutine cleanupRegion

subroutine loadRegionData(this, quantityOfInterest, filename)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use PLOT3DDescriptor_type, only : t_PLOT3DDescriptor, PLOT3D_SOLUTION_FILE

  ! <<< Enumerations >>>
  use Grid_enum
  use State_enum, only : QOI_DUMMY_FUNCTION

  ! <<< Internal modules >>>
  use ErrorHandler, only : gracefulExit, writeAndFlush
  use PLOT3DHelper
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: this
  integer, intent(in) :: quantityOfInterest
  character(len = *), intent(in) :: filename

  ! <<< Local variables >>>
  character(len = STRING_LENGTH) :: message
  logical :: success
  integer :: i, j, errorRank, procRank, ierror
  integer(kind = MPI_OFFSET_KIND) :: offset
  logical :: isSolutionFile
  real(SCALAR_KIND) :: auxiliaryData(4)

  call startTiming("loadRegionData")

  write(message, '(3A)') "Reading '", trim(filename), "'..."
  call writeAndFlush(this%comm, output_unit, message, advance = 'no')

  isSolutionFile = .false.

  do i = 1, size(this%gridCommunicators)

     success = .true.

     do j = 1, size(this%grids)
        if (this%grids(j)%index == i) then !... read one grid at a time

           offset = plot3dGetOffset(this%gridCommunicators(i), filename, i, success)
           if (.not. success) exit

           select case(quantityOfInterest)
           case (QOI_GRID, QOI_JACOBIAN, QOI_METRICS, QOI_TARGET_MOLLIFIER,                  &
                QOI_CONTROL_MOLLIFIER)
              call this%grids(j)%loadData(quantityOfInterest,                                &
                   trim(filename), offset, success)
           case default
              isSolutionFile = (quantityOfInterest /= QOI_DUMMY_FUNCTION)
              call this%states(j)%loadData(this%grids(j), quantityOfInterest,                &
                   trim(filename), offset, success)
           end select

           exit
        end if
     end do

     call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL,                               &
          MPI_LAND, this%comm, ierror)
     if (.not. success) exit
     call MPI_Barrier(this%comm, ierror)

  end do

  if (isSolutionFile) then
     auxiliaryData = real(this%states(1)%plot3dAuxiliaryData, SCALAR_KIND)
     call MPI_Bcast(auxiliaryData, 4, REAL_TYPE_MPI, 0, this%comm, ierror)
     do i = 1, 4
        this%states(:)%plot3dAuxiliaryData(i) = auxiliaryData(i)
     end do
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
     call MPI_Bcast(plot3dErrorMessage, STRING_LENGTH, MPI_CHARACTER,                        &
          errorRank, this%comm, ierror)
     call gracefulExit(this%comm, plot3dErrorMessage)
  end if

  call endTiming("loadRegionData")

end subroutine loadRegionData

subroutine saveRegionData(this, quantityOfInterest, filename)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
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
  character(len = STRING_LENGTH) :: message
  logical :: success
  integer :: i, j, nScalars, errorRank, procRank, ierror
  integer(kind = MPI_OFFSET_KIND) :: offset

  call startTiming("saveRegionData")

  write(message, '(3A)') "Writing '", trim(filename), "'..."
  call writeAndFlush(this%comm, output_unit, message, advance = 'no')

  select case(quantityOfInterest)
  case (QOI_GRID)
     call plot3dWriteSkeleton(this%comm, trim(filename),                                     &
          PLOT3D_GRID_FILE, this%globalGridSizes, success)
  case (QOI_JACOBIAN, QOI_TARGET_MOLLIFIER, QOI_CONTROL_MOLLIFIER)
     call plot3dWriteSkeleton(this%comm, trim(filename),                                     &
          PLOT3D_FUNCTION_FILE, this%globalGridSizes, success, 1)
  case (QOI_METRICS)
     call plot3dWriteSkeleton(this%comm, trim(filename),                                     &
          PLOT3D_FUNCTION_FILE, this%globalGridSizes, success,                               &
          size(this%globalGridSizes, 1) ** 2)

  case (QOI_DUMMY_FUNCTION)

     nScalars = huge(1)
     do i = 1, size(this%states)
        assert(associated(this%states(i)%dummyFunction))
        nScalars = min(nScalars, size(this%states(i)%dummyFunction, 2))
     end do
     call MPI_Allreduce(MPI_IN_PLACE, nScalars, 1,                                           &
          MPI_INTEGER, MPI_MIN, this%comm, ierror)
#ifdef DEBUG
     do i = 1, size(this%states)
        assert(size(this%states(i)%dummyFunction, 2) == nScalars)
     end do
#endif
     call plot3dWriteSkeleton(this%comm, trim(filename), PLOT3D_FUNCTION_FILE,               &
          this%globalGridSizes, success, nScalars)

  case default
     call plot3dWriteSkeleton(this%comm, trim(filename), PLOT3D_SOLUTION_FILE,               &
          this%globalGridSizes, success)
  end select

  do i = 1, size(this%gridCommunicators)

     success = .true.

     do j = 1, size(this%grids)
        if (this%grids(j)%index == i) then !... read one grid at a time.

           offset = plot3dGetOffset(this%gridCommunicators(i), filename, i, success)
           if (.not. success) exit

           select case(quantityOfInterest)
           case (QOI_GRID, QOI_JACOBIAN, QOI_METRICS, QOI_TARGET_MOLLIFIER,                  &
                QOI_CONTROL_MOLLIFIER)
              call this%grids(j)%saveData(quantityOfInterest,                                &
                   trim(filename), offset, success)
           case default
              call this%states(j)%saveData(this%grids(j), quantityOfInterest,                &
                   trim(filename), offset, success)
           end select

           exit
        end if
     end do

     call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL,                               &
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
     call MPI_Bcast(plot3dErrorMessage, STRING_LENGTH, MPI_CHARACTER,                        &
          errorRank, this%comm, ierror)
     call gracefulExit(this%comm, plot3dErrorMessage)
  end if

  call endTiming("saveRegionData")

end subroutine saveRegionData

function getCfl(this) result(cfl)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region

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
        cfl = max(cfl, this%states(i)%computeCfl(this%grids(i),                              &
             this%simulationFlags, this%solverOptions))
     end do
     call MPI_Allreduce(MPI_IN_PLACE, cfl, 1, REAL_TYPE_MPI, MPI_MAX, this%comm, ierror)
  end if

end function getCfl

function getTimeStepSize(this) result(timeStepSize)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region

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
        timeStepSize = min(timeStepSize,                                                     &
             this%states(i)%computeTimeStepSize(this%grids(i),                               &
             this%simulationFlags, this%solverOptions))
     end do
     call MPI_Allreduce(MPI_IN_PLACE, timeStepSize, 1, REAL_TYPE_MPI,                        &
          MPI_MIN, this%comm, ierror)
  else
     timeStepSize = this%solverOptions%timeStepSize
  end if

end function getTimeStepSize

subroutine computeRhs(this, mode, time)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use RhsHelper, only : computeRhsForward, computeRhsAdjoint
  use InterfaceHelper, only : exchangeInterfaceData
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: this
  integer, intent(in) :: mode
  real(SCALAR_KIND), intent(in) :: time

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j
  class(t_Patch), pointer :: patch => null()

  call startTiming("computeRhs")

  ! Semi-discrete right-hand-side operator.
  do i = 1, size(this%states)
     select case (mode)
     case (FORWARD)
        call computeRhsForward(time, this%simulationFlags, this%solverOptions,               &
             this%grids(i), this%states(i))
     case (ADJOINT)
        call computeRhsAdjoint(time, this%simulationFlags, this%solverOptions,               &
             this%grids(i), this%states(i))
     end select
  end do

  ! Multiply by Jacobian.
  do i = 1, size(this%states)
     do j = 1, this%solverOptions%nUnknowns
        this%states(i)%rightHandSide(:,j) = this%states(i)%rightHandSide(:,j) *              &
                this%grids(i)%jacobian(:,1)
     end do
  end do

  ! Update patches.
  if (allocated(this%patchFactories)) then
     do i = 1, size(this%patchFactories)
        call this%patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        do j = 1, size(this%states)
           if (patch%gridIndex /= this%grids(j)%index) cycle
           call patch%update(this%simulationFlags, this%solverOptions,                       &
                this%grids(j), this%states(j))
        end do
     end do
  end if

  ! Exchange data at block interfaces.
  call exchangeInterfaceData(this)

  ! Add patch penalties.
  if (allocated(this%patchFactories)) then
     do i = 1, size(this%patchFactories)
        call this%patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        do j = 1, size(this%states)
           if (patch%gridIndex /= this%grids(j)%index) cycle
           call patch%updateRhs(mode, this%simulationFlags, this%solverOptions,              &
                this%grids(j), this%states(j))
        end do
     end do
  end if

  ! Source terms.
  do i = 1, size(this%states)
     call this%states(i)%addSources(mode, time, this%grids(i))
  end do

  ! Zero out right-hand-side in holes.
  do i = 1, size(this%states)
     do j = 1, this%solverOptions%nUnknowns
        where (this%grids(i)%iblank == 0)
           this%states(i)%rightHandSide(:,j) = 0.0_wp
        end where
     end do
  end do

  call endTiming("computeRhs")

end subroutine computeRhs

subroutine reportGridDiagnostics(this)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Region_mod, only : t_Region

  ! <<< Internal modules >>>
  use ErrorHandler, only : writeAndFlush

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: this

  ! <<< Local variables >>>
  integer :: i, j, iGlobal, jGlobal, kGlobal, procRank, ierror
  character(len = STRING_LENGTH) :: str
  SCALAR_TYPE :: minimumJacobian, maximumJacobian

  call MPI_Comm_rank(this%comm, procRank, ierror)

  call writeAndFlush(this%comm, output_unit, "")
  write(str, '(A,I0.0,A,I2,A)') "Region is ", size(this%globalGridSizes, 1),                 &
       "D and has ", size(this%globalGridSizes, 2), " block(s):"
  call writeAndFlush(this%comm, output_unit, str)
  write(str, '(A)') repeat('-', 32)
  call writeAndFlush(this%comm, output_unit, str)

  do i = 1, size(this%globalGridSizes, 2)
     do j = 1, size(this%grids)
        if (this%grids(j)%index == i) then

           write(str, '(2X,A,I2,3(A,I4),A)') "Block ", i, ": ",                              &
                this%grids(j)%globalSize(1), " x ",                                          &
                this%grids(j)%globalSize(2), " x ",                                          &
                this%grids(j)%globalSize(3), " points"
           call writeAndFlush(this%grids(j)%comm, output_unit, str)

           call this%grids(j)%findMinimum(this%grids(j)%jacobian(:,1),                       &
                minimumJacobian, iGlobal, jGlobal, kGlobal)
           write(str, '(4X,A,(SS,ES9.2E2),3(A,I4),A)') "min. Jacobian = ",                   &
                real(minimumJacobian, SCALAR_KIND), " at (",                                 &
                iGlobal, ", ", jGlobal, ", ", kGlobal, ")"
           call writeAndFlush(this%grids(j)%comm, output_unit, str)

           call this%grids(j)%findMaximum(this%grids(j)%jacobian(:,1),                       &
                maximumJacobian, iGlobal, jGlobal, kGlobal)
           write(str, '(4X,A,(SS,ES9.2E2),3(A,I4),A)') "max. Jacobian = ",                   &
                real(maximumJacobian, SCALAR_KIND), " at (",                                 &
                iGlobal, ", ", jGlobal, ", ", kGlobal, ")"
           call writeAndFlush(this%grids(j)%comm, output_unit, str)

        end if
        call MPI_Barrier(this%grids(j)%comm, ierror)
     end do
     call MPI_Barrier(this%comm, ierror)
  end do

  call writeAndFlush(this%comm, output_unit, "")

end subroutine reportGridDiagnostics
