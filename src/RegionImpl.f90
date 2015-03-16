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
    use Region_type

    ! <<< Internal modules >>>
    use MPIHelper, only : gracefulExit, writeAndFlush
    use InputHelper, only : stripComments

    ! <<< Arguments >>>
    type(t_Region) :: this
    character(len = *), intent(in) :: filename

    ! <<< Local variables >>>
    integer, parameter :: fileUnit = 31
    integer :: i, proc, nProcs, lineNo, gridIndex, numProcsInGrid(3), istat, ierror
    character(len = STRING_LENGTH) :: line, message
    character(len = 1), parameter :: commentMarker = '#'

    call MPI_Comm_rank(this%comm, proc, ierror)
    call MPI_Comm_size(this%comm, nProcs, ierror)

    write(message, "(2A)") "Reading MPI decomposition map from '", trim(filename), "'..."
    call writeAndFlush(this%comm, output_unit, message)

    ! Check if file exists.
    if (proc == 0) then
       open(unit = fileUnit, file = trim(filename), action = 'read', status = 'old',         &
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
       open(unit = fileUnit, file = trim(filename), action = 'read', status = 'old')
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
               any(numProcsInGrid < 0) .or. sum(numProcsInGrid) > nProcs) then
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

          this%processDistributions(i,:) = numProcsInGrid(1:size(this%processDistributions,1))

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
       write(message, '(A,2(A,I0.0),A)') filename,                                           &
            ": Invalid process distribution: expected a total of ", nProcs,                  &
            " processes, got ", sum(product(this%processDistributions, dim = 1)),            &
            " processes!"
       call gracefulExit(this%comm, message)
    end if

  end subroutine readDecompositionMap

  subroutine distributeGrids(this)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env

    ! <<< Derived types >>>
    use Region_type

    ! <<< Internal modules >>>
    use MPIHelper, only : splitCommunicatorMultigrid, writeAndFlush

    ! <<< Arguments >>>
    type(t_Region) :: this

    ! <<< Local variables >>>
    integer :: nProcs, ierror
    character(len = STRING_LENGTH) :: message
    integer, allocatable :: numProcsInGrid(:)

    ! Find the size of the communicator.
    call MPI_Comm_size(this%comm, nProcs, ierror)

    write(message, "(2(A,I0.0),A)") "Distributing ", size(this%globalGridSizes, 2),          &
         " grid(s) across ", nProcs, " process(es)..."
    call writeAndFlush(this%comm, output_unit, message)

    ! Find the number of processes to be assigned to each grid: `numProcsInGrid(i)` is the
    ! number of processes assigned to grid `i`.
    allocate(numProcsInGrid(size(this%globalGridSizes, 2)))
    numProcsInGrid = sum(this%processDistributions, dim = 1)
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
    use Region_type

    ! <<< Internal modules >>>
    use MPIHelper, only : gracefulExit, writeAndFlush
    use InputHelper, only : stripComments
    use PatchDescriptor_mod, only : parsePatchType

    ! <<< Arguments >>>
    type(t_Region) :: this
    character(len = *), intent(in) :: filename

    ! <<< Local variables >>>
    integer, parameter :: fileUnit = 41
    integer :: i, proc, nPatches, lineNo, istat, ierror
    character(len = STRING_LENGTH) :: line, patchTypeString, message
    character(len = 1), parameter :: commentMarker = '#'
    integer, allocatable :: tempBuffer(:,:)

    call MPI_Comm_rank(this%comm, proc, ierror)

    ! Check if file exists.
    if (proc == 0) then
       open(unit = fileUnit, file = trim(filename), action = 'read', status = 'old',         &
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
    allocate(tempBuffer(nPatches, 9))

    ! Again, only the root process reads the file.
    if (proc == 0) then
       i = 0; lineNo = 0; istat = 0
       open(unit = fileUnit, file = trim(filename), action = 'read', status = 'old')
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
          read(line, *, iostat = istat) this%patchData(i)%name, patchTypeString,             &
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

          ! Pack patch data into tempBuffer for broadcasting.
          tempBuffer(i,1) = this%patchData(i)%gridIndex
          tempBuffer(i,2) = this%patchData(i)%normalDirection
          call parsePatchType(patchTypeString, tempBuffer(i,3)) !... validate.
          if (tempBuffer(i,3) == -1) then
             istat = -1
             write(message, "(2A,I0.0,3A)") trim(filename), ":", lineNo,                     &
                  ": Invalid type for patch '", this%patchData(i)%name, "'!"
             exit
          end if
          tempBuffer(i,4:9) = (/ this%patchData(i)%iMin, this%patchData(i)%iMax,             &
               this%patchData(i)%jMin, this%patchData(i)%jMax,                               &
               this%patchData(i)%kMin, this%patchData(i)%kMax /)

       end do
       close(fileUnit)
    end if

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
       this%patchData(i)%patchType = tempBuffer(i,3)
       this%patchData(i)%iMin = tempBuffer(i,4)
       this%patchData(i)%iMax = tempBuffer(i,5)
       this%patchData(i)%jMin = tempBuffer(i,6)
       this%patchData(i)%jMax = tempBuffer(i,7)
       this%patchData(i)%kMin = tempBuffer(i,8)
       this%patchData(i)%kMax = tempBuffer(i,9)
       call MPI_Bcast(this%patchData(i)%name, len(this%patchData(i)%name),                   &
            MPI_CHARACTER, 0, this%comm, ierror)
    end do

    SAFE_DEALLOCATE(tempBuffer)

  end subroutine readBoundaryConditions

  subroutine validatePatches(this)

    ! <<< External modules >>>
    use, intrinsic :: iso_fortran_env

    ! <<< Derived types >>>
    use Region_type

    ! <<< Internal modules >>>
    use MPIHelper, only : issueWarning, gracefulExit, writeAndFlush
    use PatchDescriptor_mod, only : validatePatchDescriptor, validatePatchesConnectivity

    ! <<< Arguments >>>
    type(t_Region) :: this

    ! <<< Local variables >>>
    integer :: i, errorCode
    character(len = STRING_LENGTH) :: message

    if (.not. allocated(this%patchData)) return

    write(message, "(A,I0.0,A)") "Validating boundary conditions..."
    call writeAndFlush(this%comm, output_unit, message)

    do i = 1, size(this%patchData)
       call validatePatchDescriptor(this%patchData(i), this%globalGridSizes,                 &
            this%simulationFlags, errorCode, message)
       if (errorCode == 1) then
          call issueWarning(this%comm, message)
       else if (errorCode == 2) then
          call gracefulExit(this%comm, message)
       end if
    end do

    call validatePatchesConnectivity(this%patchData, errorCode, message)
    if (errorCode == 2) call gracefulExit(this%comm, message)

  end subroutine validatePatches

  subroutine distributePatches(this)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_type
    use PatchDescriptor_type

    ! <<< Arguments >>>
    type(t_Region) :: this

    ! <<< Local variables >>>
    integer :: i, j, gridOffset(3), gridLocalSize(3), gridIndex, color, comm, proc, ierror
    type(t_PatchDescriptor) :: p

    if (.not. allocated(this%patchData)) return

    SAFE_DEALLOCATE(this%patchCommunicators)
    allocate(this%patchCommunicators(size(this%patchData)), source = MPI_COMM_NULL)

    do i = 1, size(this%grids)

       gridIndex = this%grids(i)%index
       gridOffset = this%grids(i)%offset
       gridLocalSize = this%grids(i)%localSize

       call MPI_Comm_rank(this%grids(i)%comm, proc, ierror)

       do j = 1, size(this%patchData)

          p = this%patchData(j)
          color = p%patchType
          if (p%gridIndex /= gridIndex .or.                                                  &
               p%iMax < gridOffset(1) + 1 .or.                                               &
               p%iMin > gridOffset(1) + gridLocalSize(1) .or.                                &
               p%jMax < gridOffset(2) + 1 .or.                                               &
               p%jMin > gridOffset(2) + gridLocalSize(2) .or.                                &
               p%kMax < gridOffset(3) + 1 .or.                                               &
               p%kMin > gridOffset(3) + gridLocalSize(3)) then
             color = MPI_UNDEFINED
          end if
          call MPI_Comm_split(this%grids(i)%comm, color, proc, comm, ierror)
          if (comm /= MPI_COMM_NULL) this%patchCommunicators(j) = comm

       end do

       call MPI_Barrier(this%grids(i)%comm, ierror)

    end do

  end subroutine distributePatches

  subroutine checkSolutionLimits(this, mode)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use State_type
    use Region_type

    ! <<< Public members >>>
    use Region_mod, only : saveRegionData

    ! <<< Internal modules >>>
    use Grid_mod, only : isVariableWithinRange
    use MPIHelper, only : gracefulExit

    ! <<< Arguments >>>
    type(t_Region) :: this
    integer, intent(in) :: mode

    ! <<< Local variables >>>
    integer :: i, iGlobal, jGlobal, kGlobal, rankReportingError, procRank, ierror
    character(len = STRING_LENGTH) :: message
    SCALAR_TYPE :: fOutsideRange

    rankReportingError = -1
    call MPI_Comm_rank(this%comm, procRank, ierror)

    do i = 1, size(this%states)

       if (.not. isVariableWithinRange(this%grids(i),                                        &
            this%states(i)%conservedVariables(:,1), fOutsideRange,                           &
            iGlobal, jGlobal, kGlobal,                                                       &
            minValue = this%solverOptions%densityRange(1),                                   &
            maxValue = this%solverOptions%densityRange(2))) then
          write(message, '(4(A,I0.0),3(A,(SS,ES9.2E2)),A)') "Density on grid ",              &
               this%grids(i)%index, " at (", iGlobal, ", ", jGlobal, ", ", kGlobal, "): ",   &
               fOutsideRange, " out of range (",                                             &
               this%solverOptions%densityRange(1), ", ",                                     &
               this%solverOptions%densityRange(2), ")!"
          rankReportingError = procRank
          exit
       end if

       if (.not. isVariableWithinRange(this%grids(i),                                        &
            this%states(i)%temperature(:,1), fOutsideRange, iGlobal, jGlobal, kGlobal,       &
            minValue = this%solverOptions%temperatureRange(1),                               &
            maxValue = this%solverOptions%temperatureRange(2))) then
          write(message, '(4(A,I0.0),3(A,(SS,ES9.2E2)),A)') "Temperature on grid ",          &
               this%grids(i)%index, " at (", iGlobal, ", ", jGlobal, ", ", kGlobal, "): ",   &
               fOutsideRange, " out of range (",                                             &
               this%solverOptions%temperatureRange(1), ", ",                                 &
               this%solverOptions%temperatureRange(2), ")!"
          rankReportingError = procRank
          exit
       end if

    end do

    call MPI_Allreduce(MPI_IN_PLACE, rankReportingError, 1,                                  &
         MPI_INTEGER, MPI_MAX, this%comm, ierror)

    if (rankReportingError /= -1) then

       if (procRank == 0 .and. rankReportingError /= 0)                                      &
            call MPI_Recv(message, STRING_LENGTH, MPI_CHARACTER, rankReportingError,         &
            rankReportingError, this%comm, MPI_STATUS_IGNORE, ierror)
       if (procRank == rankReportingError .and. rankReportingError /= 0)                     &
            call MPI_Send(message, STRING_LENGTH, MPI_CHARACTER, 0, procRank,                &
            this%comm, ierror)

       select case (mode)
       case (FORWARD)
          call saveRegionData(this, QOI_FORWARD_STATE, PROJECT_NAME // "-crashed.q")
       end select

       call gracefulExit(this%comm, message)

    end if

  end subroutine checkSolutionLimits

end module RegionImpl

subroutine setupRegion(this, comm, globalGridSizes, boundaryConditionFilename)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_type
  use PatchDescriptor_type

  ! <<< Private members >>>
  use RegionImpl

  ! <<< Public members >>>
  use Region_mod, only : cleanupRegion
  use MPITimingsHelper, only : startTiming, endTiming

  ! <<< Internal modules >>>
  use Grid_mod, only : setupGrid
  use Patch_mod, only : setupPatch, updatePatchConnectivity
  use State_mod, only : setupState
  use InputHelper, only : getRequiredOption
  use SolverOptions_mod, only : initializeSolverOptions
  use SimulationFlags_mod, only : initializeSimulationFlags

  implicit none

  ! <<< Arguments >>>
  type(t_Region) :: this
  integer, intent(in) :: comm, globalGridSizes(:,:)
  character(len = *), intent(in), optional :: boundaryConditionFilename

  ! <<< Local variables >>>
  integer :: i, j, k, nProcs, nPatches, ierror
  character(len = STRING_LENGTH) :: decompositionMapFilename
  type(t_PatchDescriptor) :: p

  call startTiming("setupRegion")

  ! Clean slate.
  call cleanupRegion(this)
  this%comm = comm
  call MPI_Comm_size(this%comm, nProcs, ierror)

  SAFE_DEALLOCATE(this%globalGridSizes)
  SAFE_DEALLOCATE(this%processDistributions)

  allocate(this%globalGridSizes(size(globalGridSizes, 1), size(globalGridSizes, 2)),         &
       source = globalGridSizes)
  allocate(this%processDistributions(size(this%globalGridSizes, 1),                          &
       size(this%globalGridSizes, 2)), source = 0)

  ! Initialize simulation flags.
  call initializeSimulationFlags(this%simulationFlags)

  ! Initialize solver options.
  call initializeSolverOptions(this%solverOptions, this%simulationFlags)

  ! Distribute the grids between available MPI processes.
  if (this%simulationFlags%manualDomainDecomp .and.                                          &
       nProcs > size(this%globalGridSizes, 2)) then
     call getRequiredOption("manual_decomposition_map_filename",                             &
          decompositionMapFilename, this%comm)
     call readDecompositionMap(this, decompositionMapFilename)
  end if
  call distributeGrids(this)
  call MPI_Barrier(this%comm, ierror)

  ! Setup grids:

  allocate(this%grids(count(this%gridCommunicators /= MPI_COMM_NULL)))
  this%grids(:)%index = 0

  do i = 1, size(this%globalGridSizes, 2)
     do j = 1, size(this%grids)
        if (this%grids(j)%index /= 0) cycle
        if (this%gridCommunicators(i) /= MPI_COMM_NULL) then
           call setupGrid(this%grids(j), i, this%globalGridSizes(:,i),                       &
                this%gridCommunicators(i), this%processDistributions(:,i),                   &
                simulationFlags = this%simulationFlags)
           exit
        end if
     end do
     call MPI_Barrier(this%comm, ierror)
  end do

  ! Setup states:

  allocate(this%states(size(this%grids)))
  do i = 1, size(this%states)
     call setupState(this%states(i), this%grids(i), this%simulationFlags, this%solverOptions)
  end do
  call MPI_Barrier(this%comm, ierror)

  ! Setup patches:

  if (present(boundaryConditionFilename)) then

     call readBoundaryConditions(this, boundaryConditionFilename)
     call validatePatches(this)
     call distributePatches(this)
     call MPI_Barrier(this%comm, ierror)

     nPatches = 0
     if (allocated(this%patchCommunicators))                                                 &
          nPatches = count(this%patchCommunicators /= MPI_COMM_NULL)
     if (nPatches > 0) then
        allocate(this%patches(nPatches))
        this%patches%index = 0
     end if

     if (allocated(this%patchData)) then

        do k = 1, size(this%grids)
           do i = 1, size(this%patchData)
              p = this%patchData(i)
              if (p%gridIndex /= this%grids(k)%index) cycle
              if (allocated(this%patches)) then
                 do j = 1, size(this%patches)
                    if (this%patches(j)%index /= 0) cycle
                    if (this%patchCommunicators(i) /= MPI_COMM_NULL) then
                       call setupPatch(this%patches(j), i, size(globalGridSizes, 1), p,      &
                            this%patchCommunicators(i), this%grids(k)%offset,                &
                            this%grids(k)%localSize, this%simulationFlags)
                       exit
                    end if
                 end do
              end if
              call MPI_Barrier(this%grids(k)%comm, ierror)
           end do
        end do

        if (allocated(this%patches)) then
           do i = 1, size(this%patches)
              call updatePatchConnectivity(this%patches(i), this%patchData)
           end do
        end if

     end if

  end if

  call endTiming("setupRegion")

end subroutine setupRegion

subroutine cleanupRegion(this)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_type

  ! <<< Internal modules >>>
  use Grid_mod, only : cleanupGrid
  use Patch_mod, only : cleanupPatch
  use State_mod, only : cleanupState

  implicit none

  ! <<< Arguments >>>
  type(t_Region) :: this

  ! <<< Local variables >>>
  integer :: i, ierror

  if (allocated(this%grids)) then
     do i = 1, size(this%grids)
        call cleanupGrid(this%grids(i))
     end do
  end if
  SAFE_DEALLOCATE(this%grids)

  if (allocated(this%states)) then
     do i = 1, size(this%states)
        call cleanupState(this%states(i))
     end do
  end if
  SAFE_DEALLOCATE(this%states)

  if (allocated(this%patches)) then
     do i = 1, size(this%patches)
        call cleanupPatch(this%patches(i))
     end do
  end if
  SAFE_DEALLOCATE(this%patches)

  if (this%comm /= MPI_COMM_NULL .and. this%comm /= MPI_COMM_WORLD)                          &
       call MPI_Comm_free(this%comm, ierror)

end subroutine cleanupRegion

subroutine loadRegionData(this, quantityOfInterest, filename)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Grid_type
  use Region_type
  use PLOT3DDescriptor_type

  ! <<< Internal modules >>>
  use Grid_mod, only : loadGridData
  use MPIHelper, only : gracefulExit, writeAndFlush
  use State_mod, only : loadStateData, getFileType
  use PLOT3DHelper
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  type(t_Region) :: this
  integer, intent(in) :: quantityOfInterest
  character(len = *), intent(in) :: filename

  ! <<< Local variables >>>
  character(len = STRING_LENGTH) :: message
  logical :: success
  integer :: i, j, errorRank, procRank, ierror
  integer(kind = MPI_OFFSET_KIND) :: offset
  real(SCALAR_KIND) :: auxiliaryData(4)

  call startTiming("loadRegionData")

  write(message, '(3A)') "Reading '", trim(filename), "'..."
  call writeAndFlush(this%comm, output_unit, message, advance = 'no')

  do i = 1, size(this%gridCommunicators)

     do j = 1, size(this%grids)
        if (this%grids(j)%index == i) then !... read one grid at a time

           offset = plot3dGetOffset(this%gridCommunicators(i), filename, i, success)
           if (.not. success) exit

           select case(quantityOfInterest)
           case (QOI_GRID, QOI_JACOBIAN, QOI_METRICS, QOI_TARGET_MOLLIFIER,                  &
                QOI_CONTROL_MOLLIFIER)
              call loadGridData(this%grids(j), quantityOfInterest,                           &
                   trim(filename), offset, success)
           case default
              call loadStateData(this%states(j), this%grids(j),                              &
                   quantityOfInterest, trim(filename), offset, success)
           end select

           exit
        end if
     end do

     call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL,                               &
          MPI_LAND, MPI_COMM_WORLD, ierror)
     if (.not. success) exit
     call MPI_Barrier(this%comm, ierror)

  end do

  if (getFileType(quantityOfInterest) == PLOT3D_SOLUTION_FILE) then
     auxiliaryData = this%states(1)%plot3dAuxiliaryData
#ifdef SCALAR_IS_COMPLEX
     call MPI_Bcast(auxiliaryData, 4, REAL_TYPE_MPI, 0, this%comm, ierror)
#else
     call MPI_Bcast(auxiliaryData, 4, SCALAR_TYPE_MPI, 0, this%comm, ierror)
#endif
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
     call MPI_Bcast(plot3dErrorMessage, STRING_LENGTH, MPI_CHARACTER, &
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
  use Grid_type
  use Region_type
  use PLOT3DDescriptor_type

  ! <<< Internal modules >>>
  use Grid_mod, only : saveGridData
  use MPIHelper, only : writeAndFlush
  use State_mod, only : saveStateData, getFileType, getNumberOfScalars
  use PLOT3DHelper
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  type(t_Region) :: this
  integer, intent(in) :: quantityOfInterest
  character(len = *), intent(in) :: filename

  ! <<< Local variables >>>
  character(len = STRING_LENGTH) :: message
  logical :: success
  integer :: i, j, fileType, errorRank, procRank, ierror
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
  case default
     fileType = getFileType(quantityOfInterest)
     if (fileType == PLOT3D_FUNCTION_FILE) then
        call plot3dWriteSkeleton(this%comm, trim(filename), fileType, this%globalGridSizes,  &
             success, getNumberOfScalars(quantityOfInterest, size(this%globalGridSizes, 1)))
     else
        call plot3dWriteSkeleton(this%comm, trim(filename), fileType,                        &
             this%globalGridSizes, success)
     end if
  end select

  do i = 1, size(this%gridCommunicators)

     do j = 1, size(this%grids)
        if (this%grids(j)%index == i) then !... read one grid at a time.

           offset = plot3dGetOffset(this%gridCommunicators(i), filename, i, success)
           if (.not. success) exit

           select case(quantityOfInterest)
           case (QOI_GRID, QOI_JACOBIAN, QOI_METRICS, QOI_TARGET_MOLLIFIER,                  &
                QOI_CONTROL_MOLLIFIER)
              call saveGridData(this%grids(j), quantityOfInterest,                           &
                   trim(filename), offset, success)
           case default
              call saveStateData(this%states(j), this%grids(j),                              &
                   quantityOfInterest, trim(filename), offset, success)
           end select

           exit
        end if
     end do

     call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL,                               &
          MPI_LAND, MPI_COMM_WORLD, ierror)
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
     call MPI_Bcast(plot3dErrorMessage, STRING_LENGTH, MPI_CHARACTER, &
          errorRank, this%comm, ierror)
     call gracefulExit(this%comm, plot3dErrorMessage)
  end if

  call endTiming("saveRegionData")

end subroutine saveRegionData

subroutine computeRhs(this, mode, time)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_type

  ! <<< Private members >>>
  use RegionImpl, only : checkSolutionLimits

  ! <<< Internal modules >>>
  use State_mod
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  type(t_Region) :: this
  integer, intent(in) :: mode
  real(SCALAR_KIND), intent(in) :: time

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, ierror
  real(wp) :: timeStepSize, cfl

  call startTiming("computeRhs")

  do i = 1, size(this%states)
     this%states(i)%rightHandSide = 0.0_wp
     call updateState(this%states(i), this%grids(i), time,                                   &
          this%simulationFlags, this%solverOptions)
  end do

  if (this%simulationFlags%useConstantCfl) then
     timeStepSize = minval(this%states(:)%timeStepSize)
     call MPI_Allreduce(MPI_IN_PLACE, timeStepSize, 1, SCALAR_TYPE_MPI, MPI_MIN,             &
          this%comm, ierror)
     this%states(:)%timeStepSize = timeStepSize
  else
     cfl = maxval(this%states(:)%cfl)
     call MPI_Allreduce(MPI_IN_PLACE, cfl, 1, SCALAR_TYPE_MPI, MPI_MAX, this%comm, ierror)
     this%states(:)%cfl = cfl
  end if

  if (this%simulationFlags%enableSolutionLimits) call checkSolutionLimits(this, mode)

  do i = 1, size(this%states)

     ! Semi-discrete right-hand-side operator.
     select case (mode)
     case (FORWARD)
        call computeRhsForward(this%states(i), this%grids(i), this%patches, time,            &
             this%simulationFlags, this%solverOptions)
     case (ADJOINT)
        call computeRhsAdjoint(this%states(i), this%grids(i), this%patches, time,            &
             this%simulationFlags, this%solverOptions)
     end select

  end do

  do i = 1, size(this%states)

     ! SAT penalties.
     select case (mode)
     case (FORWARD)
        call addPenaltiesForward(this%states(i), this%grids(i), this%patches, time,          &
             this%simulationFlags, this%solverOptions)
     case (ADJOINT)
        call addPenaltiesAdjoint(this%states(i), this%grids(i), this%patches, time,          &
             this%simulationFlags, this%solverOptions)
     end select

     ! Multiply by Jacobian and zero-out at hole points.
     do j = 1, this%states(i)%nUnknowns
        where (this%grids(i)%iblank == 0)
           this%states(i)%rightHandSide(:,j) = 0
        elsewhere
           this%states(i)%rightHandSide(:,j) = this%states(i)%rightHandSide(:,j) *           &
                this%grids(i)%jacobian(:,1)
        end where
     end do

     ! Source terms.
     select case (mode)
     case (FORWARD)
        call addSourcesForward(this%states(i), this%grids(i), this%patches, time,            &
             this%simulationFlags, this%solverOptions)
     case (ADJOINT)
        call addSourcesAdjoint(this%states(i), this%grids(i), this%patches, time,            &
             this%simulationFlags, this%solverOptions)
     end select

  end do

  call endTiming("computeRhs")

end subroutine computeRhs

subroutine reportGridDiagnostics(this)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Region_type

  ! <<< Internal modules >>>
  use Grid_mod, only : findMinimum, findMaximum
  use MPIHelper, only : writeAndFlush

  implicit none

  ! <<< Arguments >>>
  type(t_Region) :: this

  ! <<< Local variables >>>
  integer :: i, j, iGlobal, jGlobal, kGlobal, procRank, ierror
  character(len = STRING_LENGTH) :: str
  real(SCALAR_KIND) :: minimumJacobian, maximumJacobian

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

           call findMinimum(this%grids(j), this%grids(j)%jacobian(:,1),                      &
                minimumJacobian, iGlobal, jGlobal, kGlobal)
           write(str, '(4X,A,(SS,ES9.2E2),3(A,I4),A)') "min. Jacobian = ",                   &
                minimumJacobian, " at (", iGlobal, ", ", jGlobal, ", ", kGlobal, ")"
           call writeAndFlush(this%grids(j)%comm, output_unit, str)

           call findMaximum(this%grids(j), this%grids(j)%jacobian(:,1),                      &
                maximumJacobian, iGlobal, jGlobal, kGlobal)
           write(str, '(4X,A,(SS,ES9.2E2),3(A,I4),A)') "max. Jacobian = ",                   &
                maximumJacobian, " at (", iGlobal, ", ", jGlobal, ", ", kGlobal, ")"
           call writeAndFlush(this%grids(j)%comm, output_unit, str)

        end if
        call MPI_Barrier(this%grids(j)%comm, ierror)
     end do
     call MPI_Barrier(this%comm, ierror)
  end do

  call writeAndFlush(this%comm, output_unit, "")

end subroutine reportGridDiagnostics

subroutine subStepHooks(this, timestep, stage)

  ! <<< Derived types >>>
  use Region_type

  implicit none

  ! <<< Arguments >>>
  type(t_Region) :: this
  integer, intent(in) :: timestep, stage

end subroutine subStepHooks

subroutine reportResiduals(this)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env

  ! <<< Derived types >>>
  use Region_type

  ! <<< Internal modules >>>
  use Grid_mod, only : findMaximum
  use MPIHelper, only : writeAndFlush

  implicit none

  ! <<< Arguments >>>
  type(t_Region) :: this

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions, ierror
  SCALAR_TYPE, allocatable :: f(:)
  SCALAR_TYPE :: fMax
  real(wp) :: residuals(3)
  character(len = STRING_LENGTH) :: str

  do i = 1, size(this%states)

     call MPI_Cartdim_get(this%grids(i)%comm, nDimensions, ierror)

     allocate(f(size(this%states(i)%rightHandSide, 1)))

     f = abs(this%states(i)%rightHandSide(:,1))
     call findMaximum(this%grids(i), f, fMax)
     residuals(1) = real(fMax, wp)

     residuals(2) = 0.0_wp
     do j = 1, nDimensions
        f = abs(this%states(i)%rightHandSide(:,j+1))
        call findMaximum(this%grids(i), f, fMax)
        residuals(2) = max(residuals(2), real(fMax, wp))
     end do

     f = abs(this%states(i)%rightHandSide(:,nDimensions+2))
     call findMaximum(this%grids(i), f, fMax)
     residuals(3) = real(fMax, wp)

     SAFE_DEALLOCATE(f)

  end do

#ifdef SCALAR_IS_COMPLEX
  call MPI_Allreduce(MPI_IN_PLACE, residuals, 3, REAL_TYPE_MPI, MPI_MAX, this%comm, ierror)
#else
  call MPI_Allreduce(MPI_IN_PLACE, residuals, 3, SCALAR_TYPE_MPI, MPI_MAX, this%comm, ierror)
#endif
  write(str, '(2X,3(A,(ES11.4E2)))') "residuals: density = ", residuals(1),                  &
       ", momentum = ", residuals(2), ", energy = ", residuals(3)
  call writeAndFlush(this%comm, output_unit, str)

end subroutine reportResiduals
