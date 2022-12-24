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

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

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
    use Region_mod, only : t_Region
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
    class(t_Patch), pointer :: dummyPatch => null()
    integer, allocatable :: tempBuffer(:,:)

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
    integer :: i, j, gridOffset(3), gridLocalSize(3), gridIndex, color, comm,                &
         rankInGridCommunicator, rankInRegionCommunicator, rankInPatchCommunicator, ierror
    type(t_PatchDescriptor) :: p

    if (.not. allocated(this%patchData)) return

    SAFE_DEALLOCATE(this%patchCommunicators)
    SAFE_DEALLOCATE(this%patchMasterRanks)

    allocate(this%patchCommunicators(size(this%patchData)), source = MPI_COMM_NULL)
    allocate(this%patchMasterRanks(size(this%patchData)), source = -1)

    call MPI_Comm_rank(this%comm, rankInRegionCommunicator, ierror)

    do i = 1, size(this%grids)

       gridIndex = this%grids(i)%index
       gridOffset = this%grids(i)%offset
       gridLocalSize = this%grids(i)%localSize

       call MPI_Comm_rank(this%grids(i)%comm, rankInGridCommunicator, ierror)

       do j = 1, size(this%patchData)

          if (this%patchData(j)%gridIndex /= gridIndex) cycle

          p = this%patchData(j)
          color = 1
          if (p%iMax < gridOffset(1) + 1 .or.                                                &
               p%iMin > gridOffset(1) + gridLocalSize(1) .or.                                &
               p%jMax < gridOffset(2) + 1 .or.                                               &
               p%jMin > gridOffset(2) + gridLocalSize(2) .or.                                &
               p%kMax < gridOffset(3) + 1 .or.                                               &
               p%kMin > gridOffset(3) + gridLocalSize(3)) then
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
    use Region_mod, only : t_Region

    ! <<< Internal modules >>>
    use ErrorHandler, only : gracefulExit, issueWarning
    use Patch_factory, only : computeQuadratureOnPatches

    ! <<< Arguments >>>
    class(t_Region) :: this

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: mollifierNorm, timeStepFactor
    integer :: i, ierror, mpiOperator
    logical :: hasNegativeMollifier
    character(len = STRING_LENGTH) :: str, controllerNorm

    assert(allocated(this%grids))

    controllerNorm = trim(this%solverOptions%controllerNorm)
    assert_key(controllerNorm,('L1','L_Inf_with_timestep'))
    select case (controllerNorm)
    case('L1')
      mpiOperator = MPI_SUM
    case('L_Inf_with_timestep')
      mpiOperator = MPI_MAX
      timeStepFactor = sqrt( this%solverOptions%timeStepSize / this%solverOptions%controllerFactor )
    case default
      mpiOperator = -1
      timeStepFactor = - 1.0_wp
    end select

    mollifierNorm = 0.0_wp

    do i = 1, size(this%grids)

       assert(allocated(this%grids(i)%controlMollifier))

       hasNegativeMollifier = any(real(this%grids(i)%controlMollifier(:,1), wp) < 0.0_wp)
       call MPI_Allreduce(MPI_IN_PLACE, hasNegativeMollifier, 1,                             &
            MPI_LOGICAL, MPI_LOR, this%grids(i)%comm, ierror)
       if (hasNegativeMollifier) then
          write(str, '(A,I0.0,A)') "Control mollifying support function on grid ",           &
               this%grids(i)%index, " is not non-negative everywhere!"
          call gracefulExit(this%grids(i)%comm, str)
       end if

       select case(controllerNorm)
       case ('L1')
         mollifierNorm = mollifierNorm +                                                     &
              real(computeQuadratureOnPatches(this%patchFactories,                           &
              'ACTUATOR', this%grids(i), this%grids(i)%controlMollifier(:,1)), wp)
       case ('L_Inf_with_timestep')
         mollifierNorm = max(mollifierNorm,                                                  &
                            timeStepFactor * maxval(this%grids(i)%controlMollifier(:,1)) )
         call MPI_Allreduce(MPI_IN_PLACE, mollifierNorm, 1, REAL_TYPE_MPI,                   &
                            mpiOperator, this%grids(i)%comm, ierror)
       case default
         write(str, '(A)') "Solver Option 'controller_norm' is not specified!"
         call gracefulExit(this%grids(i)%comm, str)
       end select
    end do

    if (this%commGridMasters /= MPI_COMM_NULL)                                               &
         call MPI_Allreduce(MPI_IN_PLACE, mollifierNorm, 1, REAL_TYPE_MPI,                   &
         mpiOperator, this%commGridMasters, ierror)

    do i = 1, size(this%grids)
       call MPI_Bcast(mollifierNorm, 1, REAL_TYPE_MPI, 0, this%grids(i)%comm, ierror)
    end do

    if (mollifierNorm <= 0.0_wp)                                                             &
         call issueWarning(this%comm,                                                        &
         "Control mollifying support is trivial! Is an actuator patch present?")

    do i = 1, size(this%grids)
       this%grids(i)%controlMollifier = this%grids(i)%controlMollifier / mollifierNorm
    end do

  end subroutine normalizeControlMollifier

  subroutine normalizeTargetMollifier(this)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Internal modules >>>
    use ErrorHandler, only : gracefulExit, issueWarning
    use Patch_factory, only : computeQuadratureOnPatches

    ! <<< Arguments >>>
    class(t_Region) :: this

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: mollifierNorm
    integer :: i, ierror
    logical :: hasNegativeMollifier
    character(len = STRING_LENGTH) :: str

    assert(allocated(this%grids))

    mollifierNorm = 0.0_wp

    do i = 1, size(this%grids)

       assert(allocated(this%grids(i)%targetMollifier))

       hasNegativeMollifier = any(real(this%grids(i)%targetMollifier(:,1), wp) < 0.0_wp)
       call MPI_Allreduce(MPI_IN_PLACE, hasNegativeMollifier, 1,                             &
            MPI_LOGICAL, MPI_LOR, this%grids(i)%comm, ierror)
       if (hasNegativeMollifier) then
          write(str, '(A,I0.0,A)') "Target mollifying support function on grid ",            &
               this%grids(i)%index, " is not non-negative everywhere!"
          call gracefulExit(this%grids(i)%comm, str)
       end if
       mollifierNorm = mollifierNorm +                                                       &
            real(computeQuadratureOnPatches(this%patchFactories,                             &
            'COST_TARGET', this%grids(i), this%grids(i)%targetMollifier(:,1)), wp)
    end do

    if (this%commGridMasters /= MPI_COMM_NULL)                                               &
         call MPI_Allreduce(MPI_IN_PLACE, mollifierNorm, 1, REAL_TYPE_MPI,                   &
         MPI_SUM, this%commGridMasters, ierror)

    do i = 1, size(this%grids)
       call MPI_Bcast(mollifierNorm, 1, REAL_TYPE_MPI, 0, this%grids(i)%comm, ierror)
    end do
    if (mollifierNorm <= 0.0_wp)                                                             &
         call issueWarning(this%comm,                                                        &
         "Target mollifying support is trivial! Is a cost target patch present?")

    do i = 1, size(this%grids)
       this%grids(i)%targetMollifier = this%grids(i)%targetMollifier / mollifierNorm
    end do

  end subroutine normalizeTargetMollifier

  function computeRegionIntegral(region, quantityOfInterest, index) result(quadrature)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Enumerations >>>
    use State_enum, only : QOI_ADJOINT_STATE, QOI_FORWARD_STATE, QOI_DUMMY_FUNCTION

    ! <<< Arguments >>>
    class(t_Region), intent(in) :: region
    integer, intent(in), optional :: quantityOfInterest, index

    ! <<< Result >>>
    SCALAR_TYPE :: quadrature

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: quantityOfInterest_, index_
    integer :: i, ierror
    SCALAR_TYPE, allocatable :: F(:,:)

    assert(allocated(region%grids))
    assert(allocated(region%states))
    assert(size(region%grids) == size(region%states))

    quantityOfInterest_ = QOI_FORWARD_STATE
    index_ = 0
    if( present(quantityOfInterest) ) quantityOfInterest_ = quantityOfInterest
    if( present(index) ) index_ = index
    if( index_==0 ) quantityOfInterest_ = QOI_DUMMY_FUNCTION

    quadrature = 0.0_wp

    do i = 1, size(region%grids)

       assert(region%grids(i)%nGridPoints > 0)
       allocate(F(region%grids(i)%nGridPoints, 2))
       F(:,2) = 1.0_wp

       select case(quantityOfInterest_)
       case(QOI_FORWARD_STATE)
         assert(allocated(region%states(i)%conservedVariables))
         assert(size(region%states(i)%conservedVariables, 1) == region%grids(i)%nGridPoints)
         assert(size(region%states(i)%conservedVariables, 2) > 2)

         F(:,1) = region%states(i)%conservedVariables(:,index_)
       case(QOI_ADJOINT_STATE)
         assert(allocated(region%states(i)%adjointVariables))
         assert(size(region%states(i)%adjointVariables, 1) == region%grids(i)%nGridPoints)
         assert(size(region%states(i)%adjointVariables, 2) > 2)

         F(:,1) = region%states(i)%adjointVariables(:,index_)
       case(QOI_DUMMY_FUNCTION)
         F(:,1) = 1.0_wp
       end select

       quadrature = quadrature + region%grids(i)%computeInnerProduct(F(:,1),F(:,2))
       SAFE_DEALLOCATE(F)
    end do

    if (region%commGridMasters /= MPI_COMM_NULL)                                               &
         call MPI_Allreduce(MPI_IN_PLACE, quadrature, 1,                                       &
         SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

    do i = 1, size(region%grids)
       call MPI_Bcast(quadrature, 1, SCALAR_TYPE_MPI,                                          &
            0, region%grids(i)%comm, ierror)
    end do

  end function computeRegionIntegral

  function computeAdjointXmomentum(region) result(adjointXmomentum)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Arguments >>>
    class(t_Region), intent(in) :: region

    ! <<< Result >>>
    SCALAR_TYPE :: adjointXmomentum

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, nUnknowns, ierror
    SCALAR_TYPE, allocatable :: F(:,:)

    nUnknowns = region%solverOptions%nUnknowns

    assert(allocated(region%grids))
    assert(allocated(region%states))
    assert(size(region%grids) == size(region%states))

    adjointXmomentum = 0.0_wp

    do i = 1, size(region%grids)

       assert(region%grids(i)%nGridPoints > 0)
       assert(allocated(region%states(i)%adjointVariables))
       assert(size(region%states(i)%adjointVariables, 1) == region%grids(i)%nGridPoints)
       assert(size(region%states(i)%adjointVariables, 2) > 2)

       allocate(F(region%grids(i)%nGridPoints, 2))
       F(:,1) = region%states(i)%adjointVariables(:,nUnknowns)
       F(:,2) = region%states(i)%velocity(:,1)

       adjointXmomentum = adjointXmomentum + region%grids(i)%computeInnerProduct(F(:,1),F(:,2))
       SAFE_DEALLOCATE(F)
    end do

    if (region%commGridMasters /= MPI_COMM_NULL)                                               &
         call MPI_Allreduce(MPI_IN_PLACE, adjointXmomentum, 1,                                 &
         SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

    do i = 1, size(region%grids)
       call MPI_Bcast(adjointXmomentum, 1, SCALAR_TYPE_MPI,                                    &
            0, region%grids(i)%comm, ierror)
    end do

  end function computeAdjointXmomentum

  subroutine addBodyForce(region, mode, stage)

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    use, intrinsic :: iso_fortran_env, only : output_unit
    use ErrorHandler, only : writeAndFlush

    ! <<< Enumerations >>>
    use Region_enum, only : FORWARD, ADJOINT, LINEARIZED
    use State_enum, only : QOI_ADJOINT_STATE, QOI_FORWARD_STATE, QOI_DUMMY_FUNCTION

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: region
    integer, intent(in) :: mode
    integer, intent(in) :: stage

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, nDimensions
    SCALAR_TYPE :: currentXmomentum, adjointForcingFactor
    SCALAR_TYPE, allocatable :: temp(:)

    character(len=STRING_LENGTH) :: message

    call startTiming("addBodyForce")

    assert(allocated(region%states))
    assert(allocated(region%grids))
    nDimensions = region%grids(1)%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    if (stage==1) then
      currentXmomentum = computeRegionIntegral(region,QOI_FORWARD_STATE,2)
      region%momentumLossPerVolume =                                                  &
                    region%oneOverVolume / region%getTimeStepSize() *                 &
                      ( region%initialXmomentum - currentXmomentum )

      if (region%momentumLossPerVolume<0) then
        write(message,'(A,E13.6)') 'Negative body force detected: ', region%momentumLossPerVolume
        call writeAndFlush(region%comm,output_unit,message)
      end if

      if (mode==LINEARIZED) then
        region%adjointMomentumLossPerVolume =                                           &
                  - region%oneOverVolume / region%getTimeStepSize() *                   &
                    computeRegionIntegral(region,QOI_ADJOINT_STATE,2)
      end if
    end if

    select case(mode)
    case(FORWARD)
      do i = 1, size(region%states)
        region%states(i)%rightHandSide(:,2) =                                           &
                region%states(i)%rightHandSide(:,2) + region%momentumLossPerVolume

        region%states(i)%rightHandSide(:,nDimensions+2) =                               &
                             region%states(i)%rightHandSide(:,nDimensions+2) +          &
                     region%momentumLossPerVolume * region%states(i)%velocity(:,1)
      end do

    case(ADJOINT)
      if ( (stage.eq.2) .or. (stage.eq.3) ) then
        adjointForcingFactor = 2.0_wp
      else
        adjointForcingFactor = 1.0_wp
      end if
      region%adjointMomentumLossPerVolume = region%adjointMomentumLossPerVolume           &
        - adjointForcingFactor * ( computeRegionIntegral(region,QOI_ADJOINT_STATE,2)      &
                                    + computeAdjointXmomentum(region) )

      do i = 1, size(region%states)
        allocate(temp(region%grids(i)%nGridPoints))
        temp = region%momentumLossPerVolume * region%states(i)%specificVolume(:,1)                &
                                            * region%states(i)%adjointVariables(:,nDimensions+2)

        !Here the sign is minus (or plus), because magudi adjoint RK4 takes negative time step unnecessarily.
        region%states(i)%rightHandSide(:,2) = region%states(i)%rightHandSide(:,2) - temp
        region%states(i)%rightHandSide(:,1) = region%states(i)%rightHandSide(:,1) + temp          &
                                                                * region%states(i)%velocity(:,1)
        SAFE_DEALLOCATE(temp)
      end do

      if (stage.eq.1) then
        do i = 1, size(region%states)
          !Here the sign is minus, because magudi adjoint RK4 takes negative time step unnecessarily.
          region%states(i)%rightHandSide(:,2) = region%states(i)%rightHandSide(:,2)       &
           - region%adjointMomentumLossPerVolume * region%oneOverVolume / region%getTimeStepSize()
        end do
        region%adjointMomentumLossPerVolume = 0.0_wp
      end if

    case(LINEARIZED)
      do i = 1, size(region%states)
        allocate(temp(region%grids(i)%nGridPoints))
        temp = - region%states(i)%velocity(:,1) * region%states(i)%adjointVariables(:,1)&
               + region%states(i)%adjointVariables(:,2)
        temp = temp * region%states(i)%specificVolume(:,1)

        region%states(i)%rightHandSide(:,2) =                                           &
            region%states(i)%rightHandSide(:,2) + region%adjointMomentumLossPerVolume

        region%states(i)%rightHandSide(:,nDimensions+2) =                               &
            region%states(i)%rightHandSide(:,nDimensions+2)                             &
            + region%adjointMomentumLossPerVolume * region%states(i)%velocity(:,1)      &
            + region%momentumLossPerVolume * temp

        SAFE_DEALLOCATE(temp)
      end do

    end select

    call endTiming("addBodyForce")

  end subroutine addBodyForce

  subroutine computeBulkQuantities(region, bulkConservedVariables)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Arguments >>>
    class(t_Region), intent(in) :: region

    ! <<< Result >>>
    SCALAR_TYPE, intent(out) :: bulkConservedVariables(:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, ierror, nDimensions
    SCALAR_TYPE, allocatable :: F(:,:)

    assert(allocated(region%grids))
    assert(allocated(region%states))
    assert(size(region%grids) == size(region%states))

    nDimensions = size(region%globalGridSizes, 1)
    assert_key(nDimensions, (1, 2, 3))
    assert(size(bulkConservedVariables)==nDimensions+2)

    bulkConservedVariables = 0.0_wp

    do i = 1, size(region%grids)

       assert(region%grids(i)%nGridPoints > 0)
       assert(allocated(region%states(i)%conservedVariables))
       assert(size(region%states(i)%conservedVariables, 1) == region%grids(i)%nGridPoints)
       assert(size(region%states(i)%conservedVariables, 2) > 2)

       allocate(F(region%grids(i)%nGridPoints, nDimensions+3))
       F(:,1:nDimensions+2) = region%states(i)%conservedVariables
       F(:,nDimensions+3) = 1.0_wp
       do j = 1, nDimensions+2
         bulkConservedVariables(j) = bulkConservedVariables(j)                                 &
                            + region%grids(i)%computeInnerProduct(F(:,j),F(:,nDimensions+3))
       end do
       SAFE_DEALLOCATE(F)
    end do

    if (region%commGridMasters /= MPI_COMM_NULL)                                               &
         call MPI_Allreduce(MPI_IN_PLACE, bulkConservedVariables, nDimensions+2,               &
         SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

    do i = 1, size(region%grids)
       call MPI_Bcast(bulkConservedVariables, nDimensions+2, SCALAR_TYPE_MPI,                  &
            0, region%grids(i)%comm, ierror)
    end do

  end subroutine computeBulkQuantities

  subroutine normalizeKolmogorovForcing(region)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env, only : output_unit

    ! <<< Derived types >>>
    use Region_mod, only : t_Region
    use Patch_mod, only : t_Patch
    use KolmogorovForcingPatch_mod, only : t_KolmogorovForcingPatch

    ! <<< Internal modules >>>
    use Patch_factory, only : computeQuadratureOnPatches
    use ErrorHandler, only : writeAndFlush

    ! <<< Arguments >>>
    class(t_Region) :: region

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: ierror, j, k
    class(t_Patch), pointer :: patch => null()
    real(SCALAR_KIND), allocatable :: patchVolume(:), gridVolume(:), gridForcing(:)
    real(SCALAR_KIND) :: globalVolume, globalForcing
    character(len=STRING_LENGTH) :: message

    globalVolume = 0.0_wp
    globalForcing = 0.0_wp

    do j = 1, size(region%grids)
      allocate(gridVolume(region%grids(j)%nGridPoints))
      allocate(gridForcing(region%grids(j)%nGridPoints))
      gridVolume = 0.0_wp
      gridForcing = 0.0_wp
      do k = 1, size(region%patchFactories)
        call region%patchFactories(k)%connect(patch)
        if (.not. associated(patch)) cycle
        if (patch%gridIndex /= region%grids(j)%index .or. patch%nPatchPoints <= 0) cycle
        select type (patch)
        class is (t_KolmogorovForcingPatch)
          call patch%disperseAdd(patch%forcePerUnitMass, gridForcing)
          allocate(patchVolume,MOLD=patch%forcePerUnitMass)
          patchVolume = 1.0_wp
          call patch%disperseAdd(patchVolume, gridVolume)
          SAFE_DEALLOCATE(patchVolume)
        end select
      end do

      globalVolume = globalVolume +                                                                     &
                         computeQuadratureOnPatches(region%patchFactories, 'KOLMOGOROV_FORCING',        &
                                                     region%grids(j), gridVolume )
      globalForcing = globalForcing +                                                                   &
                        computeQuadratureOnPatches(region%patchFactories, 'KOLMOGOROV_FORCING',         &
                                                    region%grids(j), gridForcing )

      SAFE_DEALLOCATE(gridVolume)
      SAFE_DEALLOCATE(gridForcing)
    end do

    if (region%commGridMasters /= MPI_COMM_NULL) then
         call MPI_Allreduce(MPI_IN_PLACE, globalVolume, 1,                                          &
         SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)
         call MPI_Allreduce(MPI_IN_PLACE, globalForcing, 1,                                          &
         SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)
    end if

    do j = 1, size(region%grids)
       call MPI_Bcast(globalForcing, 1, SCALAR_TYPE_MPI,                                            &
            0, region%grids(j)%comm, ierror)
       call MPI_Bcast(globalVolume, 1, SCALAR_TYPE_MPI,                                            &
            0, region%grids(j)%comm, ierror)
    end do

    write(message,'(2(A,'//SCALAR_FORMAT//'))')                                                     &
    'Kolmogorov forcing offset: ', globalForcing, ', patch volume: ', globalVolume
    call writeAndFlush(region%comm, output_unit, message)

    do j = 1, size(region%grids)
      do k = 1, size(region%patchFactories)
        call region%patchFactories(k)%connect(patch)
        if (.not. associated(patch)) cycle
        if (patch%gridIndex /= region%grids(j)%index .or. patch%nPatchPoints <= 0) cycle
        select type (patch)
        class is (t_KolmogorovForcingPatch)
          patch%forcePerUnitMass = patch%forcePerUnitMass - globalForcing / globalVolume
        end select
      end do
    end do

  end subroutine normalizeKolmogorovForcing

end module RegionImpl

subroutine setupRegion(this, comm, globalGridSizes, simulationFlags, solverOptions, verbose)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Private members >>>
  use RegionImpl, only : readDecompositionMap, distributeGrids

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

  call startTiming("setupRegion")

  ! Clean slate.
  call this%cleanup()
  this%comm = comm
  call MPI_Comm_size(this%comm, numProcs, ierror)

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
       numProcs > size(this%globalGridSizes, 2)) then
     call getRequiredOption("decomposition_map_file",                                        &
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

  if (this%simulationFlags%enableIBM) then
    call this%levelsetFactory%cleanup()
    nullify(this%levelsetFactory)
  end if

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
  SAFE_DEALLOCATE(this%localToGlobalPatchIndex)

  if (this%comm /= MPI_COMM_NULL .and. this%comm /= MPI_COMM_WORLD)                          &
       call MPI_Comm_free(this%comm, ierror)
  this%comm = MPI_COMM_NULL

  if (this%commGridMasters /= MPI_COMM_NULL) call MPI_Comm_free(this%commGridMasters, ierror)
  this%commGridMasters = MPI_COMM_NULL

  this%timestep = 0
  this%outputOn = .true.

end subroutine cleanupRegion

subroutine setupBoundaryConditions(this, boundaryConditionFilename)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use PatchDescriptor_mod, only : t_PatchDescriptor

  ! <<< Private members >>>
  use RegionImpl

  ! <<< Internal modules >>>
  use InterfaceHelper, only : readPatchInterfaceInformation
  use MPITimingsHelper, only : startTiming, endTiming
  use Patch_factory, only : queryPatchTypeExists

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: this
  character(len = *), intent(in) :: boundaryConditionFilename

  ! <<< Local variables >>>
  integer :: i, j, k, nPatches, ierror
  type(t_PatchDescriptor) :: p
  class(t_Patch), pointer :: patch => null()

  call startTiming("setupBoundaryConditions")

  call readBoundaryConditions(this, boundaryConditionFilename)
  call readPatchInterfaceInformation(this)
  call validatePatches(this)
  call distributePatches(this)
  call MPI_Barrier(this%comm, ierror)

  nPatches = 0
  if (allocated(this%patchCommunicators))                                                    &
       nPatches = count(this%patchCommunicators /= MPI_COMM_NULL)
  if (nPatches > 0) then
     allocate(this%patchFactories(nPatches))
     allocate(this%localToGlobalPatchIndex(nPatches))
     this%localToGlobalPatchIndex = 0
  end if

  if (allocated(this%patchData)) then
     do k = 1, size(this%grids)
        do i = 1, size(this%patchData)
           p = this%patchData(i)
           if (p%gridIndex /= this%grids(k)%index) cycle
           if (allocated(this%patchFactories)) then
              do j = 1, size(this%patchFactories)
                 if (this%localToGlobalPatchIndex(j) /= 0 .or.                               &
                      this%patchCommunicators(i) == MPI_COMM_NULL) cycle
                 call this%patchFactories(j)%connect(patch, trim(p%patchType))
                 assert(associated(patch))
                 call patch%setup(i, this%patchCommunicators(i), p, this%grids(k),           &
                      this%simulationFlags, this%solverOptions)
                 this%localToGlobalPatchIndex(j) = i
                 exit
              end do
           end if
           call MPI_Barrier(this%grids(k)%comm, ierror)
        end do
     end do
  end if

  if (queryPatchTypeExists(this%patchFactories,'KOLMOGOROV_FORCING'))                        &
     call normalizeKolmogorovForcing(this)

  if (this%simulationFlags%enableController)                                                 &
     call normalizeControlMollifier(this)
  if (this%simulationFlags%enableFunctional)                                                 &
     call normalizeTargetMollifier(this)

  call endTiming("setupBoundaryConditions")

end subroutine setupBoundaryConditions

subroutine loadRegionData(this, quantityOfInterest, filename, speciesFilename)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use PLOT3DDescriptor_type, only : t_PLOT3DDescriptor, PLOT3D_SOLUTION_FILE

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
  real(wp) :: auxiliaryData(4)

  call startTiming("loadRegionData")

  isSolutionFile = .false.
  speciesFilename_ = ""
  if (present(speciesFilename)) speciesFilename_ = speciesFilename
  speciesFileOffset = int(0, MPI_OFFSET_KIND)

  select case(quantityOfInterest)
  case (QOI_GRID, QOI_JACOBIAN, QOI_TARGET_MOLLIFIER, QOI_CONTROL_MOLLIFIER,                 &
       QOI_METRICS, QOI_DUMMY_FUNCTION)
  case default
     isSolutionFile = .true.
     if (this%solverOptions%nSpecies > 0 .and. len_trim(speciesFilename_) == 0 .and.         &
          filename(len_trim(filename)-1:len_trim(filename)) /= ".q") then
        write(message, '(3A)') "Auto-detection of species filename failed: Solution file '", &
             trim(filename), "' does not have a '.q' extension!"
        call gracefulExit(this%comm, message)
     end if
     if (this%solverOptions%nSpecies > 0)                                                    &
          speciesFilename_ = filename(:len_trim(filename)-2) // ".f"
  end select

  if (present(speciesFilename)) then
     write(message, '(5A)') "Reading '", trim(filename), "', '",                             &
          trim(speciesFilename), "'..."
  else if (len_trim(speciesFilename_) > 0) then
     write(message, '(3A)') "Reading '", filename(:len_trim(filename)-2), ".q+f'..."
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
              speciesFileOffset = plot3dGetOffset(this%gridCommunicators(i),                 &
                   speciesFilename_, i, success)
              if (.not. success) exit
           end if

           select case(quantityOfInterest)
           case (QOI_GRID, QOI_JACOBIAN, QOI_METRICS, QOI_TARGET_MOLLIFIER,                  &
                QOI_CONTROL_MOLLIFIER)
              call this%grids(j)%loadData(quantityOfInterest,                                &
                   trim(filename), offset, success)
           case default
              isSolutionFile = (quantityOfInterest /= QOI_DUMMY_FUNCTION)
              call this%states(j)%loadData(this%grids(j), quantityOfInterest,                &
                   trim(filename), offset, success,                                          &
                   speciesFilename = trim(speciesFilename_),                                 &
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
        auxiliaryData = real(this%states(1)%plot3dAuxiliaryData, wp)
        call MPI_Bcast(auxiliaryData, 4, REAL_TYPE_MPI, 0, this%comm, ierror)
        do i = 1, 4
           this%states(:)%plot3dAuxiliaryData(i) = auxiliaryData(i)
        end do
        if (quantityOfInterest == QOI_FORWARD_STATE) then
           this%states(:)%time = real(auxiliaryData(4), wp)
           this%timestep = nint(real(auxiliaryData(1), wp))
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
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: message, speciesFilename
  logical :: success
  integer :: i, j, nScalars, errorRank, procRank, ierror
  integer(kind = MPI_OFFSET_KIND) :: offset, speciesFileOffset
  logical :: isSolutionFile

  if (.not. this%OutputOn) return

  call startTiming("saveRegionData")

  isSolutionFile = .false.
  speciesFilename = ""
  speciesFileOffset = int(0, MPI_OFFSET_KIND)

  select case(quantityOfInterest)
  case (QOI_GRID, QOI_JACOBIAN, QOI_TARGET_MOLLIFIER, QOI_CONTROL_MOLLIFIER,                 &
       QOI_METRICS, QOI_DUMMY_FUNCTION, QOI_NORM)
  case default
     isSolutionFile = .true.
     if (filename(len_trim(filename)-1:len_trim(filename)) /= ".q") then
        write(message, '(A)') "Solution files must have extension '.q'!"
        call gracefulExit(this%comm, message)
     end if
     if (this%solverOptions%nSpecies > 0)                                                    &
          speciesFilename = filename(:len_trim(filename)-2) // ".f"
  end select

  if (len_trim(speciesFilename) > 0) then
     write(message, '(3A)') "Writing '", filename(:len_trim(filename)-2), ".q+f'..."
  else
     write(message, '(3A)') "Writing '", trim(filename), "'..."
  end if
  call writeAndFlush(this%comm, output_unit, message, advance = 'no')

  select case(quantityOfInterest)

  case (QOI_GRID)
     call plot3dWriteSkeleton(this%comm, trim(filename),                                     &
          PLOT3D_GRID_FILE, this%globalGridSizes, success)
  case (QOI_JACOBIAN, QOI_TARGET_MOLLIFIER, QOI_CONTROL_MOLLIFIER, QOI_NORM)
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
     if (len_trim(speciesFilename) > 0)                                                      &
          call plot3dWriteSkeleton(this%comm, trim(speciesFilename), PLOT3D_FUNCTION_FILE,   &
          this%globalGridSizes, success, this%solverOptions%nSpecies)

  end select

  if (isSolutionFile) then
     if (.not. this%simulationFlags%steadyStateSimulation) then
        do i = 1, size(this%states)
           this%states(:)%plot3dAuxiliaryData(1) =                                           &
                real(this%timestep, wp)
        end do
     else
        do i = 1, size(this%states)
           this%states(:)%plot3dAuxiliaryData(1) = 0.0_wp
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
              speciesFileOffset = plot3dGetOffset(this%gridCommunicators(i),                 &
                   speciesFilename, i, success)
              if (.not. success) exit
           end if

           select case(quantityOfInterest)
           case (QOI_GRID, QOI_JACOBIAN, QOI_METRICS, QOI_TARGET_MOLLIFIER,                  &
                QOI_CONTROL_MOLLIFIER, QOI_NORM)
              call this%grids(j)%saveData(quantityOfInterest,                                &
                   trim(filename), offset, success)
           case default
              call this%states(j)%saveData(this%grids(j), quantityOfInterest,                &
                   trim(filename), offset, success, speciesFilename = trim(speciesFilename), &
                   speciesFileOffset = speciesFileOffset)
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

subroutine computeRhs(this, mode, timeStep, stage)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT, LINEARIZED

  ! <<< Internal modules >>>
  use RhsHelper, only : computeRhsForward, computeRhsAdjoint,                                 &
                        computeRhsLinearized, addInterfaceAdjointPenalty
  use InterfaceHelper, only : exchangeInterfaceData
  use MPITimingsHelper, only : startTiming, endTiming
  use RegionImpl, only : addBodyForce

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: this
  integer, intent(in) :: mode
  integer, intent(in) :: timeStep, stage

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j
  class(t_Patch), pointer :: patch => null()

  call startTiming("computeRhs")

  ! Semi-discrete right-hand-side operator.
  do i = 1, size(this%states)
     select case (mode)
     case (FORWARD)
        call computeRhsForward(this%simulationFlags, this%solverOptions,                     &
             this%grids(i), this%states(i), this%patchFactories)
     case (ADJOINT)
        call computeRhsAdjoint(this%simulationFlags, this%solverOptions,                     &
             this%grids(i), this%states(i), this%patchFactories)
     case (LINEARIZED)
        call computeRhsLinearized(this%simulationFlags, this%solverOptions,                  &
             this%grids(i), this%states(i), this%patchFactories)
     end select
  end do

  ! Exchange data at block interfaces.
  if (allocated(this%patchFactories)) then
     do i = 1, size(this%patchFactories)
        call this%patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        do j = 1, size(this%states)
           if (patch%gridIndex /= this%grids(j)%index) cycle
           select type (patch)
           class is (t_BlockInterfacePatch)
              call patch%collectInterfaceData(mode, this%simulationFlags,                    &
                   this%solverOptions, this%grids(j), this%states(j))
           end select
        end do
     end do
  end if

  call exchangeInterfaceData(this)

  ! Disperse received data at block interfaces.
  if (allocated(this%patchFactories)) then
     do i = 1, size(this%patchFactories)
        call this%patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        do j = 1, size(this%states)
           if (patch%gridIndex /= this%grids(j)%index) cycle
           select type (patch)
           class is (t_BlockInterfacePatch)
              call patch%disperseInterfaceData(mode, this%simulationFlags,                   &
                   this%solverOptions)
           end select
        end do
     end do
  end if

  ! Viscous adjoint penalties at block interfaces. NOTE!!: currently this is only discrete adjoint.
  if (mode == ADJOINT .and. this%simulationFlags%viscosityOn) then
     do i = 1, size(this%states)
        call addInterfaceAdjointPenalty(this%simulationFlags, this%solverOptions,            &
             this%grids(i), this%states(i), this%patchFactories)
     end do
  end if

  ! Multiply by Jacobian.
  do i = 1, size(this%states)
     do j = 1, this%solverOptions%nUnknowns
        this%states(i)%rightHandSide(:,j) = this%states(i)%rightHandSide(:,j) *              &
                this%grids(i)%jacobian(:,1)
     end do
  end do

  ! Update immersed boundary variables.
  if (this%simulationFlags%enableIBM) then
    call this%levelsetFactory%updateLevelset(mode, this%grids, this%states)
    do i = 1, size(this%states)
       call this%states(i)%updateIBMVariables(mode, this%grids(i), this%simulationFlags)
    end do
  end if

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

  ! Acoustic source terms.
  do i = 1, size(this%states)
     call this%states(i)%addSources(mode, this%grids(i))
  end do

  ! x-momentum conserving body force. ONLY FOR RK4
  if (this%simulationFlags%enableBodyForce)                                                  &
    call addBodyForce(this,mode,stage)

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

subroutine saveSpongeStrength(this, filename)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use SpongePatch_mod, only : t_SpongePatch

  ! <<< Enumerations >>>
  use State_enum, only : QOI_DUMMY_FUNCTION

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: this
  character(len = *), intent(in) :: filename

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  type :: t_SpongeStrengthInternal
     SCALAR_TYPE, allocatable :: buffer1(:)
     SCALAR_TYPE, pointer :: buffer2(:,:) => null()
  end type t_SpongeStrengthInternal
  type(t_SpongeStrengthInternal), allocatable :: data(:)
  integer :: i, j, ierror
  class(t_Patch), pointer :: patch => null()

  allocate(data(size(this%grids)))

  do i = 1, size(this%grids)
     allocate(data(i)%buffer1(this%grids(i)%nGridPoints))
     allocate(data(i)%buffer2(this%grids(i)%nGridPoints, 1))
     data(i)%buffer2 = 0.0_wp

     if (allocated(this%patchFactories)) then
        do j = 1, size(this%patchFactories)
           call this%patchFactories(j)%connect(patch)

           if (.not. associated(patch)) cycle
           if (patch%gridIndex /= this%grids(i)%index .or. patch%nPatchPoints <= 0) cycle

           select type (patch)
              class is (t_SpongePatch)
              call patch%disperse(patch%spongeStrength, data(i)%buffer1)
              data(i)%buffer2(:,1) = data(i)%buffer2(:,1) + data(i)%buffer1
           end select

        end do
     end if

    call MPI_Allreduce(MPI_IN_PLACE, data(i)%buffer2, size(data(i)%buffer2),                 &
         SCALAR_TYPE_MPI, MPI_SUM, this%grids(i)%comm, ierror)
     this%states(i)%dummyFunction => data(i)%buffer2
  end do

  call this%saveData(QOI_DUMMY_FUNCTION, trim(filename))

  do i = 1, size(data)
     SAFE_DEALLOCATE(data(i)%buffer1)
     if (associated(data(i)%buffer2)) deallocate(data(i)%buffer2)
     nullify(data(i)%buffer2)
  end do
  SAFE_DEALLOCATE(data)

end subroutine saveSpongeStrength

subroutine resetProbes(this)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ProbePatch_mod, only : t_ProbePatch

  ! <<< Internal modules >>>
  use InputHelper, only : getFreeUnit

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: this

  ! <<< Local variables >>>
  integer :: i, stat, fileUnit, procRank, ierror
  class(t_Patch), pointer :: patch => null()

  if (allocated(this%patchFactories)) then
     do i = 1, size(this%patchFactories)
        call this%patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        if (patch%comm == MPI_COMM_NULL) cycle
        select type (patch)
        class is (t_ProbePatch)

           call MPI_Comm_rank(patch%comm, procRank, ierror)
           if (procRank == 0) then
              open(unit = getFreeUnit(fileUnit), file = trim(patch%probeFilename),           &
                   iostat = stat, status = 'old')
              if (stat == 0) close(fileUnit, status = 'delete')
              open(unit = getFreeUnit(fileUnit), file = trim(patch%probeFilename),           &
                   action = 'write', status = 'unknown')
              close(fileUnit)
           end if
           call MPI_Barrier(patch%comm, ierror)
           patch%iProbeBuffer = 0
           patch%probeFileOffset = int(0, MPI_OFFSET_KIND)

        end select
     end do
  end if

end subroutine resetProbes

subroutine saveProbeData(this, mode, finish)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ProbePatch_mod, only : t_ProbePatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: this
  integer, intent(in) :: mode
  logical, intent(in), optional :: finish

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  logical :: finish_
  integer :: i, j
  class(t_Patch), pointer :: patch => null()

  finish_ = .false.
  if (present(finish)) finish_ = finish

  if (.not. allocated(this%patchFactories)) return

  do i = 1, size(this%patchFactories)
     call this%patchFactories(i)%connect(patch)
     if (.not. associated(patch)) cycle
     do j = 1, size(this%states)
        if (patch%gridIndex /= this%grids(j)%index .or. patch%nPatchPoints <= 0) cycle
        select type (patch)
        class is (t_ProbePatch)

           if (finish_) then
              call patch%saveData()
              patch%iProbeBuffer = 0
              cycle
           end if

           patch%iProbeBuffer = patch%iProbeBuffer + 1
           assert(patch%iProbeBuffer >= 1)
           assert(patch%iProbeBuffer <= size(patch%probeBuffer, 3))

           select case (mode)
           case (FORWARD)
              call patch%collect(this%states(j)%conservedVariables,                          &
                   patch%probeBuffer(:,:,patch%iProbeBuffer))
           case (ADJOINT)
              call patch%collect(this%states(j)%adjointVariables,                            &
                   patch%probeBuffer(:,:,patch%iProbeBuffer))
           end select

           if (patch%iProbeBuffer == size(patch%probeBuffer, 3)) then
              call patch%saveData()
              patch%iProbeBuffer = 0
           end if

        end select
     end do
  end do

end subroutine saveProbeData

subroutine connectLevelsetFactory(this)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Patch_mod, only : t_Patch
  use ImmersedBoundaryPatch_mod, only : t_ImmersedBoundaryPatch
  use LevelsetFactory_mod, only : t_LevelsetFactory
  use SinusoidalWallLevelset_mod, only : t_SinusoidalWallLevelset
  use StokesSecondWallLevelset_mod, only : t_StokesSecondWallLevelset

  ! <<< Internal modules >>>
  use InputHelper, only : getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: this

  ! <<< Local variables >>>
  character(len = STRING_LENGTH) :: levelsetType
  integer :: i, j
  class(t_Patch), pointer :: patch => null()

  ! Determine whether IBM is actually used in each state.
  do i = 1, size(this%states)
    this%states(i)%ibmPatchExists = .false.
  end do

  if (allocated(this%patchFactories)) then
    do i = 1, size(this%states)
      do j = 1, size(this%patchFactories)
        call this%patchFactories(j)%connect(patch)
        if (.not. associated(patch)) cycle
        if (patch%gridIndex /= this%grids(i)%index) cycle

        select type (patch)
        class is (t_ImmersedBoundaryPatch)
          this%states(i)%ibmPatchExists = .true.
          exit
        end select
      end do ! j = 1, size(this%patchFactories)
    end do ! i = 1, size(this%states)
  end if

  ! Determine levelset type.
  call getRequiredOption("immersed_boundary/levelset_type", levelsetType, this%comm)

  this%levelsetType = levelsetType

  select case (trim(levelsetType))

  case ('sinusoidal_wall')
    allocate(t_SinusoidalWallLevelset :: this%levelsetFactory)

  case ('stokes_second_wall')
    allocate(t_StokesSecondWallLevelset :: this%levelsetFactory)

  case default
    this%levelsetType = ""

  end select

end subroutine connectLevelsetFactory
