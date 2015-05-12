#include "config.h"

subroutine setupPatch(this, index, comm, patchDescriptor,                                    &
     grid, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use Patch_mod, only : t_Patch
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_Patch) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  character(len = STRING_LENGTH) :: key
  integer :: i, offset(3), procRank, numProcs, ierror
  integer, allocatable :: allLocalSizes(:,:), allOffsets(:,:)

  assert(index > 0)

  assert_key(grid%nDimensions, (1, 2, 3))

  assert(all(grid%offset >= 0))
  assert(all(grid%localSize > 0))

  assert(patchDescriptor%gridIndex == grid%index)
  assert(abs(patchDescriptor%normalDirection) <= grid%nDimensions)
  assert(patchDescriptor%iMin > 0)
  assert(patchDescriptor%jMin > 0)
  assert(patchDescriptor%kMin > 0)
  assert(patchDescriptor%iMax >= patchDescriptor%iMin)
  assert(patchDescriptor%jMax >= patchDescriptor%jMin)
  assert(patchDescriptor%kMax >= patchDescriptor%kMin)
  assert(len_trim(patchDescriptor%name) > 0)

  call this%cleanupBase()

  this%index = index
  this%name = trim(patchDescriptor%name)
  if (comm /= MPI_COMM_NULL) call MPI_Comm_dup(comm, this%comm, ierror)

  this%globalSize(1) = patchDescriptor%iMax - patchDescriptor%iMin + 1
  this%globalSize(2) = patchDescriptor%jMax - patchDescriptor%jMin + 1
  this%globalSize(3) = patchDescriptor%kMax - patchDescriptor%kMin + 1

  ! Zero-based index of first point on the patch belonging to the ``current'' process (this
  ! value has no meaning if the patch lies outside the part of the grid belonging to the
  ! ``current'' process).
  this%offset(1) = max(patchDescriptor%iMin, grid%offset(1) + 1)
  this%offset(2) = max(patchDescriptor%jMin, grid%offset(2) + 1)
  this%offset(3) = max(patchDescriptor%kMin, grid%offset(3) + 1)
  this%offset = min(this%offset, grid%offset + grid%localSize) - 1

  ! Extent of the patch belonging to the ``current'' process (this value has no meaning if the
  ! patch lies outside the part of the grid belonging to the ``current'' process).
  this%localSize(1) = max(patchDescriptor%iMax, grid%offset(1) + 1)
  this%localSize(2) = max(patchDescriptor%jMax, grid%offset(2) + 1)
  this%localSize(3) = max(patchDescriptor%kMax, grid%offset(3) + 1)
  this%localSize = min(this%localSize, grid%offset + grid%localSize)
  this%localSize = this%localSize - this%offset

#ifdef DEBUG
  if (comm /= MPI_COMM_NULL) then
     assert(all(this%localSize > 0))
  end if
#endif

  ! Reset size and offset if the patch lies outside the part of the grid belonging to the
  ! ``current'' process.
  if (any(this%localSize < 0)) then
     this%offset = 0
     this%localSize = 0
  end if

  this%nPatchPoints = product(this%localSize)
  this%gridIndex = grid%index
  this%normalDirection = patchDescriptor%normalDirection

  this%extent = (/ patchDescriptor%iMin, patchDescriptor%iMax, patchDescriptor%jMin,         &
       patchDescriptor%jMax, patchDescriptor%kMin, patchDescriptor%kMax /)

  this%gridLocalSize = grid%localSize
  this%gridOffset = grid%offset

  this%isCurvilinear = getOption("default/curvilinear", .true.)
  write(key, '(A,I3.3,A)') "grid", this%gridIndex, "/curvilinear"
  this%isCurvilinear = getOption(key, this%isCurvilinear)

#ifdef SCALAR_TYPE_IS_binary128_IEEE754
  if (this%comm /= MPI_COMM_NULL) then
    call MPI_Comm_size(this%comm, numProcs, ierror)
    allocate(this%mpiReduceBuffer(numProcs))
 end if
#endif

 if (this%comm /= MPI_COMM_NULL) then

    call MPI_Comm_rank(this%comm, procRank, ierror)
    call MPI_Comm_size(this%comm, numProcs, ierror)

    ! Patch derived subarray type.
    offset = this%offset - this%extent(1::2) + 1
    call MPI_Type_create_subarray(3, this%globalSize, this%localSize, offset,                &
         MPI_ORDER_FORTRAN, SCALAR_TYPE_MPI, this%mpiScalarSubarrayType, ierror)
    call MPI_Type_commit(this%mpiScalarSubarrayType, ierror)

    ! Gather local size and offset from all patch processes.
    allocate(allLocalSizes(3, numProcs), source = 0)
    allocate(allOffsets(3, numProcs), source = 0)
    call MPI_Gather(this%localSize, 3, MPI_INTEGER, allLocalSizes,                           &
         3, MPI_INTEGER, 0, this%comm, ierror)
    call MPI_Gather(this%offset, 3, MPI_INTEGER, allOffsets,                                 &
         3, MPI_INTEGER, 0, this%comm, ierror)

    if (procRank == 0) then

       ! Allocate patch subarray derived types for all processes on patch master process.
       allocate(this%mpiAllScalarSubarrayTypes(numProcs))

       ! Derived types describing subarrays on patches.
       do i = 1, numProcs
          allOffsets(:,i) = allOffsets(:,i) - this%extent(1::2) + 1
          call MPI_Type_create_subarray(3, this%globalSize, allLocalSizes(:,i),              &
               allOffsets(:,i), MPI_ORDER_FORTRAN, SCALAR_TYPE_MPI,                          &
               this%mpiAllScalarSubarrayTypes(i), ierror)
          call MPI_Type_commit(this%mpiAllScalarSubarrayTypes(i), ierror)
       end do

    end if

    SAFE_DEALLOCATE(allOffsets)
    SAFE_DEALLOCATE(allLocalSizes)

 end if

end subroutine setupPatch

subroutine cleanupPatch(this)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch

  implicit none

  ! <<< Arguments >>>
  class(t_Patch) :: this

  ! <<< Local variables >>>
  integer :: i, ierror

  if (allocated(this%mpiAllScalarSubarrayTypes)) then
     do i = 1, size(this%mpiAllScalarSubarrayTypes)
        if (this%mpiAllScalarSubarrayTypes(i) /= MPI_DATATYPE_NULL) &
             call MPI_Type_free(this%mpiAllScalarSubarrayTypes(i), ierror)
     end do
  end if
  SAFE_DEALLOCATE(this%mpiAllScalarSubarrayTypes)

  if (this%mpiScalarSubarrayType /= MPI_DATATYPE_NULL)                                       &
       call MPI_Type_free(this%mpiScalarSubarrayType, ierror)
  this%mpiScalarSubarrayType = MPI_DATATYPE_NULL

  if (this%comm /= MPI_COMM_NULL .and. this%comm /= MPI_COMM_WORLD)                          &
       call MPI_Comm_free(this%comm, ierror)
  this%comm = MPI_COMM_NULL

end subroutine cleanupPatch

subroutine collectScalarAtPatch(this, gridArray, patchArray)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  ! <<< Arguments >>>
  class(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: gridArray(:)
  SCALAR_TYPE, intent(out) :: patchArray(:)

  ! <<< Local variables >>>
  integer :: i, j, k, patchIndex, localIndex

  assert(all(this%gridLocalSize > 0) .and. size(gridArray, 1) == product(this%gridLocalSize))
  assert(all(this%localSize >= 0) .and. size(patchArray, 1) == product(this%localSize))

  call startTiming("collectAtPatch")

  do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
     do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
        do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
           patchIndex = i - this%offset(1) +                                                 &
                this%localSize(1) * (j - 1 - this%offset(2) +                                &
                this%localSize(2) * (k - 1 - this%offset(3)))
           localIndex = i - this%gridOffset(1) +                                             &
                this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                        &
                this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
           patchArray(patchIndex) = gridArray(localIndex)
        end do
     end do
  end do

  call endTiming("collectAtPatch")

end subroutine collectScalarAtPatch

subroutine collectVectorAtPatch(this, gridArray, patchArray)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  ! <<< Arguments >>>
  class(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: gridArray(:,:)
  SCALAR_TYPE, intent(out) :: patchArray(:,:)

  ! <<< Local variables >>>
  integer :: i, j, k, l, patchIndex, localIndex

  assert(all(this%gridLocalSize > 0) .and. size(gridArray, 1) == product(this%gridLocalSize))
  assert(all(this%localSize >= 0) .and. size(patchArray, 1) == product(this%localSize))
  assert(size(gridArray, 2) > 0)
  assert(size(patchArray, 2) == size(gridArray, 2))

  call startTiming("collectAtPatch")

  do l = 1, size(patchArray, 2)
     do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
        do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
           do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
              patchIndex = i - this%offset(1) +                                              &
                   this%localSize(1) * (j - 1 - this%offset(2) +                             &
                   this%localSize(2) * (k - 1 - this%offset(3)))
              localIndex = i - this%gridOffset(1) +                                          &
                   this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                     &
                   this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
              patchArray(patchIndex,l) = gridArray(localIndex,l)
           end do
        end do
     end do
  end do

  call endTiming("collectAtPatch")

end subroutine collectVectorAtPatch

subroutine collectTensorAtPatch(this, gridArray, patchArray)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  ! <<< Arguments >>>
  class(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: gridArray(:,:,:)
  SCALAR_TYPE, intent(out) :: patchArray(:,:,:)

  ! <<< Local variables >>>
  integer :: i, j, k, l, m, patchIndex, localIndex

  assert(all(this%gridLocalSize > 0) .and. size(gridArray, 1) == product(this%gridLocalSize))
  assert(all(this%localSize >= 0) .and. size(patchArray, 1) == product(this%localSize))
  assert(size(gridArray, 2) > 0)
  assert(size(patchArray, 2) == size(gridArray, 2))
  assert(size(gridArray, 3) > 0)
  assert(size(patchArray, 3) == size(gridArray, 3))

  call startTiming("collectAtPatch")

  do m = 1, size(patchArray, 3)
     do l = 1, size(patchArray, 2)
        do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
           do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
              do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                 patchIndex = i - this%offset(1) +                                           &
                      this%localSize(1) * (j - 1 - this%offset(2) +                          &
                      this%localSize(2) * (k - 1 - this%offset(3)))
                 localIndex = i - this%gridOffset(1) +                                       &
                      this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                  &
                      this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
                 patchArray(patchIndex,l,m) = gridArray(localIndex,l,m)
              end do
           end do
        end do
     end do
  end do

  call endTiming("collectAtPatch")

end subroutine collectTensorAtPatch

subroutine disperseScalarFromPatch(this, patchArray, gridArray)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  ! <<< Arguments >>>
  class(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: patchArray(:)
  SCALAR_TYPE, intent(out) :: gridArray(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, patchIndex, localIndex

  assert(all(this%localSize >= 0) .and. size(patchArray, 1) == product(this%localSize))
  assert(all(this%gridLocalSize > 0) .and. size(gridArray, 1) == product(this%gridLocalSize))

  call startTiming("disperseFromPatch")

  gridArray = 0.0_wp

  do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
     do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
        do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
           patchIndex = i - this%offset(1) +                                                 &
                this%localSize(1) * (j - 1 - this%offset(2) +                                &
                this%localSize(2) * (k - 1 - this%offset(3)))
           localIndex = i - this%gridOffset(1) +                                             &
                this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                        &
                this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
           gridArray(localIndex) = patchArray(patchIndex)
        end do
     end do
  end do

  call endTiming("disperseFromPatch")

end subroutine disperseScalarFromPatch

subroutine disperseVectorFromPatch(this, patchArray, gridArray)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  ! <<< Arguments >>>
  class(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: patchArray(:,:)
  SCALAR_TYPE, intent(out) :: gridArray(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, patchIndex, localIndex

  assert(all(this%localSize >= 0) .and. size(patchArray, 1) == product(this%localSize))
  assert(all(this%gridLocalSize > 0) .and. size(gridArray, 1) == product(this%gridLocalSize))
  assert(size(patchArray, 2) > 0)
  assert(size(gridArray, 2) == size(patchArray, 2))

  call startTiming("disperseFromPatch")

  gridArray = 0.0_wp

  do l = 1, size(patchArray, 2)
     do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
        do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
           do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
              patchIndex = i - this%offset(1) +                                              &
                   this%localSize(1) * (j - 1 - this%offset(2) +                             &
                   this%localSize(2) * (k - 1 - this%offset(3)))
              localIndex = i - this%gridOffset(1) +                                          &
                   this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                     &
                   this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
              gridArray(localIndex,l) = patchArray(patchIndex,l)
           end do
        end do
     end do
  end do

  call endTiming("disperseFromPatch")

end subroutine disperseVectorFromPatch

subroutine disperseTensorFromPatch(this, patchArray, gridArray)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  ! <<< Arguments >>>
  class(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: patchArray(:,:,:)
  SCALAR_TYPE, intent(out) :: gridArray(:,:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, m, patchIndex, localIndex

  assert(all(this%localSize >= 0) .and. size(patchArray, 1) == product(this%localSize))
  assert(all(this%gridLocalSize > 0) .and. size(gridArray, 1) == product(this%gridLocalSize))
  assert(size(patchArray, 2) > 0)
  assert(size(gridArray, 2) == size(patchArray, 2))
  assert(size(patchArray, 3) > 0)
  assert(size(gridArray, 3) == size(patchArray, 3))

  call startTiming("disperseFromPatch")

  gridArray = 0.0_wp

  do m = 1, size(patchArray, 3)
     do l = 1, size(patchArray, 2)
        do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
           do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
              do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                 patchIndex = i - this%offset(1) +                                           &
                      this%localSize(1) * (j - 1 - this%offset(2) +                          &
                      this%localSize(2) * (k - 1 - this%offset(3)))
                 localIndex = i - this%gridOffset(1) +                                       &
                      this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                  &
                      this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
                 gridArray(localIndex,l,m) = patchArray(patchIndex,l,m)
              end do
           end do
        end do
     end do
  end do

  call endTiming("disperseFromPatch")

end subroutine disperseTensorFromPatch

subroutine gatherScalarOnPatch(this, patchLocalArray, patchGlobalArray)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch

  implicit none

  ! <<< Arguments >>>
  class(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: patchLocalArray(:)
  SCALAR_TYPE, allocatable :: patchGlobalArray(:)

  ! <<< Local variables >>>
  integer :: i, mpiRequest, procRank, numProcs, ierror

  if (this%comm == MPI_COMM_NULL) return

  assert(this%nPatchPoints > 0)
  assert(size(patchLocalArray) == this%nPatchPoints)

  call MPI_Comm_rank(this%comm, procRank, ierror)
  call MPI_Comm_size(this%comm, numProcs, ierror)

#ifdef DEBUG
  if (procRank == 0) then
     assert(allocated(patchGlobalArray))
     assert(size(patchGlobalArray) == product(this%globalSize))
  end if
#endif

  call MPI_Isend(patchLocalArray, size(patchLocalArray), SCALAR_TYPE_MPI,                    &
       0, procRank, this%comm, mpiRequest, ierror)

  if (procRank == 0) then
     do i = 0, numProcs - 1
        call MPI_Recv(patchGlobalArray, 1, this%mpiAllScalarSubarrayTypes(i+1),             &
             i, i, this%comm, MPI_STATUS_IGNORE, ierror)
     end do
  end if

  call MPI_Wait(mpiRequest, MPI_STATUS_IGNORE, ierror)

end subroutine gatherScalarOnPatch

subroutine gatherVectorOnPatch(this, patchLocalArray, patchGlobalArray)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch

  implicit none

  ! <<< Arguments >>>
  class(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: patchLocalArray(:,:)
  SCALAR_TYPE, allocatable :: patchGlobalArray(:,:)

  ! <<< Local variables >>>
  integer :: i, mpiRequest, procRank, numProcs, ierror

  if (this%comm == MPI_COMM_NULL) return

  assert(this%nPatchPoints > 0)
  assert(size(patchLocalArray, 1) == this%nPatchPoints)
  assert(size(patchLocalArray, 2) > 0)

  call MPI_Comm_rank(this%comm, procRank, ierror)
  call MPI_Comm_size(this%comm, numProcs, ierror)

#ifdef DEBUG
  if (procRank == 0) then
     assert(allocated(patchGlobalArray))
     assert(size(patchGlobalArray, 1) == product(this%globalSize))
     assert(size(patchGlobalArray, 2) == size(patchLocalArray, 2))
  end if
#endif

  call MPI_Isend(patchLocalArray, size(patchLocalArray), SCALAR_TYPE_MPI, 0, procRank,       &
       this%comm, mpiRequest, ierror)

  if (procRank == 0) then
     do i = 0, numProcs - 1
        call MPI_Recv(patchGlobalArray, size(patchGlobalArray, 2),                           &
             this%mpiAllScalarSubarrayTypes(i+1), i, i, this%comm, MPI_STATUS_IGNORE, ierror)
     end do
  end if

  call MPI_Wait(mpiRequest, MPI_STATUS_IGNORE, ierror)

end subroutine gatherVectorOnPatch

subroutine gatherTensorOnPatch(this, patchLocalArray, patchGlobalArray)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch

  implicit none

  ! <<< Arguments >>>
  class(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: patchLocalArray(:,:,:)
  SCALAR_TYPE, allocatable :: patchGlobalArray(:,:,:)

  ! <<< Local variables >>>
  integer :: i, nComponents, mpiRequest, procRank, numProcs, ierror

  if (this%comm == MPI_COMM_NULL) return

  assert(this%nPatchPoints > 0)
  assert(size(patchLocalArray, 1) == this%nPatchPoints)
  assert(size(patchLocalArray, 2) > 0)
  assert(size(patchLocalArray, 3) > 0)

  call MPI_Comm_rank(this%comm, procRank, ierror)
  call MPI_Comm_size(this%comm, numProcs, ierror)

#ifdef DEBUG
  if (procRank == 0) then
     assert(allocated(patchGlobalArray))
     assert(size(patchGlobalArray, 1) == product(this%globalSize))
     assert(size(patchGlobalArray, 2) == size(patchLocalArray, 2))
     assert(size(patchGlobalArray, 3) == size(patchLocalArray, 3))
  end if
#endif

  call MPI_Isend(patchLocalArray, size(patchLocalArray), SCALAR_TYPE_MPI,                    &
       0, procRank, this%comm, mpiRequest, ierror)

  if (procRank == 0) then
     nComponents = size(patchGlobalArray, 2) * size(patchGlobalArray, 3)
     do i = 0, numProcs - 1
        call MPI_Recv(patchGlobalArray, nComponents, this%mpiAllScalarSubarrayTypes(i+1),    &
             i, i, this%comm, MPI_STATUS_IGNORE, ierror)
     end do
  end if

  call MPI_Wait(mpiRequest, MPI_STATUS_IGNORE, ierror)

end subroutine gatherTensorOnPatch

subroutine scatterScalarOnPatch(this, patchGlobalArray, patchLocalArray)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch

  implicit none

  ! <<< Arguments >>>
  class(t_Patch) :: this
  SCALAR_TYPE, intent(in), allocatable :: patchGlobalArray(:)
  SCALAR_TYPE, intent(out) :: patchLocalArray(:)

  ! <<< Local variables >>>
  integer :: i, mpiRequest, procRank, numProcs, ierror

  if (this%comm == MPI_COMM_NULL) return

  assert(this%nPatchPoints > 0)
  assert(size(patchLocalArray) == this%nPatchPoints)

  call MPI_Comm_rank(this%comm, procRank, ierror)
  call MPI_Comm_size(this%comm, numProcs, ierror)

#ifdef DEBUG
  if (procRank == 0) then
     assert(allocated(patchGlobalArray))
     assert(size(patchGlobalArray) == product(this%globalSize))
  end if
#endif

  call MPI_Irecv(patchLocalArray, size(patchLocalArray), SCALAR_TYPE_MPI,                    &
       0, procRank, this%comm, mpiRequest, ierror)

  if (procRank == 0) then
     do i = 0, numProcs - 1
        call MPI_Send(patchGlobalArray, 1, this%mpiAllScalarSubarrayTypes(i+1),              &
             i, i, this%comm, ierror)
     end do
  end if

  call MPI_Wait(mpiRequest, MPI_STATUS_IGNORE, ierror)

end subroutine scatterScalarOnPatch

subroutine scatterVectorOnPatch(this, patchGlobalArray, patchLocalArray)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch

  implicit none

  ! <<< Arguments >>>
  class(t_Patch) :: this
  SCALAR_TYPE, intent(in), allocatable :: patchGlobalArray(:,:)
  SCALAR_TYPE, intent(out) :: patchLocalArray(:,:)

  ! <<< Local variables >>>
  integer :: i, mpiRequest, procRank, numProcs, ierror

  if (this%comm == MPI_COMM_NULL) return

  assert(this%nPatchPoints > 0)
  assert(size(patchLocalArray, 1) == this%nPatchPoints)
  assert(size(patchLocalArray, 2) > 0)

  call MPI_Comm_rank(this%comm, procRank, ierror)
  call MPI_Comm_size(this%comm, numProcs, ierror)

#ifdef DEBUG
  if (procRank == 0) then
     assert(allocated(patchGlobalArray))
     assert(size(patchGlobalArray, 1) == product(this%globalSize))
     assert(size(patchGlobalArray, 2) == size(patchLocalArray, 2))
  end if
#endif

  call MPI_Irecv(patchLocalArray, size(patchLocalArray), SCALAR_TYPE_MPI,                    &
       0, procRank, this%comm, mpiRequest, ierror)

  if (procRank == 0) then
     do i = 0, numProcs - 1
        call MPI_Send(patchGlobalArray, size(patchGlobalArray, 2),                          &
             this%mpiAllScalarSubarrayTypes(i+1), i, i, this%comm, ierror)
     end do
  end if

  call MPI_Wait(mpiRequest, MPI_STATUS_IGNORE, ierror)

end subroutine scatterVectorOnPatch

subroutine scatterTensorOnPatch(this, patchGlobalArray, patchLocalArray)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch

  implicit none

  ! <<< Arguments >>>
  class(t_Patch) :: this
  SCALAR_TYPE, intent(in), allocatable :: patchGlobalArray(:,:,:)
  SCALAR_TYPE, intent(out) :: patchLocalArray(:,:,:)

  ! <<< Local variables >>>
  integer :: i, nComponents, mpiRequest, procRank, numProcs, ierror

  if (this%comm == MPI_COMM_NULL) return

  assert(this%nPatchPoints > 0)
  assert(size(patchLocalArray, 1) == this%nPatchPoints)
  assert(size(patchLocalArray, 2) > 0)

  call MPI_Comm_rank(this%comm, procRank, ierror)
  call MPI_Comm_size(this%comm, numProcs, ierror)

#ifdef DEBUG
  if (procRank == 0) then
     assert(allocated(patchGlobalArray))
     assert(size(patchGlobalArray, 1) == product(this%globalSize))
     assert(size(patchGlobalArray, 2) == size(patchLocalArray, 2))
     assert(size(patchGlobalArray, 3) == size(patchLocalArray, 3))
  end if
#endif

  call MPI_Irecv(patchLocalArray, size(patchLocalArray), SCALAR_TYPE_MPI,                    &
       0, procRank, this%comm, mpiRequest, ierror)

  if (procRank == 0) then
     nComponents = size(patchGlobalArray, 2) * size(patchGlobalArray, 3)
     do i = 0, numProcs - 1
        call MPI_Send(patchGlobalArray, nComponents, this%mpiAllScalarSubarrayTypes(i+1),   &
             i, i, this%comm, ierror)
     end do
  end if

  call MPI_Wait(mpiRequest, MPI_STATUS_IGNORE, ierror)

end subroutine scatterTensorOnPatch
