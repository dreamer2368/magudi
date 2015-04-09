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
  use InputHelper, only : getOption, getRequiredOption

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
  integer :: ierror, offset_(3)

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
  this%comm = comm
  this%nDimensions = grid%nDimensions

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
    call MPI_Comm_size(this%comm, nProcs, ierror)
    allocate(this%mpiReduceBuffer(nProcs))
#endif

  ! Derived types describing subarrays on patch.
  if (this%comm /= MPI_COMM_NULL) then
     offset_ = this%offset - this%extent(1::2) + 1
     call MPI_Type_create_subarray(3, this%globalSize, this%localSize, offset_,              &
          MPI_ORDER_FORTRAN, SCALAR_TYPE_MPI, this%mpiDerivedTypeScalarSubarray, ierror)
     call MPI_Type_commit(this%mpiDerivedTypeScalarSubarray, ierror)
     call MPI_Type_create_subarray(3, this%globalSize, this%localSize, offset_,              &
          MPI_ORDER_FORTRAN, MPI_INTEGER, this%mpiDerivedTypeIntegerSubarray, ierror)
     call MPI_Type_commit(this%mpiDerivedTypeIntegerSubarray, ierror)
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
  integer :: ierror

  if (this%mpiDerivedTypeIntegerSubarray /= MPI_DATATYPE_NULL)                               &
       call MPI_Type_free(this%mpiDerivedTypeIntegerSubarray, ierror)
  if (this%mpiDerivedTypeScalarSubarray /= MPI_DATATYPE_NULL)                                &
       call MPI_Type_free(this%mpiDerivedTypeScalarSubarray, ierror)
  if (this%comm /= MPI_COMM_NULL .and. this%comm /= MPI_COMM_WORLD)                          &
       call MPI_Comm_free(this%comm, ierror)

  this%comm = MPI_COMM_NULL
  this%mpiDerivedTypeScalarSubarray = MPI_DATATYPE_NULL
  this%mpiDerivedTypeIntegerSubarray = MPI_DATATYPE_NULL

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

subroutine gatherScalarOnPatch(this, patchLocalArray, patchGlobalArray)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch

  implicit none

  ! <<< Arguments >>>
  class(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: patchLocalArray(:)
  SCALAR_TYPE, intent(out), allocatable :: patchGlobalArray(:)

  ! <<< Local variables >>>
  integer :: i, procRank, numProcs, ierror
  integer, allocatable :: mpiRequest(:), mpiDerivedTypeScalarSubarray(:)

  if (this%comm == MPI_COMM_NULL) return

  call MPI_Comm_rank(this%comm, procRank, ierror)
  call MPI_Comm_size(this%comm, numProcs, ierror)

  if (procRank == 0)                                                                         &
       allocate(patchGlobalArray(product(this%globalSize)))
  allocate(mpiRequest(numProcs), source = MPI_REQUEST_NULL)
  allocate(mpiDerivedTypeScalarSubarray(numProcs), source = MPI_DATATYPE_NULL)

  call MPI_Allgather(this%mpiDerivedTypeScalarSubarray, 1, MPI_INTEGER,                      &
       mpiDerivedTypeScalarSubarray, 1, MPI_INTEGER, this%comm, ierror)

  if (procRank == 0) then
     do i = 0, numProcs - 1
        call MPI_Irecv(patchGlobalArray, 1, mpiDerivedTypeScalarSubarray(i+1),               &
             i, i, this%comm, mpiRequest(i+1), ierror)
     end do
  end if

  call MPI_Isend(patchLocalArray, size(patchLocalArray), SCALAR_TYPE_MPI, 0, procRank,       &
       this%comm, mpiRequest(procRank + 1), ierror)
  call MPI_Waitall(numProcs, mpiRequest, MPI_STATUSES_IGNORE, ierror)

  SAFE_DEALLOCATE(mpiDerivedTypeScalarSubarray)
  SAFE_DEALLOCATE(mpiRequest)

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
  SCALAR_TYPE, intent(out), allocatable :: patchGlobalArray(:,:)

  ! <<< Local variables >>>
  integer :: i, procRank, numProcs, ierror
  integer, allocatable :: mpiRequest(:), mpiDerivedTypeScalarSubarray(:)

  if (this%comm == MPI_COMM_NULL) return

  call MPI_Comm_rank(this%comm, procRank, ierror)
  call MPI_Comm_size(this%comm, numProcs, ierror)

  if (procRank == 0)                                                                         &
       allocate(patchGlobalArray(product(this%globalSize), size(patchLocalArray, 2)))
  allocate(mpiRequest(numProcs), source = MPI_REQUEST_NULL)
  allocate(mpiDerivedTypeScalarSubarray(numProcs), source = MPI_DATATYPE_NULL)

  call MPI_Allgather(this%mpiDerivedTypeScalarSubarray, 1, MPI_INTEGER,                      &
       mpiDerivedTypeScalarSubarray, 1, MPI_INTEGER, this%comm, ierror)

  if (procRank == 0) then
     do i = 0, numProcs - 1
        call MPI_Irecv(patchGlobalArray, size(patchGlobalArray, 2),                          &
             mpiDerivedTypeScalarSubarray(i+1), i, i, this%comm, mpiRequest(i+1), ierror)
     end do
  end if

  call MPI_Isend(patchLocalArray, size(patchLocalArray), SCALAR_TYPE_MPI, 0, procRank,       &
       this%comm, mpiRequest(procRank + 1), ierror)
  call MPI_Waitall(numProcs, mpiRequest, MPI_STATUSES_IGNORE, ierror)

  SAFE_DEALLOCATE(mpiDerivedTypeScalarSubarray)
  SAFE_DEALLOCATE(mpiRequest)

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
  SCALAR_TYPE, intent(out), allocatable :: patchGlobalArray(:,:,:)

  ! <<< Local variables >>>
  integer :: i, nComponents, procRank, numProcs, ierror
  integer, allocatable :: mpiRequest(:), mpiDerivedTypeScalarSubarray(:)

  if (this%comm == MPI_COMM_NULL) return

  call MPI_Comm_rank(this%comm, procRank, ierror)
  call MPI_Comm_size(this%comm, numProcs, ierror)

  if (procRank == 0)                                                                         &
       allocate(patchGlobalArray(product(this%globalSize), size(patchLocalArray, 2),         &
       size(patchLocalArray, 3)))
  allocate(mpiRequest(numProcs), source = MPI_REQUEST_NULL)
  allocate(mpiDerivedTypeScalarSubarray(numProcs), source = MPI_DATATYPE_NULL)

  call MPI_Allgather(this%mpiDerivedTypeScalarSubarray, 1, MPI_INTEGER,                      &
       mpiDerivedTypeScalarSubarray, 1, MPI_INTEGER, this%comm, ierror)

  if (procRank == 0) then
     nComponents = size(patchGlobalArray, 2) * size(patchGlobalArray, 3)
     do i = 0, numProcs - 1
        call MPI_Irecv(patchGlobalArray, nComponents, mpiDerivedTypeScalarSubarray(i+1),     &
             i, i, this%comm, mpiRequest(i+1), ierror)
     end do
  end if

  call MPI_Isend(patchLocalArray, size(patchLocalArray), SCALAR_TYPE_MPI, 0, procRank,       &
       this%comm, mpiRequest(procRank + 1), ierror)
  call MPI_Waitall(numProcs, mpiRequest, MPI_STATUSES_IGNORE, ierror)

  SAFE_DEALLOCATE(mpiDerivedTypeScalarSubarray)
  SAFE_DEALLOCATE(mpiRequest)

end subroutine gatherTensorOnPatch
