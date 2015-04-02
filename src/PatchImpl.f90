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
  this%normalDirection = patchDescriptor%normalDirection
  this%gridIndex = grid%index
  this%nDimensions = grid%nDimensions

  this%extent = (/ patchDescriptor%iMin, patchDescriptor%iMax,                               &
       patchDescriptor%jMin, patchDescriptor%jMax,                                           &
       patchDescriptor%kMin, patchDescriptor%kMax /)

  ! Zero-based index of first point on the patch belonging to the ``current'' process (this
  ! value has no meaning if the patch lies outside the part of the grid belonging to the
  ! ``current'' process).
  this%offset(1) = max(patchDescriptor%iMin, grid%offset(1) + 1)
  this%offset(2) = max(patchDescriptor%jMin, grid%offset(2) + 1)
  this%offset(3) = max(patchDescriptor%kMin, grid%offset(3) + 1)
  this%offset = min(this%offset, grid%offset + grid%localSize) - 1

  ! Extent of the patch belonging to the ``current'' process (this value has no meaning if the
  ! patch lies outside the part of the grid belonging to the ``current'' process).
  this%patchSize(1) = max(patchDescriptor%iMax, grid%offset(1) + 1)
  this%patchSize(2) = max(patchDescriptor%jMax, grid%offset(2) + 1)
  this%patchSize(3) = max(patchDescriptor%kMax, grid%offset(3) + 1)
  this%patchSize = min(this%patchSize, grid%offset + grid%localSize)
  this%patchSize = this%patchSize - this%offset

  ! Reset size and offset if the patch lies outside the part of the grid belonging to the
  ! ``current'' process.
  if (any(this%patchSize < 0)) then
     this%offset = 0
     this%patchSize = 0
  end if

  this%gridLocalSize = grid%localSize
  this%gridOffset = grid%offset

  this%nPatchPoints = product(this%patchSize)
  this%comm = comm

  this%isCurvilinear = getOption("default/curvilinear", .true.)
  write(key, '(A,I3.3,A)') "grid", this%gridIndex, "/curvilinear"
  this%isCurvilinear = getOption(key, this%isCurvilinear)

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

  if (this%comm /= MPI_COMM_NULL) call MPI_Comm_free(this%comm, ierror)
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
  assert(all(this%patchSize >= 0) .and. size(patchArray, 1) == product(this%patchSize)) 

  call startTiming("collectAtPatch")

  do k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
     do j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
        do i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
           patchIndex = i - this%offset(1) +                                                 &
                this%patchSize(1) * (j - 1 - this%offset(2) +                                &
                this%patchSize(2) * (k - 1 - this%offset(3)))
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
  assert(all(this%patchSize >= 0) .and. size(patchArray, 1) == product(this%patchSize))
  assert(size(gridArray, 2) > 0)
  assert(size(patchArray, 2) == size(gridArray, 2))

  call startTiming("collectAtPatch")

  do l = 1, size(patchArray, 2)
     do k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
        do j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
           do i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
              patchIndex = i - this%offset(1) +                                              &
                   this%patchSize(1) * (j - 1 - this%offset(2) +                             &
                   this%patchSize(2) * (k - 1 - this%offset(3)))
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
  assert(all(this%patchSize >= 0) .and. size(patchArray, 1) == product(this%patchSize))
  assert(size(gridArray, 2) > 0)
  assert(size(patchArray, 2) == size(gridArray, 2))
  assert(size(gridArray, 3) > 0)
  assert(size(patchArray, 3) == size(gridArray, 3))

  call startTiming("collectAtPatch")

  do m = 1, size(patchArray, 3)
     do l = 1, size(patchArray, 2)
        do k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
           do j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
              do i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
                 patchIndex = i - this%offset(1) +                                           &
                      this%patchSize(1) * (j - 1 - this%offset(2) +                          &
                      this%patchSize(2) * (k - 1 - this%offset(3)))
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
