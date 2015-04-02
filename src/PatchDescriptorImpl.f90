#include "config.h"

subroutine validatePatchDescriptor(this, globalGridSizes,                                    &
     simulationFlags, errorCode, message)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Patch_factory, only : t_PatchFactory
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  ! <<< Arguments >>>
  class(t_PatchDescriptor) :: this
  integer, intent(in) :: globalGridSizes(:,:)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  integer, intent(out) :: errorCode
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Local variables >>>
  type(t_PatchFactory) :: patchFactory
  class(t_Patch), pointer :: dummyPatch => null()
  integer :: i, extent(6)
  logical :: isExtentValid, flag, success

  assert(len_trim(this%patchType) > 0)

  errorCode = 0

  extent = (/ this%iMin, this%iMax, this%jMin, this%jMax, this%kMin, this%kMax /)

  ! Check if the grid index is valid.
  if (this%gridIndex <= 0 .or. this%gridIndex > size(globalGridSizes, 2)) then
     write(message, '(3A,I0.0,A)') "Patch '", trim(this%name),                               &
          "' has an invalid grid index: ", this%gridIndex, "!"
     errorCode = 2
     return
  end if

  ! Extents must be non-zero.
  isExtentValid = all(extent /= 0)

  ! Negative values indicate counting backwards from the end.
  do i = 1, size(globalGridSizes, 1)
     if (extent(1+2*(i-1)) < 0) extent(1+2*(i-1)) = extent(1+2*(i-1))                        &
          + globalGridSizes(i, this%gridIndex) + 1
     if (extent(2+2*(i-1)) < 0) extent(2+2*(i-1)) = extent(2+2*(i-1))                        &
          + globalGridSizes(i, this%gridIndex) + 1
  end do

  ! Check that extent describes a part of the grid.
  do i = 1, size(globalGridSizes, 1)
     isExtentValid = isExtentValid .and. (extent(2+2*(i-1)) >= extent(1+2*(i-1)))
     isExtentValid = isExtentValid .and. (extent(2+2*(i-1))                                  &
          - extent(1+2*(i-1)) + 1 <= globalGridSizes(i, this%gridIndex))
  end do
  do while (i <= 3) !... reset for direction > number of dimensions.
     extent(1+2*(i-1)) = 1
     extent(2+2*(i-1)) = 1
     i = i + 1
  end do

  ! Fail if the extent is invalid.
  if (.not. isExtentValid) then
     write(message, '(3A,I0.0,A,6(I0.0,1X))') "Patch '", trim(this%name), "' on grid ",      &
          this%gridIndex, " has an invalid extent: ",                                        &
          this%iMin, this%iMax, this%jMin, this%jMax, this%kMin, this%kMax, "!"
     errorCode = 2
     return
  end if

  ! Copy over extent to iMin, iMax, etc.
  this%iMin = extent(1)
  this%iMax = extent(2)
  this%jMin = extent(3)
  this%jMax = extent(4)
  this%kMin = extent(5)
  this%kMax = extent(6)

  call patchFactory%connect(dummyPatch, trim(this%patchType))
  flag = dummyPatch%verifyUsage(success, message)

  if (.not. success) then
     errorCode = 2
     return
  end if

  if (.not. flag) then
     write(message, '(3A,I0.0,A)') "Patch '", trim(this%name), "' on grid ",                 &
          this%gridIndex, " will not be used!"
     errorCode = 1
     return
  end if

  call patchFactory%cleanup()

end subroutine validatePatchDescriptor
