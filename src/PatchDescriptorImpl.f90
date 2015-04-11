#include "config.h"

subroutine validatePatchDescriptor(this, globalGridSizes,                                    &
     simulationFlags, solverOptions, errorCode, message)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Patch_factory, only : t_PatchFactory
  use Functional_mod, only : t_Functional
  use SolverOptions_mod, only : t_SolverOptions
  use Functional_factory, only : t_FunctionalFactory
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  ! <<< Arguments >>>
  class(t_PatchDescriptor) :: this
  integer, intent(in) :: globalGridSizes(:,:)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  integer, intent(out) :: errorCode
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Local variables >>>
  integer :: i, extent(6)
  logical :: isExtentValid, flag, success
  type(t_PatchFactory) :: patchFactory
  class(t_Patch), pointer :: dummyPatch => null()
  character(len = STRING_LENGTH) :: str
  type(t_FunctionalFactory) :: functionalFactory
  class(t_Functional), pointer :: dummyFunctional => null()

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
     write(message, '(3A,I0.0,A,6(I0.0,1X),A)') "Patch '", trim(this%name), "' on grid ",    &
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
  flag = dummyPatch%verifyUsage(this, globalGridSizes(:,this%gridIndex),                     &
       this%normalDirection, extent, simulationFlags, success, str)

  if (.not. success) then
     write(message, '(3A,I0.0,2A)') "Patch '", trim(this%name), "' on grid ",                &
          this%gridIndex, " reported: ", trim(str)
     errorCode = 2
     return
  end if

  if (.not. flag) then
     write(message, '(3A,I0.0,A)') "Patch '", trim(this%name), "' on grid ",                 &
          this%gridIndex, " will not be used!"
     errorCode = 1
     return
  end if

  select type (dummyPatch)
  class is (t_CostTargetPatch)

     call functionalFactory%connect(dummyFunctional, trim(solverOptions%costFunctionalType))
     if (associated(dummyFunctional)) then
        success = dummyFunctional%isPatchValid(this, globalGridSizes(:,this%gridIndex),      &
             this%normalDirection, extent, simulationFlags, str)
        if (.not. success) then
           write(message, '(3A,I0.0,2A)') "Patch '", trim(this%name), "' on grid ",          &
                this%gridIndex, " reported: ", trim(str)
           errorCode = 2
           return
        end if
     end if

  end select

  call patchFactory%cleanup()

end subroutine validatePatchDescriptor
