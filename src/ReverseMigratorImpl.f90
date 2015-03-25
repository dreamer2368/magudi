#include "config.h"

subroutine setupReverseMigrator(this, algorithm, region, nStages, numIntermediateStates)

  ! <<< Derived types >>>
  use Region_type
  use ReverseMigrator_type

  ! <<< Public members >>>
  use ReverseMigrator_mod, only : cleanupReverseMigrator

  implicit none

  ! <<< Arguments >>>
  type(t_ReverseMigrator) :: this
  integer, intent(in) :: algorithm
  type(t_Region), intent(in) :: region
  integer, intent(in) :: nStages, numIntermediateStates

  ! <<< Local variables >>>
  integer :: i

  assert(nStages > 0)
  assert(numIntermediateStates >= 0)
  assert(mod(numIntermediateStates, nStages) == 0)

  call cleanupReverseMigrator(this)
  
  this%algorithm = algorithm
  this%nStages = nStages
  this%numIntermediateStates = numIntermediateStates

  allocate(this%temp(size(region%states)))
  do i = 1, size(this%temp)
     assert(region%grids(i)%nGridPoints > 0)
     assert(region%states(i)%nUnknowns > 0)
     allocate(this%temp(i)%buffer(region%grids(i)%nGridPoints, &
          region%states(i)%nUnknowns, numIntermediateStates))
  end do

end subroutine setupReverseMigrator

subroutine cleanupReverseMigrator(this)

  ! <<< Derived types >>>
  use ReverseMigrator_type

  implicit none

  ! <<< Arguments >>>
  type(t_ReverseMigrator) :: this

  ! <<< Local variables >>>
  integer :: i

  if (allocated(this%temp)) then
     do i = 1, size(this%temp)
        SAFE_DEALLOCATE(this%temp(i)%buffer)
     end do
  end if

  SAFE_DEALLOCATE(this%temp)

end subroutine cleanupReverseMigrator
