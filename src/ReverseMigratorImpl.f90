#include "config.h"

subroutine setupReverseMigrator(this, algorithm, region, nStages, numIntermediateStates)

  ! <<< Derived types >>>
  use Region_type, only : t_Region
  use ReverseMigrator_type, only : t_ReverseMigrator, UNIFORM_CHECKPOINTING

  ! <<< Public members >>>
  use ReverseMigrator_mod, only : cleanupReverseMigrator

  implicit none

  ! <<< Arguments >>>
  type(t_ReverseMigrator) :: this
  character(len = *), intent(in) :: algorithm
  type(t_Region), intent(in) :: region
  integer, intent(in) :: nStages, numIntermediateStates

  ! <<< Local variables >>>
  integer :: i

  assert(nStages > 0)
  assert(numIntermediateStates >= 0)
  assert(mod(numIntermediateStates, nStages) == 0)

  call cleanupReverseMigrator(this)

  if (trim(algorithm) == "uniform checkpointing") then
     this%algorithm = UNIFORM_CHECKPOINTING
  end if

  this%nStages = nStages
  this%numIntermediateStates = numIntermediateStates

  assert(allocated(region%states))
  assert(size(region%states) > 0)

  allocate(this%temp_(size(region%states)))
  do i = 1, size(this%temp_)
     assert(region%grids(i)%nGridPoints > 0)
     assert(region%states(i)%nUnknowns > 0)
     allocate(this%temp_(i)%buffer(region%grids(i)%nGridPoints, &
          region%states(i)%nUnknowns, numIntermediateStates))
  end do

end subroutine setupReverseMigrator

subroutine cleanupReverseMigrator(this)

  ! <<< Derived types >>>
  use ReverseMigrator_type, only : t_ReverseMigrator

  implicit none

  ! <<< Arguments >>>
  type(t_ReverseMigrator) :: this

  ! <<< Local variables >>>
  integer :: i

  if (allocated(this%temp_)) then
     do i = 1, size(this%temp_)
        SAFE_DEALLOCATE(this%temp_(i)%buffer)
     end do
  end if

  SAFE_DEALLOCATE(this%temp_)

end subroutine cleanupReverseMigrator

subroutine migrateToTimestep(this, region, timestep, outputPrefix)

  ! <<< Derived types >>>
  use Region_type, only : t_Region
  use ReverseMigrator_type, only : t_ReverseMigrator

  implicit none

  ! <<< Arguments >>>
  type(t_ReverseMigrator) :: this
  type(t_Region) :: region
  integer, intent(in) :: timestep
  character(len = *), intent(in) :: outputPrefix

  assert(timestep >= 0)

end subroutine migrateToTimestep
