#include "config.h"

subroutine testComputeRhs(this,mode, timeStep, stage)

  ! <<< External modules >>>

  ! <<< Derived types >>>
  use testRegion_mod
  use Patch_mod, only : t_Patch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  implicit none

  ! <<< Arguemnts >>>
  class(t_testRegion) :: this
  integer, intent(in) :: mode
  integer, intent(in), optional :: timeStep, stage

  ! <<< Local derived types >>>
  class(t_Patch), pointer :: patch => null()

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j

  if (.not.(present(timeStep) .and. present(stage))) then
    print *, 'testRhs: timestep and stage are not provided!'
    STOP -1
  end if

  do i = 1, size(this%states)
    this%states(i)%rightHandSide = 0.0_wp
  end do

  do i = 1, size(this%states)
    select case (mode)
    case(FORWARD)
      this%tempRhsIndex = stage + 4*(timeStep-1)
      this%states(i)%rightHandSide = this%tempRhs(this%tempRhsIndex)
    end select
  end do

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
end subroutine testComputeRhs
