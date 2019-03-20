#include "config.h"

subroutine testComputeRhs(this,mode)

  ! <<< External modules >>>

  ! <<< Derived types >>>
  use testRegion_mod

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  implicit none

  ! <<< Arguemnts >>>
  class(t_testRegion) :: this
  integer, intent(in) :: mode

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i

  do i = 1, size(this%states)
    select case (mode)
    case(FORWARD)
      this%states(i)%rightHandSide = this%tempRhs(this%tempRhsIndex)
      this%tempRhsIndex = this%tempRhsIndex + 1
    case(ADJOINT)
      if(this%tempRhsIndex==1) then
        this%states(i)%rightHandSide = 0.0_wp
      else
        this%states(i)%rightHandSide = - this%states(i)%adjointForcingFactor        &
                                          * this%tempRhs(this%tempRhsIndex-1)
      end if
      this%tempRhsIndex = this%tempRhsIndex - 1
    end select
  end do
end subroutine testComputeRhs
