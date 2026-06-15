#include "config.h"

function normQuadratic(this, L2sqSum) result(val)

  ! <<< Derived types >>>
  use MSPenalty_mod, only : t_QuadraticPenalty

  implicit none

  ! <<< Arguments >>>
  class(t_QuadraticPenalty), intent(in) :: this
  SCALAR_TYPE, intent(in) :: L2sqSum

  ! <<< Result >>>
  SCALAR_TYPE :: val

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  val = 0.5_wp * L2sqSum

end function normQuadratic

function gradientFactorQuadratic(this, L2sqSum) result(val)

  ! <<< Derived types >>>
  use MSPenalty_mod, only : t_QuadraticPenalty

  implicit none

  ! <<< Arguments >>>
  class(t_QuadraticPenalty), intent(in) :: this
  SCALAR_TYPE, intent(in) :: L2sqSum

  ! <<< Result >>>
  SCALAR_TYPE :: val

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  val = 1.0_wp

end function gradientFactorQuadratic

function normHuber(this, L2sqSum) result(val)

  ! <<< Derived types >>>
  use MSPenalty_mod, only : t_HuberPenalty

  implicit none

  ! <<< Arguments >>>
  class(t_HuberPenalty), intent(in) :: this
  SCALAR_TYPE, intent(in) :: L2sqSum

  ! <<< Result >>>
  SCALAR_TYPE :: val

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp) :: L2norm

  L2norm = sqrt(real(L2sqSum, wp))
  if (L2norm >= this%threshold) then
     val = L2norm - 0.5_wp * this%threshold
  else
     val = 0.5_wp * L2sqSum / this%threshold
  end if

end function normHuber

function gradientFactorHuber(this, L2sqSum) result(val)

  ! <<< Derived types >>>
  use MSPenalty_mod, only : t_HuberPenalty

  implicit none

  ! <<< Arguments >>>
  class(t_HuberPenalty), intent(in) :: this
  SCALAR_TYPE, intent(in) :: L2sqSum

  ! <<< Result >>>
  SCALAR_TYPE :: val

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp) :: L2norm

  L2norm = sqrt(real(L2sqSum, wp))
  if (L2norm >= this%threshold) then
     val = L2norm
  else
     val = this%threshold
  end if

end function gradientFactorHuber

subroutine connectMSPenalty(penaltyPtr, penaltyType, huberThreshold)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use MSPenalty_mod, only : t_MSPenalty, t_QuadraticPenalty, t_HuberPenalty

  ! <<< Internal modules >>>
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_MSPenalty), pointer, intent(out) :: penaltyPtr
  character(len = *), intent(in) :: penaltyType
  real(SCALAR_KIND), intent(in), optional :: huberThreshold

  ! <<< Local variables >>>
  character(len = STRING_LENGTH) :: message

  nullify(penaltyPtr)

  select case (trim(penaltyType))

  case ('quadratic')
     allocate(t_QuadraticPenalty :: penaltyPtr)

  case ('huber')
     allocate(t_HuberPenalty :: penaltyPtr)
     if (present(huberThreshold)) then
        select type (penaltyPtr)
        type is (t_HuberPenalty)
           penaltyPtr%threshold = huberThreshold
        end select
     end if

  case default
     write(message, '(3A)') "Unknown time_splitting/penalty_type '", trim(penaltyType),     &
          "': expected 'quadratic' or 'huber'."
     call gracefulExit(MPI_COMM_WORLD, message)

  end select

end subroutine connectMSPenalty
