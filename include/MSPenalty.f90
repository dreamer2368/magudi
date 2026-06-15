#include "config.h"

module MSPenalty_mod

  implicit none

  type, abstract, public :: t_MSPenalty

   contains

     procedure(norm), pass, deferred :: norm
     procedure(gradientFactor), pass, deferred :: gradientFactor

  end type t_MSPenalty

  type, extends(t_MSPenalty), public :: t_QuadraticPenalty

   contains

     procedure, pass :: norm => normQuadratic
     procedure, pass :: gradientFactor => gradientFactorQuadratic

  end type t_QuadraticPenalty

  type, extends(t_MSPenalty), public :: t_HuberPenalty

     real(SCALAR_KIND) :: threshold = real(1.0e-12, SCALAR_KIND)

   contains

     procedure, pass :: norm => normHuber
     procedure, pass :: gradientFactor => gradientFactorHuber

  end type t_HuberPenalty

  abstract interface

     function norm(this, L2sqSum) result(val)

       import :: t_MSPenalty

       class(t_MSPenalty), intent(in) :: this
       SCALAR_TYPE, intent(in) :: L2sqSum

       SCALAR_TYPE :: val

     end function norm

  end interface

  abstract interface

     function gradientFactor(this, L2sqSum) result(val)

       import :: t_MSPenalty

       class(t_MSPenalty), intent(in) :: this
       SCALAR_TYPE, intent(in) :: L2sqSum

       SCALAR_TYPE :: val

     end function gradientFactor

  end interface

  interface

     function normQuadratic(this, L2sqSum) result(val)

       import :: t_QuadraticPenalty

       class(t_QuadraticPenalty), intent(in) :: this
       SCALAR_TYPE, intent(in) :: L2sqSum

       SCALAR_TYPE :: val

     end function normQuadratic

  end interface

  interface

     function gradientFactorQuadratic(this, L2sqSum) result(val)

       import :: t_QuadraticPenalty

       class(t_QuadraticPenalty), intent(in) :: this
       SCALAR_TYPE, intent(in) :: L2sqSum

       SCALAR_TYPE :: val

     end function gradientFactorQuadratic

  end interface

  interface

     function normHuber(this, L2sqSum) result(val)

       import :: t_HuberPenalty

       class(t_HuberPenalty), intent(in) :: this
       SCALAR_TYPE, intent(in) :: L2sqSum

       SCALAR_TYPE :: val

     end function normHuber

  end interface

  interface

     function gradientFactorHuber(this, L2sqSum) result(val)

       import :: t_HuberPenalty

       class(t_HuberPenalty), intent(in) :: this
       SCALAR_TYPE, intent(in) :: L2sqSum

       SCALAR_TYPE :: val

     end function gradientFactorHuber

  end interface

  public :: connectMSPenalty

  interface

     subroutine connectMSPenalty(penaltyPtr, penaltyType, huberThreshold)

       import :: t_MSPenalty

       class(t_MSPenalty), pointer, intent(out) :: penaltyPtr
       character(len = *), intent(in) :: penaltyType
       real(SCALAR_KIND), intent(in), optional :: huberThreshold

     end subroutine connectMSPenalty

  end interface

end module MSPenalty_mod
