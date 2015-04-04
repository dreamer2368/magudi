#include "config.h"

module DragCoefficient_mod

  use Functional_mod, only : t_Functional

  implicit none
  private

  type, extends(t_Functional), public :: t_DragCoefficient

     real(SCALAR_KIND) :: freeStreamPressure, direction(3)

   contains

     procedure, pass :: setup => setupDragCoefficient
     procedure, pass :: cleanup => cleanupDragCoefficient
     procedure, pass :: compute => computeDragCoefficient
     procedure, pass :: computeAdjointForcing => computeDragCoefficientAdjointForcing

  end type t_DragCoefficient

  interface

     subroutine setupDragCoefficient(this, region)

       use Region_mod, only : t_Region

       import :: t_DragCoefficient

       class(t_DragCoefficient) :: this
       class(t_Region) :: region

     end subroutine setupDragCoefficient

  end interface

  interface
     
     subroutine cleanupDragCoefficient(this)
       
       import :: t_DragCoefficient
       
       class(t_DragCoefficient) :: this
       
     end subroutine cleanupDragCoefficient
     
  end interface

  interface
     
     function computeDragCoefficient(this, time, region) result(instantaneousFunctional)

       use Region_mod, only : t_Region
       
       import :: t_DragCoefficient
       
       class(t_DragCoefficient) :: this
       real(SCALAR_KIND), intent(in) :: time
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousFunctional
       
     end function computeDragCoefficient
     
  end interface

  interface

     subroutine computeDragCoefficientAdjointForcing(this, region)

       use Region_mod, only : t_Region

       import :: t_DragCoefficient

       class(t_DragCoefficient) :: this
       class(t_Region) :: region

     end subroutine computeDragCoefficientAdjointForcing

  end interface

end module DragCoefficient_mod
