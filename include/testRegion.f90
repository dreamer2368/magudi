#include "config.h"

module testRegion_mod

  use Region_mod

  implicit none

  type, extends(t_Region), public :: t_testRegion

     SCALAR_TYPE, allocatable :: tempRhs(:)
     integer :: tempRhsIndex

   contains

     procedure, pass :: computeRhs=>testComputeRhs

  end type t_testRegion

  interface
    subroutine testComputeRhs(this,mode, timeStep, stage)
      import :: t_testRegion

      class(t_testRegion) :: this
      integer, intent(in) :: mode
      integer, intent(in), optional :: timeStep, stage
    end subroutine testComputeRhs
  end interface

end module testRegion_mod
