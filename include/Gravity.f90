#include "config.h"

module Gravity_mod

  implicit none

  type, public :: t_Gravity

   contains

     procedure, public, pass :: add => addGravity
     procedure, public, pass :: addAdjoint => addGravityAdjoint

  end type t_Gravity

  interface

     subroutine addGravity(this, iblank, density, froudeNumberInverse, rightHandSide)

       import :: t_Gravity

       class(t_Gravity) :: this
       SCALAR_TYPE, intent(in) :: density(:), froudeNumberInverse(:)
       integer, intent(in) :: iblank(:)
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine addGravity

  end interface

  interface

     subroutine addGravityAdjoint(this, iblank, adjointVariables, froudeNumberInverse,             &
          rightHandSide)

       import :: t_Gravity

       class(t_Gravity) :: this
       SCALAR_TYPE, intent(in) :: adjointVariables(:,:), froudeNumberInverse(:)
       integer, intent(in) :: iblank(:)
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine addGravityAdjoint

  end interface

end module Gravity_mod
