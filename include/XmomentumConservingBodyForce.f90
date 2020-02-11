#include "config.h"

module XmomentumConservingBodyForce_mod

  implicit none

  type, public :: t_BodyForce

     real(SCALAR_KIND) :: initialXmomentum, oneOverVolume, momentumLossPerVolume

   contains

     procedure, public, pass :: setup => setupBodyForce
     procedure, public, pass :: add => addBodyForce

  end type t_BodyForce

  interface

     subroutine setupBodyForce(this, region)

       use Region_mod, only : t_Region

       import :: t_BodyForce

       class(t_BodyForce) :: this
       class(t_Region) :: region

     end subroutine setupBodyForce

  end interface

  interface

     subroutine addBodyForce(this, region)

       use Region_mod, only : t_Region

       import :: t_BodyForce

       class(t_BodyForce) :: this
       class(t_Region) :: region

     end subroutine addBodyForce

  end interface

end module XmomentumConservingBodyForce_mod
