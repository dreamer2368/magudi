#include "config.h"

module FuelSource_mod

  implicit none
  private

  type, public :: t_FuelSource

     integer :: fuelIndex
     real(SCALAR_KIND) :: location(3), amplitude, radius

     SCALAR_TYPE, allocatable :: strength(:)

   contains

     procedure, public, pass :: setup => setupFuelSource
     procedure, public, pass :: add => addFuelSource

  end type t_FuelSource

  interface

     subroutine setupFuelSource(this, fuelIndex, location, amplitude, radius)

       use Grid_mod, only : t_Grid

       import :: t_FuelSource

       class(t_FuelSource) :: this
       integer, intent(in) :: fuelIndex
       real(SCALAR_KIND), intent(in) :: location(:), amplitude, radius

     end subroutine setupFuelSource

  end interface

  interface

     subroutine addFuelSource(this, coordinates, iblank, rightHandSide)

       import :: t_FuelSource

       class(t_FuelSource) :: this
       SCALAR_TYPE, intent(in) :: coordinates(:,:)
       integer, intent(in) :: iblank(:)
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine addFuelSource

  end interface

end module FuelSource_mod
