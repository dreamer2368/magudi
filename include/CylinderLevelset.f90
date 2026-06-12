#include "config.h"

module CylinderLevelset_mod

  use LevelsetFactory_mod, only : t_LevelsetFactory

  implicit none

  type, extends(t_LevelsetFactory), public :: t_CylinderLevelset

    real(SCALAR_KIND), allocatable :: loc(:),                & ! center of the cylinder
                                      direction(:)     ! oscillation direction.
    real(SCALAR_KIND) ::              radius,                & ! radius of the cylinder.
                                      amplitude,        & ! amplitude of oscillation.
                                      period        ! period of oscillation.

   contains

     procedure, pass :: setup => setupCylinderLevelset
     procedure, pass :: cleanup => cleanupCylinderLevelset
     procedure, pass :: updateLevelset => updateCylinderLevelset

  end type t_CylinderLevelset

  interface

     subroutine setupCylinderLevelset(this, grids, states)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State

       import :: t_CylinderLevelset

       class(t_CylinderLevelset) :: this
       class(t_Grid), intent(in) :: grids(:)
       class(t_State) :: states(:)

     end subroutine setupCylinderLevelset

  end interface

  interface

     subroutine cleanupCylinderLevelset(this)

       import :: t_CylinderLevelset

       class(t_CylinderLevelset) :: this

     end subroutine cleanupCylinderLevelset

  end interface

  interface

     subroutine updateCylinderLevelset(this, mode, grids, states)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State

       import :: t_CylinderLevelset

       class(t_CylinderLevelset) :: this
       integer, intent(in) :: mode
       class(t_Grid), intent(in) :: grids(:)
       class(t_State) :: states(:)

     end subroutine updateCylinderLevelset

  end interface

end module CylinderLevelset_mod
