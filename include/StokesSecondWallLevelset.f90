#include "config.h"

module StokesSecondWallLevelset_mod

  use LevelsetFactory_mod, only : t_LevelsetFactory

  implicit none

  type, extends(t_LevelsetFactory), public :: t_StokesSecondWallLevelset

    real(SCALAR_KIND), allocatable :: levelsetLoc(:),        & ! this point lies on level set plane.
                                      levelsetNormal(:),     & ! normal direction of the plane.
                                      levelsetDirection(:),  & ! oscillation direction.
                                      levelsetAmp,           & ! amplitude of oscillation.
                                      levelsetPeriod           ! period of oscillation.

   contains

     procedure, pass :: setup => setupStokesSecondWallLevelset
     procedure, pass :: cleanup => cleanupStokesSecondWallLevelset
     procedure, pass :: updateLevelset => updateStokesSecondWallLevelset

  end type t_StokesSecondWallLevelset

  interface

     subroutine setupStokesSecondWallLevelset(this, grids, states)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State

       import :: t_StokesSecondWallLevelset

       class(t_StokesSecondWallLevelset) :: this
       class(t_Grid), intent(in) :: grids(:)
       class(t_State) :: states(:)

     end subroutine setupStokesSecondWallLevelset

  end interface

  interface

     subroutine cleanupStokesSecondWallLevelset(this)

       import :: t_StokesSecondWallLevelset

       class(t_StokesSecondWallLevelset) :: this

     end subroutine cleanupStokesSecondWallLevelset

  end interface

  interface

     subroutine updateStokesSecondWallLevelset(this, mode, grids, states)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State

       import :: t_StokesSecondWallLevelset

       class(t_StokesSecondWallLevelset) :: this
       integer, intent(in) :: mode
       class(t_Grid), intent(in) :: grids(:)
       class(t_State) :: states(:)

     end subroutine updateStokesSecondWallLevelset

  end interface

end module StokesSecondWallLevelset_mod
