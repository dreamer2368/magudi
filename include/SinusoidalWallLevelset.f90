#include "config.h"

module SinusoidalWallLevelset_mod

  use LevelsetFactory_mod, only : t_LevelsetFactory

  implicit none

  type, private :: t_SWInternal
     SCALAR_TYPE, allocatable :: buffer(:)
  end type t_SWInternal

  type, extends(t_LevelsetFactory), public :: t_SinusoidalWallLevelset

    type(t_SWInternal), dimension(:), allocatable :: wallShapes
    real(SCALAR_KIND) :: levelsetLoc,           & ! center location of deforming part.
                         levelsetWidth,         & ! width of the deforming part.
                         levelsetHeight,        & ! height of the wall shape.
                         levelsetAmp,           & ! deformation amplitude
                         levelsetPeriod,        & ! deformation period
                         levelsetLocOffset,     & ! location offset for wall shape.
                         levelsetHeightOffset     ! height offset for wall shape.

   contains

     procedure, pass :: setup => setupSinusoidalWallLevelset
     procedure, pass :: cleanup => cleanupSinusoidalWallLevelset
     procedure, pass :: updateLevelset => updateSinusoidalWallLevelset

  end type t_SinusoidalWallLevelset

  interface

     subroutine setupSinusoidalWallLevelset(this, grids, states)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State

       import :: t_SinusoidalWallLevelset

       class(t_SinusoidalWallLevelset) :: this
       class(t_Grid), intent(in) :: grids(:)
       class(t_State) :: states(:)

     end subroutine setupSinusoidalWallLevelset

  end interface

  interface

     subroutine cleanupSinusoidalWallLevelset(this)

       import :: t_SinusoidalWallLevelset

       class(t_SinusoidalWallLevelset) :: this

     end subroutine cleanupSinusoidalWallLevelset

  end interface

  interface

     subroutine updateSinusoidalWallLevelset(this, mode, grids, states)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State

       import :: t_SinusoidalWallLevelset

       class(t_SinusoidalWallLevelset) :: this
       integer, intent(in) :: mode
       class(t_Grid), intent(in) :: grids(:)
       class(t_State) :: states(:)

     end subroutine updateSinusoidalWallLevelset

  end interface

end module SinusoidalWallLevelset_mod
