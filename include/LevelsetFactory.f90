#include "config.h"

module LevelsetFactory_mod

  implicit none

  type, abstract, public :: t_LevelsetFactory

   contains

     procedure(setup), pass, deferred :: setup
     procedure(cleanup), pass, deferred :: cleanup
     procedure(updateLevelset), pass, deferred :: updateLevelset

  end type t_LevelsetFactory

  abstract interface

     subroutine setup(this, grids, states)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State

       import :: t_LevelsetFactory

       class(t_LevelsetFactory) :: this
       class(t_Grid), intent(in) :: grids(:)
       class(t_State) :: states(:)

     end subroutine setup

  end interface

  abstract interface

     subroutine cleanup(this)

       import :: t_LevelsetFactory

       class(t_LevelsetFactory) :: this

     end subroutine cleanup

  end interface

  abstract interface

     subroutine updateLevelset(this, mode, grids, states)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State

       import :: t_LevelsetFactory

       class(t_LevelsetFactory) :: this
       integer, intent(in) :: mode
       class(t_Grid), intent(in) :: grids(:)
       class(t_State) :: states(:)

     end subroutine updateLevelset

  end interface

end module LevelsetFactory_mod
