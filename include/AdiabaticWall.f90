#include "config.h"

module AdiabaticWall_mod

  use ImpenetrableWall_mod, only : t_ImpenetrableWall

  implicit none

  type, extends(t_ImpenetrableWall), public :: t_AdiabaticWall

     real(SCALAR_KIND) :: viscousPenaltyAmounts(2)

   contains

     procedure, pass :: setup => setupAdiabaticWall
     procedure, pass :: cleanup => cleanupAdiabaticWall
     procedure, pass :: verifyUsage => verifyAdiabaticWallUsage
     procedure, pass :: updateRhs => addAdiabaticWallPenalty

  end type t_AdiabaticWall

  interface

     subroutine setupAdiabaticWall(this, index, comm, patchDescriptor,                       &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_AdiabaticWall

       class(t_AdiabaticWall) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupAdiabaticWall

  end interface

  interface

     subroutine cleanupAdiabaticWall(this)

       import :: t_AdiabaticWall

       class(t_AdiabaticWall) :: this

     end subroutine cleanupAdiabaticWall

  end interface

  interface

     function verifyAdiabaticWallUsage(this, patchDescriptor, gridSize, normalDirection,     &
          extent, simulationFlags, success, message) result(isPatchUsed)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_AdiabaticWall

       class(t_AdiabaticWall) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       logical, intent(out) :: success
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchUsed

     end function verifyAdiabaticWallUsage

  end interface

  interface

     subroutine addAdiabaticWallPenalty(this, mode, simulationFlags,                         &
          solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_AdiabaticWall

       class(t_AdiabaticWall) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine addAdiabaticWallPenalty

  end interface

end module AdiabaticWall_mod
