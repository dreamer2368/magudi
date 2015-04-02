#include "config.h"

module ImpenetrableWall_mod

  use Patch_mod, only : t_Patch

  implicit none
  private

  type, extends(t_Patch), public :: t_ImpenetrableWall

     real(SCALAR_KIND) :: inviscidPenaltyAmount

   contains

     procedure, pass :: setup => setupImpenetrableWall
     procedure, pass :: cleanup => cleanupImpenetrableWall
     procedure, pass :: update => updateImpenetrableWall
     procedure, pass :: verifyUsage => verifyImpenetrableWallUsage
     procedure, pass :: updateRhs => addImpenetrableWallPenalty

  end type t_ImpenetrableWall

  interface

     subroutine setupImpenetrableWall(this, index, comm, patchDescriptor,                    &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ImpenetrableWall

       class(t_ImpenetrableWall) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupImpenetrableWall

  end interface

  interface

     subroutine cleanupImpenetrableWall(this)

       import :: t_ImpenetrableWall

       class(t_ImpenetrableWall) :: this

     end subroutine cleanupImpenetrableWall

  end interface

  interface

     subroutine updateImpenetrableWall(this, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ImpenetrableWall

       class(t_ImpenetrableWall) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state

     end subroutine updateImpenetrableWall

  end interface

  interface

     function verifyImpenetrableWallUsage(this, patchDescriptor, gridSize, normalDirection,  &
          extent, simulationFlags, success, message) result(isPatchUsed)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ImpenetrableWall

       class(t_ImpenetrableWall) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       logical, intent(out) :: success
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchUsed

     end function verifyImpenetrableWallUsage

  end interface

  interface

     subroutine addImpenetrableWallPenalty(this, mode, simulationFlags,                      &
          solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ImpenetrableWall

       class(t_ImpenetrableWall) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine addImpenetrableWallPenalty

  end interface

end module ImpenetrableWall_mod
