#include "config.h"

module WallShapePatch_mod

  use Patch_mod, only : t_Patch

  implicit none
  private

  type, extends(t_Patch), public :: t_WallShapePatch

     integer :: nModes
     integer :: nParams
     real(SCALAR_KIND)::amplitude,                                    
     SCALAR_TYPE, allocatable :: phases(:),amplitudes(:)

   contains

     procedure, pass :: setup => setupWallShapePatch
     procedure, pass :: cleanup => cleanupWallShapePatch
     procedure, pass :: verifyUsage => verifyWallShapePatchUsage
     procedure, pass :: updateRhs => updateWallShape

  end type t_WallShapePatch

  interface

     subroutine setupWallShapePatch(this, index, comm, patchDescriptor,           &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_WallShapePatch

       class(t_WallShapePatch) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupWallShapePatch

  end interface

  interface

     subroutine cleanupWallShapePatch(this)

       import :: t_WallShapePatch

       class(t_WallShapePatch) :: this

     end subroutine cleanupWallShapePatch

  end interface

  interface

     function verifyWallShapePatchUsage(this, patchDescriptor, gridSize,          &
          normalDirection, extent, simulationFlags,                                          &
          success, message) result(isPatchUsed)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_WallShapePatch

       class(t_WallShapePatch) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       logical, intent(out) :: success
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchUsed

     end function verifyWallShapePatchUsage

  end interface

  interface

     subroutine addWallShape(this, mode, simulationFlags,                         &
          solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_WallShapePatch

       class(t_WallShapePatch) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine addWallShape

  end interface

end module WallShapePatch_mod
