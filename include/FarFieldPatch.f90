#include "config.h"

module FarFieldPatch_mod

  use Patch_mod, only : t_Patch

  implicit none
  private

  type, extends(t_Patch), public :: t_FarFieldPatch

     real(SCALAR_KIND) :: inviscidPenaltyAmount, viscousPenaltyAmount
     SCALAR_TYPE, allocatable :: viscousFluxes(:,:), targetViscousFluxes(:,:)

   contains

     procedure, pass :: setup => setupFarFieldPatch
     procedure, pass :: cleanup => cleanupFarFieldPatch
     procedure, pass :: update => updateFarFieldPatch
     procedure, pass :: verifyUsage => verifyFarFieldPatchUsage
     procedure, pass :: updateRhs => addFarFieldPenalty

  end type t_FarFieldPatch

  interface

     subroutine setupFarFieldPatch(this, index, comm, patchDescriptor,                       &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_FarFieldPatch

       class(t_FarFieldPatch) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupFarFieldPatch

  end interface

  interface

     subroutine cleanupFarFieldPatch(this)

       import :: t_FarFieldPatch

       class(t_FarFieldPatch) :: this

     end subroutine cleanupFarFieldPatch

  end interface

  interface
     
     subroutine updateFarFieldPatch(this, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags
       
       import :: t_FarFieldPatch
       
       class(t_FarFieldPatch) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       
     end subroutine updateFarFieldPatch
     
  end interface

  interface

     function verifyFarFieldPatchUsage(this, patchDescriptor, gridSize, normalDirection,       &
          extent, simulationFlags, success, message) result(isPatchUsed)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_FarFieldPatch
       
       class(t_FarFieldPatch) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       logical, intent(out) :: success
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchUsed

     end function verifyFarFieldPatchUsage

  end interface

  interface

     subroutine addFarFieldPenalty(this, mode, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_FarFieldPatch

       class(t_FarFieldPatch) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine addFarFieldPenalty

  end interface

end module FarFieldPatch_mod
