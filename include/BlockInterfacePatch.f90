#include "config.h"

module BlockInterfacePatch_mod

  use Patch_mod, only : t_Patch

  implicit none
  private

  type, extends(t_Patch), public :: t_BlockInterfacePatch

     real(SCALAR_KIND) :: inviscidPenaltyAmount, viscousPenaltyAmount
     SCALAR_TYPE, allocatable :: viscousFluxes(:,:), interfaceViscousFluxes(:,:),            &
          interfaceConservedVariables(:,:)

   contains

     procedure, pass :: setup => setupBlockInterfacePatch
     procedure, pass :: cleanup => cleanupBlockInterfacePatch
     procedure, pass :: update => updateBlockInterfacePatch
     procedure, pass :: verifyUsage => verifyBlockInterfacePatchUsage
     procedure, pass :: updateRhs => addBlockInterfacePenalty

  end type t_BlockInterfacePatch

  interface

     subroutine setupBlockInterfacePatch(this, index, comm, patchDescriptor,                 &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_BlockInterfacePatch

       class(t_BlockInterfacePatch) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupBlockInterfacePatch

  end interface

  interface

     subroutine cleanupBlockInterfacePatch(this)

       import :: t_BlockInterfacePatch

       class(t_BlockInterfacePatch) :: this

     end subroutine cleanupBlockInterfacePatch

  end interface

  interface

     subroutine updateBlockInterfacePatch(this, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_BlockInterfacePatch

       class(t_BlockInterfacePatch) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state

     end subroutine updateBlockInterfacePatch

  end interface

  interface

     function verifyBlockInterfacePatchUsage(this, patchDescriptor, gridSize,                &
          normalDirection, extent, simulationFlags,                                          &
          success, message) result(isPatchUsed)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_BlockInterfacePatch

       class(t_BlockInterfacePatch) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       logical, intent(out) :: success
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchUsed

     end function verifyBlockInterfacePatchUsage

  end interface

  interface

     subroutine addBlockInterfacePenalty(this, mode, simulationFlags,                        &
          solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_BlockInterfacePatch

       class(t_BlockInterfacePatch) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine addBlockInterfacePenalty

  end interface

end module BlockInterfacePatch_mod