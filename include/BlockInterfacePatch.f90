#include "config.h"

module BlockInterfacePatch_enum

  implicit none
  public

  integer, parameter ::  METRICS = 3081

end module BlockInterfacePatch_enum

module BlockInterfacePatch_mod

  use Patch_mod, only : t_Patch

  implicit none

  type, extends(t_Patch), public :: t_BlockInterfacePatch

     real(SCALAR_KIND) :: inviscidPenaltyAmount, viscousPenaltyAmount,                       &
                          inviscidPenaltyAmountL, inviscidPenaltyAmountR,                    &
                          viscousPenaltyAmountL, viscousPenaltyAmountR
     SCALAR_TYPE, allocatable :: conservedVariablesL(:,:), conservedVariablesR(:,:),         &
          adjointVariablesL(:,:), adjointVariablesR(:,:),                                    &
          cartesianViscousFluxesL(:,:,:), viscousFluxesL(:,:), viscousFluxesR(:,:)
     SCALAR_TYPE, allocatable :: sendBuffer(:,:), receiveBuffer(:,:)
     SCALAR_TYPE, allocatable :: metricsAlongNormalDirectionL(:,:),                          &
                                 metricsAlongNormalDirectionR(:,:)

   contains

     procedure, pass :: setup => setupBlockInterfacePatch
     procedure, pass :: cleanup => cleanupBlockInterfacePatch
     procedure, pass :: verifyUsage => verifyBlockInterfacePatchUsage
     procedure, pass :: updateRhs => addBlockInterfacePenalty
     procedure, pass :: collectInterfaceData
     procedure, pass :: disperseInterfaceData
     procedure, pass :: reshapeReceivedData

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

  interface

     subroutine collectInterfaceData(this, mode, simulationFlags, solverOptions, grid, state)

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
       class(t_State), intent(in) :: state

     end subroutine collectInterfaceData

  end interface

  interface

     subroutine disperseInterfaceData(this, mode, simulationFlags, solverOptions)

       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_BlockInterfacePatch

       class(t_BlockInterfacePatch) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine disperseInterfaceData

  end interface

  interface

     subroutine reshapeReceivedData(this, indexReordering)

       import :: t_BlockInterfacePatch

       class(t_BlockInterfacePatch) :: this
       integer, intent(in) :: indexReordering(3)

     end subroutine reshapeReceivedData

  end interface

end module BlockInterfacePatch_mod
