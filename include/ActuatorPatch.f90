#include "config.h"

module ActuatorPatch_mod

  use Patch_mod, only : t_Patch

  implicit none
  private

  type, extends(t_Patch), public :: t_ActuatorPatch

     SCALAR_TYPE, allocatable :: controlForcing(:,:)

   contains

     procedure, pass :: setup => setupActuatorPatch
     procedure, pass :: cleanup => cleanupActuatorPatch
     procedure, pass :: update => updateActuatorPatch
     procedure, pass :: verifyUsage => verifyActuatorPatchUsage
     procedure, pass :: updateRhs => addControlForcing

  end type t_ActuatorPatch

  interface

     subroutine setupActuatorPatch(this, index, comm, patchDescriptor,                       &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ActuatorPatch

       class(t_ActuatorPatch) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupActuatorPatch

  end interface

  interface

     subroutine cleanupActuatorPatch(this)

       import :: t_ActuatorPatch

       class(t_ActuatorPatch) :: this

     end subroutine cleanupActuatorPatch

  end interface

  interface

     subroutine updateActuatorPatch(this, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ActuatorPatch

       class(t_ActuatorPatch) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state

     end subroutine updateActuatorPatch

  end interface

  interface

     function verifyActuatorPatchUsage(this, patchDescriptor, gridSize, normalDirection,     &
          extent, simulationFlags, success, message) result(isPatchUsed)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ActuatorPatch

       class(t_ActuatorPatch) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       logical, intent(out) :: success
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchUsed

     end function verifyActuatorPatchUsage

  end interface

  interface

     subroutine addControlForcing(this, mode, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ActuatorPatch

       class(t_ActuatorPatch) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine addControlForcing

  end interface

end module ActuatorPatch_mod
