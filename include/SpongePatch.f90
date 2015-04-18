#include "config.h"

module SpongePatch_mod

  use Patch_mod, only : t_Patch

  implicit none
  private

  type, extends(t_Patch), public :: t_SpongePatch

     real(SCALAR_KIND) :: spongeAmount
     integer :: spongeExponent
     real(SCALAR_KIND), allocatable :: spongeStrength(:)

   contains

     procedure, pass :: setup => setupSpongePatch
     procedure, pass :: cleanup => cleanupSpongePatch
     procedure, pass :: verifyUsage => verifySpongePatchUsage
     procedure, pass :: updateRhs => addDamping

  end type t_SpongePatch

  interface

     subroutine setupSpongePatch(this, index, comm, patchDescriptor,                         &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_SpongePatch

       class(t_SpongePatch) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupSpongePatch

  end interface

  interface

     subroutine cleanupSpongePatch(this)

       import :: t_SpongePatch

       class(t_SpongePatch) :: this

     end subroutine cleanupSpongePatch

  end interface

  interface

     function verifySpongePatchUsage(this, patchDescriptor, gridSize, normalDirection,       &
          extent, simulationFlags, success, message) result(isPatchUsed)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_SpongePatch

       class(t_SpongePatch) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       logical, intent(out) :: success
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchUsed

     end function verifySpongePatchUsage

  end interface

  interface

     subroutine addDamping(this, mode, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_SpongePatch

       class(t_SpongePatch) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine addDamping

  end interface

end module SpongePatch_mod
