#include "config.h"

subroutine setupSpongePatch(this, index, comm, patchDescriptor,                              &
     grid, simulationFlags, solverOptions)

  use Grid_mod, only : t_Grid
  use SpongePatch_mod, only : t_SpongePatch
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  class(t_SpongePatch) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

end subroutine setupSpongePatch

subroutine cleanupSpongePatch(this)

  use SpongePatch_mod, only : t_SpongePatch

  class(t_SpongePatch) :: this

end subroutine cleanupSpongePatch

subroutine addDamping(this, mode, grid, state, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use Patch_mod, only : t_Patch
  use State_mod, only : t_State
  use SpongePatch_mod, only : t_SpongePatch
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_Patch) :: this
  integer, intent(in) :: mode
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

end subroutine addDamping

function verifySpongePatchUsage(this, success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use SpongePatch_mod, only : t_SpongePatch

  implicit none
  
  ! <<< Arguments >>>
  class(t_SpongePatch) :: this
  logical, intent(out) :: success
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Result >>>
  logical :: isPatchUsed

end function verifySpongePatchUsage

subroutine updateSpongePatch(this, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SpongePatch_mod, only : t_SpongePatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_SpongePatch) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state

end subroutine updateSpongePatch
