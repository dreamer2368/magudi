#include "config.h"

module SolenoidalExcitationPatch_mod

  use Patch_mod, only : t_Patch

  implicit none
  private

  integer, parameter, private :: real64 = selected_real_kind(15)

  type, extends(t_Patch), public :: t_SolenoidalExcitationPatch

     integer :: nModes
     real(SCALAR_KIND) :: location(3), speed(3), amplitude,                                  &
          gaussianFactor, frequency

     real(real64), allocatable :: angularFrequencies(:), phases(:,:)

   contains

     procedure, pass :: setup => setupSolenoidalExcitationPatch
     procedure, pass :: cleanup => cleanupSolenoidalExcitationPatch
     procedure, pass :: update => updateSolenoidalExcitationPatch
     procedure, pass :: verifyUsage => verifySolenoidalExcitationPatchUsage
     procedure, pass :: updateRhs => addSolenoidalExcitation

  end type t_SolenoidalExcitationPatch

  interface

     subroutine setupSolenoidalExcitationPatch(this, index, comm, patchDescriptor,           &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_SolenoidalExcitationPatch

       class(t_SolenoidalExcitationPatch) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupSolenoidalExcitationPatch

  end interface

  interface

     subroutine cleanupSolenoidalExcitationPatch(this)

       import :: t_SolenoidalExcitationPatch

       class(t_SolenoidalExcitationPatch) :: this

     end subroutine cleanupSolenoidalExcitationPatch

  end interface

  interface

     subroutine updateSolenoidalExcitationPatch(this, simulationFlags,                       &
          solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_SolenoidalExcitationPatch

       class(t_SolenoidalExcitationPatch) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state

     end subroutine updateSolenoidalExcitationPatch

  end interface

  interface

     function verifySolenoidalExcitationPatchUsage(this, patchDescriptor, gridSize,          &
          normalDirection, extent, simulationFlags,                                          &
          success, message) result(isPatchUsed)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_SolenoidalExcitationPatch

       class(t_SolenoidalExcitationPatch) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       logical, intent(out) :: success
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchUsed

     end function verifySolenoidalExcitationPatchUsage

  end interface

  interface

     subroutine addSolenoidalExcitation(this, mode, simulationFlags,                         &
          solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_SolenoidalExcitationPatch

       class(t_SolenoidalExcitationPatch) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine addSolenoidalExcitation

  end interface

end module SolenoidalExcitationPatch_mod
