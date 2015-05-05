#include "config.h"

module JetExcitationPatch_mod

  use SpongePatch_mod, only : t_SpongePatch

  implicit none
  private

  type, extends(t_SpongePatch), public :: t_JetExcitationPatch

     integer :: nModes
     real(SCALAR_KIND) :: amount
     real(SCALAR_KIND), allocatable :: perturbationReal(:,:,:), perturbationImag(:,:,:), &
          angularFrequencies(:)

   contains

     procedure, pass :: setup => setupJetExcitationPatch
     procedure, pass :: cleanup => cleanupJetExcitationPatch
     procedure, pass :: verifyUsage => verifyJetExcitationPatchUsage
     procedure, pass :: updateRhs => addJetExcitation

  end type t_JetExcitationPatch

  interface

     subroutine setupJetExcitationPatch(this, index, comm, patchDescriptor,                  &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_JetExcitationPatch

       class(t_JetExcitationPatch) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupJetExcitationPatch

  end interface

  interface

     subroutine cleanupJetExcitationPatch(this)

       import :: t_JetExcitationPatch

       class(t_JetExcitationPatch) :: this

     end subroutine cleanupJetExcitationPatch

  end interface

  interface

     function verifyJetExcitationPatchUsage(this, patchDescriptor, gridSize,                 &
          normalDirection, extent, simulationFlags,                                          &
          success, message) result(isPatchUsed)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_JetExcitationPatch

       class(t_JetExcitationPatch) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       logical, intent(out) :: success
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchUsed

     end function verifyJetExcitationPatchUsage

  end interface

  interface

     subroutine addJetExcitation(this, mode, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_JetExcitationPatch

       class(t_JetExcitationPatch) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine addJetExcitation

  end interface

end module JetExcitationPatch_mod
