#include "config.h"

module GaussianIgnitionPatch_mod

  use Patch_mod, only : t_Patch

  implicit none
  private

  integer, parameter, private :: real64 = selected_real_kind(15)

  type, extends(t_Patch), public :: t_GaussianIgnitionPatch

     real(SCALAR_KIND) :: origin(2), amplitude, radius,                                      &
          timeStart, timeDuration

     SCALAR_TYPE, allocatable :: strength(:)

   contains

     procedure, pass :: setup => setupGaussianIgnitionPatch
     procedure, pass :: cleanup => cleanupGaussianIgnitionPatch
     procedure, pass :: verifyUsage => verifyGaussianIgnitionPatchUsage
     procedure, pass :: updateRhs => addGaussianIgnition

  end type t_GaussianIgnitionPatch

  interface

     subroutine setupGaussianIgnitionPatch(this, index, comm, patchDescriptor,               &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_GaussianIgnitionPatch

       class(t_GaussianIgnitionPatch) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupGaussianIgnitionPatch

  end interface

  interface

     subroutine cleanupGaussianIgnitionPatch(this)

       import :: t_GaussianIgnitionPatch

       class(t_GaussianIgnitionPatch) :: this

     end subroutine cleanupGaussianIgnitionPatch

  end interface

  interface

     function verifyGaussianIgnitionPatchUsage(this, patchDescriptor, gridSize,          &
          normalDirection, extent, simulationFlags,                                          &
          success, message) result(isPatchUsed)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_GaussianIgnitionPatch

       class(t_GaussianIgnitionPatch) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       logical, intent(out) :: success
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchUsed

     end function verifyGaussianIgnitionPatchUsage

  end interface

  interface

     subroutine addGaussianIgnition(this, mode, simulationFlags,                         &
          solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_GaussianIgnitionPatch

       class(t_GaussianIgnitionPatch) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine addGaussianIgnition

  end interface

end module GaussianIgnitionPatch_mod
