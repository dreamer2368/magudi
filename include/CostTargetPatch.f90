#include "config.h"

module CostTargetPatch_mod

  use Patch_mod, only : t_Patch

  implicit none
  private

  type, extends(t_Patch), public :: t_CostTargetPatch

     real(SCALAR_KIND) :: inviscidPenaltyAmount
     SCALAR_TYPE, allocatable :: norm(:,:), adjointForcing(:,:)

   contains

     procedure, pass :: setup => setupCostTargetPatch
     procedure, pass :: cleanup => cleanupCostTargetPatch
     procedure, pass :: verifyUsage => verifyCostTargetPatchUsage
     procedure, pass :: updateRhs => addAdjointForcing
     generic :: computeInnerProduct => computeScalarInnerProductOnPatch,                     &
          computeVectorInnerProductOnPatch
     procedure, private, pass :: computeScalarInnerProductOnPatch
     procedure, private, pass :: computeVectorInnerProductOnPatch

  end type t_CostTargetPatch

  interface

     subroutine setupCostTargetPatch(this, index, comm, patchDescriptor,                     &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_CostTargetPatch

       class(t_CostTargetPatch) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupCostTargetPatch

  end interface

  interface

     subroutine cleanupCostTargetPatch(this)

       import :: t_CostTargetPatch

       class(t_CostTargetPatch) :: this

     end subroutine cleanupCostTargetPatch

  end interface

  interface

     function verifyCostTargetPatchUsage(this, patchDescriptor, gridSize, normalDirection,   &
          extent, simulationFlags, success, message) result(isPatchUsed)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_CostTargetPatch

       class(t_CostTargetPatch) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       logical, intent(out) :: success
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchUsed

     end function verifyCostTargetPatchUsage

  end interface

  interface

     subroutine addAdjointForcing(this, mode, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_CostTargetPatch

       class(t_CostTargetPatch) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine addAdjointForcing

  end interface

  interface computeInnerProduct

     function computeScalarInnerProductOnPatch(this, grid, f, g, weight) result(innerProduct)

       use Grid_mod, only : t_Grid

       import :: t_CostTargetPatch

       class(t_CostTargetPatch) :: this
       class(t_Grid) :: grid
       SCALAR_TYPE, intent(in) :: f(:), g(:)

       SCALAR_TYPE, intent(in), optional :: weight(:)

       SCALAR_TYPE :: innerProduct

     end function computeScalarInnerProductOnPatch

     function computeVectorInnerProductOnPatch(this, grid, f, g, weight) result(innerProduct)

       use Grid_mod, only : t_Grid

       import :: t_CostTargetPatch

       class(t_CostTargetPatch) :: this
       class(t_Grid) :: grid
       SCALAR_TYPE, intent(in) :: f(:,:), g(:,:)

       SCALAR_TYPE, intent(in), optional :: weight(:)

       SCALAR_TYPE :: innerProduct

     end function computeVectorInnerProductOnPatch

  end interface computeInnerProduct

end module CostTargetPatch_mod
