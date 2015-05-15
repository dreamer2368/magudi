#include "config.h"

module Patch_factory

  use Patch_mod, only : t_Patch

  implicit none
  private

  type, public :: t_PatchFactory

     class(t_Patch), pointer :: patch => null()
     character(len = STRING_LENGTH) :: patchType = ""

   contains

     procedure, pass :: connect => connectPatch
     procedure, pass :: cleanup => cleanupPatchFactory

  end type t_PatchFactory

  interface

     subroutine connectPatch(this, patchTarget, patchType, createNew)

       use Patch_mod, only : t_Patch

       import :: t_PatchFactory

       class(t_PatchFactory) :: this
       class(t_Patch), pointer, intent(out) :: patchTarget
       character(len = *), intent(in), optional :: patchType
       logical, intent(in), optional :: createNew

     end subroutine connectPatch

  end interface

  interface

     subroutine cleanupPatchFactory(this)

       import :: t_PatchFactory

       class(t_PatchFactory) :: this

     end subroutine cleanupPatchFactory

  end interface

  interface

     function queryPatchTypeExists(patchFactories, patchType,                                &
          gridIndex, normalDirection) result(patchTypeExists)

       import :: t_PatchFactory

       type(t_PatchFactory), allocatable :: patchFactories(:)
       character(len = *), intent(in) :: patchType

       integer, intent(in), optional :: gridIndex, normalDirection

       logical :: patchTypeExists

     end function queryPatchTypeExists

  end interface

  interface

     subroutine computeSpongeStrengths(patchFactories, grid)

       use Grid_mod, only : t_Grid

       import :: t_PatchFactory

       type(t_PatchFactory), allocatable :: patchFactories(:)
       class(t_Grid), intent(in) :: grid

     end subroutine computeSpongeStrengths

  end interface

  interface

     function computeQuadratureOnPatches(patchFactories, patchType,                          &
          grid, integrand) result(integral)

       use Grid_mod, only : t_Grid

       import :: t_PatchFactory

       type(t_PatchFactory), allocatable :: patchFactories(:)
       character(len = *), intent(in) :: patchType
       class(t_Grid), intent(in) :: grid
       SCALAR_TYPE, intent(in) :: integrand(:)

       SCALAR_TYPE :: integral

     end function computeQuadratureOnPatches

  end interface

  interface

     subroutine updatePatchFactories(patchFactories, simulationFlags,                        &
          solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_PatchFactory

       type(t_PatchFactory), allocatable :: patchFactories(:)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid) :: grid
       class(t_State) :: state

     end subroutine updatePatchFactories

  end interface

  interface

     subroutine computeFarFieldAdjointViscousPenalty(patchFactories, simulationFlags,        &
          solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_PatchFactory

       type(t_PatchFactory), allocatable :: patchFactories(:)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid) :: grid
       class(t_State) :: state

     end subroutine computeFarFieldAdjointViscousPenalty

  end interface

  public :: computeSpongeStrengths, computeQuadratureOnPatches, updatePatchFactories,        &
       computeFarFieldAdjointViscousPenalty, queryPatchTypeExists

end module Patch_factory
