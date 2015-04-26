#include "config.h"

module PressureDrag_mod

  use Functional_mod, only : t_Functional

  implicit none
  private

  type, extends(t_Functional), public :: t_PressureDrag

     real(SCALAR_KIND) :: direction(3)

   contains

     procedure, pass :: setup => setupPressureDrag
     procedure, pass :: cleanup => cleanupPressureDrag
     procedure, pass :: compute => computePressureDrag
     procedure, pass :: computeAdjointForcing => computePressureDragAdjointForcing
     procedure, pass :: isPatchValid => isPressureDragPatchValid

  end type t_PressureDrag

  interface

     subroutine setupPressureDrag(this, region)

       use Region_mod, only : t_Region

       import :: t_PressureDrag

       class(t_PressureDrag) :: this
       class(t_Region) :: region

     end subroutine setupPressureDrag

  end interface

  interface

     subroutine cleanupPressureDrag(this)

       import :: t_PressureDrag

       class(t_PressureDrag) :: this

     end subroutine cleanupPressureDrag

  end interface

  interface

     function computePressureDrag(this, region) result(instantaneousFunctional)

       use Region_mod, only : t_Region

       import :: t_PressureDrag

       class(t_PressureDrag) :: this
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousFunctional

     end function computePressureDrag

  end interface

  interface

     subroutine computePressureDragAdjointForcing(this, simulationFlags,                     &
          solverOptions, grid, state, patch)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use CostTargetPatch_mod, only : t_CostTargetPatch
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_PressureDrag

       class(t_PressureDrag) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       class(t_CostTargetPatch) :: patch

     end subroutine computePressureDragAdjointForcing

  end interface

  interface

     function isPressureDragPatchValid(this, patchDescriptor, gridSize, normalDirection,     &
          extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_PressureDrag

       class(t_PressureDrag) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isPressureDragPatchValid

  end interface

end module PressureDrag_mod
