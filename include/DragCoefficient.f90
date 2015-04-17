#include "config.h"

module DragCoefficient_mod

  use Functional_mod, only : t_Functional

  implicit none
  private

  type, extends(t_Functional), public :: t_DragCoefficient

     real(SCALAR_KIND) :: direction(3)

   contains

     procedure, pass :: setup => setupDragCoefficient
     procedure, pass :: cleanup => cleanupDragCoefficient
     procedure, pass :: compute => computeDragCoefficient
     procedure, pass :: computeAdjointForcing => computeDragCoefficientAdjointForcing
     procedure, pass :: isPatchValid => isDragCoefficientPatchValid

  end type t_DragCoefficient

  interface

     subroutine setupDragCoefficient(this, region)

       use Region_mod, only : t_Region

       import :: t_DragCoefficient

       class(t_DragCoefficient) :: this
       class(t_Region) :: region

     end subroutine setupDragCoefficient

  end interface

  interface

     subroutine cleanupDragCoefficient(this)

       import :: t_DragCoefficient

       class(t_DragCoefficient) :: this

     end subroutine cleanupDragCoefficient

  end interface

  interface

     function computeDragCoefficient(this, region) result(instantaneousFunctional)

       use Region_mod, only : t_Region

       import :: t_DragCoefficient

       class(t_DragCoefficient) :: this
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousFunctional

     end function computeDragCoefficient

  end interface

  interface

     subroutine computeDragCoefficientAdjointForcing(this, simulationFlags,                  &
          solverOptions, grid, state, patch)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use CostTargetPatch_mod, only : t_CostTargetPatch
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_DragCoefficient

       class(t_DragCoefficient) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       class(t_CostTargetPatch) :: patch

     end subroutine computeDragCoefficientAdjointForcing

  end interface

  interface

     function isDragCoefficientPatchValid(this, patchDescriptor, gridSize, normalDirection,  &
          extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_DragCoefficient

       class(t_DragCoefficient) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isDragCoefficientPatchValid

  end interface

end module DragCoefficient_mod
