#include "config.h"

module DragForce_mod

  use Functional_mod, only : t_Functional

  implicit none

  type, extends(t_Functional), public :: t_DragForce

     real(SCALAR_KIND) :: direction(3)

   contains

     procedure, pass :: setup => setupDragForce
     procedure, pass :: cleanup => cleanupDragForce
     procedure, pass :: compute => computeDragForce
     procedure, pass :: computeAdjointForcing => computeDragForceAdjointForcing
     procedure, pass :: isPatchValid => isDragForcePatchValid

  end type t_DragForce

  interface

     subroutine setupDragForce(this, region)

       use Region_mod, only : t_Region

       import :: t_DragForce

       class(t_DragForce) :: this
       class(t_Region) :: region

     end subroutine setupDragForce

  end interface

  interface

     subroutine cleanupDragForce(this)

       import :: t_DragForce

       class(t_DragForce) :: this

     end subroutine cleanupDragForce

  end interface

  interface

     function computeDragForce(this, region) result(instantaneousFunctional)

       use Region_mod, only : t_Region

       import :: t_DragForce

       class(t_DragForce) :: this
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousFunctional

     end function computeDragForce

  end interface

  interface

     subroutine computeDragForceAdjointForcing(this, simulationFlags,                        &
          solverOptions, grid, state, patch)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use CostTargetPatch_mod, only : t_CostTargetPatch
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_DragForce

       class(t_DragForce) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       class(t_CostTargetPatch) :: patch

     end subroutine computeDragForceAdjointForcing

  end interface

  interface

     function isDragForcePatchValid(this, patchDescriptor, gridSize, normalDirection,        &
          extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_DragForce

       class(t_DragForce) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isDragForcePatchValid

  end interface

end module DragForce_mod
