#include "config.h"

module ReactantDepletion_mod

  use Functional_mod, only : t_Functional

  implicit none

  type, extends(t_Functional), public :: t_ReactantDepletion

     integer :: reactant

   contains

     procedure, pass :: setup => setupReactantDepletion
     procedure, pass :: cleanup => cleanupReactantDepletion
     procedure, pass :: compute => computeReactantDepletion
     procedure, pass :: computeAdjointForcing => computeReactantDepletionAdjointForcing
     procedure, pass :: isPatchValid => isReactantDepletionPatchValid

  end type t_ReactantDepletion

  interface

     subroutine setupReactantDepletion(this, region)

       use Region_mod, only : t_Region

       import :: t_ReactantDepletion

       class(t_ReactantDepletion) :: this
       class(t_Region) :: region

     end subroutine setupReactantDepletion

  end interface

  interface

     subroutine cleanupReactantDepletion(this)

       import :: t_ReactantDepletion

       class(t_ReactantDepletion) :: this

     end subroutine cleanupReactantDepletion

  end interface

  interface

     function computeReactantDepletion(this, region) result(instantaneousFunctional)

       use Region_mod, only : t_Region

       import :: t_ReactantDepletion

       class(t_ReactantDepletion) :: this
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousFunctional

     end function computeReactantDepletion

  end interface

  interface

     subroutine computeReactantDepletionAdjointForcing(this, simulationFlags,                        &
          solverOptions, grid, state, patch)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use CostTargetPatch_mod, only : t_CostTargetPatch
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ReactantDepletion

       class(t_ReactantDepletion) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       class(t_CostTargetPatch) :: patch

     end subroutine computeReactantDepletionAdjointForcing

  end interface

  interface

     function isReactantDepletionPatchValid(this, patchDescriptor, gridSize,                 &
          normalDirection, extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ReactantDepletion

       class(t_ReactantDepletion) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isReactantDepletionPatchValid

  end interface

end module ReactantDepletion_mod
