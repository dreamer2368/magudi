#include "config.h"

module LighthillTensorComponent_mod

  use Functional_mod, only : t_Functional

  implicit none

  type, extends(t_Functional), public :: t_LighthillTensorComponent

     SCALAR_TYPE :: timeWindowCenter, timeWindowWidth
     logical :: useTimeWindow = .false., viscosityOn = .false.
     integer :: firstComponent = 1, secondComponent = 2

   contains

     procedure, pass :: setup => setupLighthillTensorComponent
     procedure, pass :: cleanup => cleanupLighthillTensorComponent
     procedure, pass :: compute => computeLighthillTensorComponent
     procedure, pass :: computeSpatialDistribution => computeLighthillTensorComponentSpatialDistribution
     procedure, pass :: computeAdjointForcing => computeLighthillTensorComponentAdjointForcing
     procedure, pass :: isPatchValid => isLighthillTensorComponentPatchValid

  end type t_LighthillTensorComponent

  interface

     subroutine setupLighthillTensorComponent(this, region)

       use Region_mod, only : t_Region

       import :: t_LighthillTensorComponent

       class(t_LighthillTensorComponent) :: this
       class(t_Region) :: region

     end subroutine setupLighthillTensorComponent

  end interface

  interface

     subroutine cleanupLighthillTensorComponent(this)

       import :: t_LighthillTensorComponent

       class(t_LighthillTensorComponent) :: this

     end subroutine cleanupLighthillTensorComponent

  end interface

  interface

     function computeLighthillTensorComponent(this, region) result(instantaneousFunctional)

       use Region_mod, only : t_Region

       import :: t_LighthillTensorComponent

       class(t_LighthillTensorComponent) :: this
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousFunctional

     end function computeLighthillTensorComponent

  end interface

  interface

     subroutine computeLighthillTensorComponentSpatialDistribution(this, grid, state, F)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State

       import :: t_LighthillTensorComponent

       class(t_LighthillTensorComponent) :: this
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       SCALAR_TYPE, intent(out) :: F(:,:)

     end subroutine computeLighthillTensorComponentSpatialDistribution

  end interface

  interface

     subroutine computeLighthillTensorComponentAdjointForcing(this, simulationFlags,                    &
          solverOptions, grid, state, patch)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use CostTargetPatch_mod, only : t_CostTargetPatch
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_LighthillTensorComponent

       class(t_LighthillTensorComponent) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       class(t_CostTargetPatch) :: patch

     end subroutine computeLighthillTensorComponentAdjointForcing

  end interface

  interface

     function isLighthillTensorComponentPatchValid(this, patchDescriptor, gridSize,                     &
          normalDirection, extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_LighthillTensorComponent

       class(t_LighthillTensorComponent) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isLighthillTensorComponentPatchValid

  end interface

end module LighthillTensorComponent_mod
