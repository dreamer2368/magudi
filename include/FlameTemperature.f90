#include "config.h"

module FlameTemperature_mod

  use Functional_mod, only : t_Functional

  implicit none
  private

  type, private :: t_FlameTemperatureInternal
     SCALAR_TYPE, pointer :: flameTemperature(:,:) => null()
  end type t_FlameTemperatureInternal

  type, extends(t_Functional), public :: t_FlameTemperature

     type(t_FlameTemperatureInternal), allocatable :: data_(:)

   contains

     procedure, pass :: setup => setupFlameTemperature
     procedure, pass :: cleanup => cleanupFlameTemperature
     procedure, pass :: compute => computeFlameTemperature
     procedure, pass :: computeAdjointForcing => computeFlameTemperatureAdjointForcing
     procedure, pass :: isPatchValid => isFlameTemperaturePatchValid

  end type t_FlameTemperature

  interface

     subroutine setupFlameTemperature(this, region)

       use Region_mod, only : t_Region

       import :: t_FlameTemperature

       class(t_FlameTemperature) :: this
       class(t_Region) :: region

     end subroutine setupFlameTemperature

  end interface

  interface

     subroutine cleanupFlameTemperature(this)

       import :: t_FlameTemperature

       class(t_FlameTemperature) :: this

     end subroutine cleanupFlameTemperature

  end interface

  interface

     function computeFlameTemperature(this, region) result(instantaneousFunctional)

       use Region_mod, only : t_Region

       import :: t_FlameTemperature

       class(t_FlameTemperature) :: this
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousFunctional

     end function computeFlameTemperature

  end interface

  interface

     subroutine computeFlameTemperatureAdjointForcing(this, simulationFlags,                    &
          solverOptions, grid, state, patch)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use CostTargetPatch_mod, only : t_CostTargetPatch
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_FlameTemperature

       class(t_FlameTemperature) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       class(t_CostTargetPatch) :: patch

     end subroutine computeFlameTemperatureAdjointForcing

  end interface

  interface

     function isFlameTemperaturePatchValid(this, patchDescriptor, gridSize,                     &
          normalDirection, extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_FlameTemperature

       class(t_FlameTemperature) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isFlameTemperaturePatchValid

  end interface

end module FlameTemperature_mod
