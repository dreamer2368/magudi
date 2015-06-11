#include "config.h"

module FlameIgnition_mod

  use Controller_mod, only : t_Controller

  implicit none
  private

  type, extends(t_Controller), public :: t_FlameIgnition

     character(len = STRING_LENGTH) :: sensitivityDependence
     logical :: partialSensitivity

   contains

     procedure, pass :: setup => setupFlameIgnition
     procedure, pass :: cleanup => cleanupFlameIgnition
     procedure, pass :: computeSensitivity => computeFlameIgnitionSensitivity
     procedure, pass :: updateGradient => updateFlameIgnitionGradient
     procedure, pass :: updateForcing => updateFlameIgnitionForcing
     procedure, pass :: isPatchValid => isFlameIgnitionPatchValid
     procedure, pass :: hookBeforeTimemarch => hookFlameIgnitionBeforeTimemarch
     procedure, pass :: hookAfterTimemarch => hookFlameIgnitionAfterTimemarch


  end type t_FlameIgnition

  interface

     subroutine setupFlameIgnition(this, region)

       use Region_mod, only : t_Region

       import :: t_FlameIgnition

       class(t_FlameIgnition) :: this
       class(t_Region) :: region

     end subroutine setupFlameIgnition

  end interface

  interface

     subroutine cleanupFlameIgnition(this)

       import :: t_FlameIgnition

       class(t_FlameIgnition) :: this

     end subroutine cleanupFlameIgnition

  end interface

  interface

     function computeFlameIgnitionSensitivity(this, region) result(instantaneousSensitivity)

       use Region_mod, only : t_Region

       import :: t_FlameIgnition

       class(t_FlameIgnition) :: this
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousSensitivity

     end function computeFlameIgnitionSensitivity

  end interface

  interface

     subroutine updateFlameIgnitionForcing(this, region)

       use Region_mod, only : t_Region

       import :: t_FlameIgnition

       class(t_FlameIgnition) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateFlameIgnitionForcing

  end interface

  interface

     subroutine updateFlameIgnitionGradient(this, region)

       use Region_mod, only : t_Region

       import :: t_FlameIgnition

       class(t_FlameIgnition) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateFlameIgnitionGradient

  end interface

  interface

     function isFlameIgnitionPatchValid(this, patchDescriptor, gridSize,                   &
          normalDirection, extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_FlameIgnition

       class(t_FlameIgnition) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isFlameIgnitionPatchValid

  end interface

  interface

     subroutine hookFlameIgnitionBeforeTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_FlameIgnition

       class(t_FlameIgnition) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookFlameIgnitionBeforeTimemarch

  end interface

  interface

     subroutine hookFlameIgnitionAfterTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_FlameIgnition

       class(t_FlameIgnition) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookFlameIgnitionAfterTimemarch

  end interface

end module FlameIgnition_mod
