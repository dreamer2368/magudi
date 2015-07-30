#include "config.h"

module FuelActuator_mod

  use Controller_mod, only : t_Controller

  implicit none

  type, extends(t_Controller), public :: t_FuelActuator

     integer :: fuelIndex

   contains

     procedure, pass :: setup => setupFuelActuator
     procedure, pass :: cleanup => cleanupFuelActuator
     procedure, pass :: computeSensitivity => computeFuelActuatorSensitivity
     procedure, pass :: updateForcing => updateFuelActuatorForcing
     procedure, pass :: updateGradient => updateFuelActuatorGradient
     procedure, pass :: isPatchValid => isFuelActuatorPatchValid
     procedure, pass :: hookBeforeTimemarch => hookFuelActuatorBeforeTimemarch
     procedure, pass :: hookAfterTimemarch => hookFuelActuatorAfterTimemarch


  end type t_FuelActuator

  interface

     subroutine setupFuelActuator(this, region)

       use Region_mod, only : t_Region

       import :: t_FuelActuator

       class(t_FuelActuator) :: this
       class(t_Region) :: region

     end subroutine setupFuelActuator

  end interface

  interface

     subroutine cleanupFuelActuator(this)

       import :: t_FuelActuator

       class(t_FuelActuator) :: this

     end subroutine cleanupFuelActuator

  end interface

  interface

     function computeFuelActuatorSensitivity(this, region) result(instantaneousSensitivity)

       use Region_mod, only : t_Region

       import :: t_FuelActuator

       class(t_FuelActuator) :: this
       class(t_Region) :: region

       SCALAR_TYPE :: instantaneousSensitivity

     end function computeFuelActuatorSensitivity

  end interface

  interface

     subroutine updateFuelActuatorForcing(this, region)

       use Region_mod, only : t_Region

       import :: t_FuelActuator

       class(t_FuelActuator) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateFuelActuatorForcing

  end interface

  interface

     subroutine updateFuelActuatorGradient(this, region)

       use Region_mod, only : t_Region

       import :: t_FuelActuator

       class(t_FuelActuator) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateFuelActuatorGradient

  end interface

  interface

     function isFuelActuatorPatchValid(this, patchDescriptor, gridSize,                   &
          normalDirection, extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_FuelActuator

       class(t_FuelActuator) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isFuelActuatorPatchValid

  end interface

  interface

     subroutine hookFuelActuatorBeforeTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_FuelActuator

       class(t_FuelActuator) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookFuelActuatorBeforeTimemarch

  end interface

  interface

     subroutine hookFuelActuatorAfterTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_FuelActuator

       class(t_FuelActuator) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookFuelActuatorAfterTimemarch

  end interface

end module FuelActuator_mod
