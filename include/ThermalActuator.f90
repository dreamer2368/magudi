#include "config.h"

module ThermalActuator_mod

  use Controller_mod, only : t_Controller

  implicit none
  private

  type, extends(t_Controller), public :: t_ThermalActuator

   contains

     procedure, pass :: setup => setupThermalActuator
     procedure, pass :: cleanup => cleanupThermalActuator
     procedure, pass :: computeSensitivity => computeThermalActuatorSensitivity
     procedure, pass :: computeGradient => computeThermalActuatorGradient
     procedure, pass :: updateForcing => updateThermalActuatorForcing
     procedure, pass :: updateGradient => updateThermalActuatorGradient
     procedure, pass :: isPatchValid => isThermalActuatorPatchValid
     procedure, pass :: hookBeforeTimemarch => hookThermalActuatorBeforeTimemarch
     procedure, pass :: hookAfterTimemarch => hookThermalActuatorAfterTimemarch
     procedure, pass :: addPenalty => addThermalPenalty

  end type t_ThermalActuator

  interface

     subroutine addThermalPenalty(this,cost,region, mode)

       use Region_mod, only : t_Region

       import :: t_ThermalActuator

       class(t_ThermalActuator) :: this
       class(t_Region) :: region
       SCALAR_TYPE::cost
       integer, intent(in) :: mode
     end subroutine

     subroutine setupThermalActuator(this, region)

       use Region_mod, only : t_Region

       import :: t_ThermalActuator

       class(t_ThermalActuator) :: this
       class(t_Region) :: region

     end subroutine setupThermalActuator

  end interface

  interface

     subroutine cleanupThermalActuator(this)

       import :: t_ThermalActuator

       class(t_ThermalActuator) :: this

     end subroutine cleanupThermalActuator

  end interface

interface
subroutine computeThermalActuatorSensitivity(this,timeIntegrator, region) 

use Region_mod, only : t_Region
use TimeIntegrator_mod, only : t_TimeIntegrator
import :: t_ThermalActuator

class(t_ThermalActuator) :: this
class(t_Region), intent(in) :: region
class(t_TimeIntegrator),intent(in) :: timeIntegrator

end subroutine computeThermalActuatorSensitivity
end interface

interface
subroutine computeThermalActuatorGradient(this,timeIntegrator,region) 

use Region_mod, only : t_Region
use TimeIntegrator_mod, only : t_TimeIntegrator
import :: t_ThermalActuator

class(t_ThermalActuator) :: this
class(t_Region), intent(in) :: region
class(t_TimeIntegrator),intent(in) :: timeIntegrator

end subroutine computeThermalActuatorGradient
end interface

  interface

     subroutine updateThermalActuatorForcing(this, region)

       use Region_mod, only : t_Region

       import :: t_ThermalActuator

       class(t_ThermalActuator) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateThermalActuatorForcing

  end interface

  interface

     subroutine updateThermalActuatorGradient(this, region)

       use Region_mod, only : t_Region

       import :: t_ThermalActuator

       class(t_ThermalActuator) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateThermalActuatorGradient

  end interface

  interface

     function isThermalActuatorPatchValid(this, patchDescriptor, gridSize,                   &
          normalDirection, extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ThermalActuator

       class(t_ThermalActuator) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isThermalActuatorPatchValid

  end interface

  interface

     subroutine hookThermalActuatorBeforeTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_ThermalActuator

       class(t_ThermalActuator) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookThermalActuatorBeforeTimemarch

  end interface

  interface

     subroutine hookThermalActuatorAfterTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_ThermalActuator

       class(t_ThermalActuator) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookThermalActuatorAfterTimemarch

  end interface

end module ThermalActuator_mod
