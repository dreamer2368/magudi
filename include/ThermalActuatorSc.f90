#include "config.h"

module ThermalActuatorSc_mod

  use Controller_mod, only : t_Controller

  implicit none
  private

  type, extends(t_Controller), public :: t_ThermalActuatorSc

   contains

     procedure, pass :: setup => setupThermalActuatorSc
     procedure, pass :: cleanup => cleanupThermalActuatorSc
     procedure, pass :: computeSensitivity => computeThermalActuatorScSensitivity
     procedure, pass :: updateForcing => updateThermalActuatorScForcing
     procedure, pass :: updateGradient => updateThermalActuatorScGradient
     procedure, pass :: isPatchValid => isThermalActuatorScPatchValid
     procedure, pass :: hookBeforeTimemarch => hookThermalActuatorScBeforeTimemarch
     procedure, pass :: hookAfterTimemarch => hookThermalActuatorScAfterTimemarch


  end type t_ThermalActuatorSc

  interface

     subroutine setupThermalActuatorSc(this, region)

       use Region_mod, only : t_Region

       import :: t_ThermalActuatorSc

       class(t_ThermalActuatorSc) :: this
       class(t_Region) :: region

     end subroutine setupThermalActuatorSc

  end interface

  interface

     subroutine cleanupThermalActuatorSc(this)

       import :: t_ThermalActuatorSc

       class(t_ThermalActuatorSc) :: this

     end subroutine cleanupThermalActuatorSc

  end interface

  interface

     function computeThermalActuatorScSensitivity(this, region) result(instantaneousSensitivity)

       use Region_mod, only : t_Region

       import :: t_ThermalActuatorSc

       class(t_ThermalActuatorSc) :: this
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousSensitivity

     end function computeThermalActuatorScSensitivity

  end interface

  interface

     subroutine updateThermalActuatorScForcing(this, region)

       use Region_mod, only : t_Region

       import :: t_ThermalActuatorSc

       class(t_ThermalActuatorSc) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateThermalActuatorScForcing

  end interface

  interface

     subroutine updateThermalActuatorScGradient(this, region)

       use Region_mod, only : t_Region

       import :: t_ThermalActuatorSc

       class(t_ThermalActuatorSc) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateThermalActuatorScGradient

  end interface

  interface

     function isThermalActuatorScPatchValid(this, patchDescriptor, gridSize,                   &
          normalDirection, extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ThermalActuatorSc

       class(t_ThermalActuatorSc) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isThermalActuatorScPatchValid

  end interface

  interface

     subroutine hookThermalActuatorScBeforeTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_ThermalActuatorSc

       class(t_ThermalActuatorSc) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookThermalActuatorScBeforeTimemarch

  end interface

  interface

     subroutine hookThermalActuatorScAfterTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_ThermalActuatorSc

       class(t_ThermalActuatorSc) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookThermalActuatorScAfterTimemarch

  end interface

end module ThermalActuatorSc_mod
