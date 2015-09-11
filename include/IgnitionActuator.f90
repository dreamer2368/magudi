#include "config.h"

module IgnitionActuator_mod

  use Controller_mod, only : t_Controller

  implicit none

  type, extends(t_Controller), public :: t_IgnitionActuator

     real(SCALAR_KIND), dimension(:), allocatable :: baselineValue
     character(len = STRING_LENGTH), dimension(:), allocatable :: sensitivityParameter

   contains

     procedure, pass :: setup => setupIgnitionActuator
     procedure, pass :: cleanup => cleanupIgnitionActuator
     procedure, pass :: computeSensitivity => computeIgnitionActuatorSensitivity
     procedure, pass :: updateGradient => updateIgnitionActuatorGradient
     procedure, pass :: updateForcing => updateIgnitionActuatorForcing
     procedure, pass :: isPatchValid => isIgnitionActuatorPatchValid
     procedure, pass :: hookBeforeTimemarch => hookIgnitionActuatorBeforeTimemarch
     procedure, pass :: hookAfterTimemarch => hookIgnitionActuatorAfterTimemarch


  end type t_IgnitionActuator

  interface

     subroutine setupIgnitionActuator(this, region)

       use Region_mod, only : t_Region

       import :: t_IgnitionActuator

       class(t_IgnitionActuator) :: this
       class(t_Region) :: region

     end subroutine setupIgnitionActuator

  end interface

  interface

     subroutine cleanupIgnitionActuator(this)

       import :: t_IgnitionActuator

       class(t_IgnitionActuator) :: this

     end subroutine cleanupIgnitionActuator

  end interface

  interface

     subroutine computeIgnitionActuatorSensitivity(this, region)

       use Region_mod, only : t_Region

       import :: t_IgnitionActuator

       class(t_IgnitionActuator) :: this
       class(t_Region) :: region

     end subroutine computeIgnitionActuatorSensitivity

  end interface

  interface

     subroutine updateIgnitionActuatorForcing(this, region)

       use Region_mod, only : t_Region

       import :: t_IgnitionActuator

       class(t_IgnitionActuator) :: this
       class(t_Region) :: region

     end subroutine updateIgnitionActuatorForcing

  end interface

  interface

     subroutine updateIgnitionActuatorGradient(this, region)

       use Region_mod, only : t_Region

       import :: t_IgnitionActuator

       class(t_IgnitionActuator) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateIgnitionActuatorGradient

  end interface

  interface

     function isIgnitionActuatorPatchValid(this, patchDescriptor, gridSize,                  &
          normalDirection, extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_IgnitionActuator

       class(t_IgnitionActuator) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isIgnitionActuatorPatchValid

  end interface

  interface

     subroutine hookIgnitionActuatorBeforeTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_IgnitionActuator

       class(t_IgnitionActuator) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookIgnitionActuatorBeforeTimemarch

  end interface

  interface

     subroutine hookIgnitionActuatorAfterTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_IgnitionActuator

       class(t_IgnitionActuator) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookIgnitionActuatorAfterTimemarch

  end interface

end module IgnitionActuator_mod
