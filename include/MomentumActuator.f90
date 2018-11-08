#include "config.h"

module MomentumActuator_mod

  use Controller_mod, only : t_Controller

  implicit none

  type, extends(t_Controller), public :: t_MomentumActuator

     integer :: direction, nActuatorComponents

   contains

     procedure, pass :: setup => setupMomentumActuator
     procedure, pass :: cleanup => cleanupMomentumActuator
     procedure, pass :: computeSensitivity => computeMomentumActuatorSensitivity
     procedure, pass :: updateForcing => updateMomentumActuatorForcing
     procedure, pass :: migrateToForcing => migrateToMomentumActuatorForcing
     procedure, pass :: updateGradient => updateMomentumActuatorGradient
     procedure, pass :: isPatchValid => isMomentumActuatorPatchValid
     procedure, pass :: hookBeforeTimemarch => hookMomentumActuatorBeforeTimemarch
     procedure, pass :: hookAfterTimemarch => hookMomentumActuatorAfterTimemarch

  end type t_MomentumActuator

  interface

     subroutine setupMomentumActuator(this, region)

       use Region_mod, only : t_Region

       import :: t_MomentumActuator

       class(t_MomentumActuator) :: this
       class(t_Region) :: region

     end subroutine setupMomentumActuator

  end interface

  interface

     subroutine cleanupMomentumActuator(this)

       import :: t_MomentumActuator

       class(t_MomentumActuator) :: this

     end subroutine cleanupMomentumActuator

  end interface

  interface

     function computeMomentumActuatorSensitivity(this, region) result(instantaneousSensitivity)

       use Region_mod, only : t_Region

       import :: t_MomentumActuator

       class(t_MomentumActuator) :: this
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousSensitivity

     end function computeMomentumActuatorSensitivity

  end interface

  interface

     subroutine updateMomentumActuatorForcing(this, region)

       use Region_mod, only : t_Region

       import :: t_MomentumActuator

       class(t_MomentumActuator) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateMomentumActuatorForcing

  end interface

  interface

     subroutine migrateToMomentumActuatorForcing(this, region, nTimeSteps, nStages, iTimeStep, jSubStep)

       use Region_mod, only : t_Region

       import :: t_MomentumActuator

       class(t_MomentumActuator) :: this
       class(t_Region), intent(in) :: region
       integer, intent(in) :: nTimeSteps, nStages, iTimeStep, jSubStep

     end subroutine migrateToMomentumActuatorForcing

  end interface

  interface

     subroutine updateMomentumActuatorGradient(this, region)

       use Region_mod, only : t_Region

       import :: t_MomentumActuator

       class(t_MomentumActuator) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateMomentumActuatorGradient

  end interface

  interface

     function isMomentumActuatorPatchValid(this, patchDescriptor, gridSize,                   &
          normalDirection, extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_MomentumActuator

       class(t_MomentumActuator) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isMomentumActuatorPatchValid

  end interface

  interface

     subroutine hookMomentumActuatorBeforeTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_MomentumActuator

       class(t_MomentumActuator) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookMomentumActuatorBeforeTimemarch

  end interface

  interface

     subroutine hookMomentumActuatorAfterTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_MomentumActuator

       class(t_MomentumActuator) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookMomentumActuatorAfterTimemarch

  end interface

end module MomentumActuator_mod
