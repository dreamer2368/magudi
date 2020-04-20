#include "config.h"

module GenericActuator_mod

  use Controller_mod, only : t_Controller

  implicit none

  type, extends(t_Controller), public :: t_GenericActuator

   contains

     procedure, pass :: setup => setupGenericActuator
     procedure, pass :: cleanup => cleanupGenericActuator
     procedure, pass :: computeSensitivity => computeGenericActuatorSensitivity
     procedure, pass :: updateForcing => updateGenericActuatorForcing
     procedure, pass :: updateDeltaForcing => updateGenericActuatorDeltaForcing
     procedure, pass :: migrateToForcing => migrateToGenericActuatorForcing
     procedure, pass :: updateGradient => updateGenericActuatorGradient
     procedure, pass :: isPatchValid => isGenericActuatorPatchValid
     procedure, pass :: hookBeforeTimemarch => hookGenericActuatorBeforeTimemarch
     procedure, pass :: hookAfterTimemarch => hookGenericActuatorAfterTimemarch

     procedure, pass :: collectNorm => collectGenericActuatorNorm

  end type t_GenericActuator

  interface

     subroutine setupGenericActuator(this, region)

       use Region_mod, only : t_Region

       import :: t_GenericActuator

       class(t_GenericActuator) :: this
       class(t_Region) :: region

     end subroutine setupGenericActuator

  end interface

  interface

     subroutine cleanupGenericActuator(this)

       import :: t_GenericActuator

       class(t_GenericActuator) :: this

     end subroutine cleanupGenericActuator

  end interface

  interface

     function computeGenericActuatorSensitivity(this, region) result(instantaneousSensitivity)

       use Region_mod, only : t_Region

       import :: t_GenericActuator

       class(t_GenericActuator) :: this
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousSensitivity

     end function computeGenericActuatorSensitivity

  end interface

  interface

     subroutine updateGenericActuatorForcing(this, region)

       use Region_mod, only : t_Region

       import :: t_GenericActuator

       class(t_GenericActuator) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateGenericActuatorForcing

  end interface

  interface

     subroutine updateGenericActuatorDeltaForcing(this, region)

       use Region_mod, only : t_Region

       import :: t_GenericActuator

       class(t_GenericActuator) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateGenericActuatorDeltaForcing

  end interface

  interface

     subroutine migrateToGenericActuatorForcing(this, region, startTimeStep, endTimeStep, nStages, iTimeStep, jSubStep)

       use Region_mod, only : t_Region

       import :: t_GenericActuator

       class(t_GenericActuator) :: this
       class(t_Region), intent(in) :: region
       integer, intent(in) :: startTimeStep, endTimeStep, nStages, iTimeStep, jSubStep

     end subroutine migrateToGenericActuatorForcing

  end interface

  interface

     subroutine updateGenericActuatorGradient(this, region)

       use Region_mod, only : t_Region

       import :: t_GenericActuator

       class(t_GenericActuator) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateGenericActuatorGradient

  end interface

  interface

     function isGenericActuatorPatchValid(this, patchDescriptor, gridSize,                   &
          normalDirection, extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_GenericActuator

       class(t_GenericActuator) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isGenericActuatorPatchValid

  end interface

  interface

     subroutine hookGenericActuatorBeforeTimemarch(this, region, mode, referenceTimestep)

       use Region_mod, only : t_Region

       import :: t_GenericActuator

       class(t_GenericActuator) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode
       integer, intent(in), optional :: referenceTimestep

     end subroutine hookGenericActuatorBeforeTimemarch

  end interface

  interface

     subroutine hookGenericActuatorAfterTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_GenericActuator

       class(t_GenericActuator) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookGenericActuatorAfterTimemarch

  end interface

  interface

     subroutine collectGenericActuatorNorm(this, region, timeIntegrationNorm)

       use Region_mod, only : t_Region

       import :: t_GenericActuator

       class(t_GenericActuator) :: this
       class(t_Region), intent(in) :: region
       real(SCALAR_KIND), intent(in) :: timeIntegrationNorm

     end subroutine collectGenericActuatorNorm

  end interface

end module GenericActuator_mod
