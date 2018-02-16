#include "config.h"

module Controller_mod

  implicit none

  type, abstract, public :: t_Controller

     SCALAR_TYPE :: cachedValue = real(0.0, SCALAR_KIND),                                    &
          runningTimeQuadrature = real(0.0, SCALAR_KIND)
     real(SCALAR_KIND) :: onsetTime = real(0.0, SCALAR_KIND),                                &
          duration = real(0.0, SCALAR_KIND)

   contains

     procedure, non_overridable, pass :: setupBase => setupController
     procedure, non_overridable, pass :: cleanupBase => cleanupController
     procedure, pass :: writeSensitivityToFile

     procedure(setup), pass, deferred :: setup
     procedure(cleanup), pass, deferred :: cleanup
     procedure(computeSensitivity), pass, deferred :: computeSensitivity
     procedure, non_overridable, pass :: cleanupForcing => cleanUpControlForcing
     procedure(updateForcing), pass, deferred :: updateForcing
     procedure(updateForcing), pass, deferred :: updateBaseForcing
     procedure(updateGradient), pass, deferred :: updateGradient
     procedure(isPatchValid), pass, deferred :: isPatchValid
     procedure(hookBeforeTimemarch), pass, deferred :: hookBeforeTimemarch
     procedure(hookAfterTimemarch), pass, deferred :: hookAfterTimemarch

  end type t_Controller

  abstract interface

     subroutine setup(this, region)

       use Region_mod, only : t_Region

       import :: t_Controller

       class(t_Controller) :: this
       class(t_Region) :: region

     end subroutine setup

  end interface

  abstract interface

     subroutine cleanup(this)

       import :: t_Controller

       class(t_Controller) :: this

     end subroutine cleanup

  end interface

  abstract interface

     function computeSensitivity(this, region) result(instantaneousSensitivity)

       use Region_mod, only : t_Region

       import :: t_Controller

       class(t_Controller) :: this
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousSensitivity

     end function computeSensitivity

  end interface

  abstract interface

     subroutine updateForcing(this, region)

       use Region_mod, only : t_Region

       import :: t_Controller

       class(t_Controller) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateForcing

  end interface

  abstract interface

     subroutine updateGradient(this, region)

       use Region_mod, only : t_Region

       import :: t_Controller

       class(t_Controller) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateGradient

  end interface

  abstract interface

     function isPatchValid(this, patchDescriptor, gridSize, normalDirection,                 &
          extent, simulationFlags, message)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Controller

       class(t_Controller) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isPatchValid

  end interface

  abstract interface

     subroutine hookBeforeTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_Controller

       class(t_Controller) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookBeforeTimemarch

  end interface

  abstract interface

     subroutine hookAfterTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_Controller

       class(t_Controller) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookAfterTimemarch

  end interface

  interface

     subroutine setupController(this, simulationFlags, solverOptions)

       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Controller

       class(t_Controller) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupController

  end interface

  interface

     subroutine cleanupController(this)

       import :: t_Controller

       class(t_Controller) :: this

     end subroutine cleanupController

  end interface

  interface

     subroutine writeSensitivityToFile(this, comm, filename, timestep, time, append)

       import :: t_Controller

       class(t_Controller) :: this
       integer, intent(in) :: comm
       character(len = *), intent(in) :: filename
       integer, intent(in) :: timestep
       real(SCALAR_KIND), intent(in) :: time
       logical, intent(in), optional :: append

     end subroutine writeSensitivityToFile

  end interface

  interface
     !SeungWhan: clean up control forcing (not controller!)
     subroutine cleanUpControlForcing(this,region)
     
       ! <<< Derived types >>>
       use Patch_mod, only : t_Patch
       use Region_mod, only : t_Region
       use ActuatorPatch_mod, only : t_ActuatorPatch

       import :: t_Controller
     
       implicit none
     
       class(t_Controller) :: this
       class(t_Region) :: region

     end subroutine

  end interface

end module Controller_mod
