#include "config.h"

module WallActuator_mod

  use Controller_mod, only : t_Controller

  implicit none
  private

  type, extends(t_Controller), public :: t_WallActuator
 
     integer::numP
     integer::controlIndex 
     real(SCALAR_KIND),allocatable,dimension(:)::p  !wall-shape parameterization
     real(SCALAR_KIND),allocatable,dimension(:)::instantaneousGradient    
     real(SCALAR_KIND),allocatable,dimension(:)::gradient
     !real(SCALAR_KIND)::sensitivity
 
     !this may be getting moved to their own patch
     !at this point I am using the controller to know everything about grid
     
     real(SCALAR_KIND),allocatable,dimension(:,:)::dJacobiandp    
     real(SCALAR_KIND),allocatable,dimension(:,:,:)::dMijdp
   
   contains

     procedure, pass :: setup => setupWallActuator
     procedure, pass :: cleanup => cleanupWallActuator
     procedure, pass :: computeSensitivity => computeWallActuatorSensitivity
     procedure, pass :: computeGradient => computeWallActuatorGradient
     procedure, pass :: updateForcing => updateWallActuatorForcing
     procedure, pass :: updateGradient => updateWallActuatorGradient
     procedure, pass :: isPatchValid => isWallActuatorPatchValid
     procedure, pass :: hookBeforeTimemarch => hookWallActuatorBeforeTimemarch
     procedure, pass :: hookAfterTimemarch => hookWallActuatorAfterTimemarch
   
  end type t_WallActuator

  interface

     subroutine setupWallActuator(this, region)

       use Region_mod, only : t_Region

       import :: t_WallActuator

       class(t_WallActuator) :: this
       class(t_Region) :: region

     end subroutine setupWallActuator

  end interface

  interface

     subroutine cleanupWallActuator(this)

       import :: t_WallActuator

       class(t_WallActuator) :: this

     end subroutine cleanupWallActuator

  end interface

interface

subroutine computeWallActuatorSensitivity(this,timeIntegrator, region) 

  use Region_mod, only : t_Region
  use TimeIntegrator_mod, only : t_TimeIntegrator
  import :: t_WallActuator

  class(t_WallActuator) :: this
  class(t_Region), intent(in) :: region
  class(t_TimeIntegrator),intent(in) :: timeIntegrator


end subroutine computeWallActuatorSensitivity

end interface

interface

subroutine computeWallActuatorGradient(this,timeIntegrator,region) 

use Region_mod, only : t_Region
use TimeIntegrator_mod, only : t_TimeIntegrator
import :: t_WallActuator

class(t_WallActuator) :: this  
class(t_Region), intent(in) :: region
class(t_TimeIntegrator),intent(in) :: timeIntegrator

end subroutine computeWallActuatorGradient

end interface


  interface

     subroutine updateWallActuatorForcing(this, region)

       use Region_mod, only : t_Region

       import :: t_WallActuator

       class(t_WallActuator) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateWallActuatorForcing

  end interface

  interface

     subroutine updateWallActuatorGradient(this, region)

     use Region_mod, only : t_Region
  
     import :: t_WallActuator

       class(t_WallActuator) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateWallActuatorGradient

  end interface

  interface

     function isWallActuatorPatchValid(this, patchDescriptor, gridSize,                   &
          normalDirection, extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_WallActuator

       class(t_WallActuator) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isWallActuatorPatchValid

  end interface

  interface

     subroutine hookWallActuatorBeforeTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_WallActuator

       class(t_WallActuator) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookWallActuatorBeforeTimemarch

  end interface

  interface

     subroutine hookWallActuatorAfterTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_WallActuator

       class(t_WallActuator) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookWallActuatorAfterTimemarch

  end interface

end module WallActuator_mod
