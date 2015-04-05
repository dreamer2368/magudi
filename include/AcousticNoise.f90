#include "config.h"

module AcousticNoise_mod

  use Functional_mod, only : t_Functional

  implicit none
  private

  type, private :: t_AcousticNoiseInternal
     SCALAR_TYPE, pointer :: meanPressure(:,:) => null()
  end type t_AcousticNoiseInternal

  type, extends(t_Functional), public :: t_AcousticNoise

     type(t_AcousticNoiseInternal), allocatable :: data_(:)

   contains

     procedure, pass :: setup => setupAcousticNoise
     procedure, pass :: cleanup => cleanupAcousticNoise
     procedure, pass :: compute => computeAcousticNoise
     procedure, pass :: computeAdjointForcing => computeAcousticNoiseAdjointForcing
     procedure, pass :: isPatchValid => isAcousticNoisePatchValid

  end type t_AcousticNoise

  interface

     subroutine setupAcousticNoise(this, region)

       use Region_mod, only : t_Region

       import :: t_AcousticNoise

       class(t_AcousticNoise) :: this
       class(t_Region) :: region

     end subroutine setupAcousticNoise

  end interface

  interface

     subroutine cleanupAcousticNoise(this)

       import :: t_AcousticNoise

       class(t_AcousticNoise) :: this

     end subroutine cleanupAcousticNoise

  end interface

  interface

     function computeAcousticNoise(this, time, region) result(instantaneousFunctional)

       use Region_mod, only : t_Region

       import :: t_AcousticNoise

       class(t_AcousticNoise) :: this
       real(SCALAR_KIND), intent(in) :: time
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousFunctional

     end function computeAcousticNoise

  end interface

  interface

     subroutine computeAcousticNoiseAdjointForcing(this, region)

       use Region_mod, only : t_Region

       import :: t_AcousticNoise

       class(t_AcousticNoise) :: this
       class(t_Region) :: region

     end subroutine computeAcousticNoiseAdjointForcing

  end interface

  interface

     function isAcousticNoisePatchValid(this, patchDescriptor, gridSize,                     &
          normalDirection, extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_AcousticNoise

       class(t_AcousticNoise) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isAcousticNoisePatchValid

  end interface

end module AcousticNoise_mod