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
     procedure, pass :: addAdjointForcing => addAcousticNoiseAdjointForcing

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
     
     function computeAcousticNoise(this, region) result(instantaneousFunctional)

       use Region_mod, only : t_Region
       
       import :: t_AcousticNoise
       
       class(t_AcousticNoise) :: this
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousFunctional
       
     end function computeAcousticNoise
     
  end interface

  interface

     subroutine addAcousticNoiseAdjointForcing(this, simulationFlags,                        &
          solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_AcousticNoise

       class(t_AcousticNoise) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine addAcousticNoiseAdjointForcing

  end interface

end module AcousticNoise_mod
