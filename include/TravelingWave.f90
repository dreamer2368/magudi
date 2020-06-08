#include "config.h"

module TravelingWave_mod

  use Functional_mod, only : t_Functional

  implicit none

  type, extends(t_Functional), public :: t_TravelingWave

   contains

     procedure, pass :: setup => setupTravelingWave
     procedure, pass :: cleanup => cleanupTravelingWave
     procedure, pass :: compute => computeTravelingWave
     procedure, pass :: computeSpatialDistribution => computeTravelingWaveSpatialDistribution
     procedure, pass :: computeAdjointForcing => computeTravelingWaveAdjointForcing
     procedure, pass :: isPatchValid => isTravelingWavePatchValid
     procedure, pass :: addGradient => addTravelingWaveGradient

  end type t_TravelingWave

  interface

     subroutine setupTravelingWave(this, region)

       use Region_mod, only : t_Region

       import :: t_TravelingWave

       class(t_TravelingWave) :: this
       class(t_Region) :: region

     end subroutine setupTravelingWave

  end interface

  interface

     subroutine cleanupTravelingWave(this)

       import :: t_TravelingWave

       class(t_TravelingWave) :: this

     end subroutine cleanupTravelingWave

  end interface

  interface

     function computeTravelingWave(this, region) result(instantaneousFunctional)

       use Region_mod, only : t_Region

       import :: t_TravelingWave

       class(t_TravelingWave) :: this
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousFunctional

     end function computeTravelingWave

  end interface

  interface

     subroutine computeTravelingWaveSpatialDistribution(this, grid, state, F)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State

       import :: t_TravelingWave

       class(t_TravelingWave) :: this
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       SCALAR_TYPE, intent(out) :: F(:,:)

     end subroutine computeTravelingWaveSpatialDistribution

  end interface

  interface

     subroutine computeTravelingWaveAdjointForcing(this, simulationFlags,                    &
          solverOptions, grid, state, patch)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use CostTargetPatch_mod, only : t_CostTargetPatch
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_TravelingWave

       class(t_TravelingWave) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       class(t_CostTargetPatch) :: patch

     end subroutine computeTravelingWaveAdjointForcing

  end interface

  interface

     subroutine addTravelingWaveGradient(this, simulationFlags,                    &
          solverOptions, grid, state, patch)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use CostTargetPatch_mod, only : t_CostTargetPatch
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_TravelingWave

       class(t_TravelingWave) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       class(t_CostTargetPatch) :: patch

     end subroutine addTravelingWaveGradient

  end interface

  interface

     function isTravelingWavePatchValid(this, patchDescriptor, gridSize,                     &
          normalDirection, extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_TravelingWave

       class(t_TravelingWave) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isTravelingWavePatchValid

  end interface

end module TravelingWave_mod
