#include "config.h"

module LighthillSource_mod

  use Functional_mod, only : t_Functional

  implicit none

  type, private :: t_LighthillSourceInternal
     SCALAR_TYPE, pointer :: adjointVector(:,:) => null()
  end type t_LighthillSourceInternal

  type, extends(t_Functional), public :: t_LighthillSource

     type(t_LighthillSourceInternal), allocatable :: data_(:)
     SCALAR_TYPE :: timeWindowCenter, timeWindowWidth
     logical :: useTimeWindow = .false., viscosityOn = .false.
     real(SCALAR_KIND) :: firstDirection(3), secondDirection(3)
     integer :: firstComponent, secondComponent

   contains

     procedure, pass :: setup => setupLighthillSource
     procedure, pass :: cleanup => cleanupLighthillSource
     procedure, pass :: compute => computeLighthillSource
     procedure, pass :: computeSpatialDistribution => computeLighthillSourceSpatialDistribution
     procedure, pass :: computeAdjointForcing => computeLighthillSourceAdjointForcing
     procedure, pass :: isPatchValid => isLighthillSourcePatchValid

  end type t_LighthillSource

  interface

     subroutine setupLighthillSource(this, region)

       use Region_mod, only : t_Region

       import :: t_LighthillSource

       class(t_LighthillSource) :: this
       class(t_Region) :: region

     end subroutine setupLighthillSource

  end interface

  interface

     subroutine cleanupLighthillSource(this)

       import :: t_LighthillSource

       class(t_LighthillSource) :: this

     end subroutine cleanupLighthillSource

  end interface

  interface

     function computeLighthillSource(this, region) result(instantaneousFunctional)

       use Region_mod, only : t_Region

       import :: t_LighthillSource

       class(t_LighthillSource) :: this
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousFunctional

     end function computeLighthillSource

  end interface

  interface

     subroutine computeLighthillSourceSpatialDistribution(this, grid, state, F)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State

       import :: t_LighthillSource

       class(t_LighthillSource) :: this
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       SCALAR_TYPE, intent(out) :: F(:,:)

     end subroutine computeLighthillSourceSpatialDistribution

  end interface

  interface

     subroutine computeLighthillSourceAdjointForcing(this, simulationFlags,                    &
          solverOptions, grid, state, patch)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use CostTargetPatch_mod, only : t_CostTargetPatch
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_LighthillSource

       class(t_LighthillSource) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       class(t_CostTargetPatch) :: patch

     end subroutine computeLighthillSourceAdjointForcing

  end interface

  interface

     function isLighthillSourcePatchValid(this, patchDescriptor, gridSize,                     &
          normalDirection, extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_LighthillSource

       class(t_LighthillSource) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isLighthillSourcePatchValid

  end interface

end module LighthillSource_mod
