#include "config.h"

module DensityGradient_mod

  use Functional_mod, only : t_Functional

  implicit none

  type, private :: t_DensityGradientInternal
     SCALAR_TYPE, pointer :: adjointVector(:,:) => null()
  end type t_DensityGradientInternal

  type, extends(t_Functional), public :: t_DensityGradient

     type(t_DensityGradientInternal), allocatable :: data_(:)
     SCALAR_TYPE :: timeWindowCenter, timeWindowWidth
     logical :: useTimeWindow = .false.
     integer :: secondComponent

   contains

     procedure, pass :: setup => setupDensityGradient
     procedure, pass :: cleanup => cleanupDensityGradient
     procedure, pass :: compute => computeDensityGradient
     procedure, pass :: computeSpatialDistribution => computeDensityGradientSpatialDistribution
     procedure, pass :: computeAdjointForcing => computeDensityGradientAdjointForcing
     procedure, pass :: isPatchValid => isDensityGradientPatchValid

  end type t_DensityGradient

  interface

     subroutine setupDensityGradient(this, region)

       use Region_mod, only : t_Region

       import :: t_DensityGradient

       class(t_DensityGradient) :: this
       class(t_Region) :: region

     end subroutine setupDensityGradient

  end interface

  interface

     subroutine cleanupDensityGradient(this)

       import :: t_DensityGradient

       class(t_DensityGradient) :: this

     end subroutine cleanupDensityGradient

  end interface

  interface

     function computeDensityGradient(this, region) result(instantaneousFunctional)

       use Region_mod, only : t_Region

       import :: t_DensityGradient

       class(t_DensityGradient) :: this
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousFunctional

     end function computeDensityGradient

  end interface

  interface

     subroutine computeDensityGradientSpatialDistribution(this, grid, state, F)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State

       import :: t_DensityGradient

       class(t_DensityGradient) :: this
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       SCALAR_TYPE, intent(out) :: F(:,:)

     end subroutine computeDensityGradientSpatialDistribution

  end interface

  interface

     subroutine computeDensityGradientAdjointForcing(this, simulationFlags,                    &
          solverOptions, grid, state, patch)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use CostTargetPatch_mod, only : t_CostTargetPatch
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_DensityGradient

       class(t_DensityGradient) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       class(t_CostTargetPatch) :: patch

     end subroutine computeDensityGradientAdjointForcing

  end interface

  interface

     function isDensityGradientPatchValid(this, patchDescriptor, gridSize,                     &
          normalDirection, extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_DensityGradient

       class(t_DensityGradient) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isDensityGradientPatchValid

  end interface

end module DensityGradient_mod
