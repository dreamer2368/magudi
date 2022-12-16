#include "config.h"

module ImmersedBoundaryPatch_mod

  use Patch_mod, only : t_Patch

  implicit none

  type, extends(t_Patch), public :: t_ImmersedBoundaryPatch
    ! real(SCALAR_KIND), dimension(:,:), allocatable :: levelset, levelsetNormal, levelsetCurvature,      &
    !                                                   indicatorFunction, primitiveGridNorm
    real(SCALAR_KIND) :: ratioOfSpecificHeats, ibmTemperature
    real(SCALAR_KIND) :: dti, maxSpeed, dissipationAmount, ibmEpsilon
  contains
     procedure, pass :: setup => setupImmersedBoundaryPatch
     procedure, pass :: cleanup => cleanupImmersedBoundaryPatch
     procedure, pass :: verifyUsage => verifyImmersedBoundaryPatchUsage
     procedure, pass :: updateRhs => addImmersedBoundaryPenalty
  end type t_ImmersedBoundaryPatch

  interface
    subroutine setupImmersedBoundaryPatch(this, index, comm, patchDescriptor,                            &
                                          grid, simulationFlags, solverOptions)
      use Grid_mod, only : t_Grid
      use SolverOptions_mod, only : t_SolverOptions
      use PatchDescriptor_mod, only : t_PatchDescriptor
      use SimulationFlags_mod, only : t_SimulationFlags

      import :: t_ImmersedBoundaryPatch

      class(t_ImmersedBoundaryPatch) :: this
      integer, intent(in) :: index, comm
      type(t_PatchDescriptor), intent(in) :: patchDescriptor
      class(t_Grid), intent(in) :: grid
      type(t_SimulationFlags), intent(in) :: simulationFlags
      type(t_SolverOptions), intent(in) :: solverOptions
    end subroutine setupImmersedBoundaryPatch
  end interface

  interface
    subroutine cleanupImmersedBoundaryPatch(this)
      import :: t_ImmersedBoundaryPatch
      class(t_ImmersedBoundaryPatch) :: this
    end subroutine cleanupImmersedBoundaryPatch
  end interface

  interface
    function verifyImmersedBoundaryPatchUsage(this, patchDescriptor, gridSize, normalDirection,   &
                                              extent, simulationFlags, success, message) result(isPatchUsed)
      use PatchDescriptor_mod, only : t_PatchDescriptor
      use SimulationFlags_mod, only : t_SimulationFlags
      import :: t_ImmersedBoundaryPatch
      class(t_ImmersedBoundaryPatch) :: this
      type(t_PatchDescriptor), intent(in) :: patchDescriptor
      integer, intent(in) :: gridSize(:), normalDirection, extent(6)
      type(t_SimulationFlags), intent(in) :: simulationFlags
      logical, intent(out) :: success
      character(len = STRING_LENGTH), intent(out) :: message
      logical :: isPatchUsed
    end function verifyImmersedBoundaryPatchUsage
  end interface

  interface
    subroutine addImmersedBoundaryPenalty(this, mode, simulationFlags, solverOptions, grid, state)
      use Grid_mod, only : t_Grid
      use State_mod, only : t_State
      use SolverOptions_mod, only : t_SolverOptions
      use SimulationFlags_mod, only : t_SimulationFlags
      import :: t_ImmersedBoundaryPatch
      class(t_ImmersedBoundaryPatch) :: this
      integer, intent(in) :: mode
      type(t_SimulationFlags), intent(in) :: simulationFlags
      type(t_SolverOptions), intent(in) :: solverOptions
      class(t_Grid), intent(in) :: grid
      class(t_State) :: state
    end subroutine addImmersedBoundaryPenalty
  end interface

end module ImmersedBoundaryPatch_mod
