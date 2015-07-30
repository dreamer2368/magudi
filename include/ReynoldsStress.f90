#include "config.h"

module ReynoldsStress_mod

  use Functional_mod, only : t_Functional

  implicit none

  type, private :: t_ReynoldsStressInternal
     SCALAR_TYPE, pointer :: meanVelocity(:,:) => null()
  end type t_ReynoldsStressInternal  

  type, extends(t_Functional), public :: t_ReynoldsStress

     type(t_ReynoldsStressInternal), allocatable :: data_(:)
     real(SCALAR_KIND) :: firstDirection(3), secondDirection(3)

   contains

     procedure, pass :: setup => setupReynoldsStress
     procedure, pass :: cleanup => cleanupReynoldsStress
     procedure, pass :: compute => computeReynoldsStress
     procedure, pass :: computeAdjointForcing => computeReynoldsStressAdjointForcing
     procedure, pass :: isPatchValid => isReynoldsStressPatchValid

  end type t_ReynoldsStress

  interface

     subroutine setupReynoldsStress(this, region)

       use Region_mod, only : t_Region

       import :: t_ReynoldsStress

       class(t_ReynoldsStress) :: this
       class(t_Region) :: region

     end subroutine setupReynoldsStress

  end interface

  interface

     subroutine cleanupReynoldsStress(this)

       import :: t_ReynoldsStress

       class(t_ReynoldsStress) :: this

     end subroutine cleanupReynoldsStress

  end interface

  interface

     function computeReynoldsStress(this, region) result(instantaneousFunctional)

       use Region_mod, only : t_Region

       import :: t_ReynoldsStress

       class(t_ReynoldsStress) :: this
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousFunctional

     end function computeReynoldsStress

  end interface

  interface

     subroutine computeReynoldsStressAdjointForcing(this, simulationFlags,                   &
          solverOptions, grid, state, patch)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use CostTargetPatch_mod, only : t_CostTargetPatch
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ReynoldsStress

       class(t_ReynoldsStress) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       class(t_CostTargetPatch) :: patch

     end subroutine computeReynoldsStressAdjointForcing

  end interface

  interface

     function isReynoldsStressPatchValid(this, patchDescriptor, gridSize, normalDirection,   &
          extent, simulationFlags, message) result(isPatchValid)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ReynoldsStress

       class(t_ReynoldsStress) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isReynoldsStressPatchValid

  end interface

end module ReynoldsStress_mod
