#include "config.h"

module IsothermalWall_mod

  use ImpenetrableWall_mod, only : t_ImpenetrableWall

  implicit none
  private

  type, extends(t_ImpenetrableWall), public :: t_IsothermalWall

     real(SCALAR_KIND) :: viscousPenaltyAmounts(2)
     SCALAR_TYPE, allocatable :: temperature(:), dynamicViscosity(:),                        &
          secondCoefficientOfViscosity(:), thermalDiffusivity(:)

   contains

     procedure, pass :: setup => setupIsothermalWall
     procedure, pass :: cleanup => cleanupIsothermalWall
     procedure, pass :: update => updateIsothermalWall
     procedure, pass :: verifyUsage => verifyIsothermalWallUsage
     procedure, pass :: updateRhs => addIsothermalWallPenalty

  end type t_IsothermalWall

  interface

     subroutine setupIsothermalWall(this, index, comm, patchDescriptor,                      &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_IsothermalWall

       class(t_IsothermalWall) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupIsothermalWall

  end interface

  interface

     subroutine cleanupIsothermalWall(this)

       import :: t_IsothermalWall

       class(t_IsothermalWall) :: this

     end subroutine cleanupIsothermalWall

  end interface

  interface

     subroutine updateIsothermalWall(this, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_IsothermalWall

       class(t_IsothermalWall) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state

     end subroutine updateIsothermalWall

  end interface

  interface

     function verifyIsothermalWallUsage(this, patchDescriptor, gridSize, normalDirection,    &
          extent, simulationFlags, success, message) result(isPatchUsed)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_IsothermalWall

       class(t_IsothermalWall) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       logical, intent(out) :: success
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchUsed

     end function verifyIsothermalWallUsage

  end interface

  interface

     subroutine addIsothermalWallPenalty(this, mode, simulationFlags,                        &
          solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_IsothermalWall

       class(t_IsothermalWall) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine addIsothermalWallPenalty

  end interface

end module IsothermalWall_mod
