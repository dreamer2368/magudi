#include "config.h"

module RhsHelper

  implicit none
  public

  interface

     subroutine computeRhsForward(simulationFlags, solverOptions, grid, state, patchFactories)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use Patch_factory, only : t_PatchFactory
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid) :: grid
       class(t_State) :: state
       type(t_PatchFactory), allocatable :: patchFactories(:)

     end subroutine computeRhsForward

  end interface

  interface

     subroutine computeRhsAdjoint(simulationFlags, solverOptions, grid, state, patchFactories)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use Patch_factory, only : t_PatchFactory
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid) :: grid
       class(t_State) :: state
       type(t_PatchFactory), allocatable :: patchFactories(:)

     end subroutine computeRhsAdjoint

  end interface

  interface

     subroutine computeRhsLinearized(simulationFlags, solverOptions, grid, state, patchFactories)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use Patch_factory, only : t_PatchFactory
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid) :: grid
       class(t_State) :: state
       type(t_PatchFactory), allocatable :: patchFactories(:)

     end subroutine computeRhsLinearized

  end interface

  interface

     subroutine addInterfaceAdjointPenalty(simulationFlags, solverOptions,                   &
          grid, state, patchFactories)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use Patch_factory, only : t_PatchFactory
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state
       type(t_PatchFactory), allocatable :: patchFactories(:)

     end subroutine addInterfaceAdjointPenalty

  end interface

end module RhsHelper
