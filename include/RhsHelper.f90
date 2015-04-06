#include "config.h"

module RhsHelper

  implicit none
  public

  interface

     subroutine computeRhsForward(time, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       real(SCALAR_KIND), intent(in) :: time
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid) :: grid
       class(t_State) :: state

     end subroutine computeRhsForward

  end interface

  interface

     subroutine computeRhsAdjoint(time, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       real(SCALAR_KIND), intent(in) :: time
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid) :: grid
       class(t_State) :: state

     end subroutine computeRhsAdjoint

  end interface

end module RhsHelper
