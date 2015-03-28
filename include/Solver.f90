#include "config.h"

module Solver

  implicit none
  public

  interface

     subroutine initializeSolver(region, restartFilename)

       use Region_type, only : t_Region

       type(t_Region) :: region
       character(len = *), intent(in), optional :: restartFilename

     end subroutine initializeSolver

  end interface

  interface

     subroutine solveForward(region, integrator, time, timestep, nTimesteps,                 &
          saveInterval, reportInterval, outputPrefix, costFunctional)

       use State_type, only : t_State
       use Region_type, only : t_Region
       use RK4Integrator_type, only : t_RK4Integrator

       type(t_Region) :: region
       type(t_RK4Integrator) :: integrator
       real(SCALAR_KIND), intent(inout) :: time
       integer, intent(inout) :: timestep
       integer, intent(in) :: nTimesteps
       integer, intent(in), optional :: saveInterval, reportInterval
       character(len = STRING_LENGTH), intent(in), optional :: outputPrefix
       SCALAR_TYPE, intent(out), optional :: costFunctional

     end subroutine solveForward

  end interface

  interface

     subroutine solveAdjoint(region, integrator, time, timestep, nTimesteps,                 &
          saveInterval, reportInterval, outputPrefix)

       use State_type, only : t_State
       use Region_type, only : t_Region
       use RK4Integrator_type, only : t_RK4Integrator

       type(t_Region) :: region
       type(t_RK4Integrator) :: integrator
       real(SCALAR_KIND), intent(inout) :: time
       integer, intent(inout) :: timestep
       integer, intent(in) :: nTimesteps, saveInterval
       integer, intent(in), optional :: reportInterval
       character(len = STRING_LENGTH), intent(in), optional :: outputPrefix

     end subroutine solveAdjoint

  end interface

end module Solver
