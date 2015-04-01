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
          saveInterval, outputPrefix, costFunctional)

       use State_mod, only : t_State
       use Region_type, only : t_Region
       use TimeIntegrator_mod, only : t_TimeIntegrator

       type(t_Region) :: region
       class(t_TimeIntegrator) :: integrator
       real(SCALAR_KIND), intent(inout) :: time
       integer, intent(inout) :: timestep
       integer, intent(in) :: nTimesteps
       integer, intent(in), optional :: saveInterval
       character(len = STRING_LENGTH), intent(in), optional :: outputPrefix
       SCALAR_TYPE, intent(out), optional :: costFunctional

     end subroutine solveForward

  end interface

  interface

     subroutine solveAdjoint(region, integrator, time, timestep,                             &
          nTimesteps, saveInterval, outputPrefix)

       use State_mod, only : t_State
       use Region_type, only : t_Region
       use TimeIntegrator_mod, only : t_TimeIntegrator

       type(t_Region) :: region
       class(t_TimeIntegrator) :: integrator
       real(SCALAR_KIND), intent(inout) :: time
       integer, intent(inout) :: timestep
       integer, intent(in) :: nTimesteps, saveInterval
       character(len = STRING_LENGTH), intent(in), optional :: outputPrefix

     end subroutine solveAdjoint

  end interface

end module Solver
