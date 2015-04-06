#include "config.h"

module Solver

  implicit none
  public

  interface

     subroutine initializeSolver(region, restartFilename)

       use Region_mod, only : t_Region

       class(t_Region) :: region
       character(len = *), intent(in), optional :: restartFilename

     end subroutine initializeSolver

  end interface

  interface

     subroutine solveForward(region, time, timestep, nTimesteps,                             &
          saveInterval, outputPrefix, costFunctional)

       use Region_mod, only : t_Region

       class(t_Region) :: region
       real(SCALAR_KIND), intent(inout) :: time
       integer, intent(inout) :: timestep
       integer, intent(in) :: nTimesteps
       integer, intent(in) :: saveInterval
       character(len = STRING_LENGTH), intent(in), optional :: outputPrefix
       SCALAR_TYPE, intent(out), optional :: costFunctional

     end subroutine solveForward

  end interface

  interface

     subroutine solveAdjoint(region, time, timestep, nTimesteps, saveInterval, outputPrefix)

       use Region_mod, only : t_Region

       class(t_Region) :: region
       real(SCALAR_KIND), intent(inout) :: time
       integer, intent(inout) :: timestep
       integer, intent(in) :: nTimesteps, saveInterval
       character(len = STRING_LENGTH), intent(in), optional :: outputPrefix

     end subroutine solveAdjoint

  end interface

end module Solver
