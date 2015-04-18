#include "config.h"

module Solver_mod

  use Functional_factory, only : t_FunctionalFactory
  use ResidualManager_mod, only : t_ResidualManager
  use ReverseMigrator_mod, only : t_ReverseMigrator
  use TimeIntegrator_factory, only : t_TimeIntegratorFactory

  implicit none
  private

  type, public :: t_Solver

     type(t_ResidualManager) :: residualManager
     type(t_FunctionalFactory) :: functionalFactory
     type(t_TimeIntegratorFactory) :: timeIntegratorFactory

     integer :: saveInterval, reportInterval
     character(len = STRING_LENGTH) :: outputPrefix

   contains

     procedure, pass :: setup => setupSolver
     procedure, pass :: cleanup => cleanupSolver
     procedure, pass :: runForward
     procedure, pass :: runAdjoint

  end type t_Solver

  interface

     subroutine setupSolver(this, region, restartFilename, outputPrefix)

       use Region_mod, only : t_Region

       import :: t_Solver

       class(t_Solver) :: this
       class(t_Region) :: region
       character(len = *), intent(in), optional :: restartFilename, outputPrefix

     end subroutine setupSolver

  end interface

  interface

     subroutine cleanupSolver(this)

       import :: t_Solver

       class(t_Solver) :: this

     end subroutine cleanupSolver

  end interface

  interface

     function runForward(this, region, time, timestep, nTimesteps) result(costFunctional)

       use Region_mod, only : t_Region

       import :: t_Solver

       class(t_Solver) :: this
       class(t_Region) :: region
       real(SCALAR_KIND), intent(inout) :: time
       integer, intent(inout) :: timestep
       integer, intent(in) :: nTimesteps

       SCALAR_TYPE :: costFunctional

     end function runForward

  end interface

  interface

     function runAdjoint(this, region, time, timestep, nTimesteps) result(costSensitivity)

       use Region_mod, only : t_Region

       import :: t_Solver

       class(t_Solver) :: this
       class(t_Region) :: region
       real(SCALAR_KIND), intent(inout) :: time
       integer, intent(inout) :: timestep
       integer, intent(in) :: nTimesteps

       SCALAR_TYPE :: costSensitivity

     end function runAdjoint

  end interface

end module Solver_mod
