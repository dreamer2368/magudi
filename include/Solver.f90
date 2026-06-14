#include "config.h"

module Solver_mod

  use Controller_factory, only : t_ControllerFactory
  use Functional_factory, only : t_FunctionalFactory
  use ResidualManager_mod, only : t_ResidualManager
  use ReverseMigrator_mod, only : t_ReverseMigrator
  use TimeIntegrator_factory, only : t_TimeIntegratorFactory

  implicit none

  type, public :: t_Solver

     type(t_ResidualManager) :: residualManager
     type(t_ControllerFactory) :: controllerFactory
     type(t_FunctionalFactory) :: functionalFactory
     type(t_TimeIntegratorFactory) :: timeIntegratorFactory

     integer :: nTimesteps, saveInterval, reportInterval, probeInterval
     character(len = STRING_LENGTH) :: outputPrefix

   contains

     procedure, pass :: setup => setupSolver
     procedure, pass :: cleanup => cleanupSolver
     procedure, pass :: runForward
     procedure, pass :: runAdjoint
     procedure, pass :: runLinearized
     procedure, pass :: checkGradientAccuracy

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

     function runForward(this, region, restartFilename, controlTimestepOffset) result(costFunctional)

       use Region_mod, only : t_Region

       import :: t_Solver

       class(t_Solver) :: this
       class(t_Region) :: region

       character(len = *), intent(in), optional :: restartFilename
       integer, intent(in), optional :: controlTimestepOffset

       SCALAR_TYPE :: costFunctional

     end function runForward

  end interface

  interface

     ! Multi-segment callers pass controlTimestepOffset (start timestep of the
     ! current segment, k*Nts) AND controlTotalTimesteps (the full trajectory
     ! length, Nsplit*Nts). runAdjoint derives the FORWARD-hook reference
     ! (= controlTimestepOffset) and the ADJOINT-hook reference
     ! (= controlTotalTimesteps - controlTimestepOffset - this%nTimesteps).
     ! The two arguments must be provided together or both omitted.
     function runAdjoint(this, region, controlTimestepOffset, controlTotalTimesteps,            &
                         deleteGradientFile) result(costSensitivity)

       use Region_mod, only : t_Region

       import :: t_Solver

       class(t_Solver) :: this
       class(t_Region) :: region

       integer, intent(in), optional :: controlTimestepOffset
       integer, intent(in), optional :: controlTotalTimesteps
       logical, intent(in), optional :: deleteGradientFile

       SCALAR_TYPE :: costSensitivity

     end function runAdjoint

  end interface

  interface
     ! right now dummy value does not have any meaning. just for future purpose.
     function runLinearized(this, region) result(dummyValue)

       use Region_mod, only : t_Region

       import :: t_Solver

       class(t_Solver) :: this
       class(t_Region) :: region

       SCALAR_TYPE :: dummyValue

     end function runLinearized

  end interface

  interface

     subroutine checkGradientAccuracy(this, region)

       use Region_mod, only : t_Region

       import :: t_Solver

       class(t_Solver) :: this
       class(t_Region) :: region

     end subroutine checkGradientAccuracy

  end interface

end module Solver_mod
