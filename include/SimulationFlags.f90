#include "config.h"

module SimulationFlags_mod

  implicit none

  type, public :: t_SimulationFlags

     logical :: viscosityOn           = .false., &
                enableController      = .true.,  &
                enableFunctional      = .true.,  &
                enableAdjoint         = .true.,  &
                adjointForcingSwitch  = .true.,  &
                repeatFirstDerivative = .true.,  &
                useTargetState        = .true.,  &
                dissipationOn         = .false., &
                isDomainCurvilinear   = .true.,  &
                manualDomainDecomp    = .false., &
                enableSolutionLimits  = .false., &
                useConstantCfl        = .true.,  &
                filterOn              = .false., &
                steadyStateSimulation = .false., &
                isBaselineAvailable   = .false., &
                useContinuousAdjoint  = .false., &
                compositeDissipation  = .true.,  &
                computeTimeAverage    = .false., &
                enableBodyForce       = .false., &
                checkConservation     = .false., &
                IsInitialized         = .false.

   contains

     procedure, pass :: initialize => initializeSimulationFlags
     procedure, pass :: assignSimulationFlags
     generic :: assignment(=) => assignSimulationFlags

  end type t_SimulationFlags

  interface

     subroutine initializeSimulationFlags(this)

       import :: t_SimulationFlags

       class(t_SimulationFlags) :: this

     end subroutine initializeSimulationFlags

  end interface

  interface

     subroutine assignSimulationFlags(this, simulationFlags)

       import :: t_SimulationFlags

       class(t_SimulationFlags), intent(out) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags

     end subroutine assignSimulationFlags

  end interface

end module SimulationFlags_mod
