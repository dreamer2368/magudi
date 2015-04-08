#include "config.h"

module SimulationFlags_mod

  implicit none
  private

  type, public :: t_SimulationFlags

     logical :: viscosityOn           = .false., &
                predictionOnly        = .true.,  &
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
                useContinuousAdjoint  = .false.

   contains

     procedure, pass :: initialize => initializeSimulationFlags

  end type t_SimulationFlags

  interface

     subroutine initializeSimulationFlags(this)

       import :: t_SimulationFlags

       class(t_SimulationFlags) :: this

     end subroutine initializeSimulationFlags

  end interface

end module SimulationFlags_mod
