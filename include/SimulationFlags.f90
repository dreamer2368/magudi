#include "config.h"

module SimulationFlags_type

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
                isBaselineAvailable   = .false.

  end type t_SimulationFlags

end module SimulationFlags_type

module SimulationFlags_mod

  implicit none
  public

  interface

     subroutine initializeSimulationFlags(this)

       use SimulationFlags_type

       type(t_SimulationFlags), intent(out) :: this

     end subroutine initializeSimulationFlags

  end interface

end module SimulationFlags_mod
