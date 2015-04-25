#include "config.h"

subroutine initializeSimulationFlags(this)

  ! <<< Derived types >>>
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_SimulationFlags) :: this

  this%viscosityOn           = getOption("include_viscous_terms", .false.)
  this%predictionOnly        = getOption("disable_adjoint_solver", .true.)
  this%repeatFirstDerivative = getOption("repeat_first_derivative", .true.)
  this%useTargetState        = getOption("use_target_state", .true.)
  this%dissipationOn         = getOption("add_dissipation", .false.)
  this%isDomainCurvilinear   = getOption("curvilinear_domain", .true.)
  this%manualDomainDecomp    = getOption("use_manual_domain_decomposition", .false.)
  this%enableSolutionLimits  = getOption("enable_solution_limits", .false.)
  this%useConstantCfl        = getOption("use_constant_CFL_mode", .true.)
  this%filterOn              = getOption("filter_solution", .false.)
  this%steadyStateSimulation = getOption("steady_state_simulation", .false.)
  this%isBaselineAvailable   = getOption("baseline_prediction_available", .false.)
  this%useContinuousAdjoint  = getOption("use_continuous_adjoint", .false.)
  this%compositeDissipation  = getOption("use_composite_dissipation_operator", .true.)

end subroutine initializeSimulationFlags
