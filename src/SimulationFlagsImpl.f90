#include "config.h"

subroutine initializeSimulationFlags(this)

  ! <<< Derived types >>>
  use SimulationFlags_type

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  type(t_SimulationFlags), intent(out) :: this

  this%viscosityOn           = getOption("include_viscous_terms", .false.)
  this%predictionOnly        = getOption("disable_adjoint_solver", .true.)
  this%repeatFirstDerivative = getOption("repeat_first_derivative", .true.)
  this%useTargetState        = getOption("use_target_state", .true.)
  this%baselineOnly          = getOption("disable_control", .false.)
  this%dissipationOn         = getOption("add_dissipation", .false.)
  this%isDomainCurvilinear   = getOption("curvilinear_domain", .true.)
  this%manualDomainDecomp    = getOption("use_manual_domain_decomposition", .false.)
  this%enableSolutionLimits  = getOption("enable_solution_limits", .false.)
  this%useConstantCfl        = getOption("use_constant_CFL_mode", .true.)
  this%filterOn              = getOption("filter_solution", .false.)

end subroutine initializeSimulationFlags
