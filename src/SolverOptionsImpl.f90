#include "config.h"

subroutine initializeSolverOptions(this, simulationFlags)

  ! <<< Derived types >>>
  use SolverOptions_type
  use SimulationFlags_type

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption

  implicit none

  ! <<< Arguments >>>
  type(t_SolverOptions), intent(out) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  this%ratioOfSpecificHeats = getOption("ratio_of_specific_heats", 1.4_wp)

  if (simulationFlags%viscosityOn) then
     call getRequiredOption("Reynolds_number", this%reynoldsNumberInverse)
     this%reynoldsNumberInverse = 1.0_wp / this%reynoldsNumberInverse
     this%prandtlNumberInverse = getOption("Prandtl_number", 0.72_wp)
     this%prandtlNumberInverse = 1.0_wp / this%prandtlNumberInverse
     this%powerLawExponent = getOption("viscosity_power_law_exponent", 0.666_wp)
     this%bulkViscosityRatio = getOption("bulk_viscosity_ratio", 0.6_wp)
  end if

  if (simulationFlags%enableSolutionLimits) then
     call getRequiredOption("minimum_density", this%densityRange(1))
     call getRequiredOption("maximum_density", this%densityRange(2))
     call getRequiredOption("minimum_temperature", this%temperatureRange(1))
     call getRequiredOption("maximum_temperature", this%temperatureRange(2))
  end if

  if (simulationFlags%dissipationOn) then
     call getRequiredOption("dissipation_amount", this%dissipationAmount)
  end if

  if (simulationFlags%useConstantCfl) then
     this%cfl = getOption("cfl", 0.5_wp)
  else
     call getRequiredOption("time_step_size", this%timeStepSize)
  end if

  if (simulationFlags%steadyStateSimulation) then
     call getRequiredOption("density_convergence_tolerance",  this%convergenceTolerance(1))
     call getRequiredOption("momentum_convergence_tolerance", this%convergenceTolerance(2))
     call getRequiredOption("energy_convergence_tolerance",   this%convergenceTolerance(3))
  end if

end subroutine initializeSolverOptions
