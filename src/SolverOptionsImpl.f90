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
     call getRequiredOption("Reynolds_number", this%reynoldsNumber)
     this%prandtlNumber = getOption("Prandtl_number", 0.72_wp)
     this%powerLawExponent = getOption("viscosity_power_law_exponent", 0.0_wp)
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

  this%spongeAmount = 0.0_wp
  this%spongeExponent = 2

  this%farFieldInviscidPenaltyAmount = 0.0_wp
  this%farFieldViscousPenaltyAmount = 0.0_wp

end subroutine initializeSolverOptions

subroutine updateSolverOptions(this, simulationFlags, patchData)

  ! <<< Derived types >>>
  use SolverOptions_type
  use PatchDescriptor_type
  use SimulationFlags_type

  ! <<< Internal modules >>>
  use InputHelper, only : getRequiredOption, getOption

  implicit none

  ! <<< Arguments >>>
  type(t_SolverOptions), intent(out) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_PatchDescriptor), intent(in) :: patchData(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  if (any(patchData(:)%patchType == SPONGE)) then
     call getRequiredOption("sponge_amount", this%spongeAmount)
     this%spongeExponent = getOption("sponge_exponent", 2)
  end if

  if (any(patchData(:)%patchType == SAT_FAR_FIELD)) then
     this%farFieldInviscidPenaltyAmount = &
          getOption("far_field_inviscid_penalty_amount", 1.0_wp)
     if (simulationFlags%viscosityOn) then
        this%farFieldViscousPenaltyAmount = &
             getOption("far_field_viscous_penalty_amount", 1.0_wp)
     end if
  end if

end subroutine updateSolverOptions
