#include "config.h"

module SolverOptionsImpl

  implicit none
  public

contains

  subroutine parseCostFunctionalType(costFunctionalTypeString, costFunctionalType)

    ! <<< Enumerations >>>
    use SolverOptions_enum

    ! <<< Arguments >>>
    character(len = *), intent(in) :: costFunctionalTypeString
    integer, intent(out) :: costFunctionalType

    costFunctionalType = -1

    if (trim(costFunctionalTypeString) == "sound") then
       costFunctionalType = SOUND
    else if (trim(costFunctionalTypeString) == "lift") then
       costFunctionalType = LIFT
    else if (trim(costFunctionalTypeString) == "drag") then
       costFunctionalType = DRAG
    end if

  end subroutine parseCostFunctionalType

end module SolverOptionsImpl

subroutine initializeSolverOptions(this, simulationFlags, comm)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Private members >>>
  use SolverOptionsImpl, only : parseCostFunctionalType

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_SolverOptions), intent(out) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  integer, intent(in), optional :: comm

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: comm_
  character(len = STRING_LENGTH) :: str, message

  comm_ = MPI_COMM_WORLD
  if (present(comm)) comm_ = comm

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

  if (.not. simulationFlags%predictionOnly) then
     call getRequiredOption("cost_functional_type", str)
     call parseCostFunctionalType(str, this%costFunctionalType)
     if (this%costFunctionalType == -1) then
        write(message, '(A)') "Invalid cost functional type '", trim(str), "'!"
        call gracefulExit(comm_, message)
     end if
  end if

end subroutine initializeSolverOptions
