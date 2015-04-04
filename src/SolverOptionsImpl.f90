#include "config.h"

subroutine initializeSolverOptions(this, nDimensions, simulationFlags, comm)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Functional_mod, only : t_Functional
  use SolverOptions_mod, only : t_SolverOptions
  use Functional_factory, only : t_FunctionalFactory
  use TimeIntegrator_mod, only : t_TimeIntegrator
  use SimulationFlags_mod, only : t_SimulationFlags
  use TimeIntegrator_factory, only : t_TimeIntegratorFactory

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_SolverOptions), intent(out) :: this
  integer, intent(in) :: nDimensions
  type(t_SimulationFlags), intent(in) :: simulationFlags
  integer, intent(in), optional :: comm

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: comm_
  character(len = STRING_LENGTH) :: message
  type(t_TimeIntegratorFactory) :: timeIntegratorFactory
  class(t_TimeIntegrator), pointer :: dummyTimeIntegrator => null()
  type(t_FunctionalFactory) :: functionalFactory
  class(t_Functional), pointer :: dummyFunctional => null()

  assert_key(nDimensions, (1, 2, 3))

  comm_ = MPI_COMM_WORLD
  if (present(comm)) comm_ = comm

  this%ratioOfSpecificHeats = getOption("ratio_of_specific_heats", 1.4_wp)

  this%nUnknowns = nDimensions + 2 !... extendable to reactive flows, for later.

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

  call getRequiredOption("time_integration_scheme", this%timeIntegratorType)
  call timeIntegratorFactory%connect(dummyTimeIntegrator, trim(this%timeIntegratorType))
  if (.not. associated(dummyTimeIntegrator)) then
     write(message, '(A)') "Invalid time integration scheme '",                              &
          trim(this%timeIntegratorType), "'!"
     call gracefulExit(comm_, message)
  end if

  if (.not. simulationFlags%predictionOnly) then
     call getRequiredOption("cost_functional_type", this%costFunctionalType)
     call functionalFactory%connect(dummyFunctional, trim(this%costFunctionalType))
     if (.not. associated(dummyFunctional)) then
        write(message, '(A)') "Invalid cost functional type '",                              &
             trim(this%costFunctionalType), "'!"
        call gracefulExit(comm_, message)
     end if
  end if

end subroutine initializeSolverOptions
