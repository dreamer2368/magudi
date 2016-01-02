#include "config.h"

subroutine initializeSolverOptions(this, nDimensions, simulationFlags, comm)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Controller_mod, only : t_Controller
  use Functional_mod, only : t_Functional
  use SolverOptions_mod, only : t_SolverOptions
  use Controller_factory, only : t_ControllerFactory
  use Functional_factory, only : t_FunctionalFactory
  use TimeIntegrator_mod, only : t_TimeIntegrator
  use SimulationFlags_mod, only : t_SimulationFlags
  use TimeIntegrator_factory, only : t_TimeIntegratorFactory

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit
  use ThermoChemistry, only : getMolecularWeight

  ! <<< Enumerations >>>
  use SolverOptions_enum

  implicit none

  ! <<< Arguments >>>
  class(t_SolverOptions), intent(out) :: this
  integer, intent(in) :: nDimensions
  type(t_SimulationFlags), intent(in) :: simulationFlags
  integer, intent(in), optional :: comm

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, comm_
  character(len = STRING_LENGTH) :: message, referenceSpecies, val
  real(wp) :: froudeNumberMagnitude, referenceMolecularWeight
  type(t_TimeIntegratorFactory) :: timeIntegratorFactory
  class(t_TimeIntegrator), pointer :: dummyTimeIntegrator => null()
  type(t_ControllerFactory) :: controllerFactory
  class(t_Controller), pointer :: dummyController => null()
  type(t_FunctionalFactory) :: functionalFactory
  class(t_Functional), pointer :: dummyFunctional => null()

  assert_key(nDimensions, (1, 2, 3))

  comm_ = MPI_COMM_WORLD
  if (present(comm)) comm_ = comm

  this%ratioOfSpecificHeats = getOption("ratio_of_specific_heats", 1.4_wp)

  this%nSpecies = getOption("number_of_species", 0)
  if (this%nSpecies < 0) then
     write(message, '(A)') "Number of species must be non-negative!"
     call gracefulExit(comm_, message)
  end if

  this%nUnknowns = nDimensions + 2 + this%nSpecies

  froudeNumberMagnitude = max(0.0_wp, getOption("Froude_number", 0.0_wp))
  if (froudeNumberMagnitude > 0.0_wp) then
     allocate(this%froudeNumberInverse(nDimensions))
     do i = 1, nDimensions
        write(val, '(A,I1.1)') "gravity_norm_dir", i
        this%froudeNumberInverse(i) = getOption(trim(val), 0.0_wp)
     end do
     if (abs(sqrt(sum(this%froudeNumberInverse ** 2)) - 1.0_wp) > epsilon(0.0_wp)) then
        write(message, '(A)') "Gravity norm must sum to 1!"
        call gracefulExit(comm_, message)
     end if
     do i = 1, nDimensions
        if (abs(this%froudeNumberInverse(i)) <= epsilon(0.0_wp)) then
           this%froudeNumberInverse(i) = 0.0_wp
        else
           this%froudeNumberInverse(i) = 1.0_wp / (froudeNumberMagnitude * this%froudeNumberInverse(i))
        end if
     end do
  end if

  if (simulationFlags%viscosityOn) then

     this%reynoldsNumberInverse = max(0.0_wp, getOption("Reynolds_number", 0.0_wp))
     this%prandtlNumberInverse = max(0.0_wp, getOption("Prandtl_number", 0.7_wp))
     allocate(this%schmidtNumberInverse(this%nSpecies+1))
     do i = 1, this%nSpecies
        write(message, "(A,I1.1)") "Schmidt_number_", i
        this%schmidtNumberInverse(i) = max(0.0_wp, getOption(trim(message), 0.7_wp))
     end do
     write(message, "(A)") "Schmidt_number_inert"
     this%schmidtNumberInverse(this%nSpecies+1) =                                            &
          max(0.0_wp, getOption(trim(message), 0.7_wp))

     if (this%reynoldsNumberInverse <= 0.0_wp .or. this%prandtlNumberInverse <= 0.0_wp) then
        this%powerLawExponent = 0.0_wp
        this%bulkViscosityRatio = 0.0_wp
        if (this%nSpecies > 0) this%schmidtNumberInverse = 0.0_wp
     else
        this%reynoldsNumberInverse = 1.0_wp / this%reynoldsNumberInverse
        this%prandtlNumberInverse = 1.0_wp / this%prandtlNumberInverse
        this%powerLawExponent = getOption("viscosity_power_law_exponent", 0.666_wp)
        this%bulkViscosityRatio = getOption("bulk_viscosity_ratio", 0.6_wp)
        if (this%nSpecies > 0) this%schmidtNumberInverse = 1.0_wp / this%schmidtNumberInverse
     end if

  end if

  if (simulationFlags%enableSolutionLimits) then
     call getRequiredOption("minimum_density", this%densityRange(1), comm)
     call getRequiredOption("maximum_density", this%densityRange(2), comm)
     call getRequiredOption("minimum_temperature", this%temperatureRange(1), comm)
     call getRequiredOption("maximum_temperature", this%temperatureRange(2), comm)
     if (this%nSpecies > 0) then
        call getRequiredOption("minimum_mass_fraction", this%massFractionRange(1), comm)
        call getRequiredOption("maximum_mass_fraction", this%massFractionRange(2), comm)
     end if
  end if

  if (simulationFlags%dissipationOn) then
     call getRequiredOption("dissipation_amount", this%dissipationAmount, comm)
  end if

  if (simulationFlags%useConstantCfl) then
     this%cfl = getOption("cfl", 0.5_wp)
  else
     call getRequiredOption("time_step_size", this%timeStepSize, comm)
  end if

  this%discretizationType = getOption("defaults/discretization_scheme", "SBP 4-8")

  this%timeIntegratorType = getOption("time_integration_scheme", "RK4")
  call timeIntegratorFactory%connect(dummyTimeIntegrator, trim(this%timeIntegratorType))
  if (.not. associated(dummyTimeIntegrator)) then
     write(message, '(A)') "Invalid time integration scheme '",                              &
          trim(this%timeIntegratorType), "'!"
     call gracefulExit(comm_, message)
  end if

  if (.not. simulationFlags%predictionOnly) then

     this%controllerType = getOption("controller_type", "THERMAL_ACTUATOR")
     call controllerFactory%connect(dummyController, trim(this%controllerType))
     if (.not. associated(dummyController)) then
        write(message, '(3A)') "Invalid controller type '",                                  &
             trim(this%controllerType), "'!"
        call gracefulExit(comm_, message)
     end if

     this%costFunctionalType = getOption("cost_functional_type", "SOUND")
     call functionalFactory%connect(dummyFunctional, trim(this%costFunctionalType))
     if (.not. associated(dummyFunctional)) then
        write(message, '(3A)') "Invalid cost functional type '",                             &
             trim(this%costFunctionalType), "'!"
        call gracefulExit(comm_, message)
     end if

     this%checkpointingScheme = getOption("checkpointing_scheme", "uniform checkpointing")

  end if

  val = getOption("equation_of_state", "IDEAL_GAS")
  select case (trim(val))
  case ('IDEAL_GAS')
     this%equationOfState = IDEAL_GAS

  case ('IDEAL_GAS_MIXTURE')
     if (this%nSpecies > 0) then

        this%equationOfState = IDEAL_GAS_MIXTURE

        allocate(this%speciesName(this%nSpecies+1))
        allocate(this%molecularWeightInverse(this%nSpecies+1))

        ! Get the molecular weight of the reference species.
        referenceSpecies = getOption("reference_species", "AIR")
        call getMolecularWeight(comm_, trim(referenceSpecies), referenceMolecularWeight)

        ! Get the molecular weights of each species.
        do i = 1, this%nSpecies
           write(message, "(A,I1.1)") "species_", i
           call getRequiredOption(trim(message), this%speciesName(i), comm)
           call getMolecularWeight(comm_, trim(this%speciesName(i)),                         &
                this%molecularWeightInverse(i))
           this%molecularWeightInverse(i) = referenceMolecularWeight /                       &
                this%molecularWeightInverse(i)
        end do

        ! Get the molecular weight of the inert species.
        this%speciesName(this%nSpecies+1) = getOption("inert_species", "N2")
        call getMolecularWeight(comm_, trim(this%speciesName(this%nSpecies+1)),              &
             this%molecularWeightInverse(this%nSpecies+1))
        this%molecularWeightInverse(this%nSpecies+1) = referenceMolecularWeight /            &
             this%molecularWeightInverse(this%nSpecies+1)
     else

        write(message, '(A)') "IDEAL_GAS_MIXTURE requires nSpecies > 0!"
        call gracefulExit(comm_, message)

     end if

  case default
     write(message, '(A)') "Invalid equation of state '", trim(val), "'!"
     call gracefulExit(comm_, message)
  end select

end subroutine initializeSolverOptions
