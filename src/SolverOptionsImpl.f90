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

  if (simulationFlags%viscosityOn) then

     this%reynoldsNumberInverse = max(0.0_wp, getOption("Reynolds_number", 0.0_wp))
     this%prandtlNumberInverse = max(0.0_wp, getOption("Prandtl_number", 0.72_wp))

     if (this%reynoldsNumberInverse <= 0.0_wp .or. this%prandtlNumberInverse <= 0.0_wp) then
        this%powerLawExponent = 0.0_wp
        this%bulkViscosityRatio = 0.0_wp
     else
        this%reynoldsNumberInverse = 1.0_wp / this%reynoldsNumberInverse
        this%prandtlNumberInverse = 1.0_wp / this%prandtlNumberInverse
        this%powerLawExponent = getOption("viscosity_power_law_exponent", 0.666_wp)
        this%bulkViscosityRatio = getOption("bulk_viscosity_ratio", 0.6_wp)
     end if

  end if

  if (simulationFlags%enableSolutionLimits) then
     call getRequiredOption("minimum_density", this%densityRange(1), comm)
     call getRequiredOption("maximum_density", this%densityRange(2), comm)
     call getRequiredOption("minimum_temperature", this%temperatureRange(1), comm)
     call getRequiredOption("maximum_temperature", this%temperatureRange(2), comm)
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

  if (simulationFlags%enableController) then
     this%controllerType = getOption("controller_type", "THERMAL_ACTUATOR")
     call controllerFactory%connect(dummyController, trim(this%controllerType))
     if (.not. associated(dummyController)) then
        write(message, '(3A)') "Invalid controller type '",                             &
             trim(this%controllerType), "'!"
        call gracefulExit(comm_, message)
     end if

     this%controllerNorm = getOption("controller_norm", "L1")
     this%controllerFactor = getOption("controller_factor", 12.0_wp)
  else
     this%controllerType = ""
     this%controllerNorm = ""
     this%controllerFactor = 0.0_wp
  end if

  if (simulationFlags%enableFunctional) then
     this%costFunctionalType = getOption("cost_functional_type", "SOUND")
     call functionalFactory%connect(dummyFunctional, trim(this%costFunctionalType))
     if (.not. associated(dummyFunctional)) then
        write(message, '(3A)') "Invalid cost functional type '",                             &
             trim(this%costFunctionalType), "'!"
        call gracefulExit(comm_, message)
     end if
  else
    this%costFunctionalType = ""
  end if

  if (simulationFlags%enableAdjoint) then
     this%checkpointingScheme = getOption("checkpointing_scheme", "uniform checkpointing")
  end if

  this%IsInitialized = .true.

end subroutine initializeSolverOptions

subroutine assignSolverOptions(this, solverOptions)

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

  implicit none

  ! <<< Arguments >>>
  class(t_SolverOptions), intent(out) :: this
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: comm_
  character(len = STRING_LENGTH) :: message
  type(t_TimeIntegratorFactory) :: timeIntegratorFactory
  class(t_TimeIntegrator), pointer :: dummyTimeIntegrator => null()
  type(t_ControllerFactory) :: controllerFactory
  class(t_Controller), pointer :: dummyController => null()
  type(t_FunctionalFactory) :: functionalFactory
  class(t_Functional), pointer :: dummyFunctional => null()

  comm_ = MPI_COMM_WORLD

  if (.not. solverOptions%IsInitialized) then
     write(message, '(A)') "The assigned solverOptions is not initialized!"
     call gracefulExit(comm_, message)
  end if

  this%ratioOfSpecificHeats = solverOptions%ratioOfSpecificHeats

  this%nSpecies = solverOptions%nSpecies
  if (this%nSpecies < 0) then
     write(message, '(A)') "Number of species must be non-negative!"
     call gracefulExit(comm_, message)
  end if

  this%nUnknowns = solverOptions%nUnknowns

  this%reynoldsNumberInverse = solverOptions%reynoldsNumberInverse
  this%prandtlNumberInverse = solverOptions%prandtlNumberInverse
  this%powerLawExponent = solverOptions%powerLawExponent
  this%bulkViscosityRatio = solverOptions%bulkViscosityRatio

  this%densityRange = solverOptions%densityRange
  this%temperatureRange = solverOptions%temperatureRange

  this%dissipationAmount = solverOptions%dissipationAmount

  this%cfl = solverOptions%cfl
  this%timeStepSize = solverOptions%timeStepSize

  this%discretizationType = solverOptions%discretizationType

  this%timeIntegratorType = solverOptions%timeIntegratorType
  if (len(trim(this%timeIntegratorType))>0) then
    call timeIntegratorFactory%connect(dummyTimeIntegrator, trim(this%timeIntegratorType))
    if (.not. associated(dummyTimeIntegrator)) then
       write(message, '(A)') "Invalid time integration scheme '",                              &
            trim(this%timeIntegratorType), "'!"
       call gracefulExit(comm_, message)
    end if
  end if

  this%controllerType = solverOptions%controllerType
  if ( len(trim(this%controllerType))>0 ) then
    call controllerFactory%connect(dummyController, trim(this%controllerType))
    if (.not. associated(dummyController)) then
      write(message, '(3A)') "Invalid controller type '",                             &
           trim(this%controllerType), "'!"
      call gracefulExit(comm_, message)
    end if
  end if
  this%controllerNorm = solverOptions%controllerNorm

  this%costFunctionalType = solverOptions%costFunctionalType
  if ( len(trim(this%costFunctionalType))>0 ) then
     call functionalFactory%connect(dummyFunctional, trim(this%costFunctionalType))
     if (.not. associated(dummyFunctional)) then
        write(message, '(3A)') "Invalid cost functional type '",                             &
             trim(this%costFunctionalType), "'!"
        call gracefulExit(comm_, message)
     end if
  end if

  this%checkpointingScheme = solverOptions%checkpointingScheme
  this%IsInitialized = solverOptions%IsInitialized

end subroutine assignSolverOptions
