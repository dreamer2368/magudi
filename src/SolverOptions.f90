#include "config.h"

module SolverOptions_enum

  implicit none
  public

  integer, parameter ::                                                                      &
       FORWARD   = +1,                                                                       &
       ADJOINT   = -1,                                                                       &
       UNDEFINED =  0

end module SolverOptions_enum

module SolverOptions_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  implicit none
  private


  type, public :: t_SolverOptions

     real(SCALAR_KIND) :: reynoldsNumberInverse, prandtlNumberInverse, ratioOfSpecificHeats, &
          powerLawExponent, bulkViscosityRatio, dissipationAmount, densityRange(2),          &
          temperatureRange(2), cfl, timeStepSize
     integer :: nSpecies, nUnknowns
     character(len = STRING_LENGTH) :: discretizationType, timeintegratorType,               &
          costFunctionalType, controllerType, checkpointingScheme

   contains

     procedure, pass :: initialize => initializeSolverOptions

  end type t_SolverOptions

contains

  subroutine initializeSolverOptions(this, nDimensions, simulationFlags, comm)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use SimulationFlags_mod, only : t_SimulationFlags

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

       if (this%reynoldsNumberInverse <= 0.0_wp .or.                                         &
            this%prandtlNumberInverse <= 0.0_wp) then
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

    if (.not. simulationFlags%predictionOnly) then
       this%controllerType = getOption("controller_type", "THERMAL_ACTUATOR")
       this%costFunctionalType = getOption("cost_functional_type", "ACOUSTIC_NOISE")
       this%checkpointingScheme = getOption("checkpointing_scheme", "UNIFORM_CHECKPOINTING")
    end if

  end subroutine initializeSolverOptions

end module SolverOptions_mod
