#include "config.h"

subroutine initializeSimulationFlags(this)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_SimulationFlags) :: this

  ! <<< Local variables >>>
  integer :: comm_
  character(len = STRING_LENGTH) :: message

  this%viscosityOn           = getOption("include_viscous_terms", .false.)
  this%enableController      = getOption("enable_controller", .false.)
  this%enableFunctional      = getOption("enable_functional", .false.)
  this%enableAdjoint         = getOption("enable_adjoint_solver", .false.)
  this%adjointForcingSwitch  = getOption("adjoint_forcing_switch", .true.)
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
  this%compositeDissipation  = getOption("composite_dissipation", .true.)
  this%enableBodyForce       = getOption("enable_body_force", .false.)
  this%checkConservation     = getOption("check_conservation", .false.)
  this%enableIBM             = getOption("enable_immersed_boundary", .false.)

  this%computeTimeAverage = .false.
  if (.not. this%useConstantCfl)                                                             &
       this%computeTimeAverage = getOption("compute_time_average", .false.)

  if (this%enableAdjoint) then
     this%enableController = .true.
     this%enableFunctional = .true.
  end if

  if (this%enableAdjoint .and. this%enableIBM) then
    comm_ = MPI_COMM_WORLD
    write(message, '(A)') "Adjoint mode for immersed boundary method is not implemented!"
    call gracefulExit(comm_, message)
  end if

  this%IsInitialized = .true.

end subroutine initializeSimulationFlags

subroutine assignSimulationFlags(this, simulationFlags)

  use MPI

  ! <<< Derived types >>>
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_SimulationFlags), intent(out) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags

  ! <<< Local variables >>>
  integer :: comm_
  character(len = STRING_LENGTH) :: message

  comm_ = MPI_COMM_WORLD

  if (.not. simulationFlags%IsInitialized) then
     write(message, '(A)') "The assigned simulationFlags is not initialized!"
     call gracefulExit(comm_, message)
  end if

  this%viscosityOn           = simulationFlags%viscosityOn
  this%enableController      = simulationFlags%enableController
  this%enableFunctional      = simulationFlags%enableFunctional
  this%enableAdjoint         = simulationFlags%enableAdjoint
  this%repeatFirstDerivative = simulationFlags%repeatFirstDerivative
  this%useTargetState        = simulationFlags%useTargetState
  this%dissipationOn         = simulationFlags%dissipationOn
  this%isDomainCurvilinear   = simulationFlags%isDomainCurvilinear
  this%manualDomainDecomp    = simulationFlags%manualDomainDecomp
  this%enableSolutionLimits  = simulationFlags%enableSolutionLimits
  this%useConstantCfl        = simulationFlags%useConstantCfl
  this%filterOn              = simulationFlags%filterOn
  this%steadyStateSimulation = simulationFlags%steadyStateSimulation
  this%isBaselineAvailable   = simulationFlags%isBaselineAvailable
  this%useContinuousAdjoint  = simulationFlags%useContinuousAdjoint
  this%compositeDissipation  = simulationFlags%compositeDissipation

  this%computeTimeAverage    = simulationFlags%computeTimeAverage

  this%IsInitialized         = simulationFlags%IsInitialized

end subroutine assignSimulationFlags
