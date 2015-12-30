#include "config.h"

subroutine initializeParticleOptions(this, simulationFlags, comm)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Particles_mod, only : t_ParticleOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  ! <<< Enumerations >>>
  use Particles_enum

  implicit none

  ! <<< Arguments >>>
  class(t_ParticleOptions), intent(out) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  integer, intent(in), optional :: comm

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: comm_

  if (.not. simulationFlags%particlesOn) return

  comm_ = MPI_COMM_WORLD
  if (present(comm)) comm_ = comm

  this%coefficientOfRestitution = getOption("coefficient_of_restitution", 0.85_wp)

end subroutine initializeParticleOptions

subroutine setupParticles(this, simulationFlags, comm)

  ! <<< Derived types >>>
  use Particles_mod, only : t_Particles
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit
  use MathHelper, only : pi

  ! <<< Enumerations >>>
  use Particles_enum

  implicit none

  ! <<< Arguments >>>
  class(t_Particles) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  integer, intent(in) :: comm

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  if (.not. simulationFlags%particlesOn) return

end subroutine setupParticles

subroutine cleanupParticles(this, simulationFlags)

  ! <<< Derived types >>>
  use Particles_mod, only : t_Particles
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_Particles) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags

  if (.not. simulationFlags%particlesOn) return

end subroutine cleanupParticles

subroutine loadParticleData(this, simulationFlags, filename)

  ! <<< Derived types >>>
  use Particles_mod, only : t_Particles
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_Particles) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  character(len = STRING_LENGTH), intent(in) :: filename

  if (.not. simulationFlags%particlesOn) return

end subroutine loadParticleData

subroutine particleSource(this, simulationFlags)

  ! <<< Derived types >>>
  use Particles_mod, only : t_Particles
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_Particles) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags

  if (.not. simulationFlags%particlesOn) return

end subroutine particleSource
