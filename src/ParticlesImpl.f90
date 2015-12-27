#include "config.h"

subroutine setupParticles(this, comm)

  ! <<< Derived types >>>
  use Particles_mod, only : t_Particles

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  ! <<< Enumerations >>>
  use Particles_enum

  implicit none

  ! <<< Arguments >>>
  class(t_Particles) :: this
  integer, intent(in) :: comm

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)

end subroutine setupParticles

subroutine cleanupParticles(this)

  ! <<< Derived types >>>
  use Particles_mod, only : t_Particles

  implicit none

  ! <<< Arguments >>>
  class(t_Particles) :: this

end subroutine cleanupParticles

subroutine loadParticleData(this, filename)

  ! <<< Derived types >>>
  use Particles_mod, only : t_Particles

  implicit none

  ! <<< Arguments >>>
  class(t_Particles) :: this
  character(len = STRING_LENGTH), intent(in) :: filename

end subroutine loadParticleData
