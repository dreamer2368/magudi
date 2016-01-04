#include "config.h"

subroutine setupParticles(this, particle, simulationFlags, comm)

  ! <<< Derived types >>>
  use ParticleSolver_mod, only : t_ParticleSolver
  use Particle_mod, only : t_Particle
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  ! <<< Enumerations >>>
  use ParticleSolver_enum

  implicit none

  ! <<< Arguments >>>
  class(t_ParticleSolver) :: this
  type(t_Particle), intent(in) :: particle(:)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  integer, intent(in) :: comm

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  if (.not. simulationFlags%particlesOn) return

  this%coefficientOfRestitution = getOption("coefficient_of_restitution", 0.85_wp)

end subroutine setupParticles

subroutine cleanupParticles(this, particle, simulationFlags)

  ! <<< Derived types >>>
  use ParticleSolver_mod, only : t_ParticleSolver
  use Particle_mod, only : t_Particle
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_ParticleSolver) :: this
  type(t_Particle), intent(in) :: particle(:)
  type(t_SimulationFlags), intent(in) :: simulationFlags

  if (.not. simulationFlags%particlesOn) return

end subroutine cleanupParticles

subroutine updateParticles(this, particle, simulationFlags)

  ! <<< Derived types >>>
  use ParticleSolver_mod, only : t_ParticleSolver
  use Particle_mod, only : t_Particle
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_ParticleSolver) :: this
  type(t_Particle), intent(in) :: particle(:)
  type(t_SimulationFlags), intent(in) :: simulationFlags

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  if (.not. simulationFlags%particlesOn) return

end subroutine updateParticles

subroutine computeParticleRhs(this, particle, simulationFlags)

  ! <<< Derived types >>>
  use ParticleSolver_mod, only : t_ParticleSolver
  use Particle_mod, only : t_Particle
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use ParticleSolver_enum

  implicit none

  ! <<< Arguments >>>
  class(t_ParticleSolver) :: this
  type(t_Particle), intent(in) :: particle
  type(t_SimulationFlags), intent(in) :: simulationFlags

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(WP) :: drag, tau, Rep, b1, b2, F, I, addedMass

  if (.not. simulationFlags%particlesOn) return
  
!!$  ! Interpolate the gas phase info at the droplet location
!!$  call particleIinterpolate
!!$
!!$  ! Particle mass
!!$  md = rhod*pi*myp%d**3/6.0_WP
!!$  mdi = 1.0_WP/md
!!$  
!!$  ! Particle Reynolds number (include gas volume fraction)
!!$  Rep = rhog*eg*sqrt((myp%u-ug)**2+(myp%v-vg)**2+(myp%w-wg)**2)*myp%d/mug
!!$  Rep = Rep+epsilon(1.0_WP)
!!$
!!$  ! Particle response time
!!$  tau = rhod*myp%d**2/(18.0_WP*mug)
!!$  ! Calculate drag coefficient
!!$  select case (this%dragModel)
!!$  case (STOKES)
!!$     ! Stokes drag
!!$     F = 1.0_WP
!!$  case (SCHILLER_NAUMANN)
!!$     ! Clift, Grace & Wever 1978)
!!$     F = 1.0_WP+0.15_WP*Rep**(0.687_WP)
!!$  case (TENNETI)
!!$     ! Tenneti & Subramaniam (2011)
!!$     b1 = 5.81_WP*ep/eg**3+0.48_WP*ep**(1.0_WP/3.0_WP)/eg**4
!!$     b2 = ep**3*Rep*(0.95_WP+0.61_WP*ep**3/eg**2)
!!$     F = (1.0_WP+0.15_WP*Rep**(0.687_WP))/eg**3+b1+b2
!!$     F = eg*F ! To remove mean pressure gradient forcing
!!$  case (BEETSTRA)
!!$     ! Beetstra (2006)
!!$     b1 = 10.0_WP*ep/eg**2 + eg**2 * (1.0_WP+1.5_WP*sqrt(ep))
!!$     b2 = 0.413_WP/24.0_WP*Rep/eg**2*(1.0_WP/eg+3.0_WP*eg*ep+8.4_WP*Rep**(-0.343_WP))/(1.0_WP+10.0_WP**(3.0_WP*ep)*Rep**(2.0_WP*eg-2.5_WP))
!!$     F = b1 + b2
!!$  case default
!!$     stop 'Unknown drag model'
!!$  end select
!!$  
!!$  ! Drag force over particle mass
!!$  drag = F/taud
!!$
!!$  ! Particle moment of inertia per unit mass
!!$  I = 0.1_WP*myp%d**2
!!$  
!!$  ! Gas phase source term
!!$  accx = drag*eg*(ug-myp%u)+stress(1)/rhod
!!$  accy = drag*eg*(vg-myp%v)+stress(2)/rhod
!!$  accz = drag*eg*(wg-myp%w)+stress(3)/rhod
!!$  
!!$  ! Return rhs
!!$  dxdt = myp%u
!!$  dydt = myp%v
!!$  dzdt = myp%w
!!$  dudt = accx+gravity(1)+myp%Acol(1)
!!$  dvdt = accy+gravity(2)+myp%Acol(2)
!!$  dwdt = accz+gravity(3)+myp%Acol(3)
!!$  dwxdt = myp%Tcol(1)/I
!!$  dwydt = myp%Tcol(2)/I
!!$  dwzdt = myp%Tcol(3)/I
!!$
!!$  ! Account for added mass
!!$  !madd = 0.5_WP*rhog/rhod
!!$  !accx = accx - madd * (dudt-myp%Acol(1))
!!$  !accy = accy - madd * (dvdt-myp%Acol(2))
!!$  !accz = accz - madd * (dwdt-myp%Acol(3))
!!$  !dudt = dudt - madd * (dudt-myp%Acol(1))
!!$  !dvdt = dvdt - madd * (dvdt-myp%Acol(2))
!!$  !dwdt = dwdt - madd * (dwdt-myp%Acol(3))

end subroutine computeParticleRhs

subroutine addParticleSource(this, particle, simulationFlags)

  ! <<< Derived types >>>
  use ParticleSolver_mod, only : t_ParticleSolver
  use Particle_mod, only : t_Particle
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_ParticleSolver) :: this
  type(t_Particle), intent(in) :: particle(:)
  type(t_SimulationFlags), intent(in) :: simulationFlags

  if (.not. simulationFlags%particlesOn) return

end subroutine addParticleSource
