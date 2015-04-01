#include "config.h"

subroutine createTimeIntegrator(timeIntegrator, identifier)

  ! <<< Derived types >>>
  use RK4Integrator_mod, only : t_RK4Integrator
  use TimeIntegrator_mod, only : t_TimeIntegrator
  use JamesonRK3Integrator_mod, only : t_JamesonRK3Integrator

  implicit none

  ! <<< Arguments >>>
  class(t_TimeIntegrator), pointer, intent(out) :: timeIntegrator
  character(len = *), intent(in) :: identifier

  assert_key(trim(identifier), ( \
  'RK4', \
  'JamesonRK3'))

  select case (trim(identifier))

  case ('RK4')
     allocate(t_RK4Integrator :: timeIntegrator)

  case ('JamesonRK3')
     allocate(t_JamesonRK3Integrator :: timeIntegrator)

  end select

end subroutine createTimeIntegrator

subroutine setupTimeIntegrator(this)

  ! <<< Derived types >>>
  use TimeIntegrator_mod, only : t_TimeIntegrator

  implicit none

  ! <<< Arguments >>>
  class(t_TimeIntegrator) :: this

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  call this%cleanupBase()
  
  assert(this%nStages > 0)

  allocate(this%norm(this%nStages), source = 0.0_wp)
  this%norm(this%nStages) = 1.0_wp

end subroutine setupTimeIntegrator

subroutine cleanupTimeIntegrator(this)

  ! <<< Derived types >>>
  use TimeIntegrator_mod, only : t_TimeIntegrator

  implicit none

  ! <<< Arguments >>>
  class(t_TimeIntegrator) :: this

  SAFE_DEALLOCATE(this%norm)

end subroutine cleanupTimeIntegrator
