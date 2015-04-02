#include "config.h"

subroutine connectTimeIntegrator(this, timeIntegrator, integrationScheme)

  ! <<< Derived types >>>
  use RK4Integrator_mod, only : t_RK4Integrator
  use TimeIntegrator_mod, only : t_TimeIntegrator
  use TimeIntegrator_factory, only : t_TimeIntegratorFactory
  use JamesonRK3Integrator_mod, only : t_JamesonRK3Integrator

  implicit none

  ! <<< Arguments >>>
  class(t_TimeIntegratorFactory) :: this
  class(t_TimeIntegrator), pointer, intent(out) :: timeIntegrator
  character(len = *), intent(in) :: integrationScheme

  assert_key(trim(integrationScheme), ( \
  'RK4', \
  'JamesonRK3'))

  if (.not. associated(this%timeIntegrator)) then

     select case (trim(integrationScheme))

     case ('RK4')
        allocate(t_RK4Integrator :: this%timeIntegrator)

     case ('JamesonRK3')
        allocate(t_JamesonRK3Integrator :: this%timeIntegrator)

     end select

  end if

  timeIntegrator => this%timeIntegrator

end subroutine connectTimeIntegrator

subroutine cleanupTimeIntegratorFactory(this)

  ! <<< Derived types >>>
  use TimeIntegrator_factory, only : t_TimeIntegratorFactory

  implicit none
  
  ! <<< Arguments >>>
  class(t_TimeIntegratorFactory) :: this

  if (associated(this%timeIntegrator)) then
     call this%timeIntegrator%cleanup()
     deallocate(this%timeIntegrator)
  end if
  nullify(this%timeIntegrator)
  
end subroutine cleanupTimeIntegratorFactory
