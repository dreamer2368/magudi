#include "config.h"

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
