#include "config.h"

subroutine connectTimeIntegrator(this, timeIntegratorTarget, timeIntegratorType, createNew)

  ! <<< Derived types >>>
  use RK4Integrator_mod, only : t_RK4Integrator
  use TimeIntegrator_mod, only : t_TimeIntegrator
  use TimeIntegrator_factory, only : t_TimeIntegratorFactory
  use JamesonRK3Integrator_mod, only : t_JamesonRK3Integrator

  implicit none

  ! <<< Arguments >>>
  class(t_TimeIntegratorFactory) :: this
  class(t_TimeIntegrator), pointer, intent(out) :: timeIntegratorTarget
  character(len = *), intent(in), optional :: timeIntegratorType
  logical, intent(in), optional :: createNew

  ! <<< Local variables >>>
  logical :: createNew_

  createNew_ = .false.
  if (present(createNew)) createNew_ = createNew

  if (present(timeIntegratorType) .and. .not. (associated(this%timeIntegrator) .and.         &
       .not. createNew_)) then

     if (associated(this%timeIntegrator)) deallocate(this%timeIntegrator)
     nullify(this%timeIntegrator)

     this%timeIntegratorType = timeIntegratorType

     select case (trim(timeIntegratorType))

     case ('RK4')
        allocate(t_RK4Integrator :: this%timeIntegrator)

     case ('JamesonRK3')
        allocate(t_JamesonRK3Integrator :: this%timeIntegrator)

     case default
        this%timeIntegratorType = ""

     end select

  end if

  nullify(timeIntegratorTarget)
  if (.not. associated(this%timeIntegrator)) return
  timeIntegratorTarget => this%timeIntegrator

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
