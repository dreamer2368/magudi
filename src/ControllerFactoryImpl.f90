#include "config.h"

subroutine connectController(this, controllerTarget, controllerType, createNew)

  ! <<< Derived types >>>
  use Controller_mod, only : t_Controller
  use Controller_factory, only : t_ControllerFactory
  use ThermalActuator_mod, only : t_ThermalActuator
  use WallActuator_mod, only : t_WallActuator

  implicit none

  ! <<< Arguments >>>
  class(t_ControllerFactory) :: this
  class(t_Controller), pointer, intent(out) :: controllerTarget
  character(len = *), intent(in), optional :: controllerType
  logical, intent(in), optional :: createNew

  ! <<< Local variables >>>
  logical :: createNew_

  createNew_ = .false.
  if (present(createNew)) createNew_ = createNew

  if (present(controllerType) .and. .not. (associated(this%controller) .and.                 &
       .not. createNew_)) then

     if (associated(this%controller)) deallocate(this%controller)
     nullify(this%controller)

     this%controllerType = controllerType

     select case (trim(controllerType))

     case ('THERMAL_ACTUATOR')
        allocate(t_ThermalActuator :: this%controller)

     case ('PASSIVE_WALL_ACTUATOR')
        allocate(t_WallActuator :: this%controller) 

     case default
        this%controllerType = ""

     end select

  end if

  nullify(controllerTarget)
  if (.not. associated(this%controller)) return
  controllerTarget => this%controller

end subroutine connectController

subroutine cleanupControllerFactory(this)

  ! <<< Derived types >>>
  use Controller_factory, only : t_ControllerFactory

  implicit none

  ! <<< Arguments >>>
  class(t_ControllerFactory) :: this

  if (associated(this%controller)) then
     call this%controller%cleanup()
     deallocate(this%controller)
  end if
  nullify(this%controller)

end subroutine cleanupControllerFactory
