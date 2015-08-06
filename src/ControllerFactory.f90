#include "config.h"

module Controller_factory

  use Controller_mod, only : t_Controller

  implicit none
  private

  type, public :: t_ControllerFactory

     class(t_Controller), pointer :: controller => null()
     character(len = STRING_LENGTH) :: controllerType = ""

   contains

     procedure, pass :: connect
     procedure, pass :: cleanup

  end type t_ControllerFactory

contains

  subroutine connect(this, controllerTarget, controllerType, createNew)

    ! <<< Derived types >>>
    use ThermalActuator_mod, only : t_ThermalActuator

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

    if (present(controllerType) .and. .not.                                                  &
         (associated(this%controller) .and. .not. createNew_)) then

       if (associated(this%controller)) then
          call this%controller%cleanup()
          deallocate(this%controller)
       end if

       nullify(this%controller)

       this%controllerType = controllerType

       select case (trim(controllerType))

       case ('THERMAL_ACTUATOR')
          allocate(t_ThermalActuator :: this%controller)

       case default
          this%controllerType = ""

       end select

    end if

    nullify(controllerTarget)
    if (.not. associated(this%controller)) return
    controllerTarget => this%controller

  end subroutine connect

  subroutine cleanup(this)

    class(t_ControllerFactory) :: this

    if (associated(this%controller)) then
       call this%controller%cleanup()
       deallocate(this%controller)
    end if
    nullify(this%controller)

  end subroutine cleanup

end module Controller_factory
