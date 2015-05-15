#include "config.h"

module Controller_factory

  use Controller_mod, only : t_Controller

  implicit none
  private

  type, public :: t_ControllerFactory

     class(t_Controller), pointer :: controller => null()
     character(len = STRING_LENGTH) :: controllerType = ""

   contains

     procedure, pass :: connect => connectController
     procedure, pass :: cleanup => cleanupControllerFactory

  end type t_ControllerFactory

  interface

     subroutine connectController(this, controllerTarget, controllerType, createNew)

       use Controller_mod, only : t_Controller

       import :: t_ControllerFactory

       class(t_ControllerFactory) :: this
       class(t_Controller), pointer, intent(out) :: controllerTarget
       character(len = *), intent(in), optional :: controllerType
       logical, intent(in), optional :: createNew

     end subroutine connectController

  end interface

  interface

     subroutine cleanupControllerFactory(this)

       import :: t_ControllerFactory

       class(t_ControllerFactory) :: this

     end subroutine cleanupControllerFactory

  end interface

end module Controller_factory
