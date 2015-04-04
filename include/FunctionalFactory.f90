#include "config.h"

module Functional_factory

  use Functional_mod, only : t_Functional

  implicit none
  private

  type, public :: t_FunctionalFactory

     class(t_Functional), pointer :: functional => null()
     character(len = STRING_LENGTH) :: functionalType = ""

   contains

     procedure, pass :: connect => connectFunctional
     procedure, pass :: cleanup => cleanupFunctionalFactory

  end type t_FunctionalFactory

  interface

     subroutine connectFunctional(this, functionalTarget, functionalType, createNew)

       use Functional_mod, only : t_Functional

       import :: t_FunctionalFactory

       class(t_FunctionalFactory) :: this
       class(t_Functional), pointer, intent(out) :: functionalTarget
       character(len = *), intent(in), optional :: functionalType
       logical, intent(in), optional :: createNew

     end subroutine connectFunctional

  end interface

  interface

     subroutine cleanupFunctionalFactory(this)

       import :: t_FunctionalFactory

       class(t_FunctionalFactory) :: this

     end subroutine cleanupFunctionalFactory

  end interface

end module Functional_factory
