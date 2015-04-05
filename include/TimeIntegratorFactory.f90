#include "config.h"

module TimeIntegrator_factory

  use TimeIntegrator_mod, only : t_TimeIntegrator

  implicit none
  private

  type, public :: t_TimeIntegratorFactory

     class(t_TimeIntegrator), pointer :: timeIntegrator => null()

   contains

     procedure, pass :: connect => connectTimeIntegrator
     procedure, pass :: cleanup => cleanupTimeIntegratorFactory

  end type t_TimeIntegratorFactory

  interface

     subroutine connectTimeIntegrator(this, timeIntegratorTarget, timeIntegratorType, createNew)

       use TimeIntegrator_mod, only : t_TimeIntegrator

       import :: t_TimeIntegratorFactory

       class(t_TimeIntegratorFactory) :: this
       class(t_TimeIntegrator), pointer, intent(out) :: timeIntegratorTarget
       character(len = *), intent(in), optional :: timeIntegratorType
       logical, intent(in), optional :: createNew

     end subroutine connectTimeIntegrator

  end interface

  interface

     subroutine cleanupTimeIntegratorFactory(this)

       import :: t_TimeIntegratorFactory

       class(t_TimeIntegratorFactory) :: this

     end subroutine cleanupTimeIntegratorFactory

  end interface

end module TimeIntegrator_factory
