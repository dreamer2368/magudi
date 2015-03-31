#include "config.h"

module TimeIntegrator_factory

  implicit none
  public

  interface

     subroutine createTimeIntegrator(timeIntegrator, identifier)

       use TimeIntegrator_mod, only : t_TimeIntegrator

       class(t_TimeIntegrator), pointer, intent(out) :: timeIntegrator
       character(len = *), intent(in) :: identifier

     end subroutine createTimeIntegrator

  end interface
  
end module TimeIntegrator_factory
