#include "config.h"

module MappingFunction_mod

  implicit none

  type, abstract, public :: t_MappingFunction

   contains

     procedure(compute), nopass, deferred :: compute

  end type t_MappingFunction

  interface

     subroutine compute(coordinate, gridIndex, direction, minMaxRange, isPeriodic)

       real(SCALAR_KIND), intent(out) :: coordinate(:)

       integer, intent(in), optional :: gridIndex, direction
       real(SCALAR_KIND), intent(in), optional :: minMaxRange(2)
       logical, intent(in), optional :: isPeriodic

     end subroutine compute

  end interface

end module MappingFunction_mod
