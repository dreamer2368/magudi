#include "config.h"

module UniformMap_mod

  use MappingFunction_mod, only : t_MappingFunction

  implicit none

  type, extends(t_MappingFunction), public :: t_UniformMap

   contains

     procedure, nopass :: compute => computeUniformMap

  end type t_UniformMap

  interface

     subroutine computeUniformMap(coordinate, gridIndex, direction, minMaxRange, isPeriodic)

       real(SCALAR_KIND), intent(out) :: coordinate(:)

       integer, intent(in), optional :: gridIndex, direction
       real(SCALAR_KIND), intent(in), optional :: minMaxRange(2)
       logical, intent(in), optional :: isPeriodic

     end subroutine computeUniformMap

  end interface

end module UniformMap_mod
