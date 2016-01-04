#include "config.h"

module Particle_mod

  implicit none

  type, public :: t_Particle

     integer(KIND=SCALAR_KIND) :: id
     integer, dimension(3) :: gridIndex, stop
     real(SCALAR_KIND) :: diameter, temperature
     SCALAR_TYPE, dimension(3) :: position, velocity, angularVelocity, normalCollision,      &
          tangentialCollision

  end type t_Particle

end module Particle_mod
