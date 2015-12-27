#include "config.h"

module Particles_enum

  implicit none
  public

  integer, parameter, public ::                                                              &
       STOKES           = 0,                                                                 &
       TENNETI          = 1,                                                                 &
       BEETSTRA         = 2,                                                                 &
       BELLAN           = 3

end module Particles_enum

module Particles_mod

  implicit none

  type, public :: t_Particles

     integer(KIND=SCALAR_KIND) :: id
     integer, dimension(3) :: gridIndex, stop
     real(SCALAR_KIND) :: diameter, temperature
     SCALAR_TYPE, dimension(3) :: position, velocity, angularVelocity, normalCollision,      &
          tangentialCollision

  end type t_Particles

  interface

     subroutine setupParticles(this, comm)

       import :: t_Particles

       class(t_Particles) :: this
       integer, intent(in) :: comm

     end subroutine setupParticles

  end interface

  interface

     subroutine cleanupParticles(this)

       import :: t_Particles

       class(t_Particles) :: this

     end subroutine cleanupParticles

  end interface

  interface

     subroutine loadParticleData(this, filename)

       import :: t_Particles

       class(t_Particles) :: this
       character(len = *), intent(in) :: filename

     end subroutine loadParticleData

  end interface

end module Particles_mod
