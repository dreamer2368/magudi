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

  type, public :: t_particleOptions

     integer :: nParticlesLocal, nParticlesGlobal, dragModel
     integer, dimension(:,:,:), allocatable   :: nPartInCell
     integer, dimension(:,:,:,:), allocatable :: partInCell
     SCALAR_TYPE, dimension(:,:,:), allocatable :: auxilaryGrid
     real(SCALAR_KIND) :: densityRatio, ReynoldsNumber, stokesNumber,                        &
          coefficientOfRestitution, coefficientOfFriction, collisionTime
     logical :: twoWayCoupling, collisionsOn

   contains

     procedure, pass :: initialize => initializeParticleOptions

  end type t_particleOptions

  interface

     subroutine initializeParticleOptions(this, simulationFlags, comm)

       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ParticleOptions

       class(t_ParticleOptions), intent(out) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       integer, intent(in), optional :: comm

     end subroutine initializeParticleOptions

  end interface

  interface

     subroutine setupParticles(this, simulationFlags, comm)

       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Particles

       class(t_Particles) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       integer, intent(in) :: comm

     end subroutine setupParticles

  end interface

  interface

     subroutine cleanupParticles(this, simulationFlags)

       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Particles

       class(t_Particles) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags

     end subroutine cleanupParticles

  end interface

  interface

     subroutine loadParticleData(this, simulationFlags, filename)

       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Particles

       class(t_Particles) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = *), intent(in) :: filename

     end subroutine loadParticleData

  end interface

  interface

     subroutine particleSource(this, simulationFlags)

       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Particles

       class(t_Particles) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags

     end subroutine particleSource

  end interface

end module Particles_mod
