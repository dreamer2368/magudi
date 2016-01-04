#include "config.h"

module ParticleSolver_enum

  implicit none
  public

  integer, parameter, public ::                                                              &
       STOKES           = 0,                                                                 &
       SCHILLER_NAUMANN = 1,                                                                 &
       TENNETI          = 2,                                                                 &
       BEETSTRA         = 3

end module ParticleSolver_enum

module ParticleSolver_mod

  use Particle_mod, only : t_Particle

  implicit none

  type, public :: t_ParticleSolver

     integer :: nParticlesLocal, nParticlesGlobal, dragModel
     integer, dimension(:,:,:), allocatable   :: nPartInCell
     integer, dimension(:,:,:,:), allocatable :: partInCell
     SCALAR_TYPE, dimension(:,:,:), allocatable :: auxilaryGrid
     real(SCALAR_KIND) :: densityRatio, ReynoldsNumber, stokesNumber,                        &
          coefficientOfRestitution, coefficientOfFriction, collisionTime
     logical :: twoWayCoupling, collisionsOn

   contains

     procedure, pass :: setup => setupParticles
     procedure, pass :: cleanup => cleanupParticles
     procedure, pass :: update => updateParticles
     procedure, pass :: computeRhs => computeParticleRhs
     procedure, pass :: add => addParticleSource

  end type t_ParticleSolver

  interface

     subroutine setupParticles(this, particle, simulationFlags, comm)

       use Particle_mod, only : t_Particle
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ParticleSolver

       class(t_ParticleSolver) :: this
       type(t_Particle), intent(in) :: particle(:)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       integer, intent(in) :: comm

     end subroutine setupParticles

  end interface

  interface

     subroutine cleanupParticles(this, particle, simulationFlags)

       use Particle_mod, only : t_Particle
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ParticleSolver

       class(t_ParticleSolver) :: this
       type(t_Particle), intent(in) :: particle(:)
       type(t_SimulationFlags), intent(in) :: simulationFlags

     end subroutine cleanupParticles

  end interface

  interface

     subroutine updateParticles(this, particle, simulationFlags)

       use Particle_mod, only : t_Particle
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ParticleSolver

       class(t_ParticleSolver) :: this
       type(t_Particle), intent(in) :: particle(:)
       type(t_SimulationFlags), intent(in) :: simulationFlags

     end subroutine updateParticles

  end interface

  interface

     subroutine computeParticleRhs(this, particle, simulationFlags)

       use Particle_mod, only : t_Particle
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ParticleSolver

       class(t_ParticleSolver) :: this
       type(t_Particle), intent(in) :: particle
       type(t_SimulationFlags), intent(in) :: simulationFlags

     end subroutine computeParticleRhs

  end interface

  interface

     subroutine addParticleSource(this, particle, simulationFlags)

       use Particle_mod, only : t_Particle
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ParticleSolver

       class(t_ParticleSolver) :: this
       type(t_Particle), intent(in) :: particle(:)
       type(t_SimulationFlags), intent(in) :: simulationFlags

     end subroutine addParticleSource

  end interface

end module ParticleSolver_mod
