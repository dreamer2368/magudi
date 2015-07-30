#include "config.h"

module RK4Integrator_mod

  use TimeIntegrator_mod, only : t_TimeIntegrator

  implicit none

  type, private :: t_RK4IntegratorInternal
     SCALAR_TYPE, allocatable :: buffer1(:,:), buffer2(:,:)
  end type t_RK4IntegratorInternal

  type, extends(t_TimeIntegrator), public :: t_RK4Integrator

     type(t_RK4IntegratorInternal), allocatable :: data_(:)

   contains

     procedure, pass :: setup => setupRK4Integrator
     procedure, pass :: cleanup => cleanupRK4Integrator
     procedure, pass :: substepForward => substepForwardRK4
     procedure, pass :: substepAdjoint => substepAdjointRK4

  end type t_RK4Integrator

  interface

     subroutine setupRK4Integrator(this, region)

       use Region_mod, only : t_Region

       import :: t_RK4Integrator

       class(t_RK4Integrator) :: this
       class(t_Region), intent(in) :: region

     end subroutine setupRK4Integrator

  end interface

  interface

     subroutine cleanupRK4Integrator(this)

       import :: t_RK4Integrator

       class(t_RK4Integrator) :: this

     end subroutine cleanupRK4Integrator

  end interface

  interface

     subroutine substepForwardRK4(this, region, time, timeStepSize, timestep, stage)

       use Region_mod, only : t_Region

       import :: t_RK4Integrator

       class(t_RK4Integrator) :: this
       class(t_Region) :: region
       real(SCALAR_KIND), intent(inout) :: time
       real(SCALAR_KIND), intent(in) :: timeStepSize
       integer, intent(in) :: timestep, stage

     end subroutine substepForwardRK4

  end interface

  interface

     subroutine substepAdjointRK4(this, region, time, timeStepSize, timestep, stage)

       use Region_mod, only : t_Region

       import :: t_RK4Integrator

       class(t_RK4Integrator) :: this
       class(t_Region) :: region
       real(SCALAR_KIND), intent(inout) :: time
       real(SCALAR_KIND), intent(in) :: timeStepSize
       integer, intent(in) :: timestep, stage

     end subroutine substepAdjointRK4

  end interface

end module RK4Integrator_mod
