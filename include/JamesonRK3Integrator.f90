#include "config.h"

module JamesonRK3Integrator_mod

  use TimeIntegrator_mod, only : t_TimeIntegrator

  implicit none
  private

  type, private :: t_JamesonRK3IntegratorInternal
     SCALAR_TYPE, allocatable :: buffer1(:,:), buffer2(:,:)
  end type t_JamesonRK3IntegratorInternal

  type, extends(t_TimeIntegrator), public :: t_JamesonRK3Integrator

     type(t_JamesonRK3IntegratorInternal), allocatable :: data_(:)

   contains

     procedure, pass :: setup => setupJamesonRK3Integrator
     procedure, pass :: cleanup => cleanupJamesonRK3Integrator
     procedure, pass :: substepForward => substepForwardJamesonRK3
     procedure, pass :: substepAdjoint => substepAdjointJamesonRK3

  end type t_JamesonRK3Integrator

  interface

     subroutine setupJamesonRK3Integrator(this, region)

       use Region_mod, only : t_Region

       import :: t_JamesonRK3Integrator

       class(t_JamesonRK3Integrator) :: this
       class(t_Region), intent(in) :: region

     end subroutine setupJamesonRK3Integrator

  end interface

  interface

     subroutine cleanupJamesonRK3Integrator(this)

       import :: t_JamesonRK3Integrator

       class(t_JamesonRK3Integrator) :: this

     end subroutine cleanupJamesonRK3Integrator

  end interface

  interface

     subroutine substepForwardJamesonRK3(this, region, time, timeStepSize, timestep, stage)

       use Region_mod, only : t_Region

       import :: t_JamesonRK3Integrator

       class(t_JamesonRK3Integrator) :: this
       class(t_Region) :: region
       real(SCALAR_KIND), intent(inout) :: time
       real(SCALAR_KIND), intent(in) :: timeStepSize
       integer, intent(in) :: timestep, stage

     end subroutine substepForwardJamesonRK3

  end interface

  interface

     subroutine substepAdjointJamesonRK3(this, region, time, timeStepSize, timestep, stage)

       use Region_mod, only : t_Region

       import :: t_JamesonRK3Integrator

       class(t_JamesonRK3Integrator) :: this
       class(t_Region) :: region
       real(SCALAR_KIND), intent(inout) :: time
       real(SCALAR_KIND), intent(in) :: timeStepSize
       integer, intent(in) :: timestep, stage

     end subroutine substepAdjointJamesonRK3

  end interface

end module JamesonRK3Integrator_mod
