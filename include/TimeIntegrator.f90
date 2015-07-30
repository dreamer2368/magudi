#include "config.h"

module TimeIntegrator_mod

  implicit none

  type, abstract, public :: t_TimeIntegrator

     integer :: nStages
     real(SCALAR_KIND), allocatable :: norm(:)

   contains

     procedure, non_overridable, pass :: setupBase => setupTimeIntegrator
     procedure, non_overridable, pass :: cleanupBase => cleanupTimeIntegrator

     procedure(setup), pass, deferred :: setup
     procedure(cleanup), pass, deferred :: cleanup
     procedure(subStepForward), pass, deferred :: subStepForward
     procedure(subStepAdjoint), pass, deferred :: subStepAdjoint

  end type t_TimeIntegrator

  interface

     subroutine setupTimeIntegrator(this)

       import :: t_TimeIntegrator

       class(t_TimeIntegrator) :: this

     end subroutine setupTimeIntegrator

  end interface

  interface

     subroutine cleanupTimeIntegrator(this)

       import :: t_TimeIntegrator

       class(t_TimeIntegrator) :: this

     end subroutine cleanupTimeIntegrator

  end interface

  abstract interface

     subroutine setup(this, region)

       use Region_mod, only : t_Region

       import :: t_TimeIntegrator

       class(t_TimeIntegrator) :: this
       class(t_Region), intent(in) :: region

     end subroutine setup

  end interface

  abstract interface

     subroutine cleanup(this)

       import :: t_TimeIntegrator

       class(t_TimeIntegrator) :: this

     end subroutine cleanup

  end interface

  abstract interface

     subroutine substepForward(this, region, time, timeStepSize, timestep, stage)

       use Region_mod, only : t_Region

       import :: t_TimeIntegrator

       class(t_TimeIntegrator) :: this
       class(t_Region) :: region
       real(SCALAR_KIND), intent(inout) :: time
       real(SCALAR_KIND), intent(in) :: timeStepSize
       integer, intent(in) :: timestep, stage

     end subroutine substepForward

  end interface

  abstract interface

     subroutine substepAdjoint(this, region, time, timeStepSize, timestep, stage)

       use Region_mod, only : t_Region

       import :: t_TimeIntegrator

       class(t_TimeIntegrator) :: this
       class(t_Region) :: region
       real(SCALAR_KIND), intent(inout) :: time
       real(SCALAR_KIND), intent(in) :: timeStepSize
       integer, intent(in) :: timestep, stage

     end subroutine substepAdjoint

  end interface

end module TimeIntegrator_mod
