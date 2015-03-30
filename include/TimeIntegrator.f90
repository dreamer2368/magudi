#include "config.h"

module TimeIntegrator_mod

  implicit none
  private

  type, abstract, public :: t_TimeIntegrator

     integer :: nStages
     real(SCALAR_KIND), allocatable :: norm(:)

   contains

     procedure, non_overridable, pass :: setup_ => setupTimeIntegrator
     procedure, non_overridable, pass :: cleanup_ => cleanupTimeIntegrator
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

     subroutine substepForward(this, region, time, timeStepSize, timestep, stage)

       use Region_type, only : t_Region

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

       use Region_type, only : t_Region

       import :: t_TimeIntegrator

       class(t_TimeIntegrator) :: this
       class(t_Region) :: region
       real(SCALAR_KIND), intent(inout) :: time
       real(SCALAR_KIND), intent(in) :: timeStepSize
       integer, intent(in) :: timestep, stage

     end subroutine substepAdjoint

  end interface

end module TimeIntegrator_mod
