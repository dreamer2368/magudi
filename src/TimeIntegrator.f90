#include "config.h"

module TimeIntegrator_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  implicit none
  private

  type, abstract, public :: t_TimeIntegrator

     integer :: nStages
     real(SCALAR_KIND), allocatable :: norm(:)

   contains

     procedure, non_overridable, pass :: setupBase
     procedure, non_overridable, pass :: cleanupBase

     procedure(setup), pass, deferred :: setup
     procedure(cleanup), pass, deferred :: cleanup
     procedure(subStepForward), pass, deferred :: subStepForward
     procedure(subStepAdjoint), pass, deferred :: subStepAdjoint

  end type t_TimeIntegrator

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

contains

  subroutine setupBase(this)

    implicit none

    ! <<< Arguments >>>
    class(t_TimeIntegrator) :: this

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND

    call this%cleanupBase()

    assert(this%nStages > 0)

    allocate(this%norm(this%nStages), source = 0.0_wp)
    this%norm(this%nStages) = 1.0_wp

  end subroutine setupBase

  subroutine cleanupBase(this)

    implicit none

    ! <<< Arguments >>>
    class(t_TimeIntegrator) :: this

    SAFE_DEALLOCATE(this%norm)

  end subroutine cleanupBase

end module TimeIntegrator_mod
