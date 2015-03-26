#include "config.h"

module RK4Integrator_type

  implicit none
  private

  type, private :: t_RK4Temporary

     SCALAR_TYPE, allocatable :: buffer1(:,:), buffer2(:,:)

  end type t_RK4Temporary

  type, public :: t_RK4Integrator

     type(t_RK4Temporary), allocatable :: temp_(:)

  end type t_RK4Integrator

end module RK4Integrator_type

module RK4Integrator_mod

  implicit none
  public

  interface

     subroutine setupRK4Integrator(this, region)

       use Region_type, only : t_Region
       use RK4Integrator_type, only : t_RK4Integrator

       type(t_RK4Integrator) :: this
       type(t_Region), intent(in) :: region

     end subroutine setupRK4Integrator

  end interface

  interface

     subroutine cleanupRK4Integrator(this)

       use RK4Integrator_type, only : t_RK4Integrator

       type(t_RK4Integrator) :: this

     end subroutine cleanupRK4Integrator

  end interface

  interface

     subroutine substepForward(this, region, time, timestep, stage)

       use Region_type, only : t_Region
       use RK4Integrator_type, only : t_RK4Integrator

       type(t_RK4Integrator) :: this
       type(t_Region) :: region
       real(SCALAR_KIND), intent(inout) :: time
       integer, intent(in) :: timestep, stage

     end subroutine substepForward

  end interface

  interface

     subroutine substepAdjoint(this, region, time, timestep, stage)

       use Region_type, only : t_Region
       use RK4Integrator_type, only : t_RK4Integrator

       type(t_RK4Integrator) :: this
       type(t_Region) :: region
       real(SCALAR_KIND), intent(inout) :: time
       integer, intent(in) :: timestep, stage

     end subroutine substepAdjoint

  end interface

end module RK4Integrator_mod
