#include "config.h"

module RK4Integrator_type

  implicit none
  private

  integer, parameter :: wp = SCALAR_KIND

  type, private :: t_RK4Temporary

     SCALAR_TYPE, allocatable :: buffer1(:,:), buffer2(:,:)

  end type t_RK4Temporary

  type, public :: t_RK4Integrator

     integer :: nStages = 4
     real(wp) :: norm(4) =                                                                   &
          (/ 1.0_wp / 6.0_wp, 1.0_wp / 3.0_wp, 1.0_wp / 3.0_wp, 1.0_wp / 6.0_wp /)
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

     subroutine substepForward(this, region, time, timeStepSize, timestep, stage)

       use Region_type, only : t_Region
       use RK4Integrator_type, only : t_RK4Integrator

       type(t_RK4Integrator) :: this
       type(t_Region) :: region
       real(SCALAR_KIND), intent(inout) :: time
       real(SCALAR_KIND), intent(in) :: timeStepSize
       integer, intent(in) :: timestep, stage

     end subroutine substepForward

  end interface

  interface

     subroutine substepAdjoint(this, region, time, timeStepSize, timestep, stage)

       use Region_type, only : t_Region
       use RK4Integrator_type, only : t_RK4Integrator

       type(t_RK4Integrator) :: this
       type(t_Region) :: region
       real(SCALAR_KIND), intent(inout) :: time
       real(SCALAR_KIND), intent(in) :: timeStepSize
       integer, intent(in) :: timestep, stage

     end subroutine substepAdjoint

  end interface

end module RK4Integrator_mod
