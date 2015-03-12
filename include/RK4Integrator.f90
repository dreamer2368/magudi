#include "config.h"

module RK4Integrator_type

  implicit none
  private

  type, private :: t_RK4Temporary

     SCALAR_TYPE, allocatable :: buffer1(:,:), buffer2(:,:)

  end type t_RK4Temporary

  type, public :: t_RK4Integrator

     type(t_RK4Temporary), allocatable :: temp(:)

  end type t_RK4Integrator

end module RK4Integrator_type

module RK4Integrator_mod

  implicit none

  interface

     subroutine setupRK4Integrator(this, region)

       use RK4Integrator_type
       use Region_type

       type(t_RK4Integrator) :: this
       type(t_Region), intent(in) :: region

     end subroutine setupRK4Integrator

  end interface

  interface

     subroutine cleanupRK4Integrator(this)

       use RK4Integrator_type

       type(t_RK4Integrator) :: this

     end subroutine cleanupRK4Integrator

  end interface

  interface

     subroutine stepForward(this, region, time, timestep, verbose)

       use RK4Integrator_type
       use Region_type

       type(t_RK4Integrator) :: this
       type(t_Region) :: region
       real(SCALAR_KIND), intent(inout) :: time
       integer, intent(in) :: timestep
       logical, intent(in), optional :: verbose

     end subroutine stepForward

  end interface

end module RK4Integrator_mod
