#include "config.h"

module WavywallHelperImpl

  implicit none
  public

contains

subroutine updateWallCoordinates(this,grid)
  ! <<< Derived types >>>
  use WallActuator_mod, only : t_WallActuator
  Use Grid_mod, only :t_Grid
  Use GridImpl, only :computeUnitCubeCoordinates
  implicit none

  ! <<< Arguments >>>
  class(t_Grid),intent(in) :: grid
  class(t_WallActuator) :: this

  !do nothing for now
end subroutine

subroutine computePhysicalCoordinates(this,grid)
  ! <<< Derived types >>>
  use WallActuator_mod, only : t_WallActuator
  Use Grid_mod, only :t_Grid
  Use GridImpl, only :computeUnitCubeCoordinates
  implicit none

  ! <<< Arguments >>>
  class(t_Grid),intent(in) :: grid
  class(t_WallActuator) :: this
  SCALAR_TYPE,allocatable::unitCoordinates(:,:)
 
  allocate(unitCoordinates(size(grid%coordinates(:,1)),size(grid%coordinates(1,:))))
  call computeUnitCubeCoordinates(grid,unitCoordinates)
  SAFE_DEALLOCATE(unitCoordinates)

end subroutine

subroutine compute_dMijdp(this,grid,dMijdP)
  ! <<< Derived types >>>
  use WallActuator_mod, only : t_WallActuator
  Use Grid_mod, only :t_Grid
  implicit none

  ! <<< Arguments >>>
  class(t_Grid),intent(in) :: grid
  class(t_WallActuator) :: this
  SCALAR_TYPE,allocatable::dMijdP(:,:,:,:)
  integer, parameter :: wp = SCALAR_KIND

     !assert the sizes are correct
     dMijdp(:,:,:,:)=0._wp
     call computePhysicalCoordinates(this,grid)
end subroutine

subroutine compute_Jacobian(this,grid,J)
  ! <<< Derived types >>>
  use WallActuator_mod, only : t_WallActuator
  Use Grid_mod, only :t_Grid
  implicit none

  ! <<< Arguments >>>
  class(t_Grid),intent(in) :: grid
  class(t_WallActuator) :: this
  SCALAR_TYPE,allocatable,dimension(:)::J
  integer, parameter :: wp = SCALAR_KIND

     !assert the sizes are correct 
     J(:)=0._wp
end subroutine

subroutine compute_dJdp(this,grid,dJdp)
  ! <<< Derived types >>>
  use WallActuator_mod, only : t_WallActuator
  Use Grid_mod, only :t_Grid
  implicit none

  ! <<< Arguments >>>
  class(t_Grid),intent(in) :: grid
  class(t_WallActuator) :: this
  SCALAR_TYPE,allocatable,dimension(:,:)::dJdp
  integer, parameter :: wp = SCALAR_KIND

     !assert the sizes are correct 
     dJdp(:,:)=0._wp
end subroutine

end module WavywallHelperImpl

