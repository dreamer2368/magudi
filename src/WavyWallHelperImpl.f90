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
  class(t_Grid),intent(inout) :: grid !anticpate grid changes
  class(t_WallActuator) :: this

  ! <<< Locals >>>
  SCALAR_TYPE,allocatable::unitCoordinates(:,:)
  integer:: k,j,i,index
  
  !obviously functionality only hardcoded for the two-dimensional case
  assert_key(grid%nDimensions, (2))

  allocate(unitCoordinates(size(grid%coordinates(:,1)),size(grid%coordinates(1,:))))
  call computeUnitCubeCoordinates(grid,unitCoordinates)
 
  !with an updated list of p remake the physical coordinates 
     do k = 1, grid%localSize(3)
        do j = 1, grid%localSize(2)
           do i = 1, grid%localSize(1)
               index=i + grid%localSize(1) * (j - 1 +&
                   grid%localSize(2) * (k - 1))
               grid%coordinates(index,1) = unitCoordinates(index,1)
               grid%coordinates(index,2) = unitCoordinates(index,2)+(1-unitCoordinates(index,2))*h(unitCoordinates(index,1),this%p)
           end do
        end do
     end do
  
  SAFE_DEALLOCATE(unitCoordinates)

end subroutine

real(SCALAR_KIND) function h(xi1,p)
! <<< Arguments >>>
real(SCALAR_KIND)::xi1
SCALAR_TYPE,allocatable::p(:)

! <<< Locals >>>
integer, parameter :: wp = SCALAR_KIND
real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
real(SCALAR_KIND)::shapeMollifier,amplitude

!for now we're considering a single p
assert_key(size(p), (1))

amplitude=0.005_wp
shapeMollifier=tanh(40._wp * (xi1 - 0.2_wp)) - tanh(40._wp * (xi1 - 0.8_wp))
shapeMollifier=shapeMollifier/1.2_wp
h=amplitude*shapeMollifier*cos(2._wp*pi*10._wp*xi1+p(1))

end function

subroutine computePhysicalCoordinates(this,grid)
  ! <<< Derived types >>>
  use WallActuator_mod, only : t_WallActuator
  Use Grid_mod, only :t_Grid
  Use GridImpl, only :computeUnitCubeCoordinates
  implicit none

  ! <<< Arguments >>>
  class(t_Grid),intent(in) :: grid
  class(t_WallActuator) :: this

  ! <<< Locals >>>
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

