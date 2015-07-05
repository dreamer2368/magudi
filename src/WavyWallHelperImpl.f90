#include "config.h"

!This module is specific to the following relationships for physical
!coordinate directions x1,x2,x3

!x1=f(xi1)
!x2=g(xi1,xi2,xi3(optional),p)
!x3=h(xi3)

!it also computes normalized metrics, jacobian, and their changes with respect
!to parameter list p

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
  SCALAR_TYPE::xi1,xi2,xi3
  
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
               grid%coordinates(index,1)=f(unitCoordinates(index,1))
               grid%coordinates(index,2)=g(unitCoordinates(index,1),unitCoordinates(index,2),p=this%p)
               !grid%coordinates(index,3)=h(unitCoordinates(index,3))
           end do
        end do
     end do
  
  SAFE_DEALLOCATE(unitCoordinates)

end subroutine


subroutine compute_dMijdp(this,grid,dMijdp)
  ! <<< Derived types >>>
  use WallActuator_mod, only : t_WallActuator
  Use Grid_mod, only :t_Grid
  Use GridImpl, only :computeUnitCubeCoordinates
  implicit none

  ! <<< Arguments >>>
  class(t_Grid),intent(inout) :: grid 
  class(t_WallActuator) :: this
  SCALAR_TYPE,allocatable::dMijdp(:,:,:,:)

 ! <<< Locals >>>
  SCALAR_TYPE,allocatable::dF1(:,:),dG1(:,:),dF2(:,:),dG2(:,:),dGdp_(:,:)
  SCALAR_TYPE,allocatable::dG(:,:)
  SCALAR_TYPE,allocatable::dgdp_vec(:)
  integer, parameter :: wp = SCALAR_KIND
  SCALAR_TYPE,allocatable::unitCoordinates(:,:)
  integer:: k,j,i,index
  SCALAR_TYPE::xi1,xi2,xi3

!assert sizes and dimensions
assert_key(grid%nDimensions,(2))
assert(size(dMijdp,1)==grid%ngridpoints)
assert(size(dMijdp,2)==grid%nDimensions)
assert(size(dMijdp,3)==grid%nDimensions)
assert(size(dMijdp,4)==size%this(p))

SAFE_DEALLOCATE(unitCoordinates)
SAFE_DEALLOCATE(dGdp_)
allocate(unitCoordinates(size(grid%coordinates(:,1)),size(grid%coordinates(1,:))))
allocate(dGdp_(grid%nGridPoints,size(this%p)))
allocate(dgdp_vec(size(this%p)))
call computeUnitCubeCoordinates(grid,unitCoordinates)

!with an updated list of p remake the physical coordinates 
k=1
   do j = 1, grid%localSize(2)
      do i = 1, grid%localSize(1)
          index=i + grid%localSize(1) * (j - 1 +&
              grid%localSize(2) * (k - 1))
          call dgdp(unitCoordinates(index,1),&
               unitCoordinates(index,2),p=this%p,dgdp_vec=dgdp_vec)
          dGdp_(index,:)=dgdp_vec(:)
     end do
   end do

SAFE_DEALLOCATE(dgdp_vec)
SAFE_DEALLOCATE(unitCoordinates)

SAFE_DEALLOCATE(dF1)
SAFE_DEALLOCATE(dF2)
SAFE_DEALLOCATE(dG1)
SAFE_DEALLOCATE(dG2)
SAFE_DEALLOCATE(dG)
allocate(dF1(grid%nGridPoints,1))
allocate(dF2(grid%nGridPoints,1))
allocate(dG1(grid%nGridPoints,1))
allocate(dG2(grid%nGridPoints,1))
allocate(dG(grid%nGridPoints,1))

dF1(:,1) = grid%coordinates(:,1)
call grid%firstDerivative(1)%apply(dF1, grid%localSize)
dF2(:,1) = grid%coordinates(:,1)
call grid%firstDerivative(2)%apply(dF2, grid%localSize)
dG1(:,1) = grid%coordinates(:,2)
call grid%firstDerivative(1)%apply(dG1, grid%localSize)
dG2(:,1) = grid%coordinates(:,2)
call grid%firstDerivative(2)%apply(dG2, grid%localSize)

dMijdp=0.

do i=1,size(this%p)
     dG(:,1) = dGdp_(:,i)
     call grid%firstDerivative(2)%apply(dG, grid%localSize)

     dMijdp(:,1,1,i)=-1._wp/(dF1(:,1)**2)/(dG2(:,1)**2)
     dMijdp(:,1,1,i)=dMijdp(:,1,1,i)*dG(:,1)

     dMijdp(:,2,1,i)=2._wp*dG1(:,1)*dG(:,1)/dF1(:,1)**2/dG2(:,1)**3

     dG(:,1) = dGdp_(:,i)
     call grid%firstDerivative(1)%apply(dG, grid%localSize)
     dMijdp(:,2,1,i)=dMijdp(:,2,1,i)-dG(:,1)/dF1(:,1)**2/dG2(:,1)**2

     dG(:,1) = dGdp_(:,i)
     call grid%firstDerivative(2)%apply(dG, grid%localSize)

     dMijdp(:,2,2,i)=-2._wp*dG(:,1)/dF1(:,1)/dG2(:,1)**3
end do

SAFE_DEALLOCATE(dF1)
SAFE_DEALLOCATE(dF2)
SAFE_DEALLOCATE(dG)
SAFE_DEALLOCATE(dG1)
SAFE_DEALLOCATE(dG2)
SAFE_DEALLOCATE(dGdp_)


end subroutine

subroutine compute_dJacobiandp(this,grid,dJdp)
  ! <<< Derived types >>>
  use WallActuator_mod, only : t_WallActuator
  Use Grid_mod, only :t_Grid
  Use GridImpl, only :computeUnitCubeCoordinates
  implicit none

  ! <<< Arguments >>>
  class(t_Grid),intent(inout) :: grid
  class(t_WallActuator) :: this
  SCALAR_TYPE,allocatable::dJdp(:,:)

 ! <<< Locals >>>
  SCALAR_TYPE,allocatable::dF(:,:),dG(:,:),dGdp_(:,:)
  SCALAR_TYPE,allocatable::dgdp_vec(:)
  integer, parameter :: wp = SCALAR_KIND
  SCALAR_TYPE,allocatable::unitCoordinates(:,:)
  integer:: k,j,i,index
  SCALAR_TYPE::xi1,xi2,xi3

!assert sizes and dimensions
assert_key(grid%nDimensions, (2))
assert(size(dJdp,1)==grid%ngridpoints)
assert(size(dJdp,2)==size(this%p))

SAFE_DEALLOCATE(unitCoordinates)
SAFE_DEALLOCATE(dGdp_)
allocate(unitCoordinates(size(grid%coordinates(:,1)),size(grid%coordinates(1,:))))
allocate(dGdp_(grid%nGridPoints,size(this%p)))
allocate(dgdp_vec(size(this%p)))
call computeUnitCubeCoordinates(grid,unitCoordinates)

!with an updated list of p remake the physical coordinates 
k=1
   do j = 1, grid%localSize(2)
      do i = 1, grid%localSize(1)
          index=i + grid%localSize(1) * (j - 1 +&
              grid%localSize(2) * (k - 1))
          call dgdp(unitCoordinates(index,1),&
               unitCoordinates(index,2),p=this%p,dgdp_vec=dgdp_vec)
          dGdp_(index,:)=dgdp_vec(:)
     end do
   end do

SAFE_DEALLOCATE(dgdp_vec)
SAFE_DEALLOCATE(unitCoordinates)

SAFE_DEALLOCATE(dF)
allocate(dF(grid%nGridPoints,1))
SAFE_DEALLOCATE(dG)
allocate(dG(grid%nGridPoints,1))

dF(:,1) = grid%coordinates(:,1)
call grid%firstDerivative(1)%apply(dF, grid%localSize)

dG(:,1) = grid%coordinates(:,2)
call grid%firstDerivative(2)%apply(dG, grid%localSize)

do i=1,size(this%p)
dJdp(:,i)=-1._wp/(dF(:,1)*dG(:,1)**2)
dG(:,1) = dGdp_(:,i)
call grid%firstDerivative(2)%apply(dG, grid%localSize)
dJdp(:,i)=dJdp(:,i)*dG(:,1)
end do

SAFE_DEALLOCATE(dF)
SAFE_DEALLOCATE(dG)
SAFE_DEALLOCATE(dGdp_)

end subroutine

real(SCALAR_KIND) function f(xi1)
! <<< Arguments >>>
real(SCALAR_KIND)::xi1
f=xi1 !2.*(xi1-0.5)
end function

real(SCALAR_KIND) function h(xi3)
! <<< Arguments >>>
real(SCALAR_KIND)::xi3
h=xi3
end function

real(SCALAR_KIND) function g(xi1,xi2,xi3,p)
! <<< Arguments >>>
real(SCALAR_KIND)::xi1,xi2
real(SCALAR_KIND),optional::xi3
SCALAR_TYPE,allocatable::p(:)

integer, parameter :: wp = SCALAR_KIND
SCALAR_TYPE::amplitude,shapeMollifier
SCALAR_TYPE, parameter :: pi = 4.0_wp * atan(1.0_wp)

assert(size(p) == 1)

amplitude=0.01_wp
shapeMollifier=tanh(40._wp * (xi1 - 0.2_wp)) - tanh(40._wp * (xi1 - 0.8_wp))
shapeMollifier=shapeMollifier/1.2_wp
g=xi2+(1-xi2)*amplitude*shapeMollifier*cos(2._wp*pi*10._wp*xi1+p(1))

end function

subroutine dgdp(xi1,xi2,xi3,p,dgdp_vec) 
! <<< Arguments >>>
real(SCALAR_KIND)::xi1,xi2
real(SCALAR_KIND),optional::xi3
SCALAR_TYPE,allocatable::p(:)
SCALAR_TYPE,allocatable::dgdp_vec(:)

integer::i
integer, parameter :: wp = SCALAR_KIND
SCALAR_TYPE::amplitude,shapeMollifier
SCALAR_TYPE, parameter :: pi = 4.0_wp * atan(1.0_wp)

assert(size(p) == 1)
assert(size(dgdp_vec) == size(p))

dgdp_vec=0._wp
amplitude=0.01_wp
shapeMollifier=tanh(40._wp * (xi1 - 0.2_wp)) - tanh(40._wp * (xi1 - 0.8_wp))
shapeMollifier=shapeMollifier/1.2_wp
dgdp_vec(1)=-(1-xi2)*amplitude*shapeMollifier*sin(2._wp*pi*10._wp*xi1+p(1))

!generalize here to more p

end subroutine

end module WavywallHelperImpl
