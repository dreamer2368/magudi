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
  SCALAR_TYPE,allocatable::dMijdp(:,:,:)

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
assert(size(dMijdp,2)==grid%nDimensions*grid%nDimensions)
assert(size(dMijdp,3)==size%this(p))

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
               unitCoordinates(index,2),p=this%p,dgdpVec=dgdp_vec)
          dGdp_(index,:)=dgdp_vec(:)
     end do
   end do

SAFE_DEALLOCATE(dgdp_vec)
SAFE_DEALLOCATE(unitCoordinates)

allocate(dG(grid%nGridPoints,1))

dMijdp=0.

do i=1,size(this%p)

     dG(:,1) = dGdp_(:,i)
     call grid%firstDerivative(2)%apply(dG, grid%localSize)
     dMijdp(:,1,i)=dG(:,1)
     
     dG(:,1) = dGdp_(:,i)
     call grid%firstDerivative(1)%apply(dG, grid%localSize)
     dMijdp(:,3,i)=-dG(:,1)

end do

SAFE_DEALLOCATE(dG)
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
               unitCoordinates(index,2),p=this%p,dgdpVec=dgdp_vec)
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
end do

do i=1,size(this%p)
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

integer::i,j
integer, parameter :: wp = SCALAR_KIND
SCALAR_TYPE::amplitude,shapeMollifier,gstar,width,chi
SCALAR_TYPE, parameter :: pi = 4.0_wp * atan(1.0_wp)

!assert(size(p) == 1)
assert(mod(size(p),2)==0)

width=0.2_wp
shapeMollifier=tanh(40._wp * (xi1 - 0.3_wp)) - tanh(40._wp * (xi1 - 0.7_wp))
shapeMollifier=shapeMollifier/0.8_wp

gstar=xi2

j=1
do i=2,size(p),2
gstar=gstar+(1._wp-xi2)*p(i-1)*shapeMollifier*cos(2._wp*pi*3*real(j)*xi1+p(i))
j=j+1
end do

if (xi2 .le. width) then
chi=tanh(5._wp *xi2/width)
else
chi=1._wp
end if

g=gstar-chi*(gstar-xi2)

end function

subroutine dgdp(xi1,xi2,xi3,p,dgdpVec) 
! <<< Arguments >>>
real(SCALAR_KIND)::xi1,xi2
real(SCALAR_KIND),optional::xi3
SCALAR_TYPE,allocatable::p(:)
SCALAR_TYPE,allocatable::dgdpVec(:),dgstardpVec(:)

integer::i,j
integer, parameter :: wp = SCALAR_KIND
SCALAR_TYPE::amplitude,shapeMollifier,gstar,width,chi
SCALAR_TYPE, parameter :: pi = 4.0_wp * atan(1.0_wp)

!assert(size(p) == 1)
assert(mod(size(p),2)==0)
assert(size(dgdpVec) == size(p))

allocate(dgstardpVec(size(dgdpVec)))

width=0.2_wp
dgdpVec=0._wp
shapeMollifier=tanh(40._wp * (xi1 - 0.3_wp)) - tanh(40._wp * (xi1 - 0.7_wp))
shapeMollifier=shapeMollifier/0.8_wp


j=1
do i=2,size(p),2
dgstardpVec(i-1)=(1._wp-xi2)*shapeMollifier*cos(2._wp*pi*3._wp*real(j)*xi1+p(i))
dgstardpVec(i)=-(1._wp-xi2)*p(i-1)*shapeMollifier*sin(2._wp*pi*3._wp*real(j)*xi1+p(i))
j=j+1
end do

if (xi2 .le. width) then
chi=tanh(5._wp *xi2/width)
else
chi=1._wp
end if

dgdpVec=dgstardpVec-chi*(dgstardpVec)

end subroutine

end module WavywallHelperImpl
