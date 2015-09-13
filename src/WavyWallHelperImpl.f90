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

subroutine verifyActuationAmount(this,actuationAmount,direction)
use WallActuator_mod, only : t_WallActuator
implicit none

class(t_WallActuator) :: this
integer, parameter :: wp = SCALAR_KIND
SCALAR_TYPE::sumOfSquares,actuationAmount
SCALAR_TYPE,dimension(:)::direction
integer::i
SCALAR_TYPE::dummyAmplitude

sumOfSquares=0._wp
do i=2,size(this%p),2
dummyAmplitude=this%po(i-1)+actuationAmount*direction(i-1)
sumOfSquares=sumOfSquares+dummyAmplitude**2
end do

if (sumOfSquares .gt. this%MAX_WAVY_WALL_SUM_SQUARES) then
actuationAmount=actuationAmount/sumOfSquares*this%MAX_WAVY_WALL_SUM_SQUARES
end if


end subroutine

subroutine checkAmplitudes(this)
use WallActuator_mod, only : t_WallActuator
implicit none
SCALAR_TYPE::sumOfSquares
integer, parameter :: wp = SCALAR_KIND
integer::i
class(t_WallActuator) :: this

!sumOfSquares=0._wp
!do i=2,size(this%p),2
!sumOfSquares=sumOfSquares+this%p(i-1)*this%p(i-1)
!end do
!
!if (sumOfSquares.gt.this%MAX_WAVY_WALL_SUM_SQUARES)then
!do i=2,size(this%p),2
!this%p(i-1)=this%p(i-1)/sqrt(sumOfSquares)*sqrt(this%MAX_WAVY_WALL_SUM_SQUARES)
!end do
!end if


end subroutine

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

  !call checkAmplitudes(this)
 
  !with an updated list of p remake the physical coordinates 
     do k = 1, grid%localSize(3)
        do j = 1, grid%localSize(2)
           do i = 1, grid%localSize(1)
               index=i + grid%localSize(1) * (j - 1 +&
                   grid%localSize(2) * (k - 1))
               grid%coordinates(index,1)=f(unitCoordinates(index,1))
               grid%coordinates(index,2)=g(unitCoordinates(index,1),unitCoordinates(index,2),p=this%p,beta=this%beta)
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
               unitCoordinates(index,2),p=this%p,beta=this%beta,dgdpVec=dgdp_vec)
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
               unitCoordinates(index,2),p=this%p,beta=this%beta,dgdpVec=dgdp_vec)
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
integer, parameter :: wp = SCALAR_KIND
real(SCALAR_KIND)::xmin,xmax
xmin=-25._wp; xmax=100._wp
xmin=0._wp; xmax=1._wp
f=xmin+(xmax-xmin)*xi1
end function

real(SCALAR_KIND) function h(xi3)
! <<< Arguments >>>
real(SCALAR_KIND)::xi3
h=xi3
end function

real(SCALAR_KIND) function shapeMollifier(x)
real(SCALAR_KIND)::x
integer, parameter :: wp = SCALAR_KIND
shapeMollifier=tanh(60._wp*(x - 25._wp/125._wp))-tanh(60._wp*(x-75._wp/125._wp))
shapeMollifier=tanh(60._wp*(x - 0.2_wp))-tanh(60._wp*(x-0.8_wp))
end function

real(SCALAR_KIND) function g(xi1,xi2,xi3,p,beta)
! <<< Arguments >>>
real(SCALAR_KIND)::xi1,xi2
real(SCALAR_KIND),optional::xi3
SCALAR_TYPE,allocatable::p(:),beta(:)

integer::i,j
integer, parameter :: wp = SCALAR_KIND
SCALAR_TYPE::amplitude,gstar,chi
SCALAR_TYPE, parameter :: pi = 4.0_wp * atan(1.0_wp)
SCALAR_TYPE::ymin,ymax
SCALAR_TYPE::ystar,y
SCALAR_TYPE::h

ymin=-7._wp; ymax=53._wp
ymin=0._wp; ymax=1._wp
y=ymin+(ymax-ymin)*xi2

assert(mod(size(p),2)==0)

h=0._wp
j=1
do i=2,size(p),2
h=h+p(i-1)*shapeMollifier(xi1)*cos(beta(j)*xi1+p(i))
j=j+1
end do

ystar=ymax*xi2+(1._wp-xi2)*(ymin+h)
chi=tanh(100._wp *xi2)
g=ystar-chi*(ystar-y)

end function

subroutine dgdp(xi1,xi2,xi3,p,beta,dgdpVec) 
! <<< Arguments >>>
real(SCALAR_KIND)::xi1,xi2
real(SCALAR_KIND),optional::xi3
SCALAR_TYPE,allocatable::p(:),beta(:)
SCALAR_TYPE,allocatable::dgdpVec(:),dgstardpVec(:)

integer::i,j
integer, parameter :: wp = SCALAR_KIND
SCALAR_TYPE::amplitude,gstar,width,chi
SCALAR_TYPE, parameter :: pi = 4.0_wp * atan(1.0_wp)
SCALAR_TYPE::ymin,ymax
SCALAR_TYPE::ystar,y
SCALAR_TYPE::h

!assert(size(p) == 1)
assert(mod(size(p),2)==0)
assert(size(dgdpVec) == size(p))

allocate(dgstardpVec(size(dgdpVec)))

dgdpVec=0._wp

j=1
do i=2,size(p),2
dgstardpVec(i-1)=(1._wp-xi2)*shapeMollifier(xi1)*cos(beta(j)*xi1+p(i))
dgstardpVec(i)=(-(1._wp-xi2)*p(i-1)*shapeMollifier(xi1)*sin(beta(j)*xi1+p(i)))
j=j+1
end do

chi=tanh(100._wp *xi2)
dgdpVec=dgstardpVec-chi*(dgstardpVec)

end subroutine

end module WavywallHelperImpl
