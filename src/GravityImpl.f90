#include "config.h"

subroutine addGravityForward(this, iblank, density, froudeNumberInverse, rightHandSide)

  ! <<< Derived types >>>
  use Gravity_mod, only : t_Gravity

  implicit none

  ! <<< Arguments >>>
  class(t_Gravity) :: this
  SCALAR_TYPE, intent(in) :: density(:), froudeNumberInverse(:)
  integer, intent(in) :: iblank(:)
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions

  nDimensions = size(froudeNumberInverse)
  assert_key(nDimensions, (1, 2, 3))
  assert(size(rightHandSide, 2) >= nDimensions + 2)

  do i = 1, size(rightHandSide, 1)
     if (iblank(i) == 0) cycle

     rightHandSide(i,2:nDimensions+2) = rightHandSide(i,2:nDimensions+1) +                   &
          density(i) * froudeNumberInverse(1:nDimensions)

  end do

end subroutine addGravityForward

subroutine addGravityAdjoint(this, iblank, adjointVariables, froudeNumberInverse,            &
     rightHandSide)

  ! <<< Derived types >>>
  use Gravity_mod, only : t_Gravity

  implicit none

  ! <<< Arguments >>>
  class(t_Gravity) :: this
  SCALAR_TYPE, intent(in) :: adjointVariables(:,:), froudeNumberInverse(:)
  integer, intent(in) :: iblank(:)
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions, nUnknowns, nGridPoints
  SCALAR_TYPE, allocatable :: localSourceJacobian(:,:), temp(:)

  nDimensions = size(froudeNumberInverse)
  assert_key(nDimensions, (1, 2, 3))
  assert(size(rightHandSide, 2) >= nDimensions + 2)

  nUnknowns = size(adjointVariables, 2)
  assert(nUnknowns > 0)
  assert(size(rightHandSide, 2) == nUnknowns)

  nGridPoints = size(adjointVariables,1)
  assert(nGridPoints > 0)

  allocate(localSourceJacobian(nUnknowns, nUnknowns))
  allocate(temp(nUnknowns))

  localSourceJacobian = 0.0_wp 
  localSourceJacobian(2:nDimensions+1,1) = froudeNumberInverse(1:nDimensions)

  do i = 1, nGridPoints

     if (iblank(i) == 0) cycle

     temp = matmul(transpose(localSourceJacobian), adjointVariables(i,:))

     rightHandSide(i,:) = rightHandSide(i,:) - temp

  end do !... i = 1, nGridPoints

  SAFE_DEALLOCATE(localSourceJacobian)
  SAFE_DEALLOCATE(temp)

end subroutine addGravityAdjoint
