#include "config.h"

subroutine setupFuelSource(this, fuelIndex, location, amplitude, radius)

  ! <<< Derived types >>>
  use FuelSource_mod, only : t_FuelSource

  implicit none

  ! <<< Arguments >>>
  class(t_FuelSource) :: this
  integer, intent(in) :: fuelIndex
  real(SCALAR_KIND), intent(in) :: location(:), amplitude, radius

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  assert(fuelIndex > 0)
  assert(size(location) >= 1 .and. size(location) <= 3)
  assert(radius > 0.0_wp)

  this%fuelIndex = fuelIndex

  this%location = 0.0_wp
  this%location(1:size(location)) = location

  this%amplitude = amplitude

  this%radius = radius

end subroutine setupFuelSource

subroutine addFuelSource(this, coordinates, iblank, rightHandSide)

  ! <<< Derived types >>>
  use FuelSource_mod, only : t_FuelSource

  implicit none

  ! <<< Arguments >>>
  class(t_FuelSource) :: this
  SCALAR_TYPE, intent(in) :: coordinates(:,:)
  integer, intent(in) :: iblank(:)
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  real(wp) :: r, amount, gaussianFactor

  nDimensions = size(coordinates, 2)
  assert_key(nDimensions, (1, 2, 3))
  assert(size(rightHandSide,2) >= nDimensions + 2 + this%fuelIndex)

  amount = 0.5_wp * this%amplitude / sqrt(2.0_wp * pi)

  gaussianFactor = 0.5_wp / this%radius**2

  do i = 1, size(rightHandSide, 1)

     if (iblank(i) == 0) cycle

     r = real(sum((coordinates(i,:) - this%location(1:nDimensions)) ** 2), wp)

     rightHandSide(i,nDimensions+2+this%fuelIndex) =                                         &
          rightHandSide(i,nDimensions+2+this%fuelIndex) + amount * exp(- gaussianFactor * r)

  end do

end subroutine addFuelSource
