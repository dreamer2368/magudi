#include "config.h"

subroutine setupFuelSource(this, fuelIndex, location, amplitude, frequency, radius, phase)

  ! <<< Derived types >>>
  use FuelSource_mod, only : t_FuelSource

  implicit none

  ! <<< Arguments >>>
  class(t_FuelSource) :: this
  integer, intent(in) :: fuelIndex
  real(SCALAR_KIND), intent(in) :: location(:), amplitude, frequency, radius
  real(SCALAR_KIND), intent(in), optional :: phase

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)

  assert(fuelIndex > 0)
  assert(size(location) >= 1 .and. size(location) <= 3)
  assert(radius > 0.0_wp)

  this%fuelIndex = fuelIndex

  this%location = 0.0_wp
  this%location(1:size(location)) = location

  this%amplitude = amplitude

  this%angularFrequency = 2.0_wp * pi * frequency
  this%gaussianFactor = 9.0_wp / (2.0_wp * radius ** 2)

  this%phase = 0.0_wp
  if (present(phase)) this%phase = phase

end subroutine setupFuelSource

subroutine addFuelSource(this, time, coordinates, iblank, rightHandSide)

  ! <<< Derived types >>>
  use FuelSource_mod, only : t_FuelSource

  implicit none

  ! <<< Arguments >>>
  class(t_FuelSource) :: this
  real(SCALAR_KIND), intent(in) :: time
  SCALAR_TYPE, intent(in) :: coordinates(:,:)
  integer, intent(in) :: iblank(:)
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions
  real(wp) :: a, r

  nDimensions = size(coordinates, 2)
  assert_key(nDimensions, (1, 2, 3))
  assert(size(rightHandSide,2) >= nDimensions + 2 + this%fuelIndex)

  a = this%amplitude * abs( cos(this%angularFrequency * time + this%phase) )
  do i = 1, size(rightHandSide, 1)
     if (iblank(i) == 0) cycle
     r = real(sum((coordinates(i,:) - this%location(1:nDimensions)) ** 2), wp)
     rightHandSide(i,nDimensions+2+this%fuelIndex) =                                         &
          rightHandSide(i,nDimensions+2+this%fuelIndex) + a * exp(-this%gaussianFactor * r)
  end do

end subroutine addFuelSource
