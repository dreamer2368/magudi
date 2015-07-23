#include "config.h"

subroutine setupIgnitionSource(this, location, amplitude, radius, timeStart, timeDuration)

  ! <<< Derived types >>>
  use IgnitionSource_mod, only : t_IgnitionSource

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionSource) :: this
  real(SCALAR_KIND), intent(in) :: location(:), amplitude, radius, timeStart,                &
       timeDuration

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  assert(size(location) >= 1 .and. size(location) <= 3)
  assert(radius > 0.0_wp)
  assert(timeStart > 0.0_wp)
  assert(timeDuration > 0.0_wp)

  this%location = 0.0_wp
  this%location(1:size(location)) = location

  this%amplitude = amplitude

  this%radius = radius

  this%timeStart = timeStart
  this%timeDuration = timeDuration

end subroutine setupIgnitionSource

subroutine addIgnitionSource(this, time, coordinates, iblank, ratioOfSpecificHeats,          &
     heatRelease, rightHandSide)

  ! <<< Derived types >>>
  use IgnitionSource_mod, only : t_IgnitionSource

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionSource) :: this
  real(SCALAR_KIND), intent(in) :: time, ratioOfSpecificHeats, heatRelease
  SCALAR_TYPE, intent(in) :: coordinates(:,:)
  integer, intent(in) :: iblank(:)
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  real(wp) :: r, power, timePortion, referenceTemperature, flameTemperature, gaussianFactor

  nDimensions = size(coordinates, 2)
  assert_key(nDimensions, (1, 2, 3))

  if (this%timeDuration > 0.0_wp) then
     timePortion = exp( -0.5_wp * (time - this%timeStart)**2 / this%timeDuration**2) /       &
          this%timeDuration
  else
     timePortion = 1.0_wp
  end if

  referenceTemperature = 1.0_wp / (ratioOfSpecificHeats - 1.0_wp)

  flameTemperature = referenceTemperature / (1.0_wp - heatRelease)

  power = 0.5_wp * this%amplitude * heatRelease * flameTemperature / sqrt(2.0_wp * pi)

  gaussianFactor = 0.5_wp / this%radius**2

  do i = 1, size(rightHandSide, 1)

     if (iblank(i) == 0) cycle

     r = real(sum((coordinates(i,:) - this%location(1:nDimensions)) ** 2), wp)

     rightHandSide(i,nDimensions+2) = rightHandSide(i,nDimensions+2) +                       &
          power * timePortion * exp(- gaussianFactor * r)

  end do

end subroutine addIgnitionSource
