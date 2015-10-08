#include "config.h"

subroutine setupIgnitionSource(this, ratioOfSpecificHeats, location, radius, amplitude,      &
     timeStart, timeDuration, shockMach)

  ! <<< Derived types >>>
  use IgnitionSource_mod, only : t_IgnitionSource

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionSource) :: this
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats, location(:), radius(:), amplitude, &
       timeStart, timeDuration, shockMach

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp) :: M0, M02, g, gm1, gp1

  assert(size(location) >= 1 .and. size(location) <= 3)
  assert(size(radius) >= 1 .and. size(radius) <= 3)

  this%location = 0.0_wp
  this%location(1:size(location)) = location

  this%amplitude = amplitude

  this%radius = 0.0_wp
  this%radius(1:size(radius)) = radius

  this%timeStart = timeStart
  this%timeDuration = timeDuration

  this%shockMach = shockMach
  if (this%shockMach > 0.0_wp) then
     this%vorticityDeposition = .true.
     M0 = this%shockMach
     M02 = M0 * M0
     g = ratioOfSpecificHeats
     gm1 = g - 1.0_wp
     gp1 = g + 1.0_wp
     this%vorticityCorrection = ( 2.0_wp * gm1**2 * (M02 - 1.0_wp) - gp1**3 * M02 *          &
          log(gp1 * M02) + gp1**2 * (gm1 * M02 + 2.0_wp) * log(gm1 * M02 + 2.0_wp) -         &
          8.0_wp * (1.0_wp - 3.0_wp * g) * M02 * log(M0) ) / ( 2.0_wp * gm1**2 * gp1 *       &
          M0 * (M02 - 1.0_wp) )
  end if

end subroutine setupIgnitionSource

subroutine addIgnitionSource(this, time, coordinates, iblank, density, ratioOfSpecificHeats, &
     heatRelease, rightHandSide)

  ! <<< Derived types >>>
  use IgnitionSource_mod, only : t_IgnitionSource

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionSource) :: this
  real(SCALAR_KIND), intent(in) :: time, ratioOfSpecificHeats, heatRelease
  SCALAR_TYPE, intent(in) :: coordinates(:,:), density(:)
  integer, intent(in) :: iblank(:)
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  real(wp) :: power, timePortion, referenceTemperature, flameTemperature, gaussianFactor(3), &
       vorticityContribution, r, z, M0, M02, g, gm1, gp1

  nDimensions = size(coordinates, 2)
  assert_key(nDimensions, (1, 2, 3))
  assert(size(rightHandSide, 2) >= nDimensions+2)
  assert(size(this%location) >= nDimensions)
  assert(size(this%radius) >= nDimensions)

  if (this%timeDuration > 0.0_wp) then
     timePortion = exp( -0.5_wp * (time - this%timeStart)**2 / this%timeDuration**2) /       &
          this%timeDuration / sqrt(2.0_wp * pi)
  else
     timePortion = 1.0_wp
  end if

  referenceTemperature = 1.0_wp / (ratioOfSpecificHeats - 1.0_wp)

  flameTemperature = referenceTemperature / (1.0_wp - heatRelease)

  power = 0.5_wp * this%amplitude * heatRelease * flameTemperature

  gaussianFactor = 0.0_wp
  if (minval(this%radius(1:nDimensions)) > 0.0_wp) gaussianFactor(1:nDimensions) =           &
       0.5_wp / this%radius(1:nDimensions)**2

  if (this%vorticityDeposition) then
     M0 = this%shockMach
     M02 = M0 * M0
     g = ratioOfSpecificHeats
     gm1 = g - 1.0_wp
     gp1 = g + 1.0_wp
  end if

  do i = 1, size(rightHandSide, 1)
     if (iblank(i) == 0) cycle
     rightHandSide(i,nDimensions+2) = rightHandSide(i,nDimensions+2) +                       &
          power * timePortion * exp(- sum(gaussianFactor(1:nDimensions) *                    &
          (coordinates(i,:) - this%location(1:nDimensions))**2) )

     if (this%vorticityDeposition .and. nDimensions > 1) then
        if (nDimensions == 2) then
           r = abs(coordinates(i,1) - this%location(1)) / this%radius(1)
        else
           r = sqrt( (coordinates(i,1) - this%location(1))**2 + (coordinates(i,3) -          &
                this%location(3))**2 ) / this%radius(1)
        end if

        z = abs(coordinates(i,2) - this%location(2)) / this%radius(2)

        if (r >= sqrt(M02 - 1.0_wp) .or. z >= 1.0_wp) cycle

        vorticityContribution = - this%vorticityCorrection -                                 &
             ( (-2.0_wp * gm1 * M02) / (r**2 + 1.0_wp) - 4.0_wp * g * log(r**2 + 1.0_wp) +   &
             gp1**2 * log(gm1 * M02 + 2.0_wp * (r**2 + 1.0_wp)) ) / (gm1**2 * gp1 * M0)

        rightHandSide(i,3) = rightHandSide(i,3) + timePortion * density(i) *                 &
             vorticityContribution
     end if
  end do

end subroutine addIgnitionSource
