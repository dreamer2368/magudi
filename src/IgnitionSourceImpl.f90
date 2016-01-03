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
  this%depositVorticity = .false.
  if (this%shockMach > 1.0_wp) this%depositVorticity = .true.

end subroutine setupIgnitionSource

subroutine addIgnitionSource(this, time, coordinates, iblank, density, ratioOfSpecificHeats, &
     heatRelease, rightHandSide)

  ! <<< Derived types >>>
  use IgnitionSource_mod, only : t_IgnitionSource

  ! <<< Internal modules >>>
  use MathHelper, only : pi

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionSource) :: this
  real(SCALAR_KIND), intent(in) :: time, ratioOfSpecificHeats, heatRelease
  SCALAR_TYPE, intent(in) :: coordinates(:,:), density(:)
  integer, intent(in) :: iblank(:)
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions, vorticityExponent
  real(wp) :: power, timePortion, referenceTemperature, flameTemperature, gaussianFactor(3), &
       vorticityContribution, vorticityCoefficient, vorticityLocation(3)

  nDimensions = size(coordinates, 2)
  assert_key(nDimensions, (1, 2, 3))
  assert(size(rightHandSide, 2) >= nDimensions+2)
  assert(size(this%location) >= nDimensions)
  assert(size(this%radius) >= nDimensions)

  referenceTemperature = 1.0_wp / (ratioOfSpecificHeats - 1.0_wp)

  flameTemperature = referenceTemperature / (1.0_wp - heatRelease)

  power = 0.5_wp * this%amplitude * heatRelease * flameTemperature

  if (this%timeDuration > 0.0_wp) then
     timePortion = exp( -0.5_wp * (time - this%timeStart)**2 / this%timeDuration**2)
     power = power / this%timeDuration / sqrt(2.0_wp * pi)
  else
     timePortion = 1.0_wp
  end if

  gaussianFactor = 0.0_wp
  if (minval(this%radius(1:nDimensions)) > 0.0_wp) gaussianFactor(1:nDimensions) =           &
       0.5_wp / this%radius(1:nDimensions)**2

  if (this%depositVorticity) then

     vorticityLocation = 0.0_wp
     vorticityLocation(1:nDimensions) = this%location(1:nDimensions)

     select case (nDimensions)
     case(2)
        vorticityCoefficient = 0.017_wp * (this%shockMach**4 - 1.0_wp)
        vorticityExponent = 1
     case(3)
        vorticityCoefficient = 0.022_wp * (this%shockMach**4 - 1.0_wp)
        vorticityExponent = 2
     end select
  end if

  do i = 1, size(rightHandSide, 1)
     if (iblank(i) == 0) cycle
     rightHandSide(i,nDimensions+2) = rightHandSide(i,nDimensions+2) +                       &
          power * timePortion * exp(- sum(gaussianFactor(1:nDimensions) *                    &
          (coordinates(i,:) - this%location(1:nDimensions))**2) )

     if (this%depositVorticity .and. nDimensions > 1) then

        vorticityLocation(2) = coordinates(i,2)
        vorticityContribution = vorticityCoefficient *                                       &
             (1.0_wp - 2.0_wp * sum(gaussianFactor(1:nDimensions) *                          &
             (coordinates(i,:) - vorticityLocation(1:nDimensions))**2)) *                    &
             exp(- sum(gaussianFactor(1:nDimensions) *                                       &
             (coordinates(i,:) - this%location(1:nDimensions))**2) ) ** vorticityExponent

        rightHandSide(i,3) = rightHandSide(i,3) + timePortion * density(i) *                 &
             vorticityContribution
     end if
  end do

end subroutine addIgnitionSource

subroutine addAdjointIgnitionSource(this, time, coordinates, iblank, adjointVariables,       &
     ratioOfSpecificHeats, rightHandSide)

  ! <<< Derived types >>>
  use IgnitionSource_mod, only : t_IgnitionSource

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionSource) :: this
  real(SCALAR_KIND), intent(in) :: time, ratioOfSpecificHeats
  SCALAR_TYPE, intent(in) :: coordinates(:,:), adjointVariables(:,:)
  integer, intent(in) :: iblank(:)
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions, nUnknowns, vorticityExponent
  real(wp) :: gaussianFactor(3), timePortion, vorticityContribution, vorticityCoefficient,   &
       vorticityLocation(3)
  SCALAR_TYPE, allocatable :: localSourceJacobian(:,:), temp(:)

  nDimensions = size(coordinates, 2)
  nUnknowns = size(rightHandSide, 2)
  assert_key(nDimensions, (1, 2, 3))
  assert(size(rightHandSide, 2) >= nDimensions+2)
  assert(size(this%location) >= nDimensions)
  assert(size(this%radius) >= nDimensions)

  if (.not.this%depositVorticity .or. nDimensions == 1) return

  allocate(localSourceJacobian(nUnknowns, nUnknowns))
  allocate(temp(nUnknowns))

  localSourceJacobian = 0.0_wp

  if (this%timeDuration > 0.0_wp) then
     timePortion = exp( -0.5_wp * (time - this%timeStart)**2 / this%timeDuration**2)
  else
     timePortion = 1.0_wp
  end if

  gaussianFactor = 0.0_wp
  if (minval(this%radius(1:nDimensions)) > 0.0_wp) gaussianFactor(1:nDimensions) =           &
       0.5_wp / this%radius(1:nDimensions)**2

  vorticityLocation = 0.0_wp
  vorticityLocation(1:nDimensions) = this%location(1:nDimensions)

  select case (nDimensions)
  case(2)
     vorticityCoefficient = 0.017_wp * (this%shockMach**4 - 1.0_wp)
     vorticityExponent = 1
  case(3)
     vorticityCoefficient = 0.022_wp * (this%shockMach**4 - 1.0_wp)
     vorticityExponent = 2
  end select

  do i = 1, size(rightHandSide, 1)
     if (iblank(i) == 0) cycle

     vorticityLocation(2) = coordinates(i,2)
     vorticityContribution = vorticityCoefficient *                                          &
          (1.0_wp - 2.0_wp * sum(gaussianFactor(1:nDimensions) *                             &
          (coordinates(i,:) - vorticityLocation(1:nDimensions))**2)) *                       &
          exp(- sum(gaussianFactor(1:nDimensions) *                                          &
          (coordinates(i,:) - this%location(1:nDimensions))**2) ) ** vorticityExponent

     localSourceJacobian(3,1) = timePortion * vorticityContribution

     temp = matmul(transpose(localSourceJacobian), adjointVariables(i,:))

     rightHandSide(i,:) = rightHandSide(i,:) - temp
  end do

  SAFE_DEALLOCATE(localSourceJacobian)
  SAFE_DEALLOCATE(temp)

end subroutine addAdjointIgnitionSource
