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
       vorticityContribution, vorticityCorrection, halfWidth(3), r, z, r2p1, M0, M02, g, gm1,&
       gp1, rlim

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
     M0 = this%shockMach
     M02 = M0 * M0
     g = ratioOfSpecificHeats
     gm1 = g - 1.0_wp
     gp1 = g + 1.0_wp
     rlim = sqrt(M02 - 1.0_wp)
     halfWidth = 2.0_wp * sqrt(2.0_wp * log(2.0_wp)) * this%radius
     if (nDimensions == 2) then
        vorticityCorrection =  - ( (2.0_wp * M0 * atan(rlim)) / ((g**2 - 1.0_wp) * rlim) -   &
             ((gp1 * (sqrt(4.0_wp + 2.0_wp * gm1 * M02) * atan(sqrt(2.0_wp) *                &
             sqrt((M02 - 1.0_wp) / (2.0_wp + gm1 * M02))) + rlim *                           &
             (log(gp1) + 2.0_wp * log(M0) - 2.0_wp))) / (gm1**2 * M0 * rlim)) +              &
             (8.0_wp * g * (atan(rlim) / rlim + log(M0) - 1.0_wp)) / (gm1**2 * gp1 * M0) )
     else if(nDimensions == 3) then
        vorticityCorrection = - ( 2.0_wp * gm1**2 * (M02 - 1.0_wp) - gp1**3 * M02 *          &
             log(gp1 * M02) + gp1**2 * (gm1 * M02 + 2.0_wp) * log(gm1 * M02 + 2.0_wp) -      &
             8.0_wp * (1.0_wp - 3.0_wp * g) * M02 * log(M0) ) / ( 2.0_wp * gm1**2 * gp1 *    &
             M0 * (M02 - 1.0_wp) )
     end if
  end if

  do i = 1, size(rightHandSide, 1)
     if (iblank(i) == 0) cycle
     rightHandSide(i,nDimensions+2) = rightHandSide(i,nDimensions+2) +                       &
          power * timePortion * exp(- sum(gaussianFactor(1:nDimensions) *                    &
          (coordinates(i,:) - this%location(1:nDimensions))**2) )

     if (this%depositVorticity .and. nDimensions > 1) then
        if (nDimensions == 2) then
           r = abs(coordinates(i,1) - this%location(1)) / halfWidth(1)
        else
           r = sqrt( (coordinates(i,1) - this%location(1))**2 + (coordinates(i,3) -          &
                this%location(3))**2 ) / halfWidth(1)
        end if

        z = abs(coordinates(i,2) - this%location(2)) / halfWidth(2)

        if (r >= 1.0_wp .or. z >= 1.0_wp) cycle

        r = r * rlim
        r2p1 = r * r + 1.0_wp

        vorticityContribution = vorticityCorrection +                                        &
             2.0_wp * M0 / (r2p1 * (g**2 - 1.0_wp)) -                                        &
             gp1 / (M0 * gm1**2) * log(gm1 * M02 + 2.0_wp * r2p1) +                          &
             4.0_wp * g * log(r2p1) / (gm1**2 * gp1 * M0)

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
  integer :: i, nDimensions, nUnknowns
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  real(wp) :: timePortion, vorticityContribution, vorticityCorrection, halfWidth(3), r, z,   &
       M0, M02, g, gm1, gp1, rlim, r2p1
  SCALAR_TYPE, allocatable :: localSourceJacobian(:,:), temp(:)

  if (.not.this%depositVorticity) return

  nDimensions = size(coordinates, 2)
  nUnknowns = size(rightHandSide, 2)
  assert_key(nDimensions, (1, 2, 3))
  assert(size(rightHandSide, 2) >= nDimensions+2)
  assert(size(this%location) >= nDimensions)
  assert(size(this%radius) >= nDimensions)

  allocate(localSourceJacobian(nUnknowns, nUnknowns))
  allocate(temp(nUnknowns))

  localSourceJacobian = 0.0_wp

  if (this%timeDuration > 0.0_wp) then
     timePortion = exp( -0.5_wp * (time - this%timeStart)**2 / this%timeDuration**2)
  else
     timePortion = 1.0_wp
  end if

  M0 = this%shockMach
  M02 = M0 * M0
  g = ratioOfSpecificHeats
  gm1 = g - 1.0_wp
  gp1 = g + 1.0_wp
  rlim = sqrt(M02 - 1.0_wp)
  halfWidth = 2.0_wp * sqrt(2.0_wp * log(2.0_wp)) * this%radius
  if (nDimensions == 2) then
     vorticityCorrection =  - ( (2.0_wp * M0 * atan(rlim)) / ((g**2 - 1.0_wp) * rlim) -      &
          ((gp1 * (sqrt(4.0_wp + 2.0_wp * gm1 * M02) * atan(sqrt(2.0_wp) *                   &
          sqrt((M02 - 1.0_wp) / (2.0_wp + gm1 * M02))) + rlim *                              &
          (log(gp1) + 2.0_wp * log(M0) - 2.0_wp))) / (gm1**2 * M0 * rlim)) +                 &
          (8.0_wp * g * (atan(rlim) / rlim + log(M0) - 1.0_wp)) / (gm1**2 * gp1 * M0) )
  else if(nDimensions == 3) then
     vorticityCorrection = - ( 2.0_wp * gm1**2 * (M02 - 1.0_wp) - gp1**3 * M02 *             &
          log(gp1 * M02) + gp1**2 * (gm1 * M02 + 2.0_wp) * log(gm1 * M02 + 2.0_wp) -         &
          8.0_wp * (1.0_wp - 3.0_wp * g) * M02 * log(M0) ) / ( 2.0_wp * gm1**2 * gp1 *       &
          M0 * (M02 - 1.0_wp) )
  end if

  do i = 1, size(rightHandSide, 1)
     if (iblank(i) == 0) cycle

     if (nDimensions == 2) then
        r = abs(coordinates(i,1) - this%location(1)) / halfWidth(1)
     else
        r = sqrt( (coordinates(i,1) - this%location(1))**2 + (coordinates(i,3) -             &
             this%location(3))**2 ) / halfWidth(1)
     end if

     z = abs(coordinates(i,2) - this%location(2)) / halfWidth(2)

     if (r >= 1.0_wp .or. z >= 1.0_wp) cycle

     r = r * rlim
     r2p1 = r * r + 1.0_wp

     vorticityContribution = vorticityCorrection +                                           &
          2.0_wp * M0 / (r2p1 * (g**2 - 1.0_wp)) -                                           &
          gp1 / (M0 * gm1**2) * log(gm1 * M02 + 2.0_wp * r2p1) +                             &
          4.0_wp * g * log(r2p1) / (gm1**2 * gp1 * M0)

     localSourceJacobian(3,1) = timePortion * vorticityContribution

     temp = matmul(transpose(localSourceJacobian), adjointVariables(i,:))

     rightHandSide(i,:) = rightHandSide(i,:) - temp
  end do

  SAFE_DEALLOCATE(localSourceJacobian)
  SAFE_DEALLOCATE(temp)

end subroutine addAdjointIgnitionSource
