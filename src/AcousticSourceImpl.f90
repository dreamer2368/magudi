#include "config.h"

subroutine setupAcousticSource(this, amplitude, frequency, radius, x, y, z, phase)

  ! <<< Derived types >>>
  use AcousticSource_type

  implicit none

  ! <<< Arguments >>>
  type(t_AcousticSource) :: this
  real(SCALAR_KIND), intent(in) :: amplitude, frequency, radius
  real(SCALAR_KIND), intent(in), optional :: x, y, z, phase

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)

  this%location = 0.0_wp
  if (present(x)) this%location(1) = x
  if (present(y)) this%location(2) = y
  if (present(z)) this%location(3) = z

  this%amplitude = amplitude
  this%gaussianFactor = 9.0_wp / (2.0_wp * radius ** 2)
  this%angularFrequency = 2.0_wp * pi * frequency

  if (present(phase)) this%phase = phase

end subroutine setupAcousticSource

subroutine addAcousticSource(this, time, coordinates, iblank, rightHandSide)

  ! <<< Derived types >>>
  use AcousticSource_type

  implicit none

  ! <<< Arguments >>>
  type(t_AcousticSource) :: this
  real(SCALAR_KIND), intent(in) :: time
  SCALAR_TYPE, intent(in) :: coordinates(:,:)
  integer, intent(in) :: iblank(:)
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions
  real(wp) :: a, r

  nDimensions = size(coordinates, 2)

  a = this%amplitude * cos(this%angularFrequency * time + this%phase)
  do i = 1, size(rightHandSide, 1)
     if (iblank(i) == 0) cycle
     r = real(sum((coordinates(i,:) - this%location(1:nDimensions)) ** 2), wp)
     rightHandSide(i,nDimensions+2) = rightHandSide(i,nDimensions+2) +                       &
          a * exp(-this%gaussianFactor * r)
  end do

end subroutine addAcousticSource
