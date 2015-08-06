#include "config.h"

module AcousticSource_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  implicit none

  type, public :: t_AcousticSource

     real(SCALAR_KIND) :: location(3), amplitude, gaussianFactor, angularFrequency, phase

   contains

     procedure, public, pass :: setup => setupAcousticSource
     procedure, public, pass :: add => addAcousticSource

  end type t_AcousticSource

contains

  subroutine setupAcousticSource(this, location, amplitude, frequency, radius, phase)

    implicit none

    ! <<< Arguments >>>
    class(t_AcousticSource) :: this
    real(SCALAR_KIND), intent(in) :: location(:), amplitude, frequency, radius
    real(SCALAR_KIND), intent(in), optional :: phase

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)

    assert(size(location) >= 1 .and. size(location) <= 3)
    assert(radius > 0.0_wp)

    this%location = 0.0_wp
    this%location(1:size(location)) = location

    this%amplitude = amplitude
    this%angularFrequency = 2.0_wp * pi * frequency
    this%gaussianFactor = 9.0_wp / (2.0_wp * radius ** 2)

    this%phase = 0.0_wp
    if (present(phase)) this%phase = phase

  end subroutine setupAcousticSource

  subroutine addAcousticSource(this, time, coordinates, iblank, rightHandSide)

    implicit none

    ! <<< Arguments >>>
    class(t_AcousticSource) :: this
    real(SCALAR_KIND), intent(in) :: time
    real(SCALAR_KIND), intent(in) :: coordinates(:,:)
    integer, intent(in) :: iblank(:)
    real(SCALAR_KIND), intent(inout) :: rightHandSide(:,:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, nDimensions
    real(wp) :: a, r

    nDimensions = size(coordinates, 2)
    assert_key(nDimensions, (1, 2, 3))

    a = this%amplitude * cos(this%angularFrequency * time + this%phase)
    do i = 1, size(rightHandSide, 1)
       if (iblank(i) == 0) cycle
       r = real(sum((coordinates(i,:) - this%location(1:nDimensions)) ** 2), wp)
       rightHandSide(i,nDimensions+2) = rightHandSide(i,nDimensions+2) +                     &
            a * exp(-this%gaussianFactor * r)
    end do

  end subroutine addAcousticSource

end module AcousticSource_mod
