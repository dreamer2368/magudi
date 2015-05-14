#include "config.h"

subroutine computeUniformMap(coordinate, gridIndex, direction, minMaxRange, isPeriodic)

  implicit none

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(out) :: coordinate(:)
  integer, intent(in), optional :: gridIndex, direction
  real(SCALAR_KIND), intent(in), optional :: minMaxRange(2)
  logical, intent(in), optional :: isPeriodic

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  logical :: isPeriodic_
  integer :: i, n
  real(wp) :: h

  assert(size(coordinate) >= 2)

  isPeriodic_ = .false.
  if (present(isPeriodic)) isPeriodic_ = isPeriodic

  n = size(coordinate)

  ! Find mesh spacing.
  if (isPeriodic_) then
     h = 1.0_wp / real(n, wp)
  else
     h = 1.0_wp / real(n - 1, wp)
  end if

  if (present(minMaxRange)) then
     h = h * (minMaxRange(2) - minMaxRange(1))
     coordinate(1) = minMaxRange(1)
  end if

  do i = 2, n
     coordinate(i) = coordinate(i-1) + h
  end do

end subroutine computeUniformMap
