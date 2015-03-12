#include "config.h"

module RandomNumberImpl

  implicit none

contains

  function lcg(s)

    ! <<< Arguments >>>
    integer(selected_int_kind(16)) :: s

    ! <<< Result >>>
    integer :: lcg

    ! <<< Local variables >>>
    integer, parameter :: int64 = selected_int_kind(16)

    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if

    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))

  end function lcg

end module RandomNumberImpl

subroutine initializeRandomNumberGenerator()

  ! <<< Private members >>>
  use RandomNumberImpl, only : lcg

  implicit none

  ! <<< Local variables >>>
  integer, parameter :: int64 = selected_int_kind(16), fileUnit = 11
  integer :: n, istat, i, dateTime(8)
  integer, dimension(:), allocatable :: seed1, seed2
  integer(int64) :: t

  call random_seed(size = n)

  allocate(seed1(n), seed2(n))

  call system_clock(t)
  if (t == 0) then
     call date_and_time(values = dateTime)
     t = (dateTime(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 +                            &
          dateTime(2) * 31_int64 * 24 * 60 * 60 * 1000 +                                     &
          dateTime(3) * 24_int64 * 60 * 60 * 1000 +                                          &
          dateTime(5) * 60 * 60 * 1000 +                                                     &
          dateTime(6) * 60 * 1000 +                                                          &
          dateTime(7) * 1000 +                                                               &
          dateTime(8)
  end if
  do i = 1, n
     seed1(i) = lcg(t)
  end do

  open(unit = fileUnit, file = "/dev/urandom", access = "stream",                           &
       form = "unformatted", action = "read", status = "old", iostat = istat)
  if (istat == 0) then
     read(fileUnit) seed2
     close(fileUnit)
  else
     seed2(:) = 0
  end if

  seed2 = seed2 + seed1
  call random_seed(put = seed2)

  SAFE_DEALLOCATE(seed1)
  SAFE_DEALLOCATE(seed2)

end subroutine initializeRandomNumberGenerator

function randomInteger_(a, b) result(x)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: a, b

  ! <<< Result >>>
  integer :: x

  ! <<< Local variables >>>
  real :: y

  call random_number(y)
  x = a + floor(real(b - a + 1) * y)

end function randomInteger_

function randomReal_(a, b) result(x)

  implicit none

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(in) :: a, b

  ! <<< Result >>>
  real(SCALAR_KIND) :: x

  call random_number(x)
  x = a + (b - a) * x

end function randomReal_
