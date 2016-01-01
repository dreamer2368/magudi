#include "config.h"

PURE_FUNCTION crossProduct(x, y) result(z)

  implicit none

  ! <<< Arguments >>>
  real(SCALAR_KIND), dimension(3), intent(in) :: x, y

  ! <<< Result >>>
  real(SCALAR_KIND), dimension(3) :: z

  z(1) = x(2) * y(3) - x(3) * y(2)
  z(2) = x(3) * y(1) - x(1) * y(3)
  z(3) = x(1) * y(2) - x(2) * y(1)

end function crossProduct

PURE_FUNCTION normalize(v) result(w)

  implicit none

  ! <<< Arguments >>>
  real(SCALAR_KIND), dimension(3), intent(in) :: v

  ! <<< Result >>>
  real(SCALAR_KIND), dimension(3) :: w

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp) :: norm

  norm = 1.0_wp / (sqrt(dot_product(v, v)) + epsilon(1.0_wp))
  w = v * norm

end function normalize

function besselJ0(x)

  implicit none

  ! <<< Arguments >>>
  real(SCALAR_KIND) :: besselJ0
  real(SCALAR_KIND), intent(in) :: x

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp), dimension(6), parameter :: r = (/57568490574.0_wp, -13362590354.0_wp,            &
       651619640.7_wp, -11214424.18_wp,77392.33017_wp,-184.9052456_wp/)
  real(wp), dimension(6), parameter :: s = (/57568490411.0_wp, 1029532985.0_wp,              &
       9494680.718_wp, 59272.64853_wp, 267.8532712_wp, 1.0_wp/)
  real(wp), dimension(5), parameter :: p = (/1.0_wp, -.1098628627E-2_wp, 0.2734510407E-4_wp, &
       -0.2073370639E-5_wp, 0.2093887211E-6_wp/)
  real(wp), dimension(5), parameter :: q = (/-0.1562499995E-1_wp, 0.1430488765E-3_wp,        &
       -0.6911147651E-5_wp, 0.7621095161E-6_wp, -0.934945152E-7_wp/)
  real(wp) :: ax, xx, z, y
    
  if(abs(x) < 8.0_wp) then
     y = x ** 2
     besselJ0 = (r(1) + y * (r(2) + y*(r(3) + y * (r(4) + y * (r(5) + y * r(6)))))) /          &
          (s(1) + y * (s(2) + y * (s(3) + y * (s(4) + y * (s(5) + y * s(6))))))
  else
     ax = abs(x)
     z = 8.0_wp / ax
     y = z ** 2
     xx = ax - 0.785398164_wp
     besselJ0 = sqrt(0.636619772_wp / ax) * (cos(xx) * (p(1) + y * (p(2) + y* (p(3) + y *    &
          (p(4) + y * p(5))))) - z * sin(xx) * (q(1) + y * (q(2) + y * (q(3) + y * (q(4) +   &
          y * q(5))))))
  end if

end function besselJ0

function besselJ1(x)

  implicit none

  ! <<< Arguments >>>
  real(SCALAR_KIND) :: besselJ1
  real(SCALAR_KIND), intent(in) :: x

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp), dimension(6), parameter :: r = (/72362614232.0_wp, -7895059235.0_wp,             &
       242396853.1_wp, -2972611.439_wp, 15704.48260_wp, -30.16036606_wp/)
  real(wp), dimension(6), parameter :: s = (/144725228442.0_wp, 2300535178.0_wp,             &
       18583304.74_wp, 99447.43394_wp, 376.9991397_wp, 1.0_wp/)
  real(wp), dimension(5), parameter :: p = (/1.0_wp, 0.183105E-2_wp, 0-.3516396496E-4_wp,    &
       0.2457520174E-5_wp, -0.240337019E-6_wp/)
  real(wp), dimension(5), parameter :: q = (/0.04687499995_wp, -0.2002690873E-3_wp,          &
       0.8449199096E-5_wp, -.88228987E-6_wp, 0.105787412E-6_wp/)
  real(wp) :: ax, xx, z, y
    
  if (abs(x) < 8.0_wp) then
     y = x ** 2
     besselJ1 = x * (r(1) + y * (r(2) + y * (r(3) + y * (r(4) + y * (r(5) + y * r(6)))))) /  &
          (s(1) + y * (s(2) + y * (s(3) + y * (s(4) + y * (s(5) + y * s(6))))))
  else
     ax = abs(x)
     z = 8.0_wp / ax
     y = z ** 2
     xx = ax - 2.356194491_wp
     besselJ1 = sqrt(0.636619772_wp / ax) * (cos(xx) * (p(1) + y * (p(2) + y * (p(3) + y *   &
          (p(4) + y * p(5))))) - z * sin(xx) * (q(1) + y * (q(2) + y * (q(3) + y * (q(4) +   &
          y * q(5)))))) * sign(1.0_wp, x)
  end if

end function besselJ1

function arctan(dx, dy)

  ! <<< Public members >>>
  use MathHelper, only : Pi, twoPi

  implicit none

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(in) :: dx, dy
  real(SCALAR_KIND) :: arctan

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  ! Evaluate atan.
  if (abs(dx) + abs(dy) < 1.0e-9_wp) then
     arcTan = 0.0_wp
  else
     arctan = atan(dy / dx)
  end if
    
  ! Quadrant correction.
  if (dx <= 0.0_wp) then
     arctan = Pi + arctan
  elseif (dy <= 0.0_wp .and. dx > 0.0_wp) then
     arctan = twoPi + arctan
  end if

end function arctan
