#include "config.h"

module MathHelper

  implicit none
  public

  integer, parameter, private :: wp = SCALAR_KIND

  ! Trigonometric parameters.
  real(SCALAR_KIND), parameter :: Pi    = 3.1415926535897932385_wp
  real(SCALAR_KIND), parameter :: twoPi = 6.2831853071795864770_wp

  ! Bessel first zero.
  real(SCALAR_KIND), parameter :: besselJ1Zero = 3.8317059702075123115_wp

  interface

     pure function crossProduct(x, y) result(z)

       !>  Returns cross product in 3 dimensions: z = cross(x, y).
    
       real(SCALAR_KIND), dimension(3), intent(in) :: x, y
       real(SCALAR_KIND), dimension(3)             :: z

     end function crossProduct

  end interface

  interface

     pure function normalize(v) result(w)

       !> Returns normalized vector: w = v / |v|.
    
       real(SCALAR_KIND), dimension(3), intent(in) :: v
       real(SCALAR_KIND), dimension(3)             :: w

     end function normalize

  end interface

  interface

     function gamma(xx)

       !> Returns the gamma function.
    
       real(SCALAR_KIND) :: gamma
       real(SCALAR_KIND), intent(in) :: xx
    
     end function gamma

  end interface
  
  interface

     function gammaLn(xx)

    !> Returns the log of the gamma function.
  
       real(SCALAR_KIND) :: gammaLn
       real(SCALAR_KIND), intent(in) :: xx

     end function gammaLn

  end interface
  
  interface

     function besselJ0(x)

       !> Returns the Bessel function J0(x) for any real x.
       !> [Numerical Recipes in Fortran]

       real(SCALAR_KIND) :: besselJ0
       real(SCALAR_KIND), intent(in) :: x

     end function besselJ0

  end interface

  interface

     function besselJ1(x)

       !> Returns the Bessel function J1(x) for any real x.
       !> [Numerical Recipes in Fortran]

       real(SCALAR_KIND) :: besselJ1
       real(SCALAR_KIND), intent(in) :: x

     end function besselJ1

  end interface

  interface
     function arctan(dx, dy)
    
       real(SCALAR_KIND), intent(in) :: dx, dy
       real(SCALAR_KIND) :: arctan

     end function arctan

  end interface

end module MathHelper
