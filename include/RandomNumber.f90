#include "config.h"

module RandomNumber

  interface

     subroutine initializeRandomNumberGenerator()

       !> Initializes the random number generator using a seed based on the system clock and a
       !> seed read from `/dev/urandom`.

     end subroutine initializeRandomNumberGenerator

  end interface

  interface random

     function randomInteger_(a, b) result(x)

       integer, intent(in) :: a, b

       integer :: x

     end function randomInteger_

     function randomReal_(a, b) result(x)

       real(SCALAR_KIND), intent(in) :: a, b

       real(SCALAR_KIND) :: x

     end function randomReal_

  end interface random

  private :: randomInteger_, randomReal_

end module RandomNumber
