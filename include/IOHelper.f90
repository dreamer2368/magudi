#include "config.h"

module IOHelper

  implicit none
  public

  integer :: fileIndex
  integer, dimension(128) :: iunits

  interface

     integer function iOpen()

     end function iOpen

  end interface

  interface

     integer function iClose(iu)

       integer, intent(inout) :: iu

     end function iClose

  end interface

end module IOHelper
