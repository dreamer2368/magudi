#include "config.h"

module IOHelperImpl

  implicit none
  public

end module IOHelperImpl

integer function iOpen()

  ! <<< Public members >>>
  use IOHelper, only : fileIndex, iunits

  implicit none

  integer, save :: iCall = 1
  integer :: i

  if (iCall .eq. 1) then
     fileIndex = 1
     iCall = 0
  end if
  iunits(fileIndex) = 0
  do i = 1, fileIndex
     if (iunits(i) .eq. 0) exit
  end do
  if (i .eq. fileIndex) then
     fileIndex = fileIndex + 1
     if (fileIndex .ge. 128) stop "iOpen: maximum units number exceeded!"
  end if
  iunits(i) = 1
  iOpen = i + 10

end function iOpen

integer function iClose(iu)

  ! <<< Public members >>>
  use IOHelper, only : fileIndex, iunits

  implicit none

  ! <<< Arguments >>>
  integer, intent(inout) :: iu

  iu = iu - 10
  if (iu .gt. 0 .and. iu .lt. fileIndex) then
     iunits(iu) = 0
     iClose = iu + 10
  else
     iClose = -1
  end if

end function iClose
