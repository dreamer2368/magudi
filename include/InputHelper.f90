#include "config.h"

module InputHelper

  implicit none
  public

  interface

     subroutine parseInputFile(filename, commentMarker, separator)

       character(len = *), intent(in) :: filename

       character(len = 1), intent(in), optional :: commentMarker, separator

     end subroutine parseInputFile

  end interface

  interface

     subroutine stripComments(str, commentMarker)

       character(len = *), intent(inout) :: str
       character(len = 1), intent(in) :: commentMarker

     end subroutine stripComments

  end interface

  interface getOption

     function getOptionInteger_(key, defaultValue) result(val)

       character(len = *), intent(in) :: key
       integer, intent(in) :: defaultValue

       integer :: val

     end function getOptionInteger_

     function getOptionLogical_(key, defaultValue) result(val)

       character(len = *), intent(in) :: key
       logical, intent(in) :: defaultValue

       logical :: val

     end function getOptionLogical_

     function getOptionScalar_(key, defaultValue) result(val)

       character(len = *), intent(in) :: key
       SCALAR_TYPE, intent(in) :: defaultValue

       SCALAR_TYPE :: val

     end function getOptionScalar_

     function getOptionString_(key, defaultValue) result(val)

       character(len = *), intent(in) :: key, defaultValue

       character(len = STRING_LENGTH) :: val

     end function getOptionString_

#ifdef SCALAR_IS_COMPLEX
     function getOptionReal_(key, defaultValue) result(val)

       character(len = *), intent(in) :: key
       real(SCALAR_KIND), intent(in) :: defaultValue

       real(SCALAR_KIND) :: val

     end function getOptionReal_
#endif

  end interface getOption

  interface getRequiredOption

     subroutine getRequiredOptionInteger_(key, val, comm)

       character(len = *), intent(in) :: key
       integer, intent(out) :: val

       integer, intent(in), optional :: comm

     end subroutine getRequiredOptionInteger_

     subroutine getRequiredOptionLogical_(key, val, comm)

       character(len = *), intent(in) :: key
       logical, intent(out) :: val

       integer, intent(in), optional :: comm

     end subroutine getRequiredOptionLogical_

     subroutine getRequiredOptionScalar_(key, val, comm)

       character(len = *), intent(in) :: key
       SCALAR_TYPE, intent(out) :: val

       integer, intent(in), optional :: comm

     end subroutine getRequiredOptionScalar_

     subroutine getRequiredOptionString_(key, val, comm)

       character(len = *), intent(in) :: key
       character(len = STRING_LENGTH), intent(out) :: val

       integer, intent(in), optional :: comm

     end subroutine getRequiredOptionString_

#ifdef SCALAR_IS_COMPLEX
     subroutine getRequiredOptionReal_(key, val, comm)

       character(len = *), intent(in) :: key
       real(SCALAR_KIND), intent(out) :: val

       integer, intent(in), optional :: comm

     end subroutine getRequiredOptionReal_
#endif

  end interface getRequiredOption

  private :: getOptionInteger_, getOptionLogical_, getOptionScalar_, getOptionString_,       &
       getRequiredOptionInteger_, getRequiredOptionLogical_, getRequiredOptionScalar_,       &
       getRequiredOptionString_
#ifdef SCALAR_IS_COMPLEX
  private :: getOptionReal_, getRequiredOptionReal_
#endif

end module InputHelper
