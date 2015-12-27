#include "config.h"

module Parser

  implicit none
  public

  integer, parameter :: tagLength = 128
  integer, parameter :: lineLength = 409600
  integer, dimension(128) :: iunits
  integer :: nFields, fileIndex

  type t_entryType
     character(tagLength) :: tag
     character(lineLength) :: value
  end type t_entryType

  type(t_entryType), pointer, dimension(:) :: entries

  interface parserRead

     subroutine parserReadLogical(myTag, value, Default)

       character(*), intent(in) :: myTag
       logical, intent(out) :: value
       logical, optional, intent(in) :: default

     end subroutine parserReadLogical

     subroutine parserReadInt(myTag, value, default)

       character(*), intent(in) :: myTag
       integer, intent(out) :: value
       integer, optional, intent(in) :: default

     end subroutine parserReadInt

     subroutine parserReadIntArray(myTag, value)

       character(*), intent(in) :: myTag
       integer, dimension(:), intent(out) :: value

     end subroutine parserReadIntArray

     subroutine parserReadFloat(myTag, value, default)

       character(*), intent(in) :: myTag
       real(SCALAR_KIND), intent(out) :: value
       real(SCALAR_KIND), optional, intent(in) :: default

     end subroutine parserReadFloat

     subroutine parserReadFloatArray(myTag, value)

       character(*), intent(in) :: mytag
       real(SCALAR_KIND), dimension(:), intent(out) :: value

     end subroutine parserReadFloatArray

     subroutine parserReadFloatArray2D(myTag, value)

       character(*), intent(in) :: myTag
       real(SCALAR_KIND), dimension(:,:), intent(out) :: value
    
     end subroutine parserReadFloatArray2D

     subroutine parserReadChar(myTag, value, default)

       character(*), intent(in) :: myTag
       character(len=*), intent(out) :: value
       character(len=*), optional, intent(in) :: default

     end subroutine parserReadChar
                                
     subroutine parserReadCharArray(myTag, value)

       character(*), intent(in) :: myTag
       character(*), dimension(:), intent(out) :: value

     end subroutine parserReadCharArray

  end interface parserRead

  interface

     subroutine setupParser

     end subroutine setupParser

  end interface

  interface

     subroutine parseFile(input)

       character(len = *), intent(in) :: input

     end subroutine parseFile

  end interface

  interface

     subroutine parserIsDefined(myTag, isDef)

       character(*), intent(in) :: myTag
       logical, intent(out) :: isDef

     end subroutine parserIsDefined
  
  end interface

  interface
     subroutine parserGetSize(myTag, numb)

       character(*), intent(in) :: myTag
       integer, intent(out) :: numb

     end subroutine parserGetSize

  end interface

  private :: parserReadLogical, parserReadInt, parserReadIntArray, parserReadFloat,          &
       parserReadFloatArray, parserReadFloatArray2D, parserReadChar, parserReadCharArray

end module Parser
