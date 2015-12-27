#include "config.h"

module ParserImpl

  implicit none
  public

contains

  subroutine parserPack(myField)

    ! <<< Public members >>>
    use Parser, only : t_entryType, entries, nFields

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: myField

    ! <<< Local variables >>>
    type(t_entryType), pointer, dimension(:) :: entriesNew
    
    allocate(entriesNew(nFields - 1))
    entriesNew(1:myField - 1) = entries(1:myfield - 1)
    entriesNew(myField:nFields - 1) = entries(myField + 1:nFields)
    deallocate(entries)
    nullify(entries)
    entries => entriesNew
    nFields = nFields - 1

  end subroutine parserPack

  subroutine parserSpread

    ! <<< Public members >>>
    use Parser, only : t_entryType, entries, nFields

    implicit none

    ! <<< Local variables >>>
    type(t_entryType), pointer, dimension(:) :: entriesNew
    
    if (nFields .ne. 0) then 
       allocate(entriesNew(nFields + 1))
       entriesNew(1:nFields) = entries(1:nFields)
       deallocate(entries)
       nullify(entries)
       entries => entriesNew
    else
       allocate(entries(1))
    end if
    nFields = nFields + 1

  end subroutine parserSpread
  
  subroutine parserNewEntry(myTag, myValue)

    ! <<< Public members >>>
    use Parser, only : entries, nFields

    implicit none

    ! <<< Arguments >>>
    character(*), intent(in) :: myTag
    character(*), intent(in) :: myValue

    ! <<< Local variables >>>
    integer :: iField
    logical :: isDef 
    
    call parserFieldForTag(myTag, iField, isDef)
    if (.not. isDef) then
       call parserSpread()  
       entries(nFields)%tag = myTag
       entries(nFields)%value = myValue
    else
       entries(iField)%value = myValue
    end if
    
  end subroutine parserNewEntry
  
  subroutine parserFieldForTag(myTag, myField, isDef)

    ! <<< Public members >>>
    use Parser, only : entries, nFields

    implicit none

    ! <<< Arguments >>>
    character(*), intent(in) :: myTag
    integer, intent(out) :: myField
    logical, optional, intent(out) :: isDef

    ! <<< Local variables >>>
    integer :: iField
    
    isdef = .false.
    do iField = 1, nFields
       if (entries(ifield)%tag == myTag) then
          myfield = iField
          isDef = .true.
          return
       end if
    end do

  end subroutine parserFieldForTag

end module ParserImpl

subroutine setupParser

  ! <<< Public members >>>
  use Parser, only : entries, nFields

  implicit none

  nfields = 0
  if (associated(entries)) then
     deallocate(entries)
     nullify(entries)
  end if

end subroutine setupParser

subroutine parseFile(input)

  ! <<< Public members >>>
  use Parser, only : lineLength, tagLength

  ! <<< Private members >>>
  use ParserImpl, only : parserNewEntry

  ! <<< Internal modules >>>
  use IOHelper

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: input

  ! <<< Local variables >>>
  integer :: iunit, ierror, limiter, nLines, i, j, nTags, comment
  integer, dimension(:), allocatable :: limit,line
  character(len=lineLength) :: buffer
  character(len=lineLength), dimension(:), allocatable :: file
  character(len=lineLength) :: value
  character(len=tagLength) :: tag
  
  ! Open the file.
  ierror = 0
  iunit = iOpen()
  open (iunit, file = input, form = 'formatted', status = 'old', iostat = ierror)
  if (ierror .ne. 0) stop 'Parser : unable to open the input file.'
  
  ! Count the number of lines in the file.
  ierror = 0
  nLines = 0
  do while (ierror .eq. 0)
     read(iunit, '(a)', iostat = ierror) buffer
     nLines = nLines + 1
  end do
  rewind(iunit)
  
  ! Allocate to the right size.
  allocate(file(nLines + 1), limit(nLines + 1), line(nLines + 1))
  
  ! Read everything in the buffer.
  ierror = 0
  nLines = 0
  loop: do while (ierror .eq. 0)
     read(iunit, '(a)', iostat = ierror) buffer
     if (ierror.ne.0) exit loop
     ! Remove the tabs.
     do j = 1, lineLength
        if (ichar(buffer(j:j)).EQ.9) buffer(j:j) = ' '
     end do
     ! Find comments.
     comment = scan(buffer, '!#%')
     ! Remove them.
     if (comment.NE.0) buffer(comment:) = ''
     ! Trim.
     buffer = adjustl(buffer)
     ! Add line.
     if (len_trim(buffer).NE.0) then
        nLines = nLines + 1
        file(nLines) = buffer
     end if
  end do loop
  
  ! Close the file.
  close(iunit)
  ierror = iClose(iunit)
  
  ! Get the tags.
  nTags = 0
  do i = 1, nLines
     limiter = index(file(i), ':')
     if (limiter.NE.0) then
        nTags = nTags + 1
        line(nTags) = i
        line(nTags + 1) = nLines + 1
        limit(nTags) = limiter
     end if
  end do
  
  ! Read everything now.
  do i = 1, nTags
     buffer = ''
     do j = line(i), line(i + 1) - 1
        if (j==line(i)) then
           buffer = trim(buffer) // trim(file(j))
        else
           buffer = trim(buffer) // ' ' // trim(file(j))
        end if
     end do
     read(buffer(1:limit(i) - 1),'(a)') tag
     read(buffer(limit(i) + 1:),'(a)') value
     if (len_trim(value).NE.0) then
        value = adjustl(value)
        call parserNewEntry(tag, value)
     end if
  end do

end subroutine parseFile

subroutine parserIsDefined(myTag, isDef)

  ! <<< Private members >>>
  use ParserImpl, only : parserFieldForTag

  implicit none
    
  ! <<< Arguments >>>
  character(*), intent(in) :: myTag
  logical, intent(out) :: isDef

  ! <<< Local variables >>>
  integer :: iField
    
  call parserFieldForTag(myTag, iField, isDef)

end subroutine parserIsDefined
  
subroutine parserGetSize(myTag, numb)

  ! <<< Private members >>>
  use ParserImpl, only : parserFieldForTag

  ! <<< Public members >>>
  use Parser, only : lineLength, entries

  implicit none

  ! <<< Arguments >>>
  character(*), intent(in) :: myTag
  integer, intent(out) :: numb

  ! <<< Local variables >>>
  integer :: iField
  logical :: isDef
  integer :: i
  integer, dimension(lineLength) :: counter
                                                                        
  call parserFieldForTag(myTag, iField, isDef)
  if (.not. isdef) then
     print *, 'Parser : ' // mytag // ' not defined'
     stop
  end if
                                 
  counter = 0
  do i = 1, len_trim(entries(iField)%value)
     if (entries(ifield)%value(i:i).EQ.' ') counter(i) = 1
  end do
  do i = 1 + 1, len_trim(entries(ifield)%value)
     if (counter(i).EQ.1 .AND. counter(i-1).EQ.1) counter(i-1) = 0
  end do
  numb = sum(counter) + 1

end subroutine parserGetSize

subroutine parserReadLogical(myTag, value, Default)

  ! <<< Private members >>>
  use ParserImpl, only : parserFieldForTag

  ! <<< Public members >>>
  use Parser, only : entries

  implicit none
    
  ! <<< Arguments >>>
  character(*), intent(in) :: myTag
  logical, intent(out) :: value
  logical, optional, intent(in) :: default

  ! <<< Local variables >>>
  integer :: iField, conv
  logical :: isdef
    
  call parserFieldForTag(myTag, iField, isDef)
  if (.not.isDef .AND. present(default)) then   
     value = default
  else if (.not.isDef .AND. .not.present(default)) then
     print *, 'Parser : ' // myTag // ' not defined'
     stop
  else
     if (len_trim(entries(iField)%value) .eq. 1) then          
        read(entries(ifield)%value, *) conv
        value = (conv.eq.1)
     else
        read(entries(ifield)%value, *) value
     end if
  end if

end subroutine parserReadLogical
  
subroutine parserReadInt(myTag, value, default)

  ! <<< Private members >>>
  use ParserImpl, only : parserFieldForTag

  ! <<< Public members >>>
  use Parser, only : entries

  implicit none
    
  ! <<< Arguments >>>
  character(*), intent(in) :: myTag
  integer, intent(out) :: value
  integer, optional, intent(in) :: default

  ! <<< Local variables >>>
  integer :: iField
  logical :: isDef

  call parserFieldForTag(myTag, iField, isDef)
  if (.not.isDef .AND. present(default)) then   
     value = default
  else if (.not.isDef .AND. .not.present(default)) then
     print *, 'Parser : ' // mytag // ' not defined'
     stop
  else
     read(entries(iField)%value, *) value
  end if

end subroutine parserReadInt

subroutine parserReadFloat(myTag, value, default)

  ! <<< Private members >>>
  use ParserImpl, only : parserFieldForTag

  ! <<< Public members >>>
  use Parser, only : entries

  implicit none
    
  ! <<< Arguments >>>
  character(*), intent(in) :: myTag
  real(SCALAR_KIND), intent(out) :: value
  real(SCALAR_KIND), optional, intent(in) :: default

  ! <<< Local variables >>>
  integer :: iField
  logical :: isDef
    
  call parserFieldForTag(myTag, iField, isDef)
  if (.not.isDef .AND. present(default)) then   
     value = default
  else if (.not.isDef .AND. .not.present(default)) then
     print *, 'Parser : ' // mytag // ' not defined'
     stop
  else
     read(entries(iField)%value, *) value
  end if

end subroutine parserReadFloat

subroutine parserReadChar(myTag, value, default)

  ! <<< Private members >>>
  use ParserImpl, only : parserFieldForTag

  ! <<< Public members >>>
  use Parser, only : entries

  implicit none
    
  ! <<< Arguments >>>
  character(*), intent(in) :: myTag
  character(len=*), intent(out) :: value
  character(len=*), optional, intent(in) :: default

  ! <<< Local variables >>>
  integer :: iField
  logical :: isDef

  call parserFieldForTag(myTag, iField, isDef)
  if (.not.isDef .AND. present(default)) then   
     value = default
  else if (.not.isDef .AND. .not.present(default)) then
     print *, 'Parser : ' // mytag // ' not defined'
     stop
  else
     read(entries(iField)%value, '(a)') value
  end if

end subroutine parserReadChar

subroutine parserReadIntArray(myTag, value)

  ! <<< Private members >>>
  use ParserImpl, only : parserFieldForTag

  ! <<< Public members >>>
  use Parser, only : entries

  implicit none

  ! <<< Arguments >>>
  character(*), intent(in) :: myTag
  integer, dimension(:), intent(out) :: value

  ! <<< Local variables >>>
  integer :: iField
  logical :: isDef

  ! Read them
  call parserFieldForTag(myTag, iField, isDef)
  if (.not. isDef) then
     print*,'Parser : ' // myTag // ' not defined'
     stop
  end if
  read(entries(ifield)%value, *) value

end subroutine parserReadIntArray

subroutine parserReadFloatArray(myTag, value)

  ! <<< Private members >>>
  use ParserImpl, only : parserFieldForTag

  ! <<< Public members >>>
  use Parser, only : entries

  implicit none

  ! <<< Arguments >>>
  character(*), intent(in) :: mytag
  real(SCALAR_KIND), dimension(:), intent(out) :: value

  ! <<< Local variables >>>
  integer :: iField
  logical :: isDef

  call parserFieldForTag(myTag, iField, isDef)
  if (.not. isDef) then
     print *, 'Parser : ' // mytag // ' not defined'
     stop
  end if
  read(entries(ifield)%value, *) value

end subroutine parserReadFloatArray
                                 
subroutine parserReadFloatArray2D(myTag, value)

  ! <<< Private members >>>
  use ParserImpl, only : parserFieldForTag

  ! <<< Public members >>>
  use Parser, only : entries

  implicit none

  ! <<< Arguments >>>
  character(*), intent(in) :: myTag
  real(SCALAR_KIND), dimension(:,:), intent(out) :: value

  ! <<< Local variables >>>
  integer :: iField
  logical :: isDef

  call parserFieldForTag(myTag, iField, isDef)
  if (.not. isDef) then
     print *, 'Parser : ' // myTag // ' not defined'
     stop
  end if
  read(entries(iField)%value, *) value
    
end subroutine parserReadFloatArray2D
                                
subroutine parserReadCharArray(myTag, value)

  ! <<< Private members >>>
  use ParserImpl, only : parserFieldForTag

  ! <<< Public members >>>
  use Parser, only : entries

  implicit none

  ! <<< Arguments >>>
  character(*), intent(in) :: myTag
  character(*), dimension(:), intent(out) :: value

  ! <<< Local variables >>>
  integer :: iField
  logical :: isDef

  call parserFieldForTag(myTag, iField, isDef)
  if (.not. isDef) then
     print *, 'Parser : ' // mytag // ' not defined'
     stop
  end if
  read(entries(ifield)%value,*) value

end subroutine parserReadCharArray
