#include "config.h"

module InputHelperImpl

  implicit none
  public

  type, private :: t_DictElement
     character(len = STRING_LENGTH) :: key, val
  end type t_DictElement

  type(t_DictElement), allocatable, public :: dict(:)

contains

  function split(str, leftSubstring, rightSubstring, separator) result(extraSeparatorCount)

    ! <<< Arguments >>>
    character(len = *), intent(in) :: str
    character(len = *), intent(out) :: leftSubstring, rightSubstring
    character(len = 1), intent(in) :: separator

    ! <<< Result >>>
    integer :: extraSeparatorCount

    ! <<< Local variables >>>
    integer :: i, j

    i = 1

    do j = 1, len_trim(str)
       if (str(j:j) == separator) exit
       leftSubstring(i:i) = str(j:j)
       i = i + 1
    end do

    leftSubstring(i:) = " "
    rightSubstring = str(j+1:len_trim(str))
    call trimAll(leftSubstring)
    call trimAll(rightSubstring)

    extraSeparatorCount = scan(rightSubstring, separator)

  end function split

  subroutine trimAll(str)

    ! <<< Arguments >>>
    character(len = *), intent(inout) :: str

    ! <<< Local variables >>>
    integer :: i

    do i = 1, len(str)
       if (ichar(str(i:i)) > 32) exit
       if (ichar(str(i:i)) <= 32) str(i:i) = " "
    end do

    do i = len(str), 1, -1
       if (ichar(str(i:i)) > 32) exit
       if (ichar(str(i:i)) <= 32) str(i:i) = " "
    end do

    str = trim(adjustl(str))

  end subroutine trimAll

  subroutine sort()

    ! <<< Local variables >>>
    integer :: i, j
    type(t_DictElement) :: temp

    do i = 2, size(dict)
       j = i - 1
       temp = dict(i)
       do while (dict(j)%key > temp%key)
          dict(j+1) = dict(j)
          j = j - 1
          if (j < 1) exit
       end do
       dict(j+1) = temp
    end do

  end subroutine sort

  subroutine find(key, index)

    ! <<< Arguments >>>
    character(len = *), intent(in) :: key
    integer, intent(out) :: index

    ! <<< Local variables >>>
    integer :: iStart, iEnd, iMid

    index = -1

    if (allocated(dict)) then

       iStart = 1
       iEnd = size(dict)

       do while (iEnd >= iStart) !... `dict` is sorted; use binary search.
          iMid = (iStart + iEnd) / 2
          if (dict(iMid)%key == trim(key)) then
             index = iMid
             return
          end if
          if (dict(iMid)%key > trim(key)) then
             iEnd = iMid - 1
          else
             iStart = iMid + 1
          end if
       end do

    end if

  end subroutine find

end module InputHelperImpl

subroutine getInputName(filename)

  ! <<< Internal modules >>>
  use ErrorHandler, only : writeAndFlush

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  implicit none

  ! <<< Arguments >>>
  character(len = STRING_LENGTH), intent(out) :: filename

  ! <<< Local variables >>>
  integer :: ierror
  character(len = STRING_LENGTH) :: message

  if (command_argument_count() > 1) then
     write(message, '(A)') "Usage: magudi [INPUT]"
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)') "High-performance Fortran-based adjoint optimization tool."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)')                                                                   &
          "Maximum of 1 INPUT file allowed."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     call cleanupErrorHandler()
     call MPI_Finalize(ierror)
     stop -1
  end if

  ! Get the input file name.
  if (command_argument_count() == 1) then
     call get_command_argument(1, filename)
     if (filename(1:1) == '-' .or. len_trim(filename) == 0) then
        write(message, '(A)') "No input file name was detected, using 'input'."
        call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
        filename = 'input'
     end if
  else
     write(message, '(A)') "No input file name was detected, using 'input'."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     filename = 'input'
  end if

end subroutine getInputName

function getFreeUnit(fileUnit) result(freeUnit)

  ! <<< Arguments >>>
  integer, intent(out), optional :: fileUnit

  ! <<< Result >>>
  integer :: freeUnit

  ! <<< Local variables >>>
  logical :: isOpened

  do freeUnit = 10, 10000
     inquire(unit = freeUnit, opened = isOpened)
     if (.not. isOpened) exit
  end do

  if (present(fileUnit)) fileUnit = freeUnit

end function getFreeUnit

subroutine parseInputFile(filename, commentMarker, separator)

  ! <<< External modules >>>
  use MPI

  ! <<< Private members >>>
  use InputHelperImpl, only : dict, split, sort

  ! <<< Public members >>>
  use InputHelper, only : stripComments, getFreeUnit

  ! <<< Internal modules >>>
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: filename
  character(len = 1), intent(in), optional :: commentMarker, separator

  ! <<< Local variables >>>
  character(len = 1) :: commentMarker_, separator_
  integer :: i, fileUnit, procRank, dictSize, lineNo, istat, ierror
  character(len = STRING_LENGTH) :: line, message

  assert(len_trim(filename) > 0)

  ! Assign default values to optional arguments.
  commentMarker_ = '#'
  if (present(commentMarker)) commentMarker_ = commentMarker
  separator_ = '='
  if (present(separator)) separator_ = separator

  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  ! Check if file exists.
  if (procRank == 0) then
     open(unit = getFreeUnit(fileUnit), file = trim(filename), action = 'read',              &
          status = 'old', iostat = istat)
  end if
  call MPI_Bcast(istat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  if (istat /= 0) then
     write(message, "(2A)") trim(filename), ": File not found or permission denied!"
     call gracefulExit(MPI_COMM_WORLD, message)
  end if

  ! Only the root process reads the input file.
  if (procRank == 0) then
     dictSize = 0
     do !... read once to find input dictionary size.
        read(fileUnit, '(A)', iostat = istat) line
        if (istat < 0) exit
        call stripComments(line, commentMarker_) !... skip comments.
        if (len_trim(line) == 0) cycle !... skip empty lines.
        dictSize = dictSize + 1
     end do
     close(fileUnit)
  end if

  ! Broadcast the input dictionary size to all processes.
  call MPI_Bcast(dictSize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  if (dictSize == 0) then
     write(message, "(2A)") trim(filename), ": File is empty or does not contain any input!"
     call gracefulExit(MPI_COMM_WORLD, message)
  end if

  ! Allocate memory to hold the dictionary.
  SAFE_DEALLOCATE(dict)
  allocate(dict(dictSize), stat = istat)
  if (istat /= 0) then
     write(message, "(A)") "Insufficient memory: Could not allocate storage for input!"
     call gracefulExit(MPI_COMM_WORLD, message)
  end if

  ! Again, only the root process reads the input file.
  if (procRank == 0) then
     i = 0 ; lineNo = 0 ; istat = 0
     open(unit = getFreeUnit(fileUnit), file = trim(filename),                               &
          action = 'read', status = 'old')
     do !... read again to fill input dictionary.

        read(fileUnit, '(A)', iostat = istat) line
        if (istat < 0) then
           istat = 0
           exit
        end if
        lineNo = lineNo + 1 !... need line number for reporting errors.
        call stripComments(line, commentMarker_)
        if (len_trim(line) == 0) cycle
        i = i + 1

        ! Parse in 'key <separator> value' format.
        istat = split(line, dict(i)%key, dict(i)%val, separator_)

        if (istat /= 0) then
           write(message, "(2A,I0.0,A)") trim(filename), ": Failed to parse input on line ", &
                lineNo, "!"
           exit
        end if

        if (len_trim(dict(i)%key) == 0) then
           istat = 1
           write(message, "(2A,I0.0,A)") trim(filename), ": Empty parameter key on line ",   &
                lineNo, "!"
           exit
        end if

     end do
     close(fileUnit)
  end if

  ! Broadcast istat to collectively return on error.
  call MPI_Bcast(istat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

  ! Check if an error occurred.
  if (istat /= 0) then
     call MPI_Bcast(message, len(message), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)
     call gracefulExit(MPI_COMM_WORLD, message)
  end if

  ! Sort the dictionary and broadcast it to all processes.
  call sort()
  do i = 1, dictSize
     call MPI_Bcast(dict(i)%key, STRING_LENGTH, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)
     call MPI_Bcast(dict(i)%val, STRING_LENGTH, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)
  end do

end subroutine parseInputFile

subroutine stripComments(str, commentMarker)

  ! <<< Private members >>>
  use InputHelperImpl, only : trimAll

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(inout) :: str
  character(len = 1), intent(in) :: commentMarker

  ! <<< Local variables >>>
  integer :: i

  call trimAll(str)

  do i = 1, len(str)
     if (str(i:i) == commentMarker) then
        str(i:) = " "
        exit
     end if
  end do

end subroutine stripComments

function getOptionInteger_(key, defaultValue) result(val)

  ! <<< Private members >>>
  use InputHelperImpl, only : dict, find

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  integer, intent(in) :: defaultValue

  ! <<< Result >>>
  integer :: val

  ! <<< Local variables >>>
  integer :: index, stat

  call find(key, index)

  if (index == -1) then
     val = defaultValue
     return
  end if

  read(dict(index)%val, *, iostat = stat) val
  if (stat /= 0) val = defaultValue

end function getOptionInteger_

function getOptionLogical_(key, defaultValue) result(val)

  ! <<< Private members >>>
  use InputHelperImpl, only : dict, find

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  logical, intent(in) :: defaultValue

  ! <<< Result >>>
  logical :: val

  ! <<< Local variables >>>
  integer :: index

  call find(key, index)

  if (index == -1) then
     val = defaultValue
     return
  end if

  if (trim(dict(index)%val) == "true" .or. trim(dict(index)%val) == "TRUE" .or.              &
       trim(dict(index)%val) == "True") then
     val = .true.
  else if (trim(dict(index)%val) == "false" .or. trim(dict(index)%val) == "FALSE" .or.       &
       trim(dict(index)%val) == "False") then
     val = .false.
  else
     val = defaultValue
  end if

end function getOptionLogical_

#ifdef SCALAR_IS_COMPLEX
function getOptionReal_(key, defaultValue) result(val)

  ! <<< Private members >>>
  use InputHelperImpl, only : dict, find

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  real(SCALAR_KIND), intent(in) :: defaultValue

  ! <<< Result >>>
  real(SCALAR_KIND) :: val

  ! <<< Local variables >>>
  integer :: index, stat

  call find(key, index)

  if (index == -1) then
     val = defaultValue
     return
  end if

  read(dict(index)%val, *, iostat = stat) val
  if (stat /= 0) val = defaultValue

end function getOptionReal_
#endif

function getOptionScalar_(key, defaultValue) result(val)

  ! <<< External modules >>>

  ! <<< Private members >>>
  use InputHelperImpl, only : dict, find

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  SCALAR_TYPE, intent(in) :: defaultValue

  ! <<< Result >>>
  SCALAR_TYPE :: val

  ! <<< Local variables >>>
#ifdef SCALAR_IS_COMPLEX
  real(SCALAR_KIND) :: val_
#endif
  integer :: index, stat

  call find(key, index)

  if (index == -1) then
     val = defaultValue
     return
  end if

  read(dict(index)%val, *, iostat = stat) val
  if (stat /= 0) then

#ifdef SCALAR_IS_COMPLEX
     read(dict(index)%val, *, iostat = stat) val_
     if (stat /= 0) then
        val = defaultValue
     else
        val = val_
     end if
#else
     val = defaultValue
#endif

  end if

end function getOptionScalar_

function getOptionString_(key, defaultValue) result(val)

  ! <<< Private members >>>
  use InputHelperImpl, only : dict, find

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  character(len = *), intent(in) :: defaultValue

  ! <<< Result >>>
  character(len = STRING_LENGTH) :: val

  ! <<< Local variables >>>
  integer :: index

  call find(key, index)

  if (index == -1) then
     if (len_trim(defaultValue) == 0) then
        val = ""
     else
        read(defaultValue, '(A)') val
     end if
     return
  end if

  read(dict(index)%val, *) val

end function getOptionString_

subroutine getRequiredOptionInteger_(key, val, comm)

  ! <<< External modules >>>
  use MPI

  ! <<< Private members >>>
  use InputHelperImpl, only : dict, find

  ! <<< Internal modules >>>
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  integer, intent(out) :: val
  integer, intent(in), optional :: comm

  ! <<< Local variables >>>
  integer :: comm_, index, stat
  character(len = STRING_LENGTH) :: message

  comm_ = MPI_COMM_WORLD
  if (present(comm)) comm_ = comm

  call find(key, index)
  if (index /= -1) read(dict(index)%val, *, iostat = stat) val

  if (index == -1 .or. stat /= 0) then
     write(message, "(3A)") "Required parameter '", trim(key), "' not found or invalid!"
     call gracefulExit(comm_, message)
  end if

end subroutine getRequiredOptionInteger_

subroutine getRequiredOptionLogical_(key, val, comm)

  ! <<< External modules >>>
  use MPI

  ! <<< Private members >>>
  use InputHelperImpl, only : dict, find

  ! <<< Internal modules >>>
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  logical, intent(out) :: val
  integer, intent(in), optional :: comm

  ! <<< Local variables >>>
  integer :: comm_, index, stat
  character(len = STRING_LENGTH) :: message

  comm_ = MPI_COMM_WORLD
  if (present(comm)) comm_ = comm

  call find(key, index)
  if (index /= -1) then
     if (trim(dict(index)%val) == "true" .or. trim(dict(index)%val) == "TRUE" .or.           &
          trim(dict(index)%val) == "True") then
        val = .true.
     else if (trim(dict(index)%val) == "false" .or. trim(dict(index)%val) == "FALSE" .or.    &
          trim(dict(index)%val) == "False") then
        val = .false.
     else
        stat = -1
     end if
  end if

  if (index == -1 .or. stat /= 0) then
     write(message, "(3A)") "Required parameter '", trim(key), "' not found or invalid!"
     call gracefulExit(comm, message)
  end if

end subroutine getRequiredOptionLogical_

#ifdef SCALAR_IS_COMPLEX
subroutine getRequiredOptionReal_(key, val, comm)

  ! <<< External modules >>>
  use MPI

  ! <<< Private members >>>
  use InputHelperImpl, only : dict, find

  ! <<< Internal modules >>>
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  real(SCALAR_KIND), intent(out) :: val
  integer, intent(in), optional :: comm

  ! <<< Local variables >>>
  integer :: comm_, index, stat
  character(len = STRING_LENGTH) :: message

  comm_ = MPI_COMM_WORLD
  if (present(comm)) comm_ = comm

  call find(key, index)
  if (index /= -1) read(dict(index)%val, *, iostat = stat) val

  if (index == -1 .or. stat /= 0) then
     write(message, "(3A)") "Required parameter '", trim(key), "' not found or invalid!"
     call gracefulExit(comm_, message)
  end if

end subroutine getRequiredOptionReal_
#endif

subroutine getRequiredOptionScalar_(key, val, comm)

  ! <<< External modules >>>
  use MPI

  ! <<< Private members >>>
  use InputHelperImpl, only : dict, find

  ! <<< Internal modules >>>
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  SCALAR_TYPE, intent(out) :: val
  integer, intent(in), optional :: comm

  ! <<< Local variables >>>
#ifdef SCALAR_IS_COMPLEX
  real(SCALAR_KIND) :: val_
#endif
  integer :: comm_, index, stat
  character(len = STRING_LENGTH) :: message

  comm_ = MPI_COMM_WORLD
  if (present(comm)) comm_ = comm

  call find(key, index)
  if (index /= -1) then

     read(dict(index)%val, *, iostat = stat) val
#ifdef SCALAR_IS_COMPLEX
     if (stat /= 0) then
        read(dict(index)%val, *, iostat = stat) val_
        if (stat == 0) val = val_
     end if
#endif

  end if

  if (index == -1 .or. stat /= 0) then
     write(message, "(3A)") "Required parameter '", trim(key), "' not found or invalid!"
     call gracefulExit(comm_, message)
  end if

end subroutine getRequiredOptionScalar_

subroutine getRequiredOptionString_(key, val, comm)

  ! <<< External modules >>>
  use MPI

  ! <<< Private members >>>
  use InputHelperImpl, only : dict, find

  ! <<< Internal modules >>>
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: key
  character(len = STRING_LENGTH), intent(out) :: val
  integer, intent(in), optional :: comm

  ! <<< Local variables >>>
  integer :: comm_, index
  character(len = STRING_LENGTH) :: message

  comm_ = MPI_COMM_WORLD
  if (present(comm)) comm_ = comm

  call find(key, index)
  if (index /= -1) then
     read(dict(index)%val, *) val
  else
     write(message, "(3A)") "Required parameter '", trim(key), "' not found or invalid!"
     call gracefulExit(comm_, message)
  end if

end subroutine getRequiredOptionString_
