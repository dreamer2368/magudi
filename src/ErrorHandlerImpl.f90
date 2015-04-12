#include "config.h"

module ErrorHandlerImpl

  use MPI, only : MPI_WIN_NULL

  implicit none
  public

  integer, public :: mpiWindowBase = 0, mpiWindow = MPI_WIN_NULL

end module ErrorHandlerImpl

subroutine initializeErrorHandler()

  ! <<< External modules >>>
  use MPI

  ! <<< Private members >>>
  use ErrorHandlerImpl, only : mpiWindowBase, mpiWindow

  ! <<< Public members >>>
  use ErrorHandler, only : cleanupErrorHandler

  implicit none

  ! <<< Local variables >>>
  logical :: flag
  integer :: procRank, ierror

  ! If MPI has not yet been initialized, initialize it.
  call MPI_Initialized(flag, ierror)
  if (.not. flag) then
     call MPI_Init(ierror)
  end if

  call cleanupErrorHandler()

  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  if (procRank == 0) then
     call MPI_Win_create(mpiWindowBase, int(4, MPI_ADDRESS_KIND), 4, MPI_INFO_NULL,          &
          MPI_COMM_WORLD, mpiWindow, ierror)
  else
     call MPI_Win_create(0, int(0, MPI_ADDRESS_KIND), 1, MPI_INFO_NULL,                      &
          MPI_COMM_WORLD, mpiWindow, ierror)
  end if

end subroutine initializeErrorHandler

#ifdef DEBUG

subroutine assertImpl(condition, conditionString, filename, lineNo)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env, only : error_unit

  ! <<< Private members >>>
  use ErrorHandlerImpl, only : mpiWindow

  implicit none

  ! <<< Arguments >>>
  logical, intent(in) :: condition
  character(len = *), intent(in) :: conditionString, filename
  integer, intent(in) :: lineNo

  ! <<< Local variables >>>
  integer :: i, j, one = 1, flag = 0, procRank, ierror
  character(len = len_trim(conditionString)) :: str1
  character(len = len_trim(filename)) :: str2

  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  str1(1:1) = conditionString(1:1)
  j = 2
  do i = 2, len_trim(conditionString)
     if (conditionString(i:i) /= ' ' .or. conditionString(i-1:i-1) /= ' ') then
        str1(j:j) = conditionString(i:i)
        j = j + 1
     end if
  end do
  str1(j:) = ' '

  str2 = filename
  i = index(str2, '/', .true.)
  if (i /= 0) str2(:i) = ' '
  str2 = adjustl(str2)

  if (.not. condition) then
#ifdef DEBUG
     call MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, mpiWindow, ierror)

     ! Disabling this due to poor portability (enable when MPI-3 is widely supported):

     !# call MPI_Fetch_and_op(one, flag, MPI_INTEGER, 0, int(0, MPI_ADDRESS_KIND),           &
     !#      MPI_SUM, mpiWindow, ierror)

     call MPI_Win_unlock(0, mpiWindow, ierror)
#endif
     if (flag == 0) then
        write(error_unit, '(3A,I0.0,3A)') "AssertionError at ",                              &
             trim(str2), ":", lineNo, ": ", trim(str1), "!"
        call MPI_Abort(MPI_COMM_WORLD, -1, ierror)
        call abort
     end if
  end if

end subroutine assertImpl

#endif

subroutine writeAndFlush(comm, unit, str, advance)

  ! <<< External modules >>>
  use MPI

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm, unit
  character(len = *), intent(in) :: str
  character(len = *), intent(in), optional :: advance

  ! <<< Local variables >>>
  integer :: procRank, ierror

  call MPI_Comm_rank(comm, procRank, ierror)

  if (procRank == 0) then

     ! Write `str` to `unit` from master process.
     if (present(advance)) then
        write(unit, '(A)', advance = trim(advance)) trim(str)
     else
        write(unit, '(A)') trim(str)
     end if

     ! Flush `unit` from master process.
     flush(unit)

  end if

  ! Ensure that no writes are executed from other processes in `comm` before master process
  ! has flushed `unit`.
  call MPI_Barrier(comm, ierror)

end subroutine writeAndFlush

subroutine gracefulExit(comm, errorMessage)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env, only : error_unit

  ! <<< Public members >>>
  use ErrorHandler, only : writeAndFlush

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: errorMessage

  ! <<< Local variables >>>
  integer :: comm_, result_, ierror
  logical :: flag

  ! If MPI has not yet been initialized, initialize it.
  call MPI_Initialized(flag, ierror)
  if (.not. flag) then
     call MPI_Init(ierror)
  end if

  ! Use `MPI_COMM_WORLD` if `comm` is identical/congruent/similar to it.
  comm_ = MPI_COMM_WORLD
  call MPI_Comm_compare(comm, MPI_COMM_WORLD, result_, ierror)
  if (result_ == MPI_UNEQUAL) comm_ = comm

  ! Write the error message.
  call writeAndFlush(comm_, error_unit, char(27) // "[1;31mERROR" // char(27) //             &
       "[0m: " // trim(errorMessage))

  ! Try to exit gracefully.
  if (comm_ == MPI_COMM_WORLD) then
     call MPI_Finalize(ierror)
  else
     call MPI_Abort(comm_, -1, ierror)
  end if
  stop -1

end subroutine gracefulExit

subroutine issueWarning(comm, warningMessage)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env, only : error_unit

  ! <<< Public members >>>
  use ErrorHandler, only : writeAndFlush

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: warningMessage

  ! <<< Local variables >>>
  integer :: comm_, result_, ierror
  logical :: flag

  ! Use `MPI_COMM_WORLD` if `comm` is identical/congruent/similar to it.
  comm_ = MPI_COMM_WORLD
  call MPI_Comm_compare(comm, MPI_COMM_WORLD, result_, ierror)
  if (result_ == MPI_UNEQUAL) comm_ = comm

  ! If MPI has not yet been initialized, initialize it.
  call MPI_Initialized(flag, ierror)
  if (.not. flag) then
     call MPI_Init(ierror)
  end if

  ! Write the warning message.
  call writeAndFlush(comm_, error_unit, char(27) // "[1;33mWARNING" // char(27) //           &
       "[0m: " // trim(warningMessage))

end subroutine issueWarning

subroutine cleanupErrorHandler()

  ! <<< External modules >>>
  use MPI

  ! <<< Private members >>>
  use ErrorHandlerImpl, only : mpiWindowBase, mpiWindow

  implicit none

  ! <<< Local variables >>>
  integer :: ierror

  mpiWindowBase = 0

  if (mpiWindow /= MPI_WIN_NULL) then
     call MPI_Win_free(mpiWindow, ierror)
  end if
  mpiWindow = MPI_WIN_NULL

end subroutine cleanupErrorHandler
