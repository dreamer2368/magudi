#include "config.h"

#ifdef DEBUG

subroutine assertImpl(condition, conditionString, filename, lineNo)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env, only : error_unit

  implicit none

  ! <<< Arguments >>>
  logical, intent(in) :: condition
  character(len = *), intent(in) :: conditionString, filename
  integer, intent(in) :: lineNo

  ! <<< Local variables >>>
  integer :: ierror

  if (.not. condition) then
     write(error_unit, '(3A,I0.0,3A)') "AssertionError at ", trim(filename), ":", lineNo,    &
          ": ", trim(conditionString), "!"
     call MPI_Abort(MPI_COMM_WORLD, -1, ierror)
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

  ! Use `MPI_COMM_WORLD` if `comm` is identical/congruent/similar to it.
  comm_ = MPI_COMM_WORLD
  call MPI_Comm_compare(comm, MPI_COMM_WORLD, result_, ierror)
  if (result_ == MPI_UNEQUAL) comm_ = comm

  ! If MPI has not yet been initialized, initialize it.
  call MPI_Initialized(flag, ierror)
  if (.not. flag) then
     call MPI_Init(ierror)
  end if

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
