#include "config.h"

module ErrorHandler

  implicit none
  public

#ifdef DEBUG

  interface

     !> Implements the subroutine that checks an assertion called by the C-style preprocessor
     !> macro `assert`. The current behavior is to print the details of the assertion that
     !> failed and call `MPI_Abort` from the communicator `MPI_COMM_WORLD`.

     subroutine assertImpl(condition, conditionString, filename, lineNo)

       logical, intent(in) :: condition
       character(len = *), intent(in) :: conditionString, filename
       integer, intent(in) :: lineNo

     end subroutine assertImpl

  end interface

#endif

  interface

     subroutine writeAndFlush(comm, unit, str, advance)

       !> Attempts to write the string `str` to the Fortran file unit `unit` from the master
       !> process of the MPI communicator `comm` so that it arrives at the stream in the order
       !> intended. This subroutine must be called collectively by all processes in `comm`.

       integer, intent(in) :: comm, unit
       character(len = *), intent(in) :: str

       character(len = *), intent(in), optional :: advance

     end subroutine writeAndFlush

  end interface

  interface

     subroutine gracefulExit(comm, errorMessage)

       !> Attempts to gracefully exit from an MPI application after an error has
       !> occurred. Aborts execution from all processes in `comm` if its relationship with
       !> `MPI_COMM_WORLD` is `MPI_UNEQUAL`. This subroutine must be called collectively by
       !> all processes in `comm`.

       integer, intent(in) :: comm
       character(len = *), intent(in) :: errorMessage

     end subroutine gracefulExit

  end interface

  interface

     subroutine issueWarning(comm, warningMessage)

       !> Issues a warning message and returns control to an MPI application. This subroutine
       !> must be called collectively by all processes in the MPI communicator `comm`.

       integer, intent(in) :: comm
       character(len = *), intent(in) :: warningMessage

     end subroutine issueWarning

  end interface

end module ErrorHandler
