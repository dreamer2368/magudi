#include "config.h"

program options_parser

  use MPI

  use InputHelper
  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: fileUnit, a, proc, ierror
  real(wp) :: b
  logical :: success

  call MPI_Init(ierror)

  success = .true.

  call initializeErrorHandler()

  call MPI_Comm_rank(MPI_COMM_WORLD, proc, ierror)

  if (proc == 0) then
     open(newunit = fileUnit, file = "/tmp/options.txt", action = 'write', status = 'unknown')
     write(fileUnit, '(A)') "a = 41"
     write(fileUnit, '(A)') "c = -9789"
     write(fileUnit, '(A)') "b1 = 2.0"
     write(fileUnit, '(A)') "another option = a string"
     close(fileUnit)
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierror)
  call parseInputFile("/tmp/options.txt")

  call getRequiredOption("a", a)
  success = success .and. (a == 41)
  call getRequiredOption("c", a)
  success = success .and. (a == -9789)
  success = success .and. (getOption("a", 17) == 41)
  success = success .and. (getOption("c", 99) == -9789)
  success = success .and. (getOption("b", 19) == 19)

  call getRequiredOption("b1", b)
  success = success .and. (abs(b - 2.0_wp) <= 0.0_wp)
  success = success .and. (abs(getOption("b1", -1.0_wp) - 2.0_wp) <= 0.0_wp)
  success = success .and. (abs(getOption("b2", 2.5_wp) - 2.5_wp) <= 0.0_wp)

  call cleanupErrorHandler()

  call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierror)
  call MPI_Finalize(ierror)
  if (.not. success) stop -1
  stop 0

end program options_parser
