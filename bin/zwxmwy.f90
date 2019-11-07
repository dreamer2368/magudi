#include "config.h"

program run_zwxmwy

  use ControlSpaceAdvancer, only : zWXMWY
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use ErrorHandler
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers
  use InputHelper, only : getFreeUnit

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: zFilename, WFilename, XFilename, YFilename, normFilename, message
  SCALAR_TYPE :: z
  integer :: stat, fileUnit, procRank, ierror

  ! Initialize MPI.
  call MPI_Init(ierror)
  call initializeErrorHandler()

  if (command_argument_count() < 4) then
     write(message, '(A)') "Require at least 4 arguments!"
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     call MPI_Finalize(ierror)
     stop -1
  end if

  call startTiming("total")

  z = 0.0_wp

  call get_command_argument(1, zFilename)

  call get_command_argument(2, WFilename)

  call get_command_argument(3, XFilename)

  call get_command_argument(4, YFilename)

  if (command_argument_count() == 5) then
     call get_command_argument(5, normFilename)

     z = zWXMWY(MPI_COMM_WORLD,trim(WFilename),trim(XFilename),trim(YFilename),trim(normFilename))
  elseif (command_argument_count() > 5) then
     write(message, '(A)') "Redundant arguments."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     call MPI_Finalize(ierror)
     stop -1
  else
    z = zWXMWY(MPI_COMM_WORLD,trim(WFilename),trim(XFilename),trim(YFilename))
  end if

  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
  if (procRank == 0) then
    open(unit = getFreeUnit(fileUnit), file = trim(zFilename), action='write',          &
         iostat = stat, status = 'replace')
    write(fileUnit, '(1X,SP,' // SCALAR_FORMAT // ')') z
    close(fileUnit)
  end if

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

end program run_zwxmwy
