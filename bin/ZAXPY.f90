#include "config.h"

program RunZAXPY

  use ControlSpaceAdvancer, only : ZAXPY
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use ErrorHandler
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  character(len = STRING_LENGTH) :: ZFilename, Astr, XFilename, YFilename, message
  integer :: ierror
  SCALAR_TYPE :: A

  ! Initialize MPI.
  call MPI_Init(ierror)
  call initializeErrorHandler()

  if (command_argument_count() < 3) then
     write(message, '(A)') "Require at least 3 arguments!"
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     call MPI_Finalize(ierror)
     stop -1
  end if

  call startTiming("total")

  call get_command_argument(1, ZFilename)
  ! write(message,'(2A)') "Z filename: ",trim(ZFilename)
  ! call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

  call get_command_argument(2, Astr)
  read( Astr, * ) A
  ! write(message,'(A,E13.6)') "A: ", A
  ! call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

  call get_command_argument(3, XFilename)
  ! write(message,'(2A)') "X filename: ",trim(XFilename)
  ! call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

  if (command_argument_count() == 4) then
     call get_command_argument(4, YFilename)
     ! write(message,'(2A)') "Y filename: ",trim(YFilename)
     ! call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

     call ZAXPY(MPI_COMM_WORLD,trim(ZFilename),A,trim(XFilename),trim(YFilename))
  elseif (command_argument_count() > 4) then
     write(message, '(A)') "Maximum 4 arguments are needed."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     call MPI_Finalize(ierror)
     stop -1
  else
  !   write(message, '(A)') "Y equals to 0."
  !   call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

    call ZAXPY(MPI_COMM_WORLD,trim(ZFilename),A,trim(XFilename))
  end if

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

end program RunZAXPY
