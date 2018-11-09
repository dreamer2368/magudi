#include "config.h"

program testZAXPY

  use ControlSpaceAdvancer, only : ZAXPY
  use RandomNumber, only : initializeRandomNumberGenerator
  use InputHelper, only : getFreeUnit
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use ErrorHandler
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  character(len = STRING_LENGTH) :: ZFilename, Astr, XFilename, YFilename, message
  integer :: i, stat, fileUnit, procRank, numProcs, ierror, mpiFileHandle, errorCode = 0
  SCALAR_TYPE, dimension(17) :: X, Y, Z, Zoutput
  SCALAR_TYPE :: A

  ! Initialize MPI.
  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)
  call initializeErrorHandler()

  ZFilename = 'testZ.dat'
  YFilename = 'testY.dat'
  XFilename = 'testX.dat'

  if( procRank.eq.0 ) then
    call RANDOM_NUMBER(A)
    call RANDOM_NUMBER(X)
    call RANDOM_NUMBER(Y)
    Z = A*X + Y

    open(unit = getFreeUnit(fileUnit), file = trim(XFilename), action='write',          &
          form='unformatted', access='stream', iostat = stat, status = 'replace')
    write(fileUnit) X
    close(fileUnit)
    open(unit = getFreeUnit(fileUnit), file = trim(YFilename), action='write',          &
          form='unformatted', access='stream', iostat = stat, status = 'replace')
    write(fileUnit) Y
    close(fileUnit)
  end if
  call MPI_Bcast(A, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)

  call ZAXPY(MPI_COMM_WORLD,trim(ZFilename),A,trim(XFilename),trim(YFilename))

  if( procRank.eq.0 ) then
    open(unit = getFreeUnit(fileUnit), file = trim(ZFilename), action='read',          &
          form='unformatted', access='stream', iostat = stat, status = 'old')
    read(fileUnit) Zoutput
    close(fileUnit)

    if( MAXVAL(ABS(Zoutput-Z)).ge.1.0E-15 ) then
      errorCode = -1
      print *, 'error: ', MAXVAL(ABS(Zoutput-Z))
      print *, 'error array'
      print *, Z-Zoutput
    end if
  end if

  call MPI_Bcast(errorCode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

  if( errorCode.ne.0 ) then
    stop -1
  else
    stop 0
  end if

end program testZAXPY
