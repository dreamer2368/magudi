#include "config.h"

subroutine ZAXPY(comm,ZFilename,A,XFilename,YFilename)
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  SCALAR_TYPE, intent(in) :: A
  character(len = *), intent(in) :: ZFilename, XFilename
  character(len = *), intent(in), optional :: YFilename

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, procRank, numProcs, ierror, mpiFileHandle
  character(len = STRING_LENGTH) :: message
  logical :: XFileExists, YFileExists, success
  integer(kind = MPI_OFFSET_KIND) :: XFileSize, YFileSize, ZFileSize, offset

  integer(kind = MPI_OFFSET_KIND) :: globalSize, bufferSize, sendcnt, indexQuotient
  integer(kind = MPI_OFFSET_KIND), dimension(:), allocatable :: recvcnt, displc
  SCALAR_TYPE, dimension(:), allocatable :: XBuffer, YBuffer, ZBuffer
  SCALAR_TYPE :: dummyValue

  ! Initialize MPI.
  call MPI_Comm_rank(comm, procRank, ierror)
  call MPI_Comm_size(comm, numProcs, ierror)

  call startTiming("File name and size check")

  write(message,'(2A)') "Z filename: ",trim(ZFilename)
  call writeAndFlush(comm, output_unit, message)

  write(message,'(A,E13.6)') "A: ", A
  call writeAndFlush(comm, output_unit, message)

  if (procRank == 0) then
     inquire(file = trim(XFilename), exist = XFileExists)
  end if
  call MPI_Bcast(XFileExists, 1, MPI_LOGICAL, 0, comm, ierror)
  if (.not.XFileExists) then
     write(message, '(A)') "X file does not exists!"
     call gracefulExit(comm, message)
  else
    write(message,'(2A)') "X filename: ",trim(XFilename)
    call writeAndFlush(comm, output_unit, message)

    call MPI_File_open(comm, trim(XFilename),                       &
        MPI_MODE_RDONLY, MPI_INFO_NULL, mpiFileHandle, ierror)
    call MPI_File_get_size(mpiFileHandle, XFileSize, ierror)
    call MPI_File_close(mpiFileHandle, ierror)

    globalSize = XFileSize/8
  end if
  if (PRESENT(YFilename)) then
     if (procRank == 0) then
        inquire(file = trim(YFilename), exist = YFileExists)
     end if
     call MPI_Bcast(YFileExists, 1, MPI_LOGICAL, 0, comm, ierror)
     if (.not.YFileExists) then
        write(message, '(A)') "Y file does not exists!"
        call gracefulExit(comm, message)
     else
       write(message,'(2A)') "Y filename: ",trim(YFilename)
       call writeAndFlush(comm, output_unit, message)

      call MPI_File_open(comm, trim(YFilename),                               &
                          MPI_MODE_RDONLY, MPI_INFO_NULL, mpiFileHandle, ierror)
      call MPI_File_get_size(mpiFileHandle, YFileSize, ierror)
      call MPI_File_close(mpiFileHandle, ierror)
     end if
     if (XFileSize.ne.YFileSize) then
       write(message, '(A)') "X and Y does not match in size!"
       call gracefulExit(comm, message)
     end if
  else
    YFileExists = .false.
    call MPI_Bcast(YFileExists, 1, MPI_LOGICAL, 0, comm, ierror)
    write(message, '(A)') "Y equals to 0."
    call writeAndFlush(comm, output_unit, message)
  end if

  call endTiming("File name and size check")
  call startTiming("Buffer setup")

  indexQuotient = MOD(globalSize,int(numProcs,MPI_OFFSET_KIND))
  bufferSize = globalSize/int(numProcs,MPI_OFFSET_KIND)
  allocate(recvcnt(0:numProcs-1))
  allocate(displc(0:numProcs-1))
  recvcnt(0:indexQuotient-1) = bufferSize+1
  recvcnt(indexQuotient:numProcs-1) = bufferSize
  displc = 0
  do i=1,numProcs-1
    displc(i) = SUM(recvcnt(0:i-1))
  end do
  offset = displc(procRank)*SIZEOF_SCALAR

  if( procRank<indexQuotient ) then
    sendcnt = bufferSize+1
  else
    sendcnt = bufferSize
  end if
  allocate(XBuffer(sendcnt))
  allocate(YBuffer(sendcnt))
  allocate(ZBuffer(sendcnt))
  XBuffer = 0.0_wp
  YBuffer = 0.0_wp
  ZBuffer = 0.0_wp

  call endTiming("Buffer setup")
  call startTiming("Read files")

  call MPI_Barrier(comm, ierror)

  call MPI_File_open(comm, trim(XFilename),                       &
                     MPI_MODE_RDONLY, MPI_INFO_NULL,                        &
                     mpiFileHandle, ierror)
  call MPI_FILE_SET_VIEW(mpiFileHandle,offset,MPI_DOUBLE,MPI_DOUBLE,        &
                         'native',MPI_INFO_NULL,ierror)
  call MPI_File_read_all(mpiFileHandle, XBuffer, int(sendcnt),                   &
                         SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
  call MPI_File_close(mpiFileHandle, ierror)

  call MPI_Barrier(comm, ierror)

  if (YFileExists) then
    call MPI_File_open(comm, trim(YFilename),                       &
                       MPI_MODE_RDONLY, MPI_INFO_NULL,                        &
                       mpiFileHandle, ierror)
    call MPI_FILE_SET_VIEW(mpiFileHandle,offset,MPI_DOUBLE,MPI_DOUBLE,        &
                           'native',MPI_INFO_NULL,ierror)
    call MPI_File_read_all(mpiFileHandle, YBuffer, int(sendcnt),                   &
                           SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
    call MPI_File_close(mpiFileHandle, ierror)
  end if
  call MPI_Barrier(comm, ierror)

  call endTiming("Read files")
  call startTiming("Compute z=Ax+y")

  ZBuffer = A*XBuffer + YBuffer

  call endTiming("Compute z=Ax+y")
  call startTiming("Write file")

  call MPI_Barrier(comm, ierror)

  call MPI_File_open(comm, trim(ZFilename),                       &
                     MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL,      &
                     mpiFileHandle, ierror)
  call MPI_FILE_SET_VIEW(mpiFileHandle,offset,MPI_DOUBLE,MPI_DOUBLE,        &
                         'native',MPI_INFO_NULL,ierror)
  call MPI_File_write_all(mpiFileHandle, ZBuffer, int(sendcnt),                   &
                          SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
  call MPI_File_close(mpiFileHandle, ierror)

  call MPI_Barrier(comm, ierror)

  call endTiming("Write file")

end subroutine
