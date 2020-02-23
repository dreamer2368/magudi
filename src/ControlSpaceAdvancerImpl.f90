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

    globalSize = XFileSize/SIZEOF_SCALAR
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
  call startTiming("Compute Z=aX+Y")

  ZBuffer = A*XBuffer + YBuffer

  call endTiming("Compute Z=aX+Y")
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

function zWXMWY(comm,WFilename,XFilename,YFilename, normFilename) result(z)
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: WFilename, XFilename, YFilename
  character(len = *), intent(in), optional :: normFilename

  ! <<< Result >>>
  SCALAR_TYPE :: z

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, procRank, numProcs, ierror, mpiFileHandle
  character(len = STRING_LENGTH) :: message
  logical :: WFileExists, XFileExists, YFileExists, normFileExists, success
  integer(kind = MPI_OFFSET_KIND) :: WFileSize, XFileSize, YFileSize, normFileSize, offset

  integer(kind = MPI_OFFSET_KIND) :: globalSize, bufferSize, sendcnt, indexQuotient
  integer(kind = MPI_OFFSET_KIND), dimension(:), allocatable :: recvcnt, displc
  SCALAR_TYPE, dimension(:), allocatable :: XBuffer, YBuffer, WBuffer, normBuffer
  SCALAR_TYPE :: dummyValue

  ! Initialize MPI.
  call MPI_Comm_rank(comm, procRank, ierror)
  call MPI_Comm_size(comm, numProcs, ierror)

  call startTiming("File name and size check")

  if (procRank == 0) then
     inquire(file = trim(WFilename), exist = WFileExists)
     inquire(file = trim(XFilename), exist = XFileExists)
     inquire(file = trim(YFilename), exist = YFileExists)
  end if
  call MPI_Bcast(WFileExists, 1, MPI_LOGICAL, 0, comm, ierror)
  call MPI_Bcast(XFileExists, 1, MPI_LOGICAL, 0, comm, ierror)
  call MPI_Bcast(YFileExists, 1, MPI_LOGICAL, 0, comm, ierror)

  if (.not.WFileExists) then
     write(message, '(A)') "W file does not exists!"
     call gracefulExit(comm, message)
  else
    write(message,'(2A)') "W filename: ",trim(WFilename)
    call writeAndFlush(comm, output_unit, message)

    call MPI_File_open(comm, trim(WFilename),                       &
        MPI_MODE_RDONLY, MPI_INFO_NULL, mpiFileHandle, ierror)
    call MPI_File_get_size(mpiFileHandle, WFileSize, ierror)
    call MPI_File_close(mpiFileHandle, ierror)
  end if
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
  end if
  if (.not.YFileExists) then
     write(message, '(A)') "Y file does not exists!"
     call gracefulExit(comm, message)
  else
    write(message,'(2A)') "Y filename: ",trim(YFilename)
    call writeAndFlush(comm, output_unit, message)

    call MPI_File_open(comm, trim(YFilename),                       &
        MPI_MODE_RDONLY, MPI_INFO_NULL, mpiFileHandle, ierror)
    call MPI_File_get_size(mpiFileHandle, YFileSize, ierror)
    call MPI_File_close(mpiFileHandle, ierror)
  end if

  if (PRESENT(normFilename)) then
     if (procRank == 0) then
        inquire(file = trim(normFilename), exist = normFileExists)
     end if
     call MPI_Bcast(normFileExists, 1, MPI_LOGICAL, 0, comm, ierror)
     if (.not.normFileExists) then
        write(message, '(A)') "norm file does not exists!"
        call gracefulExit(comm, message)
     else
       write(message,'(2A)') "norm filename: ",trim(normFilename)
       call writeAndFlush(comm, output_unit, message)

      call MPI_File_open(comm, trim(normFilename),                               &
                          MPI_MODE_RDONLY, MPI_INFO_NULL, mpiFileHandle, ierror)
      call MPI_File_get_size(mpiFileHandle, normFileSize, ierror)
      call MPI_File_close(mpiFileHandle, ierror)
     end if
  else
    normFileExists = .false.
    call MPI_Bcast(normFileExists, 1, MPI_LOGICAL, 0, comm, ierror)
    write(message, '(A)') "norm equals to 1 uniformly."
    call writeAndFlush(comm, output_unit, message)
  end if

  if (XFileSize.ne.WFileSize                                                    &
  .or. YFileSize.ne.WFileSize                                                   &
  .or. (normFileExists .and. normFileSize.ne.WFileSize)) then
    write(message, '(A)') "W, X, Y (and norm) do not match in size!"
    call gracefulExit(comm, message)
  end if
  globalSize = WFileSize/SIZEOF_SCALAR

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
  allocate(WBuffer(sendcnt))
  allocate(normBuffer(sendcnt))
  XBuffer = 0.0_wp
  YBuffer = 0.0_wp
  WBuffer = 0.0_wp
  normBuffer = 1.0_wp

  call endTiming("Buffer setup")
  call startTiming("Read files")

  call MPI_Barrier(comm, ierror)

  call MPI_File_open(comm, trim(WFilename),                       &
                     MPI_MODE_RDONLY, MPI_INFO_NULL,                        &
                     mpiFileHandle, ierror)
  call MPI_FILE_SET_VIEW(mpiFileHandle,offset,MPI_DOUBLE,MPI_DOUBLE,        &
                         'native',MPI_INFO_NULL,ierror)
  call MPI_File_read_all(mpiFileHandle, WBuffer, int(sendcnt),                   &
                         SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
  call MPI_File_close(mpiFileHandle, ierror)

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

  call MPI_File_open(comm, trim(YFilename),                       &
                     MPI_MODE_RDONLY, MPI_INFO_NULL,                        &
                     mpiFileHandle, ierror)
  call MPI_FILE_SET_VIEW(mpiFileHandle,offset,MPI_DOUBLE,MPI_DOUBLE,        &
                         'native',MPI_INFO_NULL,ierror)
  call MPI_File_read_all(mpiFileHandle, YBuffer, int(sendcnt),                   &
                         SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
  call MPI_File_close(mpiFileHandle, ierror)

  call MPI_Barrier(comm, ierror)

  if (normFileExists) then
    call MPI_File_open(comm, trim(normFilename),                       &
                       MPI_MODE_RDONLY, MPI_INFO_NULL,                        &
                       mpiFileHandle, ierror)
    call MPI_FILE_SET_VIEW(mpiFileHandle,offset,MPI_DOUBLE,MPI_DOUBLE,        &
                           'native',MPI_INFO_NULL,ierror)
    call MPI_File_read_all(mpiFileHandle, normBuffer, int(sendcnt),                   &
                           SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
    call MPI_File_close(mpiFileHandle, ierror)
  end if
  call MPI_Barrier(comm, ierror)

  call endTiming("Read files")
  call startTiming("Compute z=W^T*norm*(X-Y)")

  z = 0.0_wp
  z = SUM( WBuffer*normBuffer*(XBuffer - YBuffer) )

  call MPI_Allreduce(MPI_IN_PLACE, z, 1, SCALAR_TYPE_MPI, MPI_SUM, comm, ierror)
  call MPI_Bcast(z, 1, SCALAR_TYPE_MPI, 0, comm, ierror)

  call endTiming("Compute z=W^T*norm*(X-Y)")

end function

function zXdotY(comm,XFilename,YFilename,normFilename) result(z)
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: XFilename, YFilename
  character(len = *), intent(in), optional :: normFilename

  ! <<< Result >>>
  SCALAR_TYPE :: z

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, procRank, numProcs, ierror, mpiFileHandle
  character(len = STRING_LENGTH) :: message
  logical :: XFileExists, YFileExists, normFileExists, success
  integer(kind = MPI_OFFSET_KIND) :: XFileSize, YFileSize, normFileSize, offset

  integer(kind = MPI_OFFSET_KIND) :: globalSize, bufferSize, sendcnt, indexQuotient
  integer(kind = MPI_OFFSET_KIND) :: globalBufferSize, globalOffset
  integer(kind = MPI_OFFSET_KIND), dimension(:), allocatable :: recvcnt, displc
  SCALAR_TYPE, dimension(:), allocatable :: XBuffer, YBuffer, normBuffer
  SCALAR_TYPE :: dummyValue

  ! Initialize MPI.
  call MPI_Comm_rank(comm, procRank, ierror)
  call MPI_Comm_size(comm, numProcs, ierror)

  call startTiming("File name and size check")

  if (procRank == 0) then
     inquire(file = trim(XFilename), exist = XFileExists)
     inquire(file = trim(YFilename), exist = YFileExists)
  end if
  call MPI_Bcast(XFileExists, 1, MPI_LOGICAL, 0, comm, ierror)
  call MPI_Bcast(YFileExists, 1, MPI_LOGICAL, 0, comm, ierror)

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
  end if
  if (.not.YFileExists) then
     write(message, '(A)') "Y file does not exists!"
     call gracefulExit(comm, message)
  else
    write(message,'(2A)') "Y filename: ",trim(YFilename)
    call writeAndFlush(comm, output_unit, message)

    call MPI_File_open(comm, trim(YFilename),                       &
        MPI_MODE_RDONLY, MPI_INFO_NULL, mpiFileHandle, ierror)
    call MPI_File_get_size(mpiFileHandle, YFileSize, ierror)
    call MPI_File_close(mpiFileHandle, ierror)
  end if

  if (PRESENT(normFilename)) then
     if (procRank == 0) then
        inquire(file = trim(normFilename), exist = normFileExists)
     end if
     call MPI_Bcast(normFileExists, 1, MPI_LOGICAL, 0, comm, ierror)
     if (.not.normFileExists) then
        write(message, '(A)') "norm file does not exists!"
        call gracefulExit(comm, message)
     else
       write(message,'(2A)') "norm filename: ",trim(normFilename)
       call writeAndFlush(comm, output_unit, message)

      call MPI_File_open(comm, trim(normFilename),                               &
                          MPI_MODE_RDONLY, MPI_INFO_NULL, mpiFileHandle, ierror)
      call MPI_File_get_size(mpiFileHandle, normFileSize, ierror)
      call MPI_File_close(mpiFileHandle, ierror)
     end if
  else
    normFileExists = .false.
    call MPI_Bcast(normFileExists, 1, MPI_LOGICAL, 0, comm, ierror)
    write(message, '(A)') "norm equals to 1 uniformly."
    call writeAndFlush(comm, output_unit, message)
  end if

  if (XFileSize.ne.YFileSize                                                    &
  .or. (normFileExists .and. normFileSize.ne.XFileSize)) then
    write(message, '(A)') "X, Y (and norm) do not match in size!"
    call gracefulExit(comm, message)
  end if
  globalSize = XFileSize/SIZEOF_SCALAR

  call endTiming("File name and size check")

  z = 0.0_wp
  globalOffset = 0
  globalBufferSize = MERGE(globalSize, int(1e10,MPI_OFFSET_KIND),               &
                           globalSize < int(1e10,MPI_OFFSET_KIND))

  do while (globalOffset<globalSize)

    call startTiming("Buffer setup")

    if ( globalSize - globalOffset < globalBufferSize )                         &
      globalBufferSize = globalSize - globalOffset

    indexQuotient = MOD(globalBufferSize,int(numProcs,MPI_OFFSET_KIND))
    bufferSize = globalBufferSize/int(numProcs,MPI_OFFSET_KIND)
    allocate(recvcnt(0:numProcs-1))
    allocate(displc(0:numProcs-1))
    recvcnt(0:indexQuotient-1) = bufferSize+1
    recvcnt(indexQuotient:numProcs-1) = bufferSize
    displc = 0
    do i=1,numProcs-1
      displc(i) = SUM(recvcnt(0:i-1))
    end do
    offset = displc(procRank)*SIZEOF_SCALAR + globalOffset

    if( procRank<indexQuotient ) then
      sendcnt = bufferSize+1
    else
      sendcnt = bufferSize
    end if
    allocate(XBuffer(sendcnt))
    allocate(YBuffer(sendcnt))
    allocate(normBuffer(sendcnt))
    XBuffer = 0.0_wp
    YBuffer = 0.0_wp
    normBuffer = 1.0_wp

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

    call MPI_File_open(comm, trim(YFilename),                       &
                       MPI_MODE_RDONLY, MPI_INFO_NULL,                        &
                       mpiFileHandle, ierror)
    call MPI_FILE_SET_VIEW(mpiFileHandle,offset,MPI_DOUBLE,MPI_DOUBLE,        &
                           'native',MPI_INFO_NULL,ierror)
    call MPI_File_read_all(mpiFileHandle, YBuffer, int(sendcnt),                   &
                           SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
    call MPI_File_close(mpiFileHandle, ierror)

    call MPI_Barrier(comm, ierror)

    if (normFileExists) then
      call MPI_File_open(comm, trim(normFilename),                       &
                         MPI_MODE_RDONLY, MPI_INFO_NULL,                        &
                         mpiFileHandle, ierror)
      call MPI_FILE_SET_VIEW(mpiFileHandle,offset,MPI_DOUBLE,MPI_DOUBLE,        &
                             'native',MPI_INFO_NULL,ierror)
      call MPI_File_read_all(mpiFileHandle, normBuffer, int(sendcnt),                   &
                             SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
      call MPI_File_close(mpiFileHandle, ierror)
    end if
    call MPI_Barrier(comm, ierror)

    call endTiming("Read files")
    call startTiming("Compute z=X^T*norm*Y")

    dummyValue = 0.0_wp
    dummyValue = SUM( XBuffer*normBuffer*YBuffer )

    call MPI_Allreduce(MPI_IN_PLACE, dummyValue, 1, SCALAR_TYPE_MPI, MPI_SUM, comm, ierror)
    call MPI_Bcast(dummyValue, 1, SCALAR_TYPE_MPI, 0, comm, ierror)

    z = z + dummyValue

    call endTiming("Compute z=X^T*norm*Y")

    globalOffset = globalOffset + globalBufferSize

  end do

end function
