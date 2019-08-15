#include "config.h"

program extend_control_forcing

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use ErrorHandler
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  character(len = STRING_LENGTH) :: ZFilename, XFilename, message,                                    &
                                    durationStr, totalTimestepStr
  integer :: ierror
  integer :: startTimestep, duration, totalTimestep

  ! Initialize MPI.
  call MPI_Init(ierror)
  call initializeErrorHandler()

  if (command_argument_count() < 4) then
     write(message, '(A)') "Require 4 arguments!"
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     call MPI_Finalize(ierror)
     stop -1
  end if

  call startTiming("total")

  call get_command_argument(1, ZFilename)

  call get_command_argument(2, XFilename)

  call get_command_argument(3, totalTimestepStr)
  read( totalTimestepStr, * ) totalTimestep

  call get_command_argument(4, durationStr)
  read( durationStr, * ) duration

  call extendControlForcing(MPI_COMM_WORLD,trim(ZFilename),trim(XFilename),                    &
                            totalTimestep,duration)

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

contains

  subroutine extendControlForcing(comm,ZFilename,XFilename,XFileTimestep,extension)

    use MPI
    use, intrinsic :: iso_fortran_env, only : output_unit

    use InputHelper, only : parseInputFile, getOption, getRequiredOption
    use ErrorHandler
    use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: comm, XFileTimestep, extension
    character(len = *), intent(in) :: ZFilename, XFilename

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, procRank, numProcs, ierror, mpiFileHandle
    character(len = STRING_LENGTH) :: message
    logical :: XFileExists, success
    integer(kind = MPI_OFFSET_KIND) :: XFileSize, ZFileSize, offset, extensionOffset

    integer(kind = MPI_OFFSET_KIND) :: XArraySize, patchSize, extensionSize, bufferSize, sendcnt, indexQuotient
    integer(kind = MPI_OFFSET_KIND), dimension(:), allocatable :: recvcnt, displc
    SCALAR_TYPE, dimension(:), allocatable :: XBuffer, extensionBuffer
    SCALAR_TYPE :: dummyValue

    ! Initialize MPI.
    call MPI_Comm_rank(comm, procRank, ierror)
    call MPI_Comm_size(comm, numProcs, ierror)

    call startTiming("File name and size check")

    write(message,'(2A)') "Z filename: ",trim(ZFilename)
    call writeAndFlush(comm, output_unit, message)

    write(message,'(A,I8.3)') "original timestep: ", XFileTimestep
    call writeAndFlush(comm, output_unit, message)

    write(message,'(A,I8.3)') "extension duration: ", extension
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

      XArraySize = XFileSize/SIZEOF_SCALAR
    end if
    if (MOD(XArraySize,int(4*XFileTimestep,MPI_OFFSET_KIND)).ne.0) then
      write(message, '(A)') "X file size does not match with the given original time steps!"
      call gracefulExit(comm, message)
    else
      patchSize = XArraySize / 4 / XFileTimestep
      extensionSize = patchSize * 4 * extension
      extensionOffset = extensionSize * SIZEOF_SCALAR
    end if

    call endTiming("File name and size check")
    call startTiming("Extension Buffer Setup")

    indexQuotient = MOD(extensionSize,int(numProcs,MPI_OFFSET_KIND))
    bufferSize = extensionSize/int(numProcs,MPI_OFFSET_KIND)
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
    allocate(extensionBuffer(sendcnt))
    extensionBuffer = 0.0_wp

    call endTiming("Extension Buffer Setup")
    call startTiming("Write file")

    call MPI_Barrier(comm, ierror)

    call MPI_File_open(comm, trim(ZFilename),                                 &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL,      &
                       mpiFileHandle, ierror)
    call MPI_FILE_SET_VIEW(mpiFileHandle,offset,MPI_DOUBLE,MPI_DOUBLE,        &
                           'native',MPI_INFO_NULL,ierror)
    call MPI_File_write_all(mpiFileHandle, extensionBuffer, int(sendcnt),     &
                            SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
    call MPI_File_close(mpiFileHandle, ierror)

    call MPI_Barrier(comm, ierror)

    call endTiming("Write file")
    call startTiming("X Buffer setup")

    indexQuotient = MOD(XArraySize,int(numProcs,MPI_OFFSET_KIND))
    bufferSize = XArraySize/int(numProcs,MPI_OFFSET_KIND)
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
    XBuffer = 0.0_wp

    call endTiming("X Buffer setup")
    call startTiming("Read X file")

    call MPI_Barrier(comm, ierror)

    call MPI_File_open(comm, trim(XFilename),                                 &
                       MPI_MODE_RDONLY, MPI_INFO_NULL,                        &
                       mpiFileHandle, ierror)
    call MPI_FILE_SET_VIEW(mpiFileHandle,offset,MPI_DOUBLE,MPI_DOUBLE,        &
                           'native',MPI_INFO_NULL,ierror)
    call MPI_File_read_all(mpiFileHandle, XBuffer, int(sendcnt),                   &
                           SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)

    call MPI_File_close(mpiFileHandle, ierror)

    call MPI_Barrier(comm, ierror)

    call endTiming("Read X file")
    call startTiming("Write file")

    call MPI_Barrier(comm, ierror)

    call MPI_File_open(comm, trim(ZFilename),                                 &
                       MPI_MODE_RDWR, MPI_INFO_NULL,                          &
                       mpiFileHandle, ierror)
    call MPI_FILE_SET_VIEW(mpiFileHandle,offset+extensionOffset,MPI_DOUBLE,MPI_DOUBLE,        &
                           'native',MPI_INFO_NULL,ierror)
    call MPI_File_write_all(mpiFileHandle, XBuffer, int(sendcnt),                   &
                            SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
    call MPI_File_close(mpiFileHandle, ierror)

    call MPI_Barrier(comm, ierror)

    call endTiming("Write file")

  end subroutine

end program extend_control_forcing
