#include "config.h"

program paste_control_forcing

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use ErrorHandler
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  character(len = STRING_LENGTH) :: ZFilename, XFilename, message,                                    &
                                    startTimestepStr, durationStr, totalTimestepStr,                  &
                                    overlapStr
  integer :: ierror
  integer :: startTimestep, duration, totalTimestep
  logical :: overlapFile

  ! Initialize MPI.
  call MPI_Init(ierror)
  call initializeErrorHandler()

  if (command_argument_count() < 5) then
     write(message, '(A)') "Require 5 arguments!"
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     call MPI_Finalize(ierror)
     stop -1
  elseif (command_argument_count() > 6) then
    write(message, '(A)') "Too many arguments!"
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
    call MPI_Finalize(ierror)
    stop -1
  end if

  call startTiming("total")

  call get_command_argument(1, ZFilename)

  call get_command_argument(2, XFilename)

  call get_command_argument(3, totalTimestepStr)
  read( totalTimestepStr, * ) totalTimestep

  call get_command_argument(4, startTimestepStr)
  read( startTimestepStr, * ) startTimestep

  call get_command_argument(5, durationStr)
  read( durationStr, * ) duration

  overlapFile = .false.
  if (command_argument_count() == 6) then
    call get_command_argument(6, overlapStr)
    if (trim(overlapStr)=='--overlap') overlapFile = .true.
  end if

  call pasteControlForcing(MPI_COMM_WORLD,trim(ZFilename),trim(XFilename),                    &
                            totalTimestep,startTimestep,duration,overlapFile)

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

contains

  subroutine pasteControlForcing(comm,ZFilename,XFilename,totalTimestep,startTimestep,duration,overlapFile)

    use MPI
    use, intrinsic :: iso_fortran_env, only : output_unit

    use InputHelper, only : parseInputFile, getOption, getRequiredOption
    use ErrorHandler
    use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: comm, totalTimestep, startTimestep, duration
    character(len = *), intent(in) :: ZFilename, XFilename
    logical, intent(in) :: overlapFile

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, procRank, numProcs, ierror, mpiFileHandle
    character(len = STRING_LENGTH) :: message
    logical :: XFileExists, ZFileExists
    integer(kind = MPI_OFFSET_KIND) :: XFileSize, offset, globalOffset

    integer(kind = MPI_OFFSET_KIND) :: globalSize, patchSize, sliceSize, bufferSize, sendcnt, indexQuotient
    integer(kind = MPI_OFFSET_KIND), dimension(:), allocatable :: recvcnt, displc
    SCALAR_TYPE, dimension(:), allocatable :: XBuffer, ZBuffer

    ! Initialize MPI.
    call MPI_Comm_rank(comm, procRank, ierror)
    call MPI_Comm_size(comm, numProcs, ierror)

    call startTiming("File name and size check")

    write(message,'(2A)') "Z filename: ",trim(ZFilename)
    call writeAndFlush(comm, output_unit, message)

    write(message,'(A,I8.3)') "total timestep: ", totalTimestep
    call writeAndFlush(comm, output_unit, message)

    write(message,'(A,I8.3)') "start timestep: ", startTimestep
    call writeAndFlush(comm, output_unit, message)

    write(message,'(A,I8.3)') "duration: ", duration
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

      sliceSize = XFileSize/SIZEOF_SCALAR
    end if
    if (MOD(sliceSize,int(4*duration,MPI_OFFSET_KIND)).ne.0) then
      write(message, '(A)') "X file size does not match with the given duration!"
      call gracefulExit(comm, message)
    else
      patchSize = sliceSize / 4 / duration
      globalSize = patchSize * 4 * totalTimestep
      globalOffset = globalSize - patchSize * 4 * (startTimestep+duration)
      globalOffset = globalOffset * SIZEOF_SCALAR
    end if

    ! Check Z file if overlap file.
    if (overlapFile) then
      if (procRank == 0) then
         inquire(file = trim(ZFilename), exist = ZFileExists)
      end if
      call MPI_Bcast(ZFileExists, 1, MPI_LOGICAL, 0, comm, ierror)
      if (.not.ZFileExists) then
         write(message, '(A)') "Z file does not exists!"
         call gracefulExit(comm, message)
      end if
    end if

    call endTiming("File name and size check")
    call startTiming("Buffer setup")

    indexQuotient = MOD(sliceSize,int(numProcs,MPI_OFFSET_KIND))
    bufferSize = sliceSize/int(numProcs,MPI_OFFSET_KIND)
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
    allocate(ZBuffer(sendcnt))
    XBuffer = 0.0_wp
    ZBuffer = 0.0_wp

    call endTiming("Buffer setup")
    call startTiming("Read files")

    call MPI_Barrier(comm, ierror)

    ! read global file
    if (overlapFile) then
      call MPI_File_open(comm, trim(ZFilename),                                 &
                         MPI_MODE_RDONLY, MPI_INFO_NULL,                        &
                         mpiFileHandle, ierror)
      call MPI_FILE_SET_VIEW(mpiFileHandle,offset+globalOffset,MPI_DOUBLE,MPI_DOUBLE,        &
                             'native',MPI_INFO_NULL,ierror)
      call MPI_File_read_all(mpiFileHandle, ZBuffer, int(sendcnt),                   &
                             SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
      call MPI_File_close(mpiFileHandle, ierror)
    end if

    ! read slice file
    call MPI_File_open(comm, trim(XFilename),                                 &
                       MPI_MODE_RDONLY, MPI_INFO_NULL,                        &
                       mpiFileHandle, ierror)
    call MPI_FILE_SET_VIEW(mpiFileHandle,offset,MPI_DOUBLE,MPI_DOUBLE,        &
                           'native',MPI_INFO_NULL,ierror)
    call MPI_File_read_all(mpiFileHandle, XBuffer, int(sendcnt),                   &
                           SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
    call MPI_File_close(mpiFileHandle, ierror)

    call MPI_Barrier(comm, ierror)

    call endTiming("Read files")
    call startTiming("Z(slice) += X")

    ZBuffer = ZBuffer + XBuffer

    call endTiming("Z(slice) += X")
    call startTiming("Write file")

    call MPI_Barrier(comm, ierror)

    call MPI_File_open(comm, trim(ZFilename), MPI_MODE_WRONLY + MPI_MODE_CREATE,             &
                      MPI_INFO_NULL, mpiFileHandle, ierror)
    call MPI_FILE_SET_VIEW(mpiFileHandle,offset+globalOffset,MPI_DOUBLE,MPI_DOUBLE,          &
                           'native',MPI_INFO_NULL,ierror)
    call MPI_File_write_all(mpiFileHandle, ZBuffer, int(sendcnt),               &
                            SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
    call MPI_File_close(mpiFileHandle, ierror)

    call MPI_Barrier(comm, ierror)

    call endTiming("Write file")

  end subroutine

  !NOTE: this subroutine is unused and thus commented out. Should be working.
  ! subroutine createZeroControlForcing(comm,ZFilename,globalSize)
  !
  !   use MPI
  !   use, intrinsic :: iso_fortran_env, only : output_unit
  !
  !   use InputHelper, only : parseInputFile, getOption, getRequiredOption
  !   use ErrorHandler
  !   use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers
  !
  !   implicit none
  !
  !   ! <<< Arguments >>>
  !   integer, intent(in) :: comm
  !   integer(kind = MPI_OFFSET_KIND), intent(in) :: globalSize
  !   character(len = *), intent(in) :: ZFilename
  !
  !   ! <<< Local variables >>>
  !   integer, parameter :: wp = SCALAR_KIND
  !   integer :: i, procRank, numProcs, ierror, mpiFileHandle
  !   character(len = STRING_LENGTH) :: message
  !   integer(kind = MPI_OFFSET_KIND) :: offset
  !
  !   integer(kind = MPI_OFFSET_KIND) :: bufferSize, sendcnt, indexQuotient
  !   integer(kind = MPI_OFFSET_KIND), dimension(:), allocatable :: recvcnt, displc
  !   SCALAR_TYPE, dimension(:), allocatable :: ZBuffer
  !
  !   ! Initialize MPI.
  !   call MPI_Comm_rank(comm, procRank, ierror)
  !   call MPI_Comm_size(comm, numProcs, ierror)
  !
  !   call startTiming("File name and size check")
  !
  !   write(message,'(3A)') "Creating zero ", trim(ZFilename), "... "
  !   call writeAndFlush(comm, output_unit, message, advance = 'no')
  !
  !   call endTiming("File name and size check")
  !   call startTiming("Buffer setup")
  !
  !   indexQuotient = MOD(globalSize,int(numProcs,MPI_OFFSET_KIND))
  !   bufferSize = globalSize/int(numProcs,MPI_OFFSET_KIND)
  !   allocate(recvcnt(0:numProcs-1))
  !   allocate(displc(0:numProcs-1))
  !   recvcnt(0:indexQuotient-1) = bufferSize+1
  !   recvcnt(indexQuotient:numProcs-1) = bufferSize
  !   displc = 0
  !   do i=1,numProcs-1
  !     displc(i) = SUM(recvcnt(0:i-1))
  !   end do
  !   offset = displc(procRank)*SIZEOF_SCALAR
  !
  !   if( procRank<indexQuotient ) then
  !     sendcnt = bufferSize+1
  !   else
  !     sendcnt = bufferSize
  !   end if
  !   allocate(ZBuffer(sendcnt))
  !   ZBuffer = 0.0_wp
  !
  !   call endTiming("Buffer setup")
  !   call startTiming("Write file")
  !
  !   call MPI_Barrier(comm, ierror)
  !
  !   call MPI_File_open(comm, trim(ZFilename),                                 &
  !                      MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL,      &
  !                      mpiFileHandle, ierror)
  !   call MPI_FILE_SET_VIEW(mpiFileHandle,offset,MPI_DOUBLE,MPI_DOUBLE,        &
  !                          'native',MPI_INFO_NULL,ierror)
  !   call MPI_File_write_all(mpiFileHandle, ZBuffer, int(sendcnt),                   &
  !                           SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
  !   call MPI_File_close(mpiFileHandle, ierror)
  !
  !   call MPI_Barrier(comm, ierror)
  ! 
  !   call endTiming("Write file")
  !
  !   write(message, '(A)') " done!"
  !   call writeAndFlush(comm, output_unit, message)
  !
  ! end subroutine

end program paste_control_forcing
