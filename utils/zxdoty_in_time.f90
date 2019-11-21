#include "config.h"

program run_zXdotY

  use ControlSpaceAdvancer, only : zXdotY
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use ErrorHandler
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers
  use InputHelper, only : getFreeUnit

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: zFilename, XFilename, YFilename, normFilename, message
  character(len = STRING_LENGTH) :: nTimesteps_str, argument
  integer :: nTimesteps, reportInterval = 1, stat, fileUnit, procRank, ierror
  integer :: kthArgument, numberOfArguments
  logical :: lookForReportInterval = .false.

  ! Initialize MPI.
  call MPI_Init(ierror)
  call initializeErrorHandler()

  if (command_argument_count() < 5) then
     write(message, '(A)') "Require 5 arguments!"
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     call MPI_Finalize(ierror)
     stop -1
  end if

  call startTiming("total")

  call get_command_argument(1, zFilename)

  call get_command_argument(2, XFilename)

  call get_command_argument(3, YFilename)

  call get_command_argument(4, normFilename)

  call get_command_argument(5, nTimesteps_str)
  read( nTimesteps_str, * ) nTimesteps

  numberOfArguments = command_argument_count()
  do kthArgument = 6, numberOfArguments
    call get_command_argument(kthArgument,argument)
    select case(trim(adjustl(argument)))
    case("--report_interval")
      lookForReportInterval = .true.
    case default
      if (lookForReportInterval) then
        read( argument, * ) reportInterval
        assert(reportInterval>0)
        lookForReportInterval = .false.
      else
        write(message, '(3A)') "option ", trim(argument), " unknown!"
        call gracefulExit(MPI_COMM_WORLD, message)
      end if
    end select
  end do

  call zXdotYinTime(MPI_COMM_WORLD,trim(zFilename),trim(XFilename),trim(YFilename),           &
                    trim(normFilename),nTimesteps,reportInterval)

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

contains

  subroutine zXdotYinTime(comm,zFilename,XFilename,YFilename,normFilename,nTimesteps,reportInterval_)
    use MPI
    use, intrinsic :: iso_fortran_env, only : output_unit

    use InputHelper, only : parseInputFile, getOption, getRequiredOption
    use ErrorHandler
    use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: comm, nTimesteps, reportInterval_
    character(len = *), intent(in) :: zFilename, XFilename, YFilename, normFilename

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, procRank, numProcs, ierror, mpiFileHandle, timestep, stage
    character(len = STRING_LENGTH) :: message
    logical :: XFileExists, YFileExists, normFileExists, success, append
    integer(kind = MPI_OFFSET_KIND) :: XFileSize, YFileSize, normFileSize, offset, timeOffset

    integer(kind = MPI_OFFSET_KIND) :: globalSize, patchSize, bufferSize, sendcnt, indexQuotient
    integer(kind = MPI_OFFSET_KIND), dimension(:), allocatable :: recvcnt, displc
    SCALAR_TYPE, dimension(:), allocatable :: XBuffer, YBuffer, normBuffer
    SCALAR_TYPE :: z, runningTimeQuadrature, RK4norm(4)

    RK4norm = (/ 1.0_wp/6.0_wp, 1.0_wp/3.0_wp, 1.0_wp/3.0_wp, 1.0_wp/6.0_wp /)

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

    if (XFileSize.ne.YFileSize                                                    &
    .or. (normFileExists .and. normFileSize.ne.XFileSize)) then
      write(message, '(A)') "X, Y (and norm) do not match in size!"
      call gracefulExit(comm, message)
    end if
    globalSize = XFileSize/SIZEOF_SCALAR
    if (MOD(globalSize,int(4*nTimesteps,MPI_OFFSET_KIND))==0) then
      write(message, '(A)') "Number of timesteps does not match with global file size!"
      call gracefulExit(comm, message)
    end if
    patchSize = globalSize / 4 / nTimesteps

    call endTiming("File name and size check")
    call startTiming("Buffer setup")

    indexQuotient = MOD(patchSize,int(numProcs,MPI_OFFSET_KIND))
    bufferSize = patchSize/int(numProcs,MPI_OFFSET_KIND)
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
    allocate(normBuffer(sendcnt))
    XBuffer = 0.0_wp
    YBuffer = 0.0_wp
    normBuffer = 1.0_wp

    call endTiming("Buffer setup")

    runningTimeQuadrature = 0.0_wp

    do timestep = nTimesteps, 1, -1
      do stage = 4, 1, -1

        call startTiming("Read files")

        timeOffset = int( 4 * (nTimesteps-timestep) + (4-stage) , MPI_OFFSET_KIND)

        call MPI_Barrier(comm, ierror)

        call MPI_File_open(comm, trim(XFilename),                                             &
                           MPI_MODE_RDONLY, MPI_INFO_NULL,                                    &
                           mpiFileHandle, ierror)
        call MPI_FILE_SET_VIEW(mpiFileHandle,timeOffset+offset,MPI_DOUBLE,MPI_DOUBLE,         &
                               'native',MPI_INFO_NULL,ierror)
        call MPI_File_read_all(mpiFileHandle, XBuffer, int(sendcnt),                          &
                               SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
        call MPI_File_close(mpiFileHandle, ierror)

        call MPI_Barrier(comm, ierror)

        call MPI_File_open(comm, trim(YFilename),                                             &
                           MPI_MODE_RDONLY, MPI_INFO_NULL,                                    &
                           mpiFileHandle, ierror)
        call MPI_FILE_SET_VIEW(mpiFileHandle,timeOffset+offset,MPI_DOUBLE,MPI_DOUBLE,         &
                               'native',MPI_INFO_NULL,ierror)
        call MPI_File_read_all(mpiFileHandle, YBuffer, int(sendcnt),                          &
                               SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
        call MPI_File_close(mpiFileHandle, ierror)

        call MPI_Barrier(comm, ierror)

        call MPI_File_open(comm, trim(normFilename),                                          &
                           MPI_MODE_RDONLY, MPI_INFO_NULL,                                    &
                           mpiFileHandle, ierror)
        call MPI_FILE_SET_VIEW(mpiFileHandle,timeOffset+offset,MPI_DOUBLE,MPI_DOUBLE,         &
                               'native',MPI_INFO_NULL,ierror)
        call MPI_File_read_all(mpiFileHandle, normBuffer, int(sendcnt),                       &
                               SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
        call MPI_File_close(mpiFileHandle, ierror)

        call MPI_Barrier(comm, ierror)

        call endTiming("Read files")
        call startTiming("Compute z=X^T*norm*Y")

        z = 0.0_wp
        z = SUM( XBuffer*normBuffer*YBuffer )

        call MPI_Allreduce(MPI_IN_PLACE, z, 1, SCALAR_TYPE_MPI, MPI_SUM, comm, ierror)
        ! call MPI_Bcast(z, 1, SCALAR_TYPE_MPI, 0, comm, ierror)

        runningTimeQuadrature = runningTimeQuadrature + z * RK4norm(stage)

        call endTiming("Compute z=X^T*norm*Y")

      end do

      call startTiming("Write z file")
      if (MOD(nTimesteps-timestep,reportInterval)==0) then
        write(message, '(2A,I8,2(A,E13.6))') PROJECT_NAME, ": timestep = ", timestep,          &
                                             ", z = ", z,                                      &
                                             ", running-time quadrature = ", runningTimeQuadrature
        call writeAndFlush(comm, output_unit, message)

        append = (nTimesteps - timestep > 0)
        call writeQuadratureToFile(comm,trim(zFilename),timestep,z,runningTimeQuadrature,append)
      end if
      call endTiming("Write z file")

    end do

    write(message, '(A)') 'z in time collection is finished.'
    call writeAndFlush(comm, output_unit, message)

  end subroutine

  subroutine writeQuadratureToFile(comm, filename, timestep, z_, runningTimeQuadrature_, append)

    ! <<< External modules >>>
    use MPI

    ! <<< Internal modules >>>
    use InputHelper, only : getFreeUnit
    use ErrorHandler, only : gracefulExit

    ! <<< Arguments >>>
    integer, intent(in) :: comm
    character(len = *), intent(in) :: filename
    integer, intent(in) :: timestep
    SCALAR_TYPE, intent(in) :: z_, runningTimeQuadrature_
    logical, intent(in), optional :: append

    ! <<< Local variables >>>
    logical :: append_
    integer :: fileUnit, ostat, procRank, ierror
    character(len = STRING_LENGTH) :: message

    append_ = .false.
    if (present(append)) append_ = append

    call MPI_Comm_rank(comm, procRank, ierror)

    if (procRank == 0) then
       if (.not. append_) then
          open(unit = getFreeUnit(fileUnit), file = trim(filename), action = 'write',          &
               status = 'unknown', iostat = ostat)
       else
          open(unit = getFreeUnit(fileUnit), file = trim(filename), action = 'write',          &
               status = 'old', position = 'append', iostat = ostat)
       end if
    end if

    call MPI_Bcast(ostat, 1, MPI_INTEGER, 0, comm, ierror)
    if (ostat /= 0) then
       write(message, "(2A)") trim(filename), ": Failed to open file for writing!"
       call gracefulExit(comm, message)
    end if

    if (procRank == 0) then
       write(fileUnit, '(I8,1X,2(1X,SP,' // SCALAR_FORMAT // '))')                       &
            timestep, z_, runningTimeQuadrature_
    end if

    call MPI_Bcast(ostat, 1, MPI_INTEGER, 0, comm, ierror)
    if (ostat /= 0) then
       write(message, "(2A)") trim(filename), ": Error writing to file!"
       call gracefulExit(comm, message)
    end if

    if (procRank == 0) then
       flush(fileUnit)
       close(fileUnit)
    end if

    call MPI_Barrier(comm, ierror)

  end subroutine writeQuadratureToFile

end program run_zXdotY
