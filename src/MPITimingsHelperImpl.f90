#include "config.h"

module MPITimingsHelperImpl

  implicit none
  private

  integer, parameter, private :: real64 = selected_real_kind(15)

  type, public :: t_MPITimer
     character(len = STRING_LENGTH) :: name
     real(kind = real64) :: numCalls = 0.0_real64, startTimeOfLastCall = 0.0_real64, &
          accumulatedTime = 0.0_real64
     type(t_MPITimer), pointer :: left => null(), right => null()
  end type t_MPITimer

  real(kind = real64), public :: programStartTime
  type(t_MPITimer), pointer, public :: timings => null()

  public :: addNode, deleteSubtree, getNumberOfLeaves, serializeTree, sortTimings

contains

  subroutine addNode(root, node)

    ! <<< Arguments >>>
    type(t_MPITimer), pointer :: root
    type(t_MPITimer) :: node

    ! <<< Local variables >>>
    type(t_MPITimer), pointer :: current => null()

    if (.not. associated(root)) then

       allocate(root)
       root%name = trim(node%name)
       root%numCalls = node%numCalls
       root%accumulatedTime = node%accumulatedTime

    else

       current => root

       do
          if (trim(current%name) == trim(node%name)) then
             current%numCalls = current%numCalls + node%numCalls
             current%accumulatedTime = current%accumulatedTime + node%accumulatedTime
             exit
          else if (trim(node%name) > trim(current%name)) then
             if (.not. associated(current%left)) then
                allocate(current%left)
                current%left%name = trim(node%name)
                current%left%numCalls = node%numCalls
                current%left%accumulatedTime = node%accumulatedTime
                exit
             else
                current => current%left
             end if
          else
             if (.not. associated(current%right)) then
                allocate(current%right)
                current%right%name = trim(node%name)
                current%right%numCalls = node%numCalls
                current%right%accumulatedTime = node%accumulatedTime
                exit
             else
                current => current%right
             end if
          end if
       end do

    end if

  end subroutine addNode

  recursive subroutine deleteSubtree(node)

    ! <<< Arguments >>>
    type(t_MPITimer), pointer :: node

    if (.not. associated(node)) return
    if (associated(node%left)) call deleteSubtree(node%left)
    if (associated(node%right)) call deleteSubtree(node%right)
    deallocate(node)

  end subroutine deleteSubtree

  recursive function getNumberOfLeaves(node) result(nLeaves)

    ! <<< Arguments >>>
    type(t_MPITimer), pointer, intent(in) :: node

    ! <<< Result >>>
    integer :: nLeaves

    nLeaves = 0

    if (associated(node%left)) then
       nLeaves = nLeaves + getNumberOfLeaves(node%left) + 1
    end if

    if (associated(node%right)) then
       nLeaves = nLeaves + getNumberOfLeaves(node%right) + 1
    end if

  end function getNumberOfLeaves

  subroutine serializeTree(root, array)

    ! <<< Arguments >>>
    type(t_MPITimer), pointer :: root
    type(t_MPITimer), allocatable, intent(out) :: array(:)

    ! <<< Local variables >>>
    integer :: i, nNodes
    type(t_MPITimer), pointer :: current => null(), pre => null()

    nNodes = getNumberOfLeaves(root) + 1
    if (nNodes == 0) return

    allocate(array(nNodes))

    i = 0

    current => root

    do while (associated(current))

       if (.not. associated(current%left)) then

          i = i + 1
          array(i)%name = current%name
          array(i)%numCalls = current%numCalls
          array(i)%accumulatedTime = current%accumulatedTime
          current => current%right

       else

          pre => current%left

          do while (associated(pre%right))
             if (trim(pre%right%name) == trim(current%name)) exit
             pre => pre%right
          end do

          if (.not. associated(pre%right)) then

             pre%right => current
             current => current%left

          else

             nullify(pre%right)
             i = i + 1
             array(i)%name = current%name
             array(i)%numCalls = current%numCalls
             array(i)%accumulatedTime = current%accumulatedTime
             current => current%right

          end if

       end if

    end do

  end subroutine serializeTree

  subroutine sortTimings(timingsArray)

    ! <<< Arguments >>>
    type(t_MPITimer), intent(inout) :: timingsArray(:)
    type(t_MPITimer) :: temp

    ! <<< Local variables >>>
    integer :: i, j
    
    do i = 2, size(timingsArray)

       j = i - 1

       temp%name = timingsArray(i)%name
       temp%numCalls = timingsArray(i)%numCalls
       temp%accumulatedTime = timingsArray(i)%accumulatedTime

       do while (timingsArray(j)%accumulatedTime < temp%accumulatedTime)
          timingsArray(j+1)%name = timingsArray(j)%name
          timingsArray(j+1)%numCalls = timingsArray(j)%numCalls
          timingsArray(j+1)%accumulatedTime = timingsArray(j)%accumulatedTime
          j = j - 1
          if (j < 1) exit
       end do

       timingsArray(j+1)%name = temp%name
       timingsArray(j+1)%numCalls = temp%numCalls
       timingsArray(j+1)%accumulatedTime = temp%accumulatedTime

    end do
    
  end subroutine sortTimings

end module MPITimingsHelperImpl

subroutine startTiming(name)

  ! <<< External modules >>>
  use MPI

  ! <<< Private members >>>
  use MPITimingsHelperImpl, only : t_MPITimer, timings, programStartTime

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: name

  ! <<< Local variables >>>
  type(t_MPITimer), pointer :: current => null()

  if (.not. associated(timings)) then

     allocate(timings)
     timings%name = trim(name)
     timings%startTimeOfLastCall = MPI_Wtime()
     programStartTime = timings%startTimeOfLastCall

  else

     current => timings

     do
        if (trim(current%name) == trim(name)) then
           current%startTimeOfLastCall = MPI_Wtime()
           exit
        else if (trim(name) < trim(current%name)) then
           if (.not. associated(current%left)) then
              allocate(current%left)
              current%left%name = trim(name)
           end if
           current => current%left
        else
           if (.not. associated(current%right)) then
              allocate(current%right)
              current%right%name = trim(name)
           end if
           current => current%right
        end if
     end do

  end if

end subroutine startTiming

subroutine endTiming(name)

  ! <<< External modules >>>
  use MPI

  ! <<< Private members >>>
  use MPITimingsHelperImpl, only : t_MPITimer, timings

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: name

  ! <<< Local variables >>>
  integer, parameter :: real64 = selected_real_kind(15)
  type(t_MPITimer), pointer :: current => null()

  current => timings

  do while(associated(current))
     if (trim(current%name) == trim(name)) then
        current%numCalls = current%numCalls + 1.0_real64
        current%accumulatedTime = current%accumulatedTime +                                  &
             MPI_Wtime() - current%startTimeOfLastCall
        exit
     else if (trim(name) < trim(current%name)) then
        if (.not. associated(current%left)) allocate(current%left)
        current => current%left
     else
        if (.not. associated(current%right)) allocate(current%right)
        current => current%right
     end if
  end do

end subroutine endTiming

subroutine reportTimings(comm, outputUnit)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  ! <<< Private members >>>
  use MPITimingsHelperImpl

  implicit none

  ! <<< Arguments >>>
  integer, intent(in), optional :: comm, outputUnit

  ! <<< Local variables >>>
  integer, parameter :: real64 = selected_real_kind(15)
  integer :: i, j, comm_, outputUnit_, nTimings, procRank, nProcs, ierror
  real(kind = real64) :: programTotalTime
  type(t_MPITimer), allocatable :: timingsArray(:)
  integer, allocatable :: nTimingsAllProcesses(:)
  type(t_MPITimer), pointer :: globalTimings => null()
  type(t_MPITimer) :: temp

  comm_ = MPI_COMM_WORLD
  if (present(comm)) comm_ = comm

  outputUnit_ = output_unit
  if (present(outputUnit)) outputUnit_ = outputUnit

  ! Get rank and size of the MPI communicator `comm_`.
  call MPI_Comm_rank(comm_, procRank, ierror)
  call MPI_Comm_size(comm_, nProcs, ierror)

  ! Find the total time taken by the program.
  programTotalTime = (MPI_Wtime() - programStartTime) * real(nProcs, real64)

  ! Serialize the timings from each process.
  call serializeTree(timings, timingsArray)

  ! Find the number of timers on each process.
  nTimings = 0
  if (allocated(timingsArray)) nTimings = size(timingsArray)

  ! Not all processes may have timed the same number of subroutines... gather the number of
  ! timers from all processes.
  allocate(nTimingsAllProcesses(nProcs))
  call MPI_Allgather(nTimings, 1, MPI_INTEGER, nTimingsAllProcesses,                         &
       1, MPI_INTEGER, comm_, ierror)

  ! Master process adds all its timers to the `globalTimings` tree.
  if (procRank == 0) then
     do i = 1, nTimings
        call addNode(globalTimings, timingsArray(i))
     end do
  end if

  do i = 1, nProcs - 1

     if (procRank == i) then

        ! All other processes send their timings data to master.
        do j = 1, nTimings
           call MPI_Send(timingsArray(j)%name, STRING_LENGTH, MPI_CHARACTER,                 &
                0, i * 3 + 0, comm_, ierror)
           call MPI_Send(timingsArray(j)%numCalls, 1, MPI_REAL8, 0,                          &
                i * 3 + 1, comm_, ierror)
           call MPI_Send(timingsArray(j)%accumulatedTime, 1, MPI_REAL8,                      &
                0, i * 3 + 2, comm_, ierror)
        end do

     else if (procRank == 0) then

        ! Master process received the data and adds them to `globalTimings` tree.
        do j = 1, nTimingsAllProcesses(i + 1)
           call MPI_Recv(temp%name, STRING_LENGTH, MPI_CHARACTER,                            &
                i, i * 3 + 0, comm_, MPI_STATUS_IGNORE, ierror)
           call MPI_Recv(temp%numCalls, 1, MPI_REAL8,                                        &
                i, i * 3 + 1, comm_, MPI_STATUS_IGNORE, ierror)
           call MPI_Recv(temp%accumulatedTime, 1, MPI_REAL8,                                 &
                i, i * 3 + 2, comm_, MPI_STATUS_IGNORE, ierror)
           call addNode(globalTimings, temp)
        end do

     end if

  end do

  SAFE_DEALLOCATE(nTimingsAllProcesses)
  SAFE_DEALLOCATE(timingsArray)

  if (procRank == 0) then

     call serializeTree(globalTimings, timingsArray)

     if (allocated(timingsArray)) then

        call sortTimings(timingsArray)

        write(outputUnit_, '(A)') ""
        write(outputUnit_, '(2A)') PROJECT_NAME,                                             &
             " profile (timers may not be mutually exclusive):"
        write(outputUnit_, '(A)') repeat("=", 80)
        write(outputUnit_, '(A25,5X,A6,2X,A13,2X,A13,2X,A12)')                               &
             "name", "% time", "calls", "total seconds", "seconds/call"
        write(outputUnit_, '(A)') repeat("-", 80)

        do i = 1, size(timingsArray)
           write(outputUnit_, '(A25,5X,F6.2,2X,I13,2X,F13.4,2X,ES12.4)')                     &
                trim(timingsArray(i)%name),                                                  &
                timingsArray(i)%accumulatedTime / programTotalTime * 100.0_real64,           &
                nint(timingsArray(i)%numCalls), timingsArray(i)%accumulatedTime,             &
                timingsArray(i)%accumulatedTime / timingsArray(i)%numCalls
        end do

        write(outputUnit_, '(A)') repeat("=", 80)
        write(outputUnit_, '(A)') ""

        ! Flush `unit` from master process.
        flush(outputUnit_)

     end if

  end if

  ! Ensure that no writes are executed from other processes in `comm` before master process
  ! has flushed `unit`.
  call MPI_Barrier(comm_, ierror)

end subroutine reportTimings

subroutine cleanupTimers()

  ! <<< Private members >>>
  use MPITimingsHelperImpl, only : timings, deleteSubtree

  implicit none

  call deleteSubtree(timings)

end subroutine cleanupTimers
