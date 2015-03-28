#include "config.h"

program dissipation_self_adjoint

  use MPI

  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random

  implicit none

  integer :: i, direction, ierror
  logical :: success, success_

  interface

     subroutine testSelfAdjointness(identifier, direction, success, isPeriodic, tolerance)

       character(len = *), intent(in) :: identifier
       integer, intent(in) :: direction
       logical, intent(out) :: success

       logical, intent(in), optional :: isPeriodic
       real(SCALAR_KIND), intent(in), optional :: tolerance

     end subroutine testSelfAdjointness

  end interface

  call MPI_Init(ierror)

  success = .true.

  call initializeErrorHandler()

  do i = 1, 10 !... test multiple times
     do direction = 1, 3

        call testSelfAdjointness("SBP 1-2 dissipation",                                      &
             direction, success_, isPeriodic = .false.)
        success = success .and. success_
        if (.not. success) exit

        call testSelfAdjointness("SBP 1-2 dissipation",                                      &
             direction, success_, isPeriodic = .true.)
        success = success .and. success_
        if (.not. success) exit

        call testSelfAdjointness("SBP 2-4 dissipation",                                      &
             direction, success_, isPeriodic = .false.)
        success = success .and. success_
        if (.not. success) exit

        call testSelfAdjointness("SBP 2-4 dissipation",                                      &
             direction, success_, isPeriodic = .true.)
        success = success .and. success_
        if (.not. success) exit

        call testSelfAdjointness("SBP 3-6 dissipation",                                      &
             direction, success_, isPeriodic = .false.)
        success = success .and. success_
        if (.not. success) exit

        call testSelfAdjointness("SBP 3-6 dissipation",                                      &
             direction, success_, isPeriodic = .true.)
        success = success .and. success_
        if (.not. success) exit

        call testSelfAdjointness("SBP 4-8 dissipation",                                      &
             direction, success_, isPeriodic = .false.)
        success = success .and. success_
        if (.not. success) exit

        call testSelfAdjointness("SBP 4-8 dissipation",                                      &
             direction, success_, isPeriodic = .true.)
        success = success .and. success_
        if (.not. success) exit

     end do
  end do

  call cleanupErrorHandler()

  call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierror)
  call MPI_Finalize(ierror)
  if (.not. success) stop -1
  stop 0

end program dissipation_self_adjoint

subroutine testSelfAdjointness(identifier, direction, success, isPeriodic, tolerance)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use StencilOperator_type, only : t_StencilOperator

  ! <<< Internal modules >>>
  use MPIHelper, only : pigeonhole
  use RandomNumber, only : initializeRandomNumberGenerator, random
  use StencilOperator_mod, only : setupOperator, updateOperator, &
       getAdjointOperator, applyOperator, cleanupOperator

  ! <<< Arguments >>>
  character(len = *), intent(in) :: identifier
  integer, intent(in) :: direction
  logical, intent(out) :: success
  logical, intent(in), optional :: isPeriodic
  real(SCALAR_KIND), intent(in), optional :: tolerance

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  logical :: isPeriodic_(3), success_
  real(wp) :: tolerance_
  integer :: n, offset, nLocal, gridSize(3), cartesianCommunicator,                          &
       numProcesses(3), procRank, nProcs, ierror
  type(t_StencilOperator) :: A, adjointOfA
  real(wp), allocatable :: f(:,:)
  SCALAR_TYPE, allocatable :: u(:,:), v(:,:)

  if (direction < 0 .or. direction > 3) then
     success = .false.
     return
  end if

  isPeriodic_ = .false.
  if (present(isPeriodic)) isPeriodic_(direction) = isPeriodic

  tolerance_ = epsilon(0.0_wp)
  if (present(tolerance)) tolerance_ = tolerance

  call initializeRandomNumberGenerator()

  ! Find the rank and number of processes in the communicator.
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, nProcs, ierror)
#ifdef SCALAR_TYPE_IS_binary128_IEEE754
  allocate(mpiReduceBuffer(nProcs))
#endif

  numProcesses = 1
  numProcesses(direction) = nProcs

  ! Create a Cartesian communicator.
  call MPI_Cart_create(MPI_COMM_WORLD, 3, numProcesses, isPeriodic_,                         &
       .true., cartesianCommunicator, ierror)

  ! Update the operators.
  call setupOperator(A, identifier, success_)
  if (.not. success_) then
     success = .false.
     return
  end if
  call updateOperator(A, cartesianCommunicator, direction)
  call getAdjointOperator(A, adjointOfA)
  call updateOperator(adjointOfA, cartesianCommunicator, direction)

  ! Decide the size of arrays to be used for the test.
  n = random(nProcs * adjointOfA%boundaryDepth, 2 ** 16)
  call MPI_Bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

  ! Determine the offset and number of points that will be distributed to the current
  ! process.
  call pigeonhole(n, nProcs, procRank, offset, nLocal)
  gridSize = 1; gridSize(direction) = nLocal

  ! Allocate process-level data.
  allocate(f(nLocal, 1), u(nLocal, 1))

  ! Initialize `f` to random values.
  call random_number(f)

  u = f
  v = f
  call applyOperator(A, u, gridSize)
  call applyOperator(adjointOfA, v, gridSize)

  success = (maxval(abs(u - v)) / real(n, wp) <= tolerance_)

  SAFE_DEALLOCATE(f)
  SAFE_DEALLOCATE(u)
  SAFE_DEALLOCATE(v)

  call cleanupOperator(A)
  call cleanupOperator(adjointOfA)

end subroutine testSelfAdjointness
