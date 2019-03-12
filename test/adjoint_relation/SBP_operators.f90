#include "config.h"

program SBP_operators

  use MPI

  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random

  implicit none

  integer :: i, direction, ierror
  logical :: success, success_

  interface

     subroutine testAdjointRelation(identifier, direction, success, isPeriodic, tolerance)

       character(len = *), intent(in) :: identifier
       integer, intent(in) :: direction
       logical, intent(out) :: success

       logical, intent(in), optional :: isPeriodic
       real(SCALAR_KIND), intent(in), optional :: tolerance

     end subroutine testAdjointRelation

  end interface

  call MPI_Init(ierror)

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  do i = 1, 10 !... test multiple times
     do direction = 1, 3

        call testAdjointRelation("SBP 1-2 first derivative",                                 &
             direction, success_, isPeriodic = .false.)
        success = success .and. success_
        if (.not. success) exit

        call testAdjointRelation("SBP 1-2 first derivative",                                 &
             direction, success_, isPeriodic = .true.)
        success = success .and. success_
        if (.not. success) exit

        call testAdjointRelation("SBP 2-4 first derivative",                                 &
             direction, success_, isPeriodic = .false.)
        success = success .and. success_
        if (.not. success) exit

        call testAdjointRelation("SBP 2-4 first derivative",                                 &
             direction, success_, isPeriodic = .true.)
        success = success .and. success_
        if (.not. success) exit

        call testAdjointRelation("SBP 3-6 first derivative",                                 &
             direction, success_, isPeriodic = .false.)
        success = success .and. success_
        if (.not. success) exit

        call testAdjointRelation("SBP 3-6 first derivative",                                 &
             direction, success_, isPeriodic = .true.)
        success = success .and. success_
        if (.not. success) exit

        call testAdjointRelation("SBP 4-8 first derivative",                                 &
             direction, success_, isPeriodic = .false.)
        success = success .and. success_
        if (.not. success) exit

        call testAdjointRelation("SBP 4-8 first derivative",                                 &
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

end program SBP_operators

subroutine testAdjointRelation(identifier, direction, success, isPeriodic, tolerance)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use StencilOperator_mod, only : t_StencilOperator

  ! <<< Internal modules >>>
  use MPIHelper, only : pigeonhole
  use RandomNumber, only : random

  ! <<< Arguments >>>
  character(len = *), intent(in) :: identifier
  integer, intent(in) :: direction
  logical, intent(out) :: success
  logical, intent(in), optional :: isPeriodic
  real(SCALAR_KIND), intent(in), optional :: tolerance

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  logical :: isPeriodic_(3)
  real(wp) :: tolerance_
  integer :: n, offset, nLocal, gridSize(3), cartesianCommunicator,                          &
       numProcesses(3), procRank, numProcs, ierror
  type(t_StencilOperator) :: A, adjointOfA
  real(wp), allocatable :: f(:,:), g(:,:)
  SCALAR_TYPE, allocatable :: u(:,:), v(:,:), norm(:,:)
#ifdef SCALAR_TYPE_IS_binary128_IEEE754
  SCALAR_TYPE, allocatable :: mpiReduceBuffer(:)
#endif
  real(wp) :: leftHandSide, rightHandSide

  if (direction < 0 .or. direction > 3) then
     success = .false.
     return
  end if

  isPeriodic_ = .false.
  if (present(isPeriodic)) isPeriodic_(direction) = isPeriodic

  tolerance_ = epsilon(0.0_wp)
  if (present(tolerance)) tolerance_ = tolerance

  ! Find the number of processes in the communicator.
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
#ifdef SCALAR_TYPE_IS_binary128_IEEE754
  allocate(mpiReduceBuffer(numProcs))
#endif

  numProcesses = 1
  numProcesses(direction) = numProcs

  ! Create a Cartesian communicator.
  call MPI_Cart_create(MPI_COMM_WORLD, 3, numProcesses, isPeriodic_,                         &
       .true., cartesianCommunicator, ierror)

  ! Update the operators.
  call A%setup(identifier)
  call A%update(cartesianCommunicator, direction)
  call A%getAdjoint(adjointOfA)
  call adjointOfA%update(cartesianCommunicator, direction)

  ! Decide the size of arrays to be used for the test.
  n = random(numProcs * adjointOfA%boundaryDepth, 2 ** 16)
  call MPI_Bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

  ! Determine the offset and number of points that will be distributed to the current
  ! process.
  call pigeonhole(n, numProcs, procRank, offset, nLocal)
  gridSize = 1; gridSize(direction) = nLocal

  ! Allocate process-level data.
  allocate(f(nLocal, 1), g(nLocal, 1), u(nLocal, 1), v(nLocal, 1), norm(nLocal, 1))

  ! Setup the norm.
  norm = 1.0_wp
  call A%applyNorm(norm, gridSize)

  ! Initialize `f` and `g` to random values.
  call random_number(f)
  call random_number(g)

  u = f
  v = g
  call A%apply(u, gridSize)
  leftHandSide = sum(u * norm * v)
#ifdef SCALAR_TYPE_IS_binary128_IEEE754
  call MPI_Allgather(leftHandSide, 1, SCALAR_TYPE_MPI, mpiReduceBuffer,                      &
       1, SCALAR_TYPE_MPI, MPI_COMM_WORLD, ierror)
  leftHandSide = sum(mpiReduceBuffer)
#else
  call MPI_Allreduce(MPI_IN_PLACE, leftHandSide, 1, SCALAR_TYPE_MPI,                         &
       MPI_SUM, MPI_COMM_WORLD, ierror)
#endif

  u = f
  v = g
  call adjointOfA%apply(v, gridSize)
  rightHandSide = sum(u * norm * v)
#ifdef SCALAR_TYPE_IS_binary128_IEEE754
  call MPI_Allgather(rightHandSide, 1, SCALAR_TYPE_MPI, mpiReduceBuffer,                     &
       1, SCALAR_TYPE_MPI, MPI_COMM_WORLD, ierror)
  rightHandSide = sum(mpiReduceBuffer)
#else
  call MPI_Allreduce(MPI_IN_PLACE, rightHandSide, 1, SCALAR_TYPE_MPI,                        &
       MPI_SUM, MPI_COMM_WORLD, ierror)
#endif

  success = (abs(leftHandSide - rightHandSide) / real(n, wp) <= tolerance_)

  SAFE_DEALLOCATE(f)
  SAFE_DEALLOCATE(g)
  SAFE_DEALLOCATE(u)
  SAFE_DEALLOCATE(v)
  SAFE_DEALLOCATE(norm)

  call A%cleanup()
  call adjointOfA%cleanup()

end subroutine testAdjointRelation
