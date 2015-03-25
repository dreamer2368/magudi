#include "config.h"

program operator_convergence

  use MPI

  use StencilOperator_type, only : t_StencilOperator

  use StencilOperator_mod, only : setupOperator, cleanupOperator

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  type(t_StencilOperator) :: A
  integer :: direction, ierror
  logical :: success, success_

  interface

     function F1(x) result(y)
       SCALAR_TYPE, intent(in) :: x
       SCALAR_TYPE :: y
     end function F1

     function F2(x) result(y)
       SCALAR_TYPE, intent(in) :: x
       SCALAR_TYPE :: y
     end function F2

     function F3(x) result(y)
       SCALAR_TYPE, intent(in) :: x
       SCALAR_TYPE :: y
     end function F3

     function dF1(x) result(y)
       SCALAR_TYPE, intent(in) :: x
       SCALAR_TYPE :: y
     end function dF1

     function dF2(x) result(y)
       SCALAR_TYPE, intent(in) :: x
       SCALAR_TYPE :: y
     end function dF2

     function dF3(x) result(y)
       SCALAR_TYPE, intent(in) :: x
       SCALAR_TYPE :: y
     end function dF3

     subroutine testStencilOperatorConvergence(A, direction, f, g, convergenceRate, success, &
          isPeriodic, quotientExponent, startSize, nIterations, refinementFactor)

       use StencilOperator_type, only : t_StencilOperator

       type(t_StencilOperator) :: A
       integer, intent(in) :: direction

       interface

          function f(x)
            SCALAR_TYPE, intent(in) :: x
            SCALAR_TYPE :: f
          end function f

          function g(x)
            SCALAR_TYPE, intent(in) :: x
            SCALAR_TYPE :: g
          end function g

       end interface

       integer, intent(in) :: convergenceRate
       logical, intent(inout) :: success

       logical, intent(in), optional :: isPeriodic
       integer, intent(in), optional :: quotientExponent, startSize, nIterations
       real(SCALAR_KIND), intent(in), optional :: refinementFactor

     end subroutine testStencilOperatorConvergence

  end interface

  call MPI_Init(ierror)

  success = .true.

  call setupOperator(A, "SBP 1-2 first derivative", success_)
  success = success .and. success_
  if (success_) then
     do direction = 1, 3
        call testStencilOperatorConvergence(A, direction, F1, dF1, 2,                        &
             success, isPeriodic = .true., quotientExponent = 1)
        call testStencilOperatorConvergence(A, direction, F2, dF2, 1,                        &
             success, quotientExponent = 1)
        call testStencilOperatorConvergence(A, direction, F3, dF3, 1,                        &
             success, quotientExponent = 1)
     end do
  end if

  call setupOperator(A, "SBP 2-4 first derivative", success_)
  success = success .and. success_
  if (success_) then
     do direction = 1, 3
        call testStencilOperatorConvergence(A, direction, F1, dF1, 4,                        &
             success, isPeriodic = .true., quotientExponent = 1)
        call testStencilOperatorConvergence(A, direction, F2, dF2, 2,                        &
             success, quotientExponent = 1)
        call testStencilOperatorConvergence(A, direction, F3, dF3, 2,                        &
             success, quotientExponent = 1)
     end do
  end if

  call setupOperator(A, "SBP 3-6 first derivative", success_)
  success = success .and. success_
  if (success_) then
     do direction = 1, 3
        call testStencilOperatorConvergence(A, direction, F1, dF1, 6,                        &
             success, isPeriodic = .true., quotientExponent = 1)
        call testStencilOperatorConvergence(A, direction, F2, dF2, 3,                        &
             success, quotientExponent = 1)
        call testStencilOperatorConvergence(A, direction, F3, dF3, 3,                        &
             success, quotientExponent = 1)
     end do
  end if

  call setupOperator(A, "SBP 4-8 first derivative", success_)
  success = success .and. success_
  if (success_) then
     do direction = 1, 3
        call testStencilOperatorConvergence(A, direction, F1, dF1, 8, success,               &
             isPeriodic = .true., quotientExponent = 1,                                      &
             refinementFactor = 1.06_wp)
        call testStencilOperatorConvergence(A, direction, F2, dF2, 4, success,               &
             quotientExponent = 1, refinementFactor = 1.08_wp)
        call testStencilOperatorConvergence(A, direction, F3, dF3, 4, success,               &
             quotientExponent = 1, refinementFactor = 1.08_wp)
     end do
  end if

  call cleanupOperator(A)

  call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL,                                  &
       MPI_LAND, MPI_COMM_WORLD, ierror)
  call MPI_Finalize(ierror)
  if (.not. success) stop -1
  stop 0

end program operator_convergence

subroutine testStencilOperatorConvergence(A, direction, f, g, convergenceRate, success,      &
     isPeriodic, quotientExponent, startSize, nIterations, refinementFactor)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use StencilOperator_type, only : t_StencilOperator

  ! <<< Internal modules >>>
  use MPIHelper, only : pigeonhole
  use StencilOperator_mod, only : updateOperator, applyOperator

  ! <<< Arguments >>>
  type(t_StencilOperator) :: A
  integer, intent(in) :: direction
  interface
     function f(x)
       SCALAR_TYPE, intent(in) :: x
       SCALAR_TYPE :: f
     end function f
     function g(x)
       SCALAR_TYPE, intent(in) :: x
       SCALAR_TYPE :: g
     end function g
  end interface
  integer, intent(in) :: convergenceRate
  logical, intent(inout) :: success
  logical, intent(in), optional :: isPeriodic
  integer, intent(in), optional :: quotientExponent, startSize, nIterations
  real(SCALAR_KIND), intent(in), optional :: refinementFactor

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, n, startSize_, nLocal, offset, gridSize(3), nIterations_,                 &
       proc, nProcs, numProcesses(3), cartesianCommunicator, ierror
  logical :: isPeriodic_(3)
  SCALAR_TYPE :: x
  SCALAR_TYPE, allocatable :: y(:,:), yExact(:,:)
  real(wp) :: h, h_
  real(wp), allocatable :: errorHistory(:), convergenceHistory(:)

  interface

     real(SCALAR_KIND) function meanTrimmed(a)
       real(SCALAR_KIND), intent(inout) :: a(:)
     end function meanTrimmed

     subroutine sort(a)
       real(SCALAR_KIND), intent(inout) :: a(:)
     end subroutine sort

  end interface

  if (direction < 0 .or. direction > 3) then
     success = .false.
     return
  end if

  isPeriodic_ = .false.
  if (present(isPeriodic)) isPeriodic_(direction) = isPeriodic

  nIterations_ = 40
  if (present(nIterations)) nIterations_ = nIterations

  SAFE_DEALLOCATE(errorHistory)
  SAFE_DEALLOCATE(convergenceHistory)
  allocate(errorHistory(nIterations_), convergenceHistory(nIterations_ - 1))

  ! Find the rank of this process and the number of processes in the communicator.
  call MPI_Comm_rank(MPI_COMM_WORLD, proc, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, nProcs, ierror)

  numProcesses = 1
  numProcesses(direction) = nProcs

  ! Create a Cartesian communicator.
  call MPI_Cart_create(MPI_COMM_WORLD, 3, numProcesses,                                      &
       isPeriodic_, .true., cartesianCommunicator, ierror)
  call updateOperator(A, cartesianCommunicator, direction)

  if (present(startSize)) then
     startSize_ = startSize
  else
     startSize_ = max(nProcs * A%boundaryDepth, 32)
  end if
  n = startSize_

  do i = 1, nIterations_

     ! Determine the offset and number of points that will be distributed to the current
     ! process.
     call pigeonhole(n, nProcs, proc, offset, nLocal)
     gridSize(direction) = nLocal

     ! If periodic, then exclude the x = 1 plane.
     if (isPeriodic_(direction)) then
        h = 1.0_wp / real(n, wp)
     else
        h = 1.0_wp / real(n - 1, wp)
     end if

     ! Allocate process-level data.
     SAFE_DEALLOCATE(y)
     SAFE_DEALLOCATE(yExact)
     allocate(y(nLocal, 1), yExact(nLocal, 1))

     ! Evaluate the functions.
     do j = offset + 1, offset + nLocal
        x = real(j - 1, wp) * h
        y(j - offset, 1) = f(x) ; yExact(j - offset, 1) = g(x)
     end do

     ! Apply the stencil.
     gridSize = 1 ; gridSize(direction) = nLocal
     call applyOperator(A, y, gridSize)
     if (present(quotientExponent)) y = y / h ** quotientExponent

     ! Compute the error.
     errorHistory(i) = maxval(abs(y - yExact))
     call MPI_Allreduce(MPI_IN_PLACE, errorHistory(i), 1, SCALAR_TYPE_MPI,                   &
          MPI_MAX, MPI_COMM_WORLD, ierror)
     if (i > 1) then
        convergenceHistory(i - 1) = log(errorHistory(i) / errorHistory(i - 1)) / log(h / h_)
        if (convergenceHistory(i - 1) < 0.0_wp) exit
     end if
     h_ = h

     SAFE_DEALLOCATE(y)
     SAFE_DEALLOCATE(yExact)

     ! Increase the number of points for the next iteration.
     if (present(refinementFactor)) then
        n = nint(real(n, wp) * refinementFactor)
     else
        n = nint(real(n, wp) * min(1.2_wp, exp(log(2.0_wp ** 20 / startSize_) /              &
             real(nIterations_, wp))))
     end if

  end do

  i = min(i, nIterations_)
  if (i <= 2) then
     success = .false.
  else
     call sort(convergenceHistory(1:i-2))
     if (proc == 0) then
        n = nint(meanTrimmed(convergenceHistory(1:i-2)))
     end if
     success = success .and.                                                                 &
          (nint(meanTrimmed(convergenceHistory(1:i-2))) >= convergenceRate)
  end if

end subroutine testStencilOperatorConvergence

subroutine sort(a)

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(inout) :: a(:)

  ! <<< Local variables >>>
  integer :: i, j
  real(SCALAR_KIND) :: temp

  do i = 2, size(a)
     j = i - 1
     temp = a(i)
     do while (a(j) > temp)
        a(j+1) = a(j)
        j = j - 1
        if (j < 1) exit
     end do
     a(j+1) = temp
  end do

end subroutine sort

real(SCALAR_KIND) function median(a)

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(in) :: a(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: n

  n = size(a)
  if (mod(n, 2) == 1) then
     median = a((n + 1) / 2)
  else
     median = 0.5_wp * (a(n / 2) + a(n / 2 + 1))
  end if

end function median

real(SCALAR_KIND) function meanTrimmed(a)

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(inout) :: a(:)

  ! <<< Scalar variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, n
  real(wp) :: firstQuartile, thirdQuartile

  interface
     real(SCALAR_KIND) function median(a)
       real(SCALAR_KIND), intent(in) :: a(:)
     end function median
  end interface

  n = size(a)

  if (mod(n, 2) == 0) then
     firstQuartile = median(a(1:n/2))
     thirdQuartile = median(a(n/2+1:n))
  else
     firstQuartile = median(a(1:(n-1)/2))
     thirdQuartile = median(a((n+1)/2+1:n))
  end if

  meanTrimmed = 0.0_wp
  n = 0

  do i = 1, size(a)
     if (a(i) >= firstQuartile .and. a(i) <= thirdQuartile) then
        meanTrimmed = meanTrimmed + a(i)
        n = n + 1
     end if
  end do

  if (n == 0) return
  meanTrimmed = meanTrimmed / real(n, wp)

end function meanTrimmed

function F1(x) result(y)
  SCALAR_TYPE, intent(in) :: x
  SCALAR_TYPE :: y
  integer, parameter :: wp = SCALAR_KIND
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  y = sin(2.0_wp * pi * x)
end function F1

function dF1(x) result(y)
  SCALAR_TYPE, intent(in) :: x
  SCALAR_TYPE :: y
  integer, parameter :: wp = SCALAR_KIND
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  y = 2.0_wp * pi * cos(2.0_wp * pi * x)
end function dF1

function F2(x) result(y)
  SCALAR_TYPE, intent(in) :: x
  SCALAR_TYPE :: y
  integer, parameter :: wp = SCALAR_KIND
  y = sin(x + 3.0_wp) / (x + 3.0_wp)
end function F2

function dF2(x) result(y)
  SCALAR_TYPE, intent(in) :: x
  SCALAR_TYPE :: y
  integer, parameter :: wp = SCALAR_KIND
  y = cos(x + 3.0_wp) / (x + 3.0_wp) - sin(x + 3.0_wp) / (x + 3.0_wp) ** 2
end function dF2

function F3(x) result(y)
  SCALAR_TYPE, intent(in) :: x
  SCALAR_TYPE :: y
  integer, parameter :: wp = SCALAR_KIND
  y = tanh(4.0_wp * (real(x, wp) - 0.5_wp))
end function F3

function dF3(x) result(y)
  SCALAR_TYPE, intent(in) :: x
  SCALAR_TYPE :: y
  integer, parameter :: wp = SCALAR_KIND
  y = 4.0_wp * (1.0_wp - tanh(4.0_wp * (real(x, wp) - 0.5_wp)) ** 2)
end function dF3
