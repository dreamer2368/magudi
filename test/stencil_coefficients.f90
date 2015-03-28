#include "config.h"

program stencil_coefficients

  use StencilOperator_type, only : t_StencilOperator

  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use StencilOperator_mod, only : setupOperator

  implicit none

  !> Essentially, a check for typos in stencil coefficients initialized in
  !> StencilOperatorImpl.f90. Tests whether the coefficients satisfy the rules expected of
  !> them in order for the operator to converge to a certain order. Doesn't actually apply the
  !> stencil to an array.

  integer, parameter :: wp = SCALAR_KIND
  type(t_StencilOperator) :: A
  logical :: success, success_

  interface

     subroutine testStencilAccuracy(D, derivativeOrder, interiorOrderOfAccuracy,             &
          boundaryOrderOfAccuracy, success, tolerance)

       use StencilOperator_type, only : t_StencilOperator

       type(t_StencilOperator) :: D
       integer, intent(in) :: derivativeOrder, interiorOrderOfAccuracy,                      &
            boundaryOrderOfAccuracy
       logical, intent(inout) :: success

       real(SCALAR_KIND), intent(in), optional :: tolerance

     end subroutine testStencilAccuracy

  end interface

  success = .true.

  call initializeErrorHandler()

  call setupOperator(A, "SBP 1-2 first derivative", success_)
  success = success .and. success_
  if (success_) call testStencilAccuracy(A, 1, 2, 1, success)
  call setupOperator(A, "SBP 1-2 second derivative", success_)
  success = success .and. success_
  if (success_) call testStencilAccuracy(A, 2, 2, 1, success)
  call setupOperator(A, "SBP 1-2 dissipation", success_)
  success = success .and. success_
  if (success_) call testStencilAccuracy(A, 0, 2, 1, success)

  call setupOperator(A, "SBP 2-4 first derivative", success_)
  success = success .and. success_
  if (success_) call testStencilAccuracy(A, 1, 4, 2, success)
  call setupOperator(A, "SBP 2-4 second derivative", success_)
  success = success .and. success_
  if (success_) call testStencilAccuracy(A, 2, 4, 2, success)
  call setupOperator(A, "SBP 2-4 dissipation", success_)
  success = success .and. success_
  if (success_) call testStencilAccuracy(A, 0, 4, 2, success)

  call setupOperator(A, "SBP 3-6 first derivative", success_)
  success = success .and. success_
  if (success_) call testStencilAccuracy(A, 1, 6, 3, success)
  call setupOperator(A, "SBP 3-6 second derivative", success_)
  success = success .and. success_
  if (success_) call testStencilAccuracy(A, 2, 6, 3, success,                                &
       tolerance = epsilon(0.0_wp) * 2.0_wp) !... some schemes require a higher tolerance.
  call setupOperator(A, "SBP 3-6 dissipation", success_)
  success = success .and. success_
  if (success_) call testStencilAccuracy(A, 0, 6, 3, success)

  call setupOperator(A, "SBP 4-8 first derivative", success_)
  success = success .and. success_
  if (success_) call testStencilAccuracy(A, 1, 8, 4, success,                                &
       tolerance = epsilon(0.0_wp) * 300.0_wp)
  call setupOperator(A, "SBP 4-8 second derivative", success_)
  success = success .and. success_
  if (success_) call testStencilAccuracy(A, 2, 8, 4, success,                                &
       tolerance = epsilon(0.0_wp) * 4.0_wp)
  call setupOperator(A, "SBP 4-8 dissipation", success_)
  success = success .and. success_
  if (success_) call testStencilAccuracy(A, 0, 8, 4, success,                                &
       tolerance = epsilon(0.0_wp) * 20.0_wp)

  call cleanupOperator(A)

  call cleanupErrorHandler()

  if (.not. success) stop -1
  stop 0

end program stencil_coefficients

subroutine testStencilAccuracy(D, derivativeOrder, interiorOrderOfAccuracy,                  &
     boundaryOrderOfAccuracy, success, tolerance)

  ! <<< Internal modules >>>
  use StencilOperator_type, only : t_StencilOperator

  ! <<< Arguments >>>
  type(t_StencilOperator) :: D
  integer, intent(in) :: derivativeOrder, interiorOrderOfAccuracy, boundaryOrderOfAccuracy
  logical, intent(inout) :: success
  real(SCALAR_KIND), intent(in), optional :: tolerance

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i
  real(wp) :: tolerance_
  real(wp), allocatable :: boundary_stencil(:)

  interface

     function getOrderOfAccuracy(stencil, derivativeOrder, tolerance) result(orderOfAccuracy)

       real(SCALAR_KIND), allocatable, intent(in) :: stencil(:)
       integer, intent(in) :: derivativeOrder
       real(SCALAR_KIND), intent(in) :: tolerance

       integer :: orderOfAccuracy

     end function getOrderOfAccuracy

  end interface

  ! By default, use machine epsilon as tolerance. This works as long as the coefficients
  ! are typed in as fractions. If they are not available as fractions, and they are
  ! typed in decimal form, use a higher tolerance.
  tolerance_ = epsilon(0.0_wp)
  if (present(tolerance)) tolerance_ = tolerance

  success = success .and.                                                                    &
       (getOrderOfAccuracy(D%RHSinterior, derivativeOrder, tolerance_) ==                    &
       interiorOrderOfAccuracy) !... interior
  do i = 1, D%boundaryDepth
     SAFE_DEALLOCATE(boundary_stencil)
     allocate(boundary_stencil(1-i:D%boundaryWidth-i),                                       &
          source = D%RHSboundary1(:,i)) !... copy the boundary stencils.
     success = success .and.                                                                 &
          (getOrderOfAccuracy(boundary_stencil, derivativeOrder, tolerance_) >=              &
          boundaryOrderOfAccuracy)
  end do
  SAFE_DEALLOCATE(boundary_stencil)

end subroutine testStencilAccuracy

function getOrderOfAccuracy(stencil, derivativeOrder, tolerance) result(orderOfAccuracy)

  ! <<< Arguments >>>
  real(SCALAR_KIND), allocatable, intent(in) :: stencil(:)
  integer, intent(in) :: derivativeOrder
  real(SCALAR_KIND), intent(in) :: tolerance

  ! <<< Result >>>
  integer :: orderOfAccuracy

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i
  real(wp) :: error

  interface

     function factorial(n)
       integer, intent(in) :: n
       integer :: factorial
     end function factorial

  end interface

  orderOfAccuracy = 0

  ! Form the Taylor series and check that the error is below tolerance level.
  do while (orderOfAccuracy < 50)
     error = 0.0_wp
     if (derivativeOrder > 0 .and. orderOfAccuracy == derivativeOrder) error = -1.0_wp
     do i = lbound(stencil, 1), ubound(stencil, 1)
        if (orderOfAccuracy == 0) then
           error = error + stencil(i)
        else
           error = error + stencil(i) * real(i, wp) ** orderOfAccuracy /                     &
                real(factorial(orderOfAccuracy), wp)
        end if
     end do
     if (abs(error) / real(size(stencil), wp) > tolerance) exit
     orderOfAccuracy = orderOfAccuracy + 1
  end do

  if (orderOfAccuracy == 50) then
     orderOfAccuracy = huge(orderOfAccuracy)
  else
     orderOfAccuracy = orderOfAccuracy - derivativeOrder
  end if

end function getOrderOfAccuracy

function factorial(n)
  integer, intent(in) :: n
  integer :: factorial, i
  factorial = 1
  do i = 2, n
     factorial = factorial * i
  end do
end function factorial
