#include "config.h"

module StencilOperator_type

  use MPI, only : MPI_COMM_NULL

  implicit none
  private

  integer, parameter, public ::                                                              &
       SYMMETRIC      = 0,                                                                   &
       SKEW_SYMMETRIC = 1

  type, public :: t_StencilOperator

     integer :: symmetryType, interiorWidth, boundaryWidth, boundaryDepth
     real(SCALAR_KIND), allocatable :: normBoundary(:), rhsInterior(:),                      &
          rhsBoundary1(:,:), rhsBoundary2(:,:)
     logical :: hasDomainBoundary(2)
     integer :: cartesianCommunicator = MPI_COMM_NULL, direction, nGhost(2), periodicOffset(2)

  end type t_StencilOperator

end module StencilOperator_type

module StencilOperator_mod

  implicit none

  interface

     subroutine setupOperator(this, identifier, success)

       !> Sets up a stencil operator.

       use StencilOperator_type

       type(t_StencilOperator) :: this
       character(len = *), intent(in) :: identifier

       logical, intent(out), optional :: success

     end subroutine setupOperator

  end interface

  interface

     subroutine updateOperator(this, cartesianCommunicator, &
          direction, isPeriodicityOverlapping)

       !> Attaches a Cartesian communicator to a stencil operator and specifies the direction
       !> along which it will be applied using `applyOperator`.

       use StencilOperator_type

       type(t_StencilOperator) :: this
       integer, intent(in) :: cartesianCommunicator, direction

       logical, intent(in), optional :: isPeriodicityOverlapping

     end subroutine updateOperator

  end interface

  interface

     subroutine cleanupOperator(this)

       !> Deallocates stencil operator data.

       use StencilOperator_type

       type(t_StencilOperator) :: this

     end subroutine cleanupOperator

  end interface

  interface

     subroutine getAdjointOperator(this, adjointOperator)

       use StencilOperator_type

       type(t_StencilOperator) :: this, adjointOperator

     end subroutine getAdjointOperator

  end interface

  interface

     subroutine applyOperator(this, x, gridSize)

       !> Applies a stencil operator to a real/complex semidiscrete vector.

       use StencilOperator_type

       type(t_StencilOperator) :: this
       SCALAR_TYPE, intent(inout) :: x(:,:)
       integer, intent(in) :: gridSize(3)

     end subroutine applyOperator

  end interface

  interface

     subroutine applyOperatorAtInteriorPoints(this, xWithGhostPoints, x, gridSize)

       !> Applies a stencil operator to a real/complex semidiscrete vector at interior points
       !> only. This is equivalent to calling `applyOperator` if the direction is
       !> periodic. However, note the difference in call signatures.

       use StencilOperator_type

       type(t_StencilOperator) :: this
       SCALAR_TYPE, intent(in) :: xWithGhostPoints(:,:,:,:)
       SCALAR_TYPE, intent(inout) :: x(:,:)
       integer, intent(in) :: gridSize(3)

     end subroutine applyOperatorAtInteriorPoints

  end interface

  interface

     subroutine applyOperatorNorm(this, x, gridSize)

       use StencilOperator_type

       !> Multiplies a real/complex semidiscrete vector by the diagonal norm matrix in the
       !> direction along which the stencil operator is applied.

       type(t_StencilOperator) :: this
       SCALAR_TYPE, intent(inout) :: x(:,:)
       integer, intent(in) :: gridSize(3)

     end subroutine applyOperatorNorm

  end interface

  interface

     subroutine applyOperatorNormInverse(this, x, gridSize)

       use StencilOperator_type

       !> Multiplies a real/complex semidiscrete vector by the inverse of the diagonal norm
       !> matrix in the direction along which the stencil operator is applied.

       type(t_StencilOperator) :: this
       SCALAR_TYPE, intent(inout) :: x(:,:)
       integer, intent(in) :: gridSize(3)

     end subroutine applyOperatorNormInverse

  end interface

end module StencilOperator_mod
