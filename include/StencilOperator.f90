#include "config.h"

module StencilOperator_mod

  use MPI, only : MPI_COMM_NULL

  implicit none

  type, public :: t_StencilOperator

     integer :: cartesianCommunicator = MPI_COMM_NULL, direction, symmetryType,              &
          interiorWidth, boundaryWidth, boundaryDepth, nGhost(2), periodicOffset(2)
     logical :: hasDomainBoundary(2)
     real(SCALAR_KIND), allocatable :: normBoundary(:), rhsInterior(:),                      &
          rhsBoundary1(:,:), rhsBoundary2(:,:)

   contains

     procedure, pass :: setup => setupOperator
     procedure, pass :: update => updateOperator
     procedure, pass :: cleanup => cleanupOperator
     procedure, pass :: getAdjoint => getAdjointOperator
     procedure, pass :: apply => applyOperator
     procedure, pass :: applyAtInteriorPoints => applyOperatorAtInteriorPoints
     procedure, pass :: applyNorm => applyOperatorNorm
     procedure, pass :: applyNormInverse => applyOperatorNormInverse
     procedure, pass :: applyAndProjectOnBoundary => applyOperatorAndProjectOnBoundary
     procedure, pass :: projectOnBoundaryAndApply => projectOnBoundaryAndApplyOperator

  end type t_StencilOperator

  interface

     subroutine setupOperator(this, stencilScheme)

       !> Sets up a stencil operator.

       import :: t_StencilOperator

       class(t_StencilOperator) :: this
       character(len = *), intent(in) :: stencilScheme

     end subroutine setupOperator

  end interface

  interface

     subroutine updateOperator(this, cartesianCommunicator,                                  &
          direction, isPeriodicityOverlapping)

       !> Attaches a Cartesian communicator to a stencil operator and specifies the direction
       !> along which it will be applied using `applyOperator`.

       import :: t_StencilOperator

       class(t_StencilOperator) :: this
       integer, intent(in) :: cartesianCommunicator, direction

       logical, intent(in), optional :: isPeriodicityOverlapping

     end subroutine updateOperator

  end interface

  interface

     subroutine cleanupOperator(this)

       !> Deallocates stencil operator data.

       import :: t_StencilOperator

       class(t_StencilOperator) :: this

     end subroutine cleanupOperator

  end interface

  interface

     subroutine getAdjointOperator(this, adjointOperator)

       import :: t_StencilOperator

       class(t_StencilOperator), intent(in) :: this
       class(t_StencilOperator) :: adjointOperator

     end subroutine getAdjointOperator

  end interface

  interface

     subroutine applyOperator(this, x, gridSize)

       !> Applies a stencil operator to a real/complex semidiscrete vector.

       import :: t_StencilOperator

       class(t_StencilOperator), intent(in) :: this
       SCALAR_TYPE, intent(inout) :: x(:,:)
       integer, intent(in) :: gridSize(3)

     end subroutine applyOperator

  end interface

  interface

     pure subroutine applyOperatorAtInteriorPoints(this, xWithGhostPoints, x, gridSize)

       !> Applies a stencil operator to a real/complex semidiscrete vector at interior points
       !> only. This is equivalent to calling `applyOperator` if the direction is
       !> periodic. However, note the difference in call signatures.

       import :: t_StencilOperator

       class(t_StencilOperator), intent(in) :: this
       SCALAR_TYPE, intent(in) :: xWithGhostPoints(:,:,:,:)
       SCALAR_TYPE, intent(inout) :: x(:,:)
       integer, intent(in) :: gridSize(3)

     end subroutine applyOperatorAtInteriorPoints

  end interface

  interface

     pure subroutine applyOperatorAndProjectOnBoundary(this, x, gridSize, faceOrientation)

       !> Applies a stencil operator to a real/complex semidiscrete vector followed by zeroing
       !> the result at all points except points that lie on the left (right) boundary of the
       !> computational domain if `faceOrientation` is greater (lesser) than zero.

       import :: t_StencilOperator

       class(t_StencilOperator), intent(in) :: this
       SCALAR_TYPE, intent(inout) :: x(:,:)
       integer, intent(in) :: gridSize(3), faceOrientation

     end subroutine applyOperatorAndProjectOnBoundary

  end interface

  interface

     pure subroutine projectOnBoundaryAndApplyOperator(this, x, gridSize, faceOrientation)

       !> Zeros a real/complex semidiscrete vector at all points except points that lie on the
       !> left (right) boundary of the computational domain if `faceOrientation` is greater
       !> (lesser) than zero. Then, applies a stencil operator to the result.

       import :: t_StencilOperator

       class(t_StencilOperator), intent(in) :: this
       SCALAR_TYPE, intent(inout) :: x(:,:)
       integer, intent(in) :: gridSize(3), faceOrientation

     end subroutine projectOnBoundaryAndApplyOperator

  end interface

  interface

     pure subroutine applyOperatorNorm(this, x, gridSize)

       import :: t_StencilOperator

       !> Multiplies a real/complex semidiscrete vector by the diagonal norm matrix in the
       !> direction along which the stencil operator is applied.

       class(t_StencilOperator), intent(in) :: this
       SCALAR_TYPE, intent(inout) :: x(:,:)
       integer, intent(in) :: gridSize(3)

     end subroutine applyOperatorNorm

  end interface

  interface

     pure subroutine applyOperatorNormInverse(this, x, gridSize)

       import :: t_StencilOperator

       !> Multiplies a real/complex semidiscrete vector by the inverse of the diagonal norm
       !> matrix in the direction along which the stencil operator is applied.

       class(t_StencilOperator), intent(in) :: this
       SCALAR_TYPE, intent(inout) :: x(:,:)
       integer, intent(in) :: gridSize(3)

     end subroutine applyOperatorNormInverse

  end interface

end module StencilOperator_mod
