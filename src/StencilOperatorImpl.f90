#include "config.h"

module StencilOperatorImpl

  implicit none
  public

contains

  subroutine allocateData(this)

    ! <<< Derived types >>>
    use StencilOperator_type

    ! <<< Arguments >>>
    type(t_StencilOperator) :: this

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND

    allocate(this%rhsInterior(-this%interiorWidth/2:this%interiorWidth/2), source = 0.0_wp)
    allocate(this%rhsBoundary1(this%boundaryWidth, this%boundaryDepth), source = 0.0_wp)
    allocate(this%rhsBoundary2(this%boundaryWidth, this%boundaryDepth), source = 0.0_wp)
    allocate(this%normBoundary(this%boundaryDepth), source = 1.0_wp)

  end subroutine allocateData

  subroutine applyOperator_1(this, x, gridSize)

    ! <<< Derived types >>>
    use StencilOperator_type

    ! <<< Internal modules >>>
    use MPIHelper, only : fillGhostPoints

    ! <<< Arguments >>>
    type(t_StencilOperator) :: this
    SCALAR_TYPE, intent(inout) :: x(:,:)
    integer, intent(in) :: gridSize(3)

    ! <<< Local variables >>>
    integer :: i, j, k, l, m, n
    SCALAR_TYPE, allocatable :: xWithGhostPoints(:,:,:,:)

    allocate(xWithGhostPoints(gridSize(1) + sum(this%nGhost),                                &
         gridSize(2), gridSize(3), size(x, 2)))

    do l = 1, size(x, 2)
       do k = 1, gridSize(3)
          do j = 1, gridSize(2)
             do i = 1, gridSize(1)
                xWithGhostPoints(i + this%nGhost(1), j, k, l) =                              &
                     x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l)
             end do
          end do
       end do
    end do

    call fillGhostPoints(this%cartesianCommunicator, xWithGhostPoints,                       &
         1, this%nGhost, this%periodicOffset)
    call applyOperatorAtInteriorPoints_1(this, xWithGhostPoints, x, gridSize)

    n = this%boundaryWidth

    ! Left boundary.
    if (this%hasDomainBoundary(1)) then
       i = this%nGhost(1)
       do l = 1, size(x, 2)
          do k = 1, gridSize(3)
             do j = 1, gridSize(2)
                do m = 1, this%boundaryDepth
                   x(i + m - this%nGhost(1) +                                                &
                        gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) =                  &
                        sum(this%rhsBoundary1(:,m) * xWithGhostPoints(i+1:i+n,j,k,l))
                end do
             end do
          end do
       end do
    end if

    ! Right boundary.
    if (this%hasDomainBoundary(2)) then
       i = this%nGhost(1) + gridSize(1) + 1
       do l = 1, size(x, 2)
          do k = 1, gridSize(3)
             do j = 1, gridSize(2)
                do m = this%boundaryDepth, 1, -1
                   x(i - m - this%nGhost(1) +                                                &
                        gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) =                  &
                        sum(this%rhsBoundary2(:,m) * xWithGhostPoints(i-n:i-1,j,k,l))
                end do
             end do
          end do
       end do
    end if

    SAFE_DEALLOCATE(xWithGhostPoints)

  end subroutine applyOperator_1

  subroutine applyOperator_2(this, x, gridSize)

    ! <<< Derived types >>>
    use StencilOperator_type

    ! <<< Internal modules >>>
    use MPIHelper, only : fillGhostPoints

    ! <<< Arguments >>>
    type(t_StencilOperator) :: this
    SCALAR_TYPE, intent(inout) :: x(:,:)
    integer, intent(in) :: gridSize(3)

    ! <<< Local variables >>>
    integer :: i, j, k, l, m, n
    SCALAR_TYPE, allocatable :: xWithGhostPoints(:,:,:,:)

    allocate(xWithGhostPoints(gridSize(1), gridSize(2) +                                     &
         sum(this%nGhost), gridSize(3), size(x, 2)))

    do l = 1, size(x, 2)
       do k = 1, gridSize(3)
          do j = 1, gridSize(2)
             do i = 1, gridSize(1)
                xWithGhostPoints(i, j + this%nGhost(1), k, l) =                              &
                     x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l)
             end do
          end do
       end do
    end do

    call fillGhostPoints(this%cartesianCommunicator, xWithGhostPoints,                       &
         2, this%nGhost, this%periodicOffset)
    call applyOperatorAtInteriorPoints_2(this, xWithGhostPoints, x, gridSize)

    n = this%boundaryWidth

    ! Left boundary.
    if (this%hasDomainBoundary(1)) then
       j = this%nGhost(1)
       do l = 1, size(x, 2)
          do k = 1, gridSize(3)
             do i = 1, gridSize(1)
                do m = 1, this%boundaryDepth
                   x(i + gridSize(1) * (j - 1 + m - this%nGhost(1) +                         &
                        gridSize(2) * (k - 1)), l) =                                         &
                        sum(this%rhsBoundary1(:,m) * xWithGhostPoints(i,j+1:j+n,k,l))
                end do
             end do
          end do
       end do
    end if

    ! Right boundary.
    if (this%hasDomainBoundary(2)) then
       j = this%nGhost(1) + gridSize(2) + 1
       do l = 1, size(x, 2)
          do k = 1, gridSize(3)
             do i = 1, gridSize(1)
                do m = this%boundaryDepth, 1, -1
                   x(i + gridSize(1) * (j - 1 - m - this%nGhost(1) +                         &
                        gridSize(2) * (k - 1)), l) =                                         &
                        sum(this%rhsBoundary2(:,m) * xWithGhostPoints(i,j-n:j-1,k,l))
                end do
             end do
          end do
       end do
    end if

    SAFE_DEALLOCATE(xWithGhostPoints)

  end subroutine applyOperator_2

  subroutine applyOperator_3(this, x, gridSize)

    ! <<< Derived types >>>
    use StencilOperator_type

    ! <<< Internal modules >>>
    use MPIHelper, only : fillGhostPoints

    ! <<< Arguments >>>
    type(t_StencilOperator) :: this
    SCALAR_TYPE, intent(inout) :: x(:,:)
    integer, intent(in) :: gridSize(3)

    ! <<< Local variables >>>
    integer :: i, j, k, l, m, n
    SCALAR_TYPE, allocatable :: xWithGhostPoints(:,:,:,:)

    allocate(xWithGhostPoints(gridSize(1), gridSize(2),                                      &
         gridSize(3) + sum(this%nGhost), size(x, 2)))

    do l = 1, size(x, 2)
       do k = 1, gridSize(3)
          do j = 1, gridSize(2)
             do i = 1, gridSize(1)
                xWithGhostPoints(i, j, k + this%nGhost(1), l) =                              &
                     x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l)
             end do
          end do
       end do
    end do

    call fillGhostPoints(this%cartesianCommunicator, xWithGhostPoints,                       &
         3, this%nGhost, this%periodicOffset)
    call applyOperatorAtInteriorPoints_3(this, xWithGhostPoints, x, gridSize)

    n = this%boundaryWidth

    ! Left boundary.
    if (this%hasDomainBoundary(1)) then
       k = this%nGhost(1)
       do l = 1, size(x, 2)
          do j = 1, gridSize(2)
             do i = 1, gridSize(1)
                do m = 1, this%boundaryDepth
                   x(i + gridSize(1) * (j - 1 +                                              &
                        gridSize(2) * (k - 1 + m - this%nGhost(1))), l) =                    &
                        sum(this%rhsBoundary1(:,m) * xWithGhostPoints(i,j,k+1:k+n,l))
                end do
             end do
          end do
       end do
    end if

    ! Right boundary.
    if (this%hasDomainBoundary(2)) then
       k = this%nGhost(1) + gridSize(3) + 1
       do l = 1, size(x, 2)
          do j = 1, gridSize(2)
             do i = 1, gridSize(1)
                do m = this%boundaryDepth, 1, -1
                   x(i + gridSize(1) * (j - 1 +                                              &
                        gridSize(2) * (k - 1 - m - this%nGhost(1))), l) =                    &
                        sum(this%rhsBoundary2(:,m) * xWithGhostPoints(i,j,k-n:k-1,l))
                end do
             end do
          end do
       end do
    end if

    SAFE_DEALLOCATE(xWithGhostPoints)

  end subroutine applyOperator_3

  subroutine applyOperatorAtInteriorPoints_1(this, xWithGhostPoints, x, gridSize)

    ! <<< Derived types >>>
    use StencilOperator_type

    ! <<< Arguments >>>
    type(t_stencilOperator) :: this
    SCALAR_TYPE, intent(in) :: xWithGhostPoints(:,:,:,:)
    SCALAR_TYPE, intent(out) :: x(:,:)
    integer, intent(in) :: gridSize(3)

    ! <<< Local variables >>>
    integer :: i, j, k, l, n, is, ie

    n = this%interiorWidth / 2
    is = 1 + this%nGhost(1)
    ie = gridSize(1) + this%nGhost(1)
    if (this%nGhost(1) == 0) is = is + this%boundaryDepth
    if (this%nGhost(2) == 0) ie = ie - this%boundaryDepth

    select case (this%symmetryType) !... save FLOPS based on symmetry of interior stencil.

    case (SKEW_SYMMETRIC)
       do l = 1, size(x, 2)
          do k = 1, gridSize(3)
             do j = 1, gridSize(2)
                do i = is, ie
                   x(i - this%nGhost(1) + gridSize(1) * (j - 1 +                             &
                        gridSize(2) * (k - 1)), l) =                                         &
                        sum(this%rhsInterior(1:n) * (xWithGhostPoints(i+1:i+n,j,k,l) -       &
                        xWithGhostPoints(i-1:i-n:-1,j,k,l)))
                end do
             end do
          end do
       end do

    case (SYMMETRIC)
       do l = 1, size(x, 2)
          do k = 1, gridSize(3)
             do j = 1, gridSize(2)
                do i = is, ie
                   x(i - this%nGhost(1) + gridSize(1) * (j - 1 +                             &
                        gridSize(2) * (k - 1)), l) =                                         &
                        sum(this%rhsInterior(1:n) * (xWithGhostPoints(i+1:i+n,j,k,l) +       &
                        xWithGhostPoints(i-1:i-n:-1,j,k,l))) +                               &
                        this%rhsInterior(0) * xWithGhostPoints(i,j,k,l)
                end do
             end do
          end do
       end do

    end select

  end subroutine applyOperatorAtInteriorPoints_1

  subroutine applyOperatorAtInteriorPoints_2(this, xWithGhostPoints, x, gridSize)

    ! <<< Derived types >>>
    use StencilOperator_type

    ! <<< Arguments >>>
    type(t_stencilOperator) :: this
    SCALAR_TYPE, intent(in) :: xWithGhostPoints(:,:,:,:)
    SCALAR_TYPE, intent(out) :: x(:,:)
    integer, intent(in) :: gridSize(3)

    ! <<< Local variables >>>
    integer :: i, j, k, l, n, js, je

    n = this%interiorWidth / 2
    js = 1 + this%nGhost(1)
    je = gridSize(2) + this%nGhost(1)
    if (this%nGhost(1) == 0) js = js + this%boundaryDepth
    if (this%nGhost(2) == 0) je = je - this%boundaryDepth

    select case (this%symmetryType) !... save FLOPS based on symmetry of interior stencil.

    case (SKEW_SYMMETRIC)
       do l = 1, size(x, 2)
          do k = 1, gridSize(3)
             do j = js, je
                do i = 1, gridSize(1)
                   x(i + gridSize(1) * (j - 1 - this%nGhost(1) +                             &
                        gridSize(2) * (k - 1)), l) =                                         &
                        sum(this%rhsInterior(1:n) * (xWithGhostPoints(i,j+1:j+n,k,l) -       &
                        xWithGhostPoints(i,j-1:j-n:-1,k,l)))
                end do
             end do
          end do
       end do

    case (SYMMETRIC)
       do l = 1, size(x, 2)
          do k = 1, gridSize(3)
             do j = js, je
                do i = 1, gridSize(1)
                   x(i + gridSize(1) * (j - 1 - this%nGhost(1) +                             &
                        gridSize(2) * (k - 1)), l) =                                         &
                        sum(this%rhsInterior(1:n) * (xWithGhostPoints(i,j+1:j+n,k,l) +       &
                        xWithGhostPoints(i,j-1:j-n:-1,k,l))) +                               &
                        this%rhsInterior(0) * xWithGhostPoints(i,j,k,l)
                end do
             end do
          end do
       end do

    end select

  end subroutine applyOperatorAtInteriorPoints_2

  subroutine applyOperatorAtInteriorPoints_3(this, xWithGhostPoints, x, gridSize)

    ! <<< Derived types >>>
    use StencilOperator_type

    ! <<< Arguments >>>
    type(t_stencilOperator) :: this
    SCALAR_TYPE, intent(in) :: xWithGhostPoints(:,:,:,:)
    SCALAR_TYPE, intent(out) :: x(:,:)
    integer, intent(in) :: gridSize(3)

    ! <<< Local variables >>>
    integer :: i, j, k, l, n, ks, ke

    n = this%interiorWidth / 2
    ks = 1 + this%nGhost(1)
    ke = gridSize(3) + this%nGhost(1)
    if (this%nGhost(1) == 0) ks = ks + this%boundaryDepth
    if (this%nGhost(2) == 0) ke = ke - this%boundaryDepth

    select case (this%symmetryType) !... save FLOPS based on symmetry of interior stencil.

    case (SKEW_SYMMETRIC)
       do l = 1, size(x, 2)
          do k = ks, ke
             do j = 1, gridSize(2)
                do i = 1, gridSize(1)
                   x(i + gridSize(1) * (j - 1 +                                              &
                        gridSize(2) * (k - 1 - this%nGhost(1))), l) =                        &
                        sum(this%rhsInterior(1:n) * (xWithGhostPoints(i,j,k+1:k+n,l) -       &
                        xWithGhostPoints(i,j,k-1:k-n:-1,l)))
                end do
             end do
          end do
       end do

    case (SYMMETRIC)
       do l = 1, size(x, 2)
          do k = ks, ke
             do j = 1, gridSize(2)
                do i = 1, gridSize(1)
                   x(i + gridSize(1) * (j - 1 +                                              &
                        gridSize(2) * (k - 1 - this%nGhost(1))), l) =                        &
                        sum(this%rhsInterior(1:n) * (xWithGhostPoints(i,j,k+1:k+n,l) +       &
                        xWithGhostPoints(i,j,k-1:k-n:-1,l))) +                               &
                        this%rhsInterior(0) * xWithGhostPoints(i,j,k,l)
                end do
             end do
          end do
       end do

    end select

  end subroutine applyOperatorAtInteriorPoints_3

  subroutine applyOperatorNorm_1(this, x, gridSize)

    ! <<< Derived types >>>
    use StencilOperator_type

    ! <<< Arguments >>>
    type(t_StencilOperator) :: this
    SCALAR_TYPE, intent(inout) :: x(:,:)
    integer, intent(in) :: gridSize(3)

    ! <<< Local variables >>>
    integer :: i, j, k, l

    ! Left boundary.
    if (this%hasDomainBoundary(1)) then
       do l = 1, size(x, 2)
          do k = 1, gridSize(3)
             do j = 1, gridSize(2)
                do i = 1, this%boundaryDepth
                   x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) =                 &
                        this%normBoundary(i) *                                               &
                        x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l)
                end do
             end do
          end do
       end do
    end if

    ! Right boundary.
    if (this%hasDomainBoundary(2)) then
       do l = 1, size(x, 2)
          do k = 1, gridSize(3)
             do j = 1, gridSize(2)
                do i = gridSize(1) + 1 - this%boundaryDepth, gridSize(1)
                   x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) =                 &
                        this%normBoundary(gridSize(1) + 1 - i) *                             &
                        x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l)
                end do
             end do
          end do
       end do
    end if

  end subroutine applyOperatorNorm_1

  subroutine applyOperatorNorm_2(this, x, gridSize)

    ! <<< Derived types >>>
    use StencilOperator_type

    ! <<< Arguments >>>
    type(t_StencilOperator) :: this
    SCALAR_TYPE, intent(inout) :: x(:,:)
    integer, intent(in) :: gridSize(3)

    ! <<< Local variables >>>
    integer :: i, j, k, l

    ! Left boundary.
    if (this%hasDomainBoundary(1)) then
       do l = 1, size(x, 2)
          do k = 1, gridSize(3)
             do j = 1, this%boundaryDepth
                do i = 1, gridSize(1)
                   x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) =                 &
                        this%normBoundary(j) *                                               &
                        x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l)
                end do
             end do
          end do
       end do
    end if

    ! Right boundary.
    if (this%hasDomainBoundary(2)) then
       do l = 1, size(x, 2)
          do k = 1, gridSize(3)
             do j = gridSize(2) + 1 - this%boundaryDepth, gridSize(2)
                do i = 1, gridSize(1)
                   x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) =                 &
                        this%normBoundary(gridSize(2) + 1 - j) *                             &
                        x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l)
                end do
             end do
          end do
       end do
    end if

  end subroutine applyOperatorNorm_2

  subroutine applyOperatorNorm_3(this, x, gridSize)

    ! <<< Derived types >>>
    use StencilOperator_type

    ! <<< Arguments >>>
    type(t_StencilOperator) :: this
    SCALAR_TYPE, intent(inout) :: x(:,:)
    integer, intent(in) :: gridSize(3)

    ! <<< Local variables >>>
    integer :: i, j, k, l

    ! Left boundary.
    if (this%hasDomainBoundary(1)) then
       do l = 1, size(x, 2)
          do k = 1, this%boundaryDepth
             do j = 1, gridSize(2)
                do i = 1, gridSize(1)
                   x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) =                 &
                        this%normBoundary(k) *                                               &
                        x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l)
                end do
             end do
          end do
       end do
    end if

    ! Right boundary.
    if (this%hasDomainBoundary(2)) then
       do l = 1, size(x, 2)
          do k = gridSize(3) + 1 - this%boundaryDepth, gridSize(3)
             do j = 1, gridSize(2)
                do i = 1, gridSize(1)
                   x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) =                 &
                        this%normBoundary(gridSize(3) + 1 - k) *                             &
                        x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l)
                end do
             end do
          end do
       end do
    end if

  end subroutine applyOperatorNorm_3

  subroutine applyOperatorNormInverse_1(this, x, gridSize)

    ! <<< Derived types >>>
    use StencilOperator_type

    ! <<< Arguments >>>
    type(t_StencilOperator) :: this
    SCALAR_TYPE, intent(inout) :: x(:,:)
    integer, intent(in) :: gridSize(3)

    ! <<< Local variables >>>
    integer :: i, j, k, l

    ! Left boundary.
    if (this%hasDomainBoundary(1)) then
       do l = 1, size(x, 2)
          do k = 1, gridSize(3)
             do j = 1, gridSize(2)
                do i = 1, this%boundaryDepth
                   x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) =                 &
                        x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) /            &
                        this%normBoundary(i)
                end do
             end do
          end do
       end do
    end if

    ! Right boundary.
    if (this%hasDomainBoundary(2)) then
       do l = 1, size(x, 2)
          do k = 1, gridSize(3)
             do j = 1, gridSize(2)
                do i = gridSize(1) + 1 - this%boundaryDepth, gridSize(1)
                   x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) =                 &
                        x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) /            &
                        this%normBoundary(gridSize(1) + 1 - i)

                end do
             end do
          end do
       end do
    end if

  end subroutine applyOperatorNormInverse_1

  subroutine applyOperatorNormInverse_2(this, x, gridSize)

    ! <<< Derived types >>>
    use StencilOperator_type

    ! <<< Arguments >>>
    type(t_StencilOperator) :: this
    SCALAR_TYPE, intent(inout) :: x(:,:)
    integer, intent(in) :: gridSize(3)

    ! <<< Local variables >>>
    integer :: i, j, k, l

    ! Left boundary.
    if (this%hasDomainBoundary(1)) then
       do l = 1, size(x, 2)
          do k = 1, gridSize(3)
             do j = 1, this%boundaryDepth
                do i = 1, gridSize(1)
                   x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) =                 &
                        x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) /            &
                        this%normBoundary(j)
                end do
             end do
          end do
       end do
    end if

    ! Right boundary.
    if (this%hasDomainBoundary(2)) then
       do l = 1, size(x, 2)
          do k = 1, gridSize(3)
             do j = gridSize(2) + 1 - this%boundaryDepth, gridSize(2)
                do i = 1, gridSize(1)
                   x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) =                 &
                        x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) /            &
                        this%normBoundary(gridSize(2) + 1 - j)
                end do
             end do
          end do
       end do
    end if

  end subroutine applyOperatorNormInverse_2

  subroutine applyOperatorNormInverse_3(this, x, gridSize)

    ! <<< Derived types >>>
    use StencilOperator_type

    ! <<< Arguments >>>
    type(t_StencilOperator) :: this
    SCALAR_TYPE, intent(inout) :: x(:,:)
    integer, intent(in) :: gridSize(3)

    ! <<< Local variables >>>
    integer :: i, j, k, l

    ! Left boundary.
    if (this%hasDomainBoundary(1)) then
       do l = 1, size(x, 2)
          do k = 1, this%boundaryDepth
             do j = 1, gridSize(2)
                do i = 1, gridSize(1)
                   x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) =                 &
                        x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) /            &
                        this%normBoundary(k)
                end do
             end do
          end do
       end do
    end if

    ! Right boundary.
    if (this%hasDomainBoundary(2)) then
       do l = 1, size(x, 2)
          do k = gridSize(3) + 1 - this%boundaryDepth, gridSize(3)
             do j = 1, gridSize(2)
                do i = 1, gridSize(1)
                   x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) =                 &
                        x(i + gridSize(1) * (j - 1 + gridSize(2) * (k - 1)), l) /            &
                        this%normBoundary(gridSize(3) + 1 - k)
                end do
             end do
          end do
       end do
    end if

  end subroutine applyOperatorNormInverse_3

end module StencilOperatorImpl

subroutine setupOperator(this, identifier, success)

  ! <<< Derived types >>>
  use StencilOperator_type

  ! <<< Private members >>>
  use StencilOperatorImpl, only : allocateData

  ! <<< Public members >>>
  use StencilOperator_mod, only : cleanupOperator

  implicit none

  ! <<< Arguments >>>
  type(t_StencilOperator) :: this
  character(len = *), intent(in) :: identifier
  logical, intent(out), optional :: success

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp) :: x1, x2, x3

  assert_key(identifier, (     \
  'SBP 1-2 first derivative',  \
  'SBP 1-2 second derivative', \
  'SBP 1-2 dissipation',       \
  'SBP 2-4 first derivative',  \
  'SBP 2-4 second derivative', \
  'SBP 2-4 dissipation',       \
  'SBP 3-6 first derivative',  \
  'SBP 3-6 second derivative', \
  'SBP 3-6 dissipation',       \
  'SBP 4-8 first derivative',  \
  'SBP 4-8 second derivative', \
  'SBP 4-8 dissipation',       \
  'null matrix'                \
  ))

  call cleanupOperator(this)
  if (present(success)) success = .true.

  if (trim(identifier) == "SBP 1-2 first derivative") then

     this%symmetryType = SKEW_SYMMETRIC
     this%interiorWidth = 3
     this%boundaryWidth = 2
     this%boundaryDepth = 1
     call allocateData(this)

     this%rhsInterior(1:1) = (/ 1.0_wp / 2.0_wp /)
     this%rhsInterior(-1:-1:-1) = - this%rhsInterior(1:1)

     this%normBoundary = (/ 1.0_wp / 2.0_wp /)

     this%rhsBoundary1(1:2,1) = (/ -1.0_wp, 1.0_wp /)

  else if (trim(identifier) == "SBP 1-2 second derivative") then

     this%symmetryType = SYMMETRIC
     this%interiorWidth = 3
     this%boundaryWidth = 3
     this%boundaryDepth = 1
     call allocateData(this)

     this%rhsInterior(0:1) = (/ -2.0_wp, 1.0_wp /)
     this%rhsInterior(-1:-1:-1) = this%rhsInterior(1:1)

     this%normBoundary = (/ 1.0_wp / 2.0_wp /)

     this%rhsBoundary1(1:3,1) = (/ 1.0_wp, -2.0_wp, 1.0_wp /)

  else if (trim(identifier) == "SBP 1-2 dissipation") then

     this%symmetryType = SYMMETRIC
     this%interiorWidth = 3
     this%boundaryWidth = 2
     this%boundaryDepth = 1
     call allocateData(this)

     this%rhsInterior(0:1) = (/ -2.0_wp, 1.0_wp /)
     this%rhsInterior(-1:-1:-1) = this%rhsInterior(1:1)
     this%rhsInterior = this%rhsInterior / 2.0_wp

     this%normBoundary = (/ 1.0_wp / 2.0_wp /)

     this%rhsBoundary1(1:2,1) = (/ -2.0_wp, 2.0_wp /)

     this%rhsBoundary1 = this%rhsBoundary1 / 2.0_wp

  else if (trim(identifier) == "SBP 2-4 first derivative") then

     this%symmetryType = SKEW_SYMMETRIC
     this%interiorWidth = 5
     this%boundaryWidth = 6
     this%boundaryDepth = 4
     call allocateData(this)

     this%rhsInterior(1:2) = (/ 2.0_wp / 3.0_wp, -1.0_wp / 12.0_wp /)
     this%rhsInterior(-1:-2:-1) = - this%rhsInterior(1:2)

     this%normBoundary = (/ 17.0_wp / 48.0_wp, &
                            59.0_wp / 48.0_wp, &
                            43.0_wp / 48.0_wp, &
                            49.0_wp / 48.0_wp /)

     this%rhsBoundary1(1:4,1) = (/ -24.0_wp / 17.0_wp, &
                                    59.0_wp / 34.0_wp, &
                                    -4.0_wp / 17.0_wp, &
                                    -3.0_wp / 34.0_wp /)
     this%rhsBoundary1(1:3,2) = (/  -1.0_wp / 2.0_wp,  &
                                              0.0_wp,  &
                                     1.0_wp / 2.0_wp /)
     this%rhsBoundary1(1:5,3) = (/   4.0_wp / 43.0_wp, &
                                   -59.0_wp / 86.0_wp, &
                                               0.0_wp, &
                                    59.0_wp / 86.0_wp, &
                                    -4.0_wp / 43.0_wp /)
     this%rhsBoundary1(1:6,4) = (/   3.0_wp / 98.0_wp, &
                                               0.0_wp, &
                                   -59.0_wp / 98.0_wp, &
                                               0.0_wp, &
                                    32.0_wp / 49.0_wp, &
                                    -4.0_wp / 49.0_wp /)

  else if (trim(identifier) == "SBP 2-4 second derivative") then

     this%symmetryType = SYMMETRIC
     this%interiorWidth = 5
     this%boundaryWidth = 6
     this%boundaryDepth = 4
     call allocateData(this)

     this%rhsInterior(0:2) = (/ -5.0_wp / 2.0_wp, 4.0_wp / 3.0_wp, -1.0_wp / 12.0_wp /)
     this%rhsInterior(-1:-2:-1) = this%rhsInterior(1:2)

     this%normBoundary = (/ 17.0_wp / 48.0_wp, &
                            59.0_wp / 48.0_wp, &
                            43.0_wp / 48.0_wp, &
                            49.0_wp / 48.0_wp /)

     this%rhsBoundary1(1:4,1) = (/              2.0_wp, &
                                               -5.0_wp, &
                                                4.0_wp, &
                                               -1.0_wp /)
     this%rhsBoundary1(1:3,2) = (/              1.0_wp, &
                                               -2.0_wp, &
                                                1.0_wp /)
     this%rhsBoundary1(1:5,3) = (/   -4.0_wp / 43.0_wp, &
                                     59.0_wp / 43.0_wp, &
                                   -110.0_wp / 43.0_wp, &
                                     59.0_wp / 43.0_wp, &
                                     -4.0_wp / 43.0_wp /)
     this%rhsBoundary1(1:6,4) = (/   -1.0_wp / 49.0_wp, &
                                                0.0_wp, &
                                     59.0_wp / 49.0_wp, &
                                   -118.0_wp / 49.0_wp, &
                                     64.0_wp / 49.0_wp, &
                                     -4.0_wp / 49.0_wp /)

  else if (trim(identifier) == "SBP 2-4 dissipation") then

     this%symmetryType = SYMMETRIC
     this%interiorWidth = 5
     this%boundaryWidth = 6
     this%boundaryDepth = 4
     call allocateData(this)

     this%rhsInterior(0:2) = (/ -6.0_wp, 4.0_wp, -1.0_wp /)
     this%rhsInterior(-1:-2:-1) = this%rhsInterior(1:2)
     this%rhsInterior = this%rhsInterior / 16.0_wp

     this%normBoundary = (/ 17.0_wp / 48.0_wp, &
                            59.0_wp / 48.0_wp, &
                            43.0_wp / 48.0_wp, &
                            49.0_wp / 48.0_wp /)

     this%rhsBoundary1(1:3,1) = (/  -96.0_wp / 17.0_wp, &
                                    192.0_wp / 17.0_wp, &
                                    -96.0_wp / 17.0_wp /)
     this%rhsBoundary1(1:4,2) = (/  192.0_wp / 59.0_wp, &
                                   -432.0_wp / 59.0_wp, &
                                    288.0_wp / 59.0_wp, &
                                    -48.0_wp / 59.0_wp /)
     this%rhsBoundary1(1:5,3) = (/  -96.0_wp / 43.0_wp, &
                                    288.0_wp / 43.0_wp, &
                                   -336.0_wp / 43.0_wp, &
                                    192.0_wp / 43.0_wp, &
                                    -48.0_wp / 43.0_wp /)
     this%rhsBoundary1(2:6,4) = (/  -48.0_wp / 49.0_wp, &
                                    192.0_wp / 49.0_wp, &
                                   -288.0_wp / 49.0_wp, &
                                    192.0_wp / 49.0_wp, &
                                    -48.0_wp / 49.0_wp /)

     this%rhsBoundary1 = this%rhsBoundary1 / 16.0_wp

  else if (trim(identifier) == "SBP 3-6 first derivative") then

     this%symmetryType = SKEW_SYMMETRIC
     this%interiorWidth = 7
     this%boundaryWidth = 9
     this%boundaryDepth = 6
     call allocateData(this)

     this%rhsInterior(1:3) = (/ 3.0_wp / 4.0_wp, -3.0_wp / 20.0_wp, 1.0_wp / 60.0_wp /)
     this%rhsInterior(-1:-3:-1) = - this%rhsInterior(1:3)

     this%normBoundary = (/ 13649.0_wp / 43200.0_wp, &
                            12013.0_wp /  8640.0_wp, &
                             2711.0_wp /  4320.0_wp, &
                             5359.0_wp /  4320.0_wp, &
                             7877.0_wp /  8640.0_wp, &
                            43801.0_wp / 43200.0_wp /)

     this%rhsBoundary1(1:6,1) = (/  -21600.0_wp /  13649.0_wp, &
                                    104009.0_wp /  54596.0_wp, &
                                     30443.0_wp /  81894.0_wp, &
                                    -33311.0_wp /  27298.0_wp, &
                                     16863.0_wp /  27298.0_wp, &
                                    -15025.0_wp / 163788.0_wp /)
     this%rhsBoundary1(1:6,2) = (/ -104009.0_wp / 240260.0_wp, &
                                                       0.0_wp, &
                                      -311.0_wp /  72078.0_wp, &
                                     20229.0_wp /  24026.0_wp, &
                                    -24337.0_wp /  48052.0_wp, &
                                     36661.0_wp / 360390.0_wp /)
     this%rhsBoundary1(1:6,3) = (/  -30443.0_wp / 162660.0_wp, &
                                       311.0_wp /  32532.0_wp, &
                                                       0.0_wp, &
                                    -11155.0_wp /  16266.0_wp, &
                                     41287.0_wp /  32532.0_wp, &
                                    -21999.0_wp /  54220.0_wp /)
     this%rhsBoundary1(1:7,4) = (/   33311.0_wp / 107180.0_wp, &
                                    -20229.0_wp /  21436.0_wp, &
                                       485.0_wp /   1398.0_wp, &
                                                       0.0_wp, &
                                      4147.0_wp /  21436.0_wp, &
                                     25427.0_wp / 321540.0_wp, &
                                        72.0_wp /   5359.0_wp /)
     this%rhsBoundary1(1:8,5) = (/  -16863.0_wp /  78770.0_wp, &
                                     24337.0_wp /  31508.0_wp, &
                                    -41287.0_wp /  47262.0_wp, &
                                     -4147.0_wp /  15754.0_wp, &
                                                       0.0_wp, &
                                    342523.0_wp / 472620.0_wp, &
                                     -1296.0_wp /   7877.0_wp, &
                                       144.0_wp /   7877.0_wp /)
     this%rhsBoundary1(1:9,6) = (/   15025.0_wp / 525612.0_wp, &
                                    -36661.0_wp / 262806.0_wp, &
                                     21999.0_wp /  87602.0_wp, &
                                    -25427.0_wp / 262806.0_wp, &
                                   -342523.0_wp / 525612.0_wp, &
                                                       0.0_wp, &
                                     32400.0_wp /  43801.0_wp, &
                                     -6480.0_wp /  43801.0_wp, &
                                       720.0_wp /  43801.0_wp /)

  else if (trim(identifier) == "SBP 3-6 second derivative") then

     this%symmetryType = SYMMETRIC
     this%interiorWidth = 7
     this%boundaryWidth = 9
     this%boundaryDepth = 6
     call allocateData(this)

     this%rhsInterior(0:3) = (/ -49.0_wp / 18.0_wp, &
                                  3.0_wp /  2.0_wp, &
                                 -3.0_wp / 20.0_wp, &
                                  1.0_wp / 90.0_wp /)
     this%rhsInterior(-1:-3:-1) = this%rhsInterior(1:3)

     this%normBoundary = (/ 13649.0_wp / 43200.0_wp, &
                            12013.0_wp /  8640.0_wp, &
                             2711.0_wp /  4320.0_wp, &
                             5359.0_wp /  4320.0_wp, &
                             7877.0_wp /  8640.0_wp, &
                            43801.0_wp / 43200.0_wp /)

     this%rhsBoundary1(1:6,1) = (/  114170.0_wp /  40947.0_wp, &
                                   -438107.0_wp /  54596.0_wp, &
                                    336409.0_wp /  40947.0_wp, &
                                   -276997.0_wp /  81894.0_wp, &
                                      3747.0_wp /  13649.0_wp, &
                                     21035.0_wp / 163788.0_wp /)
     this%rhsBoundary1(1:6,2) = (/    6173.0_wp /   5860.0_wp, &
                                     -2066.0_wp /    879.0_wp, &
                                      3283.0_wp /   1758.0_wp, &
                                      -303.0_wp /    293.0_wp, &
                                      2111.0_wp /   3516.0_wp, &
                                      -601.0_wp /   4395.0_wp /)
     this%rhsBoundary1(1:6,3) = (/  -52391.0_wp /  81330.0_wp, &
                                    134603.0_wp /  32532.0_wp, &
                                    -21982.0_wp /   2711.0_wp, &
                                    112915.0_wp /  16266.0_wp, &
                                    -46969.0_wp /  16266.0_wp, &
                                     30409.0_wp /  54220.0_wp /)
     this%rhsBoundary1(1:7,4) = (/   68603.0_wp / 321540.0_wp, &
                                    -12423.0_wp /  10718.0_wp, &
                                    112915.0_wp /  32154.0_wp, &
                                    -75934.0_wp /  16077.0_wp, &
                                     53369.0_wp /  21436.0_wp, &
                                    -54899.0_wp / 160770.0_wp, &
                                        48.0_wp /   5359.0_wp /)
     this%rhsBoundary1(1:8,5) = (/   -7053.0_wp /  39385.0_wp, &
                                     86551.0_wp /  94524.0_wp, &
                                    -46969.0_wp /  23631.0_wp, &
                                     53369.0_wp /  15754.0_wp, &
                                    -87904.0_wp /  23631.0_wp, &
                                    820271.0_wp / 472620.0_wp, &
                                     -1296.0_wp /   7877.0_wp, &
                                        96.0_wp /   7877.0_wp /)
     this%rhsBoundary1(1:9,6) = (/   21035.0_wp / 525612.0_wp, &
                                    -24641.0_wp / 131403.0_wp, &
                                     30409.0_wp /  87602.0_wp, &
                                    -54899.0_wp / 131403.0_wp, &
                                    820271.0_wp / 525612.0_wp, &
                                   -117600.0_wp /  43801.0_wp, &
                                     64800.0_wp /  43801.0_wp, &
                                     -6480.0_wp /  43801.0_wp, &
                                       480.0_wp /  43801.0_wp /)

  else if (trim(identifier) == "SBP 3-6 dissipation") then

     this%symmetryType = SYMMETRIC
     this%interiorWidth = 7
     this%boundaryWidth = 9
     this%boundaryDepth = 6
     call allocateData(this)

     this%rhsInterior(0:3) = (/ -20.0_wp, 15.0_wp, -6.0_wp, 1.0_wp /)
     this%rhsInterior(-1:-3:-1) = this%rhsInterior(1:3)
     this%rhsInterior = this%rhsInterior / 64.0_wp

     this%normBoundary = (/ 13649.0_wp / 43200.0_wp, &
                            12013.0_wp /  8640.0_wp, &
                             2711.0_wp /  4320.0_wp, &
                             5359.0_wp /  4320.0_wp, &
                             7877.0_wp /  8640.0_wp, &
                            43801.0_wp / 43200.0_wp /)

     this%rhsBoundary1(1:4,1) = (/ -129600.0_wp / 13649.0_wp, &
                                    388800.0_wp / 13649.0_wp, &
                                   -388800.0_wp / 13649.0_wp, &
                                    129600.0_wp / 13649.0_wp /)
     this%rhsBoundary1(1:5,2) = (/   77760.0_wp / 12013.0_wp, &
                                   -241920.0_wp / 12013.0_wp, &
                                    259200.0_wp / 12013.0_wp, &
                                   -103680.0_wp / 12013.0_wp, &
                                      8640.0_wp / 12013.0_wp /)
     this%rhsBoundary1(1:6,3) = (/  -38880.0_wp /  2711.0_wp, &
                                    129600.0_wp /  2711.0_wp, &
                                   -159840.0_wp /  2711.0_wp, &
                                     90720.0_wp /  2711.0_wp, &
                                    -25920.0_wp /  2711.0_wp, &
                                      4320.0_wp /  2711.0_wp /)
     this%rhsBoundary1(1:7,4) = (/   12960.0_wp /  5359.0_wp, &
                                    -51840.0_wp /  5359.0_wp, &
                                     90720.0_wp /  5359.0_wp, &
                                    -95040.0_wp /  5359.0_wp, &
                                     64800.0_wp /  5359.0_wp, &
                                    -25920.0_wp /  5359.0_wp, &
                                      4320.0_wp /  5359.0_wp /)
     this%rhsBoundary1(2:8,5) = (/    8640.0_wp /  7877.0_wp, &
                                    -51840.0_wp /  7877.0_wp, &
                                    129600.0_wp /  7877.0_wp, &
                                   -172800.0_wp /  7877.0_wp, &
                                    129600.0_wp /  7877.0_wp, &
                                    -51840.0_wp /  7877.0_wp, &
                                      8640.0_wp /  7877.0_wp /)
     this%rhsBoundary1(3:9,6) = (/   43200.0_wp / 43801.0_wp, &
                                   -259200.0_wp / 43801.0_wp, &
                                    648000.0_wp / 43801.0_wp, &
                                   -864000.0_wp / 43801.0_wp, &
                                    648000.0_wp / 43801.0_wp, &
                                   -259200.0_wp / 43801.0_wp, &
                                     43200.0_wp / 43801.0_wp /)

     this%rhsBoundary1 = this%rhsBoundary1 / 64.0_wp

  else if (trim(identifier) == "SBP 4-8 first derivative") then

     this%symmetryType = SKEW_SYMMETRIC
     this%interiorWidth = 9
     this%boundaryWidth = 12
     this%boundaryDepth = 8
     call allocateData(this)

     this%rhsInterior(1:4) = (/  4.0_wp /   5.0_wp, &
                                -1.0_wp /   5.0_wp, &
                                 4.0_wp / 105.0_wp, &
                                -1.0_wp / 280.0_wp /)
     this%rhsInterior(-1:-4:-1) = - this%rhsInterior(1:4)

     this%normBoundary = (/ 1498139.0_wp / 5080320.0_wp, &
                            1107307.0_wp /  725760.0_wp, &
                              20761.0_wp /   80640.0_wp, &
                            1304999.0_wp /  725760.0_wp, &
                             299527.0_wp /  725760.0_wp, &
                             103097.0_wp /   80640.0_wp, &
                             670091.0_wp /  725760.0_wp, &
                            5127739.0_wp / 5080320.0_wp /)

     x1 =  541.0_wp / 1000.0_wp
     x2 = - 27.0_wp /  400.0_wp
     x3 =  187.0_wp /  250.0_wp

     this%rhsBoundary1(1,1) = -2540160.0_wp / 1498139.0_wp
     this%rhsBoundary1(2,1) = 9.0_wp * (2257920.0_wp * x1 + 11289600.0_wp * x2     &
          + 22579200.0_wp * x3 - 15849163.0_wp) / 5992556.0_wp
     this%rhsBoundary1(3,1) = 3.0_wp * (-33868800.0_wp * x1 - 162570240.0_wp * x2  &
          - 304819200.0_wp * x3 + 235236677.0_wp) / 5992556.0_wp
     this%rhsBoundary1(4,1) = (609638400.0_wp * x1 + 2743372800.0_wp * x2          &
          + 4572288000.0_wp * x3 - 3577778591.0_wp) / 17977668.0_wp
     this%rhsBoundary1(5,1) = 3.0_wp * (-16934400 * x1 - 67737600.0_wp * x2        &
          - 84672000.0_wp * x3 + 67906303.0_wp) / 1498139.0_wp
     this%rhsBoundary1(6,1) = 105.0_wp * (967680.0_wp * x1 + 2903040.0_wp * x2     &
          - 305821.0_wp) / 5992556.0_wp
     this%rhsBoundary1(7,1) = 49.0_wp * (-1244160.0_wp * x1 + 18662400.0_wp * x3   &
          - 13322233.0_wp) / 17977668.0_wp
     this%rhsBoundary1(8,1) = 3.0_wp * (-6773760.0_wp * x2 - 33868800.0_wp * x3    &
          + 24839327.0_wp) / 5992556.0_wp

     this%rhsBoundary1(1,2) = 9.0_wp * (-2257920.0_wp * x1 - 11289600.0_wp * x2    &
          - 22579200.0_wp * x3 + 15849163.0_wp) / 31004596.0_wp
     this%rhsBoundary1(2,2) = 0.0_wp
     this%rhsBoundary1(3,2) = 3.0_wp * (7257600.0_wp * x1 + 33868800.0_wp * x2     &
          + 60963840.0_wp * x3 - 47167457.0_wp) / 2214614.0_wp
     this%rhsBoundary1(4,2) = 3.0_wp * (-9676800.0_wp * x1 - 42336000.0_wp * x2    &
          - 67737600.0_wp * x3 + 53224573.0_wp) / 1107307.0_wp
     this%rhsBoundary1(5,2) = 7.0_wp * (55987200.0_wp * x1 + 217728000.0_wp * x2   &
          + 261273600.0_wp * x3 - 211102099.0_wp) / 13287684.0_wp
     this%rhsBoundary1(6,2) = 3.0_wp * (-11612160.0_wp * x1 - 33868800.0_wp * x2   &
          + 3884117.0_wp) / 2214614.0_wp
     this%rhsBoundary1(7,2) = 150.0_wp * (24192.0_wp * x1 - 338688.0_wp * x3       &
          + 240463.0_wp) / 1107307.0_wp
     this%rhsBoundary1(8,2) = (152409600.0_wp * x2 + 731566080.0_wp * x3           &
          - 536324953.0_wp) / 46506894.0_wp

     this%rhsBoundary1(1,3) = (33868800.0_wp * x1 + 162570240.0_wp * x2            &
          + 304819200.0_wp * x3 - 235236677.0_wp) / 1743924.0_wp
     this%rhsBoundary1(2,3) = (-7257600.0_wp * x1 - 33868800.0_wp * x2             &
          - 60963840.0_wp * x3 + 47167457.0_wp) / 124566.0_wp
     this%rhsBoundary1(3,3) = 0.0_wp
     this%rhsBoundary1(4,3) = (24192000.0_wp * x1 + 101606400.0_wp * x2            &
          + 152409600.0_wp * x3 - 120219461.0_wp) / 124566.0_wp
     this%rhsBoundary1(5,3) = (-72576000.0_wp * x1 - 270950400.0_wp * x2           &
          - 304819200.0_wp * x3 + 249289259.0_wp) / 249132.0_wp
     this%rhsBoundary1(6,3) = 9.0_wp * (806400.0_wp * x1 + 2257920.0_wp * x2       &
          - 290167.0_wp) / 41522.0_wp
     this%rhsBoundary1(7,3) = 6.0_wp * (-134400.0_wp * x1 + 1693440.0_wp * x3      &
          - 1191611.0_wp) / 20761.0_wp
     this%rhsBoundary1(8,3) = 5.0_wp * (-2257920.0_wp * x2 - 10160640.0_wp * x3    &
          + 7439833.0_wp) / 290654.0_wp

     this%rhsBoundary1(1,4) = (-609638400.0_wp * x1 - 2743372800.0_wp * x2         &
          - 4572288000.0_wp * x3 + 3577778591.0_wp) / 109619916.0_wp
     this%rhsBoundary1(2,4) = 3.0_wp * (9676800.0_wp * x1 + 42336000.0_wp * x2     &
          + 67737600.0_wp * x3 - 53224573.0_wp) / 1304999.0_wp
     this%rhsBoundary1(3,4) = 3.0_wp * (-24192000.0_wp * x1 - 101606400.0_wp * x2  &
          - 152409600.0_wp * x3 + 120219461.0_wp) / 2609998.0_wp
     this%rhsBoundary1(4,4) = 0.0_wp
     this%rhsBoundary1(5,4) = 9.0_wp * (16128000.0_wp * x1 + 56448000.0_wp * x2    &
          + 56448000.0_wp * x3 - 47206049.0_wp) / 5219996.0_wp
     this%rhsBoundary1(6,4) = 3.0_wp * (-19353600.0_wp * x1 - 50803200.0_wp * x2   &
          + 7628371.0_wp) / 2609998.0_wp
     this%rhsBoundary1(7,4) = 2.0_wp * (10886400.0_wp * x1 - 114307200.0_wp * x3   &
          + 79048289.0_wp) / 3914997.0_wp
     this%rhsBoundary1(8,4) = 75.0_wp * (1354752.0_wp * x2 + 5419008.0_wp * x3     &
          - 3952831.0_wp) / 18269986.0_wp

     this%rhsBoundary1(1,5) = 3.0_wp * (16934400.0_wp * x1 + 67737600.0_wp * x2 +  &
          84672000.0_wp * x3 - 67906303.0_wp) / 2096689.0_wp
     this%rhsBoundary1(2,5) = 7.0_wp * (-55987200.0_wp * x1 - 217728000.0_wp * x2  &
          - 261273600.0_wp * x3 + 211102099.0_wp) / 3594324.0_wp
     this%rhsBoundary1(3,5) = 3.0_wp * (72576000.0_wp * x1 + 270950400.0_wp * x2   &
          + 304819200.0_wp * x3 - 249289259.0_wp) / 1198108.0_wp
     this%rhsBoundary1(4,5) = 9.0_wp * (-16128000.0_wp * x1 - 56448000.0_wp * x2   &
          - 56448000.0_wp * x3 + 47206049.0_wp) / 1198108.0_wp
     this%rhsBoundary1(5,5) = 0.0_wp
     this%rhsBoundary1(6,5) = 105.0_wp * (414720.0_wp * x1 + 967680.0_wp * x2      &
          - 165527.0_wp) / 1198108.0_wp
     this%rhsBoundary1(7,5) = 15.0_wp * (-967680.0_wp * x1 + 6773760.0_wp * x3     &
          - 4472029.0_wp) / 1198108.0_wp
     this%rhsBoundary1(8,5) = (-304819200.0_wp * x2 - 914457600.0_wp * x3          &
          + 657798011.0_wp) / 25160268.0_wp
     this%rhsBoundary1(9,5) = -2592.0_wp / 299527.0_wp

     this%rhsBoundary1(1,6)  = 5.0_wp * (-967680.0_wp * x1 - 2903040.0_wp * x2     &
          + 305821.0_wp) / 1237164.0_wp
     this%rhsBoundary1(2,6)  = (11612160.0_wp * x1 + 33868800.0_wp * x2            &
          - 3884117.0_wp) / 618582.0_wp
     this%rhsBoundary1(3,6)  = 9.0_wp * (-806400.0_wp * x1 - 2257920.0_wp * x2     &
          + 290167.0_wp) / 206194.0_wp
     this%rhsBoundary1(4,6)  = (19353600.0_wp * x1 + 50803200.0_wp * x2            &
          - 7628371.0_wp) / 618582.0_wp
     this%rhsBoundary1(5,6)  = 35.0_wp * (-414720.0_wp * x1 - 967680.0_wp * x2     &
          + 165527.0_wp) / 1237164.0_wp
     this%rhsBoundary1(6,6)  = 0.0_wp
     this%rhsBoundary1(7,6)  = 80640.0_wp * x1 / 103097.0_wp
     this%rhsBoundary1(8,6)  = 80640.0_wp * x2 / 103097.0_wp
     this%rhsBoundary1(9,6)  = 3072.0_wp / 103097.0_wp
     this%rhsBoundary1(10,6) = -288.0_wp / 103097.0_wp

     this%rhsBoundary1(1,7)  = 7.0_wp * (1244160.0_wp * x1 - 18662400.0_wp * x3    &
          + 13322233.0_wp) / 8041092.0_wp
     this%rhsBoundary1(2,7)  = 150.0_wp * (-24192.0_wp * x1 + 338688.0_wp * x3     &
          - 240463.0_wp) / 670091.0_wp
     this%rhsBoundary1(3,7)  = 54.0_wp * (134400.0_wp * x1 - 1693440.0_wp * x3     &
          + 1191611.0_wp) / 670091.0_wp
     this%rhsBoundary1(4,7)  = 2.0_wp * (-10886400.0_wp * x1 + 114307200.0_wp * x3 &
          - 79048289.0_wp) / 2010273.0_wp
     this%rhsBoundary1(5,7)  = 15.0_wp * (967680.0_wp * x1 - 6773760.0_wp * x3     &
          + 4472029.0_wp) / 2680364.0_wp
     this%rhsBoundary1(6,7)  = -725760.0_wp * x1 / 670091.0_wp
     this%rhsBoundary1(7,7)  = 0.0_wp
     this%rhsBoundary1(8,7)  = 725760.0_wp * x3 / 670091.0_wp
     this%rhsBoundary1(9,7)  = -145152.0_wp / 670091.0_wp
     this%rhsBoundary1(10,7) = 27648.0_wp / 670091.0_wp
     this%rhsBoundary1(11,7) = -2592.0_wp / 670091.0_wp

     this%rhsBoundary1(1,8)  = 3.0_wp * (6773760.0_wp * x2 + 33868800.0_wp * x3    &
          - 24839327.0_wp) / 20510956.0_wp
     this%rhsBoundary1(2,8)  = (-152409600.0_wp * x2 - 731566080.0_wp * x3         &
          + 536324953.0_wp) / 30766434.0_wp
     this%rhsBoundary1(3,8)  = 45.0_wp * (2257920.0_wp * x2 + 10160640.0_wp * x3   &
          - 7439833.0_wp) / 10255478.0_wp
     this%rhsBoundary1(4,8)  = 75.0_wp * (-1354752.0_wp * x2 - 5419008.0_wp * x3   &
          + 3952831.0_wp) / 10255478.0_wp
     this%rhsBoundary1(5,8)  = (304819200.0_wp * x2 + 914457600.0_wp * x3          &
          - 657798011.0_wp) / 61532868.0_wp
     this%rhsBoundary1(6,8)  = -5080320.0_wp * x2 / 5127739.0_wp
     this%rhsBoundary1(7,8)  = -5080320.0_wp * x3 / 5127739.0_wp
     this%rhsBoundary1(8,8)  = 0.0_wp
     this%rhsBoundary1(9,8)  = 4064256.0_wp / 5127739.0_wp
     this%rhsBoundary1(10,8) = -1016064.0_wp / 5127739.0_wp
     this%rhsBoundary1(11,8) = 193536.0_wp / 5127739.0_wp
     this%rhsBoundary1(12,8) = -18144.0_wp / 5127739.0_wp

  else if (trim(identifier) == "SBP 4-8 second derivative") then

     this%symmetryType = SYMMETRIC
     this%interiorWidth = 9
     this%boundaryWidth = 12
     this%boundaryDepth = 8
     call allocateData(this)

     this%rhsInterior(0:4) = (/ -205.0_wp /  72.0_wp, &
                                   8.0_wp /   5.0_wp, &
                                  -1.0_wp /   5.0_wp, &
                                   8.0_wp / 315.0_wp, &
                                  -1.0_wp / 560.0_wp /)
     this%rhsInterior(-1:-4:-1) = this%rhsInterior(1:4)

     this%normBoundary = (/ 1498139.0_wp / 5080320.0_wp, &
                            1107307.0_wp /  725760.0_wp, &
                              20761.0_wp /   80640.0_wp, &
                            1304999.0_wp /  725760.0_wp, &
                             299527.0_wp /  725760.0_wp, &
                             103097.0_wp /   80640.0_wp, &
                             670091.0_wp /  725760.0_wp, &
                            5127739.0_wp / 5080320.0_wp /)

     this%rhsBoundary1(1:7,1)  = (/  4870382994799.0_wp / 1358976868290.0_wp, &
                                     -893640087518.0_wp /   75498714905.0_wp, &
                                      926594825119.0_wp /   60398971924.0_wp, &
                                    -1315109406200.0_wp /  135897686829.0_wp, &
                                       39126983272.0_wp /   15099742981.0_wp, &
                                       12344491342.0_wp /   75498714905.0_wp, &
                                     -451560522577.0_wp / 2717953736580.0_wp /)
     this%rhsBoundary1(1:8,2)  = (/   333806012194.0_wp /  390619153855.0_wp, &
                                     -154646272029.0_wp /  111605472530.0_wp, &
                                        1168338040.0_wp /   33481641759.0_wp, &
                                       82699112501.0_wp /  133926567036.0_wp, &
                                        -171562838.0_wp /   11160547253.0_wp, &
                                      -28244698346.0_wp /  167408208795.0_wp, &
                                       11904122576.0_wp /  167408208795.0_wp, &
                                       -2598164715.0_wp /  312495323084.0_wp /)
     this%rhsBoundary1(1:8,3)  =  (/    7838984095.0_wp /   52731029988.0_wp, &
                                        1168338040.0_wp /    5649753213.0_wp, &
                                         -88747895.0_wp /     144865467.0_wp, &
                                         423587231.0_wp /     627750357.0_wp, &
                                      -43205598281.0_wp /   22599012852.0_wp, &
                                        4876378562.0_wp /    1883251071.0_wp, &
                                       -5124426509.0_wp /    3766502142.0_wp, &
                                       10496900965.0_wp /   39548272491.0_wp /)
     this%rhsBoundary1(1:8,4)  = (/   -94978241528.0_wp /  828644350023.0_wp, &
                                       82699112501.0_wp /  157837019052.0_wp, &
                                        1270761693.0_wp /   13153084921.0_wp, &
                                     -167389605005.0_wp /  118377764289.0_wp, &
                                       48242560214.0_wp /   39459254763.0_wp, &
                                      -31673996013.0_wp /   52612339684.0_wp, &
                                       43556319241.0_wp /  118377764289.0_wp, &
                                      -44430275135.0_wp /  552429566682.0_wp /)
     this%rhsBoundary1(1:9,5)  = (/     1455067816.0_wp /   21132528431.0_wp, &
                                        -171562838.0_wp /    3018932633.0_wp, &
                                      -43205598281.0_wp /   36227191596.0_wp, &
                                       48242560214.0_wp /    9056797899.0_wp, &
                                      -52276055645.0_wp /    6037865266.0_wp, &
                                       57521587238.0_wp /    9056797899.0_wp, &
                                      -80321706377.0_wp /   36227191596.0_wp, &
                                        8078087158.0_wp /   21132528431.0_wp, &
                                             -1296.0_wp /        299527.0_wp /)
     this%rhsBoundary1(1:10,6) = (/    10881504334.0_wp /  327321118845.0_wp, &
                                      -28244698346.0_wp /  140280479505.0_wp, &
                                        4876378562.0_wp /    9352031967.0_wp, &
                                      -10557998671.0_wp /   12469375956.0_wp, &
                                       57521587238.0_wp /   28056095901.0_wp, &
                                     -278531401019.0_wp /   93520319670.0_wp, &
                                       73790130002.0_wp /   46760159835.0_wp, &
                                     -137529995233.0_wp /  785570685228.0_wp, &
                                              2048.0_wp /        103097.0_wp, &
                                              -144.0_wp /        103097.0_wp /)
     this%rhsBoundary1(1:11,7) = (/  -135555328849.0_wp / 8509847458140.0_wp, &
                                       11904122576.0_wp /  101307707835.0_wp, &
                                       -5124426509.0_wp /   13507694378.0_wp, &
                                       43556319241.0_wp /   60784624701.0_wp, &
                                      -80321706377.0_wp /   81046166268.0_wp, &
                                       73790130002.0_wp /   33769235945.0_wp, &
                                     -950494905688.0_wp /  303923123505.0_wp, &
                                      239073018673.0_wp /  141830790969.0_wp, &
                                           -145152.0_wp /        670091.0_wp, &
                                             18432.0_wp /        670091.0_wp, &
                                             -1296.0_wp /        670091.0_wp /)
     this%rhsBoundary1(2:12,8) = (/    -2598164715.0_wp /  206729925524.0_wp, &
                                       10496900965.0_wp /  155047444143.0_wp, &
                                      -44430275135.0_wp /  310094888286.0_wp, &
                                         425162482.0_wp /    2720130599.0_wp, &
                                     -137529995233.0_wp /  620189776572.0_wp, &
                                      239073018673.0_wp /  155047444143.0_wp, &
                                     -144648000000.0_wp /   51682481381.0_wp, &
                                           8128512.0_wp /       5127739.0_wp, &
                                          -1016064.0_wp /       5127739.0_wp, &
                                            129024.0_wp /       5127739.0_wp, &
                                             -9072.0_wp /       5127739.0_wp /)

  else if (trim(identifier) == "SBP 4-8 dissipation") then

     this%symmetryType = SYMMETRIC
     this%interiorWidth = 9
     this%boundaryWidth = 12
     this%boundaryDepth = 8
     call allocateData(this)

     this%rhsInterior(0:4) = (/ -70.0_wp, 56.0_wp, -28.0_wp, 8.0_wp, -1.0_wp /)
     this%rhsInterior(-1:-4:-1) = this%rhsInterior(1:4)
     this%rhsInterior = this%rhsInterior / 256.0_wp

     this%normBoundary = (/ 1498139.0_wp / 5080320.0_wp, &
                            1107307.0_wp /  725760.0_wp, &
                              20761.0_wp /   80640.0_wp, &
                            1304999.0_wp /  725760.0_wp, &
                             299527.0_wp /  725760.0_wp, &
                             103097.0_wp /   80640.0_wp, &
                             670091.0_wp /  725760.0_wp, &
                            5127739.0_wp / 5080320.0_wp /)

     this%rhsBoundary1(1:5,1)  = (/  -15240960.0_wp / 1498139.0_wp, &
                                      60963840.0_wp / 1498139.0_wp, &
                                     -91445760.0_wp / 1498139.0_wp, &
                                      60963840.0_wp / 1498139.0_wp, &
                                     -15240960.0_wp / 1498139.0_wp /)
     this%rhsBoundary1(1:6,2)  = (/    8709120.0_wp / 1107307.0_wp, &
                                     -35562240.0_wp / 1107307.0_wp, &
                                      55157760.0_wp / 1107307.0_wp, &
                                     -39191040.0_wp / 1107307.0_wp, &
                                      11612160.0_wp / 1107307.0_wp, &
                                       -725760.0_wp / 1107307.0_wp /)
     this%rhsBoundary1(1:7,3)  = (/   -1451520.0_wp /   20761.0_wp, &
                                       6128640.0_wp /   20761.0_wp, &
                                     -10080000.0_wp /   20761.0_wp, &
                                       8064000.0_wp /   20761.0_wp, &
                                      -3225600.0_wp /   20761.0_wp, &
                                        645120.0_wp /   20761.0_wp, &
                                        -80640.0_wp /   20761.0_wp /)
     this%rhsBoundary1(1:8,4)  = (/    8709120.0_wp / 1304999.0_wp, &
                                     -39191040.0_wp / 1304999.0_wp, &
                                      72576000.0_wp / 1304999.0_wp, &
                                     -73301760.0_wp / 1304999.0_wp, &
                                      46448640.0_wp / 1304999.0_wp, &
                                     -20321280.0_wp / 1304999.0_wp, &
                                       5806080.0_wp / 1304999.0_wp, &
                                       -725760.0_wp / 1304999.0_wp /)
     this%rhsBoundary1(1:9,5)  = (/   -2177280.0_wp /  299527.0_wp, &
                                      11612160.0_wp /  299527.0_wp, &
                                     -29030400.0_wp /  299527.0_wp, &
                                      46448640.0_wp /  299527.0_wp, &
                                     -52254720.0_wp /  299527.0_wp, &
                                      40642560.0_wp /  299527.0_wp, &
                                     -20321280.0_wp /  299527.0_wp, &
                                       5806080.0_wp /  299527.0_wp, &
                                       -725760.0_wp /  299527.0_wp /)
     this%rhsBoundary1(2:10,6) = (/     -80640.0_wp /  103097.0_wp, &
                                        645120.0_wp /  103097.0_wp, &
                                      -2257920.0_wp /  103097.0_wp, &
                                       4515840.0_wp /  103097.0_wp, &
                                      -5644800.0_wp /  103097.0_wp, &
                                       4515840.0_wp /  103097.0_wp, &
                                      -2257920.0_wp /  103097.0_wp, &
                                        645120.0_wp /  103097.0_wp, &
                                        -80640.0_wp /  103097.0_wp /)
     this%rhsBoundary1(3:11,7) = (/    -725760.0_wp /  670091.0_wp, &
                                       5806080.0_wp /  670091.0_wp, &
                                     -20321280.0_wp /  670091.0_wp, &
                                      40642560.0_wp /  670091.0_wp, &
                                     -50803200.0_wp /  670091.0_wp, &
                                      40642560.0_wp /  670091.0_wp, &
                                     -20321280.0_wp /  670091.0_wp, &
                                       5806080.0_wp /  670091.0_wp, &
                                       -725760.0_wp /  670091.0_wp /)
     this%rhsBoundary1(4:12,8) = (/   -5080320.0_wp / 5127739.0_wp, &
                                      40642560.0_wp / 5127739.0_wp, &
                                    -142248960.0_wp / 5127739.0_wp, &
                                     284497920.0_wp / 5127739.0_wp, &
                                    -355622400.0_wp / 5127739.0_wp, &
                                     284497920.0_wp / 5127739.0_wp, &
                                    -142248960.0_wp / 5127739.0_wp, &
                                      40642560.0_wp / 5127739.0_wp, &
                                      -5080320.0_wp / 5127739.0_wp /)

     this%rhsBoundary1 = this%rhsBoundary1 / 256.0_wp

  else if (trim(identifier) == "null matrix") then

     this%symmetryType = SYMMETRIC
     this%interiorWidth = 0
     this%boundaryWidth = 1
     this%boundaryDepth = 1
     call allocateData(this)

     this%rhsInterior(0:0) = 0.0_wp
     this%rhsBoundary1 = 0.0_wp

  else

     if (present(success)) success = .false.

  end if

  ! Fill the right-boundary coefficients.
  if (allocated(this%rhsBoundary1) .and. allocated(this%rhsBoundary2)) then
     select case (this%symmetryType)
     case (SYMMETRIC)
        this%rhsBoundary2(1:this%boundaryWidth,:) =                                          &
             +this%rhsBoundary1(this%boundaryWidth:1:-1,:)
     case (SKEW_SYMMETRIC)
        this%rhsBoundary2(1:this%boundaryWidth,:) =                                          &
             -this%rhsBoundary1(this%boundaryWidth:1:-1,:)
     end select
  end if

end subroutine setupOperator

subroutine updateOperator(this, cartesianCommunicator, direction, isPeriodicityOverlapping)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use StencilOperator_type

  implicit none

  ! <<< Arguments >>>
  type(t_StencilOperator) :: this
  integer, intent(in) :: cartesianCommunicator, direction
  logical, intent(in), optional :: isPeriodicityOverlapping

  ! <<< Local variables >>>
  integer :: nDimensions, processDistribution(3), processCoordinates(3), ierror
  logical :: isPeriodic(3)

  ! Get number of dimensions from communicator.
  call MPI_Cartdim_get(cartesianCommunicator, nDimensions, ierror)
  assert(nDimensions >= 1 .and. nDimensions <= 3)
  assert(direction >= 1 .and. direction <= nDimensions)

  ! Get process distribution, coordinates and periodicity.
  processDistribution = 1
  isPeriodic = .false.
  processCoordinates = 0
  call MPI_Cart_get(cartesianCommunicator, nDimensions, processDistribution,                 &
       isPeriodic, processCoordinates, ierror)
  assert(all(processDistribution > 0))

  ! Does this process contain the boundary of the computational domain?
  this%hasDomainBoundary(1) = (processCoordinates(direction) == 0)
  this%hasDomainBoundary(2) =                                                                &
       (processCoordinates(direction) == processDistribution(direction) - 1)
  this%hasDomainBoundary(1) = this%hasDomainBoundary(1) .and. .not. isPeriodic(direction)
  this%hasDomainBoundary(2) = this%hasDomainBoundary(2) .and. .not. isPeriodic(direction)

  ! Copy the communicator and direction.
  if (this%cartesianCommunicator /= MPI_COMM_NULL)                                           &
       call MPI_Comm_free(this%cartesianCommunicator, ierror)
  this%cartesianCommunicator = cartesianCommunicator
  this%direction = direction

  ! Set the number of ghost points.
  this%nGhost = this%interiorWidth / 2
  if (.not. isPeriodic(direction) .and. processCoordinates(direction) == 0)                  &
       this%nGhost(1) = 0
  if (.not. isPeriodic(direction) .and. processCoordinates(direction) ==                     &
       processDistribution(direction) - 1) this%nGhost(2) = 0

  ! Set the periodic offset for overlap-type periodicity.
  this%periodicOffset = 0
  if (isPeriodic(direction) .and. present(isPeriodicityOverlapping)) then
     if (isPeriodicityOverlapping .and. processCoordinates(direction) == 0)                  &
          this%periodicOffset(2) = 1
     if (isPeriodicityOverlapping .and. processCoordinates(direction) ==                     &
          processDistribution(direction) - 1) this%periodicOffset(1) = 1
  end if

end subroutine updateOperator

subroutine cleanupOperator(this)

  ! <<< Derived types >>>
  use StencilOperator_type

  implicit none

  ! <<< Arguments >>>
  type(t_StencilOperator) :: this

  SAFE_DEALLOCATE(this%rhsInterior)
  SAFE_DEALLOCATE(this%rhsBoundary1)
  SAFE_DEALLOCATE(this%rhsBoundary2)
  SAFE_DEALLOCATE(this%normBoundary)

end subroutine cleanupOperator

subroutine getAdjointOperator(this, adjointOperator)

  ! Sets up the adjoint of a stencil operator.

  ! <<< Derived types >>>
  use StencilOperator_type

  ! <<< Private members >>>
  use StencilOperatorImpl, only : allocateData

  implicit none

  ! <<< Arguments >>>
  type(t_StencilOperator) :: this, adjointOperator

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j

  assert(this%symmetryType == SYMMETRIC .or. this%symmetryType == SKEW_SYMMETRIC)
  assert(this%interiorWidth > 0)
  assert(this%boundaryWidth > 0)
  assert(this%boundaryDepth > 0)

  ! Copy basic information.
  adjointOperator%symmetryType = this%symmetryType
  adjointOperator%interiorWidth = this%interiorWidth
  adjointOperator%boundaryDepth = this%boundaryWidth
  adjointOperator%boundaryWidth = this%boundaryWidth + this%interiorWidth / 2
  call allocateData(adjointOperator)

  ! Reverse the interior stencil.
  assert(allocated(this%rhsInterior))
  assert(lbound(this%rhsInterior) == - this%interiorWidth / 2)
  assert(ubound(this%rhsInterior) == + this%interiorWidth / 2)
  do i = - this%interiorWidth / 2, this%interiorWidth / 2
     adjointOperator%rhsInterior(i) = this%rhsInterior(-i)
  end do

  ! Copy `normBoundary`.
  assert(allocated(this%normBoundary))
  assert(size(this%normBoundary) == this%boundaryDepth)
  assert(all(this%normBoundary > 0))
  SAFE_DEALLOCATE(adjointOperator%normBoundary)
  allocate(adjointOperator%normBoundary(this%boundaryDepth))
  adjointOperator%normBoundary = this%normBoundary

  ! Transpose the left-boundary coefficients.
  assert(allocated(this%rhsBoundary1))
  assert(size(this%rhsBoundary1, 1) == this%boundaryWidth)
  assert(size(this%rhsBoundary1, 2) == this%boundaryDepth)
  adjointOperator%rhsBoundary1 = 0.0_wp
  adjointOperator%rhsBoundary1(1:this%boundaryDepth,:) = transpose(this%rhsBoundary1)
  do i = this%boundaryDepth + 1, this%boundaryWidth + this%interiorWidth / 2
     do j = - this%interiorWidth / 2, this%interiorWidth / 2
        if (i + j > this%boundaryWidth) exit
        adjointOperator%rhsBoundary1(i,i+j) = this%rhsInterior(j)
     end do
  end do

  ! Pre-multiply by the inverse of the norm matrix.
  do i = 1, adjointOperator%boundaryDepth
     adjointOperator%rhsBoundary1(1:size(adjointOperator%normBoundary),i) =                  &
          adjointOperator%rhsBoundary1(1:size(adjointOperator%normBoundary),i) *             &
          adjointOperator%normBoundary
  end do

  ! Post-multiply by the norm matrix.
  do i = 1, adjointOperator%boundaryWidth
     adjointOperator%rhsBoundary1(i,1:size(adjointOperator%normBoundary)) =                  &
          adjointOperator%rhsBoundary1(i,1:size(adjointOperator%normBoundary)) /             &
          adjointOperator%normBoundary
  end do

  ! Fill the right-boundary coefficients.
  select case (adjointOperator%symmetryType)
  case (SYMMETRIC)
     adjointOperator%rhsBoundary2(1:adjointOperator%boundaryWidth,:) =                       &
          +adjointOperator%rhsBoundary1(adjointOperator%boundaryWidth:1:-1,:)
  case (SKEW_SYMMETRIC)
     adjointOperator%rhsBoundary2(1:adjointOperator%boundaryWidth,:) =                       &
          -adjointOperator%rhsBoundary1(adjointOperator%boundaryWidth:1:-1,:)
  end select

end subroutine getAdjointOperator

subroutine applyOperator(this, x, gridSize)

  ! <<< Derived types >>>
  use StencilOperator_type

  ! <<< Private members >>>
  use StencilOperatorImpl, only : applyOperator_1, applyOperator_2, applyOperator_3

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  type(t_StencilOperator) :: this
  SCALAR_TYPE, intent(inout) :: x(:,:)
  integer, intent(in) :: gridSize(3)

  call startTiming("applyOperator")

  assert(size(x, 2) > 0)
  assert(all(gridSize > 0))
  assert(size(x, 1) == product(gridSize))
  assert(this%direction >= 1 .and. this%direction <= 3)

  select case (this%direction)
  case (1)
     call applyOperator_1(this, x, gridSize)
  case (2)
     call applyOperator_2(this, x, gridSize)
  case (3)
     call applyOperator_3(this, x, gridSize)
  end select

  call endTiming("applyOperator")

end subroutine applyOperator

subroutine applyOperatorAtInteriorPoints(this, xWithGhostPoints, x, gridSize)

  ! <<< Derived types >>>
  use StencilOperator_type

  ! <<< Private members >>>
  use StencilOperatorImpl, only : applyOperatorAtInteriorPoints_1,                           &
                                  applyOperatorAtInteriorPoints_2,                           &
                                  applyOperatorAtInteriorPoints_3

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  type(t_StencilOperator) :: this
  SCALAR_TYPE, intent(in) :: xWithGhostPoints(:,:,:,:)
  SCALAR_TYPE, intent(out) :: x(:,:)
  integer, intent(in) :: gridSize(3)

  call startTiming("applyOperatorAtInteriorPoints")

  assert(size(x, 2) > 0)
  assert(size(xWithGhostPoints, 4) == size(x, 2))
  assert(all(gridSize > 0))
  assert(size(x, 1) == product(gridSize))
  assert(this%direction >= 1 .and. this%direction <= 3)

  select case (this%direction)
  case (1)
     assert(size(xWithGhostPoints, 1) >= gridSize(1))
     assert(size(xWithGhostPoints, 2) == gridSize(2))
     assert(size(xWithGhostPoints, 3) == gridSize(3))
     call applyOperatorAtInteriorPoints_1(this, xWithGhostPoints, x, gridSize)
  case (2)
     assert(size(xWithGhostPoints, 1) == gridSize(1))
     assert(size(xWithGhostPoints, 2) >= gridSize(2))
     assert(size(xWithGhostPoints, 3) == gridSize(3))
     call applyOperatorAtInteriorPoints_2(this, xWithGhostPoints, x, gridSize)
  case (3)
     assert(size(xWithGhostPoints, 1) == gridSize(1))
     assert(size(xWithGhostPoints, 2) == gridSize(2))
     assert(size(xWithGhostPoints, 3) >= gridSize(3))
     call applyOperatorAtInteriorPoints_3(this, xWithGhostPoints, x, gridSize)
  end select

  call endTiming("applyOperatorAtInteriorPoints")

end subroutine applyOperatorAtInteriorPoints

subroutine applyOperatorNorm(this, x, gridSize)

  ! <<< Derived types >>>
  use StencilOperator_type

  ! <<< Private members >>>
  use StencilOperatorImpl, only : applyOperatorNorm_1,                                       &
                                  applyOperatorNorm_2,                                       &
                                  applyOperatorNorm_3

  implicit none

  ! <<< Arguments >>>
  type(t_StencilOperator) :: this
  SCALAR_TYPE, intent(inout) :: x(:,:)
  integer, intent(in) :: gridSize(3)

  assert(size(x, 2) > 0)
  assert(all(gridSize > 0))
  assert(size(x, 1) == product(gridSize))
  assert(this%direction >= 1 .and. this%direction <= 3)

  select case (this%direction)
  case (1)
     call applyOperatorNorm_1(this, x, gridSize)
  case (2)
     call applyOperatorNorm_2(this, x, gridSize)
  case (3)
     call applyOperatorNorm_3(this, x, gridSize)
  end select

end subroutine applyOperatorNorm

subroutine applyOperatorNormInverse(this, x, gridSize)

  ! <<< Derived types >>>
  use StencilOperator_type

  ! <<< Private members >>>
  use StencilOperatorImpl, only : applyOperatorNormInverse_1,                                &
                                  applyOperatorNormInverse_2,                                &
                                  applyOperatorNormInverse_3

  implicit none

  ! <<< Arguments >>>
  type(t_StencilOperator) :: this
  SCALAR_TYPE, intent(inout) :: x(:,:)
  integer, intent(in) :: gridSize(3)

  assert(size(x, 2) > 0)
  assert(all(gridSize > 0))
  assert(size(x, 1) == product(gridSize))
  assert(this%direction >= 1 .and. this%direction <= 3)

  select case (this%direction)
  case (1)
     call applyOperatorNormInverse_1(this, x, gridSize)
  case (2)
     call applyOperatorNormInverse_2(this, x, gridSize)
  case (3)
     call applyOperatorNormInverse_3(this, x, gridSize)
  end select

end subroutine applyOperatorNormInverse
