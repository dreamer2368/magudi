#include "config.h"

subroutine setupRandomFluctuationSource(this,amplitude)

  ! <<< Derived types >>>
  use RandomFluctuationSource_mod, only : t_RandomFluctuationSource

  implicit none

  ! <<< Arguments >>>
  class(t_RandomFluctuationSource) :: this
  real(SCALAR_KIND), intent(in) :: amplitude

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  this%amplitude = amplitude

end subroutine setupRandomFluctuationSource

subroutine addRandomFluctuationSource(this, time, coordinates, iblank, rightHandSide)

  ! <<< Derived types >>>
  use RandomFluctuationSource_mod, only : t_RandomFluctuationSource

  implicit none

  ! <<< Arguments >>>
  class(t_RandomFluctuationSource) :: this
  real(SCALAR_KIND), intent(in) :: time
  SCALAR_TYPE, intent(in) :: coordinates(:,:)
  integer, intent(in) :: iblank(:)
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions
  real(wp) :: r,fluctuation

  nDimensions = size(coordinates, 2)
  assert_key(nDimensions, (1, 2, 3))

  do i = 1, size(rightHandSide, 1)
     if (iblank(i) == 0) cycle
          do j=1,size(rightHandSide,2)
          call random_number(r)
          fluctuation = this%amplitude*r
          rightHandSide(i,j) = rightHandSide(i,j)+fluctuation
          end do
  end do

end subroutine addRandomFluctuationSource
