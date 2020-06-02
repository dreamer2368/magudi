#include "config.h"

subroutine cleanupRegionVector(this)

  ! <<< Derived types >>>
  use RegionVector_mod, only : t_RegionVector

  implicit none

  ! <<< Arguments >>>
  class(t_RegionVector) :: this

  ! <<< Local variables >>>
  integer :: i

  if (allocated(this%states)) then
    do i = 1, size(this%states)
      SAFE_DEALLOCATE(this%states(i)%conservedVariables)
    end do
    deallocate(this%states)
  end if
  SAFE_DEALLOCATE(this%params)

end subroutine cleanupRegionVector

function addRegionVector(a,b) result(r)

  ! <<< Derived types >>>
  use RegionVector_mod, only : t_RegionVector

  implicit none

  ! <<< Arguments >>>
  class(t_RegionVector), intent(in) :: a, b
  type(t_RegionVector) :: r

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, bufferSize

  assert(allocated(a%states))
  assert(allocated(a%params))
  assert(allocated(b%states))
  assert(allocated(b%params))
  assert(allocated(r%states))
  assert(allocated(r%params))
  assert(size(a%states)==size(b%states))
  assert(size(a%params)==size(b%params))
  assert(size(r%states)==size(a%states))
  assert(size(r%params)==size(a%params))

  r%params = a%params + b%params

  do i = 1, size(a%states)
    assert(allocated(a%states(i)%conservedVariables))
    assert(allocated(b%states(i)%conservedVariables))
    assert(allocated(r%states(i)%conservedVariables))

    r%states(i)%conservedVariables = a%states(i)%conservedVariables             &
                                    + b%states(i)%conservedVariables
  end do

end function addRegionVector

function subtractRegionVector(a,b) result(r)

  ! <<< Derived types >>>
  use RegionVector_mod, only : t_RegionVector

  implicit none

  ! <<< Arguments >>>
  class(t_RegionVector), intent(in) :: a, b
  type(t_RegionVector) :: r

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, bufferSize

  assert(allocated(a%states))
  assert(allocated(a%params))
  assert(allocated(b%states))
  assert(allocated(b%params))
  assert(allocated(r%states))
  assert(allocated(r%params))
  assert(size(a%states)==size(b%states))
  assert(size(a%params)==size(b%params))
  assert(size(r%states)==size(a%states))
  assert(size(r%params)==size(a%params))

  r%params = a%params - b%params

  do i = 1, size(a%states)
    assert(allocated(a%states(i)%conservedVariables))
    assert(allocated(b%states(i)%conservedVariables))
    assert(allocated(r%states(i)%conservedVariables))

    r%states(i)%conservedVariables = a%states(i)%conservedVariables             &
                                    - b%states(i)%conservedVariables
  end do

end function subtractRegionVector

function multiplyRegionVector(a,b) result(r)

  ! <<< Derived types >>>
  use RegionVector_mod, only : t_RegionVector

  implicit none

  ! <<< Arguments >>>
  SCALAR_TYPE, intent(in) :: a
  class(t_RegionVector), intent(in) :: b
  type(t_RegionVector) :: r

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, bufferSize

  assert(allocated(b%states))
  assert(allocated(b%params))
  assert(allocated(r%states))
  assert(allocated(r%params))
  assert(size(r%states)==size(b%states))
  assert(size(r%params)==size(b%params))

  r%params = a * b%params

  do i = 1, size(b%states)
    assert(allocated(b%states(i)%conservedVariables))
    assert(allocated(r%states(i)%conservedVariables))

    r%states(i)%conservedVariables = a * b%states(i)%conservedVariables
  end do

end function multiplyRegionVector

function innerProduct(a,b,region) result(r)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use RegionVector_mod, only : t_RegionVector
  use Region_mod, only : t_Region

  implicit none

  ! <<< Arguments >>>
  class(t_RegionVector), intent(in) :: a, b
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: r

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, ierror

  r = 0.0_wp

  assert(allocated(a%states))
  assert(allocated(a%params))
  assert(allocated(b%states))
  assert(allocated(b%params))

  do i = 1, size(region%grids)
    assert(allocated(a%states(i)%conservedVariables))
    assert(allocated(b%states(i)%conservedVariables))

    r = r + region%grids(i)%computeInnerProduct(a%states(i)%conservedVariables, &
                                                b%states(i)%conservedVariables)
  end do

  if (region%commGridMasters /= MPI_COMM_NULL)                                  &
       call MPI_Allreduce(MPI_IN_PLACE, r, 1,                                   &
       SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(r, 1, SCALAR_TYPE_MPI, 0, region%grids(i)%comm, ierror)
  end do

  r = r + sum(a%params * b%params)

end function innerProduct
