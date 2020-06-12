#include "config.h"

module RegionVectorImpl

  implicit none
  public

contains

  subroutine saveRegionVector(this, region, mode)

    ! <<< Derived types >>>
    use Region_mod, only : t_Region
    use RegionVector_mod, only : t_RegionVector

    ! <<< Enumerations >>>
    use State_enum, only : QOI_ADJOINT_STATE, QOI_FORWARD_STATE,                &
                           QOI_RIGHT_HAND_SIDE

    implicit none

    ! <<< Arguments >>>
    class(t_RegionVector), intent(inout) :: this
    class(t_Region), intent(in) :: region
    integer, intent(in) :: mode

    ! <<< Local variables >>>
    integer :: i

    assert(size(this%states)==size(region%states))
    assert((size(this%params)==0) .or. (size(this%params)==size(region%params%buffer,1)))

    if (size(this%params)>0) this%params = region%params%buffer(:,1)

    do i = 1, size(this%states)
      select case (mode)
      case(QOI_FORWARD_STATE)
        assert(size(this%states(i)%conservedVariables,1)==size(region%states(i)%conservedVariables,1))
        assert(size(this%states(i)%conservedVariables,2)==size(region%states(i)%conservedVariables,2))

        this%states(i)%conservedVariables = region%states(i)%conservedVariables
      case(QOI_ADJOINT_STATE)
        assert(size(this%states(i)%conservedVariables,1)==size(region%states(i)%adjointVariables,1))
        assert(size(this%states(i)%conservedVariables,2)==size(region%states(i)%adjointVariables,2))

        this%states(i)%conservedVariables = region%states(i)%adjointVariables
      case(QOI_RIGHT_HAND_SIDE)
        assert(size(this%states(i)%conservedVariables,1)==size(region%states(i)%rightHandSide,1))
        assert(size(this%states(i)%conservedVariables,2)==size(region%states(i)%rightHandSide,2))

        this%states(i)%conservedVariables = region%states(i)%rightHandSide
      end select
    end do

  end subroutine saveRegionVector

  subroutine loadRegionVector(region, this, mode)

    ! <<< Derived types >>>
    use Region_mod, only : t_Region
    use RegionVector_mod, only : t_RegionVector

    ! <<< Enumerations >>>
    use State_enum, only : QOI_ADJOINT_STATE, QOI_FORWARD_STATE,                &
                           QOI_RIGHT_HAND_SIDE

    implicit none

    ! <<< Arguments >>>
    class(t_Region), intent(inout) :: region
    class(t_RegionVector), intent(in) :: this
    integer, intent(in) :: mode

    ! <<< Local variables >>>
    integer :: i

    assert(size(this%states)==size(region%states))
    assert((size(this%params)==0) .or. (size(this%params)==size(region%params%buffer,1)))

    if (size(this%params)>0) region%params%buffer(:,1) = this%params

    do i = 1, size(this%states)
      select case (mode)
      case(QOI_FORWARD_STATE)
        assert(size(this%states(i)%conservedVariables,1)==size(region%states(i)%conservedVariables,1))
        assert(size(this%states(i)%conservedVariables,2)==size(region%states(i)%conservedVariables,2))

        region%states(i)%conservedVariables = this%states(i)%conservedVariables
      case(QOI_ADJOINT_STATE)
        assert(size(this%states(i)%conservedVariables,1)==size(region%states(i)%adjointVariables,1))
        assert(size(this%states(i)%conservedVariables,2)==size(region%states(i)%adjointVariables,2))

        region%states(i)%adjointVariables = this%states(i)%conservedVariables
      case(QOI_RIGHT_HAND_SIDE)
        assert(size(this%states(i)%conservedVariables,1)==size(region%states(i)%rightHandSide,1))
        assert(size(this%states(i)%conservedVariables,2)==size(region%states(i)%rightHandSide,2))

        region%states(i)%rightHandSide = this%states(i)%conservedVariables
      end select
    end do

  end subroutine loadRegionVector

  function regionInnerProduct(a,b,region) result(r)

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

  end function regionInnerProduct

end module RegionVectorImpl

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

subroutine setFromRegionVector(this,x)

  ! <<< Derived types >>>
  use RegionVector_mod, only : t_RegionVector

  implicit none

  ! <<< Arguments >>>
  class(t_RegionVector), intent(inout) :: this
  class(t_RegionVector), intent(in) :: x

  ! <<< Local variables >>>
  integer :: i

  assert(allocated(x%states))
  assert(allocated(x%params))

  call this%cleanup()

  allocate(this%states(size(x%states)))
  allocate(this%params,MOLD=x%params)

  do i = 1, size(x%states)
    assert(allocated(x%states(i)%conservedVariables))
    allocate(this%states(i)%conservedVariables,MOLD=x%states(i)%conservedVariables)
  end do

end subroutine setFromRegionVector

subroutine setFromRegion(this,x)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use RegionVector_mod, only : t_RegionVector

  implicit none

  ! <<< Arguments >>>
  class(t_RegionVector), intent(inout) :: this
  class(t_Region), intent(in) :: x

  ! <<< Local variables >>>
  integer :: i, nParams

  assert(allocated(x%states))

  call this%cleanup()

  allocate(this%states(size(x%states)))
  nParams = 0
  if (allocated(x%params%buffer)) nParams = size(x%params%buffer,1)
  allocate(this%params(nParams))

  do i = 1, size(x%states)
    assert(allocated(x%states(i)%conservedVariables))
    allocate(this%states(i)%conservedVariables,MOLD=x%states(i)%conservedVariables)
  end do

end subroutine setFromRegion

subroutine copyRegionVector(b, a)

  ! <<< Derived types >>>
  use RegionVector_mod, only : t_RegionVector

  implicit none

  ! <<< Arguments >>>
  class(t_RegionVector), intent(out) :: b
  class(t_RegionVector), intent(in) :: a

  ! <<< Local variables >>>
  integer :: i

  call b%set(a)

  b%params = a%params

  do i = 1, size(b%states)
    assert(size(a%states(i)%conservedVariables,1)==size(b%states(i)%conservedVariables,1))
    assert(size(a%states(i)%conservedVariables,2)==size(b%states(i)%conservedVariables,2))

    b%states(i)%conservedVariables = a%states(i)%conservedVariables
  end do

end subroutine copyRegionVector

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
  assert(size(a%states)==size(b%states))
  assert(size(a%params)==size(b%params))

  call r%set(a)

  r%params = a%params + b%params

  do i = 1, size(a%states)
    assert(allocated(a%states(i)%conservedVariables))
    assert(allocated(b%states(i)%conservedVariables))

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
  assert(size(a%states)==size(b%states))
  assert(size(a%params)==size(b%params))

  call r%set(a)

  r%params = a%params - b%params

  do i = 1, size(a%states)
    assert(allocated(a%states(i)%conservedVariables))
    assert(allocated(b%states(i)%conservedVariables))

    r%states(i)%conservedVariables = a%states(i)%conservedVariables             &
                                    - b%states(i)%conservedVariables
  end do

end function subtractRegionVector

function multiplyRegionVector(a,b) result(r)

  ! <<< Derived types >>>
  use RegionVector_mod, only : t_RegionVector

  implicit none

  ! <<< Arguments >>>
  class(t_RegionVector), intent(in) :: a
  SCALAR_TYPE, intent(in) :: b
  type(t_RegionVector) :: r

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, bufferSize

  assert(allocated(a%states))
  assert(allocated(a%params))

  call r%set(a)

  r%params = b * a%params

  do i = 1, size(a%states)
    assert(allocated(a%states(i)%conservedVariables))

    r%states(i)%conservedVariables = b * a%states(i)%conservedVariables
  end do

end function multiplyRegionVector
