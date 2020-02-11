#include "config.h"

module XmomentumConservingBodyForceImpl

  implicit none
  public

contains

  function computeXmomentum(region) result(xMomentum)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! ! <<< SeungWhan: debugging >>>
    ! use, intrinsic :: iso_fortran_env, only : output_unit
    ! use ErrorHandler, only : writeAndFlush

    ! <<< Arguments >>>
    class(t_Region), intent(in) :: region

    ! <<< Result >>>
    SCALAR_TYPE :: xMomentum

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, ierror
    SCALAR_TYPE, allocatable :: F(:,:)

    ! ! <<< SeungWhan: message, timeRampFactor >>
    ! character(len=STRING_LENGTH) :: message
    ! real(wp) :: timeRampFactor

    assert(allocated(region%grids))
    assert(allocated(region%states))
    assert(size(region%grids) == size(region%states))

    xMomentum = 0.0_wp

    do i = 1, size(region%grids)

       assert(region%grids(i)%nGridPoints > 0)
       assert(allocated(region%states(i)%conservedVariables))
       assert(size(region%states(i)%conservedVariables, 1) == region%grids(i)%nGridPoints)
       assert(size(region%states(i)%conservedVariables, 2) > 2)

       allocate(F(region%grids(i)%nGridPoints, 2))
       F(:,1) = region%states(i)%conservedVariables(:,2)
       F(:,2) = 1.0_wp
       xMomentum = xMomentum + region%grids(i)%computeInnerProduct(F(:,1),F(:,2))
       SAFE_DEALLOCATE(F)
    end do

    if (region%commGridMasters /= MPI_COMM_NULL)                                               &
         call MPI_Allreduce(MPI_IN_PLACE, xMomentum, 1,                                        &
         SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

    do i = 1, size(region%grids)
       call MPI_Bcast(xMomentum, 1, SCALAR_TYPE_MPI,                                           &
            0, region%grids(i)%comm, ierror)
    end do

  end function computeXmomentum

  function computeVolume(region) result(volume)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! ! <<< SeungWhan: debugging >>>
    ! use, intrinsic :: iso_fortran_env, only : output_unit
    ! use ErrorHandler, only : writeAndFlush

    ! <<< Arguments >>>
    class(t_Region), intent(in) :: region

    ! <<< Result >>>
    SCALAR_TYPE :: volume

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, ierror
    SCALAR_TYPE, allocatable :: F(:,:)

    ! ! <<< SeungWhan: message, timeRampFactor >>
    ! character(len=STRING_LENGTH) :: message
    ! real(wp) :: timeRampFactor

    assert(allocated(region%grids))

    volume = 0.0_wp

    do i = 1, size(region%grids)

       assert(region%grids(i)%nGridPoints > 0)

       allocate(F(region%grids(i)%nGridPoints, 1))
       F(:,1) = 1.0_wp
       volume = volume + region%grids(i)%computeInnerProduct(F(:,1),F(:,1))
       SAFE_DEALLOCATE(F)
    end do

    if (region%commGridMasters /= MPI_COMM_NULL)                                               &
         call MPI_Allreduce(MPI_IN_PLACE, volume, 1,                                        &
         SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

    do i = 1, size(region%grids)
       call MPI_Bcast(volume, 1, SCALAR_TYPE_MPI,                                           &
            0, region%grids(i)%comm, ierror)
    end do

  end function computeVolume

end module XmomentumConservingBodyForceImpl

subroutine setupBodyForce(this, region)

  ! <<< Derived types >>>
  use XmomentumConservingBodyForce_mod, only : t_BodyForce
  use Region_mod, only : t_Region

  ! <<< Private members >>>
  use XmomentumConservingBodyForceImpl

  implicit none

  ! <<< Arguments >>>
  class(t_BodyForce) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)

  this%initialXmomentum = computeXmomentum(region)
  this%oneOverVolume = 1.0_wp/computeVolume(region)

end subroutine setupBodyForce

subroutine addBodyForce(this, region)

  ! <<< Derived types >>>
  use XmomentumConservingBodyForce_mod, only : t_BodyForce
  use Region_mod, only : t_Region

  ! <<< Private members >>>
  use XmomentumConservingBodyForceImpl

  use, intrinsic :: iso_fortran_env, only : output_unit
  use ErrorHandler, only : writeAndFlush

  ! <<< Internal modules >>>
  use CNSHelper
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_BodyForce) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions
  SCALAR_TYPE :: currentXmomentum
  SCALAR_TYPE, allocatable :: velocity(:)

  character(len=STRING_LENGTH) :: message

  call startTiming("addBodyForce")

  assert(allocated(region%states))
  assert(allocated(region%grids))
  nDimensions = region%grids(1)%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  currentXmomentum = computeXmomentum(region)
  this%momentumLossPerVolume = this%oneOverVolume *                             &
                                ( this%initialXmomentum - currentXmomentum )

  do i = 1, size(region%states)
    assert(region%grids(i)%nGridPoints > 0)
    allocate(velocity(region%grids(i)%nGridPoints))
    velocity = region%states(i)%conservedVariables(:,2) /                       &
                  region%states(i)%conservedVariables(:,1)

    region%states(i)%conservedVariables(:,2) =                                  &
      region%states(i)%conservedVariables(:,2) + this%momentumLossPerVolume

    region%states(i)%conservedVariables(:,nDimensions+2) =                      &
      region%states(i)%conservedVariables(:,nDimensions+2) +                    &
      this%momentumLossPerVolume * velocity

    SAFE_DEALLOCATE(velocity)
  end do

  call endTiming("addBodyForce")

end subroutine addBodyForce
