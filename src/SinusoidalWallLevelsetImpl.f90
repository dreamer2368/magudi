#include "config.h"

subroutine setupSinusoidalWallLevelset(this, grids, states)

  ! <<< Derived types >>>
  use SinusoidalWallLevelset_mod, only : t_SinusoidalWallLevelset
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State

  ! <<< Internal modules >>>
  use InputHelper, only : getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_SinusoidalWallLevelset) :: this
  class(t_Grid), intent(in) :: grids(:)
  class(t_State) :: states(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  integer :: i, n

  call this%cleanup()

  assert(size(states) == size(grids))
  allocate(this%wallShapes(size(grids)))
  do i = 1, size(grids)
    allocate(this%wallShapes(i)%buffer(grids(i)%nGridPoints))
  end do

  call getRequiredOption("immersed_boundary/location", this%levelsetLoc, grids(1)%comm)
  call getRequiredOption("immersed_boundary/width", this%levelsetWidth, grids(1)%comm)
  call getRequiredOption("immersed_boundary/amplitude", this%levelsetAmp, grids(1)%comm)
  call getRequiredOption("immersed_boundary/period", this%levelsetPeriod, grids(1)%comm)

  do n = 1, size(grids)
    this%wallShapes(n)%buffer = 0.0_wp
    do i = 1, grids(n)%nGridPoints
      if ((grids(n)%coordinates(i, 1) < (this%levelsetLoc - 0.5_wp * this%levelsetWidth)) .or.  &
          (grids(n)%coordinates(i, 1) > (this%levelsetLoc + 0.5_wp * this%levelsetWidth))) cycle
      this%wallShapes(n)%buffer(i) = 0.5_wp + 0.5_wp                                              &
         * cos(2.0_wp * pi * (grids(n)%coordinates(i, 1) - this%levelsetLoc) / this%levelsetWidth)
    end do
    this%wallShapes(n)%buffer = this%wallShapes(n)%buffer * this%levelsetAmp
  end do

end subroutine setupSinusoidalWallLevelset

subroutine cleanupSinusoidalWallLevelset(this)

  ! <<< Derived types >>>
  use SinusoidalWallLevelset_mod, only : t_SinusoidalWallLevelset

  ! <<< Arguments >>>
  class(t_SinusoidalWallLevelset) :: this

  ! <<< Local variables >>>
  integer :: i

  if (allocated(this%wallShapes)) then
    do i = 1, size(this%wallShapes)
      SAFE_DEALLOCATE(this%wallShapes(i)%buffer)
    end do
    SAFE_DEALLOCATE(this%wallShapes)
  end if

end subroutine cleanupSinusoidalWallLevelset

subroutine updateSinusoidalWallLevelset(this, mode, grids, states)

  ! <<< Derived types >>>
  use SinusoidalWallLevelset_mod, only : t_SinusoidalWallLevelset
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_SinusoidalWallLevelset) :: this
  integer, intent(in) :: mode
  class(t_Grid), intent(in) :: grids(:)
  class(t_State) :: states(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  real(wp) :: timeFactor, timeDerivativeFactor, buf, objectSpeed
  integer :: i, j

  call startTiming("updateSinusoidalWallLevelset")

  if (mode .ne. FORWARD) then
    return
  end if

  timeFactor = sin(2.0_wp * pi * states(1)%time / this%levelsetPeriod)
  timeDerivativeFactor = 2.0_wp * pi / this%levelsetPeriod                          &
                         * cos(2.0_wp * pi * states(1)%time / this%levelsetPeriod)

  do i = 1, size(states)
    if (.not. states(i)%ibmPatchExists) cycle
    
    assert(size(states(i)%levelset, 1) == size(states(i)%conservedVariables, 1))
    assert(size(states(i)%levelset, 2) == 1)
    assert(size(states(i)%levelset, 1) == size(this%wallShapes(i)%buffer))

    states(i)%levelset(:, 1) = grids(i)%coordinates(:, 2)                       &
                               - this%wallShapes(i)%buffer * timeFactor

    ! compute levelset normal.
    call grids(i)%computeGradient(states(i)%levelset(:, 1), states(i)%levelsetNormal)

    do j = 1, grids(i)%nGridPoints
      ! levelset normal magnitude
      buf = sqrt(sum(states(i)%levelsetNormal(j,:) ** 2))

      ! d/dt levelset = - wallShape * d/dt timeFactor
      objectSpeed = timeDerivativeFactor * this%wallShapes(i)%buffer(j) / buf

      ! Make the levelset normal a unit norm
      if (buf .gt. 0.0_wp) states(i)%levelsetNormal(j,:) = states(i)%levelsetNormal(j,:) / buf

      states(i)%objectVelocity(j,:) = objectSpeed * states(i)%levelsetNormal(j,:)
    end do
  end do

  call endTiming("updateSinusoidalWallLevelset")

end subroutine updateSinusoidalWallLevelset
