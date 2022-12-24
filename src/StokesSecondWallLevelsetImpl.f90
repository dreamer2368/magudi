#include "config.h"

subroutine setupStokesSecondWallLevelset(this, grids, states)

  ! <<< Derived types >>>
  use StokesSecondWallLevelset_mod, only : t_StokesSecondWallLevelset
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State

  ! <<< Internal modules >>>
  use InputHelper, only : getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_StokesSecondWallLevelset) :: this
  class(t_Grid), intent(in) :: grids(:)
  class(t_State) :: states(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  integer :: i
  real(wp) :: buf
  character(len = STRING_LENGTH) :: key

  call this%cleanup()

  assert(size(states) == size(grids))
  allocate(this%levelsetLoc(grids(1)%nDimensions))
  allocate(this%levelsetNormal(grids(1)%nDimensions))
  allocate(this%levelsetDirection(grids(1)%nDimensions))

  do i = 1, grids(1)%nDimensions
    write(key, '(A,I1.1)') "immersed_boundary/location", i
    call getRequiredOption(key, this%levelsetLoc(i), grids(1)%comm)

    write(key, '(A,I1.1)') "immersed_boundary/normal", i
    call getRequiredOption(key, this%levelsetNormal(i), grids(1)%comm)

    write(key, '(A,I1.1)') "immersed_boundary/oscillation_direction", i
    call getRequiredOption(key, this%levelsetDirection(i), grids(1)%comm)
  end do

  ! normalize normal.
  buf = 0.0_wp
  do i = 1, grids(1)%nDimensions
    buf = buf + this%levelsetNormal(i) * this%levelsetNormal(i)
  end do
  this%levelsetNormal = this%levelsetNormal / sqrt(buf)

  ! normalize direction.
  buf = 0.0_wp
  do i = 1, grids(1)%nDimensions
    buf = buf + this%levelsetDirection(i) * this%levelsetDirection(i)
  end do
  this%levelsetDirection = this%levelsetDirection / sqrt(buf)

  call getRequiredOption("immersed_boundary/amplitude", this%levelsetAmp, grids(1)%comm)
  call getRequiredOption("immersed_boundary/period", this%levelsetPeriod, grids(1)%comm)

end subroutine setupStokesSecondWallLevelset

subroutine cleanupStokesSecondWallLevelset(this)

  ! <<< Derived types >>>
  use StokesSecondWallLevelset_mod, only : t_StokesSecondWallLevelset

  ! <<< Arguments >>>
  class(t_StokesSecondWallLevelset) :: this

  SAFE_DEALLOCATE(this%levelsetLoc)
  SAFE_DEALLOCATE(this%levelsetNormal)
  SAFE_DEALLOCATE(this%levelsetDirection)

end subroutine cleanupStokesSecondWallLevelset

subroutine updateStokesSecondWallLevelset(this, mode, grids, states)

  ! <<< Derived types >>>
  use StokesSecondWallLevelset_mod, only : t_StokesSecondWallLevelset
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_StokesSecondWallLevelset) :: this
  integer, intent(in) :: mode
  class(t_Grid), intent(in) :: grids(:)
  class(t_State) :: states(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  real(wp) :: timeFactor, timeDerivativeFactor
  real(wp), dimension(grids(1)%nDimensions) :: loc, vel
  integer :: i, j

  call startTiming("updateStokesSecondWallLevelset")

  if (mode .ne. FORWARD) then
    return
  end if

  timeFactor = cos(2.0_wp * pi * states(1)%time / this%levelsetPeriod)
  timeDerivativeFactor = -2.0_wp * pi / this%levelsetPeriod                         &
                         * sin(2.0_wp * pi * states(1)%time / this%levelsetPeriod)
  loc = this%levelsetLoc + timeFactor * this%levelsetAmp * this%levelsetDirection
  vel = timeDerivativeFactor * this%levelsetAmp * this%levelsetDirection

  do i = 1, size(states)
    !NOTE: you cannot simply cycle here, since even in the same grid, some
    !processors may not have ibm patch.
    !if (.not. states(i)%ibmPatchExists) cycle

    assert(size(states(i)%levelset, 1) == size(states(i)%conservedVariables, 1))
    assert(size(states(i)%levelset, 2) == 1)

    states(i)%levelset = 0.0_wp
    do j = 1, grids(i)%nDimensions
      states(i)%levelset(:, 1) = states(i)%levelset(:, 1)                       &
        + this%levelsetNormal(j) * (grids(i)%coordinates(:, j) - loc(j))

      states(i)%levelsetNormal(:, j) = this%levelsetNormal(j)
      states(i)%objectVelocity(:, j) = vel(j)
    end do
  end do

  call endTiming("updateStokesSecondWallLevelset")

end subroutine updateStokesSecondWallLevelset
