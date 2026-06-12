#include "config.h"

subroutine setupCylinderLevelset(this, grids, states)

  ! <<< Derived types >>>
  use CylinderLevelset_mod, only : t_CylinderLevelset
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State

  ! <<< Internal modules >>>
  use InputHelper, only : getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_CylinderLevelset) :: this
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
  assert(grids(1)%nDimensions > 1)
  allocate(this%loc(grids(1)%nDimensions))
  allocate(this%direction(grids(1)%nDimensions))

  do i = 1, grids(1)%nDimensions
    write(key, '(A,I1.1)') "immersed_boundary/location", i
    call getRequiredOption(key, this%loc(i), grids(1)%comm)

    write(key, '(A,I1.1)') "immersed_boundary/oscillation_direction", i
    call getRequiredOption(key, this%direction(i), grids(1)%comm)
  end do

  ! normalize direction.
  buf = 0.0_wp
  do i = 1, grids(1)%nDimensions
    buf = buf + this%direction(i) * this%direction(i)
  end do
  this%direction = this%direction / sqrt(buf)

  call getRequiredOption("immersed_boundary/radius", this%radius, grids(1)%comm)
  call getRequiredOption("immersed_boundary/amplitude", this%amplitude, grids(1)%comm)
  call getRequiredOption("immersed_boundary/period", this%period, grids(1)%comm)

end subroutine setupCylinderLevelset

subroutine cleanupCylinderLevelset(this)

  ! <<< Derived types >>>
  use CylinderLevelset_mod, only : t_CylinderLevelset

  ! <<< Arguments >>>
  class(t_CylinderLevelset) :: this

  SAFE_DEALLOCATE(this%loc)
  SAFE_DEALLOCATE(this%direction)

end subroutine cleanupCylinderLevelset

subroutine updateCylinderLevelset(this, mode, grids, states)

  ! <<< Derived types >>>
  use CylinderLevelset_mod, only : t_CylinderLevelset
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_CylinderLevelset) :: this
  integer, intent(in) :: mode
  class(t_Grid), intent(in) :: grids(:)
  class(t_State) :: states(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  real(wp) :: timeFactor, timeDerivativeFactor, timeAccFactor, tmp
  real(wp), dimension(grids(1)%nDimensions) :: loc, vel, acc
  integer :: i, j

  call startTiming("updateCylinderLevelset")

  if (mode .ne. FORWARD) then
    return
  end if

  timeFactor = cos(2.0_wp * pi * states(1)%time / this%period)
  timeDerivativeFactor = -2.0_wp * pi / this%period                         &
                         * sin(2.0_wp * pi * states(1)%time / this%period)
  timeAccFactor = timeFactor * (-4.0_wp) * pi * pi / this%period / this%period
  loc = this%loc + timeFactor * this%amplitude * this%direction
  vel = timeDerivativeFactor * this%amplitude * this%direction
  acc = timeAccFactor * this%amplitude * this%direction

  do i = 1, size(states)
    !NOTE: you cannot simply cycle here, since even in the same grid, some
    !processors may not have ibm patch.
    !if (.not. states(i)%ibmPatchExists) cycle

    assert(size(states(i)%levelset, 1) == size(states(i)%conservedVariables, 1))
    assert(size(states(i)%levelset, 2) == 1)

    states(i)%levelset = 0.0_wp
    states(i)%levelsetNormal = 0.0_wp
    do j = 1, 2
      states(i)%levelset(:, 1) = states(i)%levelset(:, 1)                       &
        + (grids(i)%coordinates(:, j) - loc(j)) ** 2

      states(i)%levelsetNormal(:, j) = grids(i)%coordinates(:, j) - loc(j)
      states(i)%levelset(:, 1) = states(i)%levelset(:, 1) +                     &
                                   (states(i)%levelsetNormal(:, j)) ** 2
    end do

    do j = 1, grids(i)%nGridPoints
      tmp = max(sqrt(states(i)%levelset(j, 1)), 0.25_wp * this%radius)
      states(i)%levelsetNormal(j, :) = states(i)%levelsetNormal(j, :) / tmp
    end do
    states(i)%levelset = states(i)%levelset - this%radius * this%radius

    do j = 1, grids(i)%nDimensions
      states(i)%objectVelocity(:, j) = vel(j)
      states(i)%objectAcceleration(:, j) = acc(j)
    end do
  end do

  call endTiming("updateCylinderLevelset")

end subroutine updateCylinderLevelset
