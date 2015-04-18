#include "config.h"

subroutine setupSpongePatch(this, index, comm, patchDescriptor,                              &
     grid, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use SpongePatch_mod, only : t_SpongePatch
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_SpongePatch) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: key

  call this%cleanup()
  call this%setupBase(index, comm, patchDescriptor, grid, simulationFlags, solverOptions)

  if (this%nPatchPoints > 0) then
     allocate(this%spongeStrength(this%nPatchPoints), source = 0.0_wp)
  end if

  write(key, '(A)') "patches/" // trim(patchDescriptor%name) // "/"

  ! Sponge amount.
  this%spongeAmount = getOption("defaults/sponge_amount", 1.0_wp)
  this%spongeAmount = getOption(trim(key) // "sponge_amount", this%spongeAmount)

  ! Sponge exponent.
  this%spongeExponent = getOption("defaults/sponge_exponent", 2)
  this%spongeExponent = getOption(trim(key) // "sponge_exponent", this%spongeExponent)

end subroutine setupSpongePatch

subroutine cleanupSpongePatch(this)

  ! <<< Derived types >>>
  use SpongePatch_mod, only : t_SpongePatch

  implicit none

  ! <<< Arguments >>>
  class(t_SpongePatch) :: this

  call this%cleanupBase()

  SAFE_DEALLOCATE(this%spongeStrength)

end subroutine cleanupSpongePatch

subroutine addDamping(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SpongePatch_mod, only : t_SpongePatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_SpongePatch) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer :: i, j, k, l, gridIndex, patchIndex

  assert_key(mode, (FORWARD, ADJOINT))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  call startTiming("addDamping")

  select case (mode)

  case (FORWARD)
     do l = 1, solverOptions%nUnknowns
        do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
           do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
              do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                 gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                &
                      (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                  &
                      (k - 1 - this%gridOffset(3)))
                 if (grid%iblank(gridIndex) == 0) cycle
                 patchIndex = i - this%offset(1) + this%localSize(1) *                       &
                      (j - 1 - this%offset(2) + this%localSize(2) *                          &
                      (k - 1 - this%offset(3)))
                 state%rightHandSide(gridIndex, l) = state%rightHandSide(gridIndex, l) -     &
                      this%spongeStrength(patchIndex) *                                      &
                      (state%conservedVariables(gridIndex, l) -                              &
                      state%targetState(gridIndex, l))
              end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
           end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
        end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
     end do !... l = 1, size(state%rightHandSide, 2)

  case (ADJOINT)
     do l = 1, solverOptions%nUnknowns
        do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
           do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
              do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                 gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                &
                      (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                  &
                      (k - 1 - this%gridOffset(3)))
                 if (grid%iblank(gridIndex) == 0) cycle
                 patchIndex = i - this%offset(1) + this%localSize(1) *                       &
                      (j - 1 - this%offset(2) + this%localSize(2) *                          &
                      (k - 1 - this%offset(3)))
                 state%rightHandSide(gridIndex, l) = state%rightHandSide(gridIndex, l) +     &
                      this%spongeStrength(patchIndex) * state%adjointVariables(gridIndex, l)
              end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
           end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
        end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
     end do !... l = 1, solverOptions%nUnknowns

  end select

  call endTiming("addDamping")

end subroutine addDamping

function verifySpongePatchUsage(this, patchDescriptor, gridSize, normalDirection,            &
     extent, simulationFlags, success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use SpongePatch_mod, only : t_SpongePatch
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_SpongePatch) :: this
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  integer, intent(in) :: gridSize(:), normalDirection, extent(6)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  logical, intent(out) :: success
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Result >>>
  logical :: isPatchUsed

  ! <<< Local variables >>>
  integer :: i, n

  isPatchUsed = .false.

  success = .false.
  if (normalDirection > size(gridSize) .or. normalDirection == 0) then
     write(message, '(A)') "Normal direction is invalid!"
     return
  end if

  if (.not. simulationFlags%useTargetState) then
     write(message, '(A)') "No target state available for damping!"
     return
  end if

  n = size(gridSize)

  do i = 1, size(gridSize)
     if (extent((i-1)*2+1) < 0 .or. extent((i-1)*2+2) > gridSize(i) .or.                     &
          extent((i-1)*2+1) > extent((i-1)*2+2)) then
        write(message, '(A)') "Invalid extent!"
        return
     end if
     if (extent((i-1)*2+1) == extent((i-1)*2+2)) n = n - 1
  end do

  if (n /= size(gridSize)) then
     write(message, '(2(A,I0.0),A)') "Expected a ", size(gridSize),                          &
          "D patch, but extent represents a ", n, "D patch!"
     return
  end if

  i = abs(normalDirection)
  if ((normalDirection > 0 .and. extent((i-1)*2+1) /= 1) .or.                                &
       (normalDirection < 0 .and. extent((i-1)*2+2) /= gridSize(i))) then
     write(message, '(2(A,I0.0),A)') "Not aligned with a computational boundary!"
     return
  end if

  success = .true.

  isPatchUsed = .true.

end function verifySpongePatchUsage
