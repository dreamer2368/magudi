#include "config.h"

subroutine setupActuatorPatch(this, index, comm, patchDescriptor,                            &
     grid, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_ActuatorPatch) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  call this%setupBase(index, comm, patchDescriptor, grid, simulationFlags, solverOptions)

  call this%cleanup()

  this%penaltyInPhysicalCoordinates = .true.
  
  if (.not. simulationFlags%predictionOnly .and. this%nPatchPoints > 0) then
     allocate(this%controlForcing(this%nPatchPoints, solverOptions%nUnknowns))
  end if

end subroutine setupActuatorPatch

subroutine cleanupActuatorPatch(this)

  ! <<< Derived types >>>
  use ActuatorPatch_mod, only : t_ActuatorPatch

  implicit none

  ! <<< Arguments >>>
  class(t_ActuatorPatch) :: this

  call this%cleanupBase()

  SAFE_DEALLOCATE(this%controlForcing)

end subroutine cleanupActuatorPatch

subroutine addControlForcing(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  implicit none

  ! <<< Arguments >>>
  class(t_ActuatorPatch) :: this
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

  if (mode == ADJOINT) return

  call startTiming("addActuatorPenalty")

  do l = 1, solverOptions%nUnknowns
     do k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
        do j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
           do i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
              gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                   &
                   (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                     &
                   (k - 1 - this%gridOffset(3)))
              if (grid%iblank(gridIndex) == 0) cycle
              patchIndex = i - this%offset(1) + this%patchSize(1) *                          &
                   (j - 1 - this%offset(2) + this%patchSize(2) *                             &
                   (k - 1 - this%offset(3)))

              ! state%rightHandSide(gridIndex,l) = state%rightHandSide(gridIndex,l) +          &
              !      grid%controlMollifier(gridIndex,1) * this%controlForcing(patchIndex,l)

           end do !... i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
        end do !... j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
     end do !... k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
  end do !... l = 1, solverOptions%nUnknowns

  call endTiming("addActuatorPenalty")

end subroutine addControlForcing

function verifyActuatorPatchUsage(this, patchDescriptor, gridSize, normalDirection,          &
     extent, simulationFlags, success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_ActuatorPatch) :: this
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  integer, intent(in) :: gridSize(:), normalDirection, extent(6)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  logical, intent(out) :: success
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Result >>>
  logical :: isPatchUsed

  ! <<< Local variables >>>
  integer :: i

  isPatchUsed = .false.

  success = .false.

  do i = 1, size(gridSize)
     if (extent((i-1)*2+1) < 0 .or. extent((i-1)*2+2) > gridSize(i) .or.                     &
          extent((i-1)*2+1) > extent((i-1)*2+2)) then
        write(message, '(A)') "Invalid extent!"
        return
     end if
  end do

  success = .true.

  isPatchUsed = .true.
  if (simulationFlags%predictionOnly) isPatchUsed = .false.

end function verifyActuatorPatchUsage

subroutine updateActuatorPatch(this, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_ActuatorPatch) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state

  ! Nothing to be done for this patch type...

end subroutine updateActuatorPatch
