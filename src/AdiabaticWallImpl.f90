#include "config.h"

subroutine setupAdiabaticWall(this, index, comm, patchDescriptor,                            &
     grid, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use AdiabaticWall_mod, only : t_AdiabaticWall
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_AdiabaticWall) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: key
  integer :: i

  call this%cleanup()
  call this%t_ImpenetrableWall%setup(index, comm, patchDescriptor,                           &
       grid, simulationFlags, solverOptions)

  write(key, '(A,I0.0)') "patches/" // trim(patchDescriptor%name) // "/"

  ! Viscous penalty amounts.
  if (simulationFlags%viscosityOn) then
     do i = 1, 2
        write(key, '(A,I0.0)') "patches/" // trim(patchDescriptor%name) //                   &
             "/viscous_penalty_amount", i
        this%viscousPenaltyAmounts(i) = getOption(trim(key), 1.0_wp)
        this%viscousPenaltyAmounts(i) = this%viscousPenaltyAmounts(i) *                      &
             solverOptions%reynoldsNumberInverse
        this%viscousPenaltyAmounts(i) = sign(this%viscousPenaltyAmounts(i),                  &
             real(this%normalDirection, wp))
        this%viscousPenaltyAmounts(i) = this%viscousPenaltyAmounts(i) /                      &
             grid%firstDerivative(abs(this%normalDirection))%normBoundary(1)
     end do
  else
     this%viscousPenaltyAmounts = 0.0_wp
  end if

end subroutine setupAdiabaticWall

subroutine cleanupAdiabaticWall(this)

  ! <<< Derived types >>>
  use AdiabaticWall_mod, only : t_AdiabaticWall

  implicit none

  ! <<< Arguments >>>
  class(t_AdiabaticWall) :: this

  call this%t_ImpenetrableWall%cleanup()

end subroutine cleanupAdiabaticWall

subroutine addAdiabaticWallPenalty(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use AdiabaticWall_mod, only : t_AdiabaticWall
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use CNSHelper
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_AdiabaticWall) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, nUnknowns, direction, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: metricsAlongNormalDirection(:), viscousPenalties(:,:)

  assert_key(mode, (FORWARD, ADJOINT))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  call startTiming("addAdiabaticWallPenalty")

  call this%t_ImpenetrableWall%updateRhs(mode, simulationFlags, solverOptions, grid, state)

  if (.not. simulationFlags%viscosityOn) then
     call endTiming("addAdiabaticWallPenalty")
     return
  end if

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2)

  allocate(metricsAlongNormalDirection(nDimensions))
  allocate(viscousPenalties(nUnknowns - 1, 2))

  do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
     do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
        do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
           gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                      &
                (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                        &
                (k - 1 - this%gridOffset(3)))
           if (grid%iblank(gridIndex) == 0) cycle
           patchIndex = i - this%offset(1) + this%localSize(1) *                             &
                (j - 1 - this%offset(2) + this%localSize(2) *                                &
                (k - 1 - this%offset(3)))

           metricsAlongNormalDirection =                                                     &
                grid%metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

           viscousPenalties(:,:) = 0.0_wp !... TODO: implement an adiabatic wall SAT

           viscousPenalties = grid%jacobian(gridIndex, 1) * viscousPenalties

           select case (mode)

           case (FORWARD)

              state%rightHandSide(gridIndex,2:nUnknowns) =                                   &
                   state%rightHandSide(gridIndex,2:nUnknowns) -                              &
                   this%viscousPenaltyAmounts(1) * viscousPenalties(:,1) +                   &
                   this%viscousPenaltyAmounts(2) * viscousPenalties(:,2)

           case (ADJOINT)

              ! TODO: add viscous wall penalties for adjoint variables.

           end select !... mode

        end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  SAFE_DEALLOCATE(viscousPenalties)
  SAFE_DEALLOCATE(metricsAlongNormalDirection)

  call endTiming("addAdiabaticWallPenalty")

end subroutine addAdiabaticWallPenalty

function verifyAdiabaticWallUsage(this, patchDescriptor, gridSize, normalDirection,          &
     extent, simulationFlags, success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use AdiabaticWall_mod, only : t_AdiabaticWall
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_AdiabaticWall) :: this
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  integer, intent(in) :: gridSize(:), normalDirection, extent(6)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  logical, intent(out) :: success
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Result >>>
  logical :: isPatchUsed

  isPatchUsed = this%t_ImpenetrableWall%verifyUsage(patchDescriptor, gridSize,               &
       normalDirection, extent, simulationFlags, success, message)

end function verifyAdiabaticWallUsage
