#include "config.h"

subroutine setupImpenetrableWall(this, index, comm, patchDescriptor,                         &
     grid, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags
  use ImpenetrableWall_mod, only : t_ImpenetrableWall

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_ImpenetrableWall) :: this
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

  write(key, '(A)') "patches/" // trim(patchDescriptor%name) // "/"

  ! Inviscid penalty amount.
  this%inviscidPenaltyAmount = getOption(trim(key) //                                        &
       "inviscid_penalty_amount", 1.0_wp) !... default value => dual-consistent.
  this%inviscidPenaltyAmount = sign(this%inviscidPenaltyAmount,                              &
       real(this%normalDirection, wp))
  this%inviscidPenaltyAmount = this%inviscidPenaltyAmount /                                  &
       grid%firstDerivative(abs(this%normalDirection))%normBoundary(1)

end subroutine setupImpenetrableWall

subroutine cleanupImpenetrableWall(this)

  ! <<< Derived types >>>
  use ImpenetrableWall_mod, only : t_ImpenetrableWall

  implicit none

  ! <<< Arguments >>>
  class(t_ImpenetrableWall) :: this

  call this%cleanupBase()

end subroutine cleanupImpenetrableWall

subroutine addImpenetrableWallPenalty(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use ImpenetrableWall_mod, only : t_ImpenetrableWall

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use CNSHelper
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_ImpenetrableWall) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, nDimensions, nUnknowns, direction, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: localConservedVariables(:), metricsAlongNormalDirection(:),    &
       inviscidPenalty(:), deltaPressure(:), deltaInviscidPenalty(:,:)
  SCALAR_TYPE :: normalMomentum

  assert_key(mode, (FORWARD, ADJOINT))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  if (mode == ADJOINT .and. simulationFlags%useContinuousAdjoint) return

  call startTiming("addImpenetrableWallPenalty")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2)

  allocate(localConservedVariables(nUnknowns))
  allocate(metricsAlongNormalDirection(nDimensions))
  allocate(inviscidPenalty(nUnknowns))
  if (mode == ADJOINT) then
     allocate(deltaPressure(nUnknowns))
     allocate(deltaInviscidPenalty(nUnknowns, nUnknowns))
  end if

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

           localConservedVariables = state%conservedVariables(gridIndex,:)
           metricsAlongNormalDirection =                                                     &
                grid%metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

           normalMomentum = dot_product(localConservedVariables(2:nDimensions+1),            &
                metricsAlongNormalDirection)

           inviscidPenalty(1) = normalMomentum
           inviscidPenalty(2:nDimensions+1) = normalMomentum * state%velocity(gridIndex,:)
           inviscidPenalty(nDimensions+2) =                                                  &
                normalMomentum * state%specificVolume(gridIndex, 1) *                        &
                (localConservedVariables(nDimensions+2) + state%pressure(gridIndex, 1))

           select case (mode)

           case (FORWARD)

              state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -          &
                   this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) * inviscidPenalty

           case (ADJOINT)

              deltaPressure(1) = 0.5_wp * sum(state%velocity(gridIndex,:) ** 2)
              deltaPressure(2:nDimensions+1) = - state%velocity(gridIndex,:)
              deltaPressure(nDimensions+2) = 1.0_wp
              deltaPressure = deltaPressure * (solverOptions%ratioOfSpecificHeats - 1.0_wp)

              select case (nDimensions)
              case (1)
                 call computeJacobianOfInviscidFlux1D(localConservedVariables,               &
                      metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,       &
                      deltaInviscidPenalty, temperature = state%temperature(gridIndex, 1),   &
                      specificVolume = state%specificVolume(gridIndex, 1))
              case (2)
                 call computeJacobianOfInviscidFlux2D(localConservedVariables,               &
                      metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,       &
                      deltaInviscidPenalty, temperature = state%temperature(gridIndex, 1),   &
                      specificVolume = state%specificVolume(gridIndex, 1))
              case (3)
                 call computeJacobianOfInviscidFlux3D(localConservedVariables,               &
                      metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,       &
                      deltaInviscidPenalty, temperature = state%temperature(gridIndex, 1),   &
                      specificVolume = state%specificVolume(gridIndex, 1))
              end select !... nDimensions

              do l = 1, nDimensions
                 deltaInviscidPenalty(l+1,:) = deltaInviscidPenalty(l+1,:) -                 &
                      metricsAlongNormalDirection(l) * deltaPressure
              end do

              state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +          &
                   this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *                &
                   matmul(transpose(deltaInviscidPenalty),                                   &
                   state%adjointVariables(gridIndex,:))

           end select

        end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  SAFE_DEALLOCATE(deltaInviscidPenalty)
  SAFE_DEALLOCATE(deltaPressure)
  SAFE_DEALLOCATE(inviscidPenalty)
  SAFE_DEALLOCATE(metricsAlongNormalDirection)
  SAFE_DEALLOCATE(localConservedVariables)

  call endTiming("addImpenetrableWallPenalty")

end subroutine addImpenetrableWallPenalty

function verifyImpenetrableWallUsage(this, patchDescriptor, gridSize, normalDirection,       &
     extent, simulationFlags, success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags
  use ImpenetrableWall_mod, only : t_ImpenetrableWall

  implicit none

  ! <<< Arguments >>>
  class(t_ImpenetrableWall) :: this
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
  if (normalDirection > size(gridSize) .or. normalDirection == 0) then
     write(message, '(A)') "Normal direction is invalid!"
     return
  end if

  do i = 1, size(gridSize)
     if (extent((i-1)*2+1) < 0 .or. extent((i-1)*2+2) > gridSize(i) .or.                     &
          extent((i-1)*2+1) > extent((i-1)*2+2)) then
        write(message, '(A)') "Invalid extent!"
        return
     end if
  end do

  i = abs(normalDirection)
  if (extent((i-1)*2+1) /= extent((i-1)*2+2)) then
     write(message, '(A)') "Extends more than 1 grid point along normal direction!"
     return
  end if

  success = .true.

  isPatchUsed = .true.

end function verifyImpenetrableWallUsage

subroutine updateImpenetrableWall(this, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use ImpenetrableWall_mod, only : t_ImpenetrableWall

  ! <<< Internal modules >>>
  use CNSHelper, only : computeCartesianViscousFluxes

  implicit none

  ! <<< Arguments >>>
  class(t_ImpenetrableWall) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state

  ! Nothing to be done for this patch type...

end subroutine updateImpenetrableWall
