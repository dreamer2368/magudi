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

  call this%setupBase(index, comm, patchDescriptor, grid, simulationFlags, solverOptions)

  assert_key(this%nDimensions, (1, 2, 3))

  write(key, '(A)') "patches/" // trim(patchDescriptor%name) // "/"

  ! Inviscid penalty amount.
  this%inviscidPenaltyAmount = getOption(trim(key) //                                        &
       "inviscid_penalty_amount", 2.0_wp) !... default value => dual-consistent.
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
  integer :: i, j, k, l, nDimensions, nUnknowns, direction,                                  &
       incomingDirection, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: localConservedVariables(:), unitNormal(:),                     &
       metricsAlongNormalDirection(:), incomingJacobianOfInviscidFlux(:,:),                  &
       inviscidPenalty(:), viscousPenalty(:), deltaNormalMomentum(:),                        &
       deltaInviscidPenalty(:,:), deltaIncomingJacobianOfInviscidFlux(:,:,:)
  SCALAR_TYPE :: normalMomentum

  assert_key(mode, (FORWARD, ADJOINT))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  call startTiming("addImpenetrableWallPenalty")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2)

  if (mode == ADJOINT .and. simulationFlags%useContinuousAdjoint) then
     incomingDirection = -this%normalDirection
  else
     incomingDirection = this%normalDirection
  end if

  allocate(localConservedVariables(nUnknowns))
  allocate(unitNormal(nDimensions))
  allocate(metricsAlongNormalDirection(nDimensions))
  allocate(incomingJacobianOfInviscidFlux(nUnknowns, nUnknowns))
  allocate(inviscidPenalty(nUnknowns))
  if (mode == ADJOINT) then
     allocate(deltaNormalMomentum(nUnknowns))
     allocate(deltaInviscidPenalty(nUnknowns, nUnknowns))
     allocate(deltaIncomingJacobianOfInviscidFlux(nUnknowns, nUnknowns, nUnknowns))
  end if

  do k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
     do j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
        do i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
           gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                      &
                (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                        &
                (k - 1 - this%gridOffset(3)))
           if (grid%iblank(gridIndex) == 0) cycle
           patchIndex = i - this%offset(1) + this%patchSize(1) *                             &
                (j - 1 - this%offset(2) + this%patchSize(2) *                                &
                (k - 1 - this%offset(3)))

           localConservedVariables = state%conservedVariables(gridIndex,:)
           metricsAlongNormalDirection =                                                     &
                grid%metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)
           unitNormal = metricsAlongNormalDirection /                                        &
                sqrt(sum(metricsAlongNormalDirection ** 2))

           normalMomentum =                                                                  &
                dot_product(localConservedVariables(2:nDimensions+1), unitNormal)

           inviscidPenalty = 0.0_wp
           if (mode /= ADJOINT .or. .not. simulationFlags%useContinuousAdjoint) then
              normalMomentum =                                                               &
                   dot_product(localConservedVariables(2:nDimensions+1), unitNormal)
              inviscidPenalty(2:nDimensions+1) = normalMomentum * unitNormal
              inviscidPenalty(nDimensions+2) =                                               &
                   0.5_wp * state%specificVolume(gridIndex, 1) * normalMomentum ** 2
           else
              inviscidPenalty(2:nDimensions+1) =                                             &
                   - (dot_product(state%adjointVariables(gridIndex,2:nDimensions+1),         &
                   unitNormal) - (state%pressure(gridIndex,1) -                              &
                   1.0_wp / solverOptions%ratioOfSpecificHeats)) * unitNormal
           end if

           select case (mode)

           case (FORWARD)

              select case (nDimensions)
              case (1)
                 call computeIncomingJacobianOfInviscidFlux1D(localConservedVariables,       &
                      metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,       &
                      incomingDirection, incomingJacobianOfInviscidFlux,                     &
                      specificVolume = state%specificVolume(gridIndex, 1),                   &
                      temperature = state%temperature(gridIndex, 1))
              case (2)
                 call computeIncomingJacobianOfInviscidFlux2D(localConservedVariables,       &
                      metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,       &
                      incomingDirection, incomingJacobianOfInviscidFlux,                     &
                      specificVolume = state%specificVolume(gridIndex, 1),                   &
                      temperature = state%temperature(gridIndex, 1))
              case (3)
                 call computeIncomingJacobianOfInviscidFlux3D(localConservedVariables,       &
                      metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,       &
                      incomingDirection, incomingJacobianOfInviscidFlux,                     &
                      specificVolume = state%specificVolume(gridIndex, 1),                   &
                      temperature = state%temperature(gridIndex, 1))
              end select !... nDimensions

              state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -          &
                   this%inviscidPenaltyAmount * matmul(incomingJacobianOfInviscidFlux,       &
                   inviscidPenalty)

           case (ADJOINT)

              deltaNormalMomentum(1) = 0.0_wp
              deltaNormalMomentum(2:nDimensions+1) = unitNormal
              deltaNormalMomentum(nDimensions+2) = 0.0_wp

              deltaInviscidPenalty(1,:) = 0.0_wp
              do l = 1, nDimensions
                 deltaInviscidPenalty(l+1,:) = deltaNormalMomentum * unitNormal(l)
              end do
              deltaInviscidPenalty(nDimensions+2,:) =                                        &
                   state%specificVolume(gridIndex, 1) * normalMomentum * deltaNormalMomentum
              deltaInviscidPenalty(nDimensions+2,1) =                                        &
                   deltaInviscidPenalty(nDimensions+2,1) - 0.5_wp * normalMomentum ** 2 *    &
                   state%specificVolume(gridIndex, 1) ** 2

              select case (nDimensions)
              case (1)
                 call computeIncomingJacobianOfInviscidFlux1D(localConservedVariables,       &
                      metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,       &
                      incomingDirection, incomingJacobianOfInviscidFlux,                     &
                      deltaIncomingJacobianOfInviscidFlux,                                   &
                      specificVolume = state%specificVolume(gridIndex, 1),                   &
                      temperature = state%temperature(gridIndex, 1))
              case (2)
                 call computeIncomingJacobianOfInviscidFlux2D(localConservedVariables,       &
                      metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,       &
                      incomingDirection, incomingJacobianOfInviscidFlux,                     &
                      deltaIncomingJacobianOfInviscidFlux,                                   &
                      specificVolume = state%specificVolume(gridIndex, 1),                   &
                      temperature = state%temperature(gridIndex, 1))
              case (3)
                 call computeIncomingJacobianOfInviscidFlux3D(localConservedVariables,       &
                      metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,       &
                      incomingDirection, incomingJacobianOfInviscidFlux,                     &
                      deltaIncomingJacobianOfInviscidFlux,                                   &
                      specificVolume = state%specificVolume(gridIndex, 1),                   &
                      temperature = state%temperature(gridIndex, 1))
              end select !... nDimensions

              if (.not. simulationFlags%useContinuousAdjoint) then
                 do l = 1, nUnknowns
                    state%rightHandSide(gridIndex,l) = state%rightHandSide(gridIndex,l) +    &
                         this%inviscidPenaltyAmount *                                        &
                         dot_product(state%adjointVariables(gridIndex,:),                    &
                         matmul(incomingJacobianOfInviscidFlux, deltaInviscidPenalty(:,l)) + &
                         matmul(deltaIncomingJacobianOfInviscidFlux(:,:,l), inviscidPenalty))
                 end do
              else
                 state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +       &
                      this%inviscidPenaltyAmount *                                           &
                      matmul(transpose(incomingJacobianOfInviscidFlux), inviscidPenalty)
              end if

           end select

        end do !... i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)

  SAFE_DEALLOCATE(viscousPenalty)
  SAFE_DEALLOCATE(deltaIncomingJacobianOfInviscidFlux)
  SAFE_DEALLOCATE(deltaInviscidPenalty)
  SAFE_DEALLOCATE(deltaNormalMomentum)
  SAFE_DEALLOCATE(inviscidPenalty)
  SAFE_DEALLOCATE(incomingJacobianOfInviscidFlux)
  SAFE_DEALLOCATE(metricsAlongNormalDirection)
  SAFE_DEALLOCATE(unitNormal)
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
