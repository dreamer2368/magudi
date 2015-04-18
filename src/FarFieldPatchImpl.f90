#include "config.h"

subroutine setupFarFieldPatch(this, index, comm, patchDescriptor,                            &
     grid, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use FarFieldPatch_mod, only : t_FarFieldPatch
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_FarFieldPatch) :: this
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
     if (simulationFlags%viscosityOn) then
        allocate(this%viscousFluxes(this%nPatchPoints, solverOptions%nUnknowns - 1))
        allocate(this%targetViscousFluxes(this%nPatchPoints, solverOptions%nUnknowns - 1))
        allocate(this%adjointViscousPenalty(this%nPatchPoints, solverOptions%nUnknowns))
     end if
  end if

  write(key, '(A)') "patches/" // trim(patchDescriptor%name) // "/"

  ! Inviscid penalty amount.
  this%inviscidPenaltyAmount = getOption(trim(key) //                                        &
       "inviscid_penalty_amount", 1.0_wp) !... default value => dual-consistent.
  this%inviscidPenaltyAmount = sign(this%inviscidPenaltyAmount,                              &
       real(this%normalDirection, wp))
  this%inviscidPenaltyAmount = this%inviscidPenaltyAmount /                                  &
       grid%firstDerivative(abs(this%normalDirection))%normBoundary(1)

  ! Viscous penalty amount.
  if (simulationFlags%viscosityOn) then
     this%viscousPenaltyAmount = getOption(trim(key) //                                      &
          "viscous_penalty_amount", 1.0_wp)
     this%viscousPenaltyAmount = this%viscousPenaltyAmount *                                 &
          solverOptions%reynoldsNumberInverse
     this%viscousPenaltyAmount = sign(this%viscousPenaltyAmount,                             &
          real(this%normalDirection, wp))
     this%viscousPenaltyAmount = this%viscousPenaltyAmount /                                 &
          grid%firstDerivative(abs(this%normalDirection))%normBoundary(1)
  else
     this%viscousPenaltyAmount = 0.0_wp
  end if

end subroutine setupFarFieldPatch

subroutine cleanupFarFieldPatch(this)

  ! <<< Derived types >>>
  use FarFieldPatch_mod, only : t_FarFieldPatch

  implicit none

  ! <<< Arguments >>>
  class(t_FarFieldPatch) :: this

  call this%cleanupBase()

  SAFE_DEALLOCATE(this%viscousFluxes)
  SAFE_DEALLOCATE(this%targetViscousFluxes)
  SAFE_DEALLOCATE(this%adjointViscousPenalty)

end subroutine cleanupFarFieldPatch

subroutine addFarFieldPenalty(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use FarFieldPatch_mod, only : t_FarFieldPatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use CNSHelper
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_FarFieldPatch) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, nUnknowns, direction,                                     &
       incomingDirection, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: localTargetState(:), metricsAlongNormalDirection(:),           &
       incomingJacobianOfInviscidFlux(:,:)

  assert_key(mode, (FORWARD, ADJOINT))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))
  assert(allocated(state%targetState))

  call startTiming("addFarFieldPenalty")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2)

  select case (mode)

  case (FORWARD)

     if (simulationFlags%viscosityOn)                                                        &
          call this%collectViscousFluxes(simulationFlags, solverOptions, grid, state)

  case (ADJOINT)

     if (simulationFlags%viscosityOn)                                                        &
          call this%computeAdjointViscousPenalty(simulationFlags, solverOptions, grid, state)

     if (simulationFlags%useContinuousAdjoint) then
        incomingDirection = -this%normalDirection
     else
        incomingDirection = +this%normalDirection
     end if

  end select

  allocate(localTargetState(nUnknowns))
  allocate(metricsAlongNormalDirection(nDimensions))
  allocate(incomingJacobianOfInviscidFlux(nUnknowns, nUnknowns))

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

           localTargetState = state%targetState(gridIndex,:)
           metricsAlongNormalDirection =                                                     &
                grid%metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

           select case (nDimensions)
           case (1)
              call computeIncomingJacobianOfInviscidFlux1D(localTargetState,                 &
                   metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,          &
                   incomingDirection, incomingJacobianOfInviscidFlux)
           case (2)
              call computeIncomingJacobianOfInviscidFlux2D(localTargetState,                 &
                   metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,          &
                   incomingDirection, incomingJacobianOfInviscidFlux)
           case (3)
              call computeIncomingJacobianOfInviscidFlux3D(localTargetState,                 &
                   metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,          &
                   incomingDirection, incomingJacobianOfInviscidFlux)
           end select !... nDimensions

           select case (mode)

           case (FORWARD)

              state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -          &
                   this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *                &
                   matmul(incomingJacobianOfInviscidFlux,                                    &
                   state%conservedVariables(gridIndex,:) - localTargetState)

              if (simulationFlags%viscosityOn) then
                 state%rightHandSide(gridIndex,2:nUnknowns) =                                &
                      state%rightHandSide(gridIndex,2:nUnknowns) +                           &
                      this%viscousPenaltyAmount * grid%jacobian(gridIndex, 1) *              &
                      (this%viscousFluxes(patchIndex,:) -                                    &
                      this%targetViscousFluxes(patchIndex,:))
              end if

           case (ADJOINT)

              if (simulationFlags%useContinuousAdjoint) then
                 state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -       &
                      this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *             &
                      matmul(transpose(incomingJacobianOfInviscidFlux),                      &
                      state%adjointVariables(gridIndex,:))
              else
                 state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +       &
                      this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *             &
                      matmul(transpose(incomingJacobianOfInviscidFlux),                      &
                      state%adjointVariables(gridIndex,:))
              end if

              state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -          &
                   this%viscousPenaltyAmount * this%adjointViscousPenalty(patchIndex,:)

           end select !... mode

        end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  SAFE_DEALLOCATE(incomingJacobianOfInviscidFlux)
  SAFE_DEALLOCATE(metricsAlongNormalDirection)
  SAFE_DEALLOCATE(localTargetState)

  call endTiming("addFarFieldPenalty")

end subroutine addFarFieldPenalty

function verifyFarFieldPatchUsage(this, patchDescriptor, gridSize, normalDirection,          &
     extent, simulationFlags, success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use FarFieldPatch_mod, only : t_FarFieldPatch
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_FarFieldPatch) :: this
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

  if (.not. simulationFlags%useTargetState) then
     write(message, '(A)')                                                                   &
          "No target state available for enforcing far-field boundary conditions!"
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

end function verifyFarFieldPatchUsage

subroutine collectFarFieldViscousFluxes(this, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use FarFieldPatch_mod, only : t_FarFieldPatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use CNSHelper, only : computeCartesianViscousFluxes

  implicit none

  ! <<< Arguments >>>
  class(t_FarFieldPatch) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions, nUnknowns, direction
  SCALAR_TYPE, allocatable :: velocity(:,:), stressTensor(:,:),                              &
       heatFlux(:,:), viscousFluxes(:,:,:)

  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2)

  if (this%nPatchPoints > 0 .and. simulationFlags%viscosityOn) then

     allocate(velocity(this%nPatchPoints, nDimensions))
     allocate(stressTensor(this%nPatchPoints, nDimensions ** 2))
     allocate(heatFlux(this%nPatchPoints, nDimensions))
     allocate(viscousFluxes(this%nPatchPoints, nUnknowns, nDimensions))

     call this%collect(state%velocity, velocity)
     call this%collect(state%stressTensor, stressTensor)
     call this%collect(state%heatFlux, heatFlux)

     call computeCartesianViscousFluxes(nDimensions, velocity,                               &
          stressTensor, heatFlux, viscousFluxes)

     do i = 1, this%nPatchPoints
        this%viscousFluxes(i,:) = matmul(viscousFluxes(i,2:nUnknowns,:),                     &
             grid%metrics(i,1+nDimensions*(direction-1):nDimensions*direction))
     end do

     SAFE_DEALLOCATE(viscousFluxes)
     SAFE_DEALLOCATE(heatFlux)
     SAFE_DEALLOCATE(stressTensor)
     SAFE_DEALLOCATE(velocity)

  end if

end subroutine collectFarFieldViscousFluxes

subroutine computeFarFieldAdjointViscousPenalty(this,                                        &
     simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use FarFieldPatch_mod, only : t_FarFieldPatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use CNSHelper

  implicit none

  ! <<< Arguments >>>
  class(t_FarFieldPatch) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, gridIndex, patchIndex, nDimensions, nUnknowns, direction
  SCALAR_TYPE, allocatable :: temp(:,:,:), localVelocity(:), localMetricsAlongDirection1(:), &
       localMetricsAlongDirection2(:), localFluxJacobian1(:,:), localFluxJacobian2(:,:),     &
       localConservedVariables(:), localStressTensor(:), localHeatFlux(:)

  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2)

  allocate(temp(grid%nGridPoints, nUnknowns - 1, nDimensions))
  allocate(localVelocity(nDimensions))
  allocate(localMetricsAlongDirection1(nDimensions))
  allocate(localMetricsAlongDirection2(nDimensions))
  allocate(localFluxJacobian1(nUnknowns, nUnknowns))
  allocate(localFluxJacobian2(nUnknowns - 1, nUnknowns - 1))
  allocate(localConservedVariables(nUnknowns))
  allocate(localStressTensor(nDimensions ** 2))
  allocate(localHeatFlux(nDimensions))

  do i = 1, grid%nGridPoints

     localVelocity = state%velocity(i,:)
     localMetricsAlongDirection1 =                                                           &
          grid%metrics(i,1+nDimensions*(direction-1):nDimensions*direction)

     do j = 1, nDimensions

        localMetricsAlongDirection2 = grid%metrics(i,1+nDimensions*(j-1):nDimensions*j)

        select case (nDimensions)
        case (1)
           call computeSecondPartialViscousJacobian1D(localVelocity,                         &
                state%dynamicViscosity(i,1), state%secondCoefficientOfViscosity(i,1),        &
                state%thermalDiffusivity(i,1), grid%jacobian(i,1),                           &
                localMetricsAlongDirection1(1), localFluxJacobian2)
        case (2)
           call computeSecondPartialViscousJacobian2D(localVelocity,                         &
                state%dynamicViscosity(i,1), state%secondCoefficientOfViscosity(i,1),        &
                state%thermalDiffusivity(i,1), grid%jacobian(i,1),                           &
                localMetricsAlongDirection1, localMetricsAlongDirection2, localFluxJacobian2)
        case (3)
           call computeSecondPartialViscousJacobian3D(localVelocity,                         &
                state%dynamicViscosity(i,1), state%secondCoefficientOfViscosity(i,1),        &
                state%thermalDiffusivity(i,1), grid%jacobian(i,1),                           &
                localMetricsAlongDirection1, localMetricsAlongDirection2, localFluxJacobian2)
        end select

        temp(i,:,j) = matmul(transpose(localFluxJacobian2),                                  &
             state%adjointVariables(i,2:nUnknowns))

     end do !... j = 1, nDimensions

  end do !... i = 1, grid%nGridPoints

  do i = 1, nDimensions
     call grid%adjointFirstDerivative(i)%applyAtDomainBoundary(temp(:,:,i),                  &
          grid%localSize, this%normalDirection)
  end do

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

           localVelocity = state%velocity(gridIndex,:)
           localMetricsAlongDirection1 =                                                     &
                grid%metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)
           localConservedVariables = state%conservedVariables(gridIndex,:)
           localStressTensor = state%stressTensor(gridIndex,:)
           localHeatFlux = state%heatFlux(gridIndex,:)

           this%adjointViscousPenalty(patchIndex,2:nUnknowns) =                              &
                sum(temp(gridIndex,:,:), dim = 2)

           this%adjointViscousPenalty(patchIndex,nDimensions+2) =                            &
                solverOptions%ratioOfSpecificHeats * state%specificVolume(gridIndex,1) *     &
                this%adjointViscousPenalty(patchIndex,nDimensions+2)
           do l = 1, nDimensions
              this%adjointViscousPenalty(patchIndex,l+1) =                                   &
                   this%adjointViscousPenalty(patchIndex,l+1) *                              &
                   state%specificVolume(gridIndex,1) - localVelocity(l) *                    &
                   this%adjointViscousPenalty(patchIndex,nDimensions+2)
           end do

           this%adjointViscousPenalty(patchIndex,1) = - state%specificVolume(gridIndex,1) *  &
                state%conservedVariables(gridIndex,nDimensions+2) *                          &
                this%adjointViscousPenalty(patchIndex,nDimensions+2) -                       &
                sum(this%adjointViscousPenalty(patchIndex,2:nDimensions+1) * localVelocity)

           select case (nDimensions)
           case (1)
              call computeFirstPartialViscousJacobian1D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian1, state%specificVolume(gridIndex,1), localVelocity,     &
                   state%temperature(gridIndex,1))
           case (2)
              call computeFirstPartialViscousJacobian2D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian1, state%specificVolume(gridIndex,1), localVelocity,     &
                   state%temperature(gridIndex,1))
           case (3)
              call computeFirstPartialViscousJacobian3D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian1, state%specificVolume(gridIndex,1), localVelocity,     &
                   state%temperature(gridIndex,1))
           end select

           this%adjointViscousPenalty(patchIndex,:) =                                        &
                this%adjointViscousPenalty(patchIndex,:) +                                   &
                matmul(transpose(localFluxJacobian1), state%adjointVariables(gridIndex,:))
           this%adjointViscousPenalty(patchIndex,:) = grid%jacobian(gridIndex,1) *           &
                this%adjointViscousPenalty(patchIndex,:)

        end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  SAFE_DEALLOCATE(localHeatFlux)
  SAFE_DEALLOCATE(localStressTensor)
  SAFE_DEALLOCATE(localConservedVariables)
  SAFE_DEALLOCATE(localFluxJacobian2)
  SAFE_DEALLOCATE(localFluxJacobian1)
  SAFE_DEALLOCATE(localMetricsAlongDirection2)
  SAFE_DEALLOCATE(localMetricsAlongDirection1)
  SAFE_DEALLOCATE(localVelocity)
  SAFE_DEALLOCATE(temp)

end subroutine computeFarFieldAdjointViscousPenalty
