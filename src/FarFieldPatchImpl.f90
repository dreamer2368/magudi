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
  integer :: nDimensions, nUnknowns, direction

  call this%cleanup()
  call this%setupBase(index, comm, patchDescriptor, grid, simulationFlags, solverOptions)

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns >= nDimensions + 2)

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  if (this%nPatchPoints > 0) then
     allocate(this%metrics(this%nPatchPoints, nDimensions))
     if (simulationFlags%viscosityOn) then
        allocate(this%viscousFluxes(this%nPatchPoints, nUnknowns - 1))
        allocate(this%targetViscousFluxes(this%nPatchPoints, nUnknowns - 1))
        allocate(this%adjointViscousPenalty(this%nPatchPoints, nUnknowns))
     end if
  end if

  call this%collect(grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction),       &
       this%metrics)

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

  SAFE_DEALLOCATE(this%metrics)
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
  SCALAR_TYPE, allocatable :: localTargetState(:), localConservedVariables(:),               &
       localStressTensor(:), localHeatFlux(:), localSpeciesFlux(:,:), localVelocity(:),      &
       metricsAlongNormalDirection(:), incomingJacobianOfInviscidFlux(:,:),                  &
       partialJacobianOfViscousFlux(:,:)

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
  assert(nUnknowns >= nDimensions + 2)

  if (simulationFlags%viscosityOn .and. mode == FORWARD)                                     &
       call this%collectViscousFluxes(simulationFlags, solverOptions, grid, state)

  if (mode == ADJOINT .and. simulationFlags%useContinuousAdjoint) then
     incomingDirection = -this%normalDirection
  else
     incomingDirection = +this%normalDirection
  end if

  allocate(localTargetState(nUnknowns))
  allocate(metricsAlongNormalDirection(nDimensions))
  allocate(incomingJacobianOfInviscidFlux(nUnknowns, nUnknowns))
  if (mode == ADJOINT .and. simulationFlags%viscosityOn) then
     allocate(partialJacobianOfViscousFlux(nUnknowns, nUnknowns))
     allocate(localConservedVariables(nUnknowns))
     allocate(localStressTensor(nDimensions ** 2))
     allocate(localHeatFlux(nDimensions))
     allocate(localSpeciesFlux(this%nSpecies, nDimensions))
     allocate(localVelocity(nDimensions))
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

           localTargetState = state%targetState(gridIndex,:)
           metricsAlongNormalDirection = this%metrics(patchIndex,:)

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

              if (simulationFlags%viscosityOn) then

                 localConservedVariables = state%conservedVariables(gridIndex,:)
                 localStressTensor = state%stressTensor(gridIndex,:)
                 localHeatFlux = state%heatFlux(gridIndex,:)
                 localSpeciesFlux = state%speciesFlux(gridIndex,:,:)
                 localVelocity = state%velocity(gridIndex,:)

                 select case (nDimensions)
                 case (1)
                    call computeFirstPartialViscousJacobian1D(localConservedVariables,       &
                         metricsAlongNormalDirection, localStressTensor, localHeatFlux,      &
                         solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats, &
                         partialJacobianOfViscousFlux,                                       &
                         specificVolume = state%specificVolume(gridIndex, 1),                &
                         velocity = localVelocity,                                           &
                         temperature = state%temperature(gridIndex, 1))
                 case (2)
                    call computeFirstPartialViscousJacobian2D(localConservedVariables,       &
                         metricsAlongNormalDirection, localStressTensor, localHeatFlux,      &
                         solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats, &
                         partialJacobianOfViscousFlux,                                       &
                         specificVolume = state%specificVolume(gridIndex, 1),                &
                         velocity = localVelocity,                                           &
                         temperature = state%temperature(gridIndex, 1))
                 case (3)
                    call computeFirstPartialViscousJacobian3D(localConservedVariables,       &
                         metricsAlongNormalDirection, localStressTensor, localHeatFlux,      &
                         solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats, &
                         partialJacobianOfViscousFlux,                                       &
                         specificVolume = state%specificVolume(gridIndex, 1),                &
                         velocity = localVelocity,                                           &
                         temperature = state%temperature(gridIndex, 1))
                 end select !... nDimensions

                 this%adjointViscousPenalty(patchIndex,:) =                                  &
                      this%adjointViscousPenalty(patchIndex,:) -                             &
                      grid%jacobian(gridIndex,1) *                                           &
                      matmul(transpose(partialJacobianOfViscousFlux),                        &
                      state%adjointVariables(gridIndex,:))

                 state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +       &
                      this%viscousPenaltyAmount * this%adjointViscousPenalty(patchIndex,:)

              end if

           end select !... mode

        end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  SAFE_DEALLOCATE(localVelocity)
  SAFE_DEALLOCATE(localHeatFlux)
  SAFE_DEALLOCATE(localStressTensor)
  SAFE_DEALLOCATE(localConservedVariables)
  SAFE_DEALLOCATE(partialJacobianOfViscousFlux)
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
       heatFlux(:,:), speciesFlux(:,:,:), viscousFluxes(:,:,:)

  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns >= nDimensions + 2)

  if (this%nPatchPoints > 0 .and. simulationFlags%viscosityOn) then

     allocate(velocity(this%nPatchPoints, nDimensions))
     allocate(stressTensor(this%nPatchPoints, nDimensions ** 2))
     allocate(heatFlux(this%nPatchPoints, nDimensions))
     allocate(speciesFlux(this%nPatchPoints, this%nSpecies, nDimensions))
     allocate(viscousFluxes(this%nPatchPoints, nUnknowns, nDimensions))

     call this%collect(state%velocity, velocity)
     call this%collect(state%stressTensor, stressTensor)
     call this%collect(state%heatFlux, heatFlux)
     call this%collect(state%speciesFlux, speciesFlux)

     call computeCartesianViscousFluxes(nDimensions, velocity,                               &
          stressTensor, heatFlux, speciesFlux, viscousFluxes)

     do i = 1, this%nPatchPoints
        this%viscousFluxes(i,:) = matmul(viscousFluxes(i,2:nUnknowns,:), this%metrics(i,:))
     end do

     SAFE_DEALLOCATE(viscousFluxes)
     SAFE_DEALLOCATE(heatFlux)
     SAFE_DEALLOCATE(speciesFlux)
     SAFE_DEALLOCATE(stressTensor)
     SAFE_DEALLOCATE(velocity)

  end if

end subroutine collectFarFieldViscousFluxes

subroutine updateFarFieldAdjointDiffusionPenaltyTerm(this, simulationFlags, solverOptions,   &
     grid, state, direction, diffusionPenaltyTerm)

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
  integer, intent(in) :: direction
  SCALAR_TYPE, intent(inout) :: diffusionPenaltyTerm(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, gridIndex, patchIndex, nDimensions, nUnknowns
  SCALAR_TYPE, allocatable :: localVelocity(:), localMetricsAlongDirection1(:),              &
       localMetricsAlongDirection2(:), localFluxJacobian(:,:)

  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns >= nDimensions + 2)

  assert(direction >= 1 .and. direction <= nDimensions)

  assert(grid%nGridPoints > 0)
  assert(size(diffusionPenaltyTerm, 1) == grid%nGridPoints)
  assert(size(diffusionPenaltyTerm, 2) == nUnknowns - 1)

  allocate(localVelocity(nDimensions))
  allocate(localMetricsAlongDirection1(nDimensions))
  allocate(localMetricsAlongDirection2(nDimensions))
  allocate(localFluxJacobian(nUnknowns - 1, nUnknowns - 1))

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
           localMetricsAlongDirection1 = this%metrics(patchIndex,:)
           localMetricsAlongDirection2 =                                                     &
                grid%metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

           select case (nDimensions)
           case (1)
              call computeSecondPartialViscousJacobian1D(localVelocity,                      &
                   state%dynamicViscosity(gridIndex,1),                                      &
                   state%secondCoefficientOfViscosity(gridIndex,1),                          &
                   state%thermalDiffusivity(gridIndex,1), grid%jacobian(gridIndex,1),        &
                   localMetricsAlongDirection1(1), localFluxJacobian)
           case (2)
              call computeSecondPartialViscousJacobian2D(localVelocity,                      &
                   state%dynamicViscosity(gridIndex,1),                                      &
                   state%secondCoefficientOfViscosity(gridIndex,1),                          &
                   state%thermalDiffusivity(gridIndex,1), grid%jacobian(gridIndex,1),        &
                   localMetricsAlongDirection1, localMetricsAlongDirection2,                 &
                   localFluxJacobian)
           case (3)
              call computeSecondPartialViscousJacobian3D(localVelocity,                      &
                   state%dynamicViscosity(gridIndex,1),                                      &
                   state%secondCoefficientOfViscosity(gridIndex,1),                          &
                   state%thermalDiffusivity(gridIndex,1), grid%jacobian(gridIndex,1),        &
                   localMetricsAlongDirection1, localMetricsAlongDirection2,                 &
                   localFluxJacobian)
           end select

           diffusionPenaltyTerm(gridIndex,:) = diffusionPenaltyTerm(gridIndex,:) +           &
                matmul(transpose(localFluxJacobian),                                         &
                state%adjointVariables(gridIndex,2:nUnknowns))

        end do
     end do
  end do

  SAFE_DEALLOCATE(localFluxJacobian)
  SAFE_DEALLOCATE(localMetricsAlongDirection2)
  SAFE_DEALLOCATE(localMetricsAlongDirection1)
  SAFE_DEALLOCATE(localVelocity)

end subroutine updateFarFieldAdjointDiffusionPenaltyTerm
