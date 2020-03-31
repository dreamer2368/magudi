#include "config.h"

subroutine setupBlockInterfacePatch(this, index, comm, patchDescriptor,                      &
     grid, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_BlockInterfacePatch) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: key
  integer :: nDimensions, nUnknowns, nExchangedVariables, procRank, ierror

  call this%cleanup()
  call this%setupBase(index, comm, patchDescriptor, grid, simulationFlags, solverOptions)

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns >= nDimensions + 2)

  if (this%comm /= MPI_COMM_NULL) then

     assert(this%nPatchPoints > 0)

     allocate(this%metricsAlongNormalDirectionL(this%nPatchPoints,nDimensions))
     allocate(this%metricsAlongNormalDirectionR(this%nPatchPoints,nDimensions))

     allocate(this%conservedVariablesL(this%nPatchPoints, nUnknowns))
     allocate(this%conservedVariablesR(this%nPatchPoints, nUnknowns))

     if (simulationFlags%enableAdjoint) then
        allocate(this%adjointVariablesL(this%nPatchPoints, nUnknowns))
        allocate(this%adjointVariablesR(this%nPatchPoints, nUnknowns))
     end if

     if (simulationFlags%viscosityOn) then
        allocate(this%cartesianViscousFluxesL(this%nPatchPoints, nUnknowns, nDimensions))
        allocate(this%viscousFluxesL(this%nPatchPoints, nUnknowns))
        allocate(this%viscousFluxesR(this%nPatchPoints, nUnknowns))
     end if

     !SeungWhan: unclear use of predictionOnly. revisit later.
     nExchangedVariables = 3 * nUnknowns

     call MPI_Comm_rank(this%comm, procRank, ierror)
     if (procRank == 0) then
        allocate(this%sendBuffer(product(this%globalSize), nExchangedVariables))
        allocate(this%receiveBuffer(product(this%globalSize), nExchangedVariables))
     end if

  end if

  write(key, '(A)') "patches/" // trim(patchDescriptor%name) // "/"

  ! Inviscid penalty amount. Note sign is included! the opposite side has different sign, same amplitude.
  this%inviscidPenaltyAmount = getOption("defaults/inviscid_penalty_amount", 1.0_wp)
  this%inviscidPenaltyAmount = getOption(trim(key) //                                        &
       "inviscid_penalty_amount", this%inviscidPenaltyAmount)
  this%inviscidPenaltyAmount = sign(this%inviscidPenaltyAmount,                              &
       real(this%normalDirection, wp))
  this%inviscidPenaltyAmount = this%inviscidPenaltyAmount /                                  &
       grid%firstDerivative(abs(this%normalDirection))%normBoundary(1)

  ! Viscous penalty amount. Note sign is included! the opposite side has different sign, same amplitude.
  if (simulationFlags%viscosityOn) then
     this%viscousPenaltyAmount = getOption("defaults/viscous_penalty_amount", 0.5_wp)
     this%viscousPenaltyAmount = getOption(trim(key) //                                      &
          "viscous_penalty_amount", this%viscousPenaltyAmount)
     this%viscousPenaltyAmount = sign(this%viscousPenaltyAmount,                             &
          real(this%normalDirection, wp))
     this%viscousPenaltyAmount = this%viscousPenaltyAmount /                                 &
          grid%firstDerivative(abs(this%normalDirection))%normBoundary(1)
  else
     this%viscousPenaltyAmount = 0.0_wp
  end if

end subroutine setupBlockInterfacePatch

subroutine cleanupBlockInterfacePatch(this)

  ! <<< Derived types >>>
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  implicit none

  ! <<< Arguments >>>
  class(t_BlockInterfacePatch) :: this

  call this%cleanupBase()

  SAFE_DEALLOCATE(this%conservedVariablesL)
  SAFE_DEALLOCATE(this%conservedVariablesR)
  SAFE_DEALLOCATE(this%adjointVariablesL)
  SAFE_DEALLOCATE(this%adjointVariablesR)
  SAFE_DEALLOCATE(this%cartesianViscousFluxesL)
  SAFE_DEALLOCATE(this%viscousFluxesL)
  SAFE_DEALLOCATE(this%viscousFluxesR)
  SAFE_DEALLOCATE(this%sendBuffer)
  SAFE_DEALLOCATE(this%receiveBuffer)

  SAFE_DEALLOCATE(this%metricsAlongNormalDirectionL)
  SAFE_DEALLOCATE(this%metricsAlongNormalDirectionR)
end subroutine cleanupBlockInterfacePatch

subroutine addBlockInterfacePenalty(this, mode, simulationFlags, solverOptions, grid, state)
! TODO: investigate proper exchange of fluxes.. Is current implementation really correct?
  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT, LINEARIZED

  ! <<< Internal modules >>>
  use CNSHelper
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_BlockInterfacePatch) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, nDimensions, nUnknowns, direction,                                  &
             incomingDirection, gridIndex, patchIndex,                                       &
             incomingDirectionL, incomingDirectionR
  SCALAR_TYPE, allocatable :: localConservedVariablesL(:), localConservedVariablesR(:),      &
       localRoeAverage(:), localStressTensor(:), localHeatFlux(:), localVelocity(:),         &
       localMetricsAlongNormalDirection(:), incomingJacobianOfInviscidFlux(:,:),             &
       localViscousFluxJacobian(:,:), deltaIncomingJacobianOfInviscidFlux(:,:,:),            &
       deltaRoeAverage(:,:), deltaConservedVariablesL(:,:), deltaConservedVariablesR(:,:)

  assert_key(mode, (FORWARD, ADJOINT))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  call startTiming("addBlockInterfacePenalty")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns >= nDimensions + 2)

  if (mode == ADJOINT .and. simulationFlags%useContinuousAdjoint) then
     incomingDirection = -this%normalDirection
  else
     incomingDirection = +this%normalDirection
     incomingDirectionL = +this%normalDirectionL
     incomingDirectionR = +this%normalDirectionR
  end if

  allocate(localConservedVariablesL(nUnknowns))
  allocate(localConservedVariablesR(nUnknowns))
  allocate(localRoeAverage(nUnknowns))
  allocate(localMetricsAlongNormalDirection(nDimensions))
  allocate(incomingJacobianOfInviscidFlux(nUnknowns, nUnknowns))
  allocate(deltaIncomingJacobianOfInviscidFlux(nUnknowns, nUnknowns, nUnknowns))
  allocate(deltaRoeAverage(nUnknowns, nUnknowns))
  allocate(deltaConservedVariablesL(nUnknowns, nUnknowns))
  allocate(deltaConservedVariablesR(nUnknowns, nUnknowns))
  deltaConservedVariablesL = 0.0_wp
  deltaConservedVariablesR = 0.0_wp

  if (simulationFlags%viscosityOn) then
     allocate(localVelocity(nDimensions))
     allocate(localStressTensor(nDimensions ** 2))
     allocate(localHeatFlux(nDimensions))
     allocate(localViscousFluxJacobian(nUnknowns, nUnknowns))
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

           localConservedVariablesL = this%conservedVariablesL(patchIndex,:)
           localConservedVariablesR = this%conservedVariablesR(patchIndex,:)
           localMetricsAlongNormalDirection =                                                &
                grid%metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

           if (simulationFlags%viscosityOn) then
              localVelocity = state%velocity(gridIndex,:)
              localStressTensor = state%stressTensor(gridIndex,:)
              localHeatFlux = state%heatFlux(gridIndex,:)
           end if

           select case(mode)

           case(FORWARD)
             call computeRoeAverage(nDimensions, localConservedVariablesL,                  &
                  localConservedVariablesR, solverOptions%ratioOfSpecificHeats,             &
                  localRoeAverage)

             select case (nDimensions)
             case (1)
               call computeIncomingJacobianOfInviscidFlux1D(localRoeAverage,               &
                    localMetricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,  &
                    incomingDirection, incomingJacobianOfInviscidFlux)
             case (2)
               call computeIncomingJacobianOfInviscidFlux2D(localRoeAverage,               &
                    localMetricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,  &
                    incomingDirection, incomingJacobianOfInviscidFlux)
             case (3)
               call computeIncomingJacobianOfInviscidFlux3D(localRoeAverage,               &
                    localMetricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,  &
                    incomingDirection, incomingJacobianOfInviscidFlux)
             end select !... nDimensions

             state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -          &
                  this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *                &
                  matmul(incomingJacobianOfInviscidFlux,                                    &
                  localConservedVariablesL - localConservedVariablesR)

             if (simulationFlags%viscosityOn) then
                state%rightHandSide(gridIndex,2:nUnknowns) =                                &
                     state%rightHandSide(gridIndex,2:nUnknowns) +                           &
                     this%viscousPenaltyAmount * grid%jacobian(gridIndex, 1) *              &
                     (this%viscousFluxesL(patchIndex,2:nUnknowns) -                         &
                     this%viscousFluxesR(patchIndex,2:nUnknowns))
             end if

           case(ADJOINT)
             if (simulationFlags%useContinuousAdjoint) then !TODO: continuous adjoint is not verified. disabled for now.
               ! call computeRoeAverage(nDimensions, localConservedVariablesL,                  &
               !      localConservedVariablesR, solverOptions%ratioOfSpecificHeats,             &
               !      localRoeAverage)
               !
               ! select case (nDimensions)
               ! case (1)
               !    call computeIncomingJacobianOfInviscidFlux1D(localRoeAverage,               &
               !         localMetricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,  &
               !         incomingDirection, incomingJacobianOfInviscidFlux)
               ! case (2)
               !    call computeIncomingJacobianOfInviscidFlux2D(localRoeAverage,               &
               !         localMetricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,  &
               !         incomingDirection, incomingJacobianOfInviscidFlux)
               ! case (3)
               !    call computeIncomingJacobianOfInviscidFlux3D(localRoeAverage,               &
               !         localMetricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,  &
               !         incomingDirection, incomingJacobianOfInviscidFlux)
               ! end select !... nDimensions
               !
               ! state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -       &
               !      this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *             &
               !      matmul(transpose(incomingJacobianOfInviscidFlux),                      &
               !      this%adjointVariablesL(patchIndex,:) -                                 &
               !      this%adjointVariablesR(patchIndex,:))
               !
               ! if (simulationFlags%viscosityOn) then
               !
               !     select case (nDimensions)
               !     case (1)
               !        call computeFirstPartialViscousJacobian1D(localConservedVariablesL,      &
               !             localMetricsAlongNormalDirection, localStressTensor, localHeatFlux, &
               !             solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats, &
               !             localViscousFluxJacobian, state%specificVolume(gridIndex,1),        &
               !             localVelocity, state%temperature(gridIndex,1))
               !     case (2)
               !        call computeFirstPartialViscousJacobian2D(localConservedVariablesL,      &
               !             localMetricsAlongNormalDirection, localStressTensor, localHeatFlux, &
               !             solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats, &
               !             localViscousFluxJacobian, state%specificVolume(gridIndex,1),        &
               !             localVelocity, state%temperature(gridIndex,1))
               !     case (3)
               !        call computeFirstPartialViscousJacobian3D(localConservedVariablesL,      &
               !             localMetricsAlongNormalDirection, localStressTensor, localHeatFlux, &
               !             solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats, &
               !             localViscousFluxJacobian, state%specificVolume(gridIndex,1),        &
               !             localVelocity, state%temperature(gridIndex,1))
               !     end select !... nDimensions
               !
               !    state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +    &
               !         this%viscousPenaltyAmount * grid%jacobian(gridIndex, 1) *           &
               !         matmul(transpose(localViscousFluxJacobian),                         &
               !         this%adjointVariablesL(patchIndex,:) -                              &
               !         this%adjointVariablesR(patchIndex,:))
               ! end if
             else ! Discrete Adjoint

               call computeRoeAverage(nDimensions, localConservedVariablesL,                  &
                    localConservedVariablesR, solverOptions%ratioOfSpecificHeats,             &
                    localRoeAverage, deltaRoeAverage)

               localMetricsAlongNormalDirection =                                           &
                    this%metricsAlongNormalDirectionL(patchIndex,:)
               select case (nDimensions)
               case (1)
                  call computeIncomingJacobianOfInviscidFlux1D(localRoeAverage,            &
                       localMetricsAlongNormalDirection,                                   &
                       solverOptions%ratioOfSpecificHeats, incomingDirectionL,             &
                       incomingJacobianOfInviscidFlux,                                     &
                       deltaIncomingJacobianOfInviscidFlux, deltaRoeAverage)
               case (2)
                  call computeIncomingJacobianOfInviscidFlux2D(localRoeAverage,            &
                       localMetricsAlongNormalDirection,                                   &
                       solverOptions%ratioOfSpecificHeats, incomingDirectionL,             &
                       incomingJacobianOfInviscidFlux,                                     &
                       deltaIncomingJacobianOfInviscidFlux, deltaRoeAverage)
               case (3)
                  call computeIncomingJacobianOfInviscidFlux3D(localRoeAverage,            &
                       localMetricsAlongNormalDirection,                                   &
                       solverOptions%ratioOfSpecificHeats, incomingDirectionL,             &
                       incomingJacobianOfInviscidFlux,                                     &
                       deltaIncomingJacobianOfInviscidFlux, deltaRoeAverage)
               end select !... nDimensions

               state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +       &
                    this%inviscidPenaltyAmountL * grid%jacobian(gridIndex, 1) *            &
                    matmul(transpose(incomingJacobianOfInviscidFlux),                      &
                    this%adjointVariablesL(patchIndex,:))

               do l = 1, nUnknowns
                  state%rightHandSide(gridIndex,l) = state%rightHandSide(gridIndex,l) +    &
                       this%inviscidPenaltyAmountL * grid%jacobian(gridIndex, 1) *         &
                       sum(matmul(deltaIncomingJacobianOfInviscidFlux(:,:,l),              &
                       localConservedVariablesL - localConservedVariablesR) *              &
                       this%adjointVariablesL(patchIndex,:))
               end do

               if (simulationFlags%viscosityOn) then

                  select case (nDimensions)
                  case (1)
                     call computeFirstPartialViscousJacobian1D(localConservedVariablesL,      &
                          localMetricsAlongNormalDirection, localStressTensor, localHeatFlux, &
                          solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats, &
                          localViscousFluxJacobian, state%specificVolume(gridIndex,1),        &
                          localVelocity, state%temperature(gridIndex,1))
                  case (2)
                     call computeFirstPartialViscousJacobian2D(localConservedVariablesL,      &
                          localMetricsAlongNormalDirection, localStressTensor, localHeatFlux, &
                          solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats, &
                          localViscousFluxJacobian, state%specificVolume(gridIndex,1),        &
                          localVelocity, state%temperature(gridIndex,1))
                  case (3)
                     call computeFirstPartialViscousJacobian3D(localConservedVariablesL,      &
                          localMetricsAlongNormalDirection, localStressTensor, localHeatFlux, &
                          solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats, &
                          localViscousFluxJacobian, state%specificVolume(gridIndex,1),        &
                          localVelocity, state%temperature(gridIndex,1))
                  end select !... nDimensions

                  state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -    &
                       this%viscousPenaltyAmountL * grid%jacobian(gridIndex, 1) *          &
                       matmul(transpose(localViscousFluxJacobian),                         &
                       this%adjointVariablesL(patchIndex,:) )

               end if

               localMetricsAlongNormalDirection =                                           &
                     this%metricsAlongNormalDirectionR(patchIndex,:)
               select case (nDimensions)
               case (1)
                  call computeIncomingJacobianOfInviscidFlux1D(localRoeAverage,            &
                       localMetricsAlongNormalDirection,                                   &
                       solverOptions%ratioOfSpecificHeats, incomingDirectionR,             &
                       incomingJacobianOfInviscidFlux,                                     &
                       deltaIncomingJacobianOfInviscidFlux, deltaRoeAverage)
               case (2)
                  call computeIncomingJacobianOfInviscidFlux2D(localRoeAverage,            &
                       localMetricsAlongNormalDirection,                                   &
                       solverOptions%ratioOfSpecificHeats, incomingDirectionR,             &
                       incomingJacobianOfInviscidFlux,                                     &
                       deltaIncomingJacobianOfInviscidFlux, deltaRoeAverage)
               case (3)
                  call computeIncomingJacobianOfInviscidFlux3D(localRoeAverage,            &
                       localMetricsAlongNormalDirection,                                   &
                       solverOptions%ratioOfSpecificHeats, incomingDirectionR,             &
                       incomingJacobianOfInviscidFlux,                                     &
                       deltaIncomingJacobianOfInviscidFlux, deltaRoeAverage)
               end select !... nDimensions

               ! Note sign change is included in inviscidPenaltyAmount
               state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -       &
                    this%inviscidPenaltyAmountR * grid%jacobian(gridIndex, 1) *            &
                    matmul(transpose(incomingJacobianOfInviscidFlux),                      &
                    this%adjointVariablesR(patchIndex,:))

               ! Note sign change is included in inviscidPenaltyAmount
               do l = 1, nUnknowns
                  state%rightHandSide(gridIndex,l) = state%rightHandSide(gridIndex,l) -    &
                       this%inviscidPenaltyAmountR * grid%jacobian(gridIndex, 1) *         &
                       sum(matmul(deltaIncomingJacobianOfInviscidFlux(:,:,l),              &
                       localConservedVariablesL - localConservedVariablesR) *              &
                       this%adjointVariablesR(patchIndex,:))
               end do

               if (simulationFlags%viscosityOn) then

                  select case (nDimensions)
                  case (1)
                     call computeFirstPartialViscousJacobian1D(localConservedVariablesL,      &
                          localMetricsAlongNormalDirection, localStressTensor, localHeatFlux, &
                          solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats, &
                          localViscousFluxJacobian, state%specificVolume(gridIndex,1),        &
                          localVelocity, state%temperature(gridIndex,1))
                  case (2)
                     call computeFirstPartialViscousJacobian2D(localConservedVariablesL,      &
                          localMetricsAlongNormalDirection, localStressTensor, localHeatFlux, &
                          solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats, &
                          localViscousFluxJacobian, state%specificVolume(gridIndex,1),        &
                          localVelocity, state%temperature(gridIndex,1))
                  case (3)
                     call computeFirstPartialViscousJacobian3D(localConservedVariablesL,      &
                          localMetricsAlongNormalDirection, localStressTensor, localHeatFlux, &
                          solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats, &
                          localViscousFluxJacobian, state%specificVolume(gridIndex,1),        &
                          localVelocity, state%temperature(gridIndex,1))
                  end select !... nDimensions

                  state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +    &
                       this%viscousPenaltyAmountR * grid%jacobian(gridIndex, 1) *          &
                       matmul(transpose(localViscousFluxJacobian),                         &
                       this%adjointVariablesR(patchIndex,:) )

               end if
             end if

           case(LINEARIZED)

             do l = 1, nUnknowns
               deltaConservedVariablesL(l,l) = this%adjointVariablesL(patchIndex,l)
               deltaConservedVariablesR(l,l) = this%adjointVariablesR(patchIndex,l)
             end do

             call computeRoeAverage(nDimensions, localConservedVariablesL,                  &
                  localConservedVariablesR, solverOptions%ratioOfSpecificHeats,             &
                  localRoeAverage, deltaRoeAverage,                                         &
                  deltaConservedVariablesL, deltaConservedVariablesR)

             select case(nDimensions)
             case (1)
               call computeIncomingJacobianOfInviscidFlux1D(localRoeAverage,               &
                    localMetricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,  &
                    incomingDirection, incomingJacobianOfInviscidFlux,                     &
                    deltaIncomingJacobianOfInviscidFlux, deltaRoeAverage)
             case (2)
               call computeIncomingJacobianOfInviscidFlux2D(localRoeAverage,               &
                    localMetricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,  &
                    incomingDirection, incomingJacobianOfInviscidFlux,                     &
                    deltaIncomingJacobianOfInviscidFlux, deltaRoeAverage)
             case (3)
               call computeIncomingJacobianOfInviscidFlux3D(localRoeAverage,               &
                    localMetricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,  &
                    incomingDirection, incomingJacobianOfInviscidFlux,                     &
                    deltaIncomingJacobianOfInviscidFlux, deltaRoeAverage)
             end select !... nDimensions

             state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -          &
                  this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *                &
                  matmul(incomingJacobianOfInviscidFlux,                                    &
                  this%adjointVariablesL(patchIndex,:) - this%adjointVariablesR(patchIndex,:))

             state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -          &
                  this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *                &
                  matmul(sum(deltaIncomingJacobianOfInviscidFlux, dim=3),                   &
                  localConservedVariablesL - localConservedVariablesR)

           end select

        end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  SAFE_DEALLOCATE(localViscousFluxJacobian)
  SAFE_DEALLOCATE(localHeatFlux)
  SAFE_DEALLOCATE(localStressTensor)
  SAFE_DEALLOCATE(localVelocity)
  SAFE_DEALLOCATE(deltaRoeAverage)
  SAFE_DEALLOCATE(deltaIncomingJacobianOfInviscidFlux)
  SAFE_DEALLOCATE(incomingJacobianOfInviscidFlux)
  SAFE_DEALLOCATE(localMetricsAlongNormalDirection)
  SAFE_DEALLOCATE(localRoeAverage)
  SAFE_DEALLOCATE(localConservedVariablesR)
  SAFE_DEALLOCATE(localConservedVariablesL)

  call endTiming("addBlockInterfacePenalty")

end subroutine addBlockInterfacePenalty

function verifyBlockInterfacePatchUsage(this, patchDescriptor, gridSize, normalDirection,    &
     extent, simulationFlags, success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_BlockInterfacePatch) :: this
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

end function verifyBlockInterfacePatchUsage

subroutine collectInterfaceData(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT, LINEARIZED
  use BlockInterfacePatch_enum, only : METRICS

  implicit none

  ! <<< Arguments >>>
  class(t_BlockInterfacePatch) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions, nUnknowns, direction, nExchangedVariables
  SCALAR_TYPE, dimension(:,:), allocatable :: dataToBeSent

  if (this%comm == MPI_COMM_NULL) return

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns >= nDimensions + 2)

  assert(this%nPatchPoints > 0)

  select case(mode)
  case(FORWARD)
    call this%collect(state%conservedVariables, this%conservedVariablesL)

    !SeungWhan: unclear use of predictionOnly. revisit later.
    if (.not. simulationFlags%viscosityOn .and. .not. simulationFlags%enableAdjoint) then
       call this%gatherData(this%conservedVariablesL, this%sendBuffer)
       return
    end if

    if (simulationFlags%viscosityOn) then
       this%viscousFluxesL = 0.0_wp
       do j = 1, nDimensions
          do i = 2, nUnknowns
             this%viscousFluxesL(:,i) = this%viscousFluxesL(:,i) +                             &
                  this%cartesianViscousFluxesL(:,i,j) * this%metricsAlongNormalDirectionL(:,j)
          end do
       end do

       ! For viscous fluxes, metrics on the other side must be reflected on exchanging fluxes.
       this%viscousFluxesR = 0.0_wp
       do j = 1, nDimensions
          do i = 2, nUnknowns
             this%viscousFluxesR(:,i) = this%viscousFluxesR(:,i) +                             &
                  this%cartesianViscousFluxesL(:,i,j) * this%metricsAlongNormalDirectionR(:,j)
          end do
       end do
    end if
  case(ADJOINT)
    call this%collect(state%conservedVariables, this%conservedVariablesL)
    call this%collect(state%adjointVariables, this%adjointVariablesL)
  case(LINEARIZED)
    call this%collect(state%conservedVariables, this%conservedVariablesL)
    call this%collect(state%adjointVariables, this%adjointVariablesL)

    if (simulationFlags%viscosityOn) then
       this%viscousFluxesL = 0.0_wp
       do j = 1, nDimensions
          do i = 2, nUnknowns
             this%viscousFluxesL(:,i) = this%viscousFluxesL(:,i) +                             &
                  this%cartesianViscousFluxesL(:,i,j) * this%metricsAlongNormalDirectionL(:,j)
          end do
       end do

       ! For viscous fluxes, metrics on the other side must be reflected on exchanging fluxes.
       this%viscousFluxesR = 0.0_wp
       do j = 1, nDimensions
          do i = 2, nUnknowns
             this%viscousFluxesR(:,i) = this%viscousFluxesR(:,i) +                             &
                  this%cartesianViscousFluxesL(:,i,j) * this%metricsAlongNormalDirectionR(:,j)
          end do
       end do
    end if
  case(METRICS)
    call this%collect(grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction),    &
                      this%metricsAlongNormalDirectionL)
    this%inviscidPenaltyAmountL = this%inviscidPenaltyAmount
    if (simulationFlags%viscosityOn) this%viscousPenaltyAmountL = this%viscousPenaltyAmount
    this%normalDirectionL = this%normalDirection
  end select

  select case(mode)
  case(FORWARD)
    nExchangedVariables = nUnknowns
    if (simulationFlags%viscosityOn) nExchangedVariables = 2 * nUnknowns
  case(ADJOINT)
    nExchangedVariables = 2 * nUnknowns
  case(LINEARIZED)
    nExchangedVariables = 2 * nUnknowns
    if (simulationFlags%viscosityOn) nExchangedVariables = 3 * nUnknowns
  case(METRICS)
    nExchangedVariables = nUnknowns + 1
  end select

  allocate(dataToBeSent(this%nPatchPoints, nExchangedVariables))
  select case(mode)
  case(FORWARD)
    dataToBeSent(:,1:nUnknowns) = this%conservedVariablesL
    if (simulationFlags%viscosityOn) dataToBeSent(:,nUnknowns+1:) = this%viscousFluxesR
  case(ADJOINT)
    dataToBeSent(:,1:nUnknowns) = this%conservedVariablesL
    dataToBeSent(:,nUnknowns+1:) = this%adjointVariablesL
  case(LINEARIZED)
    dataToBeSent(:,1:nUnknowns) = this%conservedVariablesL
    dataToBeSent(:,nUnknowns+1:2*nUnknowns) = this%adjointVariablesL
    if (simulationFlags%viscosityOn) dataToBeSent(:,2*nUnknowns+1:) = this%viscousFluxesR
  case(METRICS)
    dataToBeSent(:,1:nDimensions) = this%metricsAlongNormalDirectionL
    dataToBeSent(:,nDimensions+1) = this%inviscidPenaltyAmountL
    if (simulationFlags%viscosityOn) dataToBeSent(:,nDimensions+2) = this%viscousPenaltyAmountL
    dataToBeSent(:,nDimensions+3) = real(this%normalDirectionL,wp)
  end select

  call this%gatherData(dataToBeSent, this%sendBuffer)

  SAFE_DEALLOCATE(dataToBeSent)

end subroutine collectInterfaceData

subroutine disperseInterfaceData(this, mode, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT, LINEARIZED
  use BlockInterfacePatch_enum, only : METRICS

  implicit none

  ! <<< Arguments >>>
  class(t_BlockInterfacePatch) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer :: nUnknowns, nExchangedVariables
  SCALAR_TYPE, dimension(:,:), allocatable :: receivedData

  if (this%comm == MPI_COMM_NULL) return

  nUnknowns = solverOptions%nUnknowns

  assert(this%nPatchPoints > 0)

  select case(mode)
  case(FORWARD)
    nExchangedVariables = nUnknowns
    if (simulationFlags%viscosityOn .or. simulationFlags%enableAdjoint)                 &
         nExchangedVariables = nExchangedVariables + nUnknowns
  case(ADJOINT)
    nExchangedVariables = 2 * nUnknowns
  case(LINEARIZED)
    nExchangedVariables = 2 * nUnknowns
    if (simulationFlags%viscosityOn) nExchangedVariables = 3 * nUnknowns
  case(METRICS)
    nExchangedVariables = nUnknowns + 1
  case default
    nExchangedVariables = 2 * nUnknowns
  end select

  allocate(receivedData(this%nPatchPoints, nExchangedVariables))

  call this%scatterData(this%receiveBuffer, receivedData)

  select case(mode)
  case(FORWARD)
    this%conservedVariablesR = receivedData(:,1:nUnknowns)
    if (simulationFlags%viscosityOn) this%viscousFluxesR = receivedData(:,nUnknowns+1:)
  case(ADJOINT)
    this%conservedVariablesR = receivedData(:,1:nUnknowns)
    this%adjointVariablesR = receivedData(:,nUnknowns+1:)
  case(LINEARIZED)
    this%conservedVariablesR = receivedData(:,1:nUnknowns)
    this%adjointVariablesR = receivedData(:,nUnknowns+1:2*nUnknowns)
    if (simulationFlags%viscosityOn) this%viscousFluxesR = receivedData(:,2*nUnknowns+1:)
  case(METRICS)
    this%metricsAlongNormalDirectionR = receivedData(:,1:nUnknowns-2)
    this%inviscidPenaltyAmountR = receivedData(1,nUnknowns-1)
    if (simulationFlags%viscosityOn) this%viscousPenaltyAmountR = receivedData(1,nUnknowns)
    this%normalDirectionR = int(receivedData(1,nUnknowns+1))
  end select

  SAFE_DEALLOCATE(receivedData)

end subroutine disperseInterfaceData

subroutine reshapeReceivedData(this, indexReordering)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_BlockInterfacePatch) :: this
  integer, intent(in) :: indexReordering(3)

  ! <<< Local variables >>>
  integer :: i, j, k, l, order(3), globalSize(3), nComponents

  if (this%comm == MPI_COMM_NULL .or. .not. allocated(this%receiveBuffer)) return

  assert(allocated(this%sendBuffer))
  assert(all(this%globalSize > 0))
  assert(size(this%sendBuffer, 1) == product(this%globalSize))
  assert(size(this%sendBuffer, 2) > 0)
  assert(all(shape(this%receiveBuffer) == shape(this%sendBuffer)))

  assert(indexReordering(3) == 3)

  if (indexReordering(1) == 1 .and. indexReordering(2) == 2) return

  call startTiming("reshapeReceivedData")

  nComponents = size(this%receiveBuffer, 2)
  globalSize = this%globalSize

  order = indexReordering

  if (abs(order(1)) == 2 .and. abs(order(2)) == 1) then

     do l = 1, nComponents
        do k = 1, globalSize(3)
           do j = 1, globalSize(2)
              do i = 1, globalSize(1)
                 this%sendBuffer(i + globalSize(1) * (j - 1 +                                &
                      globalSize(2) * (k - 1)), l) =                                         &
                      this%receiveBuffer(j + globalSize(2) * (i - 1 +                        &
                      globalSize(1) * (k - 1)), l)
              end do
           end do
        end do
     end do

     this%receiveBuffer = this%sendBuffer

     i = order(2)
     order(2) = order(1)
     order(1) = i

  end if

  if (order(1) == -1 .and. order(2) == 2) then

     do l = 1, nComponents
        do k = 1, globalSize(3)
           do j = 1, globalSize(2)
              do i = 1, globalSize(1)
                 this%sendBuffer(i + globalSize(1) * (j - 1 +                                &
                      globalSize(2) * (k - 1)), l) =                                         &
                      this%receiveBuffer(globalSize(1) + 1 - i + globalSize(1) * (j - 1 +    &
                      globalSize(2) * (k - 1)), l)
              end do
           end do
        end do
     end do

     this%receiveBuffer = this%sendBuffer

  else if (order(1) == 1 .and. order(2) == -2) then

     do l = 1, nComponents
        do k = 1, globalSize(3)
           do j = 1, globalSize(2)
              do i = 1, globalSize(1)
                 this%sendBuffer(i + globalSize(1) * (j - 1 +                                &
                      globalSize(2) * (k - 1)), l) =                                         &
                      this%receiveBuffer(i + globalSize(1) * (globalSize(2) - j +            &
                      globalSize(2) * (k - 1)), l)
              end do
           end do
        end do
     end do

     this%receiveBuffer = this%sendBuffer

  else if (order(1) == -1 .and. order(2) == -2) then

     do l = 1, nComponents
        do k = 1, globalSize(3)
           do j = 1, globalSize(2)
              do i = 1, globalSize(1)
                 this%sendBuffer(i + globalSize(1) * (j - 1 +                                &
                      globalSize(2) * (k - 1)), l) =                                         &
                      this%receiveBuffer(globalSize(1) + 1 - i +                             &
                      globalSize(1) * (globalSize(2) - j + globalSize(2) * (k - 1)), l)
              end do
           end do
        end do
     end do

     this%receiveBuffer = this%sendBuffer

  end if

  call endTiming("reshapeReceivedData")

end subroutine reshapeReceivedData
