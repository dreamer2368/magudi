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
  integer :: nDimensions, nUnknowns, nSpecies, direction

  call this%cleanup()
  call this%setupBase(index, comm, patchDescriptor, grid, simulationFlags, solverOptions)

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nSpecies = solverOptions%nSpecies
  assert(nSpecies >= 0)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2 + nSpecies)

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  if (this%nPatchPoints > 0) then
     allocate(this%metrics(this%nPatchPoints, nDimensions))
     if (simulationFlags%viscosityOn) then
        allocate(this%firstPartialViscousJacobians(this%nPatchPoints, nUnknowns, nUnknowns))
        allocate(this%secondPartialViscousJacobians(this%nPatchPoints,                       &
             nUnknowns - 1, nUnknowns - 1, nDimensions))
        allocate(this%viscousPenalty(this%nPatchPoints, nUnknowns))
     end if
  end if

  call this%collect(grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction),       &
       this%metrics)

  write(key, '(A)') "patches/" // trim(patchDescriptor%name) // "/"

  ! Inviscid penalty amount.
  this%inviscidPenaltyAmount = getOption("defaults/inviscid_penalty_amount", 1.0_wp)
  this%inviscidPenaltyAmount = getOption(trim(key) // "inviscid_penalty_amount",             &
       this%inviscidPenaltyAmount)
  this%inviscidPenaltyAmount = sign(this%inviscidPenaltyAmount,                              &
       real(this%normalDirection, wp))
  this%inviscidPenaltyAmount = this%inviscidPenaltyAmount /                                  &
       grid%firstDerivative(abs(this%normalDirection))%normBoundary(1)

  ! Viscous penalty amount.
  if (simulationFlags%viscosityOn) then
     this%viscousPenaltyAmount = getOption("defaults/viscous_penalty_amount", 1.0_wp)
     this%viscousPenaltyAmount = getOption(trim(key) // "viscous_penalty_amount",            &
          this%viscousPenaltyAmount)
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
  SAFE_DEALLOCATE(this%viscousPenalty)
  SAFE_DEALLOCATE(this%firstPartialViscousJacobians)
  SAFE_DEALLOCATE(this%secondPartialViscousJacobians)

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
  integer :: i, j, k, nDimensions, nUnknowns, nSpecies,                                      &
       direction, incomingDirection, gridIndex, patchIndex
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

  nSpecies = solverOptions%nSpecies
  assert(nSpecies >= 0)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2 + nSpecies)

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  if (mode == ADJOINT .and. simulationFlags%useContinuousAdjoint) then
     incomingDirection = -this%normalDirection
  else
     incomingDirection = +this%normalDirection
  end if

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
           metricsAlongNormalDirection = this%metrics(patchIndex,:)

           call computeIncomingJacobianOfInviscidFlux(nDimensions, nSpecies,                 &
                localTargetState, metricsAlongNormalDirection,                               &
                solverOptions%ratioOfSpecificHeats, incomingDirection,                       &
                incomingJacobianOfInviscidFlux)

           select case (mode)

           case (FORWARD)

              state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -          &
                   this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *                &
                   matmul(incomingJacobianOfInviscidFlux,                                    &
                   state%conservedVariables(gridIndex,:) - localTargetState)

              if (simulationFlags%viscosityOn)                                               &
                   state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +     &
                   this%viscousPenaltyAmount * grid%jacobian(gridIndex, 1) *                 &
                   this%viscousPenalty(patchIndex,:)

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

              if (simulationFlags%viscosityOn)                                               &
                   state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -     &
                   this%viscousPenaltyAmount * grid%jacobian(gridIndex, 1) *                 &
                   this%viscousPenalty(patchIndex,:)

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

subroutine computeFarFieldViscousJacobians(this, simulationFlags,                            &
     solverOptions, grid, state)

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
  integer :: i, j, k, l, nDimensions, nUnknowns, nSpecies, direction, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: localTargetState(:), localMetricsAlongFirstDir(:),             &
       localMetricsAlongSecondDir(:), localStressTensor(:), localHeatFlux(:),                &
       localSpeciesFlux(:,:), localVelocity(:), localFirstPartialViscousJacobian(:,:),       &
       localSecondPartialViscousJacobian(:,:)

  if (.not. simulationFlags%viscosityOn) return

  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))
  assert(allocated(state%targetState))

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nSpecies = solverOptions%nSpecies
  assert(nSpecies >= 0)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2 + nSpecies)

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  allocate(localTargetState(nUnknowns))
  allocate(localMetricsAlongFirstDir(nDimensions))
  allocate(localMetricsAlongSecondDir(nDimensions))
  allocate(localStressTensor(nDimensions ** 2))
  allocate(localHeatFlux(nDimensions))
  allocate(localSpeciesFlux(nDimensions, nSpecies))
  allocate(localVelocity(nDimensions))
  allocate(localFirstPartialViscousJacobian(nUnknowns, nUnknowns))
  allocate(localSecondPartialViscousJacobian(nUnknowns - 1, nUnknowns - 1))

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
           localMetricsAlongFirstDir = this%metrics(patchIndex,:)
           localStressTensor = state%stressTensor(gridIndex,:)
           localHeatFlux = state%heatFlux(gridIndex,:)
           localSpeciesFlux = state%speciesFlux(gridIndex,:,:)
           localVelocity = state%velocity(gridIndex,:)

           call computeFirstPartialViscousJacobian(nDimensions, nSpecies,                    &
                localTargetState, localMetricsAlongFirstDir, localStressTensor,              &
                localHeatFlux, localSpeciesFlux, solverOptions%powerLawExponent,             &
                solverOptions%ratioOfSpecificHeats, localFirstPartialViscousJacobian,        &
                specificVolume = state%specificVolume(gridIndex, 1),                         &
                velocity = localVelocity, temperature = state%temperature(gridIndex, 1),     &
                massFraction = state%massFraction(gridIndex, :))

           this%firstPartialViscousJacobians(patchIndex,:,:) =                               &
                localFirstPartialViscousJacobian

           do l = 1, nDimensions

              localMetricsAlongSecondDir =                                                   &
                   grid%metrics(gridIndex,1+nDimensions*(l-1):nDimensions*l)

              call computeSecondPartialViscousJacobian(nDimensions, nSpecies,                &
                   localVelocity, state%dynamicViscosity(gridIndex,1),                       &
                   state%secondCoefficientOfViscosity(gridIndex,1),                          &
                   state%thermalDiffusivity(gridIndex,1),                                    &
                   state%massDiffusivity(gridIndex,:), grid%jacobian(gridIndex,1),           &
                   localMetricsAlongFirstDir, localMetricsAlongSecondDir,                    &
                   localSecondPartialViscousJacobian)

              this%secondPartialViscousJacobians(patchIndex,:,:,l) =                         &
                   localSecondPartialViscousJacobian

           end do !... l = 1, nDimensions

        end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  SAFE_DEALLOCATE(localSecondPartialViscousJacobian)
  SAFE_DEALLOCATE(localFirstPartialViscousJacobian)
  SAFE_DEALLOCATE(localVelocity)
  SAFE_DEALLOCATE(localHeatFlux)
  SAFE_DEALLOCATE(localSpeciesFlux)
  SAFE_DEALLOCATE(localStressTensor)
  SAFE_DEALLOCATE(localMetricsAlongSecondDir)
  SAFE_DEALLOCATE(localMetricsAlongFirstDir)
  SAFE_DEALLOCATE(localTargetState)

end subroutine computeFarFieldViscousJacobians
