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
     this%viscousPenaltyAmounts(1) = getOption("defaults/viscous_penalty_amount", 1.0_wp)
     this%viscousPenaltyAmounts(2) = 0.0_wp
     do i = 1, 2
        write(key, '(A,I0.0)') "patches/" // trim(patchDescriptor%name) //                   &
             "/viscous_penalty_amount", i
        this%viscousPenaltyAmounts(i) = getOption(trim(key), this%viscousPenaltyAmounts(i))
        this%viscousPenaltyAmounts(i) = this%viscousPenaltyAmounts(i) /                      &
             grid%firstDerivative(abs(this%normalDirection))%normBoundary(1)
     end do
     this%viscousPenaltyAmounts(1) = this%viscousPenaltyAmounts(1) *                         &
          solverOptions%reynoldsNumberInverse
     this%viscousPenaltyAmounts(2) = sign(this%viscousPenaltyAmounts(2),                     &
          real(this%normalDirection, wp))
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
  integer :: i, j, k, nDimensions, nSpecies, nUnknowns, direction, gridIndex, patchIndex
  SCALAR_TYPE :: localTemperature
  SCALAR_TYPE, allocatable :: metricsAlongNormalDirection(:), velocity(:,:), temperature(:), &
       metrics(:,:), penaltyAtBoundary(:), penaltyNearBoundary(:,:), temp(:,:),              &
       localVelocity(:), localMassFraction(:), localFluxJacobian(:,:), adjointPenalties(:,:),&
       dynamicViscosity(:), secondCoefficientOfViscosity(:), thermalDiffusivity(:),          &
       massDiffusivity(:,:)

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

  nSpecies = solverOptions%nSpecies
  assert(nSpecies >= 0)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2 + nSpecies)

  if (this%nPatchPoints > 0) then

     allocate(temperature(this%nPatchPoints))
     allocate(dynamicViscosity(this%nPatchPoints))
     allocate(secondCoefficientOfViscosity(this%nPatchPoints))
     allocate(thermalDiffusivity(this%nPatchPoints))
     if (nSpecies > 0) allocate(massDiffusivity(this%nPatchPoints, nSpecies))

     call this%collect(state%temperature(:,1), temperature)

     call computeTransportVariables(solverOptions%nSpecies, temperature,                     &
          solverOptions%powerLawExponent, solverOptions%bulkViscosityRatio,                  &
          solverOptions%ratioOfSpecificHeats, solverOptions%reynoldsNumberInverse,           &
          solverOptions%prandtlNumberInverse, solverOptions%schmidtNumberInverse,            &
          dynamicViscosity, secondCoefficientOfViscosity, thermalDiffusivity, massDiffusivity)

     select case (mode)

     case (FORWARD)

        allocate(penaltyNearBoundary(grid%nGridPoints, nUnknowns - 1))
        allocate(temp(this%nPatchPoints, nUnknowns - 1))
        allocate(velocity(this%nPatchPoints, nDimensions))
        allocate(metrics(this%nPatchPoints, nDimensions))

        call this%collect(state%velocity, velocity)
        call this%collect(grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction), &
             metrics)

        do i = 1, nDimensions
           temp(:,i) = dynamicViscosity * velocity(:,i) * sum(metrics * metrics) +           &
                (dynamicViscosity + secondCoefficientOfViscosity) *                          &
                sum(velocity * metrics, 2) * metrics(:,i)
        end do
        temp(:,nDimensions+1:nUnknowns-1) = 0.0_wp
        if (allocated(this%hole)) then
           do i = 1, size(temp, 2)
              where (this%hole == 1) temp(:, i) = 0.0_wp
           end do
        end if

        call this%disperse(temp, penaltyNearBoundary)

        do i = 1, size(penaltyNearBoundary, 2)
           penaltyNearBoundary(:,i) = grid%jacobian(:,1) * penaltyNearBoundary(:,i)
        end do
        call grid%adjointFirstDerivative(direction)%projectOnBoundaryAndApply(               &
             penaltyNearBoundary, grid%localSize, this%normalDirection)
        call grid%firstDerivative(direction)%applyNorm(penaltyNearBoundary, grid%localSize)
        do i = 1, size(penaltyNearBoundary, 2)
           penaltyNearBoundary(:,i) = grid%jacobian(:,1) * penaltyNearBoundary(:,i)
        end do

        state%rightHandSide(:,2:nUnknowns) = state%rightHandSide(:,2:nUnknowns) -            &
             this%viscousPenaltyAmounts(2) * penaltyNearBoundary

        SAFE_DEALLOCATE(metrics)
        SAFE_DEALLOCATE(velocity)
        SAFE_DEALLOCATE(temp)
        SAFE_DEALLOCATE(penaltyNearBoundary)

     case (ADJOINT)

        allocate(temp(grid%nGridPoints, nUnknowns - 1))

        temp = state%adjointVariables(:,2:nUnknowns)
        call grid%firstDerivative(direction)%applyAndProjectOnBoundary(temp, grid%localSize, &
             this%normalDirection)

     end select

  end if

  allocate(metricsAlongNormalDirection(nDimensions))
  allocate(penaltyAtBoundary(nUnknowns))
  if (mode == ADJOINT) then
     allocate(adjointPenalties(nUnknowns, 2))
     allocate(localVelocity(nDimensions))
     allocate(localFluxJacobian(nUnknowns - 1, nUnknowns - 1))
     if (nSpecies > 0) allocate(localMassFraction(nSpecies))
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
           if (allocated(this%hole)) then
              if (this%hole(patchIndex) == 1) cycle
           end if

           metricsAlongNormalDirection =                                                     &
                grid%metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

           select case (mode)

           case (FORWARD)

              penaltyAtBoundary(1) = 0.0_wp
              penaltyAtBoundary(2:nDimensions+1) =                                           &
                   state%conservedVariables(gridIndex,2:nDimensions+1)
              penaltyAtBoundary(nDimensions+2:nUnknowns) = 0.0_wp
              penaltyAtBoundary = grid%jacobian(gridIndex, 1) * penaltyAtBoundary

              state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -          &
                   this%viscousPenaltyAmounts(1) * penaltyAtBoundary

           case (ADJOINT)

              ! TODO: add adjoint penalties for the species.

              localVelocity = state%velocity(gridIndex,:)
              localTemperature = state%temperature(gridIndex,1)
              if (nSpecies > 0) localMassFraction = state%massFraction(gridIndex,:)

              adjointPenalties(1,1) = 0.0_wp
              adjointPenalties(2:nDimensions+2,1) =                                          &
                   state%adjointVariables(gridIndex,2:nDimensions+2)

              call computeSecondPartialViscousJacobian(nDimensions, nSpecies,                &
                   solverOptions%equationOfState, localVelocity,                             &
                   dynamicViscosity(patchIndex),                                             &
                   secondCoefficientOfViscosity(patchIndex),                                 &
                   thermalDiffusivity(patchIndex), massDiffusivity(patchIndex,:),            &
                   temperature(patchIndex), solverOptions%schmidtNumberInverse,              &
                   solverOptions%molecularWeightInverse, grid%jacobian(gridIndex,1),         &
                   metricsAlongNormalDirection, metricsAlongNormalDirection,                 &
                   localFluxJacobian)

              adjointPenalties(2:nUnknowns,2) = matmul(transpose(localFluxJacobian),         &
                   temp(gridIndex,:))

              adjointPenalties(nDimensions+2,2) = solverOptions%ratioOfSpecificHeats *       &
                   state%specificVolume(gridIndex,1) * adjointPenalties(nDimensions+2,2)
              adjointPenalties(2:nDimensions+1,2) = state%specificVolume(gridIndex,1) *      &
                   adjointPenalties(2:nDimensions+1,2) - localVelocity *                     &
                   adjointPenalties(nDimensions+2,2)
              adjointPenalties(1,2) = - state%specificVolume(gridIndex,1) *                  &
                   state%conservedVariables(gridIndex,nDimensions+2) *                       &
                   adjointPenalties(nDimensions+2,2) - sum(localVelocity *                   &
                   adjointPenalties(2:nDimensions+1,2))

              adjointPenalties = grid%jacobian(gridIndex, 1) * adjointPenalties

              state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +          &
                   this%viscousPenaltyAmounts(1) * adjointPenalties(:,1) -                   &
                   this%viscousPenaltyAmounts(2) * adjointPenalties(:,2)

           end select !... mode

        end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  SAFE_DEALLOCATE(massDiffusivity)
  SAFE_DEALLOCATE(thermalDiffusivity)
  SAFE_DEALLOCATE(dynamicViscosity)
  SAFE_DEALLOCATE(secondCoefficientOfViscosity)
  SAFE_DEALLOCATE(temperature)
  SAFE_DEALLOCATE(localFluxJacobian)
  SAFE_DEALLOCATE(localMassFraction)
  SAFE_DEALLOCATE(localVelocity)
  SAFE_DEALLOCATE(adjointPenalties)
  SAFE_DEALLOCATE(metricsAlongNormalDirection)
  SAFE_DEALLOCATE(temp)

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
