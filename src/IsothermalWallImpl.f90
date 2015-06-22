#include "config.h"

subroutine setupIsothermalWall(this, index, comm, patchDescriptor,                           &
     grid, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use IsothermalWall_mod, only : t_IsothermalWall
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use CNSHelper, only : computeTransportVariables
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_IsothermalWall) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: key
  SCALAR_TYPE :: wallTemperature
  integer :: i

  call this%cleanup()
  call this%t_ImpenetrableWall%setup(index, comm, patchDescriptor,                           &
       grid, simulationFlags, solverOptions)

  write(key, '(A,I0.0)') "patches/" // trim(patchDescriptor%name) // "/"

  if (this%nPatchPoints > 0 .and. simulationFlags%viscosityOn) then

     allocate(this%temperature(this%nPatchPoints))
     allocate(this%dynamicViscosity(this%nPatchPoints))
     allocate(this%secondCoefficientOfViscosity(this%nPatchPoints))
     allocate(this%thermalDiffusivity(this%nPatchPoints))

     ! Wall temperature (used only if target state is not present).
     if (.not. simulationFlags%useTargetState) then
        wallTemperature = getOption(trim(key) // "temperature",                              &
             1.0_wp / (solverOptions%ratioOfSpecificHeats - 1.0_wp))
        this%temperature(:) = wallTemperature
        call computeTransportVariables(this%temperature, solverOptions%powerLawExponent,     &
             solverOptions%bulkViscosityRatio, solverOptions%ratioOfSpecificHeats,           &
             solverOptions%reynoldsNumberInverse, solverOptions%prandtlNumberInverse,        &
             this%dynamicViscosity, this%secondCoefficientOfViscosity,                       &
             this%thermalDiffusivity)
     end if

  end if

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

end subroutine setupIsothermalWall

subroutine cleanupIsothermalWall(this)

  ! <<< Derived types >>>
  use IsothermalWall_mod, only : t_IsothermalWall

  implicit none

  ! <<< Arguments >>>
  class(t_IsothermalWall) :: this

  call this%t_ImpenetrableWall%cleanup()

  SAFE_DEALLOCATE(this%temperature)
  SAFE_DEALLOCATE(this%dynamicViscosity)
  SAFE_DEALLOCATE(this%secondCoefficientOfViscosity)
  SAFE_DEALLOCATE(this%thermalDiffusivity)

end subroutine cleanupIsothermalWall

subroutine addIsothermalWallPenalty(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use IsothermalWall_mod, only : t_IsothermalWall
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use CNSHelper
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_IsothermalWall) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, nUnknowns, direction, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: unitNormal(:), metricsAlongNormalDirection(:), velocity(:,:),  &
       temperature(:), metrics(:,:), viscousPenalties(:,:), extraPenaltyPatch(:,:),          &
       extraPenaltyGrid(:,:), localVelocity(:), localFluxJacobian(:,:)

  assert_key(mode, (FORWARD, ADJOINT))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  if (mode == ADJOINT .and. simulationFlags%useContinuousAdjoint) return

  call startTiming("addIsothermalWallPenalty")

  call this%t_ImpenetrableWall%updateRhs(mode, simulationFlags, solverOptions, grid, state)

  if (.not. simulationFlags%viscosityOn) then
     call endTiming("addIsothermalWallPenalty")
     return
  end if

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns >= nDimensions + 2)

  if (this%nPatchPoints > 0) then

     allocate(extraPenaltyGrid(grid%nGridPoints, nUnknowns - 1))
     allocate(extraPenaltyPatch(this%nPatchPoints, nUnknowns - 1))

     select case (mode)

     case (FORWARD)

        allocate(velocity(this%nPatchPoints, nDimensions))
        allocate(temperature(this%nPatchPoints))
        allocate(metrics(this%nPatchPoints, nDimensions))

        call this%collect(state%velocity, velocity)
        call this%collect(state%temperature(:,1), temperature)
        call this%collect(grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction), &
             metrics)

        do i = 1, nDimensions
           extraPenaltyPatch(:,i) = this%dynamicViscosity * velocity(:,i) *                  &
                sum(metrics * metrics) + (this%dynamicViscosity +                            &
                this%secondCoefficientOfViscosity) * sum(velocity * metrics, 2) *            &
                metrics(:,i)
        end do
        extraPenaltyPatch(:,nDimensions+1) = this%thermalDiffusivity *                       &
             (temperature - this%temperature)

        call this%disperse(extraPenaltyPatch, extraPenaltyGrid)

        do i = 1, size(extraPenaltyGrid, 2)
           extraPenaltyGrid(:,i) = grid%jacobian(:,1) * extraPenaltyGrid(:,i)
        end do
        call grid%adjointFirstDerivative(direction)%projectOnBoundaryAndApply(               &
             extraPenaltyGrid, grid%localSize, this%normalDirection)
        call grid%firstDerivative(direction)%applyNorm(extraPenaltyGrid, grid%localSize)

        SAFE_DEALLOCATE(metrics)
        SAFE_DEALLOCATE(temperature)
        SAFE_DEALLOCATE(velocity)

     case (ADJOINT)

        extraPenaltyGrid = state%adjointVariables(:,2:nUnknowns)
        call grid%firstDerivative(direction)%applyAndProjectOnBoundary(extraPenaltyGrid,     &
             grid%localSize, this%normalDirection)

     end select

     call this%collect(extraPenaltyGrid, extraPenaltyPatch)
     SAFE_DEALLOCATE(extraPenaltyGrid)

  end if

  allocate(unitNormal(nDimensions))
  allocate(metricsAlongNormalDirection(nDimensions))
  allocate(viscousPenalties(nUnknowns, 2))
  if (mode == ADJOINT) then
     allocate(localVelocity(nDimensions))
     allocate(localFluxJacobian(nUnknowns - 1, nUnknowns - 1))
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

           metricsAlongNormalDirection =                                                     &
                grid%metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)
           unitNormal = metricsAlongNormalDirection /                                        &
                sqrt(sum(metricsAlongNormalDirection ** 2))

           select case (mode)

           case (FORWARD)

              viscousPenalties(1,:) = 0.0_wp
              viscousPenalties(2:nDimensions+2,1) =                                          &
                   state%conservedVariables(gridIndex,2:nDimensions+2)
              viscousPenalties(nDimensions+2,1) = viscousPenalties(nDimensions+2,1) -        &
                   state%conservedVariables(gridIndex,1) *                                   &
                   this%temperature(patchIndex) / solverOptions%ratioOfSpecificHeats
              viscousPenalties(2:nUnknowns,2) = extraPenaltyPatch(patchIndex,:)

              viscousPenalties = grid%jacobian(gridIndex, 1) * viscousPenalties

              state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -          &
                   this%viscousPenaltyAmounts(1) * viscousPenalties(:,1) -                   &
                   this%viscousPenaltyAmounts(2) * viscousPenalties(:,2)

           case (ADJOINT)

              localVelocity = state%velocity(gridIndex,:)

              viscousPenalties(1,1) = - state%adjointVariables(gridIndex,nDimensions+2) *    &
                   this%temperature(patchIndex) / solverOptions%ratioOfSpecificHeats
              viscousPenalties(2:nDimensions+2,1) =                                          &
                   state%adjointVariables(gridIndex,2:nDimensions+2)

              select case (nDimensions)
              case (1)
                 call computeSecondPartialViscousJacobian1D(localVelocity,                   &
                      this%dynamicViscosity(patchIndex),                                     &
                      this%secondCoefficientOfViscosity(patchIndex),                         &
                      this%thermalDiffusivity(patchIndex), grid%jacobian(gridIndex,1),       &
                      metricsAlongNormalDirection(1), localFluxJacobian)
              case (2)
                 call computeSecondPartialViscousJacobian2D(localVelocity,                   &
                      this%dynamicViscosity(patchIndex),                                     &
                      this%secondCoefficientOfViscosity(patchIndex),                         &
                      this%thermalDiffusivity(patchIndex), grid%jacobian(gridIndex,1),       &
                      metricsAlongNormalDirection, metricsAlongNormalDirection,              &
                      localFluxJacobian)
              case (3)
                 call computeSecondPartialViscousJacobian3D(localVelocity,                   &
                      this%dynamicViscosity(patchIndex),                                     &
                      this%secondCoefficientOfViscosity(patchIndex),                         &
                      this%thermalDiffusivity(patchIndex), grid%jacobian(gridIndex,1),       &
                      metricsAlongNormalDirection, metricsAlongNormalDirection,              &
                      localFluxJacobian)
              end select !... nDimensions

              viscousPenalties(2:nUnknowns,2) = matmul(transpose(localFluxJacobian),         &
                   extraPenaltyPatch(patchIndex,:))

              viscousPenalties(nDimensions+2,2) = solverOptions%ratioOfSpecificHeats *       &
                   state%specificVolume(gridIndex,1) * viscousPenalties(nDimensions+2,2)
              viscousPenalties(2:nDimensions+1,2) = state%specificVolume(gridIndex,1) *      &
                   viscousPenalties(2:nDimensions+1,2) - localVelocity *                     &
                   viscousPenalties(nDimensions+2,2)
              viscousPenalties(1,2) = - state%specificVolume(gridIndex,1) *                  &
                   state%conservedVariables(gridIndex,nDimensions+2) *                       &
                   viscousPenalties(nDimensions+2,2) - sum(localVelocity *                   &
                   viscousPenalties(2:nDimensions+1,2))

              viscousPenalties = grid%jacobian(gridIndex, 1) * viscousPenalties

              state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +          &
                   this%viscousPenaltyAmounts(1) * viscousPenalties(:,1) -                   &
                   this%viscousPenaltyAmounts(2) * viscousPenalties(:,2)

           end select !... mode

        end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  SAFE_DEALLOCATE(localFluxJacobian)
  SAFE_DEALLOCATE(localVelocity)
  SAFE_DEALLOCATE(viscousPenalties)
  SAFE_DEALLOCATE(metricsAlongNormalDirection)
  SAFE_DEALLOCATE(unitNormal)
  SAFE_DEALLOCATE(extraPenaltyPatch)

  call endTiming("addIsothermalWallPenalty")

end subroutine addIsothermalWallPenalty

function verifyIsothermalWallUsage(this, patchDescriptor, gridSize, normalDirection,         &
     extent, simulationFlags, success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use IsothermalWall_mod, only : t_IsothermalWall
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_IsothermalWall) :: this
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  integer, intent(in) :: gridSize(:), normalDirection, extent(6)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  logical, intent(out) :: success
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Result >>>
  logical :: isPatchUsed

  isPatchUsed = this%t_ImpenetrableWall%verifyUsage(patchDescriptor, gridSize,               &
       normalDirection, extent, simulationFlags, success, message)

end function verifyIsothermalWallUsage
