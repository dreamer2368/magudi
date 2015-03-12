#include "config.h"

module PatchImpl

  implicit none

contains

  subroutine allocateData(this, nDimensions, simulationFlags)

    ! <<< Derived types >>>
    use Patch_type
    use PatchDescriptor_type
    use SimulationFlags_type

    ! <<< Arguments >>>
    type(t_Patch) :: this
    integer, intent(in) :: nDimensions
    type(t_SimulationFlags), intent(in) :: simulationFlags

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND

    select case (this%patchType)

    case (SPONGE)
       allocate(this%spongeStrength(this%nPatchPoints), source = 0.0_wp)

    case (SAT_FAR_FIELD)
       allocate(this%metrics(this%nPatchPoints, nDimensions ** 2))
       if (simulationFlags%viscosityOn) then
          allocate(this%viscousFluxes(this%nPatchPoints, nDimensions + 2, nDimensions))
          allocate(this%targetViscousFluxes(this%nPatchPoints, nDimensions + 2, nDimensions))
       end if

    case (SAT_SLIP_WALL)
       allocate(this%metrics(this%nPatchPoints, nDimensions ** 2))

    case (ACTUATOR)
       allocate(this%gradient(this%nPatchPoints, 1))

    end select

  end subroutine allocateData

end module PatchImpl

subroutine parsePatchType(identifier, patchType)

  ! <<< Derived types >>>
  use PatchDescriptor_type

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: identifier
  integer, intent(out) :: patchType

  patchType = -1

  if (trim(identifier) == "SPONGE") then
     patchType = SPONGE
  else if (trim(identifier) == "ACTUATOR") then
     patchType = ACTUATOR
  else if (trim(identifier) == "SOLENOIDAL_EXCITATION_SUPPORT") then
     patchType = SOLENOIDAL_EXCITATION_SUPPORT
  else if (trim(identifier) == "SAT_FAR_FIELD") then
     patchType = SAT_FAR_FIELD
  else if (trim(identifier) == "SAT_SLIP_WALL") then
     patchType = SAT_SLIP_WALL
  else if (trim(identifier) == "SAT_ISOTHERMAL_WALL") then
     patchType = SAT_ISOTHERMAL_WALL
  else if (trim(identifier) == "SAT_BLOCK_INTERFACE") then
     patchType = SAT_BLOCK_INTERFACE
  end if

end subroutine parsePatchType

subroutine setupPatch(this, index, nDimensions, patchDescriptor,                             &
     comm, gridOffset, gridLocalSize, simulationFlags)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_type
  use PatchDescriptor_type
  use SimulationFlags_type

  ! <<< Private members >>>
  use PatchImpl, only : allocateData

  implicit none

  ! <<< Arguments >>>
  type(t_Patch) :: this
  integer, intent(in) :: index, nDimensions
  type(t_PatchDescriptor) :: patchDescriptor
  integer, intent(in) :: comm, gridOffset(3), gridLocalSize(3)
  type(t_SimulationFlags), intent(in) :: simulationFlags

  ! <<< Local variables >>>
  integer :: i

  call cleanupPatch(this)

  this%index = index
  this%normalDirection = patchDescriptor%normalDirection
  this%gridIndex = patchDescriptor%gridIndex
  this%patchType = patchDescriptor%patchType

  this%globalPatchSize(1) = patchDescriptor%iMax - patchDescriptor%iMin + 1
  this%globalPatchSize(2) = patchDescriptor%jMax - patchDescriptor%jMin + 1
  this%globalPatchSize(3) = patchDescriptor%kMax - patchDescriptor%kMin + 1

  this%patchSize(1) = max(patchDescriptor%iMax, gridOffset(1) + 1)
  this%patchSize(2) = max(patchDescriptor%jMax, gridOffset(2) + 1)
  this%patchSize(3) = max(patchDescriptor%kMax, gridOffset(3) + 1)

  this%offset(1) = max(patchDescriptor%iMin, gridOffset(1) + 1)
  this%offset(2) = max(patchDescriptor%jMin, gridOffset(2) + 1)
  this%offset(3) = max(patchDescriptor%kMin, gridOffset(3) + 1)

  do i = 1, 3
     this%offset(i) = min(this%offset(i), gridOffset(i) + gridLocalSize(i)) - 1
     this%patchSize(i) = min(this%patchSize(i), gridOffset(i) + gridLocalSize(i))            &
          - this%offset(i)
  end do

  this%gridLocalSize = gridLocalSize
  this%gridOffset = gridOffset

  this%nPatchPoints = product(this%patchSize)
  this%comm = comm

  this%isCurvilinear = simulationFlags%isDomainCurvilinear

  call allocateData(this, nDimensions, simulationFlags)

end subroutine setupPatch

subroutine cleanupPatch(this)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_type
  use PatchDescriptor_type

  implicit none

  ! <<< Arguments >>>
  type(t_Patch) :: this

  ! <<< Local variables >>>
  integer :: ierror

  SAFE_DEALLOCATE(this%viscousFluxes)
  SAFE_DEALLOCATE(this%targetViscousFluxes)
  SAFE_DEALLOCATE(this%spongeStrength)
  SAFE_DEALLOCATE(this%gradient)

  if (this%comm /= MPI_COMM_NULL) call MPI_Comm_free(this%comm, ierror)

end subroutine cleanupPatch

subroutine collectScalarAtPatch_(this, gridArray, patchArray)

  ! <<< Derived types >>>
  use Patch_type

  ! <<< Arguments >>>
  type(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: gridArray(:)
  SCALAR_TYPE, intent(out) :: patchArray(:)

  ! <<< Local variables >>>
  integer :: i, j, k, patchIndex, localIndex

  do k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
     do j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
        do i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
           patchIndex = i - this%offset(1) +                                                 &
                this%patchSize(1) * (j - 1 - this%offset(2) +                                &
                this%patchSize(2) * (k - 1 - this%offset(3)))
           localIndex = i - this%gridOffset(1) +                                             &
                this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                        &
                this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
           patchArray(patchIndex) = gridArray(localIndex)
        end do
     end do
  end do

end subroutine collectScalarAtPatch_

subroutine collectVectorAtPatch_(this, gridArray, patchArray)

  ! <<< Derived types >>>
  use Patch_type

  ! <<< Arguments >>>
  type(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: gridArray(:,:)
  SCALAR_TYPE, intent(out) :: patchArray(:,:)

  ! <<< Local variables >>>
  integer :: i, j, k, l, patchIndex, localIndex

  do l = 1, size(patchArray, 2)
     do k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
        do j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
           do i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
              patchIndex = i - this%offset(1) +                                              &
                   this%patchSize(1) * (j - 1 - this%offset(2) +                             &
                   this%patchSize(2) * (k - 1 - this%offset(3)))
              localIndex = i - this%gridOffset(1) +                                          &
                   this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                     &
                   this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
              patchArray(patchIndex,l) = gridArray(localIndex,l)
           end do
        end do
     end do
  end do

end subroutine collectVectorAtPatch_

subroutine collectTensorAtPatch_(this, gridArray, patchArray)

  ! <<< Derived types >>>
  use Patch_type

  ! <<< Arguments >>>
  type(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: gridArray(:,:,:)
  SCALAR_TYPE, intent(out) :: patchArray(:,:,:)

  ! <<< Local variables >>>
  integer :: i, j, k, l, m, patchIndex, localIndex

  do m = 1, size(patchArray, 3)
     do l = 1, size(patchArray, 2)
        do k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
           do j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
              do i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
                 patchIndex = i - this%offset(1) +                                           &
                      this%patchSize(1) * (j - 1 - this%offset(2) +                          &
                      this%patchSize(2) * (k - 1 - this%offset(3)))
                 localIndex = i - this%gridOffset(1) +                                       &
                      this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                  &
                      this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
                 patchArray(patchIndex,l,m) = gridArray(localIndex,l,m)
              end do
           end do
        end do
     end do
  end do

end subroutine collectTensorAtPatch_

subroutine addDamping(this, rightHandSide, iblank, spongeAmount,                             &
     actualVariables, targetVariables)

  ! <<< Derived types >>>
  use Patch_type

  implicit none

  ! <<< Arguments >>>
  type(t_Patch) :: this
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)
  integer, intent(in) :: iblank(:)
  real(SCALAR_KIND), intent(in) :: spongeAmount
  SCALAR_TYPE, intent(in) :: actualVariables(:,:)
  SCALAR_TYPE, intent(in), optional :: targetVariables(:,:)

  ! <<< Local variables >>>
  integer :: i, j, k, l, gridIndex, patchIndex

  if (present(targetVariables)) then

     do l = 1, size(rightHandSide, 2)
        do k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
           do j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
              do i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
                 gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                &
                      (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                  &
                      (k - 1 - this%gridOffset(3)))
                 if (iblank(gridIndex) == 0) cycle
                 patchIndex = i - this%offset(1) + this%patchSize(1) *                       &
                      (j - 1 - this%offset(2) + this%patchSize(2) *                          &
                      (k - 1 - this%offset(3)))
                 rightHandSide(gridIndex, l) = rightHandSide(gridIndex, l) +                 &
                      spongeAmount * this%spongeStrength(patchIndex) *                       &
                      (actualVariables(gridIndex, l) - targetVariables(gridIndex, l))
              end do
           end do
        end do
     end do

  else

     do l = 1, size(rightHandSide, 2)
        do k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
           do j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
              do i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
                 gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                &
                      (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                  &
                      (k - 1 - this%gridOffset(3)))
                 if (iblank(gridIndex) == 0) cycle
                 patchIndex = i - this%offset(1) + this%patchSize(1) *                       &
                      (j - 1 - this%offset(2) + this%patchSize(2) *                          &
                      (k - 1 - this%offset(3)))
                 rightHandSide(gridIndex, l) = rightHandSide(gridIndex, l) +                 &
                      spongeAmount * this%spongeStrength(patchIndex) *                       &
                      actualVariables(gridIndex, l)
              end do
           end do
        end do
     end do

  end if

end subroutine addDamping

subroutine addFarFieldPenalty(this, rightHandSide, iblank, inviscidPenaltyAmount,            &
     viscousPenaltyAmount, ratioOfSpecificHeats, nDimensions, &
     conservedVariables, targetState)

  ! <<< Derived types >>>
  use Patch_type

  ! <<< Internal modules >>>
  use CNSHelper

  implicit none

  ! <<< Arguments >>>
  type(t_Patch) :: this
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)
  integer, intent(in) :: iblank(:)
  real(SCALAR_KIND), intent(in) :: inviscidPenaltyAmount,                                    &
       viscousPenaltyAmount, ratioOfSpecificHeats
  integer, intent(in) :: nDimensions
  SCALAR_TYPE, intent(in) :: conservedVariables(:,:), targetState(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, direction, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: localTargetState(:), localMetricsAlongNormalDirection(:),      &
       localIncomingJacobianOfInviscidFluxes(:,:), viscousFluxPenalty(:,:)

  direction = abs(this%normalDirection)

  allocate(localTargetState(nDimensions + 2))
  allocate(localMetricsAlongNormalDirection(nDimensions))
  allocate(localIncomingJacobianOfInviscidFluxes(nDimensions + 2, nDimensions + 2))

  if (allocated(this%viscousFluxes) .and. allocated(this%targetViscousFluxes)) then
     allocate(viscousFluxPenalty(this%nPatchPoints, nDimensions + 1))
     viscousFluxPenalty = 0.0_wp
     do k = 1, nDimensions
        do j = 1, nDimensions + 1
           do i = 1, this%nPatchPoints
              viscousFluxPenalty(i,j) = viscousFluxPenalty(i,j) +                            &
                   this%metrics(i,k+nDimensions*(direction-1)) *                             &
                   (this%viscousFluxes(i,j+1,k) - this%targetViscousFluxes(i,j+1,k))
           end do
        end do
     end do
  end if

  do k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
     do j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
        do i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
           gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                      &
                (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                        &
                (k - 1 - this%gridOffset(3)))
           if (iblank(gridIndex) == 0) cycle
           patchIndex = i - this%offset(1) + this%patchSize(1) *                             &
                (j - 1 - this%offset(2) + this%patchSize(2) *                                &
                (k - 1 - this%offset(3)))

           localTargetState = targetState(gridIndex,:)
           localMetricsAlongNormalDirection =                                                &
                this%metrics(patchIndex,1+nDimensions*(direction-1):nDimensions*direction)

           select case (nDimensions)
           case (1)
              call computeIncomingJacobianOfInviscidFlux1D(localTargetState,                 &
                   localMetricsAlongNormalDirection, ratioOfSpecificHeats,                   &
                   sign(1, this%normalDirection), localIncomingJacobianOfInviscidFluxes)
           case (2)
              call computeIncomingJacobianOfInviscidFlux2D(localTargetState,                 &
                   localMetricsAlongNormalDirection, ratioOfSpecificHeats,                   &
                   sign(1, this%normalDirection), localIncomingJacobianOfInviscidFluxes)
           case (3)
              call computeIncomingJacobianOfInviscidFlux3D(localTargetState,                 &
                   localMetricsAlongNormalDirection, ratioOfSpecificHeats,                   &
                   sign(1, this%normalDirection), localIncomingJacobianOfInviscidFluxes)
           end select

           rightHandSide(gridIndex,:) = rightHandSide(gridIndex,:) -                         &
                sign(inviscidPenaltyAmount, real(this%normalDirection, SCALAR_KIND)) *       &
                matmul(localIncomingJacobianOfInviscidFluxes,                                &
                conservedVariables(gridIndex,:) - localTargetState)
           if (allocated(viscousFluxPenalty)) then
              rightHandSide(gridIndex,2:nDimensions+2) =                                     &
                   rightHandSide(gridIndex,2:nDimensions+2) +                                &
                   sign(viscousPenaltyAmount, real(this%normalDirection, SCALAR_KIND)) *     &
                   viscousFluxPenalty(patchIndex,:)
           end if

        end do
     end do
  end do

  SAFE_DEALLOCATE(viscousFluxPenalty)
  SAFE_DEALLOCATE(localIncomingJacobianOfInviscidFluxes)
  SAFE_DEALLOCATE(localMetricsAlongNormalDirection)
  SAFE_DEALLOCATE(localTargetState)

end subroutine addFarFieldPenalty

subroutine addWallPenalty(this, rightHandSide, iblank, inviscidPenaltyAmount,                &
     viscousPenaltyAmount, ratioOfSpecificHeats, nDimensions, conservedVariables)

  ! <<< Derived types >>>
  use Patch_type

  ! <<< Internal modules >>>
  use CNSHelper

  implicit none

  ! <<< Arguments >>>
  type(t_Patch) :: this
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)
  integer, intent(in) :: iblank(:)
  real(SCALAR_KIND), intent(in) :: inviscidPenaltyAmount,                                    &
       viscousPenaltyAmount, ratioOfSpecificHeats
  integer, intent(in) :: nDimensions
  SCALAR_TYPE, intent(in) :: conservedVariables(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, direction, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: localConservedVariables(:),                                    &
       localMetricsAlongNormalDirection(:), localIncomingJacobianOfInviscidFluxes(:,:),      &
       localInviscidPenalty(:)
  SCALAR_TYPE :: temp

  direction = abs(this%normalDirection)

  allocate(localConservedVariables(nDimensions + 2))
  allocate(localMetricsAlongNormalDirection(nDimensions))
  allocate(localIncomingJacobianOfInviscidFluxes(nDimensions + 2, nDimensions + 2))
  allocate(localInviscidPenalty(nDimensions + 2))

  do k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
     do j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
        do i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
           gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                      &
                (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                        &
                (k - 1 - this%gridOffset(3)))
           if (iblank(gridIndex) == 0) cycle
           patchIndex = i - this%offset(1) + this%patchSize(1) *                             &
                (j - 1 - this%offset(2) + this%patchSize(2) *                                &
                (k - 1 - this%offset(3)))

           localConservedVariables = conservedVariables(gridIndex,:)
           localMetricsAlongNormalDirection =                                                &
                this%metrics(patchIndex,1+nDimensions*(direction-1):nDimensions*direction)

           select case (nDimensions)
           case (1)
              call computeIncomingJacobianOfInviscidFlux1D(localConservedVariables,          &
                   localMetricsAlongNormalDirection, ratioOfSpecificHeats,                   &
                   sign(1, this%normalDirection), localIncomingJacobianOfInviscidFluxes)
           case (2)
              call computeIncomingJacobianOfInviscidFlux2D(localConservedVariables,          &
                   localMetricsAlongNormalDirection, ratioOfSpecificHeats,                   &
                   sign(1, this%normalDirection), localIncomingJacobianOfInviscidFluxes)
           case (3)
              call computeIncomingJacobianOfInviscidFlux3D(localConservedVariables,          &
                   localMetricsAlongNormalDirection, ratioOfSpecificHeats,                   &
                   sign(1, this%normalDirection), localIncomingJacobianOfInviscidFluxes)
           end select

           temp = dot_product(localConservedVariables(2:nDimensions+1),                      &
                localMetricsAlongNormalDirection)

           localInviscidPenalty(1) = 0.0_wp
           localInviscidPenalty(2:nDimensions+1) = localMetricsAlongNormalDirection * temp
           localInviscidPenalty(nDimensions+2) = &
                0.5_wp / localConservedVariables(1) * temp ** 2
           localInviscidPenalty(2:nDimensions+2) = localInviscidPenalty(2:nDimensions+2) /   &
                sum(localMetricsAlongNormalDirection ** 2)
           localInviscidPenalty = matmul(localIncomingJacobianOfInviscidFluxes,              &
                localInviscidPenalty)

           rightHandSide(gridIndex,:) = rightHandSide(gridIndex,:) -                         &
                sign(inviscidPenaltyAmount, real(this%normalDirection, SCALAR_KIND)) *       &
                localInviscidPenalty

        end do
     end do
  end do

  SAFE_DEALLOCATE(localInviscidPenalty)
  SAFE_DEALLOCATE(localIncomingJacobianOfInviscidFluxes)
  SAFE_DEALLOCATE(localMetricsAlongNormalDirection)
  SAFE_DEALLOCATE(localConservedVariables)

end subroutine addWallPenalty
