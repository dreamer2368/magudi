#include "config.h"

module PatchImpl

  implicit none
  public

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
    case (SAT_FAR_FIELD, SAT_SLIP_WALL, SAT_ISOTHERMAL_WALL, SAT_ADIABATIC_WALL)
       allocate(this%metrics(this%nPatchPoints, nDimensions ** 2))
    end select

    select case (this%patchType)

    case (SAT_FAR_FIELD)
       if (simulationFlags%viscosityOn) then
          allocate(this%viscousFluxes(this%nPatchPoints, nDimensions + 2, nDimensions))
          allocate(this%targetViscousFluxes(this%nPatchPoints, nDimensions + 2, nDimensions))
       end if

    case (SPONGE)
       allocate(this%spongeStrength(this%nPatchPoints), source = 0.0_wp)

    case (ACTUATOR)
       allocate(this%gradient(this%nPatchPoints, 1))

    end select

  end subroutine allocateData

end module PatchImpl

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

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  type(t_Patch) :: this
  integer, intent(in) :: index, nDimensions
  type(t_PatchDescriptor) :: patchDescriptor
  integer, intent(in) :: comm, gridOffset(3), gridLocalSize(3)
  type(t_SimulationFlags), intent(in) :: simulationFlags

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i
  character(len = STRING_LENGTH) :: key

  call cleanupPatch(this)

  this%index = index
  this%normalDirection = patchDescriptor%normalDirection
  this%gridIndex = patchDescriptor%gridIndex
  this%patchType = patchDescriptor%patchType

  this%extent = (/ patchDescriptor%iMin, patchDescriptor%iMax,                               &
       patchDescriptor%jMin, patchDescriptor%jMax,                                           &
       patchDescriptor%kMin, patchDescriptor%kMax /)

  ! Global patch size.
  this%globalPatchSize(1) = patchDescriptor%iMax - patchDescriptor%iMin + 1
  this%globalPatchSize(2) = patchDescriptor%jMax - patchDescriptor%jMin + 1
  this%globalPatchSize(3) = patchDescriptor%kMax - patchDescriptor%kMin + 1

  ! Zero-based index of first point on the patch belonging to the ``current'' process (this
  ! value has no meaning if the patch lies outside the part of the grid belonging to the
  ! ``current'' process).
  this%offset(1) = max(patchDescriptor%iMin, gridOffset(1) + 1)
  this%offset(2) = max(patchDescriptor%jMin, gridOffset(2) + 1)
  this%offset(3) = max(patchDescriptor%kMin, gridOffset(3) + 1)
  this%offset = min(this%offset, gridOffset + gridLocalSize) - 1

  ! Extent of the patch belonging to the ``current'' process (this value has no meaning if the
  ! patch lies outside the part of the grid belonging to the ``current'' process).
  this%patchSize(1) = max(patchDescriptor%iMax, gridOffset(1) + 1)
  this%patchSize(2) = max(patchDescriptor%jMax, gridOffset(2) + 1)
  this%patchSize(3) = max(patchDescriptor%kMax, gridOffset(3) + 1)
  this%patchSize = min(this%patchSize, gridOffset + gridLocalSize)
  this%patchSize = this%patchSize - this%offset

  ! Reset size and offset if the patch lies outside the part of the grid belonging to the
  ! ``current'' process.
  if (any(this%patchSize < 0)) then
     this%offset = 0
     this%patchSize = 0
  end if

  this%gridLocalSize = gridLocalSize
  this%gridOffset = gridOffset

  this%nPatchPoints = product(this%patchSize)
  this%comm = comm

  this%isCurvilinear = getOption("default/curvilinear", .true.)
  write(key, '(A,I3.3,A)') "grid", this%gridIndex, "/curvilinear"
  this%isCurvilinear = getOption(key, this%isCurvilinear)

  select case (this%patchType)

  case (SPONGE)

     ! Sponge amount.
     this%spongeAmount = getOption("defaults/sponge_amount", 1.0_wp)
     write(key, '(3A)') "patches/",                                                          &
          trim(patchDescriptor%name), "/sponge_amount"
     this%spongeAmount = getOption(key, this%spongeAmount)

     ! Sponge exponent.
     this%spongeExponent = getOption("defaults/sponge_exponent", 2)
     write(key, '(3A)') "patches/",                                                          &
          trim(patchDescriptor%name), "/sponge_exponent"
     this%spongeExponent = getOption(key, this%spongeExponent)

  case (SAT_FAR_FIELD, SAT_BLOCK_INTERFACE)

     ! Inviscid penalty amount.
     this%inviscidPenaltyAmount = getOption("defaults/inviscid_penalty_amount", 1.0_wp)
     write(key, '(3A)') "patches/",                                                          &
          trim(patchDescriptor%name), "/inviscid_penalty_amount"
     this%inviscidPenaltyAmount = getOption(key, this%inviscidPenaltyAmount)
     this%inviscidPenaltyAmount = sign(this%inviscidPenaltyAmount,                           &
          real(this%normalDirection, wp))

  case (SAT_SLIP_WALL, SAT_ISOTHERMAL_WALL, SAT_ADIABATIC_WALL)

     ! Inviscid penalty amount.
     this%inviscidPenaltyAmount = getOption("defaults/inviscid_penalty_amount", 2.0_wp)
     write(key, '(3A)') "patches/",                                                          &
          trim(patchDescriptor%name), "/inviscid_penalty_amount"
     this%inviscidPenaltyAmount = getOption(key, this%inviscidPenaltyAmount)
     this%inviscidPenaltyAmount = sign(this%inviscidPenaltyAmount,                           &
          real(this%normalDirection, wp))

  end select

  if (simulationFlags%viscosityOn) then

     select case (this%patchType)
     case (SAT_FAR_FIELD, SAT_ISOTHERMAL_WALL, SAT_ADIABATIC_WALL)

        ! Viscous penalty amount.
        this%viscousPenaltyAmount = getOption("defaults/viscous_penalty_amount", 1.0_wp)
        write(key, '(3A)') "patches/",                                                       &
             trim(patchDescriptor%name), "/viscous_penalty_amount"
        this%viscousPenaltyAmount = getOption(key, this%viscousPenaltyAmount)
        this%viscousPenaltyAmount = sign(this%viscousPenaltyAmount,                          &
             real(this%normalDirection, wp))

     end select

  end if

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

subroutine addDamping(this, mode, rightHandSide, iblank, solvedVariables, targetVariables)

  ! <<< Derived types >>>
  use Patch_type
  use Region_type, only : FORWARD, ADJOINT

  implicit none

  ! <<< Arguments >>>
  type(t_Patch) :: this
  integer, intent(in) :: mode
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)
  integer, intent(in) :: iblank(:)
  SCALAR_TYPE, intent(in) :: solvedVariables(:,:)
  SCALAR_TYPE, intent(in), optional :: targetVariables(:,:)

  ! <<< Local variables >>>
  integer :: i, j, k, l, gridIndex, patchIndex

  select case (mode)

  case (FORWARD)

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
                 rightHandSide(gridIndex, l) = rightHandSide(gridIndex, l) -                 &
                      this%spongeStrength(patchIndex) *                                      &
                      (solvedVariables(gridIndex, l) - targetVariables(gridIndex, l))
              end do
           end do
        end do
     end do

  case (ADJOINT)

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
                      this%spongeStrength(patchIndex) *                                      &
                      solvedVariables(gridIndex, l)
              end do
           end do
        end do
     end do

  end select

end subroutine addDamping

subroutine addFarFieldPenalty(this, mode, rightHandSide, iblank, nDimensions,                &
     ratioOfSpecificHeats, conservedVariables, targetState)

  ! <<< Derived types >>>
  use Patch_type
  use Region_type, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use CNSHelper

  implicit none

  ! <<< Arguments >>>
  type(t_Patch) :: this
  integer, intent(in) :: mode
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)
  integer, intent(in) :: iblank(:), nDimensions
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
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

  if (allocated(this%viscousFluxes)) then
     allocate(viscousFluxPenalty(this%nPatchPoints, nDimensions + 1))
     do i = 1, this%nPatchPoints
        viscousFluxPenalty(i,:) = matmul(this%viscousFluxes(i,2:nDimensions+2,:) -           &
             this%targetViscousFluxes(i,2:nDimensions+2,:),                                  &
             this%metrics(i,1+nDimensions*(direction-1):nDimensions*direction))
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
                   this%normalDirection, localIncomingJacobianOfInviscidFluxes)
           case (2)
              call computeIncomingJacobianOfInviscidFlux2D(localTargetState,                 &
                   localMetricsAlongNormalDirection, ratioOfSpecificHeats,                   &
                   this%normalDirection, localIncomingJacobianOfInviscidFluxes)
           case (3)
              call computeIncomingJacobianOfInviscidFlux3D(localTargetState,                 &
                   localMetricsAlongNormalDirection, ratioOfSpecificHeats,                   &
                   this%normalDirection, localIncomingJacobianOfInviscidFluxes)
           end select

           rightHandSide(gridIndex,:) = rightHandSide(gridIndex,:) -                         &
                this%inviscidPenaltyAmount * matmul(localIncomingJacobianOfInviscidFluxes,   &
                conservedVariables(gridIndex,:) - localTargetState)

           if (allocated(viscousFluxPenalty)) then
              rightHandSide(gridIndex,2:nDimensions+2) =                                     &
                   rightHandSide(gridIndex,2:nDimensions+2) +                                &
                   this%viscousPenaltyAmount * viscousFluxPenalty(patchIndex,:)
           end if

        end do
     end do
  end do

  SAFE_DEALLOCATE(viscousFluxPenalty)
  SAFE_DEALLOCATE(localIncomingJacobianOfInviscidFluxes)
  SAFE_DEALLOCATE(localMetricsAlongNormalDirection)
  SAFE_DEALLOCATE(localTargetState)

end subroutine addFarFieldPenalty

subroutine addWallPenalty(this, mode, rightHandSide, iblank, nDimensions,                    &
     ratioOfSpecificHeats, conservedVariables)

  ! <<< Derived types >>>
  use Patch_type
  use Region_type, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use CNSHelper

  implicit none

  ! <<< Arguments >>>
  type(t_Patch) :: this
  integer, intent(in) :: mode
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)
  integer, intent(in) :: iblank(:), nDimensions
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
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
                   this%normalDirection, localIncomingJacobianOfInviscidFluxes)
           case (2)
              call computeIncomingJacobianOfInviscidFlux2D(localConservedVariables,          &
                   localMetricsAlongNormalDirection, ratioOfSpecificHeats,                   &
                   this%normalDirection, localIncomingJacobianOfInviscidFluxes)
           case (3)
              call computeIncomingJacobianOfInviscidFlux3D(localConservedVariables,          &
                   localMetricsAlongNormalDirection, ratioOfSpecificHeats,                   &
                   this%normalDirection, localIncomingJacobianOfInviscidFluxes)
           end select

           temp = dot_product(localConservedVariables(2:nDimensions+1),                      &
                localMetricsAlongNormalDirection)

           localInviscidPenalty(1) = 0.0_wp
           localInviscidPenalty(2:nDimensions+1) = localMetricsAlongNormalDirection * temp
           localInviscidPenalty(nDimensions+2) =                                             &
                0.5_wp / localConservedVariables(1) * temp ** 2
           localInviscidPenalty(2:nDimensions+2) = localInviscidPenalty(2:nDimensions+2) /   &
                sum(localMetricsAlongNormalDirection ** 2)
           localInviscidPenalty = matmul(localIncomingJacobianOfInviscidFluxes,              &
                localInviscidPenalty)

           rightHandSide(gridIndex,:) = rightHandSide(gridIndex,:) -                         &
                this%inviscidPenaltyAmount * localInviscidPenalty

        end do
     end do
  end do

  SAFE_DEALLOCATE(localInviscidPenalty)
  SAFE_DEALLOCATE(localIncomingJacobianOfInviscidFluxes)
  SAFE_DEALLOCATE(localMetricsAlongNormalDirection)
  SAFE_DEALLOCATE(localConservedVariables)

end subroutine addWallPenalty
