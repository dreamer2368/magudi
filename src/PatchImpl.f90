#include "config.h"

module PatchImpl

  implicit none
  public

contains

  subroutine allocateData(this, nDimensions, nUnknowns, simulationFlags)

    ! <<< Derived types >>>
    use Patch_type, only : t_Patch
    use PatchDescriptor_type
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Arguments >>>
    type(t_Patch) :: this
    integer, intent(in) :: nDimensions, nUnknowns
    type(t_SimulationFlags), intent(in) :: simulationFlags

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND

    select case (this%patchType)
    case (SAT_FAR_FIELD, SAT_SLIP_WALL, SAT_ISOTHERMAL_WALL, SAT_ADIABATIC_WALL)
       allocate(this%metrics(this%nPatchPoints, nDimensions ** 2))
       if (this%isLiftMeasured .or. this%isDragMeasured)                                     &
            allocate(this%adjointSource(this%nPatchPoints, nUnknowns))
       if (this%isLiftMeasured) this%liftDirection = (/ 0.0_wp, 1.0_wp, 0.0_wp /)
       if (this%isDragMeasured) this%DragDirection = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
    end select

    select case (this%patchType)

    case (SAT_FAR_FIELD)
       if (simulationFlags%viscosityOn) then
          allocate(this%viscousFluxes(this%nPatchPoints, nUnknowns, nDimensions))
          allocate(this%targetViscousFluxes(this%nPatchPoints, nUnknowns, nDimensions))
       end if

    case (SAT_BLOCK_INTERFACE)
       allocate(this%conservedVariables(this%nPatchPoints, nUnknowns))
       allocate(this%interfaceDataBuffer1(this%nPatchPoints, 2 * nUnknowns - 1))
       allocate(this%interfaceDataBuffer2(this%nPatchPoints, 2 * nUnknowns - 1))

    case (SAT_ISOTHERMAL_WALL)
       allocate(this%wallTemperature(this%nPatchPoints))

    case (SPONGE)
       allocate(this%spongeStrength(this%nPatchPoints), source = 0.0_wp)

    case (ACTUATOR)
       allocate(this%gradient(this%nPatchPoints, 1))

    case (SOLENOIDAL_EXCITATION)
       allocate(this%solenoidalExcitationStrength(this%nPatchPoints), source = 0.0_wp)

    end select

  end subroutine allocateData

end module PatchImpl

subroutine setupPatch(this, index, nDimensions, patchDescriptor,                             &
     comm, gridOffset, gridLocalSize, nUnknowns, simulationFlags)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_type, only : t_Patch
  use PatchDescriptor_type
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Private members >>>
  use PatchImpl, only : allocateData

  ! <<< Internal modules >>>
  use Patch_mod, only : cleanupPatch
  use InputHelper, only : getOption, getRequiredOption
  use SolenoidalExcitation_mod, only : setupSolenoidalExcitation

  implicit none

  ! <<< Arguments >>>
  type(t_Patch) :: this
  integer, intent(in) :: index, nDimensions
  type(t_PatchDescriptor) :: patchDescriptor
  integer, intent(in) :: comm, gridOffset(3), gridLocalSize(3), nUnknowns
  type(t_SimulationFlags), intent(in) :: simulationFlags

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(SCALAR_KIND) :: solenoidalExcitationOrigin(3), solenoidalExcitationSpeed(3)
  character(len = STRING_LENGTH) :: key
  SCALAR_TYPE :: wallTemperature

  assert(index > 0)
  assert_key(nDimensions, (1, 2, 3))
  assert(all(gridOffset >= 0))
  assert(all(gridLocalSize > 0))
  assert(nUnknowns > 0)

  assert(patchDescriptor%gridIndex > 0)
  assert(abs(patchDescriptor%normalDirection) <= nDimensions)
  assert(patchDescriptor%iMin > 0)
  assert(patchDescriptor%jMin > 0)
  assert(patchDescriptor%kMin > 0)
  assert(patchDescriptor%iMax >= patchDescriptor%iMin)
  assert(patchDescriptor%jMax >= patchDescriptor%jMin)
  assert(patchDescriptor%kMax >= patchDescriptor%kMin)
  assert(len_trim(patchDescriptor%name) > 0)

  assert_key(patchDescriptor%patchType, ( \
  SPONGE,                \
  ACTUATOR,              \
  CONTROL_TARGET,        \
  SOLENOIDAL_EXCITATION, \
  SAT_FAR_FIELD,         \
  SAT_SLIP_WALL,         \
  SAT_ISOTHERMAL_WALL,   \
  SAT_ADIABATIC_WALL,    \
  SAT_BLOCK_INTERFACE))

  call cleanupPatch(this)

  this%index = index
  this%normalDirection = patchDescriptor%normalDirection
  this%gridIndex = patchDescriptor%gridIndex
  this%patchType = patchDescriptor%patchType

  this%extent = (/ patchDescriptor%iMin, patchDescriptor%iMax,                               &
       patchDescriptor%jMin, patchDescriptor%jMax,                                           &
       patchDescriptor%kMin, patchDescriptor%kMax /)

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

  write(key, '(3A)') "patches/", trim(patchDescriptor%name), "/"

  select case (this%patchType)

  case (SPONGE)

     ! Sponge amount.
     this%spongeAmount = getOption("defaults/sponge_amount", 1.0_wp)
     this%spongeAmount = getOption(trim(key) // "sponge_amount", this%spongeAmount)

     ! Sponge exponent.
     this%spongeExponent = getOption("defaults/sponge_exponent", 2)
     this%spongeExponent = getOption(trim(key) // "sponge_exponent", this%spongeExponent)

  case (SOLENOIDAL_EXCITATION)

     solenoidalExcitationOrigin(1) = getOption(trim(key) // "x", 0.0_wp)
     solenoidalExcitationOrigin(2) = getOption(trim(key) // "y", 0.0_wp)
     solenoidalExcitationOrigin(3) = getOption(trim(key) // "z", 0.0_wp)

     solenoidalExcitationSpeed(1) = getOption(trim(key) // "u", 0.0_wp)
     solenoidalExcitationSpeed(2) = getOption(trim(key) // "v", 0.0_wp)
     solenoidalExcitationSpeed(3) = getOption(trim(key) // "w", 0.0_wp)

     call setupSolenoidalExcitation(this%solenoidalExcitation, this%comm,                    &
          getOption(trim(key) // "number_of_modes", 1),                                      &
          solenoidalExcitationOrigin, solenoidalExcitationSpeed,                             &
          getOption(trim(key) // "amplitude", 0.01_wp),                                      &
          getOption(trim(key) // "most_unstable_frequency", 0.1_wp),                         &
          getOption(trim(key) // "radius", 1.0_wp),                                          &
          getOption(trim(key) // "seed", -1))

  case (SAT_FAR_FIELD, SAT_BLOCK_INTERFACE)

     ! Inviscid penalty amount.
     this%inviscidPenaltyAmount = getOption("defaults/inviscid_penalty_amount", 1.0_wp)
     this%inviscidPenaltyAmount = getOption(trim(key) // "inviscid_penalty_amount",          &
          this%inviscidPenaltyAmount)
     this%inviscidPenaltyAmount = sign(this%inviscidPenaltyAmount,                           &
          real(this%normalDirection, wp))

  case (SAT_SLIP_WALL, SAT_ISOTHERMAL_WALL, SAT_ADIABATIC_WALL)

     ! Inviscid penalty amount.
     this%inviscidPenaltyAmount = getOption("defaults/inviscid_penalty_amount", 2.0_wp)
     this%inviscidPenaltyAmount = getOption(trim(key) // "inviscid_penalty_amount",          &
          this%inviscidPenaltyAmount)
     this%inviscidPenaltyAmount = sign(this%inviscidPenaltyAmount,                           &
          real(this%normalDirection, wp))

     ! Measure lift/drag on this patch?
     if (.not. simulationFlags%predictionOnly) then
        this%isLiftMeasured = getOption(trim(key) // "measure_lift", .false.)
        this%isDragMeasured = getOption(trim(key) // "measure_drag", .false.)
     end if

  end select

  if (simulationFlags%viscosityOn) then

     select case (this%patchType)

     case (SAT_FAR_FIELD, SAT_ISOTHERMAL_WALL, SAT_ADIABATIC_WALL)

        ! Viscous penalty amount.
        this%viscousPenaltyAmount = getOption("defaults/viscous_penalty_amount", 1.0_wp)
        this%viscousPenaltyAmount = getOption(trim(key) // "viscous_penalty_amount",         &
             this%viscousPenaltyAmount)
        this%viscousPenaltyAmount = sign(this%viscousPenaltyAmount,                          &
             real(this%normalDirection, wp))

     end select

  end if

  if (this%nPatchPoints > 0) call allocateData(this, nDimensions, nUnknowns, simulationFlags)

  if (allocated(this%wallTemperature) .and. .not. simulationFlags%useTargetState) then
     call getRequiredOption(trim(key) // "temperature", wallTemperature)
     this%wallTemperature = wallTemperature
  end if

end subroutine setupPatch

subroutine cleanupPatch(this)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_type, only : t_Patch
  use PatchDescriptor_type

  ! <<< Internal modules >>>
  use SolenoidalExcitation_mod, only : cleanupSolenoidalExcitation

  implicit none

  ! <<< Arguments >>>
  type(t_Patch) :: this

  ! <<< Local variables >>>
  integer :: ierror

  call cleanupSolenoidalExcitation(this%solenoidalExcitation)

  SAFE_DEALLOCATE(this%viscousFluxes)
  SAFE_DEALLOCATE(this%targetViscousFluxes)
  SAFE_DEALLOCATE(this%wallTemperature)
  SAFE_DEALLOCATE(this%adjointSource)
  SAFE_DEALLOCATE(this%metrics)
  SAFE_DEALLOCATE(this%spongeStrength)
  SAFE_DEALLOCATE(this%gradient)
  SAFE_DEALLOCATE(this%solenoidalExcitationStrength)

  if (this%comm /= MPI_COMM_NULL) call MPI_Comm_free(this%comm, ierror)
  this%comm = MPI_COMM_NULL
  this%commOfConformingPatch = MPI_COMM_NULL

end subroutine cleanupPatch

subroutine collectScalarAtPatch_(this, gridArray, patchArray)

  ! <<< Derived types >>>
  use Patch_type, only : t_Patch

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  ! <<< Arguments >>>
  type(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: gridArray(:)
  SCALAR_TYPE, intent(out) :: patchArray(:)

  ! <<< Local variables >>>
  integer :: i, j, k, patchIndex, localIndex

  assert(all(this%gridLocalSize > 0) .and. size(gridArray, 1) == product(this%gridLocalSize))
  assert(all(this%patchSize >= 0) .and. size(patchArray, 1) == product(this%patchSize)) 

  call startTiming("collectAtPatch")

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

  call endTiming("collectAtPatch")

end subroutine collectScalarAtPatch_

subroutine collectVectorAtPatch_(this, gridArray, patchArray)

  ! <<< Derived types >>>
  use Patch_type, only : t_Patch

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  ! <<< Arguments >>>
  type(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: gridArray(:,:)
  SCALAR_TYPE, intent(out) :: patchArray(:,:)

  ! <<< Local variables >>>
  integer :: i, j, k, l, patchIndex, localIndex

  assert(all(this%gridLocalSize > 0) .and. size(gridArray, 1) == product(this%gridLocalSize))
  assert(all(this%patchSize >= 0) .and. size(patchArray, 1) == product(this%patchSize))
  assert(size(gridArray, 2) > 0)
  assert(size(patchArray, 2) == size(gridArray, 2))

  call startTiming("collectAtPatch")

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

  call endTiming("collectAtPatch")

end subroutine collectVectorAtPatch_

subroutine collectTensorAtPatch_(this, gridArray, patchArray)

  ! <<< Derived types >>>
  use Patch_type, only : t_Patch

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  ! <<< Arguments >>>
  type(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: gridArray(:,:,:)
  SCALAR_TYPE, intent(out) :: patchArray(:,:,:)

  ! <<< Local variables >>>
  integer :: i, j, k, l, m, patchIndex, localIndex

  assert(all(this%gridLocalSize > 0) .and. size(gridArray, 1) == product(this%gridLocalSize))
  assert(all(this%patchSize >= 0) .and. size(patchArray, 1) == product(this%patchSize))
  assert(size(gridArray, 2) > 0)
  assert(size(patchArray, 2) == size(gridArray, 2))
  assert(size(gridArray, 3) > 0)
  assert(size(patchArray, 3) == size(gridArray, 3))

  call startTiming("collectAtPatch")

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

  call endTiming("collectAtPatch")

end subroutine collectTensorAtPatch_

subroutine addDamping(this, mode, rightHandSide, iblank, solvedVariables, targetVariables)

  ! <<< Derived types >>>
  use Patch_type, only : t_Patch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

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

  assert(all(this%gridLocalSize > 0))
  assert(size(rightHandSide, 1) == product(this%gridLocalSize))
  assert(size(iblank) == product(this%gridLocalSize))
  assert(size(solvedVariables, 1) == product(this%gridLocalSize))
  assert(size(rightHandSide, 2) > 0)
  assert(size(solvedVariables, 2) == size(rightHandSide, 2))

  assert_key(mode, (FORWARD, ADJOINT))

  call startTiming("addDamping")

  select case (mode)

  case (FORWARD)

     assert(present(targetVariables))
     assert(size(targetVariables, 1) == product(this%gridLocalSize))
     assert(size(targetVariables, 2) == size(rightHandSide, 2))

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
              end do !... i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
           end do !... j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
        end do !... k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
     end do !... l = 1, size(rightHandSide, 2)

  case (ADJOINT)

     assert(.not. present(targetVariables))

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
              end do !... i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
           end do !... j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
        end do !... k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
     end do !... l = 1, size(rightHandSide, 2)

  end select

  call endTiming("addDamping")

end subroutine addDamping

subroutine addFarFieldPenalty(this, mode, rightHandSide, iblank, nDimensions,                &
     ratioOfSpecificHeats, conservedVariables, targetState, adjointVariables)

  ! <<< Derived types >>>
  use Patch_type, only : t_Patch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use CNSHelper
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  type(t_Patch) :: this
  integer, intent(in) :: mode
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)
  integer, intent(in) :: iblank(:), nDimensions
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
  SCALAR_TYPE, intent(in) :: conservedVariables(:,:), targetState(:,:)
  SCALAR_TYPE, intent(in), optional :: adjointVariables(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nUnknowns, direction, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: localTargetState(:), metricsAlongNormalDirection(:),           &
       incomingJacobianOfInviscidFlux(:,:), viscousFluxPenalty(:,:)

  assert(all(this%gridLocalSize > 0))
  assert(size(rightHandSide, 1) == product(this%gridLocalSize))
  assert(size(iblank) == product(this%gridLocalSize))
  assert_key(nDimensions, (1, 2, 3))
  assert(ratioOfSpecificHeats > 1.0_wp)
  assert(size(conservedVariables, 1) == product(this%gridLocalSize))
  assert(size(targetState, 1) == product(this%gridLocalSize))
  assert(size(rightHandSide, 2) > 0)
  assert(size(conservedVariables, 2) == size(rightHandSide, 2))
  assert(size(targetState, 2) == size(rightHandSide, 2))

  assert_key(mode, (FORWARD, ADJOINT))

#ifdef DEBUG
  select case (mode)
  case (FORWARD)
     assert(.not. present(adjointVariables))
  case (ADJOINT)
     assert(present(adjointVariables))
     assert(size(adjointVariables, 1) == product(this%gridLocalSize))
     assert(size(adjointVariables, 2) == size(rightHandSide, 2))
  end select
#endif

  call startTiming("addFarFieldPenalty")

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  nUnknowns = size(rightHandSide, 2)
  assert(nUnknowns == nDimensions + 2)

  allocate(localTargetState(nUnknowns))
  allocate(metricsAlongNormalDirection(nDimensions))
  allocate(incomingJacobianOfInviscidFlux(nUnknowns, nUnknowns))

  if (allocated(this%viscousFluxes)) then
     allocate(viscousFluxPenalty(this%nPatchPoints, nUnknowns - 1))
     do i = 1, this%nPatchPoints
        viscousFluxPenalty(i,:) = matmul(this%viscousFluxes(i,2:nUnknowns,:) -               &
             this%targetViscousFluxes(i,2:nUnknowns,:),                                      &
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
           metricsAlongNormalDirection =                                                     &
                this%metrics(patchIndex,1+nDimensions*(direction-1):nDimensions*direction)

           select case (nDimensions)
           case (1)
              call computeIncomingJacobianOfInviscidFlux1D(localTargetState,                 &
                   metricsAlongNormalDirection, ratioOfSpecificHeats,                        &
                   this%normalDirection, incomingJacobianOfInviscidFlux)
           case (2)
              call computeIncomingJacobianOfInviscidFlux2D(localTargetState,                 &
                   metricsAlongNormalDirection, ratioOfSpecificHeats,                        &
                   this%normalDirection, incomingJacobianOfInviscidFlux)
           case (3)
              call computeIncomingJacobianOfInviscidFlux3D(localTargetState,                 &
                   metricsAlongNormalDirection, ratioOfSpecificHeats,                        &
                   this%normalDirection, incomingJacobianOfInviscidFlux)
           end select !... nDimensions

           select case (mode)
           case (FORWARD)

              rightHandSide(gridIndex,:) = rightHandSide(gridIndex,:) -                      &
                   this%inviscidPenaltyAmount * matmul(incomingJacobianOfInviscidFlux,       &
                   conservedVariables(gridIndex,:) - localTargetState)

              if (allocated(viscousFluxPenalty)) then
                 rightHandSide(gridIndex,2:nUnknowns) =                                      &
                      rightHandSide(gridIndex,2:nUnknowns) +                                 &
                      this%viscousPenaltyAmount * viscousFluxPenalty(patchIndex,:)
              end if

           case (ADJOINT)

              rightHandSide(gridIndex,:) = rightHandSide(gridIndex,:) +                      &
                   this%inviscidPenaltyAmount *                                              &
                   matmul(transpose(incomingJacobianOfInviscidFlux),                         &
                   adjointVariables(gridIndex,:))

              ! TODO: add viscous far-field penalties for adjoint variables.

           end select !... mode

        end do !... i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)

  SAFE_DEALLOCATE(viscousFluxPenalty)
  SAFE_DEALLOCATE(incomingJacobianOfInviscidFlux)
  SAFE_DEALLOCATE(metricsAlongNormalDirection)
  SAFE_DEALLOCATE(localTargetState)

  call endTiming("addFarFieldPenalty")

end subroutine addFarFieldPenalty

subroutine addWallPenalty(this, mode, rightHandSide, iblank, nDimensions,                    &
     ratioOfSpecificHeats, conservedVariables, adjointVariables)

  ! <<< Derived types >>>
  use Patch_type, only : t_Patch
  use PatchDescriptor_type, only : SAT_ISOTHERMAL_WALL, SAT_ADIABATIC_WALL

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use CNSHelper
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  type(t_Patch) :: this
  integer, intent(in) :: mode
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)
  integer, intent(in) :: iblank(:), nDimensions
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
  SCALAR_TYPE, intent(in) :: conservedVariables(:,:)
  SCALAR_TYPE, intent(in), optional :: adjointVariables(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, nUnknowns, direction, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: localConservedVariables(:),                                    &
       metricsAlongNormalDirection(:), incomingJacobianOfInviscidFlux(:,:),                  &
       inviscidPenalty(:), viscousPenalty(:), deltaNormalVelocity(:),                        &
       deltaInviscidPenalty(:,:), deltaIncomingJacobianOfInviscidFlux(:,:,:)
  SCALAR_TYPE :: normalVelocity

  assert(all(this%gridLocalSize > 0))
  assert(size(rightHandSide, 1) == product(this%gridLocalSize))
  assert(size(iblank) == product(this%gridLocalSize))
  assert_key(nDimensions, (1, 2, 3))
  assert(ratioOfSpecificHeats > 1.0_wp)
  assert(size(conservedVariables, 1) == product(this%gridLocalSize))
  assert(size(rightHandSide, 2) > 0)
  assert(size(conservedVariables, 2) == size(rightHandSide, 2))

  assert_key(mode, (FORWARD, ADJOINT))

#ifdef DEBUG
  select case (mode)
  case (FORWARD)
     assert(.not. present(adjointVariables))
  case (ADJOINT)
     assert(present(adjointVariables))
     assert(size(adjointVariables, 1) == product(this%gridLocalSize))
     assert(size(adjointVariables, 2) == size(rightHandSide, 2))
  end select
#endif

  call startTiming("addWallPenalty")

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  nUnknowns = size(rightHandSide, 2)
  assert(nUnknowns == nDimensions + 2)

  allocate(localConservedVariables(nUnknowns))
  allocate(metricsAlongNormalDirection(nDimensions))
  allocate(incomingJacobianOfInviscidFlux(nUnknowns, nUnknowns))
  allocate(inviscidPenalty(nUnknowns))
  if (mode == ADJOINT) then
     allocate(deltaNormalVelocity(nUnknowns))
     allocate(deltaInviscidPenalty(nUnknowns, nUnknowns))
     allocate(deltaIncomingJacobianOfInviscidFlux(nUnknowns, nUnknowns, nUnknowns))
  end if

  select case (this%patchType)
  case (SAT_ISOTHERMAL_WALL, SAT_ADIABATIC_WALL)
     allocate(viscousPenalty(nUnknowns - 1))
  end select

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
           metricsAlongNormalDirection =                                                     &
                this%metrics(patchIndex,1+nDimensions*(direction-1):nDimensions*direction)

           normalVelocity = dot_product(localConservedVariables(2:nDimensions+1) /           &
                localConservedVariables(1), metricsAlongNormalDirection /                    &
                sqrt(sum(metricsAlongNormalDirection ** 2)))
           inviscidPenalty(1) = 0.0_wp
           inviscidPenalty(2:nDimensions+1) = metricsAlongNormalDirection * normalVelocity / &
                sqrt(sum(metricsAlongNormalDirection ** 2))
           inviscidPenalty(nDimensions+2) =                                                  &
                0.5_wp * localConservedVariables(1) * normalVelocity ** 2

           select case (this%patchType)
           case (SAT_ISOTHERMAL_WALL)
              viscousPenalty(1:nDimensions+1) = localConservedVariables(2:nDimensions+2)
              viscousPenalty(nDimensions+1) = viscousPenalty(nDimensions+1) -                &
                   localConservedVariables(1) *                                              &
                   this%wallTemperature(patchIndex) / ratioOfSpecificHeats
           end select

           select case (mode)

           case (FORWARD)

              select case (nDimensions)
              case (1)
                 call computeIncomingJacobianOfInviscidFlux1D(localConservedVariables,       &
                      metricsAlongNormalDirection, ratioOfSpecificHeats,                     &
                      this%normalDirection, incomingJacobianOfInviscidFlux)
              case (2)
                 call computeIncomingJacobianOfInviscidFlux2D(localConservedVariables,       &
                      metricsAlongNormalDirection, ratioOfSpecificHeats,                     &
                      this%normalDirection, incomingJacobianOfInviscidFlux)
              case (3)
                 call computeIncomingJacobianOfInviscidFlux3D(localConservedVariables,       &
                      metricsAlongNormalDirection, ratioOfSpecificHeats,                     &
                      this%normalDirection, incomingJacobianOfInviscidFlux)
              end select !... nDimensions

           case (ADJOINT)

              deltaNormalVelocity(1) = - normalVelocity / localConservedVariables(1)
              deltaNormalVelocity(2:nDimensions+1) = metricsAlongNormalDirection /           &
                sqrt(sum(metricsAlongNormalDirection ** 2)) / localConservedVariables(1)
              deltaNormalVelocity(nDimensions+2) = 0.0_wp
              deltaInviscidPenalty(1,:) = 0.0_wp
              do l = 1, nDimensions
                 deltaInviscidPenalty(l+1,:) =                                               &
                      metricsAlongNormalDirection(l) * deltaNormalVelocity /                 &
                      sqrt(sum(metricsAlongNormalDirection ** 2))
              end do
              deltaInviscidPenalty(nDimensions+2,:) =                                        &
                   localConservedVariables(1) * normalVelocity * deltaNormalVelocity
              deltaInviscidPenalty(nDimensions+2,1) =                                        &
                   deltaInviscidPenalty(nDimensions+2,1) + 0.5_wp * normalVelocity ** 2

              select case (nDimensions)
              case (1)
                 call computeIncomingJacobianOfInviscidFlux1D(localConservedVariables,       &
                      metricsAlongNormalDirection, ratioOfSpecificHeats,                     &
                      this%normalDirection, incomingJacobianOfInviscidFlux,                  &
                      deltaIncomingJacobianOfInviscidFlux)
              case (2)
                 call computeIncomingJacobianOfInviscidFlux2D(localConservedVariables,       &
                      metricsAlongNormalDirection, ratioOfSpecificHeats,                     &
                      this%normalDirection, incomingJacobianOfInviscidFlux,                  &
                      deltaIncomingJacobianOfInviscidFlux)
              case (3)
                 call computeIncomingJacobianOfInviscidFlux3D(localConservedVariables,       &
                      metricsAlongNormalDirection, ratioOfSpecificHeats,                     &
                      this%normalDirection, incomingJacobianOfInviscidFlux,                  &
                      deltaIncomingJacobianOfInviscidFlux)
              end select !... nDimensions

           end select !... mode

           select case (mode)

           case (FORWARD)

              rightHandSide(gridIndex,:) = rightHandSide(gridIndex,:) -                      &
                   this%inviscidPenaltyAmount * matmul(incomingJacobianOfInviscidFlux,       &
                   inviscidPenalty)

              select case (this%patchType)
              case (SAT_ISOTHERMAL_WALL, SAT_ADIABATIC_WALL)
                 rightHandSide(gridIndex,2:nUnknowns) =                                      &
                      rightHandSide(gridIndex,2:nUnknowns) -                                 &
                      this%viscousPenaltyAmount * viscousPenalty
              end select

           case (ADJOINT)

              do l = 1, nUnknowns
                 rightHandSide(gridIndex,l) = rightHandSide(gridIndex,l) +                   &
                      this%inviscidPenaltyAmount *                                           &
                      dot_product(adjointVariables(gridIndex,:),                             &
                      matmul(deltaIncomingJacobianOfInviscidFlux(:,:,l), inviscidPenalty) +  &
                      matmul(incomingJacobianOfInviscidFlux, deltaInviscidPenalty(:,l)))
              end do

           end select

        end do !... i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)

  SAFE_DEALLOCATE(viscousPenalty)
  SAFE_DEALLOCATE(deltaIncomingJacobianOfInviscidFlux)
  SAFE_DEALLOCATE(deltaInviscidPenalty)
  SAFE_DEALLOCATE(deltaNormalVelocity)
  SAFE_DEALLOCATE(inviscidPenalty)
  SAFE_DEALLOCATE(incomingJacobianOfInviscidFlux)
  SAFE_DEALLOCATE(metricsAlongNormalDirection)
  SAFE_DEALLOCATE(localConservedVariables)

  call endTiming("addWallPenalty")

end subroutine addWallPenalty

subroutine updateSolenoidalExcitationStrength(this, coordinates, iblank)

  ! <<< Derived types >>>
  use Patch_type, only : t_Patch

  implicit none

  ! <<< Arguments >>>
  type(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: coordinates(:,:)
  integer, intent(in) :: iblank(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, gridIndex, patchIndex

  assert(all(this%gridLocalSize > 0))
  assert(size(coordinates, 1) == product(this%gridLocalSize))
  assert(size(iblank) == product(this%gridLocalSize))
  assert(this%solenoidalExcitation%gaussianFactor >= 0.0_wp)

  nDimensions = size(coordinates, 2)
  assert_key(nDimensions, (1, 2, 3))

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

           this%solenoidalExcitationStrength(patchIndex) =                                   &
                this%solenoidalExcitation%amplitude *                                        &
                exp(- this%solenoidalExcitation%gaussianFactor *                             &
                sum((real(coordinates(gridIndex,:), wp) -                                    &
                this%solenoidalExcitation%location(1:nDimensions)) ** 2))

        end do !... i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)

end subroutine updateSolenoidalExcitationStrength

subroutine addSolenoidalExcitation(this, coordinates, iblank, time, rightHandSide)

  ! <<< Derived types >>>
  use Patch_type, only : t_Patch

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  type(t_Patch) :: this
  SCALAR_TYPE, intent(in) :: coordinates(:,:)
  integer, intent(in) :: iblank(:)
  real(SCALAR_KIND), intent(in) :: time
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, gridIndex, patchIndex, nDimensions
  real(wp) :: location(3), speed(3), gaussianFactor, localStrength, localOrigin(3), temp(4)
  real(wp), allocatable :: angularFrequencies(:), phases(:,:)

  call startTiming("addSolenoidalExcitation")

  nDimensions = size(coordinates, 2)

  gaussianFactor = this%solenoidalExcitation%gaussianFactor
  location = this%solenoidalExcitation%location
  speed = this%solenoidalExcitation%speed

  allocate(angularFrequencies(this%solenoidalExcitation%nModes))
  allocate(phases(this%solenoidalExcitation%nModes, nDimensions))

  angularFrequencies = real(this%solenoidalExcitation%angularFrequencies, wp)
  phases = real(this%solenoidalExcitation%phases(:,1:nDimensions), wp)

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

           localStrength = this%solenoidalExcitationStrength(patchIndex)
           localOrigin(1:nDimensions) = real(coordinates(gridIndex,:), wp) -                 &
                location(1:nDimensions) - speed(1:nDimensions) * time

           do l = 1, this%solenoidalExcitation%nModes

              temp(1) = sin(angularFrequencies(l) * localOrigin(1) + phases(l,1))
              temp(2) = cos(angularFrequencies(l) * localOrigin(1) + phases(l,1))
              temp(3) = sin(angularFrequencies(l) * localOrigin(2) + phases(l,2))
              temp(4) = cos(angularFrequencies(l) * localOrigin(2) + phases(l,2))

              rightHandSide(gridIndex,2) = rightHandSide(gridIndex,2) +                      &
                   localStrength * temp(1) *                                                 &
                   (angularFrequencies(l) * temp(4) - 2.0_wp * gaussianFactor *              &
                   (coordinates(gridIndex,2) - location(2)) * temp(3))
              rightHandSide(gridIndex,3) = rightHandSide(gridIndex,3) +                      &
                   localStrength * temp(3) *                                                 &
                   (angularFrequencies(l) * temp(2) - 2.0_wp * gaussianFactor *              &
                   (coordinates(gridIndex,1) - location(1)) * temp(1))

           end do !... l = 1, this%solenoidalExcitation%nModes

        end do !... i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)

  SAFE_DEALLOCATE(phases)
  SAFE_DEALLOCATE(angularFrequencies)

  call endTiming("addSolenoidalExcitation")

end subroutine addSolenoidalExcitation

subroutine updatePatchConnectivity(this, patchData)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_type, only : t_Patch
  use PatchDescriptor_type

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  type(t_Patch) :: this
  type(t_PatchDescriptor), intent(in) :: patchData(:)

  ! <<< Local variables >>>
  character(len = STRING_LENGTH) :: key, val
  integer :: i

  select case (this%patchType)

  case (SAT_BLOCK_INTERFACE)

     this%indexOfConformingPatch = 0
     this%commOfConformingPatch = MPI_COMM_NULL

     write(key, '(3A)') "patches/", trim(patchData(this%index)%name), "/conforms_with"
     val = getOption(key, "")

     if (len_trim(val) == 0) then

        do i = 1, size(patchData)
           if (i == this%index) cycle
           write(key, '(3A)') "patches/", trim(patchData(i)%name), "/conforms_with"
           val = getOption(key, "")
           if (trim(val) == trim(patchData(this%index)%name) .and.                           &
                patchData(i)%patchType == SAT_BLOCK_INTERFACE .and.                          &
                patchData(i)%normalDirection == -this%normalDirection) then
              this%indexOfConformingPatch = i
              exit
           end if
        end do !... i = 1, size(patchData)

     else

        do i = 1, size(patchData)
           if (i == this%index) cycle
           if (trim(val) == trim(patchData(i)%name) .and.                                    &
                patchData(i)%patchType == SAT_BLOCK_INTERFACE .and.                          &
                patchData(i)%normalDirection == -this%normalDirection) then
              this%indexOfConformingPatch = i
              exit
           end if
        end do !... i = 1, size(patchData)

     end if

  end select

end subroutine updatePatchConnectivity
