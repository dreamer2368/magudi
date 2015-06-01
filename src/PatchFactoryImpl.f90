#include "config.h"

subroutine connectPatch(this, patchTarget, patchType, createNew)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Patch_factory, only : t_PatchFactory
  use SpongePatch_mod, only : t_SpongePatch
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use AdiabaticWall_mod, only : t_AdiabaticWall
  use FarFieldPatch_mod, only : t_FarFieldPatch
  use IsothermalWall_mod, only : t_IsothermalWall
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use ImpenetrableWall_mod, only : t_ImpenetrableWall
  use JetExcitationPatch_mod, only : t_JetExcitationPatch
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch
  use SolenoidalExcitationPatch_mod, only : t_SolenoidalExcitationPatch
  use GaussianIgnitionPatch_mod, only : t_GaussianIgnitionPatch

  implicit none

  ! <<< Arguments >>>
  class(t_PatchFactory) :: this
  class(t_Patch), pointer, intent(out) :: patchTarget
  character(len = *), intent(in), optional :: patchType
  logical, intent(in), optional :: createNew

  ! <<< Local variables >>>
  logical :: createNew_

  createNew_ = .false.
  if (present(createNew)) createNew_ = createNew

  if (present(patchType) .and. .not. (associated(this%patch) .and. .not. createNew_)) then

     if (associated(this%patch)) deallocate(this%patch)
     nullify(this%patch)

     this%patchType = patchType

     select case (trim(patchType))

     case ('SPONGE')
        allocate(t_SpongePatch :: this%patch)

     case ('ACTUATOR')
        allocate(t_ActuatorPatch :: this%patch)

     case ('SAT_ADIABATIC_WALL')
        allocate(t_AdiabaticWall :: this%patch)

     case ('SAT_FAR_FIELD')
        allocate(t_FarFieldPatch :: this%patch)

     case ('SAT_ISOTHERMAL_WALL')
        allocate(t_IsothermalWall :: this%patch)

     case ('COST_TARGET')
        allocate(t_CostTargetPatch :: this%patch)

     case ('SAT_SLIP_WALL')
        allocate(t_ImpenetrableWall :: this%patch)

     case ('SAT_BLOCK_INTERFACE')
        allocate(t_BlockInterfacePatch :: this%patch)

     case ('SOLENOIDAL_EXCITATION')
        allocate(t_SolenoidalExcitationPatch :: this%patch)

     case ('GAUSSIAN_IGNITION')
        allocate(t_GaussianIgnitionPatch :: this%patch)

     case ('JET_EXCITATION')
        allocate(t_JetExcitationPatch :: this%patch)

     case default
        this%patchType = ""

     end select

  end if

  nullify(patchTarget)
  if (.not. associated(this%patch)) return
  patchTarget => this%patch

end subroutine connectPatch

subroutine cleanupPatchFactory(this)

  use Patch_factory, only : t_PatchFactory

  class(t_PatchFactory) :: this

  if (associated(this%patch)) then
     call this%patch%cleanup()
     deallocate(this%patch)
  end if
  nullify(this%patch)

end subroutine cleanupPatchFactory

function queryPatchTypeExists(patchFactories, patchType,                                     &
     gridIndex, normalDirection) result(patchTypeExists)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Patch_factory, only : t_PatchFactory

  implicit none

  ! <<< Arguments >>>
  type(t_PatchFactory), allocatable :: patchFactories(:)
  character(len = *), intent(in) :: patchType
  integer, intent(in), optional :: gridIndex, normalDirection

  ! <<< Result >>>
  logical :: patchTypeExists

  ! <<< Local variables >>>
  integer :: i
  class(t_Patch), pointer :: patch => null()

#ifdef DEBUG
  if (present(normalDirection)) then
     assert_key(normalDirection, (1, 2, 3))
  end if
#endif

  patchTypeExists = .false.

  if (allocated(patchFactories)) then
     do i = 1, size(patchFactories)
        if (trim(patchType) /= trim(patchFactories(i)%patchType)) cycle

        call patchFactories(i)%connect(patch, patchType)

        if (associated(patch)) then
           if (present(gridIndex)) then
              if (patch%gridIndex /= gridIndex) cycle
           end if
           if (present(normalDirection)) then
              if (abs(patch%normalDirection) /= normalDirection) cycle
           end if
           patchTypeExists = .true.
        end if

     end do
  end if

end function queryPatchTypeExists

subroutine computeSpongeStrengths(patchFactories, grid)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use Patch_mod, only : t_Patch
  use Patch_factory, only : t_PatchFactory
  use SpongePatch_mod, only : t_SpongePatch

  ! <<< Internal modules >>>
  use MPIHelper, only : gatherAlongDirection
  use Patch_factory, only : queryPatchTypeExists

  implicit none

  ! <<< Arguments >>>
  type(t_PatchFactory), allocatable :: patchFactories(:)
  class(t_Grid), intent(in) :: grid

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, iMin, iMax, jMin, jMax, kMin, kMax, direction, nDimensions, ierror
  logical :: spongesExistAlongDirection
  class(t_Patch), pointer :: patch => null()
  SCALAR_TYPE, dimension(:,:), allocatable :: coordinateDerivatives, arcLength,              &
       globalArcLengthsAlongDirection
  real(wp), allocatable :: curveLengthIntegrand(:)

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  assert(grid%nGridPoints > 0)
  assert(all(grid%localSize > 0) .and. product(grid%localSize) == grid%nGridPoints)
  assert(all(grid%offset >= 0))
  assert(all(grid%globalSize >= grid%localSize))

  do direction =  1, nDimensions

     ! Check if there are sponge patches along direction `direction`.
     spongesExistAlongDirection = queryPatchTypeExists(patchFactories,                       &
          'SPONGE', grid%index, direction)
     call MPI_Allreduce(MPI_IN_PLACE, spongesExistAlongDirection, 1, MPI_LOGICAL,            &
          MPI_LOR, grid%comm, ierror) !... reduce across grid-level processes.
     if (.not. spongesExistAlongDirection) cycle

     ! Compute local arc length.
     allocate(arcLength(grid%nGridPoints, 1))
     allocate(coordinateDerivatives(grid%nGridPoints, nDimensions))
     call grid%computeCoordinateDerivatives(direction, coordinateDerivatives)
     arcLength(:,1) = sqrt(sum(coordinateDerivatives ** 2, dim = 2))
     SAFE_DEALLOCATE(coordinateDerivatives) !... no longer needed.

     ! Gather arc length along direction `direction`.
     allocate(globalArcLengthsAlongDirection(grid%nGridPoints / grid%localSize(direction) *  &
          grid%globalSize(direction), 1))
     call gatherAlongDirection(grid%comm, arcLength, grid%localSize,                         &
          direction, grid%offset(direction), globalArcLengthsAlongDirection)
     SAFE_DEALLOCATE(arcLength) !... no longer needed.

     allocate(curveLengthIntegrand(grid%globalSize(direction)))

     if (allocated(patchFactories)) then
        do l = 1, size(patchFactories)

           call patchFactories(l)%connect(patch)
           if (.not. associated(patch)) cycle
           if (patch%gridIndex /= grid%index .or. patch%nPatchPoints <= 0 .or.               &
                abs(patch%normalDirection) /= direction) cycle
           select type (patch)
           class is (t_SpongePatch)

              select case (direction)

              case (1)

                 iMin = patch%extent(1)
                 iMax = patch%extent(2)

                 do k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)
                    do j = patch%offset(2) + 1,                                              &
                         patch%offset(2) + patch%localSize(2)

                       do i = 1, grid%globalSize(1)
                          curveLengthIntegrand(i) =                                          &
                               real(globalArcLengthsAlongDirection(i +                       &
                               grid%globalSize(1) * (j - 1 - grid%offset(2) +                &
                               grid%localSize(2) * (k - 1 - grid%offset(3))), 1), wp)
                       end do

                       if (patch%normalDirection > 0) then
                          do i = patch%offset(1) + 1,                                        &
                               patch%offset(1) + patch%localSize(1)
                             patch%spongeStrength(i - patch%offset(1) +                      &
                                  patch%localSize(1) * (j - 1 - patch%offset(2) +            &
                                  patch%localSize(2) * (k - 1 - patch%offset(3)))) =         &
                                  sum(curveLengthIntegrand(iMin : i - 1)) /                  &
                                  sum(curveLengthIntegrand(iMin : iMax - 1))
                          end do
                       else
                          do i = patch%offset(1) + 1,                                        &
                               patch%offset(1) + patch%localSize(1)
                             patch%spongeStrength(i - patch%offset(1) +                      &
                                  patch%localSize(1) * (j - 1 - patch%offset(2) +            &
                                  patch%localSize(2) * (k - 1 - patch%offset(3)))) =         &
                                  sum(curveLengthIntegrand(i + 1 : iMax)) /                  &
                                  sum(curveLengthIntegrand(iMin + 1 : iMax))
                          end do
                       end if

                    end do !... j = patch%offset(2) + 1,                                     &
                    !...       patch%offset(2) + patch%localSize(2)
                 end do !... k = patch%offset(3) + 1,                                        &
                 !...       patch%offset(3) + patch%localSize(3)

              case (2)

                 jMin = patch%extent(3)
                 jMax = patch%extent(4)

                 do k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)
                    do i = patch%offset(1) + 1,                                              &
                         patch%offset(1) + patch%localSize(1)

                       do j = 1, grid%globalSize(2)
                          curveLengthIntegrand(j) =                                          &
                               real(globalArcLengthsAlongDirection(i - grid%offset(1) +      &
                               grid%localSize(1) * (j - 1 +                                  &
                               grid%globalSize(2) * (k - 1 - grid%offset(3))), 1), wp)
                       end do

                       if (patch%normalDirection > 0) then
                          do j = patch%offset(2) + 1,                                        &
                               patch%offset(2) + patch%localSize(2)
                             patch%spongeStrength(i - patch%offset(1) +                      &
                                  patch%localSize(1) * (j - 1 - patch%offset(2) +            &
                                  patch%localSize(2) * (k - 1 - patch%offset(3)))) =         &
                                  sum(curveLengthIntegrand(jMin : j - 1)) /                  &
                                  sum(curveLengthIntegrand(jMin : jMax - 1))
                          end do
                       else
                          do j = patch%offset(2) + 1,                                        &
                               patch%offset(2) + patch%localSize(2)
                             patch%spongeStrength(i - patch%offset(1) +                      &
                                  patch%localSize(1) * (j - 1 - patch%offset(2) +            &
                                  patch%localSize(2) * (k - 1 - patch%offset(3)))) =         &
                                  sum(curveLengthIntegrand(j + 1 : jMax)) /                  &
                                  sum(curveLengthIntegrand(jMin + 1 : jMax))
                          end do
                       end if

                    end do !... i = patch%offset(1) + 1,                                     &
                    !...       patch%offset(1) + patch%localSize(1)
                 end do !... k = patch%offset(3) + 1,                                        &
                 !...       patch%offset(3) + patch%localSize(3)

              case (3)

                 kMin = patch%extent(5)
                 kMax = patch%extent(6)

                 do j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
                    do i = patch%offset(1) + 1,                                              &
                         patch%offset(1) + patch%localSize(1)

                       do k = 1, grid%globalSize(3)
                          curveLengthIntegrand(k) =                                          &
                               real(globalArcLengthsAlongDirection(i - grid%offset(1) +      &
                               grid%localSize(1) * (j - 1 - grid%offset(2) +                 &
                               grid%localSize(2) * (k - 1)), 1), wp)
                       end do

                       if (patch%normalDirection > 0) then
                          do k = patch%offset(3) + 1,                                        &
                               patch%offset(3) + patch%localSize(3)
                             patch%spongeStrength(i - patch%offset(1) +                      &
                                  patch%localSize(1) * (j - 1 - patch%offset(2) +            &
                                  patch%localSize(2) * (k - 1 - patch%offset(3)))) =         &
                                  sum(curveLengthIntegrand(kMin : k - 1)) /                  &
                                  sum(curveLengthIntegrand(kMin : kMax - 1))
                          end do
                       else
                          do k = patch%offset(3) + 1,                                        &
                               patch%offset(3) + patch%localSize(3)
                             patch%spongeStrength(i - patch%offset(1) +                      &
                                  patch%localSize(1) * (j - 1 - patch%offset(2) +            &
                                  patch%localSize(2) * (k - 1 - patch%offset(3)))) =         &
                                  sum(curveLengthIntegrand(k + 1 : kMax)) /                  &
                                  sum(curveLengthIntegrand(kMin + 1 : kMax))
                          end do
                       end if

                    end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
                 end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)

              end select !... select case (direction)

              patch%spongeStrength = patch%spongeAmount *                                    &
                   (1.0_wp - patch%spongeStrength) ** real(patch%spongeExponent, wp)

           end select !... select type (patch)

        end do !... l = 1, size(patchFactories)
     end if !... allocated(patchFactories)

     SAFE_DEALLOCATE(curveLengthIntegrand)
     SAFE_DEALLOCATE(globalArcLengthsAlongDirection)

  end do

end subroutine computeSpongeStrengths

function computeQuadratureOnPatches(patchFactories, patchType,                               &
     grid, integrand) result(integral)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use Patch_mod, only : t_Patch
  use Patch_factory, only : t_PatchFactory

  implicit none

  ! <<< Arguments >>>
  type(t_PatchFactory), allocatable :: patchFactories(:)
  character(len = *), intent(in) :: patchType
  class(t_Grid), intent(in) :: grid
  SCALAR_TYPE, intent(in) :: integrand(:)

  ! <<< Result >>>
  SCALAR_TYPE :: integral

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l
  SCALAR_TYPE, allocatable :: mask(:)
  class(t_Patch), pointer :: patch => null()

  assert(grid%nGridPoints > 0)
  assert(all(grid%localSize > 0) .and. product(grid%localSize) == grid%nGridPoints)
  assert(all(grid%offset >= 0))

  assert(size(integrand) == grid%nGridPoints)

  assert(allocated(grid%iblank))
  assert(size(grid%iblank) == grid%nGridPoints)

  allocate(mask(grid%nGridPoints))
  mask = 0.0_wp

  if (allocated(patchFactories)) then
     do l = 1, size(patchFactories)

        if (patchFactories(l)%patchType /= trim(patchType)) cycle
        call patchFactories(l)%connect(patch)
        if (.not. associated(patch)) cycle
        if (patch%gridIndex /= grid%index) cycle

        assert(all(patch%gridLocalSize == grid%localSize))
        assert(all(patch%gridOffset == grid%offset))

        do k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)
           do j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
              do i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
                 mask(i - grid%offset(1) + grid%localSize(1) * (j - 1 - grid%offset(2) +     &
                      grid%localSize(2) * (k - 1 - grid%offset(3)))) = 1.0_wp
              end do
           end do
        end do

     end do !... l = 1, size(patchFactories)
  end if !... allocated(patchFactories)

  where (grid%iblank == 0)
     mask = 0.0_wp
  end where

  integral = grid%computeInnerProduct(mask, integrand)

  SAFE_DEALLOCATE(mask)

end function computeQuadratureOnPatches

subroutine updatePatchFactories(patchFactories, simulationFlags, solverOptions, grid, state)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use Patch_mod, only : t_Patch
  use State_mod, only : t_State
  use Patch_factory, only : t_PatchFactory
  use FarFieldPatch_mod, only : t_FarFieldPatch
  use SolverOptions_mod, only : t_SolverOptions
  use IsothermalWall_mod, only : t_IsothermalWall
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use CNSHelper
  use Patch_factory, only : queryPatchTypeExists
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  type(t_PatchFactory), allocatable :: patchFactories(:)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions, ierror
  logical :: flag
  class(t_Patch), pointer :: patch => null()
  SCALAR_TYPE, allocatable :: gridNorm(:,:), targetTemperature(:), targetViscousFluxes(:,:,:)

  call startTiming("updatePatches")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  if (allocated(patchFactories)) then
     do i = 1, size(patchFactories)

        call patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        if (patch%gridIndex /= grid%index .or. patch%nPatchPoints <= 0) cycle

        allocate(gridNorm(grid%nGridPoints, 1))

        select type (patch)
        class is (t_CostTargetPatch)
           gridNorm = 1.0_wp
           do j = 1, nDimensions
              if (j /= abs(patch%normalDirection))                                           &
                   call grid%firstDerivative(j)%applyNorm(gridNorm, grid%localSize)
           end do
           call patch%collect(gridNorm, patch%norm)
        end select

        SAFE_DEALLOCATE(gridNorm)

     end do
  end if

  if (simulationFlags%viscosityOn .and. simulationFlags%useTargetState) then

     flag = queryPatchTypeExists(patchFactories, 'SAT_ISOTHERMAL_WALL', grid%index)
     call MPI_Allreduce(MPI_IN_PLACE, flag, 1, MPI_LOGICAL, MPI_LOR, grid%comm, ierror)
     if (flag) then
        allocate(targetTemperature(grid%nGridPoints))
        call computeDependentVariables(nDimensions, solverOptions%nSpecies,                  &
             state%targetState, solverOptions%ratioOfSpecificHeats,                          &
             temperature = targetTemperature)
     end if

     if (allocated(patchFactories)) then
        do i = 1, size(patchFactories)

           call patchFactories(i)%connect(patch)
           if (.not. associated(patch)) cycle
           if (patch%gridIndex /= grid%index .or. patch%nPatchPoints <= 0) cycle

           select type (patch)
           class is (t_IsothermalWall)
              call patch%collect(targetTemperature, patch%temperature)
              call computeTransportVariables(solverOptions%nSpecies, patch%temperature,      &
                   solverOptions%powerLawExponent, solverOptions%bulkViscosityRatio,         &
                   solverOptions%ratioOfSpecificHeats, solverOptions%reynoldsNumberInverse,  &
                   solverOptions%prandtlNumberInverse, solverOptions%schmidtNumberInverse,   &
                   patch%dynamicViscosity, patch%secondCoefficientOfViscosity,               &
                   patch%thermalDiffusivity, patch%massDiffusivity)
           end select

        end do
     end if

     SAFE_DEALLOCATE(targetTemperature)

     flag = queryPatchTypeExists(patchFactories, 'SAT_FAR_FIELD', grid%index)
     call MPI_Allreduce(MPI_IN_PLACE, flag, 1, MPI_LOGICAL, MPI_LOR, grid%comm, ierror)
     if (flag) then
        allocate(targetViscousFluxes(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))
        call state%update(grid, simulationFlags, solverOptions, state%targetState)
        call computeCartesianViscousFluxes(nDimensions, solverOptions%nSpecies,              &
             state%velocity, state%massFraction, state%stressTensor, state%heatFlux,         &
             state%speciesFlux, targetViscousFluxes)
     end if

     if (allocated(patchFactories)) then
        do i = 1, size(patchFactories)

           call patchFactories(i)%connect(patch)
           if (.not. associated(patch)) cycle
           if (patch%gridIndex /= grid%index .or. patch%nPatchPoints <= 0) cycle

           select type (patch)
           class is (t_FarFieldPatch)
              call patch%collect(targetViscousFluxes, patch%viscousFluxes)
           end select

        end do
     end if

     SAFE_DEALLOCATE(targetViscousFluxes)

  end if

  call endTiming("updatePatches")

end subroutine updatePatchFactories
