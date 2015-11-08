#include "config.h"

subroutine setupFlameTemperature(this, region)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use FlameTemperature_mod, only : t_FlameTemperature

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables
  use InputHelper, only : getOption, getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_FlameTemperature) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, ierror
  character(len = STRING_LENGTH) :: key

  assert(allocated(region%states))
  assert(size(region%states) > 0)
  assert(region%solverOptions%nSpecies > 0)
  assert(region%states(1)%combustion%nReactions > 0)

  call this%cleanup()

  call this%setupBase(region%simulationFlags, region%solverOptions)

  allocate(this%data_(size(region%globalGridSizes, 2)))

  do i = 1, size(this%data_)
     do j = 1, size(region%grids)
        if (region%grids(j)%index == i) then
           allocate(this%data_(i)%flameTemperature(1, 1))
           region%states(j)%dummyFunction => this%data_(i)%flameTemperature
           exit
        end if
     end do
     call MPI_Barrier(region%comm, ierror)
  end do

  do i = 1, size(this%data_)
     this%data_(i)%flameTemperature = 1.0_wp / ( (region%solverOptions%ratioOfSpecificHeats -&
          1.0_wp) * (1.0_wp - region%states(1)%combustion%heatRelease) )
  end do

  write(key, '(A)') "flame_temperature/"
  this%weightBurnRegion = getOption(trim(key) // "weight_burn_region", .false.)
  if (this%weightBurnRegion) call getRequiredOption(trim(key) // "burn_radius",              &
       this%burnRadius, region%comm)

  this%useTimeRamp = getOption(trim(key) // "use_time_ramp", .false.)
  if (this%useTimeRamp) then
     call getRequiredOption(trim(key) // "ramp_width", this%rampWidthInverse)
     this%rampWidthInverse = 1.0_wp / this%rampWidthInverse
     call getRequiredOption(trim(key) // "ramp_peak", this%rampPeak)
  end if

end subroutine setupFlameTemperature

subroutine cleanupFlameTemperature(this)

  ! <<< Derived types >>>
  use FlameTemperature_mod, only : t_FlameTemperature

  implicit none

  ! <<< Arguments >>>
  class(t_FlameTemperature) :: this

  ! <<< Local variables >>>
  integer :: i

  call this%cleanupBase()

  if (allocated(this%data_)) then
     do i = 1, size(this%data_)
        if (associated(this%data_(i)%flameTemperature)) deallocate(this%data_(i)%flameTemperature)
        nullify(this%data_(i)%flameTemperature)
     end do
  end if
  SAFE_DEALLOCATE(this%data_)

end subroutine cleanupFlameTemperature

function computeFlameTemperature(this, region) result(instantaneousFunctional)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use FlameTemperature_mod, only : t_FlameTemperature
  use Patch_mod, only : t_Patch
  use CostTargetPatch_mod, only : t_CostTargetPatch

  ! <<< Arguments >>>
  class(t_FlameTemperature) :: this
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, H2, O2, ierror
  SCALAR_TYPE :: YF0, YO0, s, Z, Zst, gaussianFactor, referenceTemperature, timeRampFactor
  SCALAR_TYPE, allocatable :: F(:,:), W(:)

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(allocated(region%patchFactories))
  assert(size(region%grids) == size(region%states))
  assert(region%solverOptions%nSpecies >= 2)

  instantaneousFunctional = 0.0_wp
  referenceTemperature = 1.0_wp / (region%solverOptions%ratioOfSpecificHeats - 1.0_wp)

  timeRampFactor = 1.0_wp
  if (this%useTimeRamp)                                                                      &
       timeRampFactor = exp(-0.5_wp * (region%states(1)%time - this%rampPeak)**2 *           &
       this%rampWidthInverse**2)

  do i = 1, size(region%grids)

     assert(region%grids(i)%nGridPoints > 0)
     assert(allocated(region%grids(i)%targetMollifier))
     assert(size(region%grids(i)%targetMollifier, 1) == region%grids(i)%nGridPoints)
     assert(size(region%grids(i)%targetMollifier, 2) == 1)
     assert(allocated(region%states(i)%massFraction))
     assert(size(region%states(i)%massFraction, 1) == region%grids(i)%nGridPoints)
     assert(size(region%states(i)%massFraction, 2) == region%solverOptions%nSpecies)
     assert(allocated(region%states(i)%temperature))
     assert(size(region%states(i)%temperature, 1) == region%grids(i)%nGridPoints)
     assert(size(region%states(i)%temperature, 2) == 1)

     j = region%grids(i)%index

     assert(associated(this%data_(j)%flameTemperature))
     assert(size(this%data_(j)%flameTemperature, 1) == 1)
     assert(size(this%data_(j)%flameTemperature, 2) == 1)

     allocate(W(region%grids(i)%nGridPoints))

     if (this%weightBurnRegion) then

        H2 = region%states(i)%combustion%H2
        O2 = region%states(i)%combustion%O2
        YF0 = region%states(i)%combustion%Y0(H2)
        YO0 = region%states(i)%combustion%Y0(O2)
        s = region%states(i)%combustion%stoichiometricRatio
        Zst = 1.0_wp / (1.0_wp + s * YF0 / YO0)
        gaussianFactor = -0.5_wp / this%burnRadius**2

        do k = 1, region%grids(i)%nGridPoints
           Z = ( region%states(i)%massFraction(k, H2) -                                      &
                region%states(i)%massFraction(k, O2) / s  + YO0 / s ) / (YF0 + YO0 / s)
           W(k) = region%grids(i)%targetMollifier(k, 1) * exp(gaussianFactor * (Z - Zst) **2)
        end do

     else

        W = region%grids(i)%targetMollifier(:, 1)

     end if

     allocate(F(region%grids(i)%nGridPoints, 1))

     F = (region%states(i)%temperature - referenceTemperature) /                             &
          (this%data_(j)%flameTemperature(1,1) - referenceTemperature)

     instantaneousFunctional = instantaneousFunctional +                                     &
          region%grids(i)%computeInnerProduct(F, F, W)

     SAFE_DEALLOCATE(W)
     SAFE_DEALLOCATE(F)

  end do

  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE, instantaneousFunctional, 1,                          &
       SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(instantaneousFunctional, 1, SCALAR_TYPE_MPI,                             &
          0, region%grids(i)%comm, ierror)
  end do

  this%auxilaryFunctional = instantaneousFunctional
  instantaneousFunctional = instantaneousFunctional * timeRampFactor
  this%cachedValue = instantaneousFunctional

end function computeFlameTemperature

subroutine computeFlameTemperatureAdjointForcing(this, simulationFlags, solverOptions,       &
     grid, state, patch)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use FlameTemperature_mod, only : t_FlameTemperature
  use SolverOptions_mod, only : t_SolverOptions
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use CNSHelper

  ! <<< Enumerations >>>
  use SolverOptions_enum

  implicit none

  ! <<< Arguments >>>
  class(t_FlameTemperature) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  class(t_CostTargetPatch) :: patch

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, nSpecies, gridIndex, patchIndex, H2, O2
  SCALAR_TYPE :: flameTemperature, referenceTemperature, F, W, Z, Zst, s, YF0, YO0,          &
       gaussianFactor, timeRampFactor, timeStep
  SCALAR_TYPE, dimension(:), allocatable :: deltaTemperature

  nDimensions = grid%nDimensions
  nSpecies = solverOptions%nSpecies
  assert_key(nDimensions, (1, 2, 3))
  assert(grid%nGridPoints > 0)
  assert(nSpecies >= 2)

  i = grid%index
  flameTemperature = this%data_(i)%flameTemperature(1,1)
  referenceTemperature = 1.0_wp / (solverOptions%ratioOfSpecificHeats - 1.0_wp)

  timeRampFactor = 1.0_wp
  if (this%useTimeRamp) then
     timeStep = state%computeTimeStepSize(grid, simulationFlags, solverOptions)
     timeRampFactor = exp(-0.5_wp * (state%time - 0.5_wp * timeStep - this%rampPeak)**2 *    &
          this%rampWidthInverse**2)
  end if

  if (this%weightBurnRegion) then

     assert(allocated(state%massFraction))
     assert(size(state%massFraction, 1) == grid%nGridPoints)
     assert(size(state%massFraction, 2) == nSpecies)

     gaussianFactor = -0.5_wp / this%burnRadius**2
     H2 = state%combustion%H2
     O2 = state%combustion%O2
     YF0 = state%combustion%Y0(H2)
     YO0 = state%combustion%Y0(O2)
     s = state%combustion%stoichiometricRatio
     Zst = 1.0_wp / (1.0_wp + s * YF0 / YO0)

  end if

  allocate(deltaTemperature(nDimensions + nSpecies + 2))

  do k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)
     do j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
        do i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
           gridIndex = i - patch%gridOffset(1) + patch%gridLocalSize(1) *                    &
                (j - 1 - patch%gridOffset(2) + patch%gridLocalSize(2) *                      &
                (k - 1 - patch%gridOffset(3)))

           if (grid%iblank(gridIndex) == 0) cycle

           patchIndex = i - patch%offset(1) + patch%localSize(1) *                           &
                (j - 1 - patch%offset(2) + patch%localSize(2) *                              &
                (k - 1 - patch%offset(3)))

           patch%adjointForcing(patchIndex,:) = 0.0_wp

           if (this%weightBurnRegion) then

              Z = ( state%massFraction(gridIndex, H2) -                                      &
                   state%massFraction(gridIndex, O2) / s  + YO0 / s ) / (YF0 + YO0 / s)
              W = grid%targetMollifier(gridIndex, 1) * exp(gaussianFactor * (Z - Zst) **2)

              ! First apply -2*W*(T-T0)/(Tf-T0)^2*dT/dQ.
              call computeDeltaVariables(nDimensions, nSpecies,                              &
                   state%conservedVariables(gridIndex,:), solverOptions%equationOfState,     &
                   solverOptions%ratioOfSpecificHeats, solverOptions%molecularWeightInverse, &
                   deltaTemperature = deltaTemperature)

              F = - 2.0_wp * W * timeRampFactor *                                            &
                   (state%temperature(gridIndex, 1) - referenceTemperature) /                &
                   (flameTemperature - referenceTemperature)**2

              patch%adjointForcing(patchIndex,:) = F * deltaTemperature

              ! Now apply -((T-Tf)/(Tf-T0))^2*dW/dQ.
              F = W * timeRampFactor *                                                       &
                   ( (state%temperature(gridIndex, 1) - referenceTemperature) /              &
                   (flameTemperature - referenceTemperature) )**2 *                          &
                   ( (Z - Zst) / this%burnRadius**2 ) * state%specificVolume(gridIndex,1) /  &
                   (YF0 + YO0 / s)

              patch%adjointForcing(patchIndex,1) = patch%adjointForcing(patchIndex,1) +      &
                   (state%massFraction(gridIndex, O2) / s -                                  &
                   state%massFraction(gridIndex, H2)) * F
              patch%adjointForcing(patchIndex,nDimensions+2+H2) = F
              patch%adjointForcing(patchIndex,nDimensions+2+O2) = - F / s

           else

              call computeDeltaVariables(nDimensions, nSpecies,                              &
                   state%conservedVariables(gridIndex,:), solverOptions%equationOfState,     &
                   solverOptions%ratioOfSpecificHeats, solverOptions%molecularWeightInverse, &
                   deltaTemperature = deltaTemperature)

              F = - 2.0_wp * grid%targetMollifier(gridIndex, 1) * timeRampFactor *           &
                   (state%temperature(gridIndex, 1) - referenceTemperature) /                &
                   (flameTemperature - referenceTemperature)**2

              patch%adjointForcing(patchIndex,:) = F * deltaTemperature

           end if

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
     end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

  SAFE_DEALLOCATE(deltaTemperature)

end subroutine computeFlameTemperatureAdjointForcing

function isFlameTemperaturePatchValid(this, patchDescriptor, gridSize, normalDirection,         &
     extent, simulationFlags, message) result(isPatchValid)

  ! <<< Derived types >>>
  use FlameTemperature_mod, only : t_FlameTemperature
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_FlameTemperature) :: this
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  integer, intent(in) :: gridSize(:), normalDirection, extent(6)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Result >>>
  logical :: isPatchValid

  ! <<< Local variables >>>
  integer :: i, n

  isPatchValid = .false.

  n = size(gridSize)

  do i = 1, size(gridSize)
     if (extent((i-1)*2+1) == extent((i-1)*2+2)) n = n - 1
  end do

  if (n /= size(gridSize)) then
     write(message, '(2(A,I0.0),A)') "Expected a ", size(gridSize),                          &
          "D patch, but extent represents a ", n, "D patch!"
     return
  end if

  isPatchValid = .true.

end function isFlameTemperaturePatchValid
