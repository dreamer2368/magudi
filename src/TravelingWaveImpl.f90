#include "config.h"

subroutine setupTravelingWave(this, region)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use TravelingWave_mod, only : t_TravelingWave

  ! <<< Enumerations >>>
  use State_enum, only : QOI_DUMMY_FUNCTION

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables
  use InputHelper, only : getOption, getRequiredOption


  implicit none

  ! <<< Arguments >>>
  class(t_TravelingWave) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  call this%cleanup()

  call this%setupBase(region%simulationFlags, region%solverOptions)

  SAFE_DEALLOCATE(region%params%buffer)
  allocate(region%params%buffer(1,1))
  region%params%buffer = getOption("traveling_wave/initial_speed",0.0_wp)

end subroutine setupTravelingWave

subroutine cleanupTravelingWave(this)

  ! <<< Derived types >>>
  use TravelingWave_mod, only : t_TravelingWave

  implicit none

  ! <<< Arguments >>>
  class(t_TravelingWave) :: this

  ! <<< Local variables >>>
  integer :: i

  call this%cleanupBase()

end subroutine cleanupTravelingWave

subroutine computeTravelingWaveSpatialDistribution(this, grid, state, F)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State

  use TravelingWave_mod, only : t_TravelingWave

  ! <<< Internal modules >>>
  use CNSHelper, only: transformFluxes

  implicit none

  ! <<< Arguments >>>
  class(t_TravelingWave) :: this
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  SCALAR_TYPE, intent(out) :: F(:,:)
  SCALAR_TYPE, allocatable :: fluxes1(:,:,:), fluxes2(:,:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  assert(allocated(F))
  assert(size(F,1)==grids%nGridPoints)
  assert(size(F,2)==nDimensions+2)
  assert(size(state%conservedVariables,1)==grids%nGridPoints)
  assert(size(state%conservedVariables,2)==nDimensions+2)

  allocate(fluxes1(grid%nGridPoints, nDimensions+2, nDimensions))
  allocate(fluxes2(grid%nGridPoints, nDimensions+2, nDimensions))
  fluxes1 = 0.0_wp
  fluxes1(:,:,1) = state%conservedVariables

  call transformFluxes(nDimensions, fluxes1, grid%metrics, fluxes2, grid%isCurvilinear)

  SAFE_DEALLOCATE(fluxes1)

  do i = 1, nDimensions
     call grid%firstDerivative(i)%apply(fluxes2(:,:,i), grid%localSize)
  end do

  F = sum(fluxes2, dim = 3)
  do i = 1, nDimensions+2
    F(:,i) = F(:,i) * grid%jacobian(:,1)
  end do
  ! F = - this%params(1) * F - state%rightHandSide

  SAFE_DEALLOCATE(fluxes2)

end subroutine computeTravelingWaveSpatialDistribution

function computeTravelingWave(this, region) result(instantaneousFunctional)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use TravelingWave_mod, only : t_TravelingWave

  ! <<< Internal modules >>>
  use Patch_factory, only : computeQuadratureOnPatches

  ! <<< SeungWhan: debugging >>>
  use, intrinsic :: iso_fortran_env, only : output_unit
  use ErrorHandler, only : writeAndFlush

  ! <<< Arguments >>>
  class(t_TravelingWave) :: this
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, ierror
  SCALAR_TYPE, allocatable :: F(:,:)

  ! <<< SeungWhan: message, timeRampFactor >>
  character(len=STRING_LENGTH) :: message

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(size(region%grids) == size(region%states))

  instantaneousFunctional = 0.0_wp

  do i = 1, size(region%grids)

     assert(region%grids(i)%nGridPoints > 0)
     assert(allocated(region%grids(i)%targetMollifier))
     assert(size(region%grids(i)%targetMollifier, 1) == region%grids(i)%nGridPoints)
     assert(size(region%grids(i)%targetMollifier, 2) == 1)

     allocate(F(region%grids(i)%nGridPoints, region%solverOptions%nUnknowns))
     call this%computeSpatialDistribution(region%grids(i),region%states(i),F)
     instantaneousFunctional = instantaneousFunctional + 0.5_wp *                            &
       ! computeQuadratureOnPatches( region%patchFactories, 'COST_TARGET', region%grids(i),    &
       !   sum(( - region%params%buffer(1,1)*F - region%states(i)%rightHandSide)**2, dim=2) *  &
       !                                                 region%grids(i)%targetMollifier(:,1) )
       computeQuadratureOnPatches( region%patchFactories, 'COST_TARGET', region%grids(i),    &
         sum(region%states(i)%rightHandSide**2, dim=2) * region%grids(i)%targetMollifier(:,1) )
     SAFE_DEALLOCATE(F)
  end do

  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE, instantaneousFunctional, 1,                          &
       SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(instantaneousFunctional, 1, SCALAR_TYPE_MPI,                             &
          0, region%grids(i)%comm, ierror)
  end do

  this%cachedValue = instantaneousFunctional

end function computeTravelingWave

subroutine computeTravelingWaveAdjointForcing(this, simulationFlags, solverOptions,          &
     grid, state, patch)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use TravelingWave_mod, only : t_TravelingWave
  use SolverOptions_mod, only : t_SolverOptions
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_TravelingWave) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  class(t_CostTargetPatch) :: patch

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: F(:,:)

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  assert(grid%nGridPoints > 0)
  assert(allocated(grid%targetMollifier))
  assert(size(grid%targetMollifier, 1) == grid%nGridPoints)
  assert(size(grid%targetMollifier, 2) == 1)

  allocate(F(grid%nGridPoints, nDimensions+2))
  call this%computeSpatialDistribution(grid,state,F)
  ! state%adjointVariables = F
  SAFE_DEALLOCATE(F)


  !
  ! allocate(meanPressure(patch%nPatchPoints))
  ! i = grid%index
  ! call patch%collect(this%data_(i)%meanPressure(:,1), meanPressure)
  !
  ! timeRampFactor = 1.0_wp
  ! if (this%useTimeWindow)                                                                    &
  !     timeRampFactor = timeRampFactor                                                        &
  !         * exp( -(state%time-this%timeWindowCenter)**2/2.0_wp/this%timeWindowWidth/this%timeWindowWidth )
  ! ! ideal_mean_pressure = 1.0_wp/solverOptions%ratioOfSpecificHeats
  !
  ! do k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)
  !    do j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  !       do i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
  !          gridIndex = i - patch%gridOffset(1) + patch%gridLocalSize(1) *                    &
  !               (j - 1 - patch%gridOffset(2) + patch%gridLocalSize(2) *                      &
  !               (k - 1 - patch%gridOffset(3)))
  !          if (grid%iblank(gridIndex) == 0) cycle
  !          patchIndex = i - patch%offset(1) + patch%localSize(1) *                           &
  !               (j - 1 - patch%offset(2) + patch%localSize(2) *                              &
  !               (k - 1 - patch%offset(3)))
  !
  !          F = - 2.0_wp * grid%targetMollifier(gridIndex, 1) *                               &
  !               timeRampFactor *                                                             &
  !               (solverOptions%ratioOfSpecificHeats - 1.0_wp) *                              &
  !               (state%pressure(gridIndex, 1) - meanPressure(patchIndex))
  !               ! (state%pressure(gridIndex, 1) - ideal_mean_pressure)
  !
  !          patch%adjointForcing(patchIndex,nDimensions+2) = F
  !          patch%adjointForcing(patchIndex,2:nDimensions+1) =                                &
  !               - state%velocity(gridIndex,:) * F
  !          patch%adjointForcing(patchIndex,1) =                                              &
  !               0.5_wp * sum(state%velocity(gridIndex,:) ** 2) * F
  !
  !       end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
  !    end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  ! end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)
  !
  ! SAFE_DEALLOCATE(meanPressure)

end subroutine computeTravelingWaveAdjointForcing

subroutine addTravelingWaveGradient(this, simulationFlags, solverOptions,          &
                                     grid, state, patch)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use TravelingWave_mod, only : t_TravelingWave
  use SolverOptions_mod, only : t_SolverOptions
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_TravelingWave) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  class(t_CostTargetPatch) :: patch

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, gridIndex, patchIndex

end subroutine addTravelingWaveGradient

function isTravelingWavePatchValid(this, patchDescriptor, gridSize, normalDirection,         &
     extent, simulationFlags, message) result(isPatchValid)

  ! <<< Derived types >>>
  use TravelingWave_mod, only : t_TravelingWave
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_TravelingWave) :: this
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

end function isTravelingWavePatchValid
