#include "config.h"

subroutine setupDensityGradient(this, region)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use DensityGradient_mod, only : t_DensityGradient

  ! <<< Enumerations >>>
  use State_enum, only : QOI_DUMMY_FUNCTION

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables
  use InputHelper, only : getOption, getRequiredOption


  implicit none

  ! <<< Arguments >>>
  class(t_DensityGradient) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, ierror

  assert(allocated(region%states))
  assert(size(region%states) > 0)

  nDimensions = region%grids(1)%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  call this%cleanup()

  call this%setupBase(region%simulationFlags, region%solverOptions)

  this%useTimeWindow = getOption("functional/use_time_window",.false.)
  if (this%useTimeWindow) then
    call getRequiredOption("functional/time_window_center",this%timeWindowCenter)
    call getRequiredOption("functional/time_window_width",this%timeWindowWidth)
  end if

  this%secondComponent = getOption("density_gradient_component2", 1)

  allocate(this%data_(size(region%globalGridSizes, 2)))
  do i = 1, size(this%data_)
     do j = 1, size(region%grids)
        if (region%grids(j)%index == i) then
           allocate(this%data_(i)%adjointVector(region%grids(j)%nGridPoints,1))
           this%data_(i)%adjointVector = region%grids(j)%targetMollifier

           ! NOTE: adjoint derivative is transposed derivative PRE-MULTIPLIED by norm.
           k = this%secondComponent
           call region%grids(j)%adjointFirstDerivative(k)%apply(this%data_(i)%adjointVector,  &
                                                                region%grids(j)%localSize)
           exit
        end if
     end do
     call MPI_Barrier(region%comm, ierror)
  end do

end subroutine setupDensityGradient

subroutine cleanupDensityGradient(this)

  ! <<< Derived types >>>
  use DensityGradient_mod, only : t_DensityGradient

  implicit none

  ! <<< Arguments >>>
  class(t_DensityGradient) :: this

  ! <<< Local variables >>>
  integer :: i

  call this%cleanupBase()

  if (allocated(this%data_)) then
     do i = 1, size(this%data_)
        if (associated(this%data_(i)%adjointVector)) deallocate(this%data_(i)%adjointVector)
        nullify(this%data_(i)%adjointVector)
     end do
  end if
  SAFE_DEALLOCATE(this%data_)

end subroutine cleanupDensityGradient

subroutine computeDensityGradientSpatialDistribution(this, grid, state, F)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use DensityGradient_mod, only : t_DensityGradient

  ! <<< Internal modules >>>
  use CNSHelper, only : computeCartesianInviscidFluxes, computeCartesianViscousFluxes, transformFluxes

  ! <<< Arguments >>>
  class(t_DensityGradient) :: this
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  SCALAR_TYPE, intent(out) :: F(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: k, nDimensions, nUnknowns

  assert(size(F, 1) == grid%nGridPoints)
  assert(size(F, 2) == 1)

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))
  nUnknowns = nDimensions + 2

  assert(grid%nGridPoints > 0)
  assert(allocated(grid%targetMollifier))
  assert(size(grid%targetMollifier, 1) == grid%nGridPoints)
  assert(size(grid%targetMollifier, 2) == 1)

  F(:,1) = state%conservedVariables(:,1)
  k = this%secondComponent
  call grid%firstDerivative(k)%apply(F, grid%localSize)

end subroutine

function computeDensityGradient(this, region) result(instantaneousFunctional)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use DensityGradient_mod, only : t_DensityGradient

  ! <<< Internal modules >>>
  use Patch_factory, only : computeQuadratureOnPatches

  ! <<< SeungWhan: debugging >>>
  use, intrinsic :: iso_fortran_env, only : output_unit
  use ErrorHandler, only : writeAndFlush

  ! <<< Arguments >>>
  class(t_DensityGradient) :: this
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, ierror
  SCALAR_TYPE, allocatable :: F(:,:)

  real(wp) :: timeRampFactor

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(size(region%grids) == size(region%states))

  instantaneousFunctional = 0.0_wp

  ! SeungWhan: compute timeRampFactor
  timeRampFactor = 1.0_wp
  if (this%useTimeWindow)                                                                    &
      timeRampFactor = timeRampFactor                                                        &
          * exp( -(region%states(1)%time-this%timeWindowCenter)**2/2.0_wp/this%timeWindowWidth/this%timeWindowWidth )

  do i = 1, size(region%grids)

     allocate(F(region%grids(i)%nGridPoints, 1))

     call this%computeSpatialDistribution(region%grids(i),region%states(i),F)

     instantaneousFunctional = instantaneousFunctional +                                     &
          timeRampFactor *                                                                   &
          computeQuadratureOnPatches(region%patchFactories, 'COST_TARGET',                   &
                  region%grids(i), F(:,1) * region%grids(i)%targetMollifier(:,1))

  end do

  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE, instantaneousFunctional, 1,                          &
       SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(instantaneousFunctional, 1, SCALAR_TYPE_MPI,                             &
          0, region%grids(i)%comm, ierror)
  end do

  this%cachedValue = instantaneousFunctional

end function computeDensityGradient

subroutine computeDensityGradientAdjointForcing(this, simulationFlags, solverOptions,          &
     grid, state, patch)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use DensityGradient_mod, only : t_DensityGradient
  use SolverOptions_mod, only : t_SolverOptions
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use CNSHelper!, only : transformFluxes

  implicit none

  ! <<< Arguments >>>
  class(t_DensityGradient) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  class(t_CostTargetPatch) :: patch

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, nUnknowns, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: adjointVector(:)
  real(wp) :: timeRampFactor

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))
  nUnknowns = nDimensions + 2

  allocate(adjointVector(patch%nPatchPoints))
  i = grid%index
  call patch%collect(this%data_(i)%adjointVector(:,1), adjointVector)

  timeRampFactor = 1.0_wp
  if (this%useTimeWindow)                                                                    &
      timeRampFactor = timeRampFactor                                                        &
          * exp( -(state%time-this%timeWindowCenter)**2/2.0_wp/this%timeWindowWidth/this%timeWindowWidth )

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

           patch%adjointForcing(patchIndex,1) = - adjointVector(patchIndex)                  &
                                                   * timeRampFactor

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
     end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

  SAFE_DEALLOCATE(adjointVector)

end subroutine computeDensityGradientAdjointForcing

function isDensityGradientPatchValid(this, patchDescriptor, gridSize, normalDirection,         &
     extent, simulationFlags, message) result(isPatchValid)

  ! <<< Derived types >>>
  use DensityGradient_mod, only : t_DensityGradient
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_DensityGradient) :: this
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

end function isDensityGradientPatchValid
