#include "config.h"

subroutine setupLighthillTensorComponent(this, region)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use LighthillTensorComponent_mod, only : t_LighthillTensorComponent

  ! <<< Enumerations >>>
  use State_enum, only : QOI_DUMMY_FUNCTION

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables
  use InputHelper, only : getOption, getRequiredOption


  implicit none

  ! <<< Arguments >>>
  class(t_LighthillTensorComponent) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix, message

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

  ! this%viscosityOn = region%simulationFlags%viscosityOn
  this%viscosityOn = .false.

  this%firstComponent = getOption("lighthill_tensor_component/first", 1)
  this%secondComponent = getOption("lighthill_tensor_component/second", 1)
  assert_key(this%firstComponent,(1,2,3))
  assert_key(this%secondComponent,(1,2,3))

end subroutine setupLighthillTensorComponent

subroutine cleanupLighthillTensorComponent(this)

  ! <<< Derived types >>>
  use LighthillTensorComponent_mod, only : t_LighthillTensorComponent

  implicit none

  ! <<< Arguments >>>
  class(t_LighthillTensorComponent) :: this

  ! <<< Local variables >>>
  integer :: i

  call this%cleanupBase()

end subroutine cleanupLighthillTensorComponent

subroutine computeLighthillTensorComponentSpatialDistribution(this, grid, state, F)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use LighthillTensorComponent_mod, only : t_LighthillTensorComponent

  ! <<< Internal modules >>>
  use CNSHelper, only : computeCartesianInviscidFluxes, computeCartesianViscousFluxes, transformFluxes

  ! <<< Arguments >>>
  class(t_LighthillTensorComponent) :: this
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  SCALAR_TYPE, intent(out) :: F(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, nUnknowns, ierror
  logical :: computeViscousFlux = .false.
  SCALAR_TYPE, allocatable :: fluxes1(:,:,:), fluxes2(:,:,:)

  computeViscousFlux = this%viscosityOn

  assert(size(F, 1) == grid%nGridPoints)
  assert(size(F, 2) == 1)

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))
  nUnknowns = nDimensions + 2

  assert(grid%nGridPoints > 0)
  assert(allocated(grid%targetMollifier))
  assert(size(grid%targetMollifier, 1) == grid%nGridPoints)
  assert(size(grid%targetMollifier, 2) == 1)
  assert(allocated(state%pressure))
  assert(size(state%pressure, 1) == grid%nGridPoints)
  assert(size(state%pressure, 2) == 1)
  assert(allocated(state%velocity))
  assert(size(state%velocity, 1) == grid%nGridPoints)
  assert(size(state%velocity, 2) == nDimensions)
  assert(allocated(state%stressTensor))
  assert(size(state%stressTensor, 1) == grid%nGridPoints)
  assert(size(state%stressTensor, 2) == nDimensions ** 2)
  assert(allocated(state%heatFlux))
  assert(size(state%heatFlux, 1) == grid%nGridPoints)
  assert(size(state%heatFlux, 2) == nDimensions)

  allocate(fluxes1(grid%nGridPoints, nUnknowns, nDimensions))
  allocate(fluxes2(grid%nGridPoints, nUnknowns, nDimensions))

  ! Compute Cartesian form of inviscid fluxes.
  call computeCartesianInviscidFluxes(nDimensions, state%conservedVariables,                 &
       state%velocity, state%pressure(:,1), fluxes1)

  ! Compute Cartesian form of viscous fluxes if viscous terms are included and computed using
  ! repeated first derivatives.
  if (computeViscousFlux) then
     call computeCartesianViscousFluxes(nDimensions, state%velocity,                         &
          state%stressTensor, state%heatFlux, fluxes2)
     fluxes1 = fluxes1 - fluxes2 !... Cartesian form of total fluxes.
  end if
  do k = 1, nDimensions
    fluxes1(:,k+1,k) = fluxes1(:,k+1,k) - state%conservedVariables(:,1)
  end do

  F(:,1) = fluxes1(:,this%firstComponent+1,this%secondComponent)
  SAFE_DEALLOCATE(fluxes1)
  SAFE_DEALLOCATE(fluxes2)

end subroutine

function computeLighthillTensorComponent(this, region) result(instantaneousFunctional)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use LighthillTensorComponent_mod, only : t_LighthillTensorComponent

  ! <<< Internal modules >>>
  use Patch_factory, only : computeQuadratureOnPatches

  ! <<< SeungWhan: debugging >>>
  use, intrinsic :: iso_fortran_env, only : output_unit
  use ErrorHandler, only : writeAndFlush

  ! <<< Arguments >>>
  class(t_LighthillTensorComponent) :: this
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, nUnknowns, ierror
  SCALAR_TYPE, allocatable :: F(:,:)

  character(len=STRING_LENGTH) :: message
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

end function computeLighthillTensorComponent

subroutine computeLighthillTensorComponentAdjointForcing(this, simulationFlags, solverOptions,          &
     grid, state, patch)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use LighthillTensorComponent_mod, only : t_LighthillTensorComponent
  use SolverOptions_mod, only : t_SolverOptions
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use CNSHelper

  implicit none

  ! <<< Arguments >>>
  class(t_LighthillTensorComponent) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  class(t_CostTargetPatch) :: patch

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, nDimensions, nUnknowns, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: temp1(:,:,:), temp2(:,:,:), thirdPartialViscousJacobian(:,:),  &
       localFluxJacobian1(:,:), localFluxJacobian2(:,:), localConservedVariables(:),         &
       localVelocity(:), unitVector1(:), unitVector2(:),                                     &
       localStressTensor(:), localHeatFlux(:), localAdjointDiffusion(:,:)
  SCALAR_TYPE, allocatable :: F(:)
  real(wp) :: timeRampFactor

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))
  nUnknowns = nDimensions + 2

  allocate(localFluxJacobian1(nUnknowns, nUnknowns))
  allocate(localConservedVariables(nUnknowns))
  allocate(localVelocity(nDimensions))
  allocate(unitVector2(nDimensions))
  allocate(F(nUnknowns))
  unitVector2 = 0.0_wp
  unitVector2(this%secondComponent) = 1.0_wp

  timeRampFactor = 1.0_wp
  if (this%useTimeWindow)                                                                    &
      timeRampFactor = timeRampFactor                                                        &
          * exp( -(state%time-this%timeWindowCenter)**2/2.0_wp/this%timeWindowWidth/this%timeWindowWidth )
  ! ideal_mean_pressure = 1.0_wp/solverOptions%ratioOfSpecificHeats

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

           localConservedVariables = state%conservedVariables(gridIndex,:)
           localVelocity = state%velocity(gridIndex,:)

           select case (nDimensions)
           case (1)
              call computeJacobianOfInviscidFlux1D(localConservedVariables,                             &
                   unitVector2, solverOptions%ratioOfSpecificHeats,                                     &
                   localFluxJacobian1, specificVolume = state%specificVolume(gridIndex,1),              &
                   velocity = localVelocity, temperature = state%temperature(gridIndex,1))
           case (2)
              call computeJacobianOfInviscidFlux2D(localConservedVariables,                             &
                   unitVector2, solverOptions%ratioOfSpecificHeats,                                     &
                   localFluxJacobian1, specificVolume = state%specificVolume(gridIndex,1),              &
                   velocity = localVelocity, temperature = state%temperature(gridIndex,1))
           case (3)
              call computeJacobianOfInviscidFlux3D(localConservedVariables,                             &
                   unitVector2, solverOptions%ratioOfSpecificHeats,                                     &
                   localFluxJacobian1, specificVolume = state%specificVolume(gridIndex,1),              &
                   velocity = localVelocity, temperature = state%temperature(gridIndex,1))
           end select !... nDimensions

           ! For density identity tensor jacobian.
           localFluxJacobian1(this%secondComponent,1) = localFluxJacobian1(this%secondComponent,1)      &
                                                            - 1.0_wp

           F = - localFluxJacobian1(this%firstComponent+1,:) * timeRampFactor

           patch%adjointForcing(patchIndex,:) = F

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
     end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

  SAFE_DEALLOCATE(F)

end subroutine computeLighthillTensorComponentAdjointForcing

function isLighthillTensorComponentPatchValid(this, patchDescriptor, gridSize, normalDirection,         &
     extent, simulationFlags, message) result(isPatchValid)

  ! <<< Derived types >>>
  use LighthillTensorComponent_mod, only : t_LighthillTensorComponent
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_LighthillTensorComponent) :: this
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

end function isLighthillTensorComponentPatchValid
