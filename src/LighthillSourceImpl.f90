#include "config.h"

subroutine setupLighthillSource(this, region)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use LighthillSource_mod, only : t_LighthillSource

  ! <<< Enumerations >>>
  use State_enum, only : QOI_DUMMY_FUNCTION

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables
  use InputHelper, only : getOption, getRequiredOption


  implicit none

  ! <<< Arguments >>>
  class(t_LighthillSource) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, nDimensions, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix, message
  SCALAR_TYPE, allocatable :: temp1(:,:), temp2(:,:)

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

  if( region%simulationFlags%enableAdjoint ) then
    allocate(this%data_(size(region%globalGridSizes, 2)))
    do i = 1, size(this%data_)
       do j = 1, size(region%grids)
          if (region%grids(j)%index == i) then
             allocate(this%data_(i)%adjointTensor(region%grids(j)%nGridPoints,nDimensions,nDimensions))

             allocate(temp1(region%grids(j)%nGridPoints,1))
             allocate(temp2(region%grids(j)%nGridPoints,nDimensions))

             temp2 = 0.0_wp

             ! NOTE: adjoint derivative is transposed derivative PRE-MULTIPLIED by norm.
             do k = 1, nDimensions
               temp1 = region%grids(j)%targetMollifier
               call region%grids(j)%adjointFirstDerivative(k)%apply(temp1, region%grids(j)%localSize)
               do l = 1, nDimensions
                 temp2(:,l) = temp2(:,l) + temp1(:,1) * region%grids(j)%metrics(:,l+nDimensions*(k-1))          &
                                                      * region%grids(j)%jacobian(:,1)
               end do
             end do
             do k = 1, nDimensions
               this%data_(i)%adjointTensor(:,:,k) = temp2
               call region%grids(j)%adjointFirstDerivative(k)%apply(this%data_(i)%adjointTensor(:,:,k),         &
                                                                    region%grids(j)%localSize)
             end do

             SAFE_DEALLOCATE(temp1)
             SAFE_DEALLOCATE(temp2)
             exit
          end if
       end do
       call MPI_Barrier(region%comm, ierror)
    end do
  end if

end subroutine setupLighthillSource

subroutine cleanupLighthillSource(this)

  ! <<< Derived types >>>
  use LighthillSource_mod, only : t_LighthillSource

  implicit none

  ! <<< Arguments >>>
  class(t_LighthillSource) :: this

  ! <<< Local variables >>>
  integer :: i

  call this%cleanupBase()

end subroutine cleanupLighthillSource

subroutine computeLighthillSourceSpatialDistribution(this, grid, state, F)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use LighthillSource_mod, only : t_LighthillSource

  ! <<< Internal modules >>>
  use CNSHelper, only : computeCartesianInviscidFluxes, computeCartesianViscousFluxes, transformFluxes

  ! <<< Arguments >>>
  class(t_LighthillSource) :: this
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  SCALAR_TYPE, intent(out) :: F(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, nUnknowns, ierror
  logical :: computeViscousFlux = .false.
  SCALAR_TYPE, allocatable :: fluxes1(:,:,:), fluxes2(:,:,:), F1(:,:,:), F2(:,:,:)

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
  allocate(F1(grid%nGridPoints, 1, nDimensions))
  allocate(F2(grid%nGridPoints, 1, nDimensions))

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

  ! Transform fluxes from Cartesian to contravariant form: `fluxes1` has the Cartesian form of
  ! total fluxes... upon return, `fluxes2` has the contravariant form.
  call transformFluxes(nDimensions, fluxes1, grid%metrics, fluxes2, grid%isCurvilinear)

  ! Take derivatives of fluxes.
  do k = 1, nDimensions
      call grid%firstDerivative(k)%apply(fluxes2(:,2:nDimensions+1,k),         &
                                                    grid%localSize)
  end do
  F1(:,1,:) = sum(fluxes2(:,2:nDimensions+1,:), dim = 3)
  do k = 1, nDimensions
     F1(:,1,k) = F1(:,1,k) * grid%jacobian(:,1)
  end do

  ! Transform fluxes from Cartesian to contravariant form: `F1` has the Cartesian form of
  ! total fluxes... upon return, `F2` has the contravariant form.
  call transformFluxes(nDimensions, F1, grid%metrics, F2, grid%isCurvilinear)

  ! Take derivatives of fluxes.
  do j = 1, nDimensions
      call grid%firstDerivative(j)%apply(F2(:,:,j), grid%localSize)
  end do

  F(:,1) = sum(F2(:,1,:), dim = 2) * grid%jacobian(:,1)
  SAFE_DEALLOCATE(F1)
  SAFE_DEALLOCATE(F2)
  SAFE_DEALLOCATE(fluxes1)
  SAFE_DEALLOCATE(fluxes2)

end subroutine

function computeLighthillSource(this, region) result(instantaneousFunctional)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use LighthillSource_mod, only : t_LighthillSource

  ! <<< Internal modules >>>
  use Patch_factory, only : computeQuadratureOnPatches

  ! <<< SeungWhan: debugging >>>
  use, intrinsic :: iso_fortran_env, only : output_unit
  use ErrorHandler, only : writeAndFlush

  ! <<< Arguments >>>
  class(t_LighthillSource) :: this
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

end function computeLighthillSource

subroutine computeLighthillSourceAdjointForcing(this, simulationFlags, solverOptions,          &
     grid, state, patch)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use LighthillSource_mod, only : t_LighthillSource
  use SolverOptions_mod, only : t_SolverOptions
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use CNSHelper!, only : transformFluxes

  implicit none

  ! <<< Arguments >>>
  class(t_LighthillSource) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  class(t_CostTargetPatch) :: patch

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, nDimensions, nUnknowns, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: localFluxJacobian1(:,:), localConservedVariables(:),         &
                              localVelocity(:), localMetricsAlongDirection1(:)
  SCALAR_TYPE, allocatable :: F(:)
  SCALAR_TYPE, allocatable :: adjointTensor(:,:,:), jacobian(:)
  real(wp) :: timeRampFactor

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))
  nUnknowns = nDimensions + 2

  allocate(adjointTensor(patch%nPatchPoints,nDimensions,nDimensions))
  allocate(jacobian(patch%nPatchPoints))
  i = grid%index
  call patch%collect(this%data_(i)%adjointTensor, adjointTensor)
  call patch%collect(grid%jacobian(:,1), jacobian)

  allocate(localFluxJacobian1(nUnknowns, nUnknowns))
  allocate(localConservedVariables(nUnknowns))
  allocate(localVelocity(nDimensions))
  allocate(localMetricsAlongDirection1(nDimensions))

  allocate(F(nUnknowns))

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

           F = 0.0_wp

           localConservedVariables = state%conservedVariables(gridIndex,:)
           localVelocity = state%velocity(gridIndex,:)

           do l = 1, nDimensions
             localMetricsAlongDirection1 = grid%metrics(gridIndex,1+nDimensions*(l-1):nDimensions*l)

             select case (nDimensions)
             case (1)
                call computeJacobianOfInviscidFlux1D(localConservedVariables,                     &
                     localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                     localFluxJacobian1, specificVolume = state%specificVolume(gridIndex,1),              &
                     velocity = localVelocity, temperature = state%temperature(gridIndex,1))
             case (2)
                call computeJacobianOfInviscidFlux2D(localConservedVariables,                     &
                     localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                     localFluxJacobian1, specificVolume = state%specificVolume(gridIndex,1),              &
                     velocity = localVelocity, temperature = state%temperature(gridIndex,1))
             case (3)
                call computeJacobianOfInviscidFlux3D(localConservedVariables,                     &
                     localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                     localFluxJacobian1, specificVolume = state%specificVolume(gridIndex,1),              &
                     velocity = localVelocity, temperature = state%temperature(gridIndex,1))
             end select !... nDimensions

             localFluxJacobian1(2:nDimensions+1,1) = localFluxJacobian1(2:nDimensions+1,1)                  &
                                                       - localMetricsAlongDirection1
             F = F + matmul( transpose(localFluxJacobian1(2:nDimensions+1,:)),                              &
                             adjointTensor(patchIndex,:,l) )

           end do

           patch%adjointForcing(patchIndex,:) = - F * timeRampFactor                                        &
                                                    * jacobian(patchIndex)

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
     end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

  SAFE_DEALLOCATE(adjointTensor)
  SAFE_DEALLOCATE(jacobian)
  SAFE_DEALLOCATE(F)

end subroutine computeLighthillSourceAdjointForcing

function isLighthillSourcePatchValid(this, patchDescriptor, gridSize, normalDirection,         &
     extent, simulationFlags, message) result(isPatchValid)

  ! <<< Derived types >>>
  use LighthillSource_mod, only : t_LighthillSource
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_LighthillSource) :: this
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

end function isLighthillSourcePatchValid
