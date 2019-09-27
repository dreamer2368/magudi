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
  integer :: i, j, k, nDimensions, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix, message
  SCALAR_TYPE, allocatable :: temp(:,:)

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

  this%firstDirection  = 0.0_wp
  this%secondDirection = 0.0_wp

  this%firstDirection(1) = getOption("lighthill_tensor_direction1_x", 1.0_wp)
  if (nDimensions >= 2)                                                                      &
       this%firstDirection(2) = getOption("lighthill_tensor_direction1_y", 0.0_wp)
  if (nDimensions == 3)                                                                      &
       this%firstDirection(3) = getOption("lighthill_tensor_direction1_z", 0.0_wp)

  this%secondDirection(1) = getOption("lighthill_tensor_direction2_x", 1.0_wp)
  if (nDimensions >= 2)                                                                      &
       this%secondDirection(2) = getOption("lighthill_tensor_direction2_y", 0.0_wp)
  if (nDimensions == 3)                                                                      &
       this%secondDirection(3) = getOption("lighthill_tensor_direction2_z", 0.0_wp)

  if (sum(this%firstDirection ** 2) <= epsilon(0.0_wp) .or.                                  &
       sum(this%secondDirection ** 2) <= epsilon(0.0_wp)) then
     write(message, '(A)')                                                                   &
          "Unable to determine unit vectors for computing Lighthill tensor!"
     call gracefulExit(region%comm, message)
  end if

  this%firstDirection = this%firstDirection / sqrt(sum(this%firstDirection ** 2))
  this%secondDirection = this%secondDirection / sqrt(sum(this%secondDirection ** 2))

  this%firstComponent = getOption("lighthill_tensor_component1", 1)
  assert_key(this%firstComponent, (1,2,3))
  this%secondComponent = getOption("lighthill_tensor_component2", 1)
  assert_key(this%secondComponent, (1,2,3))

  if( region%simulationFlags%enableAdjoint ) then
    allocate(this%data_(size(region%globalGridSizes, 2)))
    do i = 1, size(this%data_)
       do j = 1, size(region%grids)
          if (region%grids(j)%index == i) then
             allocate(this%data_(i)%adjointVector(region%grids(j)%nGridPoints,nDimensions))

             allocate(temp(region%grids(j)%nGridPoints,1))

             ! NOTE: adjoint derivative is transposed derivative PRE-MULTIPLIED by norm.
             do k = 1, nDimensions
               temp = region%grids(j)%targetMollifier * region%grids(j)%jacobian
               call region%grids(j)%adjointFirstDerivative(k)%apply(temp, region%grids(j)%localSize)
               this%data_(i)%adjointVector(:,k) = temp(:,1)
             end do

             SAFE_DEALLOCATE(temp)
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
  ! do k = 1, nDimensions
  !   fluxes1(:,k+1,k) = fluxes1(:,k+1,k) - state%conservedVariables(:,1)
  ! end do

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

  ! do i = 1, grid%nGridPoints
  !   F(i,1) = dot_product( F1(i,1,:), this%firstDirection(1:nDimensions) )
  ! end do

  ! F(:,1) = fluxes2(:,this%firstComponent+1,this%secondComponent)
  F(:,1) = F1(:,1,this%firstComponent)

  ! Verified part.
  ! ! Transform fluxes from Cartesian to contravariant form: `F1` has the Cartesian form of
  ! ! total fluxes... upon return, `F2` has the contravariant form.
  ! call transformFluxes(nDimensions, F1, grid%metrics, F2, grid%isCurvilinear)
  !
  ! ! Take derivatives of fluxes.
  ! do j = 1, nDimensions
  !     call grid%firstDerivative(j)%apply(F2(:,:,j), grid%localSize)
  ! end do
  !
  ! F(:,1) = sum(F2(:,1,:), dim = 2) * grid%jacobian(:,1)
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
  SCALAR_TYPE, allocatable :: temp1(:,:,:), temp2(:,:,:), thirdPartialViscousJacobian(:,:),  &
       localFluxJacobian1(:,:), localFluxJacobian2(:,:), localConservedVariables(:),         &
       localVelocity(:), localMetricsAlongDirection1(:),                                     &
       localStressTensor(:), localHeatFlux(:), localAdjointDiffusion(:,:)
  SCALAR_TYPE, allocatable :: divTensor(:,:,:), divTensor2(:,:), F(:)
  SCALAR_TYPE, allocatable :: adjointVector(:,:), jacobian(:)
  real(wp) :: temp, timeRampFactor

  ! <<< Internal interface >>>
  interface

     subroutine computeLighthillDivTensor(simulationFlags, solverOptions, grid, state, patch, secondComponent, divTensor)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags
       use CostTargetPatch_mod, only : t_CostTargetPatch

       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid) :: grid
       class(t_State) :: state
       class(t_CostTargetPatch) :: patch
       integer, intent(in) :: secondComponent
       SCALAR_TYPE, intent(out) :: divTensor(:,:,:)

     end subroutine computeLighthillDivTensor

  end interface

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))
  nUnknowns = nDimensions + 2

  allocate(adjointVector(patch%nPatchPoints,nDimensions))
  ! allocate(jacobian(patch%nPatchPoints))
  i = grid%index
  call patch%collect(this%data_(i)%adjointVector, adjointVector)
  ! call patch%collect(grid%jacobian(:,1), jacobian)

  allocate(localFluxJacobian1(nUnknowns, nUnknowns))
  allocate(localConservedVariables(nUnknowns))
  allocate(localVelocity(nDimensions))
  allocate(localMetricsAlongDirection1(nDimensions))

  ! allocate(divTensor(grid%nGridPoints,nUnknowns,nDimensions))
  allocate(divTensor(grid%nGridPoints,nDimensions,nUnknowns))
  allocate(divTensor2(nDimensions,nUnknowns))
  ! allocate(F(grid%nGridPoints,nUnknowns))
  allocate(F(nUnknowns))
  ! call computeLighthillDivTensor(simulationFlags,solverOptions,grid,state,patch,this%secondComponent,divTensor)

  ! ! Transform fluxes from Cartesian to contravariant form: `fluxes1` has the Cartesian form of
  ! ! total fluxes... upon return, `fluxes2` has the contravariant form.
  ! call transformFluxes(nDimensions, divTensor, grid%metrics, divTensor2, grid%isCurvilinear)
  !
  ! ! Take derivatives of fluxes.
  ! do k = 1, nDimensions
  !     call grid%firstDerivative(k)%apply(divTensor2(:,:,k), grid%localSize)
  ! end do
  ! F = sum(divTensor2, dim = 3)

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

           ! F = matmul( transpose(divTensor(gridIndex,:,:)), this%firstDirection(1:nDimensions) )
           ! F = divTensor(gridIndex,this%firstComponent,:)
           !
           ! patch%adjointForcing(patchIndex,:) = - F * timeRampFactor                                       &
           !                                          * grid%targetMollifier(gridIndex, 1)

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

             ! F = F + matmul( transpose(localFluxJacobian1(2:nDimensions+1,:)),                              &
             !                 this%firstDirection(1:nDimensions) )                                           &
             !           * this%secondDirection(l)
             ! F(1) = F(1) - dot_product( localMetricsAlongDirection1, this%firstDirection(1:nDimensions) )   &
             !                  * this%secondDirection(l)
             F = F + localFluxJacobian1(this%firstComponent+1,:) * adjointVector(patchIndex,l)
             ! F(1) = F(1) - localMetricsAlongDirection1(this%firstComponent)

           end do

           patch%adjointForcing(patchIndex,:) = - F * timeRampFactor!                                        &
                                                    ! * jacobian(patchIndex)

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
     end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

  SAFE_DEALLOCATE(adjointVector)
  ! SAFE_DEALLOCATE(jacobian)
  SAFE_DEALLOCATE(divTensor)
  SAFE_DEALLOCATE(divTensor2)
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

subroutine computeLighthillDivTensor(simulationFlags, solverOptions, grid, state, patch, secondComponent, divTensor)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use CostTargetPatch_mod, only : t_CostTargetPatch

  ! <<< Enumerations >>>
  use Region_enum, only : ADJOINT

  ! <<< Internal modules >>>
  use CNSHelper

  implicit none

  ! <<< Arguments >>>
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid) :: grid
  class(t_State) :: state
  class(t_CostTargetPatch) :: patch
  integer, intent(in) :: secondComponent
  SCALAR_TYPE, intent(out) :: divTensor(:,:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, m, gridIndex, patchIndex, nDimensions, nUnknowns
  SCALAR_TYPE, allocatable :: temp1(:,:,:), temp2(:,:,:), thirdPartialViscousJacobian(:,:),  &
       localFluxJacobian1(:,:), localFluxJacobian2(:,:), localConservedVariables(:),         &
       localVelocity(:), localMetricsAlongDirection1(:), localMetricsAlongDirection2(:),     &
       localStressTensor(:), localHeatFlux(:), localAdjointDiffusion(:,:)

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns >= nDimensions + 2)

  assert(size(secondDirection) == nDimensions)
  assert(size(divTensor,1) == grid%nGridPoints)
  assert(size(divTensor,2) == nUnknowns)
  assert(size(divTensor,3) == nDimensions)
  divTensor = 0.0_wp

  ! allocate(temp1(grid%nGridPoints, nUnknowns, nDimensions))
  allocate(temp1(grid%nGridPoints, nDimensions, nUnknowns))

  allocate(localFluxJacobian1(nUnknowns, nUnknowns))
  allocate(localConservedVariables(nUnknowns))
  allocate(localVelocity(nDimensions))
  allocate(localMetricsAlongDirection1(nDimensions))

  ! do l = 1, nDimensions
  l = secondComponent

    temp1 = 0.0_wp

    do gridIndex = 1, grid%nGridPoints
    ! do k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)
    !    do j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
    !       do i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
    !          gridIndex = i - patch%gridOffset(1) + patch%gridLocalSize(1) *                    &
    !               (j - 1 - patch%gridOffset(2) + patch%gridLocalSize(2) *                      &
    !               (k - 1 - patch%gridOffset(3)))
             if (grid%iblank(gridIndex) == 0) cycle

             localConservedVariables = state%conservedVariables(gridIndex,:)
             localVelocity = state%velocity(gridIndex,:)

            localMetricsAlongDirection1 = grid%metrics(gridIndex,1+nDimensions*(l-1):nDimensions*l)

            select case (nDimensions)
            case (1)
               call computeJacobianOfInviscidFlux1D(localConservedVariables,                             &
                    localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,                     &
                    localFluxJacobian1, specificVolume = state%specificVolume(gridIndex,1),              &
                    velocity = localVelocity, temperature = state%temperature(gridIndex,1))
            case (2)
               call computeJacobianOfInviscidFlux2D(localConservedVariables,                             &
                    localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,                     &
                    localFluxJacobian1, specificVolume = state%specificVolume(gridIndex,1),              &
                    velocity = localVelocity, temperature = state%temperature(gridIndex,1))
            case (3)
               call computeJacobianOfInviscidFlux3D(localConservedVariables,                             &
                    localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,                     &
                    localFluxJacobian1, specificVolume = state%specificVolume(gridIndex,1),              &
                    velocity = localVelocity, temperature = state%temperature(gridIndex,1))
            end select !... nDimensions

            ! For density identity tensor jacobian.
            localFluxJacobian1(2:nDimensions+1,1) = localFluxJacobian1(2:nDimensions+1,1)         &
                                                             - localMetricsAlongDirection1

            temp1(gridIndex,:,:) = localFluxJacobian1(2:nDimensions+1,:)

   !       end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
   !    end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
   ! end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)
   end do

   do m = 1, nUnknowns
     print *, 'adjoint, ',m,' component, before derivative:'
     print *, temp1(:,1,m)
     call grid%firstDerivative(l)%apply(temp1(:,:,m), grid%localSize)
     print *, 'adjoint, ',m,' component, after derivative:'
     print *, temp1(:,1,m)
   end do

   divTensor = divTensor + temp1

 ! end do !... l = 1, nDimensions

 do j = 1, size(divTensor,3)
   do i = 1, size(divTensor,2)
     divTensor(:,i,j) = divTensor(:,i,j) * grid%jacobian(:,1)
   end do
 end do

  SAFE_DEALLOCATE(localConservedVariables)
  SAFE_DEALLOCATE(localFluxJacobian1)
  SAFE_DEALLOCATE(localHeatFlux)
  SAFE_DEALLOCATE(localStressTensor)
  SAFE_DEALLOCATE(localFluxJacobian2)

  SAFE_DEALLOCATE(localMetricsAlongDirection1)
  SAFE_DEALLOCATE(localVelocity)
  SAFE_DEALLOCATE(temp1)

end subroutine computeLighthillDivTensor
