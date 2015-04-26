#include "config.h"

program drag_adjoint_forcing

  use MPI
  use, intrinsic :: iso_fortran_env, only : real64

  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random

  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use Functional_mod, only : t_Functional
  use Functional_factory, only : t_FunctionalFactory
  use PressureDrag_mod, only : t_PressureDrag

  use Region_enum, only : ADJOINT

  use CNSHelper, only : computeDependentVariables, computeJacobianOfInviscidFlux2D

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  real(real64), parameter :: testDuration = 3.0_real64
  character(len = STRING_LENGTH), parameter :: discretizationTypes(4) =                      &
       (/ "SBP 1-2", "SBP 2-4", "SBP 3-6", "SBP 4-8" /)
  real(SCALAR_KIND), parameter :: tolerance = sqrt(epsilon(0.0_wp))
  integer :: i, direction, nDimensions, gridSize(3,1), extent(6), numProcs, ierror
  logical :: success, success_
  real(real64) :: startTime
  character(len = STRING_LENGTH) :: discretizationType
  logical :: isDomainCurvilinear
  real(SCALAR_KIND) :: dragDirection(3)

  type(t_Region) :: region(2)
  type(t_FunctionalFactory) :: functionalFactory(2)
  class(t_Functional), pointer :: dragCoefficient => null()

  interface

     subroutine setupTestRegion(region, discretizationType, gridSize, patchExtent,           &
          normalDirection, useContinuousAdjoint, isDomainCurvilinear, success)

       use Region_mod, only : t_Region
       use SolverOptions_mod, only : t_SolverOptions

       class(t_Region) :: region
       integer, intent(in) :: gridSize(:,:), patchExtent(6), normalDirection
       character(len = *), intent(in) :: discretizationType
       logical, intent(in) :: useContinuousAdjoint, isDomainCurvilinear
       logical, intent(out) :: success

     end subroutine setupTestRegion

     subroutine randomizeTestRegionData(grid, state, patchFactories,                         &
          ratioOfSpecificHeats, dragDirection, success)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use Patch_factory, only : t_PatchFactory

       class(t_Grid) :: grid
       class(t_State) :: state
       type(t_PatchFactory), intent(in) :: patchFactories(2)
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats, dragDirection(3)
       logical, intent(out) :: success

     end subroutine randomizeTestRegionData

     subroutine zeroOutRhsOnOtherPatches(nDimensions, direction, gridSize, rightHandSide)

       integer, intent(in) :: nDimensions, direction, gridSize(3)
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine zeroOutRhsOnOtherPatches

  end interface

  call MPI_Init(ierror)

  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)
  if (numProcs /= 1) then
     call MPI_Finalize(ierror)
     stop -1
  end if

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  startTime = MPI_Wtime()

  do

     do nDimensions = 1, 3

        direction = random(1, nDimensions)

        gridSize(:,:) = 1
        do i = 1, nDimensions
           gridSize(i,1) = random(20, 60)
        end do

        extent(1::2) = 1
        extent(2::2) = gridSize(:,1)
        if (random(0, 1) == 0) then
           extent(1+2*(direction-1)) = 1
           extent(2+2*(direction-1)) = 1
        else
           extent(1+2*(direction-1)) = gridSize(direction, 1)
           extent(2+2*(direction-1)) = gridSize(direction, 1)
           direction = -direction
        end if

        discretizationType = discretizationTypes(random(1, size(discretizationTypes)))
        isDomainCurvilinear = (random(0, 1) == 0)

        call setupTestRegion(region(1), trim(discretizationType), gridSize(1:nDimensions,:), &
             extent, direction, .true., isDomainCurvilinear, success_)
        success = success .and. success_
        if (.not. success) exit

        call setupTestRegion(region(2), trim(discretizationType), gridSize(1:nDimensions,:), &
             extent, direction, .false., isDomainCurvilinear, success_)
        success = success .and. success_
        if (.not. success) exit

        dragDirection = 0.0_wp
        do i = 1, nDimensions
           dragDirection(i) = random(0.0_wp, 1.0_wp)
        end do

        ! Generate random data that satisfies the boundary conditions.
        call randomizeTestRegionData(region(1)%grids(1), region(1)%states(1),                &
             region(1)%patchFactories, region(1)%solverOptions%ratioOfSpecificHeats,         &
             dragDirection, success_)
        success = success .and. success_
        if (.not. success) exit

        ! Copy data from `region(1)` to `region(2)`
        region(2)%grids(1)%jacobian = region(1)%grids(1)%jacobian
        region(2)%grids(1)%metrics = region(1)%grids(1)%metrics
        region(2)%states(1)%conservedVariables = region(1)%states(1)%conservedVariables
        region(2)%states(1)%adjointVariables = region(1)%states(1)%adjointVariables
        call computeDependentVariables(nDimensions, region(2)%states(1)%conservedVariables,  &
             region(2)%solverOptions%ratioOfSpecificHeats,                                   &
             region(2)%states(1)%specificVolume(:,1), region(2)%states(1)%velocity,          &
             region(2)%states(1)%pressure(:,1), region(2)%states(1)%temperature(:,1))

        do i = 1, 2

           ! Create a functional for computing drag coefficient.
           call functionalFactory(i)%connect(dragCoefficient, "PRESSURE_DRAG", .true.)
           success = success .and. associated(dragCoefficient)
           if (.not. success) exit

           ! Setup the drag coefficient functional:

           select type (dragCoefficient)

           class is (t_PressureDrag)
              call dragCoefficient%setup(region(i))
              dragCoefficient%direction = dragDirection
              call dragCoefficient%updateAdjointForcing(region(i))

           class default
              success = .false.
              exit

           end select

           ! Compute adjoint RHS.
           call region(i)%computeRhs(ADJOINT)
           call zeroOutRhsOnOtherPatches(nDimensions, direction, gridSize(:,1),              &
                region(i)%states(1)%rightHandSide)

        end do !... i = 1, 2

        if (.not. success) exit

        success = success .and. (maxval(abs(region(1)%states(1)%rightHandSide -              &
             region(2)%states(1)%rightHandSide)) < tolerance)
        if (.not. success) exit

     end do !... nDimensions = 1, 3

     if (.not. success) exit
     if (MPI_Wtime() - startTime > testDuration) exit !... run for `testDuration` seconds

  end do

  do i = 1, 2
     call functionalFactory(i)%cleanup()
     call region(i)%cleanup()
  end do
  call cleanupErrorHandler()

  call MPI_Finalize(ierror)

  if (.not. success) stop -1
  stop 0

end program drag_adjoint_forcing

subroutine setupTestRegion(region, discretizationType, gridSize, patchExtent,                &
     normalDirection, useContinuousAdjoint, isDomainCurvilinear, success)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use SolverOptions_mod, only : t_SolverOptions
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SimulationFlags_mod, only : t_SimulationFlags
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use ImpenetrableWall_mod, only : t_ImpenetrableWall

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region
  integer, intent(in) :: gridSize(:,:), patchExtent(6), normalDirection
  character(len = *), intent(in) :: discretizationType
  logical, intent(in) :: useContinuousAdjoint, isDomainCurvilinear
  logical, intent(out) :: success

  ! <<< Local variables >>>
  integer :: nDimensions, errorCode
  type(t_SimulationFlags) :: simulationFlags
  type(t_SolverOptions) :: solverOptions
  type(t_PatchDescriptor) :: patchDescriptor
  character(len = STRING_LENGTH) :: str
  class(t_Patch), pointer :: patch => null()

  assert(size(gridSize, 2) == 1)

  nDimensions = size(gridSize, 1)
  assert_key(nDimensions, (1, 2, 3))

  success = .false.

  call simulationFlags%initialize()
  simulationFlags%predictionOnly = .false.
  simulationFlags%isDomainCurvilinear = isDomainCurvilinear
  simulationFlags%useContinuousAdjoint = useContinuousAdjoint

  call solverOptions%initialize(nDimensions, simulationFlags, MPI_COMM_WORLD)
  solverOptions%costFunctionalType = "PRESSURE_DRAG"
  solverOptions%discretizationType = trim(discretizationType)
  call region%setup(MPI_COMM_WORLD, gridSize, simulationFlags = simulationFlags,             &
       solverOptions = solverOptions, verbose = .false.)

  allocate(region%patchFactories(2))

  ! Setup wall patch:

  patchDescriptor = t_PatchDescriptor("impenetrableWall", "SAT_SLIP_WALL", 1,                &
       normalDirection, patchExtent(1), patchExtent(2), patchExtent(3),                      &
       patchExtent(4), patchExtent(5), patchExtent(6))

  call patchDescriptor%validate(gridSize, simulationFlags, solverOptions, errorCode, str)
  if (errorCode /= 0) return

  call region%patchFactories(1)%connect(patch, trim(patchDescriptor%patchType), .true.)
  if (.not. associated(patch)) return

  select type (patch)
  class is (t_ImpenetrableWall)
     call patch%setup(1, region%comm, patchDescriptor, region%grids(1),                      &
          simulationFlags, solverOptions)
  class default
     return
  end select

  patchDescriptor = t_PatchDescriptor("targetRegion", "COST_TARGET", 1,                      &
       normalDirection, patchExtent(1), patchExtent(2), patchExtent(3),                      &
       patchExtent(4), patchExtent(5), patchExtent(6))

  call patchDescriptor%validate(gridSize, simulationFlags, solverOptions, errorCode, str)
  if (errorCode /= 0) return

  call region%patchFactories(2)%connect(patch, trim(patchDescriptor%patchType), .true.)
  if (.not. associated(patch)) return

  select type (patch)
  class is (t_CostTargetPatch)
     call patch%setup(1, region%comm, patchDescriptor, region%grids(1),                      &
          simulationFlags, solverOptions)
  class default
     return
  end select

  success = .true.

end subroutine setupTestRegion

subroutine randomizeTestRegionData(grid, state, patchFactories, ratioOfSpecificHeats,        &
     dragDirection, success)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use Patch_mod, only : t_Patch
  use State_mod, only : t_State
  use Patch_factory, only : t_PatchFactory
  use ImpenetrableWall_mod, only : t_ImpenetrableWall

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables
  use RandomNumber, only : random

  implicit none

  ! <<< Arguments >>>
  class(t_Grid) :: grid
  class(t_State) :: state
  type(t_PatchFactory), intent(in) :: patchFactories(2)
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats, dragDirection(3)
  logical, intent(out) :: success

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions
  real(SCALAR_KIND), allocatable :: F(:,:)
  class(t_Patch), pointer :: patch => null()
  class(t_ImpenetrableWall), pointer :: impenetrableWall => null()

  interface

     subroutine applyForwardBoundaryConditions(patch, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use ImpenetrableWall_mod, only : t_ImpenetrableWall

       class(t_ImpenetrableWall), intent(in) :: patch
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine applyForwardBoundaryConditions

     subroutine applyAdjointBoundaryConditions(patch, grid, state, dragDirection)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use ImpenetrableWall_mod, only : t_ImpenetrableWall

       class(t_ImpenetrableWall), intent(in) :: patch
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state
       real(SCALAR_KIND), intent(in) :: dragDirection(3)

     end subroutine applyAdjointBoundaryConditions

  end interface

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  success = .false.

  ! Jacobian.
  do i = 1, grid%nGridPoints
     grid%jacobian(i,1) = random(1e-4_wp, 1e4_wp)
  end do

  ! Normalized metrics.
  allocate(F(grid%nGridPoints, nDimensions ** 2))
  call random_number(F)
  grid%metrics = F
  SAFE_DEALLOCATE(F)

  ! Conserved variables.
  do i = 1, grid%nGridPoints
     state%conservedVariables(i,1) = random(0.01_wp, 10.0_wp)
     do j = 1, nDimensions
        state%conservedVariables(i,j+1) =                                                    &
             state%conservedVariables(i,1) * random(-10.0_wp, 10.0_wp)
     end do
     state%conservedVariables(i, nDimensions + 2) = state%conservedVariables(i,1) *          &
          random(0.01_wp, 10.0_wp) / ratioOfSpecificHeats +                                  &
          0.5_wp / state%conservedVariables(i,1) *                                           &
          sum(state%conservedVariables(i,2:nDimensions+1) ** 2)
  end do

  assert(all(state%conservedVariables(:,1) > 0.0_wp))

  ! Adjoint variables.
  allocate(F(grid%nGridPoints, size(state%adjointVariables, 2)))
  call random_number(F)
  state%adjointVariables = F
  ! state%adjointVariables = 0.0_wp
  SAFE_DEALLOCATE(F)

  call patchFactories(1)%connect(patch)
  select type (patch)
  class is (t_ImpenetrableWall)
     impenetrableWall => patch
  class default
     return
  end select

  ! Apply forward boundary conditions.
  call applyForwardBoundaryConditions(impenetrableWall, grid, state)

  ! Compute dependent variables.
  call computeDependentVariables(nDimensions, state%conservedVariables,                      &
       ratioOfSpecificHeats, state%specificVolume(:,1), state%velocity,                      &
       state%pressure(:,1), state%temperature(:,1))

  assert(all(state%specificVolume(:,1) > 0.0_wp))
  assert(all(state%temperature(:,1) > 0.0_wp))

  ! Apply adjoint boundary conditions.
  call applyAdjointBoundaryConditions(impenetrableWall, grid, state, dragDirection)

  success = .true.

end subroutine randomizeTestRegionData

subroutine applyForwardBoundaryConditions(patch, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use ImpenetrableWall_mod, only : t_ImpenetrableWall

  implicit none

  ! <<< Arguments >>>
  class(t_ImpenetrableWall), intent(in) :: patch
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, gridIndex, patchIndex, direction, nDimensions
  SCALAR_TYPE, allocatable :: unitNormal(:), localVelocity(:)

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = patch%normalDirection
  assert(abs(direction) >= 1 .and. abs(direction) <= nDimensions)

  allocate(unitNormal(nDimensions))
  allocate(localVelocity(nDimensions))

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

           unitNormal = grid%metrics(gridIndex, 1 + nDimensions * (abs(direction) - 1) :     &
                nDimensions * abs(direction))
           unitNormal = unitNormal / sqrt(sum(unitNormal ** 2))

           localVelocity = state%conservedVariables(gridIndex,2:nDimensions+1) /             &
                state%conservedVariables(gridIndex,1)

           state%conservedVariables(gridIndex,2:nDimensions+1) =                             &
                localVelocity - unitNormal * dot_product(localVelocity, unitNormal)

           state%conservedVariables(gridIndex,nDimensions+2) =                               &
                state%conservedVariables(gridIndex,nDimensions+2) +                          &
                0.5_wp * state%conservedVariables(gridIndex,1) *                             &
                (sum(state%conservedVariables(gridIndex,2:nDimensions+1) ** 2 -              &
                localVelocity ** 2))

           state%conservedVariables(gridIndex,2:nDimensions+1) =                             &
                state%conservedVariables(gridIndex,1) *                                      &
                state%conservedVariables(gridIndex,2:nDimensions+1)

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
     end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

  SAFE_DEALLOCATE(unitNormal)
  SAFE_DEALLOCATE(localVelocity)

end subroutine applyForwardBoundaryConditions

subroutine applyAdjointBoundaryConditions(patch, grid, state, dragDirection)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use ImpenetrableWall_mod, only : t_ImpenetrableWall

  implicit none

  ! <<< Arguments >>>
  class(t_ImpenetrableWall), intent(in) :: patch
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state
  real(SCALAR_KIND), intent(in) :: dragDirection(3)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, gridIndex, patchIndex, direction, nDimensions
  SCALAR_TYPE, allocatable :: unitNormal(:)

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = patch%normalDirection
  assert(abs(direction) >= 1 .and. abs(direction) <= nDimensions)

  allocate(unitNormal(nDimensions))

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

           unitNormal = grid%metrics(gridIndex, 1 + nDimensions * (abs(direction) - 1) :     &
                nDimensions * abs(direction))
           unitNormal = unitNormal / sqrt(sum(unitNormal ** 2))

           state%adjointVariables(gridIndex,2:nDimensions+1) =                               &
                state%adjointVariables(gridIndex,2:nDimensions+1) - unitNormal *             &
                dot_product(state%adjointVariables(gridIndex,2:nDimensions+1) -              &
                dragDirection(1:nDimensions), unitNormal)

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
     end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

  SAFE_DEALLOCATE(unitNormal)

end subroutine applyAdjointBoundaryConditions

subroutine zeroOutRhsOnOtherPatches(nDimensions, direction, gridSize, rightHandSide)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: nDimensions, direction, gridSize(3)
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, extent(6)

  assert_key(nDimensions, (1, 2, 3))
  assert(size(rightHandSide, 1) > 0)
  assert(size(rightHandSide, 2) == nDimensions + 2)

  assert(direction /= 0 .and. abs(direction) <= nDimensions)
  assert(all(gridSize > 0))

#ifdef DEBUG
  if (nDimensions < 3) then
     assert(all(gridSize(nDimensions+1:3) == 1))
  end if
#endif

  do l = 1, nDimensions

     if (l == direction) cycle

     extent(1::2) = 1
     extent(2::2) = gridSize

     extent(1+2*(l-1)) = 1
     extent(2+2*(l-1)) = 1

     do i = extent(1), extent(2)
        do j = extent(3), extent(4)
           do k = extent(5), extent(6)
              rightHandSide(i+gridSize(1)*(j-1+gridSize(2)*(k-1)),:) = 0.0_wp
           end do
        end do
     end do

  end do

  do l = 1, nDimensions

     if (l == -direction) cycle

     extent(1::2) = 1
     extent(2::2) = gridSize

     extent(1+2*(l-1)) = gridSize(l)
     extent(2+2*(l-1)) = gridSize(l)

     do i = extent(1), extent(2)
        do j = extent(3), extent(4)
           do k = extent(5), extent(6)
              rightHandSide(i+gridSize(1)*(j-1+gridSize(2)*(k-1)),:) = 0.0_wp
           end do
        end do
     end do

  end do

end subroutine zeroOutRhsOnOtherPatches
