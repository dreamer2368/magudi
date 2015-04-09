#include "config.h"

program isothermal_wall_SAT

  use MPI

  use, intrinsic :: iso_fortran_env, only : error_unit

  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random

  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use IsothermalWall_mod, only : t_IsothermalWall
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  use Region_enum, only : FORWARD, ADJOINT

  use CNSHelper, only : computeDependentVariables

  implicit none

  integer, parameter :: wp = SCALAR_KIND, n = 1
  real(SCALAR_KIND), parameter :: tolerance = sqrt(epsilon(0.0_wp))
  logical :: success
  integer :: i, j, k, direction, nDimensions, gridSize(3, 1),                                &
       extent(6), errorCode, numProcs, ierror
  character(len = STRING_LENGTH) :: message
  real(SCALAR_KIND), allocatable :: F(:,:)

  type(t_SimulationFlags) :: simulationFlags
  type(t_SolverOptions) :: solverOptions
  type(t_Grid) :: grid
  type(t_State) :: state
  type(t_PatchDescriptor) :: patchDescriptor
  type(t_IsothermalWall) :: patch

  interface

     subroutine applyForwardBoundaryConditions(patch, grid, state, ratioOfSpecificHeats)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use IsothermalWall_mod, only : t_IsothermalWall

       ! <<< Arguments >>>
       type(t_IsothermalWall), intent(in) :: patch
       type(t_Grid), intent(in) :: grid
       type(t_State) :: state
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats

     end subroutine applyForwardBoundaryConditions

     subroutine applyAdjointBoundaryConditions(patch, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use IsothermalWall_mod, only : t_IsothermalWall

       ! <<< Arguments >>>
       type(t_IsothermalWall), intent(in) :: patch
       type(t_Grid), intent(in) :: grid
       type(t_State) :: state

     end subroutine applyAdjointBoundaryConditions

     subroutine addSurfaceIntegralContribution(patch, grid, state, solverOptions)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use IsothermalWall_mod, only : t_IsothermalWall

       ! <<< Arguments >>>
       type(t_IsothermalWall), intent(in) :: patch
       type(t_Grid), intent(in) :: grid
       type(t_State) :: state
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine addSurfaceIntegralContribution

  end interface

  ! Initialize MPI.
  call MPI_Init(ierror)

  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)
  if (numProcs /= 1) then
     call MPI_Finalize(ierror)
     stop -1
  end if

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  call simulationFlags%initialize()
  simulationFlags%predictionOnly = .false.
  simulationFlags%viscosityOn = .true.

  do nDimensions = 2, 2!1, 3

     do k = 1, n

        direction = random(1, nDimensions)

        gridSize(:,:) = 1
        do i = 1, nDimensions
           gridSize(i,1) = random(1, 60)
        end do

        simulationFlags%isDomainCurvilinear = (random(0, 2) == 0)

        call solverOptions%initialize(nDimensions, simulationFlags)
        solverOptions%reynoldsNumberInverse = random(1e-6_wp, 1e-2_wp)
        solverOptions%prandtlNumberInverse = random(0.01_wp, 1e2_wp)
        solverOptions%bulkViscosityRatio = random(0.0_wp, 10.0_wp)
        solverOptions%powerLawExponent = random(0.0_wp, 2.0_wp)

        call grid%setup(1, gridSize(1:nDimensions,1), MPI_COMM_WORLD,                        &
             simulationFlags = simulationFlags)
        call grid%setupSpatialDiscretization(simulationFlags)

        call state%setup(grid, simulationFlags, solverOptions)

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

        patchDescriptor = t_PatchDescriptor("testPatch",                                     &
             "SAT_ISOTHERMAL_WALL", 1, direction,                                            &
             extent(1), extent(2), extent(3), extent(4), extent(5), extent(6))

        call patchDescriptor%validate(gridSize, simulationFlags,                             &
             solverOptions, errorCode, message)
        success = success .and. (errorCode == 0)
        if (.not. success) then
           write(error_unit, '(A)') trim(message)
           exit
        end if

        call patch%setup(1, MPI_COMM_WORLD, patchDescriptor,                                 &
             grid, simulationFlags, solverOptions)

        ! Randomize grid metrics.
        grid%jacobian = 1.0_wp
        allocate(F(grid%nGridPoints, nDimensions ** 2))
        call random_number(F)
        grid%metrics = F
        SAFE_DEALLOCATE(F)

        ! Randomize conserved variables.
        do i = 1, grid%nGridPoints
           state%conservedVariables(i,1) = random(0.01_wp, 10.0_wp)
           do j = 1, nDimensions
              state%conservedVariables(i,j+1) =                                              &
                   state%conservedVariables(i,1) * random(-10.0_wp, 10.0_wp)
           end do
           state%conservedVariables(i, nDimensions + 2) = state%conservedVariables(i,1) *    &
                random(0.01_wp, 10.0_wp) / solverOptions%ratioOfSpecificHeats +              &
                0.5_wp / state%conservedVariables(i,1) *                                     &
                sum(state%conservedVariables(i,2:nDimensions+1) ** 2)
        end do

        assert(all(state%conservedVariables(:,1) > 0.0_wp))

        ! Randomize wall temperature on patch.
        do i = 1, patch%nPatchPoints
           patch%temperature(i) = random(0.01_wp, 100.0_wp)
        end do

        ! Apply forward boundary conditions.
        call applyForwardBoundaryConditions(patch, grid, state,                              &
             solverOptions%ratioOfSpecificHeats)
        call computeDependentVariables(nDimensions, state%conservedVariables,                &
          solverOptions%ratioOfSpecificHeats, state%specificVolume(:,1), state%velocity,     &
          state%pressure(:,1), state%temperature(:,1))
        assert(all(state%specificVolume(:,1) > 0.0_wp))
        assert(all(state%temperature(:,1) > 0.0_wp))

        ! Check forward SAT consistency.
        state%rightHandSide = 0.0_wp
        call patch%updateRhs(FORWARD, simulationFlags, solverOptions, grid, state)

        success = success .and. all(abs(state%rightHandSide) < tolerance)
        if (.not. success) exit

        ! Randomize adjoint variables.
        allocate(F(grid%nGridPoints, solverOptions%nUnknowns))
        call random_number(F)
        state%adjointVariables = F
        SAFE_DEALLOCATE(F)

        ! Check adjoint SAT consistency.
        call applyAdjointBoundaryConditions(patch, grid, state)
        state%rightHandSide = 0.0_wp
        call patch%updateRhs(ADJOINT, simulationFlags, solverOptions, grid, state)
        call addSurfaceIntegralContribution(patch, grid, state, solverOptions)

        success = success .and. any(abs(state%rightHandSide) < tolerance)
        if (.not. success) exit

        call patch%cleanup()
        call state%cleanup()
        call grid%cleanup()

     end do !... k = 1, n

     if (.not. success) exit

  end do !... nDimensions = 1, 3

  call cleanupErrorHandler()

  call MPI_Finalize(ierror)

  if (.not. success) stop -1
  stop 0

end program isothermal_wall_SAT

subroutine applyForwardBoundaryConditions(patch, grid, state, ratioOfSpecificHeats)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use IsothermalWall_mod, only : t_IsothermalWall

  implicit none

  ! <<< Arguments >>>
  type(t_IsothermalWall), intent(in) :: patch
  type(t_Grid), intent(in) :: grid
  type(t_State) :: state
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, gridIndex, patchIndex, direction, nDimensions

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = patch%normalDirection
  assert(abs(direction) >= 1 .and. abs(direction) <= nDimensions)

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

           state%conservedVariables(gridIndex,2:nDimensions+1) = 0.0_wp

           state%conservedVariables(gridIndex,nDimensions+2) =                               &
                state%conservedVariables(gridIndex,1) * patch%temperature(patchIndex) /      &
                ratioOfSpecificHeats

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
     end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

end subroutine applyForwardBoundaryConditions

subroutine applyAdjointBoundaryConditions(patch, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use IsothermalWall_mod, only : t_IsothermalWall

  implicit none

  ! <<< Arguments >>>
  type(t_IsothermalWall), intent(in) :: patch
  type(t_Grid), intent(in) :: grid
  type(t_State) :: state

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

           localVelocity = state%adjointVariables(gridIndex,2:nDimensions+1)

           state%adjointVariables(gridIndex,2:nDimensions+1) =                               &
                localVelocity - unitNormal * dot_product(localVelocity, unitNormal)

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
     end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

  SAFE_DEALLOCATE(unitNormal)
  SAFE_DEALLOCATE(localVelocity)

end subroutine applyAdjointBoundaryConditions

subroutine addSurfaceIntegralContribution(patch, grid, state, solverOptions)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use IsothermalWall_mod, only : t_IsothermalWall

  ! <<< Internal modules >>>
  use CNSHelper

  implicit none

  ! <<< Arguments >>>
  type(t_IsothermalWall), intent(in) :: patch
  type(t_Grid), intent(in) :: grid
  type(t_State) :: state
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, gridIndex, patchIndex, direction, nDimensions, nUnknowns
  SCALAR_TYPE, allocatable :: localConservedVariables(:), localMetricsAlongDirection(:),     &
       localFluxJacobian(:,:)

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = patch%normalDirection
  assert(abs(direction) >= 1 .and. abs(direction) <= nDimensions)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2)

  allocate(localConservedVariables(nUnknowns))
  allocate(localMetricsAlongDirection(nDimensions))
  allocate(localFluxJacobian(nUnknowns, nUnknowns))

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
           localMetricsAlongDirection = grid%metrics(gridIndex, 1 + nDimensions *            &
                (abs(direction) - 1) : nDimensions * abs(direction))

           select case (nDimensions)
           case (1)
              call computeJacobianOfInviscidFlux1D(localConservedVariables,                  &
                   localMetricsAlongDirection, solverOptions%ratioOfSpecificHeats,           &
                   localFluxJacobian)
           case (2)
              call computeJacobianOfInviscidFlux2D(localConservedVariables,                  &
                   localMetricsAlongDirection, solverOptions%ratioOfSpecificHeats,           &
                   localFluxJacobian)
           case (3)
              call computeJacobianOfInviscidFlux3D(localConservedVariables,                  &
                   localMetricsAlongDirection, solverOptions%ratioOfSpecificHeats,           &
                   localFluxJacobian)
           end select !... nDimensions

           state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -             &
                sign(1.0_wp / grid%firstDerivative(abs(direction))%normBoundary(1),          &
                real(direction, wp)) * matmul(transpose(localFluxJacobian),                  &
                state%adjointVariables(gridIndex,:))

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
     end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

  SAFE_DEALLOCATE(localFluxJacobian)
  SAFE_DEALLOCATE(localMetricsAlongDirection)
  SAFE_DEALLOCATE(localConservedVariables)

end subroutine addSurfaceIntegralContribution
