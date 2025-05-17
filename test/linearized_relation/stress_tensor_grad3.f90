#include "config.h"

program stress_tensor_grad3

  use MPI

  use CNSHelper
  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random

  implicit none

  logical :: success, success_, isPeriodic
  integer :: i, j, nDimensions, direction, ierror
  integer :: procRank
  character(len = STRING_LENGTH), parameter :: discretizationTypes(4) =                      &
       (/ "SBP 1-2", "SBP 2-4", "SBP 3-6", "SBP 4-8" /)

  interface

     real(SCALAR_KIND) function meanTrimmed(a)
       real(SCALAR_KIND), intent(inout) :: a(:)
     end function meanTrimmed

     subroutine sort(a)
       real(SCALAR_KIND), intent(inout) :: a(:)
     end subroutine sort

     subroutine testAdjointRelation(identifier, nDimensions, success, isPeriodic, direction, tolerance)

       character(len = *), intent(in) :: identifier
       integer, intent(in) :: nDimensions
       logical, intent(out) :: success

       logical, intent(in) :: isPeriodic
       integer, intent(in) :: direction
       real(SCALAR_KIND), intent(in), optional :: tolerance

     end subroutine testAdjointRelation

  end interface

  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  do nDimensions = 1, 3
    do j = 1, 4!... for each discretizationTypes
      success = .true.
      do i = 1, 10 !... test multiple times
        ! Didn't test periodic grid yet!!
        ! isPeriodic = .true.
        ! do direction = 1, nDimensions
        !   call testAdjointRelation(discretizationTypes(j), nDimensions,           &
        !                            success_, isPeriodic, direction)
        !   success = success .and. success_
        ! end do
        ! if( .not. success_) then
        !   if( procRank == 0 ) then
        !     print *, 'Failed, ', trim(discretizationTypes(j))
        !     print *, 'dimension: ', nDimensions
        !     print *, 'periodicity: ', isPeriodic
        !     print *, 'direction: ', direction
        !   end if
        !   exit
        ! end if

        isPeriodic = .false.
        do direction = 1, nDimensions
          call testAdjointRelation(discretizationTypes(j), nDimensions,           &
                                   success_, isPeriodic, direction)
          success = success .and. success_
        end do
        if( .not. success_) then
          if( procRank == 0 ) then
            print *, 'Failed, ', trim(discretizationTypes(j))
            print *, 'dimension: ', nDimensions
            print *, 'periodicity: ', isPeriodic
            print *, 'direction: ', direction
          end if
          exit
        end if
      end do
      if( procRank == 0 .and. success ) then
        print *, 'Success, ', trim(discretizationTypes(j))
        print *, 'dimension: ', nDimensions
      end if
    end do
  end do

  call cleanupErrorHandler()

  call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierror)
  call MPI_Finalize(ierror)
  if (.not. success) stop -1
  stop 0

end program stress_tensor_grad3

subroutine sort(a)

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(inout) :: a(:)

  ! <<< Local variables >>>
  integer :: i, j
  real(SCALAR_KIND) :: temp

  do i = 2, size(a)
     j = i - 1
     temp = a(i)
     do while (a(j) > temp)
        a(j+1) = a(j)
        j = j - 1
        if (j < 1) exit
     end do
     a(j+1) = temp
  end do

end subroutine sort

real(SCALAR_KIND) function median(a)

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(in) :: a(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: n

  n = size(a)
  if (mod(n, 2) == 1) then
     median = a((n + 1) / 2)
  else
     median = 0.5_wp * (a(n / 2) + a(n / 2 + 1))
  end if

end function median

real(SCALAR_KIND) function meanTrimmed(a)

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(inout) :: a(:)

  ! <<< Scalar variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, n
  real(wp) :: firstQuartile, thirdQuartile

  interface
     real(SCALAR_KIND) function median(a)
       real(SCALAR_KIND), intent(in) :: a(:)
     end function median
  end interface

  n = size(a)

  if (mod(n, 2) == 0) then
     firstQuartile = median(a(1:n/2))
     thirdQuartile = median(a(n/2+1:n))
  else
     firstQuartile = median(a(1:(n-1)/2))
     thirdQuartile = median(a((n+1)/2+1:n))
  end if

  meanTrimmed = 0.0_wp
  n = 0

  do i = 1, size(a)
     if (a(i) >= firstQuartile .and. a(i) <= thirdQuartile) then
        meanTrimmed = meanTrimmed + a(i)
        n = n + 1
     end if
  end do

  if (n == 0) return
  meanTrimmed = meanTrimmed / real(n, wp)

end function meanTrimmed

subroutine testAdjointRelation(identifier, nDimensions, success, isPeriodic, direction, tolerance)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use StencilOperator_mod, only : t_StencilOperator
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  use PatchDescriptor_mod, only : t_PatchDescriptor
  use Patch_factory, only : t_PatchFactory, updatePatchFactories
  use Patch_mod, only : t_Patch
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use DragForce_mod, only : t_DragForce

  use Region_enum, only : FORWARD, ADJOINT

  use CNSHelper
  use RhsHelperImpl, only : addDragForceAdjointPenalty

  ! <<< Internal modules >>>
  use MPIHelper, only : pigeonhole
  use RandomNumber, only : random
  use ErrorHandler

  ! <<< Arguments >>>
  character(len = *), intent(in) :: identifier
  integer, intent(in) :: nDimensions
  logical, intent(out) :: success
  logical, intent(in) :: isPeriodic
  integer, intent(in) :: direction
  real(SCALAR_KIND), intent(in), optional :: tolerance

  ! <<< interface >>>
  interface
     real(SCALAR_KIND) function meanTrimmed(a)
       real(SCALAR_KIND), intent(inout) :: a(:)
     end function meanTrimmed

     subroutine sort(a)
       real(SCALAR_KIND), intent(inout) :: a(:)
     end subroutine sort
  end interface

  ! <<< Local derived type variables >>>
  type(t_SimulationFlags) :: simulationFlags
  type(t_SolverOptions) :: solverOptions
  type(t_Grid) :: grid
  type(t_State) :: state0, state1
  type(t_PatchDescriptor) :: patchDescriptor
  type(t_PatchFactory), allocatable :: patchFactories(:)
  class(t_Patch), pointer :: patch => null()
  ! type(t_FarFieldPatch) :: patch

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp) :: scalar1, scalar2, J0, J1, tolerance_, normBoundaryFactor,           &
              stepSizes(32), errorHistory(32), convergenceHistory(31), scalarHistory(32)
  integer :: i, j, k, l, m, gridSize(nDimensions, 1), nUnknowns, normalDirection, &
              patchIndex, gridIndex, forceDirection, sgn, errorCode, extent(6)
  real(SCALAR_KIND), allocatable :: F(:,:), stress0(:,:), stress1(:,:),           &
                                    temp1(:,:), temp2(:,:),                       &
                                    deltaConservedVariables(:,:), deltaPrimitiveVariables(:,:),&
                                    ones(:)
  SCALAR_TYPE, dimension(nDimensions) :: h, gridPerturbation
  character(len = STRING_LENGTH) :: errorMessage

  tolerance_ = 1.0E-11
  if( present(tolerance) ) tolerance_ = tolerance

  success = .true.

  ! set up simulation flags
  call simulationFlags%initialize()
  simulationFlags%enableController = .true.
  simulationFlags%enableFunctional = .true.
  simulationFlags%enableAdjoint = .true.
  ! randomize curvilinear domain
  ! simulationFlags%isDomainCurvilinear = (random(0, 2) == 0)
  simulationFlags%isDomainCurvilinear = .false.
  simulationFlags%viscosityOn = .true.
  simulationFlags%repeatFirstDerivative = .true. ! this is default value.
  simulationFlags%useTargetState = .true.
  simulationFlags%storeVelocityGradient = .true.

  ! randomize grid size
  ! Note that too small grid size will yield null matrix for stencil operators.
  gridSize(:,:) = 1
  do i = 1, nDimensions
     gridSize(i,1) = random(20, 40)
  end do

  ! initialize solver option, grid, and states
  call solverOptions%initialize(nDimensions, simulationFlags)
  solverOptions%costFunctionalType = "DRAG"
  solverOptions%discretizationType = trim(identifier)
  solverOptions%reynoldsNumberInverse = 1.0_wp / random(5.0_wp,100.0_wp)
  solverOptions%prandtlNumberInverse = 1.0_wp / random(0.1_wp,5.0_wp)
  solverOptions%powerLawExponent = random(0.1_wp, 5.0_wp)
  solverOptions%bulkViscosityRatio = random(0.1_wp, 5.0_wp)
  call grid%setup(1, gridSize(1:nDimensions,1), MPI_COMM_WORLD,                        &
       simulationFlags = simulationFlags)
  call grid%setupSpatialDiscretization(simulationFlags, solverOptions)

  ! randomize grid coordinates (didn't reflect periodicity)
  h = 1.0_wp / real(grid%globalSize(1:nDimensions)-1,wp)
  do i = 1, grid%nGridPoints
    call random_number(gridPerturbation)
    gridPerturbation = (2.0_wp * gridPerturbation - 1.0_wp) * 0.05_wp * h
    grid%coordinates(i,:) = grid%coordinates(i,:) + gridPerturbation
  end do

  ! grid is not randomized: do not put any argument in updateGrid!!
  call grid%update()
  ! print *, 'Jacobian range: (', minval(grid%jacobian), ', ', maxval(grid%jacobian), ')'
  ! print *, 'Metrics range: (', minval(grid%metrics), ', ', maxval(grid%metrics), ')'

  call state0%setup(grid, simulationFlags, solverOptions)
  call state1%setup(grid, simulationFlags, solverOptions)

  ! setup patch
  allocate(patchFactories(1))
  ! normalDirection = min(2, nDimensions)
  ! forceDirection = 1
  normalDirection = direction
  forceDirection = random(1, nDimensions)
  if (random(0, 1) == 0) forceDirection = - forceDirection
  extent(1::2) = 1
  extent(2::2) = -1
  if (random(0, 1) == 0) then
     extent(1+2*(normalDirection-1)) = 1
     extent(2+2*(normalDirection-1)) = 1
  else
     extent(1+2*(normalDirection-1)) = gridSize(normalDirection, 1)
     extent(2+2*(normalDirection-1)) = gridSize(normalDirection, 1)
     normalDirection = -normalDirection
  end if
  patchDescriptor = t_PatchDescriptor("testPatch", "COST_TARGET", 1, normalDirection,     &
  ! patchDescriptor = t_PatchDescriptor("testPatch", "COST_TARGET", 1, 0,     &
                      extent(1), extent(2), extent(3), extent(4), extent(5), extent(6))
  ! patchDescriptor = t_PatchDescriptor("testPatch", "COST_TARGET", 1, 0,     &
  !                                     1, -1, 1, -1, 1, -1)
  call patchDescriptor%validate(gridSize, simulationFlags,                             &
                                solverOptions, errorCode, errorMessage)
  if (errorCode .ne. 0) then
    print *, trim(errorMessage)
    success = .false.
    return
  end if
  call patchFactories(1)%connect(patch,trim(patchDescriptor%patchType))
  call patch%setup(1, MPI_COMM_WORLD, patchDescriptor,                                 &
       grid, simulationFlags, solverOptions)

  ! Randomize conserved variables.
  do i = 1, grid%nGridPoints
     state0%conservedVariables(i,1) = random(0.01_wp, 10.0_wp)
     do j = 1, nDimensions
        state0%conservedVariables(i,j+1) =                                              &
             state0%conservedVariables(i,1) * random(-10.0_wp, 10.0_wp)
     end do
     state0%conservedVariables(i, nDimensions + 2) = state0%conservedVariables(i,1) *   &
          random(0.01_wp, 10.0_wp) / solverOptions%ratioOfSpecificHeats +               &
          0.5_wp / state0%conservedVariables(i,1) *                                     &
          sum(state0%conservedVariables(i,2:nDimensions+1) ** 2)
  end do
  assert(all(state0%conservedVariables(:,1) > 0.0_wp))

  ! Compute dependent variables.
  call state0%update(grid, simulationFlags, solverOptions)

  ! Update Patch Factories.
  call updatePatchFactories(patchFactories, simulationFlags,                &
                            solverOptions, grid, state0)

  ! Randomize delta conserved variables.
  allocate(deltaConservedVariables(grid%nGridPoints, solverOptions%nUnknowns))
  allocate(deltaPrimitiveVariables(grid%nGridPoints, solverOptions%nUnknowns))
  do i = 1, grid%nGridPoints
    do j = 1, nDimensions + 2
       deltaPrimitiveVariables(i,j) = random(-1.0_wp, 1.0_wp)
    end do
  end do
  deltaConservedVariables(:,1) = deltaPrimitiveVariables(:,1)
  do j = 1, nDimensions
     deltaConservedVariables(:,j+1) = state0%conservedVariables(:,j+1) /                     &
          state0%conservedVariables(:,1) * deltaPrimitiveVariables(:,1) +                    &
          state0%conservedVariables(:,1) * deltaPrimitiveVariables(:,j+1)
  end do
  deltaConservedVariables(:,nDimensions+2) = state0%conservedVariables(:,nDimensions+2) /    &
       state0%conservedVariables(:,1) * deltaPrimitiveVariables(:,1) +                       &
       sum(state0%conservedVariables(:,2:nDimensions+1) *                                    &
       deltaPrimitiveVariables(:,2:nDimensions+1), dim = 2) +                          &
       state0%conservedVariables(:,1) / solverOptions%ratioOfSpecificHeats *                               &
       deltaPrimitiveVariables(:,nDimensions+2)

  ! Compute baseline stress tensor
  ! stress tensor needs to be transformed, to be represented on the computational surface.
  allocate(stress0(grid%nGridPoints, 1))
  allocate(stress1(grid%nGridPoints, 1))
  allocate(ones(grid%nGridPoints))
  ones = 1.0_wp
  J0 = 0.0_wp
  stress0 = 0.0_wp
  i = abs(normalDirection)
  j = abs(forceDirection)
  sgn = SIGN(1, forceDirection * normalDirection)
  stress0(:, 1) = sgn * grid%metrics(:,i+nDimensions*(i-1)) *                     &
                  state0%dynamicViscosity(:,1) *                                  &
                  state0%velocityGradient(:,i+nDimensions*(j-1))
  select type (patch)
    class is (t_CostTargetPatch)
      J0 = patch%computeInnerProduct(grid, stress0(:,1), ones)
  end select

  ! Compute adjoint rhs for viscous flux
  nUnknowns = solverOptions%nUnknowns
  allocate(temp1(grid%nGridPoints, nUnknowns))
  allocate(temp2(grid%nGridPoints, 1))
  temp1 = 0.0_wp
  temp2 = 0.0_wp

  i = abs(normalDirection)
  j = abs(forceDirection)
  sgn = SIGN(1, forceDirection * normalDirection)
  ! normBoundaryFactor = 1.0_wp
  normBoundaryFactor = 1.0_wp / grid%firstDerivative(abs(patch%normalDirection))%normBoundary(1)
  temp2(:, 1) = - sgn * normBoundaryFactor *                                     &
                     grid%jacobian(:, 1) *                                          &
                     grid%metrics(:,i+nDimensions*(i-1)) *                          &
                     state0%dynamicViscosity(:,1) *                                  &
                     state0%velocityGradient(:,i+nDimensions*(j-1)) *                &
                     solverOptions%powerLawExponent /                               &
                     state0%temperature(:, 1)

  temp1(:, nUnknowns) = temp2(:, 1) *                                   &
                     solverOptions%ratioOfSpecificHeats *            &
                     state0%specificVolume(:, 1)
  do k = 1, nDimensions
    temp1(:, k+1) = - state0%velocity(:, k) * temp1(:, nUnknowns)
  end do
  temp1(:, 1) = - solverOptions%ratioOfSpecificHeats *                              &
                    state0%specificVolume(:, 1) * state0%specificVolume(:, 1) *      &
                    state0%conservedVariables(:, nUnknowns)
  temp1(:, 1) = temp1(:, 1) * temp2(:, 1) -                                         &
                sum(state0%velocity * temp1(:, 2:nDimensions+1), dim=2)

  select type (patch)
    class is (t_CostTargetPatch)
      call patch%collect(temp1, patch%adjointForcing)
  end select
  
  state0%rightHandSide = 0.0_wp
  ! (3) Add patch penalties
  select type (patch)
    class is (t_CostTargetPatch)
      call patch%updateRhs(ADJOINT, simulationFlags, solverOptions, grid, state0)
  end select

  !! This part is moved to addDragForceAdjointPenalty.
  l = abs(normalDirection)
  m = abs(forceDirection)
  sgn = SIGN(1, forceDirection * normalDirection)
  normBoundaryFactor = 1.0_wp / grid%firstDerivative(abs(patch%normalDirection))%normBoundary(1)
  
  temp2 = 0.0_wp
  do k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)
    do j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
        do i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
          gridIndex = i - patch%gridOffset(1) + patch%gridLocalSize(1) *         &
              (j - 1 - patch%gridOffset(2) + patch%gridLocalSize(2) *           &
              (k - 1 - patch%gridOffset(3)))
          if (grid%iblank(gridIndex) == 0) cycle
          patchIndex = i - patch%offset(1) + patch%localSize(1) *                &
              (j - 1 - patch%offset(2) + patch%localSize(2) *                   &
              (k - 1 - patch%offset(3)))

          !   temp2(:, 1) = - sgn * grid%targetMollifier(:, 1) *                 &
          temp2(gridIndex, 1) = temp2(gridIndex, 1) -                               &
                                  sgn * normBoundaryFactor *                        &
                                  state0%dynamicViscosity(gridIndex, 1) *            &
                                  grid%metrics(gridIndex, l+nDimensions*(l-1))**2 * &
                                  grid%jacobian(gridIndex, 1)

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
    end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

  call grid%adjointFirstDerivative(l)%apply(temp2, grid%localSize)
  temp2(:, 1) = temp2(:, 1) * grid%jacobian(:, 1)! * normBoundaryFactor

  state0%rightHandSide(:, 1) = state0%rightHandSide(:, 1) -                               &
                        state0%velocity(:, m) * state0%specificVolume(:, 1) * temp2(:, 1)
  state0%rightHandSide(:, m+1) = state0%rightHandSide(:, m+1) +                           &
                                  state0%specificVolume(:, 1) * temp2(:, 1)

  SAFE_DEALLOCATE(temp1)
  SAFE_DEALLOCATE(temp2)

  ! <R^{\dagger}u, \delta v>
  scalar1 = grid%computeInnerProduct(state0%rightHandSide,deltaConservedVariables)

  ! <u, \delta R(v)>
  ! Prepare step sizes
  stepSizes(1) = 0.01_wp
  do k = 2, size(stepSizes)
     stepSizes(k) = stepSizes(k-1) * 10.0_wp**(-0.25_wp)
  end do
  do k = 1, size(stepSizes)
    !(1) finite difference on conserved variables
    state1%conservedVariables = state0%conservedVariables + stepSizes(k) * deltaConservedVariables
    ! assert(all(state1%conservedVariables(:,1) > 0.0_wp))
    ! Compute dependent variables.
    call state1%update(grid, simulationFlags, solverOptions)
    
    ! ! No deviation from second jacobian of viscous flux
    ! state1%velocityGradient = state0%velocityGradient
    ! ! No deviation from dynamic viscosity
    ! state1%dynamicViscosity = state0%dynamicViscosity

    ! (2) Compute baseline SAT far-field
    stress1 = 0.0_wp
    i = abs(normalDirection)
    j = abs(forceDirection)
    sgn = SIGN(1, forceDirection * normalDirection)
    stress1(:, 1) = sgn * grid%metrics(:,i+nDimensions*(i-1)) *                     &
                    state1%dynamicViscosity(:,1) *                                  &
                    state1%velocityGradient(:,i+nDimensions*(j-1))
    select type (patch)
      class is (t_CostTargetPatch)
        J1 = patch%computeInnerProduct(grid, stress1(:,1), ones)
    end select

    ! (3) <u, \delta R(v)>
    scalar2 = J1 - J0

    errorHistory(k) = 0.0_wp
    scalarHistory(k) = scalar2 / stepSizes(k)
    if( abs(scalar1)>0.0_wp ) errorHistory(k) = abs( (scalar2/stepSizes(k) + scalar1)/scalar1 )

    if (k > 1) then
       convergenceHistory(k-1) = log(errorHistory(k) / errorHistory(k-1)) /              &
            log(stepSizes(k) / stepSizes(k-1))
       if (k > 5) then
         if(sum(convergenceHistory(k-3:k-1))/3.0_wp < 0.0_wp) exit
       end if
    end if
  end do

  if (maxval(errorHistory).le.tolerance_) then
    success = .true.
  elseif (k > 3) then
     call sort(convergenceHistory(:k-2))
     success = success .and. nint(meanTrimmed(convergenceHistory(:k-2))).ge.1
  else
     success = .false.
     print *, 'minimum error: ', minval(errorHistory)
  end if

  if (.not. success) then
    do i = 1, min(k, size(stepSizes))
      write(errorMessage,'(E8.3,3(3X,E32.15))') stepSizes(i), scalar1,                         &
                                            scalarHistory(i), errorHistory(i)
      call writeAndFlush(MPI_COMM_WORLD,output_unit,errorMessage)
    end do
  end if

  SAFE_DEALLOCATE(deltaConservedVariables)
  SAFE_DEALLOCATE(deltaPrimitiveVariables)

  call patch%cleanup()
  call patchFactories(1)%cleanup()
  SAFE_DEALLOCATE(patchFactories)
  call state0%cleanup()
  call state1%cleanup()
  call grid%cleanup()

end subroutine testAdjointRelation
