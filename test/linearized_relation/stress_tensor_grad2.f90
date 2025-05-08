#include "config.h"

program stress_tensor_grad

  use MPI

  use CNSHelper
  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random

  implicit none

  logical :: success, success_, isPeriodic
  integer :: i, j, nDimensions, ierror
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

     subroutine testLinearizedRelation(identifier, nDimensions, success, isPeriodic, tolerance)

       character(len = *), intent(in) :: identifier
       integer, intent(in) :: nDimensions
       logical, intent(out) :: success

       logical, intent(in), optional :: isPeriodic
       real(SCALAR_KIND), intent(in), optional :: tolerance

     end subroutine testLinearizedRelation

  end interface

  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  do nDimensions = 2, 3
    do j = 1, 4 !... for each discretizationTypes
      success = .true.
      do i = 1, 10 !... test multiple times
        ! Didn't test periodic grid yet!!
        ! isPeriodic = .true.
        ! call testLinearizedRelation(discretizationTypes(j), nDimensions,           &
        !                          success_, isPeriodic)
        ! success = success .and. success_
        ! if( .not. success) then
        !   if( procRank == 0 ) then
        !     print *, 'Failed, ', trim(discretizationTypes(j))
        !     print *, 'dimension: ', nDimensions
        !     print *, 'periodicity: ', isPeriodic
        !   end if
        !   exit
        ! end if

        isPeriodic = .false.
        call testLinearizedRelation(discretizationTypes(j), nDimensions,           &
                                 success_, isPeriodic)
        success = success .and. success_
        if( .not. success_) then
          if( procRank == 0 ) then
            print *, 'Failed, ', trim(discretizationTypes(j))
            print *, 'dimension: ', nDimensions
            print *, 'periodicity: ', isPeriodic
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

end program stress_tensor_grad

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

subroutine testLinearizedRelation(identifier, nDimensions, success, isPeriodic, tolerance)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use StencilOperator_mod, only : t_StencilOperator
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  use Region_enum, only : FORWARD, ADJOINT

  use CNSHelper

  ! <<< Internal modules >>>
  use MPIHelper, only : pigeonhole
  use RandomNumber, only : random
  use ErrorHandler

  ! <<< Arguments >>>
  character(len = *), intent(in) :: identifier
  integer, intent(in) :: nDimensions
  logical, intent(out) :: success
  logical, intent(in), optional :: isPeriodic
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

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp) :: scalar1, scalar2, J0, J1, scalarHistory(32),                            &
              stepSizes(32), errorHistory(32), convergenceHistory(31)
  integer :: i, j, k, l, gridSize(nDimensions, 1), nUnknowns, &
             normalDirection, forceDirection
  real(SCALAR_KIND), allocatable :: F(:,:), stress0(:,:), stress1(:,:),           &
                                    temp1(:,:), temp2(:,:),                       &
                                    deltaConservedVariables(:,:), deltaPrimitiveVariables(:,:),&
                                    ones(:)
  SCALAR_TYPE :: tmp
  SCALAR_TYPE, dimension(nDimensions) :: h, gridPerturbation
  real(SCALAR_KIND), allocatable :: patchNorm(:,:)
  character(len = STRING_LENGTH) :: errorMessage

  success = .true.

  ! decide stress direction
  call random_number(tmp)
  normalDirection = 1 + FLOOR(nDimensions * tmp)
  ! select force direction perpendicular to normal direction.
  forceDirection = normalDirection
  do while (forceDirection .eq. normalDirection)
    call random_number(tmp)
    forceDirection = 1 + FLOOR(nDimensions * tmp)
  end do

  ! set up simulation flags
  call simulationFlags%initialize()
  simulationFlags%enableController = .true.
  simulationFlags%enableFunctional = .true.
  simulationFlags%enableAdjoint = .true.
  ! we only consider non-curvilinear case.
  simulationFlags%isDomainCurvilinear = .false.
  simulationFlags%viscosityOn = .true.
  simulationFlags%repeatFirstDerivative = .true. ! this is default value.
  simulationFlags%storeVelocityGradient = .true.

  ! randomize grid size
  ! Note that too small grid size will yield null matrix for stencil operators.
  gridSize(:,:) = 1
  do i = 1, nDimensions
     gridSize(i,1) = random(20, 40)
  end do

  ! initialize solver option, grid, and states
  call solverOptions%initialize(nDimensions, simulationFlags)
  solverOptions%discretizationType = trim(identifier)
  solverOptions%reynoldsNumberInverse = 1.0_wp / 5.0_wp
  solverOptions%prandtlNumberInverse = 1.0_wp / 0.72_wp
  solverOptions%powerLawExponent = 0.666_wp
  solverOptions%bulkViscosityRatio = 0.6_wp
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

  ! Randomize conserved variables.
  do i = 1, grid%nGridPoints
     state0%conservedVariables(i,1) = random(0.01_wp, 10.0_wp)
     do j = 1, nDimensions
        state0%conservedVariables(i,j+1) =                                              &
             state0%conservedVariables(i,1) * random(-10.0_wp, 10.0_wp)
     end do
     state0%conservedVariables(i, nDimensions + 2) = state0%conservedVariables(i,1) *    &
          random(0.01_wp, 10.0_wp) / solverOptions%ratioOfSpecificHeats +              &
          0.5_wp / state0%conservedVariables(i,1) *                                     &
          sum(state0%conservedVariables(i,2:nDimensions+1) ** 2)
  end do
  assert(all(state0%conservedVariables(:,1) > 0.0_wp))

  ! Compute dependent variables.
  call state0%update(grid,simulationFlags,solverOptions)

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

  ! norm on computational domain. stress is integrated on here.
  ! On production run, this norm can be found from t_CostTargetPatch.
  allocate(patchNorm(grid%nGridPoints, 1))
  patchNorm = 1.0_wp
  ! volumetric integral for now.
  do i = 1, nDimensions
    call grid%firstDerivative(i)%applyNorm(patchNorm, grid%localSize)
  end do

  ! Compute baseline stress tensor
  ! stress tensor needs to be transformed, to be represented on the computational surface.
  allocate(stress0(grid%nGridPoints, 1))
  allocate(stress1(grid%nGridPoints, 1))
  allocate(ones(grid%nGridPoints))
  ones = 1.0_wp
  J0 = 0.0_wp
  stress0 = 0.0_wp
  i = normalDirection
  j = forceDirection
  stress0(:, 1) = grid%metrics(:,i+nDimensions*(i-1)) *                           &
                  state0%dynamicViscosity(:,1) *                                  &
                  state0%velocityGradient(:,i+nDimensions*(j-1))
  J0 = sum(stress0(:,1) * patchNorm(:,1) * ones)

  ! Compute jacobian for stress tensor
  nUnknowns = solverOptions%nUnknowns

  ! (2) A * dQ
  allocate(temp1(grid%nGridPoints, solverOptions%nUnknowns))
  temp1 = 0.0_wp

  temp1(:, nUnknowns) = solverOptions%ratioOfSpecificHeats * state0%specificVolume(:, 1)
  do i = 1, nDimensions
    temp1(:, i+1) = - state0%velocity(:, i) * temp1(:, nUnknowns)
  end do
  temp1(:, 1) = - solverOptions%ratioOfSpecificHeats *                              &
                  state0%specificVolume(:, 1) * state0%specificVolume(:, 1) *       &
                  state0%conservedVariables(:, nUnknowns)
  temp1(:, 1) = temp1(:, 1) - sum(state0%velocity * temp1(:, 2:nDimensions+1), dim=2)

  i = normalDirection
  j = forceDirection
  do k = 1, nUnknowns
    temp1(:, k) = temp1(:, k) *                                                     &
                    grid%metrics(:,i+nDimensions*(i-1)) *                           &
                    state0%dynamicViscosity(:,1) *                                  &
                    state0%velocityGradient(:,i+nDimensions*(j-1)) *                &
                    solverOptions%powerLawExponent /                                &
                    state0%temperature(:, 1)
  end do

  ! (3) B * D * C * dQ
  allocate(temp2(grid%nGridPoints, 1))
  i = normalDirection
  j = forceDirection

  temp2(:, 1) = state0%dynamicViscosity(:, 1) *                     &
                  grid%metrics(:, i+nDimensions*(i-1))**2 *         &
                  grid%jacobian(:, 1)
  call grid%adjointFirstDerivative(i)%apply(temp2, grid%localSize)

  temp1(:, 1) = temp1(:, 1) - state0%velocity(:, j) * state0%specificVolume(:, 1) * temp2(:, 1)
  temp1(:, j+1) = temp1(:, j+1) + state0%specificVolume(:, 1) * temp2(:, 1)

  ! <u, \partial R\delta v>
  scalar1 = 0.0_wp
  do k = 1, size(temp1, 2)
    scalar1 = scalar1 + sum(temp1(:,k) * patchNorm(:,1) * deltaConservedVariables(:,k))
  end do

  SAFE_DEALLOCATE(temp1)
  SAFE_DEALLOCATE(temp2)

  ! <u, \delta R(v)>
  ! Prepare step sizes
  stepSizes(1) = 0.001_wp
  do k = 2, size(stepSizes)
     stepSizes(k) = stepSizes(k-1) * 10.0_wp**(-0.25_wp)
  end do
  do k = 1, size(stepSizes)
    !(1) finite difference on conserved variables
    state1%conservedVariables = state0%conservedVariables + stepSizes(k) * deltaConservedVariables
    assert(all(state1%conservedVariables(:,1) > 0.0_wp))
    ! Compute dependent variables.
    !! update all dependent/transport variables.
    call state1%update(grid,simulationFlags,solverOptions)
    !! update dependent variable after derivative only.
    ! call computeDependentVariables(nDimensions, state1%conservedVariables,                    &
    !      solverOptions%ratioOfSpecificHeats, state1%specificVolume(:,1), state1%velocity,       &
    !      state1%pressure(:,1), state1%temperature(:,1))
    ! call computeTransportVariables(state0%temperature(:,1), solverOptions%powerLawExponent,   &
    !      solverOptions%bulkViscosityRatio, solverOptions%ratioOfSpecificHeats,              &
    !      solverOptions%reynoldsNumberInverse, solverOptions%prandtlNumberInverse,           &
    !      state1%dynamicViscosity(:,1), state1%secondCoefficientOfViscosity(:,1),                &
    !      state1%thermalDiffusivity(:,1))
    ! call grid%computeGradient(state1%velocity, state1%stressTensor)
    ! call computeStressTensor(nDimensions, state1%stressTensor, state1%dynamicViscosity(:,1), &
    !      state1%secondCoefficientOfViscosity(:,1))
    
    i = normalDirection
    j = forceDirection
    stress1(:, 1) = grid%metrics(:,i+nDimensions*(i-1)) *                           &
                    state1%dynamicViscosity(:,1) *                                  &
                    state1%velocityGradient(:,i+nDimensions*(j-1))
    J1 = sum(stress1(:,1) * patchNorm(:,1) * ones)

    scalar2 = J1 - J0

    scalarHistory(k) = scalar2/stepSizes(k)
    errorHistory(k) = abs( (scalar2/stepSizes(k) - scalar1)/scalar1 )

    if (k > 1) then
       convergenceHistory(k-1) = log(errorHistory(k) / errorHistory(k-1)) /              &
            log(stepSizes(k) / stepSizes(k-1))
       if (k > 5) then
         if ( sum(convergenceHistory(k-3:k-1))/3.0_wp < 0.0_wp) exit
       end if
    end if
  end do

  if (k > 2) then
     call sort(convergenceHistory(:k-2))
     success = success .and. nint(meanTrimmed(convergenceHistory(:k-2))).ge.1
  else
     success = .false.
  end if

  if (.not. success) then
    do i = 1, k
      write(errorMessage,'(E8.3,3(3X,E32.15))') stepSizes(i), scalar1,                         &
                                            scalarHistory(i), errorHistory(i)
      call writeAndFlush(MPI_COMM_WORLD,output_unit,errorMessage)
    end do
  end if

  SAFE_DEALLOCATE(deltaConservedVariables)
  SAFE_DEALLOCATE(deltaPrimitiveVariables)

  call state0%cleanup()
  call state1%cleanup()
  call grid%cleanup()

end subroutine testLinearizedRelation
