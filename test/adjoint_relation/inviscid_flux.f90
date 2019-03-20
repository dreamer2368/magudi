#include "config.h"

program inviscid_flux

  use MPI

  use CNSHelper
  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random

  implicit none

  logical :: success, success_, isPeriodic
  integer :: i, j, k, nDimensions, ierror
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

     subroutine testAdjointRelation(identifier, nDimensions, success, isPeriodic, tolerance)

       character(len = *), intent(in) :: identifier
       integer, intent(in) :: nDimensions
       logical, intent(out) :: success

       logical, intent(in) :: isPeriodic
       real(SCALAR_KIND), intent(in), optional :: tolerance

     end subroutine testAdjointRelation

  end interface

  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  do nDimensions = 1, 3
    do j = 1, 4 !... for each discretizationTypes
      success = .true.
      do i = 1, 10 !... test multiple times
        ! Didn't test periodic grid yet!!
        ! isPeriodic = .true.
        ! call testAdjointRelation(discretizationTypes(j), nDimensions,           &
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
        call testAdjointRelation(discretizationTypes(j), nDimensions,           &
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

end program inviscid_flux

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

subroutine testAdjointRelation(identifier, nDimensions, success, isPeriodic, tolerance)

  ! <<< External modules >>>
  use MPI

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

  ! <<< Arguments >>>
  character(len = *), intent(in) :: identifier
  integer, intent(in) :: nDimensions
  logical, intent(out) :: success
  logical, intent(in) :: isPeriodic
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
  logical :: isPeriodic_(3), hasNegativeJacobian
  real(wp) :: scalar1, scalar2, tolerance_,&
              stepSizes(32), errorHistory(32), convergenceHistory(31)
  integer :: i, j, k, gridSize(nDimensions, 1), nUnknowns
  real(SCALAR_KIND), allocatable :: F(:,:), fluxes1(:,:,:), fluxes2(:,:,:),           &
                                    temp1(:,:,:), localFluxJacobian1(:,:),            &
                                    localConservedVariables(:), localVelocity(:),     &
                                    localMetricsAlongDirection1(:),                   &
                                    adjointRightHandSide(:,:),                        &
                                    deltaConservedVariables(:,:), deltaPrimitiveVariables(:,:)
  SCALAR_TYPE, dimension(nDimensions) :: h, gridPerturbation
  character(len = STRING_LENGTH) :: errorMessage

  success = .true.

  ! set up simulation flags
  call simulationFlags%initialize()
  simulationFlags%enableController = .true.
  simulationFlags%enableFunctional = .true.
  simulationFlags%enableAdjoint = .true.
  simulationFlags%isDomainCurvilinear = .true.

  ! randomize grid size
  ! Note that too small grid size will yield null matrix for stencil operators.
  gridSize(:,:) = 1
  do i = 1, nDimensions
     gridSize(i,1) = random(20, 40)
  end do

  ! initialize solver option, grid, and states
  call solverOptions%initialize(nDimensions, simulationFlags)
  solverOptions%discretizationType = trim(identifier)
  call grid%setup(1, gridSize(:,1), MPI_COMM_WORLD,                        &
       simulationFlags = simulationFlags)
  call grid%setupSpatialDiscretization(simulationFlags, solverOptions)

  ! randomize grid coordinates (didn't reflect periodicity)
  h = 1.0_wp / real(grid%globalSize(1:nDimensions)-1,wp)
  do i = 1, grid%nGridPoints
    call random_number(gridPerturbation)
    gridPerturbation = (2.0_wp * gridPerturbation - 1.0_wp) * 0.13_wp * h
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
  call computeDependentVariables(nDimensions, state0%conservedVariables,                &
    solverOptions%ratioOfSpecificHeats, state0%specificVolume(:,1), state0%velocity,     &
    state0%pressure(:,1), state0%temperature(:,1))
  assert(all(state0%specificVolume(:,1) > 0.0_wp))
  assert(all(state0%temperature(:,1) > 0.0_wp))

  ! Randomize adjoint variables.
  allocate(F(grid%nGridPoints, solverOptions%nUnknowns))
  call random_number(F)
  state0%adjointVariables = F
  SAFE_DEALLOCATE(F)

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

  ! Compute baseline inviscid flux
  ! (1) Cartesian form
  allocate(fluxes1(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))
  allocate(fluxes2(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))
  call computeCartesianInviscidFluxes(nDimensions, state0%conservedVariables,                 &
       state0%velocity, state0%pressure(:,1), fluxes1)
  ! (2) Transform fluxes
  call transformFluxes(nDimensions, fluxes1, grid%metrics, fluxes2, grid%isCurvilinear)
  ! (3) Take derivatives of fluxes
  do i = 1, nDimensions
     call grid%firstDerivative(i)%apply(fluxes2(:,:,i), grid%localSize)
  end do
  state0%rightHandSide = - sum(fluxes2, dim = 3) !... update RHS.
  SAFE_DEALLOCATE(fluxes1)
  SAFE_DEALLOCATE(fluxes2)
  ! (4) Multiply by Jacobian
  do j = 1, solverOptions%nUnknowns
     state0%rightHandSide(:,j) = state0%rightHandSide(:,j) * grid%jacobian(:,1)
  end do

  ! Compute adjoint rhs for inviscid flux
  nUnknowns = solverOptions%nUnknowns
  allocate(adjointRightHandSide(grid%nGridPoints, solverOptions%nUnknowns))
  allocate(temp1(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))
  ! (1) Partial derivatives of adjoint variables w.r.t. *computational* coordinates.
  adjointRightHandSide = 0.0_wp
  do i = 1, nDimensions
     temp1(:,:,i) = state0%adjointVariables
     call grid%adjointFirstDerivative(i)%apply(temp1(:,:,i), grid%localSize)
  end do
  ! (2) Compute jacobian of inviscid flux dotted with adjoint derivatives
  allocate(localFluxJacobian1(nUnknowns, nUnknowns))
  allocate(localConservedVariables(nUnknowns))
  allocate(localVelocity(nDimensions))
  allocate(localMetricsAlongDirection1(nDimensions))
  do j = 1, grid%nGridPoints
     localConservedVariables = state0%conservedVariables(j,:)
     localVelocity = state0%velocity(j,:)
     do i = 1, nDimensions
        localMetricsAlongDirection1 = grid%metrics(j,1+nDimensions*(i-1):nDimensions*i)
        select case (nDimensions)
        case (1)
           call computeJacobianOfInviscidFlux1D(localConservedVariables,                     &
                localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                localFluxJacobian1, specificVolume = state0%specificVolume(j,1),              &
                velocity = localVelocity, temperature = state0%temperature(j,1))
        case (2)
           call computeJacobianOfInviscidFlux2D(localConservedVariables,                     &
                localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                localFluxJacobian1, specificVolume = state0%specificVolume(j,1),              &
                velocity = localVelocity, temperature = state0%temperature(j,1))
        case (3)
           call computeJacobianOfInviscidFlux3D(localConservedVariables,                     &
                localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                localFluxJacobian1, specificVolume = state0%specificVolume(j,1),              &
                velocity = localVelocity, temperature = state0%temperature(j,1))
        end select !... nDimensions

        adjointRightHandSide(j,:) = adjointRightHandSide(j,:) +                                &
             matmul(transpose(localFluxJacobian1), temp1(j,:,i))
     end do !... i = 1, nDimensions
  end do !... j = 1, grid%nGridPoints
  ! (3) Multiply by Jacobian
  do j = 1, solverOptions%nUnknowns
     adjointRightHandSide(:,j) = adjointRightHandSide(:,j) * grid%jacobian(:,1)
  end do

  SAFE_DEALLOCATE(localFluxJacobian1)
  SAFE_DEALLOCATE(localConservedVariables)
  SAFE_DEALLOCATE(localVelocity)
  SAFE_DEALLOCATE(localMetricsAlongDirection1)
  SAFE_DEALLOCATE(temp1)

  ! <R^{\dagger}u, \delta v>
  scalar1 = grid%computeInnerProduct(adjointRightHandSide,deltaConservedVariables)
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
    call computeDependentVariables(nDimensions, state1%conservedVariables,                &
      solverOptions%ratioOfSpecificHeats, state1%specificVolume(:,1), state1%velocity,     &
      state1%pressure(:,1), state1%temperature(:,1))
    assert(all(state1%specificVolume(:,1) > 0.0_wp))
    assert(all(state1%temperature(:,1) > 0.0_wp))

    ! (2)Compute baseline inviscid flux
    ! (2-1) Cartesian form
    allocate(fluxes1(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))
    allocate(fluxes2(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))
    ! print *, size(fluxes1,1),size(fluxes1,2),size(fluxes1,3)
    call computeCartesianInviscidFluxes(nDimensions, state1%conservedVariables,                 &
         state1%velocity, state1%pressure(:,1), fluxes1)
    ! (2-2) Transform fluxes
    call transformFluxes(nDimensions, fluxes1, grid%metrics, fluxes2, grid%isCurvilinear)
    ! (2-3) Take derivatives of fluxes
    do i = 1, nDimensions
       call grid%firstDerivative(i)%apply(fluxes2(:,:,i), grid%localSize)
    end do
    state1%rightHandSide = - sum(fluxes2, dim = 3) !... update RHS.
    SAFE_DEALLOCATE(fluxes1)
    SAFE_DEALLOCATE(fluxes2)
    ! (2-4) Multiply by Jacobian
    do j = 1, solverOptions%nUnknowns
       state1%rightHandSide(:,j) = state1%rightHandSide(:,j) * grid%jacobian(:,1)
    end do

    ! (3) <u, \delta R(v)>
    scalar2 = grid%computeInnerProduct(state0%adjointVariables,                             &
                                        state1%rightHandSide - state0%rightHandSide)

    errorHistory(k) = abs( (scalar2/stepSizes(k) + scalar1)/scalar1 )
    ! print *, stepSizes(k), -scalar1, scalar2/stepSizes(k), errorHistory(k)

    if (k > 1) then
       convergenceHistory(k-1) = log(errorHistory(k) / errorHistory(k-1)) /              &
            log(stepSizes(k) / stepSizes(k-1))
       if (k > 5) then
         if (sum(convergenceHistory(k-3:k-1))/3.0_wp < 0.0_wp) exit
       end if
    end if
  end do

  if (k > 2) then
     call sort(convergenceHistory(:k-2))
     success = success .and. nint(meanTrimmed(convergenceHistory(:k-2))).ge.1
  else
     success = .false.
  end if

  SAFE_DEALLOCATE(adjointRightHandSide)
  SAFE_DEALLOCATE(deltaConservedVariables)
  SAFE_DEALLOCATE(deltaPrimitiveVariables)

  call state0%cleanup()
  call state1%cleanup()
  call grid%cleanup()

end subroutine testAdjointRelation
