#include "config.h"

program SAT_farfield_linearized

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

     subroutine testLinearizedRelation(identifier, nDimensions, success, isPeriodic, direction, tolerance)

       character(len = *), intent(in) :: identifier
       integer, intent(in) :: nDimensions
       logical, intent(out) :: success

       logical, intent(in) :: isPeriodic
       integer, intent(in) :: direction
       real(SCALAR_KIND), intent(in), optional :: tolerance

     end subroutine testLinearizedRelation

  end interface

  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  do nDimensions = 1, 3
    do j = 1, 4!... for each discretizationTypes
      do direction = 1, nDimensions
        success = .true.
        isPeriodic = .false.
        do i = 1, 10 !... test multiple times

          call testLinearizedRelation(discretizationTypes(j), nDimensions,           &
                                   success_, isPeriodic, direction)
          success = success .and. success_
          if( .not. success_) then
            exit
          end if
        end do
        if( .not. success) then
          if( procRank == 0 ) then
            print *, 'Failed, ', trim(discretizationTypes(j))
            print *, 'dimension: ', nDimensions
            print *, 'periodicity: ', isPeriodic
            print *, 'direction: ', direction
          end if
        end if
      end do !... for each direction
      if( procRank == 0 .and. success ) then
        print *, 'Success, ', trim(discretizationTypes(j))
        print *, 'dimension: ', nDimensions
      end if
    end do !... for each discretizationTypes
  end do

  call cleanupErrorHandler()

  call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierror)
  call MPI_Finalize(ierror)
  if (.not. success) stop -1
  stop 0

end program SAT_farfield_linearized

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

subroutine testLinearizedRelation(identifier, nDimensions, success, isPeriodic, direction, tolerance)

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
  use Patch_factory, only : t_PatchFactory
  use Patch_mod, only : t_Patch
  use FarFieldPatch_mod, only : t_FarFieldPatch

  use Region_enum, only : FORWARD, ADJOINT, LINEARIZED

  use CNSHelper
  use RhsHelperImpl, only : addFarFieldAdjointPenalty

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
  real(wp) :: scalar1, scalar2, tolerance_, scalarHistory(32),                             &
              stepSizes(32), errorHistory(32), convergenceHistory(31)
  integer :: i, j, k, gridSize(nDimensions, 1), nUnknowns, direction_, errorCode, extent(6)
  real(SCALAR_KIND), allocatable :: F(:,:), fluxes2(:,:,:),                           &
                                    temp1(:,:,:), localFluxJacobian1(:,:),            &
                                    localConservedVariables(:), localVelocity(:),     &
                                    localMetricsAlongDirection1(:),                   &
                                    localMetricsAlongDirection2(:),                   &
                                    linearizedRightHandSide(:,:),                        &
                                    deltaConservedVariables(:,:), deltaPrimitiveVariables(:,:),&
                                    localFluxJacobian2(:,:), localStressTensor(:),    &
                                    localHeatFlux(:), localLinearizedDiffusion(:,:),     &
                                    temp2(:,:,:), temp3(:,:), targetViscousFluxes(:,:,:)
  SCALAR_TYPE, dimension(nDimensions) :: h, gridPerturbation
  character(len = STRING_LENGTH) :: errorMessage

  tolerance_ = 1.0E-9
  if( present(tolerance) ) tolerance_ = tolerance

  success = .true.

  ! set up simulation flags
  call simulationFlags%initialize()
  simulationFlags%enableController = .true.
  simulationFlags%enableFunctional = .true.
  simulationFlags%enableAdjoint = .true.
  ! randomize curvilinear domain
  ! simulationFlags%isDomainCurvilinear = (random(0, 2) == 0)
  simulationFlags%isDomainCurvilinear = .true.
  simulationFlags%viscosityOn = .true.
  simulationFlags%repeatFirstDerivative = .true. ! this is default value.
  simulationFlags%useTargetState = .true.

  ! randomize grid size
  ! Note that too small grid size will yield null matrix for stencil operators.
  gridSize(:,:) = 1
  do i = 1, nDimensions
     gridSize(i,1) = random(20, 40)
  end do

  ! initialize solver option, grid, and states
  call solverOptions%initialize(nDimensions, simulationFlags)
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

  call state0%setup(grid, simulationFlags, solverOptions)
  call state1%setup(grid, simulationFlags, solverOptions)

  ! setup patch
  allocate(patchFactories(1))
  direction_ = direction
  extent(1::2) = 1
  extent(2::2) = -1
  if (random(0, 1) == 0) then
     extent(1+2*(direction_-1)) = 1
     extent(2+2*(direction_-1)) = 1
  else
     extent(1+2*(direction_-1)) = gridSize(direction_, 1)
     extent(2+2*(direction_-1)) = gridSize(direction_, 1)
     direction_ = -direction_
  end if
  patchDescriptor = t_PatchDescriptor("testPatch", "SAT_FAR_FIELD", 1, direction_,     &
       extent(1), extent(2), extent(3), extent(4), extent(5), extent(6))
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
  select type(patch)
    class is (t_FarFieldPatch)
      patch%inviscidPenaltyAmount = sign(random(0.01_wp,10.0_wp),                       &
                                      real(patch%normalDirection,wp))                   &
                                      /grid%firstDerivative(abs(direction_))%normBoundary(1)
      if (simulationFlags%viscosityOn)                                                  &
        patch%viscousPenaltyAmount = sign(random(0.01_wp,10.0_wp),                      &
                                     real(patch%normalDirection,wp))                    &
                                     /grid%firstDerivative(abs(direction_))%normBoundary(1)
  end select

  ! Randomize target state.
  do i = 1, grid%nGridPoints
    state0%targetState(i,1) = random(0.01_wp, 10.0_wp)
    do j = 1, nDimensions
       state0%targetState(i,j+1) =                                              &
            state0%targetState(i,1) * random(-10.0_wp, 10.0_wp)
    end do
    state0%targetState(i, nDimensions + 2) = state0%targetState(i,1) *   &
         random(0.01_wp, 10.0_wp) / solverOptions%ratioOfSpecificHeats +               &
         0.5_wp / state0%targetState(i,1) *                                     &
         sum(state0%targetState(i,2:nDimensions+1) ** 2)
  end do
  assert(all(state0%targetState(:,1) > 0.0_wp))
  state1%targetState = state0%targetState
  if( simulationFlags%viscosityOn ) then
    allocate(targetViscousFluxes(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))
    call state0%update(grid, simulationFlags, solverOptions, state0%targetState)
    call computeCartesianViscousFluxes(nDimensions, state0%velocity,                      &
         state0%stressTensor, state0%heatFlux, targetViscousFluxes)
    call patchFactories(1)%connect(patch)
    select type (patch)
    class is (t_FarFieldPatch)
      call patch%collect(targetViscousFluxes, patch%targetViscousFluxes)
    end select
  end if

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

  ! Compute baseline SAT far-field
  ! (1) Cartesian form
  allocate(fluxes2(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))
  if (simulationFlags%viscosityOn .and. simulationFlags%repeatFirstDerivative) then
     call computeCartesianViscousFluxes(nDimensions, state0%velocity,                         &
          state0%stressTensor, state0%heatFlux, fluxes2)
  end if
  ! (2) Send viscous fluxes to patch
  if (simulationFlags%viscosityOn) then
    call patchFactories(1)%connect(patch)
    select type (patch)
      class is (t_FarFieldPatch)
        call patch%collect(fluxes2, patch%viscousFluxes)
    end select
  end if
  ! (3) Add patch penalty
  state0%rightHandSide = 0.0_wp
  select type (patch)
    class is (t_FarFieldPatch)
      call patch%updateRhs(FORWARD, simulationFlags, solverOptions, grid, state0)
  end select
  SAFE_DEALLOCATE(fluxes2)

  ! Compute adjoint rhs for viscous flux
  nUnknowns = solverOptions%nUnknowns
  allocate(linearizedRightHandSide(grid%nGridPoints, solverOptions%nUnknowns))
  allocate(temp3(grid%nGridPoints, solverOptions%nUnknowns))
  temp3 = state0%rightHandSide
  state0%rightHandSide = 0.0_wp
  state0%adjointVariables = deltaConservedVariables
  ! (1) Collect delta viscous fluxes on far-field patches: This is the same as linearized viscous flux.
  ! (1-1) D * C * dQ
  allocate(temp2(grid%nGridPoints, solverOptions%nUnknowns-1, nDimensions))
  temp2 = 0.0_wp
  ! temp2(:,1,1) = deltaConservedVariables(:,1)
  do i = 1, nDimensions
    temp2(:,i,1) = - state0%velocity(:,i) * deltaConservedVariables(:,1)        &
                  + deltaConservedVariables(:,i+1)
  end do

  temp2(:,nUnknowns-1,1) = - state0%specificVolume(:,1)                         &
                           * state0%conservedVariables(:,nUnknowns)             &
                           * deltaConservedVariables(:,1)
  temp2(:,nUnknowns-1,1) = temp2(:,nUnknowns-1,1)                               &
                  - sum(state0%velocity * temp2(:,1:nDimensions,1), dim=2)      &
                  + deltaConservedVariables(:,nUnknowns)

  temp2(:,nUnknowns-1,1) = temp2(:,nUnknowns-1,1)                                 &
                             * solverOptions%ratioOfSpecificHeats
  do i = 1, nUnknowns-1
    temp2(:,i,1) = temp2(:,i,1) * state0%specificVolume(:,1)
  end do

  do i = 2, nDimensions
    temp2(:,:,i) = temp2(:,:,1)
  end do
  do i = 1, nDimensions
    call grid%firstDerivative(i)%apply(temp2(:,:,i), grid%localSize)
  end do
  ! (1-2) A * dQ
  allocate(temp1(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))
  temp1 = 0.0_wp

  allocate(localFluxJacobian1(nUnknowns, nUnknowns))
  allocate(localConservedVariables(nUnknowns))
  allocate(localVelocity(nDimensions))
  allocate(localMetricsAlongDirection1(nDimensions))
  if (simulationFlags%viscosityOn) then
     allocate(localFluxJacobian2(nUnknowns, nUnknowns))
     allocate(localStressTensor(nDimensions ** 2))
     allocate(localHeatFlux(nDimensions))
  end if
  do j = 1, grid%nGridPoints
     localConservedVariables = state0%conservedVariables(j,:)
     localVelocity = state0%velocity(j,:)
     if (simulationFlags%viscosityOn) then
        localStressTensor = state0%stressTensor(j,:)
        localHeatFlux = state0%heatFlux(j,:)
     end if
     do i = 1, nDimensions
        localMetricsAlongDirection1 = grid%metrics(j,1+nDimensions*(i-1):nDimensions*i)
        if (simulationFlags%viscosityOn) then
           select case (nDimensions)
           case (1)
              call computeFirstPartialViscousJacobian1D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian2, specificVolume = state0%specificVolume(j,1),           &
                   velocity = localVelocity, temperature = state0%temperature(j,1))
           case (2)
              call computeFirstPartialViscousJacobian2D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian2, specificVolume = state0%specificVolume(j,1),           &
                   velocity = localVelocity, temperature = state0%temperature(j,1))
           case (3)
              call computeFirstPartialViscousJacobian3D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian2, specificVolume = state0%specificVolume(j,1),           &
                   velocity = localVelocity, temperature = state0%temperature(j,1))
           end select
           ! This sign change is for combining with inviscid flux.
           ! For collecting viscous flux, we should not change the sign.
           ! localFluxJacobian1 = - localFluxJacobian2
           localFluxJacobian1 = localFluxJacobian2
        end if

        temp1(j,:,i) = temp1(j,:,i) +                                                         &
                        matmul(localFluxJacobian1, deltaConservedVariables(j,:))
     end do !... i = 1, nDimensions
  end do !... j = 1, grid%nGridPoints

  SAFE_DEALLOCATE(localFluxJacobian1)
  SAFE_DEALLOCATE(localConservedVariables)
  SAFE_DEALLOCATE(localHeatFlux)
  SAFE_DEALLOCATE(localStressTensor)
  SAFE_DEALLOCATE(localFluxJacobian2)

  if (simulationFlags%viscosityOn) then

     allocate(localMetricsAlongDirection2(nDimensions))
     allocate(localFluxJacobian2(nUnknowns - 1, nUnknowns - 1))
     allocate(localLinearizedDiffusion(nUnknowns - 1, nDimensions))

     do k = 1, grid%nGridPoints

        localVelocity = state0%velocity(k,:)
        localLinearizedDiffusion = 0.0_wp

        do i = 1, nDimensions

           localMetricsAlongDirection1 = grid%metrics(k,1+nDimensions*(i-1):nDimensions*i)

           do j = 1, nDimensions

              localMetricsAlongDirection2 = grid%metrics(k,1+nDimensions*(j-1):nDimensions*j)

              select case (nDimensions)
              case (1)
                 call computeSecondPartialViscousJacobian1D(localVelocity,                   &
                      state0%dynamicViscosity(k,1), state0%secondCoefficientOfViscosity(k,1),  &
                      state0%thermalDiffusivity(k,1), grid%jacobian(k,1),                     &
                      localMetricsAlongDirection1(1), localFluxJacobian2)
              case (2)
                 call computeSecondPartialViscousJacobian2D(localVelocity,                   &
                      state0%dynamicViscosity(k,1), state0%secondCoefficientOfViscosity(k,1),  &
                      state0%thermalDiffusivity(k,1), grid%jacobian(k,1),                     &
                      localMetricsAlongDirection1, localMetricsAlongDirection2,              &
                      localFluxJacobian2)
              case (3)
                 call computeSecondPartialViscousJacobian3D(localVelocity,                   &
                      state0%dynamicViscosity(k,1), state0%secondCoefficientOfViscosity(k,1),  &
                      state0%thermalDiffusivity(k,1), grid%jacobian(k,1),                     &
                      localMetricsAlongDirection1, localMetricsAlongDirection2,              &
                      localFluxJacobian2)
              end select !... nDimensions

              localLinearizedDiffusion(:,i) = localLinearizedDiffusion(:,i) +                      &
                                         matmul(localFluxJacobian2,temp2(k,:,j))

           end do !... j = 1, nDimensions

        end do !... i = 1, nDimensions

        do i = 1, nDimensions
          ! This sign change is for combining with inviscid flux.
          ! For collecting viscous flux, we should not change the sign.
           ! temp1(k,2:nUnknowns,i) = temp1(k,2:nUnknowns,i) - localLinearizedDiffusion(:,i)
           temp1(k,2:nUnknowns,i) = temp1(k,2:nUnknowns,i) + localLinearizedDiffusion(:,i)
        end do

     end do !... k = 1, grid%nGridPoints

     SAFE_DEALLOCATE(temp2)

     SAFE_DEALLOCATE(localLinearizedDiffusion)
     SAFE_DEALLOCATE(localFluxJacobian2)
     SAFE_DEALLOCATE(localMetricsAlongDirection2)

  end if !... simulationFlags%viscosityOn
  ! (2) Send viscous fluxes to patch
  if (simulationFlags%viscosityOn) then
    call patchFactories(1)%connect(patch)
    select type (patch)
      class is (t_FarFieldPatch)
        ! patch%viscousFluxes = 0.0_wp
        ! call patch%collect(temp1(:,:,abs(patch%normalDirection)),                               &
        !                     patch%viscousFluxes(:,:,abs(patch%normalDirection)) )
        call patch%collect(temp1, patch%viscousFluxes)
    end select
  end if
  SAFE_DEALLOCATE(temp1)

  ! (2) Add patch penalties
  select type (patch)
    class is (t_FarFieldPatch)
      call patch%updateRhs(LINEARIZED, simulationFlags, solverOptions, grid, state0)
  end select
  linearizedRightHandSide = state0%rightHandSide
  state0%rightHandSide = temp3
  SAFE_DEALLOCATE(temp3)

  ! Randomize adjoint variables.
  allocate(F(grid%nGridPoints, solverOptions%nUnknowns))
  call random_number(F)
  state0%adjointVariables = F
  SAFE_DEALLOCATE(F)

  ! <R^{\dagger}u, \delta v>
  scalar1 = grid%computeInnerProduct(state0%adjointVariables,linearizedRightHandSide)

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

    ! ! Update state for only first partial Jacobian.
    ! call computeDependentVariables(nDimensions, state1%conservedVariables,                    &
    !      solverOptions%ratioOfSpecificHeats, state1%specificVolume(:,1), state1%velocity,       &
    !      state1%pressure(:,1), state1%temperature(:,1))
    ! call computeTransportVariables(state1%temperature(:,1), solverOptions%powerLawExponent,   &
    !      solverOptions%bulkViscosityRatio, solverOptions%ratioOfSpecificHeats,              &
    !      solverOptions%reynoldsNumberInverse, solverOptions%prandtlNumberInverse,           &
    !      state1%dynamicViscosity(:,1), state1%secondCoefficientOfViscosity(:,1),                &
    !      state1%thermalDiffusivity(:,1))
    ! call grid%computeGradient(state0%velocity, state1%stressTensor)
    ! call computeStressTensor(nDimensions, state1%stressTensor, state1%dynamicViscosity(:,1), &
    !      state1%secondCoefficientOfViscosity(:,1))
    !
    ! call grid%computeGradient(state0%temperature(:,1), state1%heatFlux)
    ! do i = 1, nDimensions
    !    state1%heatFlux(:,i) = - state1%thermalDiffusivity(:,1) * state1%heatFlux(:,i)
    ! end do

    ! (2) Compute baseline SAT far-field
    ! (2-1) Cartesian form
    allocate(fluxes2(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))
    if (simulationFlags%viscosityOn .and. simulationFlags%repeatFirstDerivative) then
       call computeCartesianViscousFluxes(nDimensions, state1%velocity,                         &
            state1%stressTensor, state1%heatFlux, fluxes2)
    end if
    ! (2-2) Send viscous fluxes to patch
    if (simulationFlags%viscosityOn) then
      call patchFactories(1)%connect(patch)
      select type (patch)
        class is (t_FarFieldPatch)
          call patch%collect(fluxes2, patch%viscousFluxes)
      end select
    end if
    ! (2-3) Add patch penalty
    state1%rightHandSide = 0.0_wp
    select type (patch)
      class is (t_FarFieldPatch)
        call patch%updateRhs(FORWARD, simulationFlags, solverOptions, grid, state1)
    end select
    SAFE_DEALLOCATE(fluxes2)

    ! (3) <u, \delta R(v)>
    scalar2 = grid%computeInnerProduct(state0%adjointVariables,                             &
                                        state1%rightHandSide - state0%rightHandSide)

    scalarHistory(k) = scalar2/stepSizes(k)
    errorHistory(k) = 0.0_wp
    if( abs(scalar1)>0.0_wp ) errorHistory(k) = abs( (scalar2/stepSizes(k) - scalar1)/scalar1 )

    if (k > 1) then
       convergenceHistory(k-1) = log(errorHistory(k) / errorHistory(k-1)) /              &
            log(stepSizes(k) / stepSizes(k-1))
       if (k > 5) then
         if(sum(convergenceHistory(k-3:k-1))/3.0_wp < 0.0_wp) exit
       end if
    end if
  end do

  if (k > 3) then
     call sort(convergenceHistory(:k-2))
     success = nint(meanTrimmed(convergenceHistory(:k-2))).ge.1
     if ((.not. success) .and. (maxval(errorHistory).le.tolerance_)) then
       write(errorMessage,'(A,E8.3)')                                                     &
       'Error does not follow first-order slope, but is smaller than the tolerance ',     &
       tolerance_
       call writeAndFlush(MPI_COMM_WORLD,output_unit,errorMessage)
       success = .true.
     end if
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

  SAFE_DEALLOCATE(targetViscousFluxes)
  SAFE_DEALLOCATE(linearizedRightHandSide)
  SAFE_DEALLOCATE(deltaConservedVariables)
  SAFE_DEALLOCATE(deltaPrimitiveVariables)

  call patch%cleanup()
  call patchFactories(1)%cleanup()
  SAFE_DEALLOCATE(patchFactories)
  call state0%cleanup()
  call state1%cleanup()
  call grid%cleanup()

end subroutine testLinearizedRelation
