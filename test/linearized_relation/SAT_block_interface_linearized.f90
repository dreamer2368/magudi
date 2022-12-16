#include "config.h"

program SAT_block_interface_linearized

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

     subroutine testLinearizedRelation(identifier, nDimensions, success, direction, tolerance)

       character(len = *), intent(in) :: identifier
       integer, intent(in) :: nDimensions
       logical, intent(out) :: success

       integer, intent(in) :: direction
       real(SCALAR_KIND), intent(in), optional :: tolerance

     end subroutine testLinearizedRelation

  end interface

  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  do nDimensions = 1,3
    do j = 1, 4!... for each discretizationTypes
      success = .true.
      do i = 1, 1 !... test multiple times
        ! Didn't test periodic grid yet!!
        ! isPeriodic = .true.
        ! do direction = 1, nDimensions
        !   call testLinearizedRelation(discretizationTypes(j), nDimensions,           &
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
        do direction = 1,1
          call testLinearizedRelation(discretizationTypes(j), nDimensions,           &
                                   success_, direction)
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

end program SAT_block_interface_linearized

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

subroutine testLinearizedRelation(identifier, nDimensions, success, direction, tolerance)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use StencilOperator_mod, only : t_StencilOperator
  use Region_mod, only : t_Region
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  use PatchDescriptor_mod, only : t_PatchDescriptor
  use Patch_factory, only : t_PatchFactory
  use Patch_mod, only : t_Patch
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  use Region_enum, only : FORWARD, ADJOINT, LINEARIZED
  use BlockInterfacePatch_enum, only : METRICS

  use CNSHelper
  use RhsHelper, only : addInterfaceAdjointPenalty

  ! <<< Internal modules >>>
  use MPIHelper, only : pigeonhole
  use RandomNumber, only : random
  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use Patch_factory, only : updatePatchFactories
  use InterfaceHelper, only : checkFunctionContinuityAtInterfaces, exchangeInterfaceData
  use ErrorHandler

  ! <<< Arguments >>>
  character(len = *), intent(in) :: identifier
  integer, intent(in) :: nDimensions
  logical, intent(out) :: success
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
  type(t_Region) :: region
  type(t_Grid) :: grid(2)
  type(t_State) :: state0(2), state1(2), deltaState(2)
  class(t_Patch), pointer :: patch => null()
  ! type(t_FarFieldPatch) :: patch

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp) :: scalar1, scalar2, tolerance_, hx, scalarHistory(32),                    &
              stepSizes(32), errorHistory(32), convergenceHistory(31)
  integer :: i, j, k, l, gridSize(nDimensions, 2), gridIndex,                         &
             nUnknowns, ierror, normalDirection
  real(SCALAR_KIND), allocatable :: F(:,:), fluxes2(:,:,:),                           &
                                    deltaConservedVariables(:,:), deltaPrimitiveVariables(:,:),&
                                    temp1(:,:,:), localFluxJacobian1(:,:),            &
                                    localConservedVariables(:), localVelocity(:),     &
                                    localMetricsAlongDirection1(:),                   &
                                    localMetricsAlongDirection2(:),                   &
                                    localFluxJacobian2(:,:), localStressTensor(:),    &
                                    localHeatFlux(:), localLinearizedDiffusion(:,:),  &
                                    temp2(:,:,:)
  SCALAR_TYPE :: inviscidPenaltyAmount, viscousPenaltyAmount
  SCALAR_TYPE, dimension(nDimensions) :: h, gridPerturbation
  character(len = STRING_LENGTH) :: errorMessage
  character(len = STRING_LENGTH) :: filename

  tolerance_ = 1.0E-11
  if( present(tolerance) ) tolerance_ = tolerance

  success = .true.

  ! set up simulation flags
  call simulationFlags%initialize()
  simulationFlags%enableController = .false.
  simulationFlags%enableFunctional = .false.
  simulationFlags%enableAdjoint = .true.
  simulationFlags%isDomainCurvilinear = .true.
  simulationFlags%viscosityOn = .true.
  simulationFlags%repeatFirstDerivative = .true. ! this is default value.
  simulationFlags%useTargetState = .true.

  ! randomize grid size
  ! Note that too small grid size will yield null matrix for stencil operators.
  gridSize(:,:) = 1
  do i = 1, nDimensions
     gridSize(i,1) = random(20, 40)
     gridSize(i,2) = gridSize(i,1)
     if (i==direction) gridSize(i,2) = random(20, 40)
  end do

  ! initialize solver option, grid, and states
  call solverOptions%initialize(nDimensions, simulationFlags)
  solverOptions%discretizationType = trim(identifier)
  solverOptions%reynoldsNumberInverse = 1.0_wp / random(5.0_wp,100.0_wp)
  solverOptions%prandtlNumberInverse = 1.0_wp / random(0.1_wp,5.0_wp)
  solverOptions%powerLawExponent = random(0.1_wp, 1.0_wp)
  solverOptions%bulkViscosityRatio = random(0.1_wp, 1.0_wp)

  call region%setup(MPI_COMM_WORLD, gridSize, simulationFlags, solverOptions)

  ! randomize grid 1 coordinates (!!!This is not parallelized.)
  hx = 1.0E-5
  region%grids(1)%coordinates(:,direction) = region%grids(1)%coordinates(:,direction)*1.0E-5
  h = 1.0_wp
  h(direction) = 1.0E-5
  h = h / real(region%grids(1)%globalSize(1:nDimensions)-1,wp)
  do k = 1, region%grids(1)%globalSize(3)
    do j = 1, region%grids(1)%globalSize(2)
      do i = 1, region%grids(1)%globalSize(1)
        call random_number(gridPerturbation)
        gridPerturbation = (2.0_wp * gridPerturbation - 1.0_wp) * 0.05_wp * h

        gridIndex = i + region%grids(1)%globalSize(1)*( j-1 + region%grids(1)%globalSize(2)*( k-1 ) )
        region%grids(1)%coordinates(gridIndex,:) = region%grids(1)%coordinates(gridIndex,:) + gridPerturbation

        select case(direction)
        case(1)
          if (i==region%grids(1)%globalSize(1)) then
            gridIndex = 1 + region%grids(2)%globalSize(1)*( j-1 + region%grids(2)%globalSize(2)*( k-1 ) )
            region%grids(2)%coordinates(gridIndex,:) = region%grids(2)%coordinates(gridIndex,:) + gridPerturbation
          end if
        case(2)
          if (j==region%grids(1)%globalSize(2)) then
            gridIndex = i + region%grids(2)%globalSize(1)*( 0 + region%grids(2)%globalSize(2)*( k-1 ) )
            region%grids(2)%coordinates(gridIndex,:) = region%grids(2)%coordinates(gridIndex,:) + gridPerturbation
          end if
        case(3)
          if (k==region%grids(1)%globalSize(3)) then
            gridIndex = i + region%grids(2)%globalSize(1)*( j-1 + region%grids(2)%globalSize(2)*( 0 ) )
            region%grids(2)%coordinates(gridIndex,:) = region%grids(2)%coordinates(gridIndex,:) + gridPerturbation
          end if
        end select
      end do
    end do
  end do
  ! randomize grid 2 coordinates (!!!This is not parallelized.)
  h = 1.0_wp / real(region%grids(2)%globalSize(1:nDimensions)-1,wp)
  do k = 1, region%grids(2)%globalSize(3)
    do j = 1, region%grids(2)%globalSize(2)
      do i = 1, region%grids(2)%globalSize(1)
        select case(direction)
        case(1)
          if (i==1) cycle
        case(2)
          if (j==1) cycle
        case(3)
          if (k==1) cycle
        end select

        call random_number(gridPerturbation)
        gridPerturbation = (2.0_wp * gridPerturbation - 1.0_wp) * 0.05_wp * h

        gridIndex = i + region%grids(2)%globalSize(1)*( j-1 + region%grids(2)%globalSize(2)*( k-1 ) )
        region%grids(2)%coordinates(gridIndex,:) = region%grids(2)%coordinates(gridIndex,:) + gridPerturbation
      end do
    end do
  end do
  ! shift grid 1 coordinates (!!!This is not parallelized.)
  region%grids(2)%coordinates(:,direction) = region%grids(2)%coordinates(:,direction) + 1.0E-5

  do i = 1, size(region%grids)
     call region%grids(i)%update()
  end do

  ! initialize states
  do i = 1, size(state0)
     call state0(i)%setup(region%grids(i), region%simulationFlags, region%solverOptions)
     call state1(i)%setup(region%grids(i), region%simulationFlags, region%solverOptions)
     call deltaState(i)%setup(region%grids(i), region%simulationFlags, region%solverOptions)
  end do
  call MPI_Barrier(region%comm, ierror)

  ! Randomize conserved variables.
  do k = 1, size(state0)
    do i = 1, region%grids(k)%nGridPoints
       state0(k)%conservedVariables(i,1) = random(0.01_wp, 10.0_wp)
       state0(k)%targetState(i,1) = random(0.01_wp, 10.0_wp)
       do j = 1, nDimensions
          state0(k)%conservedVariables(i,j+1) =                                              &
               state0(k)%conservedVariables(i,1) * random(-10.0_wp, 10.0_wp)
          state0(k)%targetState(i,j+1) =                                              &
                state0(k)%targetState(i,1) * random(-10.0_wp, 10.0_wp)
       end do
       state0(k)%conservedVariables(i, nDimensions + 2) = state0(k)%conservedVariables(i,1) *    &
            random(0.01_wp, 10.0_wp) / region%solverOptions%ratioOfSpecificHeats +              &
            0.5_wp / state0(k)%conservedVariables(i,1) *                                     &
            sum(state0(k)%conservedVariables(i,2:nDimensions+1) ** 2)
       state0(k)%targetState(i, nDimensions + 2) = state0(k)%targetState(i,1) *   &
             random(0.01_wp, 10.0_wp) / solverOptions%ratioOfSpecificHeats +               &
             0.5_wp / state0(k)%targetState(i,1) *                                     &
             sum(state0(k)%targetState(i,2:nDimensions+1) ** 2)
    end do
  end do
  do i = 1, size(state0)
    assert(all(state0(i)%conservedVariables(:,1) > 0.0_wp))
    assert(all(state0(i)%targetState(:,1) > 0.0_wp))
  end do

  ! Compute dependent variables.
  do i = 1, size(state0)
    call state0(i)%update(region%grids(i),region%simulationFlags,region%solverOptions)
  end do

  ! Randomize adjoint variables.
  do i = 1, size(state0)
    allocate(F(region%grids(i)%nGridPoints, region%solverOptions%nUnknowns))
    call random_number(F)
    state0(i)%adjointVariables = F
    SAFE_DEALLOCATE(F)
  end do

  ! Randomize delta conserved variables.
  do k = 1, size(state0)
    allocate(deltaPrimitiveVariables(region%grids(k)%nGridPoints, region%solverOptions%nUnknowns))
    do i = 1, region%grids(k)%nGridPoints
      do j = 1, nDimensions + 2
         deltaPrimitiveVariables(i,j) = random(-1.0_wp, 1.0_wp)
      end do
    end do
    deltaState(k)%conservedVariables(:,1) = deltaPrimitiveVariables(:,1)
    do j = 1, nDimensions
       deltaState(k)%conservedVariables(:,j+1) = state0(k)%conservedVariables(:,j+1) /                     &
            state0(k)%conservedVariables(:,1) * deltaPrimitiveVariables(:,1) +                    &
            state0(k)%conservedVariables(:,1) * deltaPrimitiveVariables(:,j+1)
    end do
    deltaState(k)%conservedVariables(:,nDimensions+2) = state0(k)%conservedVariables(:,nDimensions+2) /    &
         state0(k)%conservedVariables(:,1) * deltaPrimitiveVariables(:,1) +                       &
         sum(state0(k)%conservedVariables(:,2:nDimensions+1) *                                    &
         deltaPrimitiveVariables(:,2:nDimensions+1), dim = 2) +                          &
         state0(k)%conservedVariables(:,1) / region%solverOptions%ratioOfSpecificHeats *                               &
         deltaPrimitiveVariables(:,nDimensions+2)
    SAFE_DEALLOCATE(deltaPrimitiveVariables)
  end do

  ! Set region to state0
  do i = 1, size(state0)
    region%states(i)%conservedVariables = state0(i)%conservedVariables
    region%states(i)%targetState = state0(i)%targetState
    call region%states(i)%update(region%grids(i), region%simulationFlags, region%solverOptions)

    region%states(i)%adjointVariables = deltaState(i)%conservedVariables
  end do

  ! Parse options from the input file.
  filename = PROJECT_NAME // ".inp"
  call parseInputFile(filename)

  ! Write out some useful information.
  ! call region%reportGridDiagnostics()

  ! Setup boundary conditions.
  call getRequiredOption("boundary_condition_file", filename)
  call region%setupBoundaryConditions(filename)

  ! Check continuity at block interfaces.
  ! call checkFunctionContinuityAtInterfaces(region, epsilon(0.0_wp))

  inviscidPenaltyAmount = random(0.01_wp,10.0_wp)
  viscousPenaltyAmount = random(0.01_wp,10.0_wp)
  do i = 1, size(region%patchFactories)
     call region%patchFactories(i)%connect(patch)
     if (.not. associated(patch)) cycle
     do j = 1, size(region%states)
       if (patch%gridIndex /= region%grids(j)%index) cycle
       select type (patch)
       class is (t_BlockInterfacePatch)
         patch%inviscidPenaltyAmount = sign(inviscidPenaltyAmount,                              &
                                        real(patch%normalDirection, wp))                        &
                                        /region%grids(j)%firstDerivative(abs(patch%normalDirection))%normBoundary(1)
         patch%viscousPenaltyAmount = sign(viscousPenaltyAmount,                             &
                                        real(patch%normalDirection, wp))                      &
                                        /region%grids(j)%firstDerivative(abs(patch%normalDirection))%normBoundary(1)
       end select
     end do
  end do

  ! Update patches.
  do i = 1, size(region%grids)
     call updatePatchFactories(region%patchFactories, region%simulationFlags,                &
          region%solverOptions, region%grids(i), region%states(i))
  end do

  ! Exchange jacobian at block interfaces.
   do i = 1, size(region%patchFactories)
      call region%patchFactories(i)%connect(patch)
      if (.not. associated(patch)) cycle
      do j = 1, size(region%states)
         if (patch%gridIndex /= region%grids(j)%index) cycle
         select type (patch)
         class is (t_BlockInterfacePatch)
            call patch%collectInterfaceData(METRICS, simulationFlags,                    &
                 solverOptions, region%grids(j), region%states(j))
         end select
      end do
   end do
  call exchangeInterfaceData(region)
  ! (4) Disperse received data at block interfaces.
   do i = 1, size(region%patchFactories)
      call region%patchFactories(i)%connect(patch)
      if (.not. associated(patch)) cycle
      do j = 1, size(region%states)
         if (patch%gridIndex /= region%grids(j)%index) cycle
         select type (patch)
         class is (t_BlockInterfacePatch)
            call patch%disperseInterfaceData(METRICS, region%simulationFlags,                   &
                 region%solverOptions)
         end select
      end do
   end do

  ! Compute baseline SAT block interface
  do k = 1, 2
    ! (1) Cartesian form
    allocate(fluxes2(region%grids(k)%nGridPoints, solverOptions%nUnknowns, nDimensions))
    if (simulationFlags%viscosityOn .and. simulationFlags%repeatFirstDerivative) then
       call computeCartesianViscousFluxes(nDimensions, region%states(k)%velocity,                         &
            region%states(k)%stressTensor, region%states(k)%heatFlux, fluxes2)
    end if
    ! (2) Send viscous fluxes to patch
    if (simulationFlags%viscosityOn) then
      do i = 1, size(region%patchFactories)
        call region%patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        if (patch%gridIndex /= region%grids(k)%index) cycle
        select type (patch)
        class is (t_BlockInterfacePatch)
           call patch%collect(fluxes2, patch%cartesianViscousFluxesL)
        end select
      end do
    end if
    SAFE_DEALLOCATE(fluxes2)
  end do
  ! (3) Exchange data at block interfaces.
   do i = 1, size(region%patchFactories)
      call region%patchFactories(i)%connect(patch)
      if (.not. associated(patch)) cycle
      do j = 1, size(region%states)
         if (patch%gridIndex /= region%grids(j)%index) cycle
         select type (patch)
         class is (t_BlockInterfacePatch)
            call patch%collectInterfaceData(FORWARD, simulationFlags,                    &
                 solverOptions, region%grids(j), region%states(j))
         end select
      end do
   end do
  call exchangeInterfaceData(region)
  ! (4) Disperse received data at block interfaces.
   do i = 1, size(region%patchFactories)
      call region%patchFactories(i)%connect(patch)
      if (.not. associated(patch)) cycle
      do j = 1, size(region%states)
         if (patch%gridIndex /= region%grids(j)%index) cycle
         select type (patch)
         class is (t_BlockInterfacePatch)
            call patch%disperseInterfaceData(FORWARD, region%simulationFlags,                   &
                 region%solverOptions)
         end select
      end do
   end do
  ! (5) Add patch penalty
  do i = 1, size(region%states)
    region%states(i)%rightHandSide = 0.0_wp
  end do
  do i = 1, size(region%patchFactories)
     call region%patchFactories(i)%connect(patch)
     if (.not. associated(patch)) cycle
     do j = 1, size(region%states)
        if (patch%gridIndex /= region%grids(j)%index) cycle
        call patch%updateRhs(FORWARD, simulationFlags, solverOptions,              &
             region%grids(j), region%states(j))
     end do
  end do

  do i = 1, size(state0)
    state0(i)%rightHandSide = region%states(i)%rightHandSide
  end do

  ! Compute adjoint rhs for SAT block interface
  nUnknowns = solverOptions%nUnknowns
  region%states(1)%rightHandSide = 0.0_wp
  region%states(2)%rightHandSide = 0.0_wp
  ! (1) Collect delta viscous fluxes on far-field patches: This is the same as linearized viscous flux.
  if (simulationFlags%viscosityOn) then
    do l = 1, size(region%grids)
      ! (1-1) D * C * dQ
      allocate(temp2(region%grids(l)%nGridPoints, solverOptions%nUnknowns-1, nDimensions))
      temp2 = 0.0_wp
      do i = 1, nDimensions
        temp2(:,i,1) = - region%states(l)%velocity(:,i) * region%states(l)%adjointVariables(:,1)        &
                      + region%states(l)%adjointVariables(:,i+1)
      end do

      temp2(:,nUnknowns-1,1) = - region%states(l)%specificVolume(:,1)                         &
                               * region%states(l)%conservedVariables(:,nUnknowns)             &
                               * region%states(l)%adjointVariables(:,1)
      temp2(:,nUnknowns-1,1) = temp2(:,nUnknowns-1,1)                               &
                      - sum(region%states(l)%velocity * temp2(:,1:nDimensions,1), dim=2)      &
                      + region%states(l)%adjointVariables(:,nUnknowns)

      temp2(:,nUnknowns-1,1) = temp2(:,nUnknowns-1,1)                                 &
                                 * solverOptions%ratioOfSpecificHeats
      do i = 1, nUnknowns-1
        temp2(:,i,1) = temp2(:,i,1) * region%states(l)%specificVolume(:,1)
      end do

      do i = 2, nDimensions
        temp2(:,:,i) = temp2(:,:,1)
      end do
      do i = 1, nDimensions
        call region%grids(l)%firstDerivative(i)%apply(temp2(:,:,i), region%grids(l)%localSize)
      end do

      ! (1-2) A * dQ
      allocate(temp1(region%grids(l)%nGridPoints, solverOptions%nUnknowns, nDimensions))
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
      do j = 1, region%grids(l)%nGridPoints
         localConservedVariables = region%states(l)%conservedVariables(j,:)
         localVelocity = region%states(l)%velocity(j,:)
         if (simulationFlags%viscosityOn) then
            localStressTensor = region%states(l)%stressTensor(j,:)
            localHeatFlux = region%states(l)%heatFlux(j,:)
         end if
         do i = 1, nDimensions
            localMetricsAlongDirection1 = region%grids(l)%metrics(j,1+nDimensions*(i-1):nDimensions*i)
            if (simulationFlags%viscosityOn) then
               select case (nDimensions)
               case (1)
                  call computeFirstPartialViscousJacobian1D(localConservedVariables,             &
                       localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                       solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                       localFluxJacobian2, specificVolume = region%states(l)%specificVolume(j,1),           &
                       velocity = localVelocity, temperature = region%states(l)%temperature(j,1))
               case (2)
                  call computeFirstPartialViscousJacobian2D(localConservedVariables,             &
                       localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                       solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                       localFluxJacobian2, specificVolume = region%states(l)%specificVolume(j,1),           &
                       velocity = localVelocity, temperature = region%states(l)%temperature(j,1))
               case (3)
                  call computeFirstPartialViscousJacobian3D(localConservedVariables,             &
                       localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                       solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                       localFluxJacobian2, specificVolume = region%states(l)%specificVolume(j,1),           &
                       velocity = localVelocity, temperature = region%states(l)%temperature(j,1))
               end select
               ! This sign change is for combining with inviscid flux.
               ! For collecting viscous flux, we should not change the sign.
               ! localFluxJacobian1 = - localFluxJacobian2
               localFluxJacobian1 = localFluxJacobian2
            end if

            temp1(j,:,i) = temp1(j,:,i) +                                                         &
                            matmul(localFluxJacobian1, region%states(l)%adjointVariables(j,:))
         end do !... i = 1, nDimensions
      end do !... j = 1, grid%nGridPoints

      SAFE_DEALLOCATE(localFluxJacobian1)
      SAFE_DEALLOCATE(localConservedVariables)
      SAFE_DEALLOCATE(localHeatFlux)
      SAFE_DEALLOCATE(localStressTensor)
      SAFE_DEALLOCATE(localFluxJacobian2)

      allocate(localMetricsAlongDirection2(nDimensions))
      allocate(localFluxJacobian2(nUnknowns - 1, nUnknowns - 1))
      allocate(localLinearizedDiffusion(nUnknowns - 1, nDimensions))

      do k = 1, region%grids(l)%nGridPoints

         localVelocity = region%states(l)%velocity(k,:)
         localLinearizedDiffusion = 0.0_wp

         do i = 1, nDimensions

            localMetricsAlongDirection1 = region%grids(l)%metrics(k,1+nDimensions*(i-1):nDimensions*i)

            do j = 1, nDimensions

               localMetricsAlongDirection2 = region%grids(l)%metrics(k,1+nDimensions*(j-1):nDimensions*j)

               select case (nDimensions)
               case (1)
                  call computeSecondPartialViscousJacobian1D(localVelocity,                   &
                       region%states(l)%dynamicViscosity(k,1), region%states(l)%secondCoefficientOfViscosity(k,1),  &
                       region%states(l)%thermalDiffusivity(k,1), region%grids(l)%jacobian(k,1),                     &
                       localMetricsAlongDirection1(1), localFluxJacobian2)
               case (2)
                  call computeSecondPartialViscousJacobian2D(localVelocity,                   &
                       region%states(l)%dynamicViscosity(k,1), region%states(l)%secondCoefficientOfViscosity(k,1),  &
                       region%states(l)%thermalDiffusivity(k,1), region%grids(l)%jacobian(k,1),                     &
                       localMetricsAlongDirection1, localMetricsAlongDirection2,              &
                       localFluxJacobian2)
               case (3)
                  call computeSecondPartialViscousJacobian3D(localVelocity,                   &
                       region%states(l)%dynamicViscosity(k,1), region%states(l)%secondCoefficientOfViscosity(k,1),  &
                       region%states(l)%thermalDiffusivity(k,1), region%grids(l)%jacobian(k,1),                     &
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

      SAFE_DEALLOCATE(localVelocity)
      SAFE_DEALLOCATE(localMetricsAlongDirection1)
      SAFE_DEALLOCATE(localLinearizedDiffusion)
      SAFE_DEALLOCATE(localFluxJacobian2)
      SAFE_DEALLOCATE(localMetricsAlongDirection2)
      SAFE_DEALLOCATE(temp2)

      ! Send viscous fluxes to patches that will use it.
       do i = 1, size(region%patchFactories)
          call region%patchFactories(i)%connect(patch)
          if (.not. associated(patch)) cycle
          if (patch%gridIndex /= region%grids(l)%index .or. patch%nPatchPoints <= 0) cycle

          select type (patch)
          class is (t_BlockInterfacePatch)
             normalDirection = abs(patch%normalDirection)
             call patch%collect(temp1(:,:,normalDirection), patch%viscousFluxesL)
          end select

       end do
       SAFE_DEALLOCATE(temp1)
    end do !... l = 1, size(region%grids)
  end if

  ! (1) Exchange data at block interfaces.
  do i = 1, size(region%patchFactories)
    call region%patchFactories(i)%connect(patch)
    if (.not. associated(patch)) cycle
    do j = 1, size(region%states)
       if (patch%gridIndex /= region%grids(j)%index) cycle
       select type (patch)
       class is (t_BlockInterfacePatch)
          call patch%collectInterfaceData(LINEARIZED, simulationFlags,                    &
               solverOptions, region%grids(j), region%states(j))
       end select
    end do
  end do
  call exchangeInterfaceData(region)
  ! (2) Disperse received data at block interfaces.
  do i = 1, size(region%patchFactories)
    call region%patchFactories(i)%connect(patch)
    if (.not. associated(patch)) cycle
    do j = 1, size(region%states)
       if (patch%gridIndex /= region%grids(j)%index) cycle
       select type (patch)
       class is (t_BlockInterfacePatch)
          call patch%disperseInterfaceData(LINEARIZED, simulationFlags,                   &
               solverOptions)
       end select
    end do
  end do
  ! (3) Add patch penalties.
  do i = 1, size(region%patchFactories)
    call region%patchFactories(i)%connect(patch)
    if (.not. associated(patch)) cycle
    do j = 1, size(region%states)
       if (patch%gridIndex /= region%grids(j)%index) cycle
       call patch%updateRhs(LINEARIZED, region%simulationFlags, region%solverOptions,              &
            region%grids(j), region%states(j))
    end do
  end do

  ! <R^{\dagger}u, \delta v>
  scalar1 = 0.0_wp
  do i = 1, size(region%states)
    scalar1 = scalar1 + region%grids(i)%computeInnerProduct(state0(i)%adjointVariables,           &
                                                            region%states(i)%rightHandSide)
  end do
  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE, scalar1, 1,                          &
       SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(scalar1, 1, SCALAR_TYPE_MPI,                             &
          0, region%grids(i)%comm, ierror)
  end do

  ! <u, \delta R(v)>
  ! Prepare step sizes
  stepSizes(1) = 0.01_wp
  do k = 2, size(stepSizes)
     stepSizes(k) = stepSizes(k-1) * 10.0_wp**(-0.25_wp)
  end do
  do k = 1, size(stepSizes)
    !(1) finite difference on conserved variables
    do i = 1, size(state0)
      region%states(i)%conservedVariables = state0(i)%conservedVariables + stepSizes(k) * deltaState(i)%conservedVariables
      assert(all(region%states(i)%conservedVariables(:,1) > 0.0_wp))

      ! Compute dependent variables.
      call region%states(i)%update(region%grids(i),region%simulationFlags,region%solverOptions)
      assert(all(region%states(i)%specificVolume(:,1) > 0.0_wp))
      assert(all(region%states(i)%temperature(:,1) > 0.0_wp))

      ! ! Update state for only first partial Jacobian.
      ! call computeDependentVariables(nDimensions, region%states(i)%conservedVariables,                    &
      !      solverOptions%ratioOfSpecificHeats, region%states(i)%specificVolume(:,1), region%states(i)%velocity,       &
      !      region%states(i)%pressure(:,1), region%states(i)%temperature(:,1))
      ! call computeTransportVariables(region%states(i)%temperature(:,1), solverOptions%powerLawExponent,   &
      !      solverOptions%bulkViscosityRatio, solverOptions%ratioOfSpecificHeats,              &
      !      solverOptions%reynoldsNumberInverse, solverOptions%prandtlNumberInverse,           &
      !      region%states(i)%dynamicViscosity(:,1), region%states(i)%secondCoefficientOfViscosity(:,1),                &
      !      region%states(i)%thermalDiffusivity(:,1))
      ! call region%grids(i)%computeGradient(state0(i)%velocity, region%states(i)%stressTensor)
      ! call computeStressTensor(nDimensions, region%states(i)%stressTensor, region%states(i)%dynamicViscosity(:,1), &
      !      region%states(i)%secondCoefficientOfViscosity(:,1))
      !
      ! call region%grids(i)%computeGradient(state0(i)%temperature(:,1), region%states(i)%heatFlux)
      ! do j = 1, nDimensions
      !    region%states(i)%heatFlux(:,j) = - region%states(i)%thermalDiffusivity(:,1) * region%states(i)%heatFlux(:,j)
      ! end do

    end do

    ! (2) Compute baseline SAT block interface
    do j = 1, 2
      ! (2-1) Cartesian form
      allocate(fluxes2(region%grids(j)%nGridPoints, solverOptions%nUnknowns, nDimensions))
      if (simulationFlags%viscosityOn .and. simulationFlags%repeatFirstDerivative) then
         call computeCartesianViscousFluxes(nDimensions, region%states(j)%velocity,                         &
              region%states(j)%stressTensor, region%states(j)%heatFlux, fluxes2)
      end if
      ! (2) Send viscous fluxes to patch
      if (simulationFlags%viscosityOn) then
        do i = 1, size(region%patchFactories)
          call region%patchFactories(i)%connect(patch)
          if (.not. associated(patch)) cycle
          if (patch%gridIndex /= region%grids(j)%index) cycle
          select type (patch)
          class is (t_BlockInterfacePatch)
             call patch%collect(fluxes2, patch%cartesianViscousFluxesL)
          end select
        end do
      end if
      SAFE_DEALLOCATE(fluxes2)
    end do
    ! (3) Exchange data at block interfaces.
     do i = 1, size(region%patchFactories)
        call region%patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        do j = 1, size(region%states)
           if (patch%gridIndex /= region%grids(j)%index) cycle
           select type (patch)
           class is (t_BlockInterfacePatch)
              call patch%collectInterfaceData(FORWARD, simulationFlags,                    &
                   solverOptions, region%grids(j), region%states(j))
           end select
        end do
     end do
    call exchangeInterfaceData(region)
    ! (4) Disperse received data at block interfaces.
     do i = 1, size(region%patchFactories)
        call region%patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        do j = 1, size(region%states)
           if (patch%gridIndex /= region%grids(j)%index) cycle
           select type (patch)
           class is (t_BlockInterfacePatch)
              call patch%disperseInterfaceData(FORWARD, region%simulationFlags,                   &
                   region%solverOptions)
           end select
        end do
     end do
    ! (5) Add patch penalty
    do i = 1, size(region%states)
      region%states(i)%rightHandSide = 0.0_wp
    end do
    do i = 1, size(region%patchFactories)
       call region%patchFactories(i)%connect(patch)
       if (.not. associated(patch)) cycle
       do j = 1, size(region%states)
          if (patch%gridIndex /= region%grids(j)%index) cycle
          call patch%updateRhs(FORWARD, simulationFlags, solverOptions,              &
               region%grids(j), region%states(j))
       end do
    end do

    ! (3) <u, \delta R(v)>
    scalar2 = 0.0_wp
    do i = 1, size(state0)
      scalar2 = scalar2 + region%grids(i)%computeInnerProduct(state0(i)%adjointVariables,                             &
                                                    region%states(i)%rightHandSide - state0(i)%rightHandSide)
    end do
    if (region%commGridMasters /= MPI_COMM_NULL)                                               &
         call MPI_Allreduce(MPI_IN_PLACE, scalar2, 1,                          &
         SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

    do i = 1, size(region%grids)
       call MPI_Bcast(scalar2, 1, SCALAR_TYPE_MPI,                             &
            0, region%grids(i)%comm, ierror)
    end do

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

  SAFE_DEALLOCATE(deltaConservedVariables)
  SAFE_DEALLOCATE(deltaPrimitiveVariables)

  ! call patch%cleanup()
  ! call patchFactories(1)%cleanup()
  ! SAFE_DEALLOCATE(patchFactories)
  do i = 1, 2
    call state0(i)%cleanup()
    call state1(i)%cleanup()
    call deltaState(i)%cleanup()
    call grid(i)%cleanup()
  end do
  call region%cleanup()

end subroutine testLinearizedRelation
