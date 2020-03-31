#include "config.h"

module RhsHelperImpl

  implicit none
  public

contains

  subroutine addDissipation(mode, simulationFlags, solverOptions, grid, state)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod, only : t_State
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Enumerations >>>
    use Region_enum, only : FORWARD, ADJOINT, LINEARIZED

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    integer, intent(in) :: mode
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid), intent(in) :: grid
    class(t_State) :: state

    ! <<< Local variables >>>
    integer :: i, j, nDimensions, nUnknowns
    SCALAR_TYPE, allocatable :: dissipationTerm(:,:)
    real(SCALAR_KIND) :: dissipationAmount

    if (.not. simulationFlags%dissipationOn) return

    call startTiming("addDissipation")

    assert_key(mode, (FORWARD, ADJOINT))

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    nUnknowns = solverOptions%nUnknowns

    allocate(dissipationTerm(grid%nGridPoints, nUnknowns))

    select case (mode)
    case (FORWARD)
       dissipationAmount = + solverOptions%dissipationAmount
    case (ADJOINT)
       dissipationAmount = - solverOptions%dissipationAmount
    case (LINEARIZED)
       dissipationAmount = + solverOptions%dissipationAmount
    end select

    do i = 1, nDimensions

       select case (mode)
       case (FORWARD)
          dissipationTerm = state%conservedVariables
       case (ADJOINT)
          dissipationTerm = state%adjointVariables
       case (LINEARIZED)
          dissipationTerm = state%adjointVariables
       end select

       call grid%dissipation(i)%apply(dissipationTerm, grid%localSize)

       if (.not. simulationFlags%compositeDissipation) then
          do j = 1, nUnknowns
             dissipationTerm(:,j) = - grid%arcLengths(:,i) * dissipationTerm(:,j)
          end do
          call grid%dissipationTranspose(i)%apply(dissipationTerm, grid%localSize)
          call grid%firstDerivative(i)%applyNormInverse(dissipationTerm, grid%localSize)
       end if

       state%rightHandSide = state%rightHandSide + dissipationAmount * dissipationTerm

    end do

    SAFE_DEALLOCATE(dissipationTerm)

    call endTiming("addDissipation")

  end subroutine addDissipation

  subroutine addFarFieldAdjointPenalty(simulationFlags, solverOptions,                       &
       grid, state, patchFactories)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod, only : t_State
    use Patch_mod, only : t_Patch
    use Patch_factory, only : t_PatchFactory
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags
    use FarFieldPatch_mod, only : t_FarFieldPatch

    ! <<< Internal modules >>>
    use CNSHelper
    use Patch_factory, only : queryPatchTypeExists

    implicit none

    ! <<< Arguments >>>
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid), intent(in) :: grid
    class(t_State) :: state
    type(t_PatchFactory), allocatable :: patchFactories(:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    logical :: farFieldPatchesExist
    integer :: i, j, k, l, iPatchFactory, nDimensions, nUnknowns,                            &
         direction, gridIndex, patchIndex, ierror
    SCALAR_TYPE, allocatable :: localConservedVariables(:), localVelocity(:),                &
         localMetricsAlongDirection1(:), localMetricsAlongDirection2(:),                     &
         localViscousFluxJacobian(:,:), temp1(:,:,:), temp2(:,:)
    class(t_Patch), pointer :: patch => null()

    farFieldPatchesExist = queryPatchTypeExists(patchFactories,                              &
         'SAT_FAR_FIELD', grid%index)
    call MPI_Allreduce(MPI_IN_PLACE, farFieldPatchesExist, 1, MPI_LOGICAL,                   &
         MPI_LOR, grid%comm, ierror) !... reduce across grid-level processes.
    if (.not. farFieldPatchesExist) return

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns >= nDimensions + 2)

    allocate(localConservedVariables(nUnknowns))
    allocate(localVelocity(nDimensions))
    allocate(localMetricsAlongDirection1(nDimensions))
    allocate(localMetricsAlongDirection2(nDimensions))
    allocate(localViscousFluxJacobian(nUnknowns - 1, nUnknowns - 1))
    allocate(temp1(grid%nGridPoints, nUnknowns - 1, nDimensions))
    allocate(temp2(grid%nGridPoints, nUnknowns - 1))

    temp1 = 0.0_wp

    if (allocated(patchFactories)) then
       do iPatchFactory = 1, size(patchFactories)
          call patchFactories(iPatchFactory)%connect(patch)
          if (.not. associated(patch)) cycle
          if (patch%gridIndex /= grid%index .or. patch%nPatchPoints <= 0) cycle

          assert(all(grid%offset == patch%gridOffset))
          assert(all(grid%localSize == patch%gridLocalSize))

          direction = abs(patch%normalDirection)

          select type (patch)
             class is (t_FarFieldPatch)
             assert(direction >= 1 .and. direction <= nDimensions)

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

                      localConservedVariables = state%conservedVariables(gridIndex,:)
                      localVelocity = state%velocity(gridIndex,:)
                      localMetricsAlongDirection1 = grid%metrics(gridIndex,                  &
                           1 + nDimensions * (direction - 1) : nDimensions * direction)

                      do l = 1, nDimensions

                         localMetricsAlongDirection2 =                                       &
                              grid%metrics(gridIndex,1+nDimensions*(l-1):nDimensions*l)

                         select case (nDimensions)
                         case (1)
                            call computeSecondPartialViscousJacobian1D(localVelocity,        &
                                 state%dynamicViscosity(gridIndex,1),                        &
                                 state%secondCoefficientOfViscosity(gridIndex,1),            &
                                 state%thermalDiffusivity(gridIndex,1),                      &
                                 grid%jacobian(gridIndex,1), localMetricsAlongDirection1(1), &
                                 localViscousFluxJacobian)
                         case (2)
                            call computeSecondPartialViscousJacobian2D(localVelocity,        &
                                 state%dynamicViscosity(gridIndex,1),                        &
                                 state%secondCoefficientOfViscosity(gridIndex,1),            &
                                 state%thermalDiffusivity(gridIndex,1),                      &
                                 grid%jacobian(gridIndex,1), localMetricsAlongDirection1,    &
                                 localMetricsAlongDirection2, localViscousFluxJacobian)
                         case (3)
                            call computeSecondPartialViscousJacobian3D(localVelocity,        &
                                 state%dynamicViscosity(gridIndex,1),                        &
                                 state%secondCoefficientOfViscosity(gridIndex,1),            &
                                 state%thermalDiffusivity(gridIndex,1),                      &
                                 grid%jacobian(gridIndex,1), localMetricsAlongDirection1,    &
                                 localMetricsAlongDirection2, localViscousFluxJacobian)
                         end select !... nDimensions

                         temp1(gridIndex,:,l) = temp1(gridIndex,:,l) -                       &
                              patch%viscousPenaltyAmount *                                   &
                              matmul(transpose(localViscousFluxJacobian),                    &
                              state%adjointVariables(gridIndex,2:nUnknowns))

                      end do !... l = 1, nDimensions

                   end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
                end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
             end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

          end select !... type (patch)

       end do !... iPatchFactory = 1, size(patchFactories)
    end if !... allocated(patchFactories)

    do i = 1, nDimensions
       call grid%adjointFirstDerivative(i)%apply(temp1(:,:,i), grid%localSize)
    end do
    temp2 = sum(temp1, dim = 3)

    temp2(:,nDimensions+1) = solverOptions%ratioOfSpecificHeats *                            &
         state%specificVolume(:,1) * temp2(:,nDimensions+1)
    do i = 1, nDimensions
       temp2(:,i) = state%specificVolume(:,1) * temp2(:,i) -                                 &
            state%velocity(:,i) * temp2(:,nDimensions+1)
    end do

    state%rightHandSide(:,2:nUnknowns) = state%rightHandSide(:,2:nUnknowns) + temp2
    state%rightHandSide(:,1) = state%rightHandSide(:,1) -                                    &
         state%specificVolume(:,1) * state%conservedVariables(:,nDimensions+2) *             &
         temp2(:,nDimensions+1) - sum(state%velocity * temp2(:,1:nDimensions), dim = 2)

    SAFE_DEALLOCATE(temp2)
    SAFE_DEALLOCATE(temp1)
    SAFE_DEALLOCATE(localViscousFluxJacobian)
    SAFE_DEALLOCATE(localMetricsAlongDirection2)
    SAFE_DEALLOCATE(localMetricsAlongDirection1)
    SAFE_DEALLOCATE(localVelocity)
    SAFE_DEALLOCATE(localConservedVariables)

  end subroutine addFarFieldAdjointPenalty

end module RhsHelperImpl

subroutine computeRhsForward(simulationFlags, solverOptions, grid, state, patchFactories)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use Patch_mod, only : t_Patch
  use State_mod, only : t_State
  use Patch_factory, only : t_PatchFactory
  use FarFieldPatch_mod, only : t_FarFieldPatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD

  ! <<< Private members >>>
  use RhsHelperImpl, only : addDissipation

  ! <<< Internal modules >>>
  use CNSHelper
  use Patch_factory, only : queryPatchTypeExists
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid) :: grid
  class(t_State) :: state
  type(t_PatchFactory), allocatable :: patchFactories(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions
  SCALAR_TYPE, allocatable :: fluxes1(:,:,:), fluxes2(:,:,:)
  class(t_Patch), pointer :: patch => null()

  call startTiming("computeRhsForward")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  allocate(fluxes1(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))
  allocate(fluxes2(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))

  state%rightHandSide = 0.0_wp

  ! Compute Cartesian form of inviscid fluxes.
  call computeCartesianInviscidFluxes(nDimensions, state%conservedVariables,                 &
       state%velocity, state%pressure(:,1), fluxes1)

  ! Compute Cartesian form of viscous fluxes if viscous terms are included and computed using
  ! repeated first derivatives.
  if (simulationFlags%viscosityOn .and. simulationFlags%repeatFirstDerivative) then
     call computeCartesianViscousFluxes(nDimensions, state%velocity,                         &
          state%stressTensor, state%heatFlux, fluxes2)
     fluxes1 = fluxes1 - fluxes2 !... Cartesian form of total fluxes.
  end if

  ! Send viscous fluxes to patches that will use it.
  if (simulationFlags%viscosityOn .and. allocated(patchFactories)) then
     do i = 1, size(patchFactories)
        call patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        if (patch%gridIndex /= grid%index .or. patch%nPatchPoints <= 0) cycle

        select type (patch)
        class is (t_FarFieldPatch)
           call patch%collect(fluxes2, patch%viscousFluxes)
        class is (t_BlockInterfacePatch)
           call patch%collect(fluxes2, patch%cartesianViscousFluxesL)
        end select

     end do
  end if

  ! Transform fluxes from Cartesian to contravariant form: `fluxes1` has the Cartesian form of
  ! total fluxes... upon return, `fluxes2` has the contravariant form.
  call transformFluxes(nDimensions, fluxes1, grid%metrics, fluxes2, grid%isCurvilinear)

  SAFE_DEALLOCATE(fluxes1) !... no longer needed.

  ! Take derivatives of fluxes.
  do i = 1, nDimensions
     call grid%firstDerivative(i)%apply(fluxes2(:,:,i), grid%localSize)
  end do
  state%rightHandSide = state%rightHandSide - sum(fluxes2, dim = 3) !... update RHS.

  SAFE_DEALLOCATE(fluxes2) !... no longer needed

  ! Add dissipation if required.
  if (simulationFlags%dissipationOn)                                                         &
       call addDissipation(FORWARD, simulationFlags, solverOptions, grid, state)

  call endTiming("computeRhsForward")

end subroutine computeRhsForward

subroutine computeRhsAdjoint(simulationFlags, solverOptions, grid, state, patchFactories)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use Patch_factory, only : t_PatchFactory
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  ! <<< Enumerations >>>
  use Region_enum, only : ADJOINT

  ! <<< Private members >>>
  use RhsHelperImpl, only : addDissipation, addFarFieldAdjointPenalty

  ! <<< Internal modules >>>
  use CNSHelper
  use Patch_factory, only : queryPatchTypeExists
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid) :: grid
  class(t_State) :: state
  type(t_PatchFactory), allocatable :: patchFactories(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, nUnknowns
  SCALAR_TYPE, allocatable :: temp1(:,:,:), temp2(:,:),                                      &
       localFluxJacobian1(:,:), localFluxJacobian2(:,:), localConservedVariables(:),         &
       localVelocity(:), localMetricsAlongDirection1(:), localMetricsAlongDirection2(:),     &
       localStressTensor(:), localHeatFlux(:), localAdjointDiffusion(:,:)

  call startTiming("computeRhsAdjoint")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns >= nDimensions + 2)

  allocate(temp1(grid%nGridPoints, nUnknowns, nDimensions))

  state%rightHandSide = 0.0_wp

  ! Partial derivatives of adjoint variables w.r.t. *computational* coordinates.
  do i = 1, nDimensions
     temp1(:,:,i) = state%adjointVariables
     call grid%adjointFirstDerivative(i)%apply(temp1(:,:,i), grid%localSize)
  end do

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

     localConservedVariables = state%conservedVariables(j,:)
     localVelocity = state%velocity(j,:)
     if (simulationFlags%viscosityOn) then
        localStressTensor = state%stressTensor(j,:)
        localHeatFlux = state%heatFlux(j,:)
     end if

     do i = 1, nDimensions

        localMetricsAlongDirection1 = grid%metrics(j,1+nDimensions*(i-1):nDimensions*i)

        select case (nDimensions)
        case (1)
           call computeJacobianOfInviscidFlux1D(localConservedVariables,                     &
                localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                localFluxJacobian1, specificVolume = state%specificVolume(j,1),              &
                velocity = localVelocity, temperature = state%temperature(j,1))
        case (2)
           call computeJacobianOfInviscidFlux2D(localConservedVariables,                     &
                localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                localFluxJacobian1, specificVolume = state%specificVolume(j,1),              &
                velocity = localVelocity, temperature = state%temperature(j,1))
        case (3)
           call computeJacobianOfInviscidFlux3D(localConservedVariables,                     &
                localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                localFluxJacobian1, specificVolume = state%specificVolume(j,1),              &
                velocity = localVelocity, temperature = state%temperature(j,1))
        end select !... nDimensions

        if (simulationFlags%viscosityOn) then
           select case (nDimensions)
           case (1)
              call computeFirstPartialViscousJacobian1D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian2, specificVolume = state%specificVolume(j,1),           &
                   velocity = localVelocity, temperature = state%temperature(j,1))
           case (2)
              call computeFirstPartialViscousJacobian2D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian2, specificVolume = state%specificVolume(j,1),           &
                   velocity = localVelocity, temperature = state%temperature(j,1))
           case (3)
              call computeFirstPartialViscousJacobian3D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian2, specificVolume = state%specificVolume(j,1),           &
                   velocity = localVelocity, temperature = state%temperature(j,1))
           end select
           localFluxJacobian1 = localFluxJacobian1 - localFluxJacobian2
        end if

        state%rightHandSide(j,:) = state%rightHandSide(j,:) +                                &
             matmul(transpose(localFluxJacobian1), temp1(j,:,i))

     end do !... i = 1, nDimensions

  end do !... j = 1, grid%nGridPoints

  SAFE_DEALLOCATE(localConservedVariables)
  SAFE_DEALLOCATE(localFluxJacobian1)
  SAFE_DEALLOCATE(localHeatFlux)
  SAFE_DEALLOCATE(localStressTensor)
  SAFE_DEALLOCATE(localFluxJacobian2)

  if (simulationFlags%viscosityOn) then

     allocate(temp2(grid%nGridPoints, nUnknowns - 1))

     allocate(localMetricsAlongDirection2(nDimensions))
     allocate(localFluxJacobian2(nUnknowns - 1, nUnknowns - 1))
     allocate(localAdjointDiffusion(nUnknowns - 1, nDimensions))

     do k = 1, grid%nGridPoints

        localVelocity = state%velocity(k,:)
        localAdjointDiffusion = 0.0_wp

        do j = 1, nDimensions

           localMetricsAlongDirection2 = grid%metrics(k,1+nDimensions*(j-1):nDimensions*j)

           do i = 1, nDimensions

              localMetricsAlongDirection1 = grid%metrics(k,1+nDimensions*(i-1):nDimensions*i)

              select case (nDimensions)
              case (1)
                 call computeSecondPartialViscousJacobian1D(localVelocity,                   &
                      state%dynamicViscosity(k,1), state%secondCoefficientOfViscosity(k,1),  &
                      state%thermalDiffusivity(k,1), grid%jacobian(k,1),                     &
                      localMetricsAlongDirection1(1), localFluxJacobian2)
              case (2)
                 call computeSecondPartialViscousJacobian2D(localVelocity,                   &
                      state%dynamicViscosity(k,1), state%secondCoefficientOfViscosity(k,1),  &
                      state%thermalDiffusivity(k,1), grid%jacobian(k,1),                     &
                      localMetricsAlongDirection1, localMetricsAlongDirection2,              &
                      localFluxJacobian2)
              case (3)
                 call computeSecondPartialViscousJacobian3D(localVelocity,                   &
                      state%dynamicViscosity(k,1), state%secondCoefficientOfViscosity(k,1),  &
                      state%thermalDiffusivity(k,1), grid%jacobian(k,1),                     &
                      localMetricsAlongDirection1, localMetricsAlongDirection2,              &
                      localFluxJacobian2)
              end select !... nDimensions

              localAdjointDiffusion(:,j) = localAdjointDiffusion(:,j) +                      &
                   matmul(transpose(localFluxJacobian2), temp1(k,2:nUnknowns,i))

           end do !... i = 1, nDimensions

        end do !... j = 1, nDimensions

        do j = 1, nDimensions
           temp1(k,2:nUnknowns,j) = localAdjointDiffusion(:,j)
        end do

     end do !... k = 1, grid%nGridPoints

     do j = 1, nDimensions
        call grid%adjointFirstDerivative(j)%apply(temp1(:,2:nUnknowns,j), grid%localSize)
     end do
     temp2 = sum(temp1(:,2:nUnknowns,:), dim = 3)

     temp2(:,nDimensions+1) = solverOptions%ratioOfSpecificHeats *                           &
          state%specificVolume(:,1) * temp2(:,nDimensions+1)
     do i = 1, nDimensions
        temp2(:,i) = state%specificVolume(:,1) * temp2(:,i) -                                &
             state%velocity(:,i) * temp2(:,nDimensions+1)
     end do

     state%rightHandSide(:,2:nUnknowns) = state%rightHandSide(:,2:nUnknowns) - temp2
     state%rightHandSide(:,1) = state%rightHandSide(:,1) +                                   &
          state%specificVolume(:,1) * state%conservedVariables(:,nDimensions+2) *            &
          temp2(:,nDimensions+1) + sum(state%velocity * temp2(:,1:nDimensions), dim = 2)

     SAFE_DEALLOCATE(temp2)

     SAFE_DEALLOCATE(localAdjointDiffusion)
     SAFE_DEALLOCATE(localFluxJacobian2)
     SAFE_DEALLOCATE(localMetricsAlongDirection2)

  end if !... simulationFlags%viscosityOn

  SAFE_DEALLOCATE(localMetricsAlongDirection1)
  SAFE_DEALLOCATE(localVelocity)
  SAFE_DEALLOCATE(temp1)

  ! Add dissipation if required.
  if (simulationFlags%dissipationOn)                                                         &
       call addDissipation(ADJOINT, simulationFlags, solverOptions, grid, state)

  ! Viscous penalties on far-field patches.
  if (simulationFlags%viscosityOn)                                                           &
       call addFarFieldAdjointPenalty(simulationFlags, solverOptions,                        &
       grid, state, patchFactories)

  call endTiming("computeRhsAdjoint")

end subroutine computeRhsAdjoint

subroutine computeRhsLinearized(simulationFlags, solverOptions, grid, state, patchFactories)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use Patch_mod, only : t_Patch
  use State_mod, only : t_State
  use Patch_factory, only : t_PatchFactory
  use FarFieldPatch_mod, only : t_FarFieldPatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  ! <<< Enumerations >>>
  use Region_enum, only : LINEARIZED

  ! <<< Private members >>>
  use RhsHelperImpl, only : addDissipation

  ! <<< Internal modules >>>
  use CNSHelper
  use Patch_factory, only : queryPatchTypeExists
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid) :: grid
  class(t_State) :: state
  type(t_PatchFactory), allocatable :: patchFactories(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, nUnknowns, normalDirection
  SCALAR_TYPE, allocatable :: fluxes1(:,:,:), fluxes2(:,:,:),                   &
              localFluxJacobian1(:,:), localConservedVariables(:),              &
              localVelocity(:), localMetricsAlongDirection1(:),                 &
              localFluxJacobian2(:,:), localStressTensor(:),                    &
              localHeatFlux(:), localLinearizedDiffusion(:,:),                  &
              localMetricsAlongDirection2(:), temp(:,:,:)
  class(t_Patch), pointer :: patch => null()

  call startTiming("computeRhsLinearized")

  nDimensions = grid%nDimensions
  nUnknowns = solverOptions%nUnknowns
  assert_key(nDimensions, (1, 2, 3))

  allocate(fluxes1(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))
  allocate(fluxes2(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))

  state%rightHandSide = 0.0_wp

  ! Compute Cartesian form of delta inviscid fluxes.
  fluxes1 = 0.0_wp

  allocate(localFluxJacobian1(nUnknowns, nUnknowns))
  allocate(localConservedVariables(nUnknowns))
  allocate(localVelocity(nDimensions))
  allocate(localMetricsAlongDirection1(nDimensions))
  do j = 1, grid%nGridPoints
     localConservedVariables = state%conservedVariables(j,:)
     localVelocity = state%velocity(j,:)
     do i = 1, nDimensions
        localMetricsAlongDirection1 = grid%metrics(j,1+nDimensions*(i-1):nDimensions*i)
        select case (nDimensions)
        case (1)
           call computeJacobianOfInviscidFlux1D(localConservedVariables,                     &
                localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                localFluxJacobian1, specificVolume = state%specificVolume(j,1),              &
                velocity = localVelocity, temperature = state%temperature(j,1))
        case (2)
           call computeJacobianOfInviscidFlux2D(localConservedVariables,                     &
                localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                localFluxJacobian1, specificVolume = state%specificVolume(j,1),              &
                velocity = localVelocity, temperature = state%temperature(j,1))
        case (3)
           call computeJacobianOfInviscidFlux3D(localConservedVariables,                     &
                localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                localFluxJacobian1, specificVolume = state%specificVolume(j,1),              &
                velocity = localVelocity, temperature = state%temperature(j,1))
        end select !... nDimensions

        fluxes1(j,:,i) = fluxes1(j,:,i) +                                                     &
                        matmul(localFluxJacobian1, state%adjointVariables(j,:))

     end do !... i = 1, nDimensions
  end do !... j = 1, grid%nGridPoints

  ! Compute Cartesian form of delta viscous fluxes if viscous terms are included and computed using
  ! repeated first derivatives.
  fluxes2 = 0.0_wp
  if (simulationFlags%viscosityOn .and. simulationFlags%repeatFirstDerivative) then
    ! for second partial viscous terms
    allocate(temp(grid%nGridPoints, nUnknowns-1, nDimensions))
    temp = 0.0_wp
    do i = 1, nDimensions
      temp(:,i,1) = - state%velocity(:,i) * state%adjointVariables(:,1)         &
                     + state%adjointVariables(:,i+1)
    end do
    temp(:,nUnknowns-1,1) = - state%specificVolume(:,1)                         &
                             * state%conservedVariables(:,nUnknowns)            &
                             * state%adjointVariables(:,1)
    temp(:,nUnknowns-1,1) = temp(:,nUnknowns-1,1)                               &
                    - sum(state%velocity * temp(:,1:nDimensions,1), dim=2)      &
                    + state%adjointVariables(:,nUnknowns)
    temp(:,nUnknowns-1,1) = temp(:,nUnknowns-1,1)                               &
                               * solverOptions%ratioOfSpecificHeats
    do i = 1, nUnknowns-1
      temp(:,i,1) = temp(:,i,1) * state%specificVolume(:,1)
    end do
    do i = 2, nDimensions
      temp(:,:,i) = temp(:,:,1)
    end do
    do i = 1, nDimensions
      call grid%firstDerivative(i)%apply(temp(:,:,i), grid%localSize)
    end do

    allocate(localFluxJacobian2(nUnknowns-1, nUnknowns-1))
    allocate(localStressTensor(nDimensions ** 2))
    allocate(localHeatFlux(nDimensions))
    allocate(localMetricsAlongDirection2(nDimensions))
    allocate(localLinearizedDiffusion(nUnknowns-1, nDimensions))

    do k = 1, grid%nGridPoints
       localConservedVariables = state%conservedVariables(k,:)
       localVelocity = state%velocity(k,:)
       localStressTensor = state%stressTensor(k,:)
       localHeatFlux = state%heatFlux(k,:)
       localLinearizedDiffusion = 0.0_wp
       do i = 1, nDimensions
          localMetricsAlongDirection1 = grid%metrics(k,1+nDimensions*(i-1):nDimensions*i)
           select case (nDimensions)
           case (1)
              call computeFirstPartialViscousJacobian1D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian1, specificVolume = state%specificVolume(k,1),           &
                   velocity = localVelocity, temperature = state%temperature(k,1))
           case (2)
              call computeFirstPartialViscousJacobian2D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian1, specificVolume = state%specificVolume(k,1),           &
                   velocity = localVelocity, temperature = state%temperature(k,1))
           case (3)
              call computeFirstPartialViscousJacobian3D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian1, specificVolume = state%specificVolume(k,1),           &
                   velocity = localVelocity, temperature = state%temperature(k,1))
           end select

          fluxes2(k,:,i) = fluxes2(k,:,i) +                                                    &
                          matmul(localFluxJacobian1, state%adjointVariables(k,:))

          do j = 1, nDimensions

             localMetricsAlongDirection2 = grid%metrics(k,1+nDimensions*(j-1):nDimensions*j)

             select case (nDimensions)
             case (1)
                call computeSecondPartialViscousJacobian1D(localVelocity,                   &
                     state%dynamicViscosity(k,1), state%secondCoefficientOfViscosity(k,1),  &
                     state%thermalDiffusivity(k,1), grid%jacobian(k,1),                     &
                     localMetricsAlongDirection1(1), localFluxJacobian2)
             case (2)
                call computeSecondPartialViscousJacobian2D(localVelocity,                   &
                     state%dynamicViscosity(k,1), state%secondCoefficientOfViscosity(k,1),  &
                     state%thermalDiffusivity(k,1), grid%jacobian(k,1),                     &
                     localMetricsAlongDirection1, localMetricsAlongDirection2,              &
                     localFluxJacobian2)
             case (3)
                call computeSecondPartialViscousJacobian3D(localVelocity,                   &
                     state%dynamicViscosity(k,1), state%secondCoefficientOfViscosity(k,1),  &
                     state%thermalDiffusivity(k,1), grid%jacobian(k,1),                     &
                     localMetricsAlongDirection1, localMetricsAlongDirection2,              &
                     localFluxJacobian2)
             end select !... nDimensions

             localLinearizedDiffusion(:,i) = localLinearizedDiffusion(:,i) +                      &
                                        matmul(localFluxJacobian2,temp(k,:,j))

          end do !... j = 1, nDimensions

          fluxes2(k,2:nUnknowns,i) = fluxes2(k,2:nUnknowns,i) + localLinearizedDiffusion(:,i)
       end do !... i = 1, nDimensions
    end do !... k = 1, grid%nGridPoints

  end if

  ! Send viscous fluxes to patches that will use it.
  if (simulationFlags%viscosityOn .and. allocated(patchFactories)) then
     do i = 1, size(patchFactories)
        call patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        if (patch%gridIndex /= grid%index .or. patch%nPatchPoints <= 0) cycle

        select type (patch)
        class is (t_FarFieldPatch)
           call patch%collect(fluxes2, patch%viscousFluxes)
        class is (t_BlockInterfacePatch)
          normalDirection = abs(patch%normalDirection)
           call patch%collect(fluxes2(:,:,normalDirection), patch%viscousFluxesL)
        end select

     end do
  end if

  fluxes1 = fluxes1 - fluxes2

  SAFE_DEALLOCATE(fluxes2) !... no longer needed.

  ! Take derivatives of fluxes.
  do i = 1, nDimensions
     call grid%firstDerivative(i)%apply(fluxes1(:,:,i), grid%localSize)
  end do
  state%rightHandSide = state%rightHandSide - sum(fluxes1, dim = 3) !... update RHS.

  SAFE_DEALLOCATE(fluxes1) !... no longer needed

  ! Add dissipation if required.
  if (simulationFlags%dissipationOn)                                                         &
       call addDissipation(LINEARIZED, simulationFlags, solverOptions, grid, state)

  call endTiming("computeRhsLinearized")

end subroutine computeRhsLinearized

subroutine addInterfaceAdjointPenalty(simulationFlags, solverOptions,                        &
     grid, state, patchFactories)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use Patch_mod, only : t_Patch
  use Patch_factory, only : t_PatchFactory
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  ! <<< Internal modules >>>
  use CNSHelper
  use Patch_factory, only : queryPatchTypeExists

  implicit none

  ! <<< Arguments >>>
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state
  type(t_PatchFactory), allocatable :: patchFactories(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  logical :: blockInterfacesExist
  integer :: i, j, k, l, iPatchFactory, nDimensions, nUnknowns,                              &
       direction, gridIndex, patchIndex, ierror
  SCALAR_TYPE, allocatable :: localConservedVariablesL(:), localVelocity(:),                 &
       localMetricsAlongDirection1(:), localMetricsAlongDirection2(:),                       &
       localViscousFluxJacobian(:,:), temp1(:,:,:), temp2(:,:)
  class(t_Patch), pointer :: patch => null()

  blockInterfacesExist = queryPatchTypeExists(patchFactories,                                &
       'SAT_BLOCK_INTERFACE', grid%index)
  call MPI_Allreduce(MPI_IN_PLACE, blockInterfacesExist, 1, MPI_LOGICAL,                     &
       MPI_LOR, grid%comm, ierror) !... reduce across grid-level processes.
  if (.not. blockInterfacesExist) return

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns >= nDimensions + 2)

  allocate(localConservedVariablesL(nUnknowns))
  allocate(localVelocity(nDimensions))
  allocate(localMetricsAlongDirection1(nDimensions))
  allocate(localMetricsAlongDirection2(nDimensions))
  allocate(localViscousFluxJacobian(nUnknowns - 1, nUnknowns - 1))
  allocate(temp1(grid%nGridPoints, nUnknowns - 1, nDimensions))
  allocate(temp2(grid%nGridPoints, nUnknowns - 1))

  temp1 = 0.0_wp

  if (allocated(patchFactories)) then
     do iPatchFactory = 1, size(patchFactories)
        call patchFactories(iPatchFactory)%connect(patch)
        if (.not. associated(patch)) cycle
        if (patch%gridIndex /= grid%index .or. patch%nPatchPoints <= 0) cycle

        assert(all(grid%offset == patch%gridOffset))
        assert(all(grid%localSize == patch%gridLocalSize))

        select type (patch)
        class is (t_BlockInterfacePatch)

          direction = abs(patch%normalDirection)
          assert(direction >= 1 .and. direction <= nDimensions)

           do k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)
              do j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
                 do i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
                    gridIndex = i - patch%gridOffset(1) + patch%gridLocalSize(1) *           &
                         (j - 1 - patch%gridOffset(2) + patch%gridLocalSize(2) *             &
                         (k - 1 - patch%gridOffset(3)))
                    if (grid%iblank(gridIndex) == 0) cycle
                    patchIndex = i - patch%offset(1) + patch%localSize(1) *                  &
                         (j - 1 - patch%offset(2) + patch%localSize(2) *                     &
                         (k - 1 - patch%offset(3)))

                    localConservedVariablesL = patch%conservedVariablesL(patchIndex,:)
                    localVelocity = state%velocity(gridIndex,:)
                    ! localMetricsAlongDirection1 = grid%metrics(gridIndex,                    &
                    !      1 + nDimensions * (direction - 1) : nDimensions * direction)

                    ! This is only discrete adjoint!!
                    do l = 1, nDimensions

                       localMetricsAlongDirection2 =                                         &
                            grid%metrics(gridIndex,1+nDimensions*(l-1):nDimensions*l)

                       localMetricsAlongDirection1 = patch%metricsAlongNormalDirectionL(patchIndex,:)

                       select case (nDimensions)
                       case (1)
                          call computeSecondPartialViscousJacobian1D(localVelocity,          &
                               state%dynamicViscosity(gridIndex,1),                          &
                               state%secondCoefficientOfViscosity(gridIndex,1),              &
                               state%thermalDiffusivity(gridIndex,1),                        &
                               grid%jacobian(gridIndex,1), localMetricsAlongDirection1(1),   &
                               localViscousFluxJacobian)
                       case (2)
                          call computeSecondPartialViscousJacobian2D(localVelocity,          &
                               state%dynamicViscosity(gridIndex,1),                          &
                               state%secondCoefficientOfViscosity(gridIndex,1),              &
                               state%thermalDiffusivity(gridIndex,1),                        &
                               grid%jacobian(gridIndex,1), localMetricsAlongDirection1,      &
                               localMetricsAlongDirection2, localViscousFluxJacobian)
                       case (3)
                          call computeSecondPartialViscousJacobian3D(localVelocity,          &
                               state%dynamicViscosity(gridIndex,1),                          &
                               state%secondCoefficientOfViscosity(gridIndex,1),              &
                               state%thermalDiffusivity(gridIndex,1),                        &
                               grid%jacobian(gridIndex,1), localMetricsAlongDirection1,      &
                               localMetricsAlongDirection2, localViscousFluxJacobian)
                       end select !... nDimensions

                       ! Note sign change is included in viscousPenaltyAmount
                       temp1(gridIndex,:,l) = temp1(gridIndex,:,l) -                                        &
                            patch%viscousPenaltyAmountL * matmul(transpose(localViscousFluxJacobian),       &
                            patch%adjointVariablesL(patchIndex,2:nUnknowns) )

                       localMetricsAlongDirection1 = patch%metricsAlongNormalDirectionR(patchIndex,:)

                       select case (nDimensions)
                       case (1)
                           call computeSecondPartialViscousJacobian1D(localVelocity,          &
                                state%dynamicViscosity(gridIndex,1),                          &
                                state%secondCoefficientOfViscosity(gridIndex,1),              &
                                state%thermalDiffusivity(gridIndex,1),                        &
                                grid%jacobian(gridIndex,1), localMetricsAlongDirection1(1),   &
                                localViscousFluxJacobian)
                       case (2)
                           call computeSecondPartialViscousJacobian2D(localVelocity,          &
                                state%dynamicViscosity(gridIndex,1),                          &
                                state%secondCoefficientOfViscosity(gridIndex,1),              &
                                state%thermalDiffusivity(gridIndex,1),                        &
                                grid%jacobian(gridIndex,1), localMetricsAlongDirection1,      &
                                localMetricsAlongDirection2, localViscousFluxJacobian)
                       case (3)
                           call computeSecondPartialViscousJacobian3D(localVelocity,          &
                                state%dynamicViscosity(gridIndex,1),                          &
                                state%secondCoefficientOfViscosity(gridIndex,1),              &
                                state%thermalDiffusivity(gridIndex,1),                        &
                                grid%jacobian(gridIndex,1), localMetricsAlongDirection1,      &
                                localMetricsAlongDirection2, localViscousFluxJacobian)
                       end select !... nDimensions

                       ! Note sign change is included in viscousPenaltyAmount
                       temp1(gridIndex,:,l) = temp1(gridIndex,:,l) +                                        &
                            patch%viscousPenaltyAmountR * matmul(transpose(localViscousFluxJacobian),       &
                            patch%adjointVariablesR(patchIndex,2:nUnknowns) )

                    end do !... l = 1, nDimensions

                 end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
              end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
           end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

        end select !... type (patch)

     end do !... iPatchFactory = 1, size(patchFactories)
  end if !... allocated(patchFactories)

  do i = 1, nDimensions
     call grid%adjointFirstDerivative(i)%apply(temp1(:,:,i), grid%localSize)
  end do
  temp2 = sum(temp1, dim = 3)

  temp2(:,nDimensions+1) = solverOptions%ratioOfSpecificHeats *                              &
       state%specificVolume(:,1) * temp2(:,nDimensions+1)
  do i = 1, nDimensions
     temp2(:,i) = state%specificVolume(:,1) * temp2(:,i) -                                   &
          state%velocity(:,i) * temp2(:,nDimensions+1)
  end do

  state%rightHandSide(:,2:nUnknowns) = state%rightHandSide(:,2:nUnknowns) + temp2
  state%rightHandSide(:,1) = state%rightHandSide(:,1) -                                      &
       state%specificVolume(:,1) * state%conservedVariables(:,nDimensions+2) *               &
       temp2(:,nDimensions+1) - sum(state%velocity * temp2(:,1:nDimensions), dim = 2)

  SAFE_DEALLOCATE(temp2)
  SAFE_DEALLOCATE(temp1)
  SAFE_DEALLOCATE(localViscousFluxJacobian)
  SAFE_DEALLOCATE(localMetricsAlongDirection2)
  SAFE_DEALLOCATE(localMetricsAlongDirection1)
  SAFE_DEALLOCATE(localVelocity)
  SAFE_DEALLOCATE(localConservedVariablesL)

end subroutine addInterfaceAdjointPenalty
