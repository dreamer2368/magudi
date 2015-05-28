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
    use Region_enum, only : FORWARD, ADJOINT

    ! <<< Internal modules >>>
    use CNSHelper, only : computeSpectralRadius
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    integer, intent(in) :: mode
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid), intent(in) :: grid
    class(t_State) :: state

    ! <<< Local variables >>>
    integer :: i, j, nDimensions, nUnknowns
    SCALAR_TYPE, allocatable :: dissipationTerm(:,:), spectralRadius(:,:)
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
    end select

    if (simulationFlags%compositeDissipation) then

       do i = 1, nDimensions

          select case (mode)
          case (FORWARD)
             dissipationTerm = state%conservedVariables
          case (ADJOINT)
             dissipationTerm = state%adjointVariables
          end select

          call grid%dissipation(i)%apply(dissipationTerm, grid%localSize)
          state%rightHandSide = state%rightHandSide + dissipationAmount * dissipationTerm

       end do

    else

       allocate(spectralRadius(grid%nGridPoints, nDimensions))
       call computeSpectralRadius(nDimensions, solverOptions%ratioOfSpecificHeats,           &
            state%velocity, state%temperature(:,1), grid%metrics,                            &
            spectralRadius, grid%isCurvilinear)

       do i = 1, nDimensions

          select case (mode)
          case (FORWARD)
             dissipationTerm = state%conservedVariables
          case (ADJOINT)
             dissipationTerm = state%adjointVariables
          end select

          call grid%dissipation(i)%apply(dissipationTerm, grid%localSize)

          do j = 1, nUnknowns
             dissipationTerm(:,j) = spectralRadius(:,i) * dissipationTerm(:,j)
          end do

          call grid%dissipationTranspose(i)%apply(dissipationTerm, grid%localSize)
          call grid%firstDerivative(i)%applyNormInverse(dissipationTerm, grid%localSize)

          state%rightHandSide = state%rightHandSide - dissipationAmount * dissipationTerm

       end do

       SAFE_DEALLOCATE(spectralRadius)

    end if

    SAFE_DEALLOCATE(dissipationTerm)

    call endTiming("addDissipation")

  end subroutine addDissipation

  subroutine computeForwardFarFieldViscousPenalty(simulationFlags,                           &
       solverOptions, grid, state, patchFactories)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use Patch_mod, only : t_Patch
    use State_mod, only : t_State
    use Patch_factory, only : t_PatchFactory
    use FarFieldPatch_mod, only : t_FarFieldPatch
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Arguments >>>
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid) :: grid
    class(t_State) :: state
    type(t_PatchFactory), allocatable :: patchFactories(:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, nDimensions, nUnknowns, nSpecies
    class(t_Patch), pointer :: patch => null()
    SCALAR_TYPE, allocatable :: patchConservedVariables(:,:), patchTargetState(:,:),         &
         viscousPenaltyAlongDirection(:,:), temp(:,:)

    if (.not. simulationFlags%viscosityOn) return

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    nSpecies = solverOptions%nSpecies
    assert(nSpecies >= 0)

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns == nDimensions + nSpecies + 2)

    assert(grid%nGridPoints > 0)
    assert(allocated(state%targetState))
    assert(size(state%targetState, 1) == grid%nGridPoints)
    assert(size(state%targetState, 2) == nUnknowns)

    allocate(temp(grid%nGridPoints, nUnknowns - 1))

    do i = 1, nDimensions

       temp = 0.0_wp

       if (allocated(patchFactories)) then
          do j = 1, size(patchFactories)
             call patchFactories(j)%connect(patch)
             if (.not. associated(patch)) cycle
             if (patch%gridIndex /= grid%index .or. patch%nPatchPoints <= 0) cycle
             select type (patch)
             class is (t_FarFieldPatch)

                allocate(patchConservedVariables(patch%nPatchPoints, nUnknowns))
                allocate(patchTargetState(patch%nPatchPoints, nUnknowns))
                call patch%collect(state%conservedVariables, patchConservedVariables)
                call patch%collect(state%targetState, patchTargetState)

                patchConservedVariables = patchConservedVariables - patchTargetState
                do k = 1, patch%nPatchPoints
                   patchConservedVariables(k,2:nDimensions+1) =                              &
                        (patchConservedVariables(k,2:nDimensions+1) -                        &
                        patchTargetState(k,2:nDimensions+1) * patchConservedVariables(k,1) / &
                        patchTargetState(k,1)) / patchTargetState(k,1)
                   patchConservedVariables(k,nDimensions+2) =                                &
                        solverOptions%ratioOfSpecificHeats *                                 &
                        (patchConservedVariables(k,nDimensions+2) -                          &
                        sum(patchConservedVariables(k,2:nDimensions+1) *                     &
                        patchTargetState(k,2:nDimensions+1)) -                               &
                        patchTargetState(k,nDimensions+2) * patchConservedVariables(k,1) /   &
                        patchTargetState(k,1)) / patchTargetState(k,1)
                end do

                call patch%disperseAdd(patchConservedVariables(:,2:nUnknowns), temp)

                SAFE_DEALLOCATE(patchTargetState)
                SAFE_DEALLOCATE(patchConservedVariables)

             end select !... type(patch)
          end do !... j = 1, size(patchFactories)
       end if !... allocated(patchFactories)

       call grid%firstDerivative(i)%apply(temp, grid%localSize)

       if (allocated(patchFactories)) then
          do j = 1, size(patchFactories)
             call patchFactories(j)%connect(patch)
             if (.not. associated(patch)) cycle
             if (patch%gridIndex /= grid%index .or. patch%nPatchPoints <= 0) cycle
             select type (patch)
             class is (t_FarFieldPatch)

                allocate(viscousPenaltyAlongDirection(patch%nPatchPoints, nUnknowns - 1))

                patch%viscousPenalty(:,1) = 0.0_wp
                call patch%collect(temp, viscousPenaltyAlongDirection)
                do k = 1, patch%nPatchPoints
                   viscousPenaltyAlongDirection(k,:) =                                       &
                        matmul(patch%secondPartialViscousJacobians(k,:,:,i),                 &
                        viscousPenaltyAlongDirection(k,:))
                end do

                if (i == 1) then
                   patch%viscousPenalty(:,1) = 0.0_wp
                   patch%viscousPenalty(:,2:nUnknowns) = viscousPenaltyAlongDirection
                else
                   patch%viscousPenalty(:,2:nUnknowns) =                                     &
                        patch%viscousPenalty(:,2:nUnknowns) + viscousPenaltyAlongDirection
                end if

                SAFE_DEALLOCATE(viscousPenaltyAlongDirection)

                if (i == nDimensions) then

                   allocate(patchConservedVariables(patch%nPatchPoints, nUnknowns))
                   allocate(patchTargetState(patch%nPatchPoints, nUnknowns))
                   call patch%collect(state%conservedVariables, patchConservedVariables)
                   call patch%collect(state%targetState, patchTargetState)

                   do k = 1, patch%nPatchPoints
                      patch%viscousPenalty(k,:) = patch%viscousPenalty(k,:) +                &
                           matmul(patch%firstPartialViscousJacobians(k,:,:),                 &
                           patchConservedVariables(k,:) - patchTargetState(k,:))
                   end do

                   SAFE_DEALLOCATE(patchTargetState)
                   SAFE_DEALLOCATE(patchConservedVariables)

                end if

             end select !... type(patch)
          end do !... j = 1, size(patchFactories)
       end if !... allocated(patchFactories)

    end do !... i = 1, nDimensions

    SAFE_DEALLOCATE(temp)

  end subroutine computeForwardFarFieldViscousPenalty

  subroutine computeAdjointFarFieldViscousPenalty(simulationFlags,                           &
       solverOptions, grid, state, patchFactories)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use Patch_mod, only : t_Patch
    use State_mod, only : t_State
    use Patch_factory, only : t_PatchFactory
    use FarFieldPatch_mod, only : t_FarFieldPatch
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Arguments >>>
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid) :: grid
    class(t_State) :: state
    type(t_PatchFactory), allocatable :: patchFactories(:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, nDimensions, nUnknowns, nSpecies
    class(t_Patch), pointer :: patch => null()
    SCALAR_TYPE, allocatable :: patchAdjointVariables(:,:), patchTargetState(:,:),           &
         viscousPenaltyAlongDirection(:,:), temp(:,:)

    if (.not. simulationFlags%viscosityOn) return

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    nSpecies = solverOptions%nSpecies
    assert(nSpecies >= 0)

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns == nDimensions + nSpecies + 2)

    assert(grid%nGridPoints > 0)
    assert(allocated(state%targetState))
    assert(size(state%targetState, 1) == grid%nGridPoints)
    assert(size(state%targetState, 2) == nUnknowns)

    allocate(temp(grid%nGridPoints, nUnknowns - 1))

    do i = 1, nDimensions

       temp = 0.0_wp

       if (allocated(patchFactories)) then
          do j = 1, size(patchFactories)
             call patchFactories(j)%connect(patch)
             if (.not. associated(patch)) cycle
             if (patch%gridIndex /= grid%index .or. patch%nPatchPoints <= 0) cycle
             select type (patch)
             class is (t_FarFieldPatch)

                allocate(patchAdjointVariables(patch%nPatchPoints, nUnknowns))
                call patch%collect(state%adjointVariables, patchAdjointVariables)

                do k = 1, patch%nPatchPoints
                   patchAdjointVariables(k,1) = 0.0_wp
                   patchAdjointVariables(k,2:nUnknowns) =                                    &
                        matmul(transpose(patch%secondPartialViscousJacobians(k,:,:,i)),      &
                        patchAdjointVariables(k,2:nUnknowns))
                end do

                call patch%disperseAdd(patchAdjointVariables(:,2:nUnknowns), temp)

                SAFE_DEALLOCATE(patchAdjointVariables)

             end select !... type(patch)
          end do !... j = 1, size(patchFactories)
       end if !... allocated(patchFactories)

       call grid%adjointFirstDerivative(i)%apply(temp, grid%localSize)

       if (allocated(patchFactories)) then
          do j = 1, size(patchFactories)
             call patchFactories(j)%connect(patch)
             if (.not. associated(patch)) cycle
             if (patch%gridIndex /= grid%index .or. patch%nPatchPoints <= 0) cycle
             select type (patch)
             class is (t_FarFieldPatch)

                allocate(viscousPenaltyAlongDirection(patch%nPatchPoints, nUnknowns))
                allocate(patchTargetState(patch%nPatchPoints, nUnknowns))
                call patch%collect(state%targetState, patchTargetState)

                viscousPenaltyAlongDirection(:,1) = 0.0_wp
                call patch%collect(temp, viscousPenaltyAlongDirection(:,2:nUnknowns))

                viscousPenaltyAlongDirection(:,nDimensions+2) =                              &
                     solverOptions%ratioOfSpecificHeats *                                    &
                     viscousPenaltyAlongDirection(:,nDimensions+2) / patchTargetState(:,1)
                do k = 1, nDimensions
                   viscousPenaltyAlongDirection(:,k+1) =                                     &
                        (viscousPenaltyAlongDirection(:,k+1) - patchTargetState(:,k+1) *     &
                        viscousPenaltyAlongDirection(:,nDimensions+2)) / patchTargetState(:,1)
                end do
                viscousPenaltyAlongDirection(:,1) =                                          &
                     - sum(patchTargetState(:,2:nDimensions+2) *                             &
                     viscousPenaltyAlongDirection(:,2:nDimensions+2), dim = 2) /             &
                     patchTargetState(:,1)

                if (i == 1) then
                   patch%viscousPenalty = viscousPenaltyAlongDirection
                else
                   patch%viscousPenalty = patch%viscousPenalty + viscousPenaltyAlongDirection
                end if

                SAFE_DEALLOCATE(patchTargetState)
                SAFE_DEALLOCATE(viscousPenaltyAlongDirection)

                if (i == nDimensions) then

                   allocate(patchAdjointVariables(patch%nPatchPoints, nUnknowns))
                   call patch%collect(state%adjointVariables, patchAdjointVariables)

                   do k = 1, patch%nPatchPoints
                      patch%viscousPenalty(k,:) = patch%viscousPenalty(k,:) +                &
                           matmul(transpose(patch%firstPartialViscousJacobians(k,:,:)),      &
                           patchAdjointVariables(k,:))
                   end do

                   SAFE_DEALLOCATE(patchAdjointVariables)

                end if

             end select !... type(patch)
          end do !... j = 1, size(patchFactories)
       end if !... allocated(patchFactories)

    end do !... i = 1, nDimensions

    SAFE_DEALLOCATE(temp)

  end subroutine computeAdjointFarFieldViscousPenalty

end module RhsHelperImpl

subroutine computeRhsForward(simulationFlags, solverOptions, grid, state, patchFactories)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use Patch_factory, only : t_PatchFactory
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD

  ! <<< Private members >>>
  use RhsHelperImpl, only : addDissipation, computeForwardFarFieldViscousPenalty

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
  integer :: i, nDimensions, nSpecies, ierror
  SCALAR_TYPE, allocatable :: fluxes1(:,:,:), fluxes2(:,:,:)
  logical :: flag

  call startTiming("computeRhsForward")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nSpecies = solverOptions%nSpecies
  assert(nSpecies >= 0)

  allocate(fluxes1(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))
  allocate(fluxes2(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))

  state%rightHandSide = 0.0_wp

  ! Compute Cartesian form of inviscid fluxes.
  call computeCartesianInvsicidFluxes(nDimensions, nSpecies, state%conservedVariables,       &
       state%velocity, state%pressure(:,1), fluxes1)

  ! Compute Cartesian form of viscous fluxes if viscous terms are included and computed using
  ! repeated first derivatives.
  if (simulationFlags%viscosityOn .and. simulationFlags%repeatFirstDerivative) then
     call computeCartesianViscousFluxes(nDimensions, nSpecies, state%velocity,               &
          state%massFraction, state%stressTensor, state%heatFlux, state%speciesFlux,         &
          fluxes2)
     fluxes1 = fluxes1 - fluxes2 !... Cartesian form of total fluxes.
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

  ! Update penalties on patches that require grid-level derivatives.
  flag = queryPatchTypeExists(patchFactories, 'SAT_FAR_FIELD', grid%index)
  call MPI_Allreduce(MPI_IN_PLACE, flag, 1, MPI_LOGICAL, MPI_LOR, grid%comm, ierror)
  if (flag) call computeForwardFarFieldViscousPenalty(simulationFlags,                       &
       solverOptions, grid, state, patchFactories)

  call endTiming("computeRhsForward")

end subroutine computeRhsForward

subroutine computeRhsAdjoint(simulationFlags, solverOptions, combustion, grid, state,        &
     patchFactories)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use Patch_factory, only : t_PatchFactory
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use Combustion_mod, only : t_Combustion

  ! <<< Enumerations >>>
  use Region_enum, only : ADJOINT

  ! <<< Private members >>>
  use RhsHelperImpl, only : addDissipation, computeAdjointFarFieldViscousPenalty

  ! <<< Internal modules >>>
  use CNSHelper
  use Patch_factory, only : queryPatchTypeExists
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  type(t_Combustion), intent(in) :: combustion
  class(t_Grid) :: grid
  class(t_State) :: state
  type(t_PatchFactory), allocatable :: patchFactories(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, nSpecies, nUnknowns, ierror
  SCALAR_TYPE, allocatable :: temp1(:,:,:), temp2(:,:),                                      &
       localFluxJacobian1(:,:), localFluxJacobian2(:,:), localConservedVariables(:),         &
       localVelocity(:), localMassFraction(:), localMetricsAlongDirection1(:),               &
       localMetricsAlongDirection2(:), localStressTensor(:), localHeatFlux(:),               &
       localAdjointDiffusion(:,:), localSpeciesFlux(:,:), localSourceJacobian(:,:)
  logical :: flag

  call startTiming("computeRhsAdjoint")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))
  nspecies = solverOptions%nSpecies
  assert(nSpecies >= 0)
  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2 + nSpecies)

  allocate(temp1(grid%nGridPoints, nUnknowns, nDimensions))

  state%rightHandSide = 0.0_wp

  ! Partial derivatives of adjoint variables w.r.t. *computational* coordinates.
  do i = 1, nDimensions
     temp1(:,:,i) = state%adjointVariables
     call grid%adjointFirstDerivative(i)%apply(temp1(:,:,i), grid%localSize)
  end do

  allocate(localFluxJacobian1(nUnknowns, nUnknowns))
  allocate(localConservedVariables(solverOptions%nUnknowns))
  allocate(localVelocity(nDimensions))
  allocate(localMassFraction(nSpecies))
  allocate(localMetricsAlongDirection1(nDimensions))

  if (simulationFlags%viscosityOn) then
     allocate(localFluxJacobian2(solverOptions%nUnknowns, solverOptions%nUnknowns))
     allocate(localStressTensor(nDimensions ** 2))
     allocate(localHeatFlux(nDimensions))
     allocate(localSpeciesFlux(nSpecies,nDimensions))
  end if

  if (nSpecies > 0 .and. combustion%nReactions > 0) then
     allocate(localSourceJacobian(nUnknowns, nUnknowns))
  end if

  do j = 1, grid%nGridPoints

     localConservedVariables = state%conservedVariables(j,:)
     localVelocity = state%velocity(j,:)
     if (simulationFlags%viscosityOn) then
        localStressTensor = state%stressTensor(j,:)
        localHeatFlux = state%heatFlux(j,:)
        localSpeciesFlux = state%speciesFlux(j,:,:)
     end if

     if (nSpecies > 0) then
        localMassFraction = state%massFraction(j,:)
     end if

     do i = 1, nDimensions

        localMetricsAlongDirection1 = grid%metrics(j,1+nDimensions*(i-1):nDimensions*i)

        call computeJacobianOfInviscidFlux(nDimensions, nSpecies,                            &
             localConservedVariables, localMetricsAlongDirection1,                           &
             solverOptions%ratioOfSpecificHeats, localFluxJacobian1,                         &
             specificVolume = state%specificVolume(j,1), velocity = localVelocity,           &
             temperature = state%temperature(j,1), massFraction = localMassFraction)

        if (simulationFlags%viscosityOn) then
           call computeFirstPartialViscousJacobian(nDimensions, nSpecies,                    &
                localConservedVariables, localMetricsAlongDirection1, localStressTensor,     &
                localHeatFlux, localSpeciesFlux, solverOptions%powerLawExponent,             &
                solverOptions%ratioOfSpecificHeats, localFluxJacobian2,                      &
                specificVolume = state%specificVolume(j,1), velocity = localVelocity,        &
                temperature = state%temperature(j,1), massFraction = localMassFraction)
           localFluxJacobian1 = localFluxJacobian1 - localFluxJacobian2
        end if

        if (nSpecies > 0 .and. combustion%nReactions > 0) then
           call computeJacobianOfSource(nDimensions, nSpecies,                               &
                localConservedVariables, localMetricsAlongDirection1,                        &
                solverOptions%ratioOfSpecificHeats, combustion,                              &
                localSourceJacobian, specificVolume = state%specificVolume(j,1),             &
                velocity = localVelocity, temperature = state%temperature(j,1),              &
                massFraction = localMassFraction)
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

              call computeSecondPartialViscousJacobian(nDimensions, nSpecies,                &
                   localVelocity, state%dynamicViscosity(k,1),                               &
                   state%secondCoefficientOfViscosity(k,1),                                  &
                   state%thermalDiffusivity(k,1), state%massDiffusivity(k,:),                &
                   grid%jacobian(k,1), localMetricsAlongDirection1,                          &
                   localMetricsAlongDirection2, localFluxJacobian2)

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

  if (nSpecies > 0 .and. combustion%nReactions > 0) then
     SAFE_DEALLOCATE(localSourceJacobian)
  end if

  SAFE_DEALLOCATE(localMetricsAlongDirection1)
  SAFE_DEALLOCATE(localVelocity)
  SAFE_DEALLOCATE(temp1)

  ! Add dissipation if required.
  if (simulationFlags%dissipationOn)                                                         &
       call addDissipation(ADJOINT, simulationFlags, solverOptions, grid, state)

  ! Update penalties on patches that require grid-level derivatives.
  flag = queryPatchTypeExists(patchFactories, 'SAT_FAR_FIELD', grid%index)
  call MPI_Allreduce(MPI_IN_PLACE, flag, 1, MPI_LOGICAL, MPI_LOR, grid%comm, ierror)
  if (flag) call computeAdjointFarFieldViscousPenalty(simulationFlags,                       &
       solverOptions, grid, state, patchFactories)

  call endTiming("computeRhsAdjoint")

end subroutine computeRhsAdjoint
