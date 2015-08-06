#include "config.h"

module RhsHelper

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  implicit none
  private

  public :: computeRhsForward, computeRhsAdjoint

contains

  subroutine computeRhsForward(simulationFlags, solverOptions, grid, state, patchFactories)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use Patch_mod, only : t_Patch
    use State_mod, only : t_State
    use Patch_factory, only : t_PatchFactory
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags
    use InflowOutflowPatch_mod, only : t_InflowOutflowPatch
    use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : FORWARD

    ! <<< Internal modules >>>
    use CNSHelper
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
    real(SCALAR_KIND), allocatable :: fluxes1(:,:,:), fluxes2(:,:,:)
    class(t_Patch), pointer :: patch => null()

    call startTiming("computeRhsForward")

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    allocate(fluxes1(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))
    allocate(fluxes2(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))

    state%rightHandSide = 0.0_wp

    ! Compute Cartesian form of inviscid fluxes.
    call computeCartesianInvsicidFluxes(nDimensions, state%conservedVariables,               &
         state%velocity, state%pressure(:,1), fluxes1)

    ! Compute Cartesian form of viscous fluxes if viscous terms are included and computed
    ! using repeated first derivatives.
    if (simulationFlags%viscosityOn .and. simulationFlags%repeatFirstDerivative) then
       call computeCartesianViscousFluxes(nDimensions, state%velocity,                       &
            state%stressTensor, state%heatFlux, fluxes2)
       fluxes1 = fluxes1 - fluxes2 !... Cartesian form of total fluxes.
    end if

    ! Send viscous fluxes to patches that will use it.
    if (simulationFlags%viscosityOn .and. allocated(patchFactories)) then
       do i = 1, size(patchFactories)
          call patchFactories(i)%connect(patch)
          if (.not. associated(patch)) cycle
          if (patch%gridIndex /= grid%index) cycle
          select type (patch)
          class is (t_InflowOutflowPatch)
             call patch%collect(fluxes2, patch%viscousFluxes)
          class is (t_BlockInterfacePatch)
             call patch%collect(fluxes2, patch%cartesianViscousFluxes)
          end select
       end do
    end if

    ! Transform fluxes from Cartesian to contravariant form: `fluxes1` has the Cartesian form
    ! of total fluxes... upon return, `fluxes2` has the contravariant form.
    call transformFluxes(nDimensions, fluxes1, grid%metrics, fluxes2, grid%isCurvilinear)

    SAFE_DEALLOCATE(fluxes1) !... no longer needed.

    ! Take derivatives of fluxes.
    do i = 1, nDimensions
       call grid%firstDerivative(i)%apply(fluxes2(:,:,i), grid%localSize)
    end do
    state%rightHandSide = state%rightHandSide - sum(fluxes2, dim = 3) !... update RHS.

    SAFE_DEALLOCATE(fluxes2) !... no longer needed

    ! Add dissipation if required.
    if (simulationFlags%dissipationOn)                                                       &
         call addDissipation(FORWARD, simulationFlags, solverOptions, grid, state)

    call endTiming("computeRhsForward")

  end subroutine computeRhsForward

  subroutine computeRhsAdjoint(simulationFlags, solverOptions, grid, state)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod, only : t_State
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : ADJOINT

    ! <<< Internal modules >>>
    use CNSHelper
    use MPITimingsHelper, only : startTiming, endTiming

    implicit none

    ! <<< Arguments >>>
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid) :: grid
    class(t_State) :: state

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, nDimensions, nUnknowns
    real(SCALAR_KIND), allocatable :: adjointGradient(:,:,:), fluxJacobian1(:,:),            &
         fluxJacobian2(:,:), localAdjointDiffusion(:,:), adjointDiffusion(:,:)

    call startTiming("computeRhsAdjoint")

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns >= nDimensions + 2)

    allocate(adjointGradient(grid%nGridPoints, nUnknowns, nDimensions))
    allocate(fluxJacobian1(nUnknowns, nUnknowns))
    if (simulationFlags%viscosityOn) allocate(fluxJacobian2(nUnknowns, nUnknowns))

    state%rightHandSide = 0.0_wp

    ! Partial derivatives of adjoint variables w.r.t. *computational* coordinates.
    do i = 1, nDimensions
       adjointGradient(:,:,i) = state%adjointVariables
       call grid%adjointFirstDerivative(i)%apply(adjointGradient(:,:,i), grid%localSize)
    end do

    do j = 1, grid%nGridPoints

       do i = 1, nDimensions

          call computeInviscidJacobian(nDimensions, j, state%conservedVariables,             &
               grid%metrics(:,1+nDimensions*(i-1):nDimensions*i),                            &
               solverOptions%ratioOfSpecificHeats, fluxJacobian1,                            &
               specificVolume = state%specificVolume(:,1), velocity = state%velocity,        &
               temperature = state%temperature(:,1))

          if (simulationFlags%viscosityOn) then
             call computeFirstPartialViscousJacobian(nDimensions, j,                         &
                  state%conservedVariables,                                                  &
                  grid%metrics(:,1+nDimensions*(i-1):nDimensions*i), state%stressTensor,     &
                  state%heatFlux, solverOptions%powerLawExponent,                            &
                  solverOptions%ratioOfSpecificHeats, fluxJacobian2,                         &
                  specificVolume = state%specificVolume(:,1), velocity = state%velocity,     &
                  temperature = state%temperature(:,1))
             fluxJacobian1 = fluxJacobian1 - fluxJacobian2
          end if

          state%rightHandSide(j,:) = state%rightHandSide(j,:) +                              &
               matmul(transpose(fluxJacobian1), adjointGradient(j,:,i))

       end do !... i = 1, nDimensions

    end do !... j = 1, grid%nGridPoints

    SAFE_DEALLOCATE(fluxJacobian1)
    SAFE_DEALLOCATE(fluxJacobian2)

    if (simulationFlags%viscosityOn) then

       allocate(fluxJacobian2(nUnknowns - 1, nUnknowns - 1))
       allocate(localAdjointDiffusion(nUnknowns - 1, nDimensions))

       do k = 1, grid%nGridPoints

          do j = 1, nDimensions
             do i = 1, nDimensions
                call computeSecondPartialViscousJacobian(nDimensions, k, state%velocity,     &
                     state%dynamicViscosity(:,1), state%secondCoefficientOfViscosity(:,1),   &
                     state%thermalDiffusivity(:,1), grid%jacobian(:,1),                      &
                     grid%metrics(:,1+nDimensions*(i-1):nDimensions*i),                      &
                     grid%metrics(:,1+nDimensions*(j-1):nDimensions*j), fluxJacobian2)
                localAdjointDiffusion(:,j) = localAdjointDiffusion(:,j) +                    &
                     matmul(transpose(fluxJacobian2), adjointGradient(k,2:nUnknowns,i))
             end do !... i = 1, nDimensions
          end do !... j = 1, nDimensions

          do j = 1, nDimensions
             adjointGradient(k,2:nUnknowns,j) = localAdjointDiffusion(:,j) !... re-use memory.
          end do

       end do !... k = 1, grid%nGridPoints

       do j = 1, nDimensions
          call grid%adjointFirstDerivative(j)%apply(adjointGradient(:,2:nUnknowns,j),        &
               grid%localSize)
       end do

       allocate(adjointDiffusion(grid%nGridPoints, nUnknowns - 1))
       adjointDiffusion = sum(adjointGradient(:,2:nUnknowns,:), dim = 3)

       adjointDiffusion(:,nDimensions+1) = solverOptions%ratioOfSpecificHeats *              &
            state%specificVolume(:,1) * adjointDiffusion(:,nDimensions+1)
       do i = 1, nDimensions
          adjointDiffusion(:,i) = state%specificVolume(:,1) * adjointDiffusion(:,i) -        &
               state%velocity(:,i) * adjointDiffusion(:,nDimensions+1)
       end do

       state%rightHandSide(:,2:nUnknowns) = state%rightHandSide(:,2:nUnknowns) -             &
            adjointDiffusion
       state%rightHandSide(:,1) = state%rightHandSide(:,1) + state%specificVolume(:,1) *     &
            state%conservedVariables(:,nDimensions+2) * adjointDiffusion(:,nDimensions+1) +  &
            sum(state%velocity * adjointDiffusion(:,1:nDimensions), dim = 2)

    end if !... simulationFlags%viscosityOn

    SAFE_DEALLOCATE(adjointGradient)
    SAFE_DEALLOCATE(fluxJacobian1)
    SAFE_DEALLOCATE(fluxJacobian2)
    SAFE_DEALLOCATE(localAdjointDiffusion)
    SAFE_DEALLOCATE(adjointDiffusion)

    ! Add dissipation if required.
    if (simulationFlags%dissipationOn)                                                       &
         call addDissipation(ADJOINT, simulationFlags, solverOptions, grid, state)

    call endTiming("computeRhsAdjoint")

  end subroutine computeRhsAdjoint

  subroutine addDissipation(mode, simulationFlags, solverOptions, grid, state)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod, only : t_State
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : FORWARD, ADJOINT

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    integer, intent(in) :: mode
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid), intent(in) :: grid
    class(t_State) :: state

    ! <<< Local variables >>>
    integer :: i, j, nDimensions
    real(SCALAR_KIND), allocatable :: dissipationTerm(:,:)
    real(SCALAR_KIND) :: dissipationAmount

    if (.not. simulationFlags%dissipationOn) return

    call startTiming("addDissipation")

    assert_key(mode, (FORWARD, ADJOINT))

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    allocate(dissipationTerm(grid%nGridPoints, solverOptions%nUnknowns))

    select case (mode)
    case (FORWARD)
       dissipationAmount = + solverOptions%dissipationAmount
    case (ADJOINT)
       dissipationAmount = - solverOptions%dissipationAmount
    end select

    do i = 1, nDimensions

       select case (mode)
       case (FORWARD)
          dissipationTerm = state%conservedVariables
       case (ADJOINT)
          dissipationTerm = state%adjointVariables
       end select

       call grid%dissipation(i)%apply(dissipationTerm, grid%localSize)

       if (.not. simulationFlags%compositeDissipation) then
          do j = 1, solverOptions%nUnknowns
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

end module RhsHelper
