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

    call endTiming("endDissipation")

  end subroutine addDissipation

end module RhsHelperImpl

subroutine computeRhsForward(simulationFlags, solverOptions, grid, state)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD

  ! <<< Private members >>>
  use RhsHelperImpl, only : addDissipation

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
  integer :: i, nDimensions
  SCALAR_TYPE, allocatable :: fluxes1(:,:,:), fluxes2(:,:,:)

  call startTiming("computeRhsForward")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  allocate(fluxes1(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))
  allocate(fluxes2(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))

  state%rightHandSide = 0.0_wp

  ! Compute Cartesian form of inviscid fluxes.
  call computeCartesianInvsicidFluxes(nDimensions, state%conservedVariables,                 &
       state%velocity, state%pressure(:,1), fluxes1)

  ! Compute Cartesian form of viscous fluxes if viscous terms are included and computed using
  ! repeated first derivatives.
  if (simulationFlags%viscosityOn .and. simulationFlags%repeatFirstDerivative) then
     call computeCartesianViscousFluxes(nDimensions, state%velocity,                         &
          state%stressTensor, state%heatFlux, fluxes2)
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

  call endTiming("computeRhsForward")

end subroutine computeRhsForward

subroutine computeRhsAdjoint(simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : ADJOINT

  ! <<< Private members >>>
  use RhsHelperImpl, only : addDissipation

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
  integer :: i, j, k, nDimensions
  SCALAR_TYPE, allocatable :: temp1(:,:,:), temp2(:,:),                                      &
       localFluxJacobian1(:,:), localFluxJacobian2(:,:), localConservedVariables(:),         &
       localVelocity(:), localMetricsAlongDirection1(:), localMetricsAlongDirection2(:),     &
       localStressTensor(:), localHeatFlux(:), localAdjointDiffusion(:,:)

  call startTiming("computeRhsAdjoint")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  allocate(temp1(grid%nGridPoints, solverOptions%nUnknowns, nDimensions))

  state%rightHandSide = 0.0_wp

  ! Partial derivatives of adjoint variables w.r.t. *computational* coordinates.
  do i = 1, nDimensions
     temp1(:,:,i) = state%adjointVariables
     call grid%adjointFirstDerivative(i)%apply(temp1(:,:,i), grid%localSize)
  end do

  allocate(localFluxJacobian1(solverOptions%nUnknowns, solverOptions%nUnknowns))
  allocate(localConservedVariables(solverOptions%nUnknowns))
  allocate(localVelocity(nDimensions))
  allocate(localMetricsAlongDirection1(nDimensions))

  if (simulationFlags%viscosityOn) then
     allocate(localFluxJacobian2(solverOptions%nUnknowns, solverOptions%nUnknowns))
     allocate(localStressTensor(nDimensions ** 2))
     allocate(localHeatFlux(nDimensions))
  end if

  do i = 1, grid%nGridPoints

     localConservedVariables = state%conservedVariables(i,:)
     localVelocity = state%velocity(i,:)
     if (simulationFlags%viscosityOn) then
        localStressTensor = state%stressTensor(i,:)
        localHeatFlux = state%heatFlux(i,:)
     end if

     do j = 1, nDimensions

        localMetricsAlongDirection1 = grid%metrics(i,1+nDimensions*(j-1):nDimensions*j)

        select case (nDimensions)
        case (1)
           call computeJacobianOfInviscidFlux1D(localConservedVariables,                     &
                localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                localFluxJacobian1, specificVolume = state%specificVolume(i,1),              &
                velocity = localVelocity, temperature = state%temperature(i,1))
        case (2)
           call computeJacobianOfInviscidFlux2D(localConservedVariables,                     &
                localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                localFluxJacobian1, specificVolume = state%specificVolume(i,1),              &
                velocity = localVelocity, temperature = state%temperature(i,1))
        case (3)
           call computeJacobianOfInviscidFlux3D(localConservedVariables,                     &
                localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                localFluxJacobian1, specificVolume = state%specificVolume(i,1),              &
                velocity = localVelocity, temperature = state%temperature(i,1))
        end select !... nDimensions

        if (simulationFlags%viscosityOn) then
           select case (nDimensions)
           case (1)
              call computeFirstPartialViscousJacobian1D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian2, specificVolume = state%specificVolume(i,1),           &
                   velocity = localVelocity, temperature = state%temperature(i,1))
           case (2)
              call computeFirstPartialViscousJacobian2D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian2, specificVolume = state%specificVolume(i,1),           &
                   velocity = localVelocity, temperature = state%temperature(i,1))
           case (3)
              call computeFirstPartialViscousJacobian3D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian2, specificVolume = state%specificVolume(i,1),           &
                   velocity = localVelocity, temperature = state%temperature(i,1))
           end select
           localFluxJacobian1 = localFluxJacobian1 - localFluxJacobian2
        end if

        state%rightHandSide(i,:) = state%rightHandSide(i,:) +                                &
             matmul(transpose(localFluxJacobian1), temp1(i,:,j))

     end do !... j = 1, nDimensions

  end do !... i = 1, grid%nGridPoints

  SAFE_DEALLOCATE(localConservedVariables)
  SAFE_DEALLOCATE(localFluxJacobian1)

  if (simulationFlags%viscosityOn) then

     SAFE_DEALLOCATE(localHeatFlux)
     SAFE_DEALLOCATE(localStressTensor)
     SAFE_DEALLOCATE(localFluxJacobian2)

     allocate(temp2(grid%nGridPoints, solverOptions%nUnknowns - 1))

     allocate(localMetricsAlongDirection2(nDimensions))
     allocate(localFluxJacobian2(solverOptions%nUnknowns - 1, solverOptions%nUnknowns - 1))
     allocate(localAdjointDiffusion(solverOptions%nUnknowns - 1, nDimensions))

     do i = 1, grid%nGridPoints

        localVelocity = state%velocity(i,:)
        localAdjointDiffusion = 0.0_wp

        do j = 1, nDimensions

           localMetricsAlongDirection1 = grid%metrics(i,1+nDimensions*(j-1):nDimensions*j)

           do k = 1, nDimensions

              localMetricsAlongDirection2 = grid%metrics(i,1+nDimensions*(k-1):nDimensions*k)

              select case (nDimensions)
              case (1)
                 call computeSecondPartialViscousJacobian1D(localVelocity,                   &
                      state%dynamicViscosity(i,1), state%secondCoefficientOfViscosity(i,1),  &
                      state%thermalDiffusivity(i,1), grid%jacobian(i,1),                     &
                      localMetricsAlongDirection1(1), localFluxJacobian2)
              case (2)
                 call computeSecondPartialViscousJacobian2D(localVelocity,                   &
                      state%dynamicViscosity(i,1), state%secondCoefficientOfViscosity(i,1),  &
                      state%thermalDiffusivity(i,1), grid%jacobian(i,1),                     &
                      localMetricsAlongDirection2, localMetricsAlongDirection1,              &
                      localFluxJacobian2)
              case (3)
                 call computeSecondPartialViscousJacobian3D(localVelocity,                   &
                      state%dynamicViscosity(i,1), state%secondCoefficientOfViscosity(i,1),  &
                      state%thermalDiffusivity(i,1), grid%jacobian(i,1),                     &
                      localMetricsAlongDirection2, localMetricsAlongDirection1,              &
                      localFluxJacobian2)
              end select !... nDimensions

              localAdjointDiffusion(:,j) = localAdjointDiffusion(:,j) +                      &
                   matmul(transpose(localFluxJacobian2), temp1(i,2:solverOptions%nUnknowns,k))

           end do !... k = 1, nDimensions

        end do !... j = 1, nDimensions

        do j = 1, nDimensions
           temp1(i,2:solverOptions%nUnknowns,j) = localAdjointDiffusion(:,j)
        end do

     end do !... i = 1, grid%nGridPoints

     do i = 1, nDimensions
        call grid%adjointFirstDerivative(i)%apply(temp1(:,2:solverOptions%nUnknowns,i),      &
             grid%localSize)
     end do
     temp2 = sum(temp1(:,2:solverOptions%nUnknowns,:), dim = 3)

     temp2(:,nDimensions+1) = solverOptions%ratioOfSpecificHeats *                           &
          state%specificVolume(:,1) * temp2(:,nDimensions+1)
     do i = 1, nDimensions
        temp2(:,i) = state%specificVolume(:,1) * temp2(:,i) -                                &
             state%velocity(:,i) * temp2(:,nDimensions+1)
     end do

     state%rightHandSide(:,2:solverOptions%nUnknowns) =                                      &
          state%rightHandSide(:,2:solverOptions%nUnknowns) - temp2
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

  call endTiming("computeRhsAdjoint")

end subroutine computeRhsAdjoint
