#include "config.h"

module InflowOutflowPatch_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  use Patch_mod, only : t_Patch

  implicit none
  private

  type, extends(t_Patch), public :: t_InflowOutflowPatch

     real(SCALAR_KIND) :: inviscidPenaltyAmount, viscousPenaltyAmount
     real(SCALAR_KIND), allocatable :: viscousFluxes(:,:,:), targetViscousFluxes(:,:,:)

   contains

     procedure, pass :: setup
     procedure, pass :: cleanup
     procedure, pass :: updateRhs

  end type t_InflowOutflowPatch

contains

  subroutine setup(this, name, comm, grid, state, extent,                                    &
       normalDirection, simulationFlags, solverOptions)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod, only : t_State
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Internal modules >>>
    use CNSHelper, only : computeCartesianViscousFluxes
    use InputHelper, only : getOption
    use ErrorHandler, only : gracefulExit

    implicit none

    ! <<< Arguments >>>
    class(t_InflowOutflowPatch) :: this
    character(len = *), intent(in) :: name
    integer, intent(in) :: comm
    class(t_Grid), intent(in) :: grid
    class(t_State), intent(in) :: state
    integer, intent(in) :: extent(6), normalDirection
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    character(len = STRING_LENGTH) :: key, message
    integer :: nDimensions, nUnknowns, direction
    real(SCALAR_KIND), allocatable :: velocity(:,:), stressTensor(:,:), heatFlux(:,:)

    call this%cleanup()
    call this%setupBase(name, comm, grid, extent, normalDirection)

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns >= nDimensions + 2)

    direction = abs(this%normalDirection)

    if (direction <= 0 .or. direction > nDimensions) then
       write(message, '(3A)') "Invalid normal direction for patch '", trim(name), "'!"
       call gracefulExit(grid%comm, message)
    end if

    if (.not. simulationFlags%useTargetState) then
       write(message, '(A)') "Patch '", trim(name),                                          &
            "' requires a target state for enforcing far-field boundary conditions!"
       call gracefulExit(grid%comm, message)
    end if

    if (extent((direction-1)*2+1) /= extent((direction-1)*2+2)) then
       write(message, '(A)') "Patch '", trim(name),                                          &
            "' is not allowed to extend more than 1 grid point along normal direction!"
       call gracefulExit(grid%comm, message)
    end if

    write(key, '(A)') "patches/" // trim(name) // "/"

    ! Inviscid penalty amount.
    this%inviscidPenaltyAmount = getOption("defaults/inviscid_penalty_amount", 1.0_wp)
    this%inviscidPenaltyAmount = getOption(trim(key) // "inviscid_penalty_amount",           &
         this%inviscidPenaltyAmount)
    this%inviscidPenaltyAmount = sign(this%inviscidPenaltyAmount,                            &
         real(this%normalDirection, wp))
    this%inviscidPenaltyAmount = this%inviscidPenaltyAmount /                                &
         grid%firstDerivative(direction)%normBoundary(1)

    ! Viscous penalty amount.
    if (simulationFlags%viscosityOn) then
       this%viscousPenaltyAmount = getOption("defaults/viscous_penalty_amount", 1.0_wp)
       this%viscousPenaltyAmount = getOption(trim(key) // "viscous_penalty_amount",          &
            this%viscousPenaltyAmount)
       this%viscousPenaltyAmount = sign(this%viscousPenaltyAmount,                           &
            real(this%normalDirection, wp))
       this%viscousPenaltyAmount = this%viscousPenaltyAmount /                               &
            grid%firstDerivative(direction)%normBoundary(1)
    else
       this%viscousPenaltyAmount = 0.0_wp
    end if

    if (simulationFlags%viscosityOn .and. allocated(state%targetState)) then

       call state%update(grid, simulationFlags, solverOptions, useTargetState = .true.)

       if (this%nPatchPoints > 0) then

          allocate(this%viscousFluxes(this%nPatchPoints, nUnknowns, nDimensions))
          allocate(this%targetViscousFluxes(this%nPatchPoints, nUnknowns, nDimensions))

          allocate(velocity(this%nPatchPoints, nDimensions))
          allocate(stressTensor(this%nPatchPoints, nDimensions ** 2))
          allocate(heatFlux(this%nPatchPoints, nDimensions))

          call this%collect(state%velocity, velocity)
          call this%collect(state%stressTensor, stressTensor)
          call this%collect(state%heatFlux, heatFlux)
          call computeCartesianViscousFluxes(nDimensions, velocity,                          &
               stressTensor, heatFlux, this%targetViscousFluxes)

          SAFE_DEALLOCATE(heatFlux)
          SAFE_DEALLOCATE(stressTensor)
          SAFE_DEALLOCATE(velocity)

       end if

    end if

  end subroutine setup

  subroutine cleanup(this)

    implicit none

    ! <<< Arguments >>>
    class(t_InflowOutflowPatch) :: this

    call this%cleanupBase()

    SAFE_DEALLOCATE(this%viscousFluxes)
    SAFE_DEALLOCATE(this%targetViscousFluxes)

  end subroutine cleanup

  subroutine updateRhs(this, mode, simulationFlags, solverOptions, grid, state)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod, only : t_State
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : FORWARD, ADJOINT

    ! <<< Internal modules >>>
    use CNSHelper
    use MPITimingsHelper, only : startTiming, endTiming

    implicit none

    ! <<< Arguments >>>
    class(t_InflowOutflowPatch) :: this
    integer, intent(in) :: mode
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid), intent(in) :: grid
    class(t_State) :: state

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, l, nDimensions, nUnknowns, direction, incomingDirection,             &
         gridIndex, patchIndex, depthIndex, depthBlockOffset(3), depthBlockSize(3)
    real(SCALAR_KIND), allocatable :: fluxJacobian(:,:), temp1(:,:), temp2(:,:)

    assert_key(mode, (FORWARD, ADJOINT))
    assert(allocated(state%targetState))

    call startTiming("addInflowOutflowPenalty")

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    direction = abs(this%normalDirection)
    assert(direction >= 1 .and. direction <= nDimensions)

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns >= nDimensions + 2)

    if (mode == ADJOINT .and. simulationFlags%useContinuousAdjoint) then
       incomingDirection = -this%normalDirection
    else
       incomingDirection = +this%normalDirection
    end if

    allocate(fluxJacobian(nUnknowns, nUnknowns))

    select case (mode)

    case (FORWARD)

       do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
          do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
             do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                 &
                     (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                   &
                     (k - 1 - this%gridOffset(3)))
                if (grid%iblank(gridIndex) == 0) cycle
                patchIndex = i - this%offset(1) + this%localSize(1) *                        &
                     (j - 1 - this%offset(2) + this%localSize(2) *                           &
                     (k - 1 - this%offset(3)))

                call computeIncomingInviscidJacobian(nDimensions, gridIndex,                 &
                     state%targetState,                                                      &
                     grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction),      &
                     solverOptions%ratioOfSpecificHeats, incomingDirection, fluxJacobian)

                state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -        &
                     this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *              &
                     matmul(fluxJacobian, state%conservedVariables(gridIndex,:) -            &
                     state%targetState(gridIndex,:))

                if (simulationFlags%viscosityOn) then
                   state%rightHandSide(gridIndex,2:nUnknowns) =                              &
                        state%rightHandSide(gridIndex,2:nUnknowns) +                         &
                        this%viscousPenaltyAmount * grid%jacobian(gridIndex, 1) *            &
                        matmul(this%viscousFluxes(patchIndex,2:nUnknowns,:) -                &
                        this%targetViscousFluxes(patchIndex,2:nUnknowns,:),                  &
                        grid%metrics(gridIndex,                                              &
                        1+nDimensions*(direction-1):nDimensions*direction))
                end if

             end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
          end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
       end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

    case (ADJOINT)

       do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
          do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
             do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                 &
                     (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                   &
                     (k - 1 - this%gridOffset(3)))
                if (grid%iblank(gridIndex) == 0) cycle
                patchIndex = i - this%offset(1) + this%localSize(1) *                        &
                     (j - 1 - this%offset(2) + this%localSize(2) *                           &
                     (k - 1 - this%offset(3)))

                call computeIncomingInviscidJacobian(nDimensions, gridIndex,                 &
                     state%targetState,                                                      &
                     grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction),      &
                     solverOptions%ratioOfSpecificHeats, incomingDirection, fluxJacobian)

                state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +        &
                     sign(this%inviscidPenaltyAmount, real(incomingDirection, wp)) *         &
                     grid%jacobian(gridIndex, 1) * matmul(transpose(fluxJacobian),           &
                     state%adjointVariables(gridIndex,:))

                if (simulationFlags%viscosityOn) then

                   call computeFirstPartialViscousJacobian(nDimensions, gridIndex,           &
                        state%conservedVariables,                                            &
                        grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction),   &
                        state%stressTensor, state%heatFlux, solverOptions%powerLawExponent,  &
                        solverOptions%ratioOfSpecificHeats, fluxJacobian,                    &
                        state%specificVolume(:,1), state%velocity, state%temperature(:,1))

                   state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +     &
                        sign(this%viscousPenaltyAmount, real(incomingDirection, wp)) *       &
                        grid%jacobian(gridIndex, 1) * matmul(transpose(fluxJacobian),        &
                        state%adjointVariables(gridIndex,:))

                end if

             end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
          end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
       end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

    end select

    SAFE_DEALLOCATE(fluxJacobian)

    if (mode == ADJOINT .and. simulationFlags%viscosityOn .and. this%nPatchPoints > 0) then

       allocate(fluxJacobian(nUnknowns - 1, nUnknowns - 1))

       depthBlockSize = this%localSize
       depthBlockSize(direction) = grid%adjointFirstDerivative(direction)%boundaryDepth

       depthBlockOffset = this%offset
       if (this%normalDirection < 0) depthBlockOffset(direction) =                           &
            depthBlockOffset(direction) + 1 - depthBlockSize(direction)

       allocate(temp1(product(depthBlockSize), nUnknowns - 1), source = 0.0_wp)
       allocate(temp2(product(depthBlockSize), nUnknowns), source = 0.0_wp)

       do l = 1, nDimensions

          do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
             do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
                do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                   gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *              &
                        (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                &
                        (k - 1 - this%gridOffset(3)))
                   if (grid%iblank(gridIndex) == 0) cycle
                   depthIndex = i - depthBlockOffset(1) + depthBlockSize(1) *                &
                        (j - 1 - depthBlockOffset(2) + depthBlockSize(2) *                   &
                        (k - 1 - depthBlockOffset(3)))

                   call computeSecondPartialViscousJacobian(nDimensions, gridIndex,          &
                        state%velocity, state%dynamicViscosity(:,1),                         &
                        state%secondCoefficientOfViscosity(:,1),                             &
                        state%thermalDiffusivity(:,1), grid%jacobian(:,1),                   &
                        grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction),   &
                        grid%metrics(:,1+nDimensions*(l-1):nDimensions*l), fluxJacobian)

                   temp1(depthIndex,:) = matmul(transpose(fluxJacobian),                     &
                        state%adjointVariables(gridIndex,2:nUnknowns))

                end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
             end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
          end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

          call grid%adjointFirstDerivative(l)%projectOnBoundaryAndApply(temp1,               &
               depthBlockSize, this%normalDirection)
          temp2(:,2:nUnknowns) = temp2(:,2:nUnknowns) - temp1

       end do !... l = 1, nDimensions

       do k = depthBlockOffset(3) + 1, depthBlockOffset(3) + depthBlockSize(3)
          do j = depthBlockOffset(2) + 1, depthBlockOffset(2) + depthBlockSize(2)
             do i = depthBlockOffset(1) + 1, depthBlockOffset(1) + depthBlockSize(1)
                gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                 &
                     (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                   &
                     (k - 1 - this%gridOffset(3)))
                if (grid%iblank(gridIndex) == 0) cycle
                depthIndex = i - depthBlockOffset(1) + depthBlockSize(1) *                   &
                     (j - 1 - depthBlockOffset(2) + depthBlockSize(2) *                      &
                     (k - 1 - depthBlockOffset(3)))

                temp2(depthIndex,nDimensions+2) = solverOptions%ratioOfSpecificHeats *       &
                     state%specificVolume(gridIndex,1) * temp2(depthIndex,nDimensions+2)
                temp2(depthIndex,2:nDimensions+1) = state%specificVolume(gridIndex,1) *      &
                     temp2(depthIndex,2:nDimensions+1) - state%velocity(gridIndex,:) *       &
                     temp2(depthIndex,nDimensions+2)
                temp2(depthIndex,1) = - state%specificVolume(gridIndex,1) *                  &
                     state%conservedVariables(gridIndex,nDimensions+2) *                     &
                     temp2(depthIndex,nDimensions+2) + sum(state%velocity(gridIndex,:) *     &
                     temp2(depthIndex,2:nDimensions+1))

                state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +        &
                     sign(this%viscousPenaltyAmount, real(incomingDirection, wp)) *          &
                     grid%jacobian(gridIndex, 1) * temp2(depthIndex,:)

             end do !... i = this%offset(1) + 1, this%offset(1) + depthBlockSize(1)
          end do !... j = this%offset(2) + 1, this%offset(2) + depthBlockSize(2)
       end do !... k = this%offset(3) + 1, this%offset(3) + depthBlockSize(3)

       SAFE_DEALLOCATE(temp2)
       SAFE_DEALLOCATE(temp1)
       SAFE_DEALLOCATE(fluxJacobian)

    end if

    call endTiming("addInflowOutflowPenalty")

  end subroutine updateRhs

end module InflowOutflowPatch_mod
