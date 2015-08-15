#include "config.h"

module BlockInterfacePatch_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  use Patch_mod, only : t_Patch

  implicit none
  private

  type, extends(t_Patch), public :: t_BlockInterfacePatch

     real(SCALAR_KIND) :: inviscidPenaltyAmount, viscousPenaltyAmount
     real(SCALAR_KIND), allocatable :: conservedVariablesR(:,:), adjointVariablesR(:,:),     &
          cartesianViscousFluxes(:,:,:), viscousFluxesL(:,:), viscousFluxesR(:,:),           &
          sendBuffer(:,:), receiveBuffer(:,:)

   contains

     procedure, pass :: setup
     procedure, pass :: cleanup
     procedure, pass :: updateRhs
     procedure, pass :: packSendBuffer
     procedure, pass :: unpackReceiveBuffer
     procedure, pass :: reshapeReceiveBuffer

  end type t_BlockInterfacePatch

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
    class(t_BlockInterfacePatch) :: this
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

    if (extent((direction-1)*2+1) /= extent((direction-1)*2+2)) then
       write(message, '(A)') "Patch '", trim(name),                                          &
            "' is not allowed to extend more than 1 grid point along normal direction!"
       call gracefulExit(grid%comm, message)
    end if

    write(key, '(A)') "patches/" // trim(name) // "/"

    ! Inviscid penalty amount.
    this%inviscidPenaltyAmount = getOption("interfaces/inviscid_penalty_amount", 1.0_wp)
    this%inviscidPenaltyAmount = getOption(trim(key) // "inviscid_penalty_amount",           &
         this%inviscidPenaltyAmount)
    this%inviscidPenaltyAmount = sign(this%inviscidPenaltyAmount,                            &
         real(this%normalDirection, wp))
    this%inviscidPenaltyAmount = this%inviscidPenaltyAmount /                                &
         grid%firstDerivative(direction)%normBoundary(1)

    ! Viscous penalty amount.
    if (simulationFlags%viscosityOn) then
       this%viscousPenaltyAmount = getOption("interfaces/viscous_penalty_amount", 0.5_wp)
       this%viscousPenaltyAmount = getOption(trim(key) // "viscous_penalty_amount",          &
            this%viscousPenaltyAmount)
       this%viscousPenaltyAmount = sign(this%viscousPenaltyAmount,                           &
            real(this%normalDirection, wp))
       this%viscousPenaltyAmount = this%viscousPenaltyAmount /                               &
            grid%firstDerivative(direction)%normBoundary(1)
    else
       this%viscousPenaltyAmount = 0.0_wp
    end if

    if (this%nPatchPoints > 0) then

       allocate(this%conservedVariablesR(this%nPatchPoints, nUnknowns))
       if (simulationFlags%viscosityOn)                                                      &
            allocate(this%cartesianViscousFluxes(this%nPatchPoints, nUnknowns, nDimensions))

       if (.not. simulationFlags%predictionOnly) then
          allocate(this%adjointVariablesR(this%nPatchPoints, nUnknowns))
          if (simulationFlags%viscosityOn) then
             allocate(this%viscousFluxesL(this%nPatchPoints, nUnknowns - 1))
             allocate(this%viscousFluxesR(this%nPatchPoints, nUnknowns - 1))
          end if
       end if

    end if

  end subroutine setup

  subroutine cleanup(this)

    implicit none

    ! <<< Arguments >>>
    class(t_BlockInterfacePatch) :: this

    call this%cleanupBase()

    SAFE_DEALLOCATE(this%conservedVariablesR)
    SAFE_DEALLOCATE(this%adjointVariablesR)
    SAFE_DEALLOCATE(this%cartesianViscousFluxes)
    SAFE_DEALLOCATE(this%viscousFluxesL)
    SAFE_DEALLOCATE(this%viscousFluxesR)
    SAFE_DEALLOCATE(this%sendBuffer)
    SAFE_DEALLOCATE(this%receiveBuffer)

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
    class(t_BlockInterfacePatch) :: this
    integer, intent(in) :: mode
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid), intent(in) :: grid
    class(t_State) :: state

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, l, nDimensions, nUnknowns, direction,                                &
         incomingDirection, gridIndex, patchIndex, depthIndex, depthBlockSize(3)
    real(SCALAR_KIND), allocatable :: fluxJacobian(:,:), deltaFluxJacobian(:,:,:),           &
         temp1(:,:), temp2(:,:)

    assert_key(mode, (FORWARD, ADJOINT))

    call startTiming("Interface penalty")

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
    if (mode == ADJOINT .and. .not. simulationFlags%useContinuousAdjoint)                    &
         allocate(deltaFluxJacobian(nUnknowns, nUnknowns, nUnknowns))

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

                call computeIncomingInviscidJacobian(nDimensions, gridIndex, patchIndex,     &
                     state%conservedVariables, this%conservedVariablesR,                     &
                     grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction),      &
                     solverOptions%ratioOfSpecificHeats, incomingDirection, fluxJacobian)

                state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -        &
                     this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *              &
                     matmul(fluxJacobian, state%conservedVariables(gridIndex,:) -            &
                     this%conservedVariablesR(patchIndex,:))

                if (simulationFlags%viscosityOn) then
                   state%rightHandSide(gridIndex,2:nUnknowns) =                              &
                        state%rightHandSide(gridIndex,2:nUnknowns) -                         &
                        this%viscousPenaltyAmount * grid%jacobian(gridIndex, 1) *            &
                        (this%viscousFluxesL(patchIndex,:) -                                 &
                        this%viscousFluxesR(patchIndex,:))
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

                if (simulationFlags%useContinuousAdjoint) then
                   call computeIncomingInviscidJacobian(nDimensions, gridIndex, patchIndex,  &
                        state%conservedVariables, this%conservedVariablesR,                  &
                        grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction),   &
                        solverOptions%ratioOfSpecificHeats, +incomingDirection,              &
                        fluxJacobian)
                else
                   call computeIncomingInviscidJacobian(nDimensions, gridIndex, patchIndex,  &
                        state%conservedVariables, this%conservedVariablesR,                  &
                        grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction),   &
                        solverOptions%ratioOfSpecificHeats, +incomingDirection,              &
                        fluxJacobian, deltaFluxJacobian)
                end if

                state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +        &
                     sign(this%inviscidPenaltyAmount, real(incomingDirection, wp)) *         &
                     grid%jacobian(gridIndex, 1) * matmul(transpose(fluxJacobian),           &
                     state%adjointVariables(gridIndex,:))

                if (.not. simulationFlags%useContinuousAdjoint) then
                   do l = 1, nUnknowns
                      state%rightHandSide(gridIndex,l) = state%rightHandSide(gridIndex,l) +  &
                           this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *        &
                           sum(matmul(deltaFluxJacobian(:,:,l),                              &
                           state%conservedVariables(gridIndex,:) -                           &
                           this%conservedVariablesR(patchIndex,:)) *                         &
                           state%adjointVariables(gridIndex,:))
                   end do
                end if

                if (simulationFlags%useContinuousAdjoint) then
                   call computeIncomingInviscidJacobian(nDimensions, gridIndex, patchIndex,  &
                        state%conservedVariables, this%conservedVariablesR,                  &
                        grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction),   &
                        solverOptions%ratioOfSpecificHeats, -incomingDirection,              &
                        fluxJacobian)
                else
                   call computeIncomingInviscidJacobian(nDimensions, gridIndex, patchIndex,  &
                        state%conservedVariables, this%conservedVariablesR,                  &
                        grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction),   &
                        solverOptions%ratioOfSpecificHeats, -incomingDirection,              &
                        fluxJacobian, deltaFluxJacobian)
                end if

                state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +        &
                     sign(this%inviscidPenaltyAmount, real(incomingDirection, wp)) *         &
                     grid%jacobian(gridIndex, 1) * matmul(transpose(fluxJacobian),           &
                     this%adjointVariablesR(patchIndex,:))

                if (.not. simulationFlags%useContinuousAdjoint) then
                   do l = 1, nUnknowns
                      state%rightHandSide(gridIndex,l) = state%rightHandSide(gridIndex,l) +  &
                           this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *        &
                           sum(matmul(deltaFluxJacobian(:,:,l),                              &
                           state%conservedVariables(gridIndex,:) -                           &
                           this%conservedVariablesR(patchIndex,:)) *                         &
                           this%adjointVariablesR(patchIndex,:))
                   end do
                end if

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
                        state%adjointVariables(gridIndex,:) +                                &
                        this%adjointVariablesR(patchIndex,:))

                end if

             end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
          end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
       end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

    end select

    SAFE_DEALLOCATE(fluxJacobian)
    SAFE_DEALLOCATE(deltaFluxJacobian)

    if (mode == ADJOINT .and. simulationFlags%viscosityOn) then

       allocate(fluxJacobian(nUnknowns - 1, nUnknowns - 1))

       depthBlockSize = this%localSize
       depthBlockSize(direction) = grid%adjointFirstDerivative(direction)%boundaryDepth

       if (this%nPatchPoints > 0) then
          allocate(temp1(product(depthBlockSize), nUnknowns - 1), source = 0.0_wp)
          allocate(temp2(product(depthBlockSize), nUnknowns), source = 0.0_wp)
       end if

       do l = 1, nDimensions

          do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
             do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
                do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                   gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *              &
                        (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                &
                        (k - 1 - this%gridOffset(3)))
                   if (grid%iblank(gridIndex) == 0) cycle
                   depthIndex = i - this%offset(1) + depthBlockSize(1) *                     &
                        (j - 1 - this%offset(2) + depthBlockSize(2) *                        &
                        (k - 1 - this%offset(3)))

                   call computeSecondPartialViscousJacobian(nDimensions, gridIndex,          &
                        state%velocity, state%dynamicViscosity(:,1),                         &
                        state%secondCoefficientOfViscosity(:,1),                             &
                        state%thermalDiffusivity(:,1), grid%jacobian(:,1),                   &
                        grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction),   &
                        grid%metrics(:,1+nDimensions*(l-1):nDimensions*l), fluxJacobian)

                   temp1(depthIndex,:) = matmul(transpose(fluxJacobian),                     &
                        state%adjointVariables(gridIndex,2:nUnknowns) +                      &
                        this%adjointVariablesR(patchIndex,2:nUnknowns))

                end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
             end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
          end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

          call grid%adjointFirstDerivative(l)%projectOnBoundaryAndApply(temp1,               &
               depthBlockSize, this%normalDirection)
          if (this%nPatchPoints > 0) temp2(:,2:nUnknowns) = temp2(:,2:nUnknowns) - temp1

       end do !... l = 1, nDimensions

       do k = this%offset(3) + 1, this%offset(3) + depthBlockSize(3)
          do j = this%offset(2) + 1, this%offset(2) + depthBlockSize(2)
             do i = this%offset(1) + 1, this%offset(1) + depthBlockSize(1)
                gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                 &
                     (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                   &
                     (k - 1 - this%gridOffset(3)))
                if (grid%iblank(gridIndex) == 0) cycle
                depthIndex = i - this%offset(1) + depthBlockSize(1) *                        &
                     (j - 1 - this%offset(2) + depthBlockSize(2) *                           &
                     (k - 1 - this%offset(3)))

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

    call endTiming("Interface penalty")

  end subroutine updateRhs

  subroutine packSendBuffer(this, mode, simulationFlags, solverOptions, grid, state)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod, only : t_State
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : FORWARD, ADJOINT

    implicit none

    ! <<< Arguments >>>
    class(t_BlockInterfacePatch) :: this
    integer, intent(in) :: mode
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid), intent(in) :: grid
    class(t_State), intent(in) :: state

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, nDimensions, nUnknowns, nExchangedComponents, direction, procRank, ierror
    real(SCALAR_KIND), allocatable :: sendBuffer(:,:), metricsAlongNormalDirection(:,:)

    if (this%comm == MPI_COMM_NULL) return

    assert(this%nPatchPoints > 0)

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    direction = abs(this%normalDirection)
    assert(direction >= 1 .and. direction <= nDimensions)

    nUnknowns = solverOptions%nUnknowns

    nExchangedComponents = nUnknowns
    if (mode == ADJOINT .and. .not. simulationFlags%predictionOnly)                          &
         nExchangedComponents = nExchangedComponents + nUnknowns
    if (mode == FORWARD .and. simulationFlags%viscosityOn)                                   &
         nExchangedComponents = nExchangedComponents + nUnknowns - 1

    call MPI_Comm_rank(this%comm, procRank, ierror)
    if (procRank == 0)                                                                       &
         allocate(this%sendBuffer(product(this%globalSize), nExchangedComponents))
    call MPI_Barrier(this%comm, ierror)

    if (mode == FORWARD .and. simulationFlags%viscosityOn) then

       allocate(metricsAlongNormalDirection(this%nPatchPoints, nDimensions))
       call this%collect(grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction),  &
            metricsAlongNormalDirection)

       this%viscousFluxesL = 0.0_wp
       do j = 1, nDimensions
          do i = 2, nUnknowns
             this%viscousFluxesL(:,i-1) = this%viscousFluxesL(:,i-1) +                       &
                  this%cartesianViscousFluxes(:,i,j) * metricsAlongNormalDirection(:,j)
          end do
       end do

    end if

    allocate(sendBuffer(this%nPatchPoints, nExchangedComponents))

    call this%collect(state%conservedVariables, sendBuffer(:,1:nUnknowns))
    if (mode == ADJOINT .and. .not. simulationFlags%predictionOnly)                          &
         call this%collect(state%adjointVariables, sendBuffer(:,nUnknowns+1:))
    if (mode == FORWARD .and. simulationFlags%viscosityOn)                                   &
         sendBuffer(:,nUnknowns+1:) = this%viscousFluxesL

    call this%gather(sendBuffer, this%sendBuffer)

    SAFE_DEALLOCATE(sendBuffer)
    SAFE_DEALLOCATE(metricsAlongNormalDirection)

  end subroutine packSendBuffer

  subroutine unpackReceiveBuffer(this, mode, simulationFlags, solverOptions)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : FORWARD, ADJOINT

    implicit none

    ! <<< Arguments >>>
    class(t_BlockInterfacePatch) :: this
    integer, intent(in) :: mode
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: nUnknowns, nExchangedComponents
    real(SCALAR_KIND), allocatable :: receiveBuffer(:,:)

    if (this%comm == MPI_COMM_NULL) return

    assert(this%nPatchPoints > 0)

    nUnknowns = solverOptions%nUnknowns

    nExchangedComponents = nUnknowns
    if (mode == ADJOINT .and. .not. simulationFlags%predictionOnly)                          &
         nExchangedComponents = nExchangedComponents + nUnknowns
    if (mode == FORWARD .and. simulationFlags%viscosityOn)                                   &
         nExchangedComponents = nExchangedComponents + nUnknowns - 1

    allocate(receiveBuffer(this%nPatchPoints, nExchangedComponents))

    call this%scatter(this%receiveBuffer, receiveBuffer)

    this%conservedVariablesR = receiveBuffer(:,1:nUnknowns)
    if (mode == ADJOINT .and. .not. simulationFlags%predictionOnly)                          &
         this%adjointVariablesR = receiveBuffer(:,nUnknowns+1:)
    if (mode == FORWARD .and. simulationFlags%viscosityOn)                                   &
         this%viscousFluxesR = receiveBuffer(:,nUnknowns+1:)

    SAFE_DEALLOCATE(receiveBuffer)

  end subroutine unpackReceiveBuffer

  subroutine reshapeReceiveBuffer(this, indexReordering)

    ! <<< External modules >>>
    use MPI

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    implicit none

    ! <<< Arguments >>>
    class(t_BlockInterfacePatch) :: this
    integer, intent(in) :: indexReordering(3)

    ! <<< Local variables >>>
    integer :: i, j, k, l, order(3), globalSize(3), nComponents, ierror

    if (this%comm == MPI_COMM_NULL) return

    assert(indexReordering(3) == 3)
    if (indexReordering(1) == 1 .and. indexReordering(2) == 2) return

    call startTiming("Reshape interface buffer")

    if (allocated(this%receiveBuffer)) then

       assert(allocated(this%sendBuffer))
       assert(all(this%globalSize > 0))
       assert(size(this%sendBuffer, 1) == product(this%globalSize))
       assert(size(this%sendBuffer, 2) > 0)
       assert(all(shape(this%receiveBuffer) == shape(this%sendBuffer)))

       nComponents = size(this%receiveBuffer, 2)
       globalSize = this%globalSize

       order = indexReordering

       if (abs(order(1)) == 2 .and. abs(order(2)) == 1) then

          do l = 1, nComponents
             do k = 1, globalSize(3)
                do j = 1, globalSize(2)
                   do i = 1, globalSize(1)
                      this%sendBuffer(i + globalSize(1) * (j - 1 +                           &
                           globalSize(2) * (k - 1)), l) =                                    &
                           this%receiveBuffer(j + globalSize(2) * (i - 1 +                   &
                           globalSize(1) * (k - 1)), l)
                   end do
                end do
             end do
          end do

          this%receiveBuffer = this%sendBuffer

          i = order(2)
          order(2) = order(1)
          order(1) = i

       end if

       if (order(1) == -1 .and. order(2) == 2) then

          do l = 1, nComponents
             do k = 1, globalSize(3)
                do j = 1, globalSize(2)
                   do i = 1, globalSize(1)
                      this%sendBuffer(i + globalSize(1) * (j - 1 +                           &
                           globalSize(2) * (k - 1)), l) =                                    &
                           this%receiveBuffer(globalSize(1) + 1 - i +                        &
                           globalSize(1) * (j - 1 + globalSize(2) * (k - 1)), l)
                   end do
                end do
             end do
          end do

          this%receiveBuffer = this%sendBuffer

       else if (order(1) == 1 .and. order(2) == -2) then

          do l = 1, nComponents
             do k = 1, globalSize(3)
                do j = 1, globalSize(2)
                   do i = 1, globalSize(1)
                      this%sendBuffer(i + globalSize(1) * (j - 1 +                           &
                           globalSize(2) * (k - 1)), l) =                                    &
                           this%receiveBuffer(i + globalSize(1) * (globalSize(2) - j +       &
                           globalSize(2) * (k - 1)), l)
                   end do
                end do
             end do
          end do

          this%receiveBuffer = this%sendBuffer

       else if (order(1) == -1 .and. order(2) == -2) then

          do l = 1, nComponents
             do k = 1, globalSize(3)
                do j = 1, globalSize(2)
                   do i = 1, globalSize(1)
                      this%sendBuffer(i + globalSize(1) * (j - 1 +                           &
                           globalSize(2) * (k - 1)), l) =                                    &
                           this%receiveBuffer(globalSize(1) + 1 - i +                        &
                           globalSize(1) * (globalSize(2) - j + globalSize(2) * (k - 1)), l)
                   end do
                end do
             end do
          end do

          this%receiveBuffer = this%sendBuffer

       end if

    end if

    call MPI_Barrier(this%comm, ierror)

    call endTiming("Reshape interface buffer")

  end subroutine reshapeReceiveBuffer

end module BlockInterfacePatch_mod
