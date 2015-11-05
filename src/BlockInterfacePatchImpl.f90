#include "config.h"

subroutine setupBlockInterfacePatch(this, index, comm, patchDescriptor,                      &
     grid, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_BlockInterfacePatch) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: key
  integer :: nExchangedVariables, procRank, ierror

  call this%cleanup()
  call this%setupBase(index, comm, patchDescriptor, grid, simulationFlags, solverOptions)

  if (this%comm /= MPI_COMM_NULL) then

     assert(this%nPatchPoints > 0)

     nExchangedVariables = solverOptions%nUnknowns

     if (simulationFlags%viscosityOn) then
        allocate(this%viscousFluxes(this%nPatchPoints, solverOptions%nUnknowns - 1))
        nExchangedVariables = nExchangedVariables + solverOptions%nUnknowns - 1
     end if

     call MPI_Comm_rank(this%comm, procRank, ierror)
     if (procRank == 0) then
        allocate(this%sendBuffer(product(this%globalSize), nExchangedVariables))
        allocate(this%receiveBuffer(product(this%globalSize), nExchangedVariables))
     end if

  end if

  write(key, '(A)') "patches/" // trim(patchDescriptor%name) // "/"

  ! Inviscid penalty amount.
  this%inviscidPenaltyAmount = getOption(trim(key) //                                        &
       "inviscid_penalty_amount", 1.0_wp) !... default value => dual-consistent.
  this%inviscidPenaltyAmount = sign(this%inviscidPenaltyAmount,                              &
       real(this%normalDirection, wp))
  this%inviscidPenaltyAmount = this%inviscidPenaltyAmount /                                  &
       grid%firstDerivative(abs(this%normalDirection))%normBoundary(1)

  ! Viscous penalty amount.
  if (simulationFlags%viscosityOn) then
     this%viscousPenaltyAmount = getOption(trim(key) //                                      &
          "viscous_penalty_amount", 0.5_wp)
     this%viscousPenaltyAmount = sign(this%viscousPenaltyAmount,                             &
          real(this%normalDirection, wp))
     this%viscousPenaltyAmount = this%viscousPenaltyAmount /                                 &
          grid%firstDerivative(abs(this%normalDirection))%normBoundary(1)
  else
     this%viscousPenaltyAmount = 0.0_wp
  end if

end subroutine setupBlockInterfacePatch

subroutine cleanupBlockInterfacePatch(this)

  ! <<< Derived types >>>
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  implicit none

  ! <<< Arguments >>>
  class(t_BlockInterfacePatch) :: this

  call this%cleanupBase()

  SAFE_DEALLOCATE(this%viscousFluxes)
  SAFE_DEALLOCATE(this%sendBuffer)
  SAFE_DEALLOCATE(this%receiveBuffer)

end subroutine cleanupBlockInterfacePatch

subroutine addBlockInterfacePenalty(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

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
  integer :: i, j, k, nDimensions, nUnknowns, nSpecies, direction,                           &
       incomingDirection, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: receivedData(:,:), interfaceConservedVariables(:,:),           &
       interfaceAdjointVariables(:,:), interfaceViscousFluxes(:,:),                          &
       localConservedVariables(:), localInterfaceConservedVariables(:), localRoeAverage(:),  &
       localMetricsAlongNormalDirection(:), incomingJacobianOfInviscidFlux(:,:)

  assert_key(mode, (FORWARD, ADJOINT))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  call startTiming("addBlockInterfacePenalty")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nSpecies = solverOptions%nSpecies
  assert(nSpecies >= 0)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2 + nSpecies)

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  if (this%comm /= MPI_COMM_NULL) then

     allocate(interfaceConservedVariables(this%nPatchPoints, nUnknowns))
     allocate(interfaceViscousFluxes(this%nPatchPoints, nUnknowns - 1))

     if (mode == ADJOINT)                                                                    &
          allocate(interfaceAdjointVariables(this%nPatchPoints, nUnknowns))

     if (simulationFlags%viscosityOn) then
        if (mode == ADJOINT) then
           allocate(receivedData(this%nPatchPoints, 3 * nUnknowns - 1))
        else
           allocate(receivedData(this%nPatchPoints, 2 * nUnknowns - 1))
        end if
     else
        if (mode == ADJOINT) then
           allocate(receivedData(this%nPatchPoints, 2 * nUnknowns))
        else
           allocate(receivedData(this%nPatchPoints, nUnknowns))
        end if
     end if

  end if

  call this%scatterData(this%receiveBuffer, receivedData)

  interfaceConservedVariables = receivedData(:,1:nUnknowns)
  if (mode == ADJOINT) interfaceAdjointVariables = receivedData(:,nUnknowns+1:2*nUnknowns)
  if (simulationFlags%viscosityOn) then
     if (mode == ADJOINT) then
        interfaceViscousFluxes = receivedData(:,2*nUnknowns+1:)
     else
        interfaceViscousFluxes = receivedData(:,nUnknowns+1:)
     end if
  end if

  SAFE_DEALLOCATE(receivedData)

  if (mode == ADJOINT .and. simulationFlags%useContinuousAdjoint) then
     incomingDirection = -this%normalDirection
  else
     incomingDirection = +this%normalDirection
  end if

  allocate(localConservedVariables(nUnknowns))
  allocate(localInterfaceConservedVariables(nUnknowns))
  allocate(localRoeAverage(nUnknowns))
  allocate(localMetricsAlongNormalDirection(nDimensions))
  allocate(incomingJacobianOfInviscidFlux(nUnknowns, nUnknowns))

  do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
     do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
        do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
           gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                      &
                (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                        &
                (k - 1 - this%gridOffset(3)))
           if (grid%iblank(gridIndex) == 0) cycle
           patchIndex = i - this%offset(1) + this%localSize(1) *                             &
                (j - 1 - this%offset(2) + this%localSize(2) *                                &
                (k - 1 - this%offset(3)))

           localConservedVariables = state%conservedVariables(gridIndex,:)
           localInterfaceConservedVariables = interfaceConservedVariables(patchIndex,:)
           localMetricsAlongNormalDirection =                                                &
                grid%metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

           call computeRoeAverage(nDimensions, localConservedVariables,                      &
                localInterfaceConservedVariables, solverOptions%ratioOfSpecificHeats,        &
                localRoeAverage)

           call computeIncomingJacobianOfInviscidFlux(nDimensions, nSpecies, localRoeAverage,&
                localMetricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,        &
                incomingDirection, incomingJacobianOfInviscidFlux)

           select case (mode)

           case (FORWARD)

              state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -          &
                   this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *                &
                   matmul(incomingJacobianOfInviscidFlux,                                    &
                   localConservedVariables - localInterfaceConservedVariables)

              if (simulationFlags%viscosityOn) then
                 state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -       &
                      this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *             &
                      matmul(transpose(incomingJacobianOfInviscidFlux),                      &
                      state%adjointVariables(gridIndex,:) -                                  &
                      interfaceAdjointVariables(patchIndex,:))
              else
                 state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +       &
                      this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *             &
                      matmul(transpose(incomingJacobianOfInviscidFlux),                      &
                      state%adjointVariables(gridIndex,:))
              end if

           case (ADJOINT)

              if (simulationFlags%useContinuousAdjoint) then
                 state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -       &
                      this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *             &
                      matmul(transpose(incomingJacobianOfInviscidFlux),                      &
                      state%adjointVariables(gridIndex,:))
              end if

              ! TODO: viscous penalties.

           end select !... mode

        end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  SAFE_DEALLOCATE(incomingJacobianOfInviscidFlux)
  SAFE_DEALLOCATE(localMetricsAlongNormalDirection)
  SAFE_DEALLOCATE(localRoeAverage)
  SAFE_DEALLOCATE(localInterfaceConservedVariables)
  SAFE_DEALLOCATE(localConservedVariables)
  SAFE_DEALLOCATE(interfaceViscousFluxes)
  SAFE_DEALLOCATE(interfaceConservedVariables)

  call endTiming("addBlockInterfacePenalty")

end subroutine addBlockInterfacePenalty

function verifyBlockInterfacePatchUsage(this, patchDescriptor, gridSize, normalDirection,    &
     extent, simulationFlags, success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_BlockInterfacePatch) :: this
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  integer, intent(in) :: gridSize(:), normalDirection, extent(6)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  logical, intent(out) :: success
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Result >>>
  logical :: isPatchUsed

  ! <<< Local variables >>>
  integer :: i

  isPatchUsed = .false.

  success = .false.
  if (normalDirection > size(gridSize) .or. normalDirection == 0) then
     write(message, '(A)') "Normal direction is invalid!"
     return
  end if

  do i = 1, size(gridSize)
     if (extent((i-1)*2+1) < 0 .or. extent((i-1)*2+2) > gridSize(i) .or.                     &
          extent((i-1)*2+1) > extent((i-1)*2+2)) then
        write(message, '(A)') "Invalid extent!"
        return
     end if
  end do

  i = abs(normalDirection)
  if (extent((i-1)*2+1) /= extent((i-1)*2+2)) then
     write(message, '(A)') "Extends more than 1 grid point along normal direction!"
     return
  end if

  success = .true.

  isPatchUsed = .true.

end function verifyBlockInterfacePatchUsage

subroutine collectInterfaceData(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : ADJOINT
  use SolverOptions_enum

  ! <<< Internal modules >>>
  use CNSHelper, only : computeCartesianViscousFluxes

  implicit none

  ! <<< Arguments >>>
  class(t_BlockInterfacePatch) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state

  ! <<< Local variables >>>
  integer :: i, direction, nDimensions, nUnknowns, nSpecies, procRank, ierror
  SCALAR_TYPE, dimension(:,:), allocatable :: velocity, massFraction, stressTensor,          &
       heatFlux, enthalpyFlux, metricsAlongNormalDirection, dataToBeSent
  SCALAR_TYPE, dimension(:,:,:), allocatable :: speciesFlux, viscousFluxes

  if (this%comm == MPI_COMM_NULL) return

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nSpecies = solverOptions%nSpecies
  assert(nSpecies >= 0)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2 + nSpecies)

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  assert(this%nPatchPoints > 0)

  if (simulationFlags%viscosityOn) then
     if (mode == ADJOINT) then
        allocate(dataToBeSent(this%nPatchPoints, 3 * nUnknowns - 1))
     else
        allocate(dataToBeSent(this%nPatchPoints, 2 * nUnknowns - 1))
     end if
  else
     if (mode == ADJOINT) then
        allocate(dataToBeSent(this%nPatchPoints, 2 * nUnknowns))
     else
        allocate(dataToBeSent(this%nPatchPoints, nUnknowns))
     end if
  end if

  call MPI_Comm_rank(this%comm, procRank, ierror)

  if (procRank == 0 .and. size(dataToBeSent, 2) /= size(this%sendBuffer, 2)) then
     SAFE_DEALLOCATE(this%sendBuffer)
     SAFE_DEALLOCATE(this%receiveBuffer)     
     allocate(this%sendBuffer(product(this%globalSize), size(dataToBeSent, 2)))
     allocate(this%receiveBuffer(product(this%globalSize), size(dataToBeSent, 2)))
  end if

  call MPI_Barrier(this%comm, ierror)

  call this%collect(state%conservedVariables, dataToBeSent(:,1:nUnknowns))

  if (mode == ADJOINT)                                                                       &
       call this%collect(state%adjointVariables, dataToBeSent(:,nUnknowns+1:2*nUnknowns))

  if (simulationFlags%viscosityOn) then

     allocate(velocity(this%nPatchPoints, nDimensions))
     allocate(massFraction(this%nPatchPoints, nSpecies))
     allocate(stressTensor(this%nPatchPoints, nDimensions ** 2))
     allocate(heatFlux(this%nPatchPoints, nDimensions))
     allocate(viscousFluxes(this%nPatchPoints, nUnknowns, nDimensions))
     allocate(metricsAlongNormalDirection(this%nPatchPoints, nDimensions))
     if (nSpecies > 0) then
        allocate(speciesFlux(this%nPatchPoints, nSpecies, nDimensions))
        if (solverOptions%equationOfState == IDEAL_GAS_MIXTURE)                              &
             allocate(enthalpyFlux(this%nPatchPoints, nDimensions))
     end if

     call this%collect(state%velocity, velocity)
     call this%collect(state%massFraction, massFraction)
     call this%collect(state%stressTensor, stressTensor)
     call this%collect(state%heatFlux, heatFlux)
     call this%collect(state%speciesFlux, speciesFlux)
     if (solverOptions%equationOfState == IDEAL_GAS_MIXTURE)                                 &
          call this%collect(state%enthalpyFlux, enthalpyFlux)
     call this%collect(grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction),    &
          metricsAlongNormalDirection)

     if (solverOptions%equationOfState == IDEAL_GAS) then
        call computeCartesianViscousFluxes(nDimensions, nSpecies, velocity, stressTensor,    &
             heatFlux, viscousFluxes, massFraction = massFraction, speciesFlux = speciesFlux)
     else
        call computeCartesianViscousFluxes(nDimensions, nSpecies, velocity, stressTensor,    &
             heatFlux, viscousFluxes, massFraction, speciesFlux, enthalpyFlux)
     end if

     do i = 1, this%nPatchPoints
        this%viscousFluxes(i,:) = matmul(viscousFluxes(i,2:nUnknowns,:),                     &
             metricsAlongNormalDirection(i,:))
     end do

     if (mode == ADJOINT) then
        dataToBeSent(:,2*nUnknowns+1:) = this%viscousFluxes
     else
        dataToBeSent(:,nUnknowns+1:) = this%viscousFluxes
     end if

     SAFE_DEALLOCATE(viscousFluxes)
     SAFE_DEALLOCATE(metricsAlongNormalDirection)
     SAFE_DEALLOCATE(heatFlux)
     SAFE_DEALLOCATE(speciesFlux)
     SAFE_DEALLOCATE(enthalpyFlux)
     SAFE_DEALLOCATE(stressTensor)
     SAFE_DEALLOCATE(velocity)
     SAFE_DEALLOCATE(massFraction)

  end if

  call this%gatherData(dataToBeSent, this%sendBuffer)
  SAFE_DEALLOCATE(dataToBeSent)

end subroutine collectInterfaceData

subroutine reshapeReceivedData(this, indexReordering)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_BlockInterfacePatch) :: this
  integer, intent(in) :: indexReordering(3)

  ! <<< Local variables >>>
  integer :: i, j, k, l, order(3), globalSize(3), nComponents

  if (this%comm == MPI_COMM_NULL .or. .not. allocated(this%receiveBuffer)) return

  assert(allocated(this%sendBuffer))
  assert(all(this%globalSize > 0))
  assert(size(this%sendBuffer, 1) == product(this%globalSize))
  assert(size(this%sendBuffer, 2) > 0)
  assert(all(shape(this%receiveBuffer) == shape(this%sendBuffer)))

  assert(indexReordering(3) == 3)

  if (indexReordering(1) == 1 .and. indexReordering(2) == 2) return

  call startTiming("reshapeReceivedData")

  nComponents = size(this%receiveBuffer, 2)
  globalSize = this%globalSize

  order = indexReordering

  if (abs(order(1)) == 2 .and. abs(order(2)) == 1) then

     do l = 1, nComponents
        do k = 1, globalSize(3)
           do j = 1, globalSize(2)
              do i = 1, globalSize(1)
                 this%sendBuffer(i + globalSize(1) * (j - 1 +                                &
                      globalSize(2) * (k - 1)), l) =                                         &
                      this%receiveBuffer(j + globalSize(2) * (i - 1 +                        &
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
                 this%sendBuffer(i + globalSize(1) * (j - 1 +                                &
                      globalSize(2) * (k - 1)), l) =                                         &
                      this%receiveBuffer(globalSize(1) + 1 - i + globalSize(1) * (j - 1 +    &
                      globalSize(2) * (k - 1)), l)
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
                 this%sendBuffer(i + globalSize(1) * (j - 1 +                                &
                      globalSize(2) * (k - 1)), l) =                                         &
                      this%receiveBuffer(i + globalSize(1) * (globalSize(2) - j +            &
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
                 this%sendBuffer(i + globalSize(1) * (j - 1 +                                &
                      globalSize(2) * (k - 1)), l) =                                         &
                      this%receiveBuffer(globalSize(1) + 1 - i +                             &
                      globalSize(1) * (globalSize(2) - j + globalSize(2) * (k - 1)), l)
              end do
           end do
        end do
     end do

     this%receiveBuffer = this%sendBuffer

  end if

  call endTiming("reshapeReceivedData")

end subroutine reshapeReceivedData
