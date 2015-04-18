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

     allocate(this%interfaceConservedVariables(this%nPatchPoints, solverOptions%nUnknowns))
     nExchangedVariables = solverOptions%nUnknowns

     if (simulationFlags%viscosityOn) then
        allocate(this%viscousFluxes(this%nPatchPoints, solverOptions%nUnknowns - 1))
        allocate(this%interfaceViscousFluxes(this%nPatchPoints, solverOptions%nUnknowns - 1))
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
          "viscous_penalty_amount", 1.0_wp)
     this%viscousPenaltyAmount = this%viscousPenaltyAmount *                                 &
          solverOptions%reynoldsNumberInverse
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
  SAFE_DEALLOCATE(this%interfaceViscousFluxes)
  SAFE_DEALLOCATE(this%interfaceConservedVariables)
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
  integer :: i, j, k, nDimensions, nUnknowns, direction, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: receivedData(:,:), localConservedVariables(:),                 &
       localInterfaceConservedVariables(:), roeAverage(:), metricsAlongNormalDirection(:),   &
       incomingJacobianOfInviscidFlux(:,:)

  assert_key(mode, (FORWARD, ADJOINT))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  call startTiming("addBlockInterfacePenalty")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2)

  if (this%comm /= MPI_COMM_NULL) then
     if (simulationFlags%viscosityOn) then
        allocate(receivedData(this%nPatchPoints, 2 * nUnknowns - 1))
     else
        allocate(receivedData(this%nPatchPoints, nUnknowns))
     end if
     call this%scatterData(this%receiveBuffer, receivedData)
     this%interfaceConservedVariables = receivedData(:,1:nUnknowns)
     if (simulationFlags%viscosityOn)                                                        &
          this%interfaceViscousFluxes = receivedData(:,nUnknowns+1:)
     SAFE_DEALLOCATE(receivedData)
  end if

  allocate(localConservedVariables(nUnknowns))
  allocate(localInterfaceConservedVariables(nUnknowns))
  allocate(roeAverage(nUnknowns))
  allocate(metricsAlongNormalDirection(nDimensions))
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
           localInterfaceConservedVariables = this%interfaceConservedVariables(patchIndex,:)
           metricsAlongNormalDirection =                                                     &
                grid%metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

           call computeRoeAverage(nDimensions, localConservedVariables,                      &
                localInterfaceConservedVariables, solverOptions%ratioOfSpecificHeats,        &
                roeAverage)

           select case (nDimensions)
           case (1)
              call computeIncomingJacobianOfInviscidFlux1D(roeAverage,                       &
                   metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,          &
                   this%normalDirection, incomingJacobianOfInviscidFlux)
           case (2)
              call computeIncomingJacobianOfInviscidFlux2D(roeAverage,                       &
                   metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,          &
                   this%normalDirection, incomingJacobianOfInviscidFlux)
           case (3)
              call computeIncomingJacobianOfInviscidFlux3D(roeAverage,                       &
                   metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,          &
                   this%normalDirection, incomingJacobianOfInviscidFlux)
           end select !... nDimensions

           select case (mode)
           case (FORWARD)

              state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -          &
                   this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *                &
                   matmul(incomingJacobianOfInviscidFlux,                                    &
                   localConservedVariables - this%interfaceConservedVariables(patchIndex,:))

              if (simulationFlags%viscosityOn) then
                 state%rightHandSide(gridIndex,2:nUnknowns) =                                &
                      state%rightHandSide(gridIndex,2:nUnknowns) +                           &
                      this%viscousPenaltyAmount * grid%jacobian(gridIndex, 1) *              &
                      (this%viscousFluxes(patchIndex,:) -                                    &
                      this%interfaceViscousFluxes(patchIndex,:))
              end if

           case (ADJOINT)

              ! TODO.

           end select !... mode

        end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  SAFE_DEALLOCATE(incomingJacobianOfInviscidFlux)
  SAFE_DEALLOCATE(metricsAlongNormalDirection)
  SAFE_DEALLOCATE(roeAverage)
  SAFE_DEALLOCATE(localInterfaceConservedVariables)
  SAFE_DEALLOCATE(localConservedVariables)

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
  integer :: i, direction, nDimensions, nUnknowns
  SCALAR_TYPE, dimension(:,:), allocatable :: velocity,                                      &
       stressTensor, heatFlux, dataToBeSent
  SCALAR_TYPE, allocatable :: viscousFluxes(:,:,:)

  if (this%comm == MPI_COMM_NULL) return

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  nUnknowns = solverOptions%nUnknowns

  assert(this%nPatchPoints > 0)

  if (simulationFlags%viscosityOn) then

     allocate(velocity(this%nPatchPoints, nDimensions))
     allocate(stressTensor(this%nPatchPoints, nDimensions ** 2))
     allocate(heatFlux(this%nPatchPoints, nDimensions))
     allocate(viscousFluxes(this%nPatchPoints, nUnknowns, nDimensions))

     call this%collect(state%velocity, velocity)
     call this%collect(state%stressTensor, stressTensor)
     call this%collect(state%heatFlux, heatFlux)

     call computeCartesianViscousFluxes(nDimensions, velocity,                               &
          stressTensor, heatFlux, viscousFluxes)

     do i = 1, this%nPatchPoints
        this%viscousFluxes(i,:) = matmul(viscousFluxes(i,2:nUnknowns,:),                     &
             grid%metrics(i,1+nDimensions*(direction-1):nDimensions*direction))
     end do

     SAFE_DEALLOCATE(viscousFluxes)
     SAFE_DEALLOCATE(heatFlux)
     SAFE_DEALLOCATE(stressTensor)
     SAFE_DEALLOCATE(velocity)

  end if

  if (simulationFlags%viscosityOn) then
     allocate(dataToBeSent(this%nPatchPoints, 2 * nUnknowns - 1))
  else
     allocate(dataToBeSent(this%nPatchPoints, nUnknowns))
  end if

  call this%collect(state%conservedVariables, dataToBeSent(:,1:nUnknowns))
  if (simulationFlags%viscosityOn) dataToBeSent(:,nUnknowns+1:) = this%viscousFluxes
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
