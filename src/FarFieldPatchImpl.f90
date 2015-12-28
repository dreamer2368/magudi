#include "config.h"

subroutine setupFarFieldPatch(this, index, comm, patchDescriptor,                            &
     grid, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use FarFieldPatch_mod, only : t_FarFieldPatch
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_FarFieldPatch) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: nDimensions, nUnknowns, nSpecies, direction, i, j, k, gridIndex, patchIndex
  real(SCALAR_KIND) :: radius, holeRadius, holePosition(3)
  logical :: holeInsideShape
  character(len = STRING_LENGTH) :: key, message, holeShape

  call this%cleanup()
  call this%setupBase(index, comm, patchDescriptor, grid, simulationFlags, solverOptions)

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nSpecies = solverOptions%nSpecies
  assert(nSpecies >= 0)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2 + nSpecies)

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  if (this%nPatchPoints > 0) then
     if (simulationFlags%viscosityOn) then
        allocate(this%viscousFluxes(this%nPatchPoints, nUnknowns, nDimensions))
        allocate(this%targetViscousFluxes(this%nPatchPoints, nUnknowns, nDimensions))
     end if
  end if


  write(key, '(A)') "patches/" // trim(patchDescriptor%name) // "/"

  ! Inviscid penalty amount.
  this%inviscidPenaltyAmount = getOption("defaults/inviscid_penalty_amount", 1.0_wp)
  this%inviscidPenaltyAmount = getOption(trim(key) // "inviscid_penalty_amount",             &
       this%inviscidPenaltyAmount)
  this%inviscidPenaltyAmount = sign(this%inviscidPenaltyAmount,                              &
       real(this%normalDirection, wp))
  this%inviscidPenaltyAmount = this%inviscidPenaltyAmount /                                  &
       grid%firstDerivative(direction)%normBoundary(1)

  ! Viscous penalty amount.
  if (simulationFlags%viscosityOn) then
     this%viscousPenaltyAmount = getOption("defaults/viscous_penalty_amount", 1.0_wp)
     this%viscousPenaltyAmount = getOption(trim(key) // "viscous_penalty_amount",            &
          this%viscousPenaltyAmount)
     this%viscousPenaltyAmount = sign(this%viscousPenaltyAmount,                             &
          real(this%normalDirection, wp))
     this%viscousPenaltyAmount = this%viscousPenaltyAmount /                                 &
          grid%firstDerivative(direction)%normBoundary(1)
  else
     this%viscousPenaltyAmount = 0.0_wp
  end if

  ! Add a hole to the patch.
  write(key, '(A,I0.0)') "patches/" // trim(patchDescriptor%name) // "/"
  if (getOption(trim(key) // "include_hole", .false.)) then

     holeInsideShape = getOption(trim(key) // "hole_inside_shape", .true.)
     call getRequiredOption(trim(key) // "hole_shape", holeShape, comm)

     select case(trim(holeShape))

     case("CIRCLE")

        call getRequiredOption(trim(key) // "hole_radius", holeRadius, comm)
        holePosition = 0.0_wp
        do i = 1, grid%nDimensions
           write(message, '(A,I0.0)') trim(key) // "hole_position_", i
           call getRequiredOption(trim(message), holePosition(i), comm)
        end do

        allocate(this%hole(this%nPatchPoints))
        this%hole = 0
        do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
           do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
              do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                 gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                &
                      (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                  &
                      (k - 1 - this%gridOffset(3)))
                 if (grid%iblank(gridIndex) == 0) cycle
                 patchIndex = i - this%offset(1) + this%localSize(1) *                       &
                      (j - 1 - this%offset(2) + this%localSize(2) *                          &
                      (k - 1 - this%offset(3)))

                 radius = sqrt(sum((grid%coordinates(gridIndex, 1:grid%nDimensions) -        &
                      holePosition(1:grid%nDimensions)) ** 2))
                 if (holeInsideShape .and. radius <= holeRadius) then
                    this%hole(patchIndex) = 1
                 else if (.not.holeInsideShape .and. radius > holeRadius) then
                    this%hole(patchIndex) = 1
                 end if

              end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
           end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
        end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

        case default

        write(message, '(A)') "Unknown hole shape on patch " //                              &
             trim(patchDescriptor%name) // " `" // trim(holeShape) // "'"
        call gracefulExit(comm, message)

     end select

  end if

end subroutine setupFarFieldPatch

subroutine cleanupFarFieldPatch(this)

  ! <<< Derived types >>>
  use FarFieldPatch_mod, only : t_FarFieldPatch

  implicit none

  ! <<< Arguments >>>
  class(t_FarFieldPatch) :: this

  call this%cleanupBase()

  SAFE_DEALLOCATE(this%viscousFluxes)
  SAFE_DEALLOCATE(this%targetViscousFluxes)

end subroutine cleanupFarFieldPatch

subroutine addFarFieldPenalty(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use FarFieldPatch_mod, only : t_FarFieldPatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use CNSHelper
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_FarFieldPatch) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, nUnknowns, nSpecies,                                      &
       direction, incomingDirection, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: localTargetState(:), metricsAlongNormalDirection(:),           &
       incomingJacobianOfInviscidFlux(:,:)

  assert_key(mode, (FORWARD, ADJOINT))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))
  assert(allocated(state%targetState))

  call startTiming("addFarFieldPenalty")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nSpecies = solverOptions%nSpecies
  assert(nSpecies >= 0)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2 + nSpecies)

  direction = abs(this%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  if (mode == ADJOINT .and. simulationFlags%useContinuousAdjoint) then
     incomingDirection = -this%normalDirection
  else
     incomingDirection = +this%normalDirection
  end if

  allocate(localTargetState(nUnknowns))
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
           if (allocated(this%hole)) then
              if (this%hole(patchIndex) == 1) cycle
           end if

           localTargetState = state%targetState(gridIndex,:)
           metricsAlongNormalDirection =                                                     &
                grid%metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

           call computeIncomingJacobianOfInviscidFlux(nDimensions, nSpecies,                 &
                localTargetState, metricsAlongNormalDirection,                               &
                solverOptions%ratioOfSpecificHeats, incomingDirection,                       &
                incomingJacobianOfInviscidFlux)

           select case (mode)

           case (FORWARD)

              state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -          &
                   this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *                &
                   matmul(incomingJacobianOfInviscidFlux,                                    &
                   state%conservedVariables(gridIndex,:) - localTargetState)

              if (simulationFlags%viscosityOn)                                               &
                   state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +     &
                   this%viscousPenaltyAmount * grid%jacobian(gridIndex, 1) *                 &
                   matmul(this%viscousFluxes(patchIndex,:,:) -                               &
                   this%targetViscousFluxes(patchIndex,:,:), metricsAlongNormalDirection)

           case (ADJOINT)

              if (simulationFlags%useContinuousAdjoint) then
                 state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -       &
                      this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *             &
                      matmul(transpose(incomingJacobianOfInviscidFlux),                      &
                      state%adjointVariables(gridIndex,:))
              else
                 state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +       &
                      this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *             &
                      matmul(transpose(incomingJacobianOfInviscidFlux),                      &
                      state%adjointVariables(gridIndex,:))
              end if

           end select !... mode

        end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  SAFE_DEALLOCATE(incomingJacobianOfInviscidFlux)
  SAFE_DEALLOCATE(metricsAlongNormalDirection)
  SAFE_DEALLOCATE(localTargetState)

  call endTiming("addFarFieldPenalty")

end subroutine addFarFieldPenalty

function verifyFarFieldPatchUsage(this, patchDescriptor, gridSize, normalDirection,          &
     extent, simulationFlags, success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use FarFieldPatch_mod, only : t_FarFieldPatch
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_FarFieldPatch) :: this
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

  if (.not. simulationFlags%useTargetState) then
     write(message, '(A)')                                                                   &
          "No target state available for enforcing far-field boundary conditions!"
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

end function verifyFarFieldPatchUsage
