#include "config.h"

subroutine setupJetExcitationPatch(this, index, comm, patchDescriptor,                       &
     grid, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags
  use JetExcitationPatch_mod, only : t_JetExcitationPatch

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit
  use PLOT3DHelper

  implicit none

  ! <<< Arguments >>>
  class(t_JetExcitationPatch) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer :: i, nDimensions, procRank, ierror
  character(len = STRING_LENGTH) :: outputPrefix, key, filename, message
  logical :: success
  integer, allocatable :: globalPatchSizes(:,:)
  integer(kind = MPI_OFFSET_KIND) :: offset

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (2, 3))

  call this%cleanup()
  call this%t_SpongePatch%setup(index, comm, patchDescriptor,                                &
       grid, simulationFlags, solverOptions)

  outputPrefix = getOption("output_prefix", PROJECT_NAME)

  write(key, '(A)') "patches/" // trim(patchDescriptor%name) // "/"

  call getRequiredOption(trim(key) // "number_of_modes", this%nModes, this%comm)
  this%nModes = min(max(0, this%nModes), 99)

  if (this%comm /= MPI_COMM_NULL) call MPI_Comm_rank(this%comm, procRank, ierror)

  if (this%nModes > 0 .and. this%nPatchPoints > 0) then

     allocate(this%angularFrequencies(this%nModes))
     allocate(this%perturbationReal(this%nPatchPoints, solverOptions%nUnknowns, this%nModes))
     allocate(this%perturbationImag(this%nPatchPoints, solverOptions%nUnknowns, this%nModes))
     
     do i = 1, this%nModes

        write(key, '(A,I2,A)') "patches/" // trim(patchDescriptor%name) // "/mode", i, "/"
        call getRequiredOption(trim(key) // "angular_frequency",                             &
             this%angularFrequencies(i), this%comm)
        
        write(filename, '(2A,I2,A)') outputPrefix, "-", i, ".eigenmodes_real.q"
        call plot3dDetectFormat(this%comm, filename, success,                                &
             globalGridSizes = globalPatchSizes)
        if (.not. success) call gracefulExit(this%comm, plot3dErrorMessage)

        if (size(globalPatchSizes, 1) /= grid%nDimensions .or.                               &
             size(globalPatchSizes, 2) < grid%index) then
           write(message, '(2A)') trim(filename), ": mismatch in grid dimensions."
           call gracefulExit(this%comm, message)
        end if

        offset = plot3dGetOffset(this%comm, filename, grid%index, success)
        if (.not. success) call gracefulExit(this%comm, plot3dErrorMessage)
        call plot3dReadSingleSolution(this%comm, trim(filename), offset,                     &
             this%mpiAllScalarSubarrayTypes(procRank + 1), this%globalSize,                  &
             this%perturbationReal(:,:,i), success)
        if (.not. success) call gracefulExit(this%comm, plot3dErrorMessage)

        write(filename, '(2A,I2,A)') outputPrefix, "-", i, ".eigenmodes_imag.q"
        call plot3dDetectFormat(this%comm, filename, success,                                &
             globalGridSizes = globalPatchSizes)
        if (.not. success) call gracefulExit(this%comm, plot3dErrorMessage)

        if (size(globalPatchSizes, 1) /= grid%nDimensions .or.                               &
             size(globalPatchSizes, 2) < grid%index) then
           write(message, '(2A)') trim(filename), ": mismatch in grid dimensions."
           call gracefulExit(this%comm, message)
        end if

        offset = plot3dGetOffset(this%comm, filename, grid%index, success)
        if (.not. success) call gracefulExit(this%comm, plot3dErrorMessage)
        call plot3dReadSingleSolution(this%comm, trim(filename), offset,                     &
             this%mpiAllScalarSubarrayTypes(procRank + 1), this%globalSize,                  &
             this%perturbationImag(:,:,i), success)
        if (.not. success) call gracefulExit(this%comm, plot3dErrorMessage)

     end do !... i = 1, this%nModes

  end if !... this%nModes > 0 .and. this%nPatchPoints > 0

end subroutine setupJetExcitationPatch

subroutine cleanupJetExcitationPatch(this)

  ! <<< Derived types >>>
  use JetExcitationPatch_mod, only : t_JetExcitationPatch

  implicit none

  ! <<< Arguments >>>
  class(t_JetExcitationPatch) :: this

  call this%t_SpongePatch%cleanup()

  SAFE_DEALLOCATE(this%angularFrequencies)
  SAFE_DEALLOCATE(this%perturbationReal)
  SAFE_DEALLOCATE(this%perturbationImag)

end subroutine cleanupJetExcitationPatch

subroutine addJetExcitation(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use JetExcitationPatch_mod, only : t_JetExcitationPatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_JetExcitationPatch) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, nDimensions, gridIndex, patchIndex

  assert_key(mode, (FORWARD, ADJOINT))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  if (mode == ADJOINT) return

  call startTiming("addJetExcitation")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (2, 3))

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
           do l = 1, this%nModes
              state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -          &
                   this%spongeStrength(patchIndex) *                                         &
                   (this%perturbationReal(patchIndex,:,l) *                                  &
                   cos(this%angularFrequencies(l) * state%time) -                            &
                   this%perturbationImag(patchIndex,:,l) *                                   &
                   sin(this%angularFrequencies(l) * state%time))
           end do
        end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  call endTiming("addJetExcitation")

end subroutine addJetExcitation

function verifyJetExcitationPatchUsage(this, patchDescriptor, gridSize,                      &
     normalDirection, extent, simulationFlags,                                               &
     success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags
  use JetExcitationPatch_mod, only : t_JetExcitationPatch

  implicit none

  ! <<< Arguments >>>
  class(t_JetExcitationPatch) :: this
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  integer, intent(in) :: gridSize(:), normalDirection, extent(6)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  logical, intent(out) :: success
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Result >>>
  logical :: isPatchUsed

  isPatchUsed = this%t_SpongePatch%verifyUsage(patchDescriptor, gridSize, normalDirection,   &
       extent, simulationFlags, success, message)

end function verifyJetExcitationPatchUsage
