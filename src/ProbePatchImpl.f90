#include "config.h"

subroutine setupProbePatch(this, index, comm, patchDescriptor,                               &
     grid, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use ProbePatch_mod, only : t_ProbePatch
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_ProbePatch) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: outputPrefix
  integer :: probeBufferSize

  call this%cleanup()
  call this%setupBase(index, comm, patchDescriptor, grid, simulationFlags, solverOptions)

  outputPrefix = getOption("output_prefix", PROJECT_NAME)
  write(this%probeFilename, '(4A)') trim(outputPrefix), ".probe_",                           &
       trim(patchDescriptor%name), ".dat"

  probeBufferSize = getOption("probe_buffer_size", 1)
  if (this%nPatchPoints > 0) then
     allocate(this%probeBuffer(this%nPatchPoints, solverOptions%nUnknowns, probeBufferSize))
     this%probeBuffer = 0.0_wp
  end if

end subroutine setupProbePatch

subroutine cleanupProbePatch(this)

  ! <<< External modules >>>
  use MPI, only : MPI_OFFSET_KIND

  ! <<< Derived types >>>
  use ProbePatch_mod, only : t_ProbePatch

  implicit none

  ! <<< Arguments >>>
  class(t_ProbePatch) :: this

  call this%cleanupBase()

  SAFE_DEALLOCATE(this%probeBuffer)

  this%iProbeBuffer = 0
  this%probeFileOffset = int(0, MPI_OFFSET_KIND)

end subroutine cleanupProbePatch

subroutine updateProbePatch(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use ProbePatch_mod, only : t_ProbePatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_ProbePatch) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

end subroutine updateProbePatch

function verifyProbePatchUsage(this, patchDescriptor, gridSize, normalDirection,             &
     extent, simulationFlags, success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use ProbePatch_mod, only : t_ProbePatch
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_ProbePatch) :: this
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

  do i = 1, size(gridSize)
     if (extent((i-1)*2+1) < 0 .or. extent((i-1)*2+2) > gridSize(i) .or.                     &
          extent((i-1)*2+1) > extent((i-1)*2+2)) then
        write(message, '(A)') "Invalid extent!"
        return
     end if
  end do

  success = .true.

  isPatchUsed = .true.

end function verifyProbePatchUsage

subroutine saveSolutionOnProbe(this)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use ProbePatch_mod, only : t_ProbePatch

  ! <<< Internal modules >>>
  use InputHelper, only : getFreeUnit

  ! <<< Arguments >>>
  class(t_ProbePatch) :: this

  ! <<< Local variables >>>
  integer :: arrayOfSizes(5), arrayOfSubsizes(5), arrayOfStarts(5), mpiScalarSubarrayType,   &
       mpiFileHandle, dataSize, ierror

  if (this%comm == MPI_COMM_NULL .or. this%iProbeBuffer == 0) return

  arrayOfSizes(1:3) = this%globalSize
  arrayOfSizes(4) = size(this%probeBuffer, 2)
  arrayOfSizes(5) = this%iProbeBuffer
  arrayOfSubsizes(1:3) = this%localSize
  arrayOfSubsizes(4) = size(this%probeBuffer, 2)
  arrayOfSubsizes(5) = this%iProbeBuffer
  arrayOfStarts(1:3) = this%offset - this%extent(1::2) + 1
  arrayOfStarts(4:5) = 0
  call MPI_Type_create_subarray(5, arrayOfSizes, arrayOfSubsizes, arrayOfStarts,             &
       MPI_ORDER_FORTRAN, SCALAR_TYPE_MPI, mpiScalarSubarrayType, ierror)
  call MPI_Type_commit(mpiScalarSubarrayType, ierror)

  call MPI_File_open(this%comm, trim(this%probeFilename) // char(0), MPI_MODE_WRONLY,        &
       MPI_INFO_NULL, mpiFileHandle, ierror)

  assert(this%probeFileOffset >= 0)

  call MPI_File_set_view(mpiFileHandle, this%probeFileOffset, SCALAR_TYPE_MPI,               &
       mpiScalarSubarrayType, "native", MPI_INFO_NULL, ierror)

  dataSize = this%nPatchPoints * size(this%probeBuffer, 2) * this%iProbeBuffer
  call MPI_File_write_all(mpiFileHandle, this%probeBuffer, dataSize,                         &
       SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)

  this%probeFileOffset = this%probeFileOffset +                                              &
       SIZEOF_SCALAR * product(int(this%globalSize, MPI_OFFSET_KIND)) *                      &
       size(this%probeBuffer, 2) * this%iProbeBuffer

  call MPI_File_close(mpiFileHandle, ierror)

  call MPI_Type_free(mpiScalarSubarrayType, ierror)

end subroutine saveSolutionOnProbe
