#include "config.h"

module ProbePatch_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  use MPI, only : MPI_OFFSET_KIND
  use Patch_mod, only : t_Patch

  implicit none
  private

  type, extends(t_Patch), public :: t_ProbePatch

     integer :: iProbeBuffer = 0
     integer(kind = MPI_OFFSET_KIND) :: probeFileOffset = int(0, MPI_OFFSET_KIND)
     character(len = STRING_LENGTH) :: probeFilename
     real(SCALAR_KIND), allocatable :: probeBuffer(:,:,:)
     integer :: saveMode, loadMode

   contains

     procedure, pass :: setup
     procedure, pass :: cleanup
     procedure, pass :: updateRhs
     procedure, pass :: update
     procedure, pass :: loadData
     procedure, pass :: saveData
     procedure, pass :: reset
     procedure, pass :: seekToEOF

  end type t_ProbePatch

contains

  subroutine setup(this, name, comm, grid, state, extent,                                    &
       normalDirection, simulationFlags, solverOptions)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod, only : t_State
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : FORWARD, UNDEFINED

    ! <<< Internal modules >>>
    use InputHelper, only : getOption
    use ErrorHandler, only : issueWarning

    implicit none

    ! <<< Arguments >>>
    class(t_ProbePatch) :: this
    character(len = *), intent(in) :: name
    integer, intent(in) :: comm
    class(t_Grid), intent(in) :: grid
    class(t_State), intent(in) :: state
    integer, intent(in) :: extent(6), normalDirection
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions

    ! <<< Local variables >>>
    character(len = STRING_LENGTH) :: outputPrefix
    integer :: nDimensions, nUnknowns

    call this%cleanup()
    call this%setupBase(name, comm, grid, extent, normalDirection)

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns >= nDimensions + 2)

    outputPrefix = getOption("output_prefix", PROJECT_NAME)
    write(this%probeFilename, '(4A)') trim(outputPrefix), ".probe_",                         &
         trim(name), ".dat"

    this%saveMode = FORWARD
    this%loadMode = UNDEFINED

  end subroutine setup

  subroutine cleanup(this)

    implicit none

    ! <<< Arguments >>>
    class(t_ProbePatch) :: this

    call this%cleanupBase()

    SAFE_DEALLOCATE(this%probeBuffer)

    this%iProbeBuffer = 0
    this%probeFileOffset = int(0, MPI_OFFSET_KIND)

  end subroutine cleanup

  subroutine updateRhs(this, mode, simulationFlags, solverOptions, grid, state)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod, only : t_State
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

  end subroutine updateRhs

  subroutine update(this, mode, simulationFlags, solverOptions, grid, state)

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
    class(t_ProbePatch) :: this
    integer, intent(in) :: mode
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid), intent(in) :: grid
    class(t_State) :: state

    ! <<< Local variables >>>
    integer :: nDimensions, nUnknowns

    if (this%comm == MPI_COMM_NULL) return

    assert_key(mode, (FORWARD, ADJOINT))
    if (mode /= this%saveMode) return

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns >= nDimensions + 2)

    assert(allocated(this%probeBuffer))
    assert(this%iProbeBuffer >= 1 .and. this%iProbeBuffer <= size(this%probeBuffer, 3))
    assert(size(this%probeBuffer, 2) == nUnknowns)

    select case (mode)
    case (FORWARD)
       call this%collect(state%conservedVariables, this%probeBuffer(:,:,this%iProbeBuffer))
    case (ADJOINT)
       call this%collect(state%adjointVariables, this%probeBuffer(:,:,this%iProbeBuffer))
    end select

  end subroutine update

  subroutine loadData(this)

    ! <<< External modules >>>
    use MPI

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : UNDEFINED

    implicit none

    ! <<< Arguments >>>
    class(t_ProbePatch) :: this

    ! <<< Local variables >>>
    integer :: arrayOfSizes(5), arrayOfSubsizes(5), arrayOfStarts(5),                        &
         mpiScalarSubarrayType, mpiFileHandle, dataSize, ierror
    integer(kind = MPI_OFFSET_KIND) :: nBytesToRead, fileSize

    if (this%comm == MPI_COMM_NULL) return
    if (this%saveMode == UNDEFINED .or. this%loadMode == UNDEFINED) return

    call MPI_File_open(this%comm, trim(this%probeFilename) // char(0), MPI_MODE_RDONLY,      &
         MPI_INFO_NULL, mpiFileHandle, ierror)
    call MPI_File_get_size(mpiFileHandle, fileSize, ierror)

    nBytesToRead = SIZEOF_SCALAR * product(int(this%globalSize, MPI_OFFSET_KIND)) *          &
         size(this%probeBuffer, 2) * this%iProbeBuffer

    if (this%saveMode /= this%loadMode) then
       if (this%probeFileOffset - nBytesToRead <= int(0, MPI_OFFSET_KIND)) then
          this%iProbeBuffer = int(this%probeFileOffset / (SIZEOF_SCALAR *                    &
               product(int(this%globalSize, MPI_OFFSET_KIND)) * size(this%probeBuffer, 2)))
          this%probeFileOffset = 0
       else
          this%probeFileOffset = this%probeFileOffset - nBytesToRead
       end if
    else
       if (this%probeFileOffset + nBytesToRead >= fileSize) then
          this%iProbeBuffer = int((fileSize - this%probeFileOffset) / (SIZEOF_SCALAR *       &
               product(int(this%globalSize, MPI_OFFSET_KIND)) * size(this%probeBuffer, 2)))
       end if
    end if

    assert(this%probeFileOffset >= 0 .and. this%probeFileOffset < fileSize)

    arrayOfSizes(1:3) = this%globalSize
    arrayOfSizes(4) = size(this%probeBuffer, 2)
    arrayOfSizes(5) = this%iProbeBuffer
    arrayOfSubsizes(1:3) = this%localSize
    arrayOfSubsizes(4) = size(this%probeBuffer, 2)
    arrayOfSubsizes(5) = this%iProbeBuffer
    arrayOfStarts(1) = this%offset(1) - this%iMin + 1
    arrayOfStarts(2) = this%offset(2) - this%jMin + 1
    arrayOfStarts(3) = this%offset(3) - this%kMin + 1
    arrayOfStarts(4:5) = 0
    call MPI_Type_create_subarray(5, arrayOfSizes, arrayOfSubsizes, arrayOfStarts,           &
         MPI_ORDER_FORTRAN, SCALAR_TYPE_MPI, mpiScalarSubarrayType, ierror)
    call MPI_Type_commit(mpiScalarSubarrayType, ierror)

    call MPI_File_set_view(mpiFileHandle, this%probeFileOffset, SCALAR_TYPE_MPI,             &
         mpiScalarSubarrayType, "native", MPI_INFO_NULL, ierror)

    dataSize = this%nPatchPoints * size(this%probeBuffer, 2) * this%iProbeBuffer
    call MPI_File_read_all(mpiFileHandle, this%probeBuffer, dataSize,                        &
         SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)

    call MPI_File_close(mpiFileHandle, ierror)

    call MPI_Type_free(mpiScalarSubarrayType, ierror)

    if (this%saveMode == this%loadMode)                                                      &
         this%probeFileOffset = max(fileSize, this%probeFileOffset + nBytesToRead)

  end subroutine loadData

  subroutine saveData(this)

    ! <<< External modules >>>
    use MPI

    ! <<< Arguments >>>
    class(t_ProbePatch) :: this

    ! <<< Local variables >>>
    integer :: arrayOfSizes(5), arrayOfSubsizes(5), arrayOfStarts(5),                        &
         mpiScalarSubarrayType, mpiFileHandle, dataSize, ierror

    if (this%comm == MPI_COMM_NULL .or. this%iProbeBuffer == 0) return

    arrayOfSizes(1:3) = this%globalSize
    arrayOfSizes(4) = size(this%probeBuffer, 2)
    arrayOfSizes(5) = this%iProbeBuffer
    arrayOfSubsizes(1:3) = this%localSize
    arrayOfSubsizes(4) = size(this%probeBuffer, 2)
    arrayOfSubsizes(5) = this%iProbeBuffer
    arrayOfStarts(1) = this%offset(1) - this%iMin + 1
    arrayOfStarts(2) = this%offset(2) - this%jMin + 1
    arrayOfStarts(3) = this%offset(3) - this%kMin + 1
    arrayOfStarts(4:5) = 0
    call MPI_Type_create_subarray(5, arrayOfSizes, arrayOfSubsizes, arrayOfStarts,           &
         MPI_ORDER_FORTRAN, SCALAR_TYPE_MPI, mpiScalarSubarrayType, ierror)
    call MPI_Type_commit(mpiScalarSubarrayType, ierror)

    call MPI_File_open(this%comm, trim(this%probeFilename) // char(0), MPI_MODE_WRONLY,      &
         MPI_INFO_NULL, mpiFileHandle, ierror)

    assert(this%probeFileOffset >= 0)

    call MPI_File_set_view(mpiFileHandle, this%probeFileOffset, SCALAR_TYPE_MPI,             &
         mpiScalarSubarrayType, "native", MPI_INFO_NULL, ierror)

    dataSize = this%nPatchPoints * size(this%probeBuffer, 2) * this%iProbeBuffer
    call MPI_File_write_all(mpiFileHandle, this%probeBuffer, dataSize,                       &
         SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)

    this%probeFileOffset = this%probeFileOffset +                                            &
         SIZEOF_SCALAR * product(int(this%globalSize, MPI_OFFSET_KIND)) *                    &
         size(this%probeBuffer, 2) * this%iProbeBuffer

    call MPI_File_close(mpiFileHandle, ierror)

    call MPI_Type_free(mpiScalarSubarrayType, ierror)

  end subroutine saveData

  subroutine reset(this)

    ! <<< External modules >>>
    use MPI

    ! <<< Internal modules >>>
    use InputHelper, only : getFreeUnit

    implicit none

    ! <<< Arguments >>>
    class(t_ProbePatch) :: this

    ! <<< Local variables >>>
    integer :: fileUnit, stat, procRank, ierror

    if (this%comm == MPI_COMM_NULL) return

    call MPI_Comm_rank(this%comm, procRank, ierror)

    if (procRank == 0) then
       open(unit = getFreeUnit(fileUnit), file = trim(this%probeFilename),                   &
            iostat = stat, status = 'old')
       if (stat == 0) close(fileUnit, status = 'delete')
       open(unit = getFreeUnit(fileUnit), file = trim(this%probeFilename),                   &
            action = 'write', status = 'unknown')
       close(fileUnit)
    end if
    call MPI_Barrier(this%comm, ierror)
    this%iProbeBuffer = 0
    this%probeFileOffset = int(0, MPI_OFFSET_KIND)

  end subroutine reset

  subroutine seekToEOF(this)

    ! <<< External modules >>>
    use MPI

    ! <<< Internal modules >>>
    use InputHelper, only : getFreeUnit
    use ErrorHandler, only : gracefulExit

    implicit none

    ! <<< Arguments >>>
    class(t_ProbePatch) :: this

    ! <<< Local variables >>>
    integer :: mpiFileHandle, procRank, ierror
    logical :: fileExists
    character(len = STRING_LENGTH) :: message

    if (this%comm == MPI_COMM_NULL) return

    call MPI_Comm_rank(this%comm, procRank, ierror)

    if (procRank == 0) inquire(file = trim(this%probeFilename), exist = fileExists)
    call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, this%comm, ierror)

    if (.not. fileExists) then
       write(message, '(3A,I0.0,A)') "Probe file not available for patch '",                 &
            trim(this%name), "' on grid ", this%gridIndex, "!"
       call gracefulExit(this%comm, message)
    end if

    this%iProbeBuffer = size(this%probeBuffer, 3) + 1

    call MPI_File_open(this%comm, trim(this%probeFilename) // char(0),                       &
         MPI_MODE_WRONLY, MPI_INFO_NULL, mpiFileHandle, ierror)
    call MPI_File_get_size(mpiFileHandle, this%probeFileOffset, ierror)
    call MPI_File_close(mpiFileHandle, ierror)

  end subroutine seekToEOF

end module ProbePatch_mod
