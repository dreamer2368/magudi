#include "config.h"

module ActuatorPatchImpl

  implicit none
  public

contains

  subroutine loadActuatorGradient(this, gradientBuffer, mpiScalarSubarrayType)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use ActuatorPatch_mod, only : t_ActuatorPatch

    ! <<< Arguments >>>
    class(t_ActuatorPatch) :: this
    SCALAR_TYPE, intent(in), optional :: gradientBuffer(:,:,:)
    integer, intent(in), optional :: mpiScalarSubarrayType

    ! <<< Local variables >>>
    integer :: mpiScalarSubarrayType_, mpiFileHandle, ierror

    mpiScalarSubarrayType_ = this%mpiScalarSubarrayType
    if (present(mpiScalarSubarrayType)) mpiScalarSubarrayType_ = mpiScalarSubarrayType

    call MPI_File_open(this%comm, trim(this%gradientFilename) // char(0), MPI_MODE_RDONLY,   &
         MPI_INFO_NULL, mpiFileHandle, ierror)

    if (present(gradientBuffer)) then
       this%gradientFileOffset = this%gradientFileOffset -                                   &
            SIZEOF_SCALAR * product(int(this%globalSize, MPI_OFFSET_KIND)) *                 &
            size(gradientBuffer, 2) * size(gradientBuffer, 3)
    else
       this%gradientFileOffset = this%gradientFileOffset -                                   &
            SIZEOF_SCALAR * product(int(this%globalSize, MPI_OFFSET_KIND)) *                 &
            size(this%gradientBuffer, 2) * size(this%gradientBuffer, 3)
    end if
    assert(this%gradientFileOffset >= 0)

    call MPI_File_set_view(mpiFileHandle, this%gradientFileOffset, SCALAR_TYPE_MPI,          &
         mpiScalarSubarrayType_, "native", MPI_INFO_NULL, ierror)

    if (present(gradientBuffer)) then
       call MPI_File_read_all(mpiFileHandle, gradientBuffer, size(gradientBuffer),           &
            SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
    else
       call MPI_File_read_all(mpiFileHandle, this%gradientBuffer, size(this%gradientBuffer), &
            SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
    end if

    call MPI_File_close(mpiFileHandle, ierror)

  end subroutine loadActuatorGradient

  subroutine saveActuatorGradient(this, gradientBuffer, mpiScalarSubarrayType)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use ActuatorPatch_mod, only : t_ActuatorPatch

    ! <<< Arguments >>>
    class(t_ActuatorPatch) :: this
    SCALAR_TYPE, intent(in), optional :: gradientBuffer(:,:,:)
    integer, intent(in), optional :: mpiScalarSubarrayType

    ! <<< Local variables >>>
    integer :: mpiScalarSubarrayType_, mpiFileHandle, ierror

    mpiScalarSubarrayType_ = this%mpiScalarSubarrayType
    if (present(mpiScalarSubarrayType)) mpiScalarSubarrayType_ = mpiScalarSubarrayType

    call MPI_File_open(this%comm, trim(this%gradientFilename) // char(0), MPI_MODE_WRONLY,   &
         MPI_INFO_NULL, mpiFileHandle, ierror)
    call MPI_File_set_view(mpiFileHandle, this%gradientFileOffset, SCALAR_TYPE_MPI,          &
         mpiScalarSubarrayType_, "native", MPI_INFO_NULL, ierror)
    if (present(gradientBuffer)) then
       call MPI_File_write_all(mpiFileHandle, gradientBuffer, size(gradientBuffer),          &
            SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
       this%gradientFileOffset = this%gradientFileOffset +                                   &
            SIZEOF_SCALAR * product(int(this%globalSize, MPI_OFFSET_KIND)) *                 &
            size(gradientBuffer, 2) * size(gradientBuffer, 3)
    else
       call MPI_File_write_all(mpiFileHandle, this%gradientBuffer,                           &
            size(this%gradientBuffer), SCALAR_TYPE_MPI,                                      &
            MPI_STATUS_IGNORE, ierror)
       this%gradientFileOffset = this%gradientFileOffset +                                   &
            SIZEOF_SCALAR * product(int(this%globalSize, MPI_OFFSET_KIND)) *                 &
            size(this%gradientBuffer, 2) * size(this%gradientBuffer, 3)
    end if
    call MPI_File_close(mpiFileHandle, ierror)

  end subroutine saveActuatorGradient

end module ActuatorPatchImpl

subroutine setupActuatorPatch(this, index, comm, patchDescriptor,                            &
     grid, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_ActuatorPatch) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  character(len = STRING_LENGTH) :: key, outputPrefix
  integer :: gradientBufferSize, arrayOfSizes(5), arrayOfSubsizes(5), arrayOfStarts(5), ierror

#ifdef DEBUG
  if (.not. simulationFlags%predictionOnly) then
     assert(grid%nGridPoints > 0)
     assert(allocated(grid%controlMollifier))
     assert(size(grid%controlMollifier, 1) == grid%nGridPoints)
     assert(size(grid%controlMollifier, 2) == 1)
  end if
#endif

  call this%cleanup()
  call this%setupBase(index, comm, patchDescriptor, grid, simulationFlags, solverOptions)

  outputPrefix = getOption("output_prefix", PROJECT_NAME)
  write(this%gradientFilename, '(4A)') trim(outputPrefix), "-",                              &
       trim(patchDescriptor%name), "_gradient.dat"

  write(key, '(A)') "patches/" // trim(patchDescriptor%name) // "/"

  ! Interval for saving the gradient.
  gradientBufferSize = getOption("defaults/gradient_buffer_size", 1)
  gradientBufferSize = getOption(trim(key) // "gradient_buffer_size", gradientBufferSize)

  if (.not. simulationFlags%predictionOnly .and. this%nPatchPoints > 0) then

     allocate(this%mollifier(this%nPatchPoints))
     call this%collect(grid%controlMollifier(:,1), this%mollifier)
     allocate(this%controlForcing(this%nPatchPoints, solverOptions%nUnknowns))
     allocate(this%gradientBuffer(this%nPatchPoints, solverOptions%nUnknowns,                &
          gradientBufferSize))

     assert(this%comm /= MPI_COMM_NULL)

     if (this%mpiScalarSubarrayType /= MPI_DATATYPE_NULL)                                    &
          call MPI_Type_free(this%mpiScalarSubarrayType, ierror)
     this%mpiScalarSubarrayType = MPI_DATATYPE_NULL

     arrayOfSizes(1:3) = this%globalSize
     arrayOfSizes(4) = size(this%gradientBuffer, 2)
     arrayOfSizes(5) = size(this%gradientBuffer, 3)
     arrayOfSubsizes(1:3) = this%localSize
     arrayOfSubsizes(4) = size(this%gradientBuffer, 2)
     arrayOfSubsizes(5) = size(this%gradientBuffer, 3)
     arrayOfStarts(1:3) = this%offset - this%extent(1::2) + 1
     arrayOfStarts(4:5) = 0
     call MPI_Type_create_subarray(5, arrayOfSizes, arrayOfSubsizes, arrayOfStarts,          &
          MPI_ORDER_FORTRAN, SCALAR_TYPE_MPI, this%mpiScalarSubarrayType, ierror)
     call MPI_Type_commit(this%mpiScalarSubarrayType, ierror)

  end if

end subroutine setupActuatorPatch

subroutine cleanupActuatorPatch(this)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use ActuatorPatch_mod, only : t_ActuatorPatch

  implicit none

  ! <<< Arguments >>>
  class(t_ActuatorPatch) :: this

  call this%cleanupBase()

  SAFE_DEALLOCATE(this%mollifier)
  SAFE_DEALLOCATE(this%controlForcing)
  SAFE_DEALLOCATE(this%gradientBuffer)

  this%iGradientBuffer = 0
  this%gradientFileOffset = int(0, MPI_OFFSET_KIND)

end subroutine cleanupActuatorPatch

subroutine updateActuatorPatch(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Private members >>>
  use ActuatorPatchImpl, only : loadActuatorGradient, saveActuatorGradient

  implicit none

  ! <<< Arguments >>>
  class(t_ActuatorPatch) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, nDimensions, nUnknowns, gridIndex, patchIndex

  assert_key(mode, (FORWARD, ADJOINT))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  if (mode == FORWARD .and. abs(state%actuationAmount) <= 0.0_wp) return

  call startTiming("updateActuatorPatch")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2)

  select case (mode)

  case (FORWARD)

     assert(this%iGradientBuffer > 0)
     assert(this%iGradientBuffer <= size(this%gradientBuffer, 3))

     if (this%iGradientBuffer == size(this%gradientBuffer, 3)) call loadActuatorGradient(this)

     ! Temporarily hard-coded to be an energy actuation.
     this%controlForcing = 0.0_wp
     this%controlForcing(:,nDimensions+2) = - state%actuationAmount *                        &
          this%gradientBuffer(:,nDimensions+2,this%iGradientBuffer)

     this%iGradientBuffer = this%iGradientBuffer - 1
     if (this%iGradientBuffer == 0) this%iGradientBuffer = size(this%gradientBuffer, 3)

     do l = 1, nUnknowns
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

                 state%rightHandSide(gridIndex,l) = state%rightHandSide(gridIndex,l) +       &
                      grid%controlMollifier(gridIndex,1) * this%controlForcing(patchIndex,l)

              end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
           end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
        end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
     end do !... l = 1, nUnknowns

  case (ADJOINT)

     this%iGradientBuffer = this%iGradientBuffer + 1
     if (this%iGradientBuffer == size(this%gradientBuffer, 3) + 1) then
        call saveActuatorGradient(this)
        this%iGradientBuffer = 1
     end if

     assert(this%iGradientBuffer > 0)
     assert(this%iGradientBuffer <= size(this%gradientBuffer, 3))

     call this%collect(state%adjointVariables, this%gradientBuffer(:,:,this%iGradientBuffer))
     do l = 1, nUnknowns
        this%gradientBuffer(:,l,this%iGradientBuffer) = this%mollifier *                     &
             this%gradientBuffer(:,l,this%iGradientBuffer)
     end do

  end select

  call endTiming("updateActuatorPatch")

end subroutine updateActuatorPatch

function verifyActuatorPatchUsage(this, patchDescriptor, gridSize, normalDirection,          &
     extent, simulationFlags, success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_ActuatorPatch) :: this
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
  if (simulationFlags%predictionOnly) isPatchUsed = .false.

end function verifyActuatorPatchUsage

subroutine startBufferedActuatorGradientIO(this, mode)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use ActuatorPatch_mod, only : t_ActuatorPatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use InputHelper, only : getFreeUnit
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_ActuatorPatch) :: this
  integer, intent(in) :: mode

  ! <<< Local variables >>>
  character(len = STRING_LENGTH) :: message
  integer :: fileUnit, mpiFileHandle, procRank, ierror
  logical :: fileExists

  if (this%comm /= MPI_COMM_NULL) then

     call MPI_Comm_rank(this%comm, procRank, ierror)

     select case (mode)

     case (FORWARD)

        if (procRank == 0) inquire(file = trim(this%gradientFilename), exist = fileExists)
        call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, this%comm, ierror)
        if (.not. fileExists) then
           write(message, '(3A,I0.0,A)') "No gradient information available for patch '",    &
                trim(this%name), "' on grid ", this%gridIndex, "!"
           call gracefulExit(this%comm, message)
        end if
        this%iGradientBuffer = size(this%gradientBuffer, 3)
        call MPI_File_open(this%comm, trim(this%gradientFilename) // char(0),                &
             MPI_MODE_WRONLY, MPI_INFO_NULL, mpiFileHandle, ierror)
        call MPI_File_get_size(mpiFileHandle, this%gradientFileOffset, ierror)
        call MPI_File_close(mpiFileHandle, ierror)

     case (ADJOINT)

        if (procRank == 0) then
           open(unit = getFreeUnit(fileUnit), file = trim(this%gradientFilename),            &
                action = 'write', status = 'unknown')
           close(fileUnit)
        end if
        call MPI_Barrier(this%comm, ierror)
        this%iGradientBuffer = 0
        this%gradientFileOffset = int(0, MPI_OFFSET_KIND)

     end select

  end if

end subroutine startBufferedActuatorGradientIO

subroutine finishBufferedActuatorGradientIO(this, mode)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use ActuatorPatch_mod, only : t_ActuatorPatch

  ! <<< Enumerations >>>
  use Region_enum, only : ADJOINT

  ! <<< Private members >>>
  use ActuatorPatchImpl, only : saveActuatorGradient

  implicit none

  ! <<< Arguments >>>
  class(t_ActuatorPatch) :: this
  integer, intent(in) :: mode

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: arrayOfSizes(5), arrayOfSubsizes(5), arrayOfStarts(5),                          &
       mpiScalarSubarrayType, ierror

  select case (mode)

  case (ADJOINT)

     if (this%iGradientBuffer == size(this%gradientBuffer, 3)) then
        call saveActuatorGradient(this)
     else
        arrayOfSizes(1:3) = this%globalSize
        arrayOfSizes(4) = size(this%gradientBuffer, 2)
        arrayOfSizes(5) = this%iGradientBuffer
        arrayOfSubsizes(1:3) = this%localSize
        arrayOfSubsizes(4) = size(this%gradientBuffer, 2)
        arrayOfSubsizes(5) = this%iGradientBuffer
        arrayOfStarts(1:3) = this%offset - this%extent(1::2) + 1
        arrayOfStarts(4:5) = 0
        call MPI_Type_create_subarray(5, arrayOfSizes, arrayOfSubsizes, arrayOfStarts,       &
             MPI_ORDER_FORTRAN, SCALAR_TYPE_MPI, mpiScalarSubarrayType, ierror)
        call MPI_Type_commit(mpiScalarSubarrayType, ierror)
        call saveActuatorGradient(this, this%gradientBuffer(:,:,:this%iGradientBuffer),      &
             mpiScalarSubarrayType)
        call MPI_Type_free(mpiScalarSubarrayType, ierror)
        mpiScalarSubarrayType = MPI_DATATYPE_NULL
     end if

  end select

end subroutine finishBufferedActuatorGradientIO
