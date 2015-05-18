#include "config.h"

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
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: outputPrefix

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
  write(this%gradientFilename, '(4A)') trim(outputPrefix), ".gradient_",                     &
       trim(patchDescriptor%name), ".dat"

  if (.not. simulationFlags%predictionOnly .and. this%nPatchPoints > 0) then
     allocate(this%controlForcing(this%nPatchPoints, solverOptions%nUnknowns))
     this%controlForcing = 0.0_wp
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
  integer :: i, j, k, l, nDimensions, nUnknowns, nSpecies, gridIndex, patchIndex

  assert_key(mode, (FORWARD, ADJOINT))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  if (mode == FORWARD .and. abs(state%actuationAmount) <= 0.0_wp) return

  call startTiming("updateActuatorPatch")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nSpecies = solverOptions%nSpecies
  assert(nSpecies >= 0)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2 + nSpecies)

  select case (mode)

  case (FORWARD)

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

                 state%rightHandSide(gridIndex, l) = state%rightHandSide(gridIndex, l) +     &
                      grid%controlMollifier(gridIndex, 1) * this%controlForcing(patchIndex, l)

              end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
           end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
        end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
     end do !... l = 1, nUnknowns

  ! case (ADJOINT)


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

subroutine loadActuatorGradient(this)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use ActuatorPatch_mod, only : t_ActuatorPatch

  ! <<< Arguments >>>
  class(t_ActuatorPatch) :: this

  ! <<< Local variables >>>
  integer :: arrayOfSizes(5), arrayOfSubsizes(5), arrayOfStarts(5),                          &
       mpiScalarSubarrayType, mpiFileHandle, dataSize, ierror
  integer(kind = MPI_OFFSET_KIND) :: nBytesToRead

  if (this%comm == MPI_COMM_NULL) return

  arrayOfSizes(1:3) = this%globalSize
  arrayOfSizes(4) = size(this%gradientBuffer, 2)
  arrayOfSizes(5) = this%iGradientBuffer
  arrayOfSubsizes(1:3) = this%localSize
  arrayOfSubsizes(4) = size(this%gradientBuffer, 2)
  arrayOfSubsizes(5) = this%iGradientBuffer
  arrayOfStarts(1:3) = this%offset - this%extent(1::2) + 1
  arrayOfStarts(4:5) = 0
  call MPI_Type_create_subarray(5, arrayOfSizes, arrayOfSubsizes, arrayOfStarts,             &
       MPI_ORDER_FORTRAN, SCALAR_TYPE_MPI, mpiScalarSubarrayType, ierror)
  call MPI_Type_commit(mpiScalarSubarrayType, ierror)

  call MPI_File_open(this%comm, trim(this%gradientFilename) // char(0), MPI_MODE_RDONLY,     &
       MPI_INFO_NULL, mpiFileHandle, ierror)

  nBytesToRead = SIZEOF_SCALAR * product(int(this%globalSize, MPI_OFFSET_KIND)) *            &
       size(this%gradientBuffer, 2) * this%iGradientBuffer
  if (this%gradientFileOffset - nBytesToRead <= int(0, MPI_OFFSET_KIND)) then
     this%iGradientBuffer = int(this%gradientFileOffset / (SIZEOF_SCALAR * &
          product(int(this%globalSize, MPI_OFFSET_KIND)) * size(this%gradientBuffer, 2)))
     this%gradientFileOffset = 0
  else
     this%gradientFileOffset = this%gradientFileOffset - nBytesToRead
  end if

  assert(this%gradientFileOffset >= 0)

  call MPI_File_set_view(mpiFileHandle, this%gradientFileOffset, SCALAR_TYPE_MPI,            &
       mpiScalarSubarrayType, "native", MPI_INFO_NULL, ierror)

  dataSize = this%nPatchPoints * size(this%gradientBuffer, 2) * this%iGradientBuffer
  call MPI_File_read_all(mpiFileHandle, this%gradientBuffer, dataSize,                       &
       SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)

  call MPI_File_close(mpiFileHandle, ierror)

  call MPI_Type_free(mpiScalarSubarrayType, ierror)

end subroutine loadActuatorGradient

subroutine saveActuatorGradient(this)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use ActuatorPatch_mod, only : t_ActuatorPatch

  ! <<< Arguments >>>
  class(t_ActuatorPatch) :: this

  ! <<< Local variables >>>
  integer :: arrayOfSizes(5), arrayOfSubsizes(5), arrayOfStarts(5),                          &
       mpiScalarSubarrayType, mpiFileHandle, dataSize, ierror

  if (this%comm == MPI_COMM_NULL .or. this%iGradientBuffer == 0) return

  arrayOfSizes(1:3) = this%globalSize
  arrayOfSizes(4) = size(this%gradientBuffer, 2)
  arrayOfSizes(5) = this%iGradientBuffer
  arrayOfSubsizes(1:3) = this%localSize
  arrayOfSubsizes(4) = size(this%gradientBuffer, 2)
  arrayOfSubsizes(5) = this%iGradientBuffer
  arrayOfStarts(1:3) = this%offset - this%extent(1::2) + 1
  arrayOfStarts(4:5) = 0
  call MPI_Type_create_subarray(5, arrayOfSizes, arrayOfSubsizes, arrayOfStarts,             &
       MPI_ORDER_FORTRAN, SCALAR_TYPE_MPI, mpiScalarSubarrayType, ierror)
  call MPI_Type_commit(mpiScalarSubarrayType, ierror)

  call MPI_File_open(this%comm, trim(this%gradientFilename) // char(0), MPI_MODE_WRONLY,     &
       MPI_INFO_NULL, mpiFileHandle, ierror)

  assert(this%gradientFileOffset >= 0)

  call MPI_File_set_view(mpiFileHandle, this%gradientFileOffset, SCALAR_TYPE_MPI,            &
       mpiScalarSubarrayType, "native", MPI_INFO_NULL, ierror)

  dataSize = this%nPatchPoints * size(this%gradientBuffer, 2) * this%iGradientBuffer
  call MPI_File_write_all(mpiFileHandle, this%gradientBuffer, dataSize,                      &
       SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)

  this%gradientFileOffset = this%gradientFileOffset +                                        &
       SIZEOF_SCALAR * product(int(this%globalSize, MPI_OFFSET_KIND)) *                      &
       size(this%gradientBuffer, 2) * this%iGradientBuffer

  call MPI_File_close(mpiFileHandle, ierror)

  call MPI_Type_free(mpiScalarSubarrayType, ierror)

end subroutine saveActuatorGradient
