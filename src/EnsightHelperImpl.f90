#include "config.h"

module EnsightHelperImpl

  implicit none
  public

contains

  pure subroutine swapIntegerEndianness_(a)

    ! <<< Arguments >>>
    integer, intent(inout) :: a

    ! <<< Local variables >>>
    integer :: b

    call mvbits(a,  0, 8, b, 24)
    call mvbits(a,  8, 8, b, 16)
    call mvbits(a, 16, 8, b,  8)
    call mvbits(a, 24, 8, b,  0)

    a = b

  end subroutine swapIntegerEndianness_

  pure subroutine swapScalarEndianness_(a)

    ! <<< Arguments >>>
    SCALAR_TYPE, intent(inout) :: a

    ! <<< Local variables >>>
    integer :: i
    integer(kind = 1) :: b(SIZEOF_SCALAR), c(SIZEOF_SCALAR)

    b = transfer(a, b)
    do i = 1, SIZEOF_SCALAR
       c(i) = b(SIZEOF_SCALAR+1-i)
    end do
    a = transfer(c, a)

  end subroutine swapScalarEndianness_

end module EnsightHelperImpl

subroutine plot3dDetectFormat(comm, filename, success, descriptor,                           &
     globalGridSizes, includeFunctionFiles)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_c_binding

  ! <<< Derived types >>>
  use PLOT3DDescriptor_type

  ! <<< Public members >>>
  use PLOT3DHelper, only : plot3dErrorMessage

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: filename
  logical, intent(out) :: success
  type(t_PLOT3DDescriptor), intent(out), optional :: descriptor
  integer, allocatable, intent(out), optional :: globalGridSizes(:,:)
  logical, intent(in), optional :: includeFunctionFiles

  interface

     function detect_format(filename, includeFunctionFiles, nGrids, nDimensions, nScalars,   &
          hasIblank, isEndiannessNative, fileType, globalGridSizes)                          &
          bind(C, name = "detect_format")
       use, intrinsic :: iso_c_binding
       integer(kind = C_INT) :: detect_format
       character(kind = C_CHAR) :: filename(STRING_LENGTH + 1)
       integer(kind = C_INT), value :: includeFunctionFiles
       integer(kind = C_INT32_T) :: nGrids
       integer(kind = C_INT) :: nDimensions, nScalars, isEndiannessNative, hasIblank, fileType
       type(C_PTR) :: globalGridSizes
     end function detect_format

     subroutine free_memory(globalGridSizes) bind(C, name = "free_memory")
       use, intrinsic :: iso_c_binding
       type(C_PTR) :: globalGridSizes
     end subroutine free_memory

  end interface

  ! <<< Local variables >>>
  integer :: i, j, procRank, ierror
  integer(C_INT) :: nDimensions, nScalars, errorCode, includeFunctionFiles_,                 &
       isEndiannessNative, hasIblank, file_type
  character(C_CHAR) :: filename_(STRING_LENGTH + 1)
  type(C_PTR) :: globalGridSizes_ptr
  integer(kind = C_INT32_T) :: nGrids
  integer(C_INT), pointer :: globalGridSizes_(:)

  ! By default, exclude function files.
  includeFunctionFiles_ = 0
  if (present(includeFunctionFiles)) then
     if (includeFunctionFiles) includeFunctionFiles_ = 1
  end if

  ! Get my rank in `comm`.
  call MPI_Comm_rank(comm, procRank, ierror)

  ! Only the ``master'' process in `comm` calls the function that implements the detection
  ! algorithm and broadcasts an error code to other processes.
  if (procRank == 0) then
     do i = 1, len_trim(filename)
        filename_(i) = filename(i:i)
     end do
     filename_(len_trim(filename)+1) = char(0)
     errorCode = detect_format(filename_, includeFunctionFiles_, nGrids, nDimensions,        &
          nScalars, isEndiannessNative, hasIblank, file_type, globalGridSizes_ptr)
     if (.not. present(globalGridSizes)) call free_memory(globalGridSizes_ptr)
  end if
  call MPI_Bcast(errorCode, 1, MPI_INT, 0, comm, ierror)

  ! If an error occured, update the static `plot3dErrorMessage` variable and return.
  if (errorCode /= 0) then
     select case (errorCode)
     case (-1)
        write(plot3dErrorMessage, '(2A)') trim(filename),                                    &
             ": File not found or permission denied."
     case (-2)
        write(plot3dErrorMessage, '(2A)') trim(filename), ": Unexpected end of file."
     case (-3)
        write(plot3dErrorMessage, '(2A)') trim(filename),                                    &
             ": Not a valid multi-block whole-format PLOT3D file."
     case (-4)
        write(plot3dErrorMessage, '(2A)') trim(filename), ": Inconsistent record markers."
     end select
     success = .false.
     return
  end if

  ! If not, broadcast the number of dimensions and number of grids to all processes.
  call MPI_Bcast(nDimensions, 1, MPI_INT, 0, comm, ierror)
  call MPI_Bcast(nGrids, 1, MPI_INT32_T, 0, comm, ierror)

#ifdef DEBUG
  if (nDimensions < 0 .or. nDimensions > 3) then
     write(plot3dErrorMessage, '(3A,2(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": detect_format returned an invalid number of dimensions: ",                      &
          nDimensions, "!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (nGrids < 0) then
     write(plot3dErrorMessage, '(3A,2(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": detect_format returned an invalid number of grids: ",                           &
          nGrids, "!"
     success = .false.
     return
  end if
#endif

  if (present(descriptor)) then

     ! Fill the descriptor array with useful information.
     descriptor%nDimensions = nDimensions
     descriptor%nGrids = int(nGrids, kind = kind(descriptor%nGrids))
     call MPI_Bcast(isEndiannessNative, 1, MPI_INT, 0, comm, ierror)
     descriptor%isEndiannessNative = (isEndiannessNative /= 0)
     call MPI_Bcast(hasIblank, 1, MPI_INT, 0, comm, ierror)
     descriptor%hasIblank = (hasIblank /= 0)
     call MPI_Bcast(file_type, 1, MPI_INT, 0, comm, ierror)

     ! Character identifier for the file type.
     select case (file_type)
     case (0)
        descriptor%fileType = PLOT3D_GRID_FILE
     case (1)
        descriptor%fileType = PLOT3D_SOLUTION_FILE
     case (2)
        descriptor%fileType = PLOT3D_FUNCTION_FILE
     case default
#ifdef DEBUG
        write(plot3dErrorMessage, '(3A,2(I0.0,A))') "In ", __FILE__, ":", __LINE__,          &
             ": detect_format returned an invalid file type: ",                              &
             file_type, "!"
        success = .false.
        return
#endif
     end select

     ! If this is a function file, broadcast the number of components to all processes.
     if (descriptor%fileType == PLOT3D_FUNCTION_FILE) then
        call MPI_Bcast(nScalars, 1, MPI_INT, 0, comm, ierror)
#ifdef DEBUG
        if (nScalars < 0) then
           write(plot3dErrorMessage, '(3A,2(I0.0,A))') "In ", __FILE__, ":", __LINE__,       &
                ": detect_format returned an invalid number of scalars: ",                   &
                nScalars, "!"
           success = .false.
           return
        end if
#endif
        descriptor%nScalars = nScalars
     else
        descriptor%nScalars = 0
     end if

  end if

  ! If requested, read and broadcast the grid dimensions:

  if (present(globalGridSizes)) then

     SAFE_DEALLOCATE(globalGridSizes)
     allocate(globalGridSizes(nDimensions, nGrids))

     if (procRank == 0) then
        call C_F_POINTER(globalGridSizes_ptr, globalGridSizes_, [nDimensions * nGrids])
        do j = 1, size(globalGridSizes, 2)
           do i = 1, size(globalGridSizes, 1)
              globalGridSizes(i,j) = globalGridSizes_(i + nDimensions * (j - 1))
           end do
        end do
        call free_memory(globalGridSizes_ptr)
     end if

     call MPI_Bcast(globalGridSizes, size(globalGridSizes), MPI_INTEGER, 0, comm, ierror)

  end if

  ! Successful return code.
  success = .true.

end subroutine plot3dDetectFormat

function plot3dGetOffsetFromGridSizes_(fileType, globalGridSizes, gridIndex,                 &
     hasIblank, success, nScalars) result(offset)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use PLOT3DDescriptor_type

  ! <<< Public members >>>
  use PLOT3DHelper, only : plot3dErrorMessage

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: fileType, globalGridSizes(:,:), gridIndex
  logical, intent(in) :: hasIblank
  logical, intent(out) :: success
  integer, intent(in), optional :: nScalars

  ! <<< Result >>>
  integer(kind = MPI_OFFSET_KIND) :: offset

  ! <<< Local variables >>>
  integer :: nScalars_

  offset = 0

  nScalars_ = 1
  if (present(nScalars)) nScalars_ = nScalars
#ifdef DEBUG
  if (nScalars_ < 0) then
     write(plot3dErrorMessage, '(3A,2(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value for argument nScalars: ", nScalars_, "!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (size(globalGridSizes, 2) <= 0 .or. size(globalGridSizes, 1) <= 0 .or.                  &
       size(globalGridSizes, 1) > 3) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid shape of array argument globalGridSizes: (",                            &
          size(globalGridSizes, 1), ", ",                                                    &
          size(globalGridSizes, 2), ")!"
     success = .false.
     return
  end if
  if (any(globalGridSizes <= 0)) then
     write(plot3dErrorMessage, '(3A,I0.0,A)') "In ", __FILE__, ":", __LINE__,                &
          ": Invalid array argument globalGridSizes (some values were non-negative)!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (gridIndex < 0 .or. gridIndex > size(globalGridSizes, 2)) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value for argument gridIndex: ", gridIndex,                             &
          ". Expected a positive integer <= ",                                               &
          size(globalGridSizes, 2), "!"
     success = .false.
     return
  end if
#endif

  select case (fileType)

  case (PLOT3D_GRID_FILE)
     if (hasIblank) then
        offset = 4 + 5 * SIZEOF_PLOT3D_OFF +                                                 &
             12 * int(size(globalGridSizes, 2), MPI_OFFSET_KIND) +                           &
             2 * SIZEOF_PLOT3D_OFF * (int(gridIndex, MPI_OFFSET_KIND) - 1) +                 &
             (3 * SIZEOF_SCALAR + 4) *                                                       &
             sum(product(int(globalGridSizes(:,:gridIndex-1), MPI_OFFSET_KIND), dim = 1))
     else
        offset = 4 + 5 * SIZEOF_PLOT3D_OFF +                                                 &
             12 * int(size(globalGridSizes, 2), MPI_OFFSET_KIND) +                           &
             2 * SIZEOF_PLOT3D_OFF * (int(gridIndex, MPI_OFFSET_KIND) - 1) +                 &
             3 * SIZEOF_SCALAR *                                                             &
             sum(product(int(globalGridSizes(:,:gridIndex-1), MPI_OFFSET_KIND), dim = 1))
     end if

  case (PLOT3D_SOLUTION_FILE)
     offset = 4 + 5 * SIZEOF_PLOT3D_OFF +                                                    &
          12 * int(size(globalGridSizes, 2), MPI_OFFSET_KIND) +                              &
          4 * SIZEOF_PLOT3D_OFF * (int(gridIndex, MPI_OFFSET_KIND) - 1) +                    &
          (5 * SIZEOF_SCALAR) *                                                              &
          sum(product(int(globalGridSizes(:,:gridIndex-1), MPI_OFFSET_KIND), dim = 1)) +     &
          4 * SIZEOF_SCALAR * (int(gridIndex, MPI_OFFSET_KIND) - 1)

  case (PLOT3D_FUNCTION_FILE)
     offset = 4 + 5 * SIZEOF_PLOT3D_OFF +                                                    &
          16 * int(size(globalGridSizes, 2), MPI_OFFSET_KIND) +                              &
          2 * SIZEOF_PLOT3D_OFF * (int(gridIndex, MPI_OFFSET_KIND) - 1) +                    &
          (int(nScalars_, MPI_OFFSET_KIND) * SIZEOF_SCALAR) *                                &
          sum(product(int(globalGridSizes(:,:gridIndex-1), MPI_OFFSET_KIND), dim = 1))

  case default
#ifdef DEBUG
     write(plot3dErrorMessage, '(3A,2(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value for argument fileType: ", fileType, "!"
     success = .false.
     return
#endif

  end select

  ! Successful return code.
  success = .true.

end function plot3dGetOffsetFromGridSizes_

function plot3dGetOffsetFromFile_(comm, filename, gridIndex, success) result(offset)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use PLOT3DDescriptor_type

  ! <<< Public members >>>
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dGetOffset, plot3dErrorMessage

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: filename
  integer, intent(in) :: gridIndex
  logical, intent(out) :: success

  ! <<< Result >>>
  integer(kind = MPI_OFFSET_KIND) :: offset

  ! <<< Local variables >>>
  integer, allocatable :: globalGridSizes(:,:)
  type(t_PLOT3DDescriptor) :: descriptor

  offset = 0

  call plot3dDetectFormat(comm, filename, success, descriptor,                               &
       globalGridSizes, includeFunctionFiles = .true.)
  if (.not. success) then
     return
  end if

#ifdef DEBUG
  if (gridIndex < 0 .or. gridIndex > size(globalGridSizes, 2)) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value for argument gridIndex: ", gridIndex,                             &
          ". Expected a positive integer <= ",                                               &
          size(globalGridSizes, 2), "!"
     success = .false.
     return
  end if
#endif

  offset = plot3dGetOffset(descriptor%fileType, globalGridSizes, gridIndex,                  &
       descriptor%hasIblank, success, descriptor%nScalars)
  SAFE_DEALLOCATE(globalGridSizes)

  ! Successful return code.
  success = .true.

end function plot3dGetOffsetFromFile_

subroutine plot3dWriteSkeleton(comm, filename, fileType, globalGridSizes, success, nScalars)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_c_binding

  ! <<< Derived types >>>
  use PLOT3DDescriptor_type

  ! <<< Public members >>>
  use PLOT3DHelper, only : plot3dErrorMessage

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: filename
  integer, intent(in) :: fileType, globalGridSizes(:,:)
  logical, intent(out) :: success
  integer, intent(in), optional :: nScalars

  interface
     function write_skeleton(filename, nGrids, nComp, fileType, globalGridSizes)             &
          bind(C, name = "write_skeleton")
       use, intrinsic :: iso_c_binding
       integer(kind = C_INT) :: write_skeleton
       character(kind = C_CHAR) :: filename(STRING_LENGTH + 1)
       integer(kind = C_INT32_T), value :: nGrids, nComp
       integer(kind = C_INT), value :: fileType
       integer(kind = C_INT32_T) :: globalGridSizes(*)
     end function write_skeleton
  end interface

  ! <<< Local variables >>>
  integer :: i, idim, iGrid, procRank, ierror
  integer(C_INT) :: errorCode, fileType_
  character(C_CHAR) :: filename_(STRING_LENGTH + 1)
  integer(C_INT32_T) :: nGrids, nScalars_
  integer(C_INT32_T), pointer :: globalGridSizes_(:)

  call MPI_Comm_rank(comm, procRank, ierror)

#ifdef DEBUG
  if (fileType == PLOT3D_FUNCTION_FILE .and. .not. present(nScalars)) then
     write(plot3dErrorMessage, '(3A,I0.0,2A)') "In ", __FILE__, ":", __LINE__,               &
          ": Unable to write PLOT3D function file skeleton",                                 &
          ": required argument nScalars not specified!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (fileType == PLOT3D_FUNCTION_FILE .and. present(nScalars)) then
     if (nScalars <= 0) then
        write(plot3dErrorMessage, '(3A,2(I0.0,A))') "In ", __FILE__, ":", __LINE__,          &
             ": Invalid value for argument nScalars: ", nScalars, "!"
        success = .false.
        return
     end if
  end if
#endif

#ifdef DEBUG
  if (size(globalGridSizes, 2) <= 0 .or. size(globalGridSizes, 1) <= 0 .or.                  &
       size(globalGridSizes, 1) > 3) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid shape of array argument globalGridSizes: (",                            &
          size(globalGridSizes, 1), ", ",                                                    &
          size(globalGridSizes, 2), ")!"
     success = .false.
     return
  end if
  if (any(globalGridSizes <= 0)) then
     write(plot3dErrorMessage, '(3A,I0.0,A)') "In ", __FILE__, ":", __LINE__,                &
          ": Invalid array argument globalGridSizes (some values were non-negative)!"
     success = .false.
     return
  end if
#endif

  select case (fileType)
  case (PLOT3D_GRID_FILE)
     fileType_ = 0
  case (PLOT3D_SOLUTION_FILE)
     fileType_ = 1
  case (PLOT3D_FUNCTION_FILE)
     fileType_ = 2
  case default
#ifdef DEBUG
     write(plot3dErrorMessage, '(3A,2(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value for argument fileType: ", fileType, "!"
     success = .false.
     return
#endif
  end select

  ! Only the ``master'' process in `comm` calls the function that writes the skeleton of the
  ! PLOT3D file and broadcasts an error code to other processes.
  if (procRank == 0) then

     nGrids = size(globalGridSizes, 2)
     nScalars_ = 0
     if (fileType == PLOT3D_FUNCTION_FILE) nScalars_ = nScalars

     do i = 1, len_trim(filename)
        filename_(i) = filename(i:i)
     end do
     filename_(len_trim(filename)+1) = char(0)

     allocate(globalGridSizes_(3 * nGrids))
     globalGridSizes_ = 1

     do iGrid = 1, size(globalGridSizes, 2)
        do idim = 1, size(globalGridSizes, 1)
           globalGridSizes_(idim + 3 * (iGrid - 1)) = globalGridSizes(idim, iGrid)
        end do
     end do

     errorCode = write_skeleton(filename_, nGrids, nScalars_, fileType_, globalGridSizes_)
     deallocate(globalGridSizes_)

  end if
  call MPI_Bcast(errorCode, 1, MPI_INT, 0, comm, ierror)

  ! If an error occured, print an error message and abort.
  if (errorCode /= 0) then
     select case (errorCode)
     case (-1)
        write(plot3dErrorMessage, '(2A)') trim(filename), ": Could not open file for writing."
     case default
        write(plot3dErrorMessage, '(2A)') trim(filename), ": Failed to write to file."
     end select
     success = .false.
     return
  end if

  call MPI_Barrier(comm, ierror)

  ! Successful return code.
  success = .true.

end subroutine plot3dWriteSkeleton

subroutine plot3dWriteSingleGrid(comm, filename, offset, mpiDerivedTypeScalarSubarray,       &
     mpiDerivedTypeIntegerSubarray, globalGridSize, coordinates, iblank, success)

  ! <<< External modules >>>
  use MPI

  ! <<< Public members >>>
  use PLOT3DHelper, only : plot3dErrorMessage

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: filename
  integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
  integer, intent(in) :: mpiDerivedTypeScalarSubarray,                                       &
       mpiDerivedTypeIntegerSubarray,                                                        &
       globalGridSize(3)
  real(SCALAR_KIND), intent(in) :: coordinates(:,:)
  integer, intent(in) :: iblank(:)
  logical, intent(out) :: success

  ! <<< Local variables >>>
  integer :: i, mpiFileHandle, ierror

#ifdef DEBUG
  if (offset < 0) then
     write(plot3dErrorMessage, '(3A,2(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value of argument offset: ", offset,                                    &
          ". Expected a non-negative value!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (size(coordinates, 1) <= 0 .or. size(coordinates, 2) <= 0 .or.                          &
       size(coordinates, 2) > 3) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid shape of array argument coordinates: (",                                &
          size(coordinates, 1), ", ",                                                        &
          size(coordinates, 2), ")!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (size(iblank) /= size(coordinates, 1)) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid shape of array argument iblank: (", size(iblank),                       &
          "). Expected (", size(coordinates, 1), ")!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (any(globalGridSize <= 0)) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value of argument globalGridSize: (/ ",                                 &
          globalGridSize(1), ", ", globalGridSize(2), ", ",                                  &
          globalGridSize(3), " /)!"
     success = .false.
     return
  end if
#endif

  call MPI_File_open(comm, trim(filename) // char(0), MPI_MODE_WRONLY,                       &
       MPI_INFO_NULL, mpiFileHandle, ierror)

  do i = 1, size(coordinates, 2)
     call MPI_File_set_view(mpiFileHandle, offset, SCALAR_TYPE_MPI,                          &
          mpiDerivedTypeScalarSubarray, "native",                                            &
          MPI_INFO_NULL, ierror)
     call MPI_File_write_all(mpiFileHandle, coordinates(:,i), size(coordinates, 1),          &
          SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
     offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
  end do

  do while (i <= 3)
     offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
     i = i + 1
  end do

  call MPI_File_set_view(mpiFileHandle, offset, MPI_INTEGER, mpiDerivedTypeIntegerSubarray,  &
       "native", MPI_INFO_NULL, ierror)
  call MPI_File_write_all(mpiFileHandle, iblank, size(iblank),                               &
       MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
  offset = offset + 4 * product(int(globalGridSize, MPI_OFFSET_KIND)) + 2 * SIZEOF_PLOT3D_OFF

  call MPI_File_close(mpiFileHandle, ierror)

  ! Successful return code.
  success = .true.

end subroutine plot3dWriteSingleGrid

subroutine plot3dWriteSingleAuxiliarySolutionData(comm, filename,                            &
     offset, auxiliarySolutionData, success)

  ! <<< External modules >>>
  use MPI

  ! <<< Public members >>>
  use PLOT3DHelper, only : plot3dErrorMessage

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: filename
  integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
  real(SCALAR_KIND), intent(in) :: auxiliarySolutionData(4)
  logical, intent(out) :: success

  ! <<< Local variables >>>
  integer :: mpiFileHandle, ierror

#ifdef DEBUG
  if (offset < 0) then
     write(plot3dErrorMessage, '(3A,2(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value of argument offset: ", offset,                                    &
          ". Expected a non-negative value!"
     success = .false.
     return
  end if
#endif

  call MPI_File_open(comm, trim(filename) // char(0), MPI_MODE_WRONLY,                       &
       MPI_INFO_NULL, mpiFileHandle, ierror)

  call MPI_File_set_view(mpiFileHandle, offset, SCALAR_TYPE_MPI,                             &
       SCALAR_TYPE_MPI, "native", MPI_INFO_NULL, ierror)
  call MPI_File_write_all(mpiFileHandle, auxiliarySolutionData, 4,                           &
       SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
  offset = offset + 4 * SIZEOF_SCALAR + 2 * SIZEOF_PLOT3D_OFF

  call MPI_File_close(mpiFileHandle, ierror)

  ! Successful return code.
  success = .true.

end subroutine plot3dWriteSingleAuxiliarySolutionData

subroutine plot3dWriteSingleSolution(comm, filename, offset, mpiDerivedTypeScalarSubarray,   &
     globalGridSize, solutionVector, success)

  ! <<< External modules >>>
  use MPI

  ! <<< Public members >>>
  use PLOT3DHelper, only : plot3dErrorMessage

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: filename
  integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
  integer, intent(in) :: mpiDerivedTypeScalarSubarray, globalGridSize(3)
  real(SCALAR_KIND), intent(in) :: solutionVector(:,:)
  logical, intent(out) :: success

  ! <<< Local variables >>>
  integer :: i, mpiFileHandle, ierror

#ifdef DEBUG
  if (offset < 0) then
     write(plot3dErrorMessage, '(3A,2(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value of argument offset: ", offset,                                    &
          ". Expected a non-negative value!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (size(solutionVector, 1) <= 0 .or. size(solutionVector, 2) < 3 .or.                     &
       size(solutionVector, 2) > 5) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid shape of array argument solutionVector: (",                             &
          size(solutionVector, 1), ", ",                                                     &
          size(solutionVector, 2), ")!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (any(globalGridSize <= 0)) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value of argument globalGridSize: (/ ",                                 &
          globalGridSize(1), ", ", globalGridSize(2), ", ",                                  &
          globalGridSize(3), " /)!"
     success = .false.
     return
  end if
#endif

  call MPI_File_open(comm, trim(filename) // char(0), MPI_MODE_WRONLY,                       &
       MPI_INFO_NULL, mpiFileHandle, ierror)

  do i = 1, size(solutionVector, 2) - 1
     call MPI_File_set_view(mpiFileHandle, offset, SCALAR_TYPE_MPI,                          &
          mpiDerivedTypeScalarSubarray, "native",                                            &
          MPI_INFO_NULL, ierror)
     call MPI_File_write_all(mpiFileHandle, solutionVector(:,i), size(solutionVector, 1),    &
          SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
     offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
  end do

  do while (i <= 4)
     offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
     i = i + 1
  end do

  i = size(solutionVector, 2)
  call MPI_File_set_view(mpiFileHandle, offset, SCALAR_TYPE_MPI,                             &
       mpiDerivedTypeScalarSubarray, "native", MPI_INFO_NULL, ierror)
  call MPI_File_write_all(mpiFileHandle, solutionVector(:,i), size(solutionVector, 1),       &
       SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
  offset = offset + SIZEOF_SCALAR *                                                          &
       product(int(globalGridSize, MPI_OFFSET_KIND)) + 2 * SIZEOF_PLOT3D_OFF

  call MPI_File_close(mpiFileHandle, ierror)

  ! Successful return code.
  success = .true.

end subroutine plot3dWriteSingleSolution

subroutine plot3dWriteSingleFunction(comm, filename, offset, mpiDerivedTypeScalarSubarray,   &
     globalGridSize, functionVector, success)

  ! <<< External modules >>>
  use MPI

  ! <<< Public members >>>
  use PLOT3DHelper, only : plot3dErrorMessage

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: filename
  integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
  integer, intent(in) :: mpiDerivedTypeScalarSubarray, globalGridSize(3)
  real(SCALAR_KIND), intent(in) :: functionVector(:,:)
  logical, intent(out) :: success

  ! <<< Local variables >>>
  integer :: i, mpiFileHandle, ierror

#ifdef DEBUG
  if (offset < 0) then
     write(plot3dErrorMessage, '(3A,2(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value of argument offset: ", offset,                                    &
          ". Expected a non-negative value!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (size(functionVector, 1) <= 0 .or. size(functionVector, 2) <= 0) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid shape of array argument functionVector: (",                             &
          size(functionVector, 1), ", ",                                                     &
          size(functionVector, 2), ")!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (any(globalGridSize <= 0)) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value of argument globalGridSize: (/ ",                                 &
          globalGridSize(1), ", ", globalGridSize(2), ", ",                                  &
          globalGridSize(3), " /)!"
     success = .false.
     return
  end if
#endif

  call MPI_File_open(comm, trim(filename) // char(0), MPI_MODE_WRONLY,                       &
       MPI_INFO_NULL, mpiFileHandle, ierror)

  do i = 1, size(functionVector, 2)
     call MPI_File_set_view(mpiFileHandle, offset, SCALAR_TYPE_MPI,                          &
          mpiDerivedTypeScalarSubarray, "native",                                            &
          MPI_INFO_NULL, ierror)
     call MPI_File_write_all(mpiFileHandle, functionVector(:,i), size(functionVector, 1),    &
          SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
     offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
  end do
  offset = offset + 2 * SIZEOF_PLOT3D_OFF

  call MPI_File_close(mpiFileHandle, ierror)

  ! Successful return code.
  success = .true.

end subroutine plot3dWriteSingleFunction

subroutine plot3dReadSingleGrid(comm, filename, offset, mpiDerivedTypeScalarSubarray,        &
     mpiDerivedTypeIntegerSubarray, globalGridSize, coordinates, iblank, success)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use PLOT3DDescriptor_type

  ! <<< Private members >>>
  use PLOT3DHelperImpl, only : swapEndianness

  ! <<< Public members >>>
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: filename
  integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
  integer, intent(in) :: mpiDerivedTypeScalarSubarray,                                       &
       mpiDerivedTypeIntegerSubarray,                                                        &
       globalGridSize(3)
  real(SCALAR_KIND), intent(out) :: coordinates(:,:)
  integer, intent(out) :: iblank(:)
  logical, intent(out) :: success

  ! <<< Local variables >>>
  integer :: i, idim, procRank, mpiFileHandle, ierror
  type(t_PLOT3DDescriptor) :: descriptor

#ifdef DEBUG
  if (offset < 0) then
     write(plot3dErrorMessage, '(3A,2(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value of argument offset: ", offset,                                    &
          ". Expected a non-negative value!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (size(coordinates, 1) <= 0 .or. size(coordinates, 2) <= 0 .or.                          &
       size(coordinates, 2) > 3) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid shape of array argument coordinates: (",                                &
          size(coordinates, 1), ", ",                                                        &
          size(coordinates, 2), ")!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (size(iblank) /= size(coordinates, 1)) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid shape of array argument iblank: (", size(iblank),                       &
          "). Expected (", size(coordinates, 1), ")!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (any(globalGridSize <= 0)) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value of argument globalGridSize: (/ ",                                 &
          globalGridSize(1), ", ", globalGridSize(2), ", ",                                  &
          globalGridSize(3), " /)!"
     success = .false.
     return
  end if
#endif

  call MPI_Comm_rank(comm, procRank, ierror)
  call plot3dDetectFormat(comm, filename, success, descriptor)
  if (.not. success) return

  if (descriptor%fileType /= PLOT3D_GRID_FILE) then
     write(plot3dErrorMessage, '(3A)') "'", trim(filename), "': Not a valid PLOT3D grid file."
     success = .false.
     return
  end if

  call MPI_File_open(comm, trim(filename) // char(0), MPI_MODE_RDONLY,                       &
       MPI_INFO_NULL, mpiFileHandle, ierror)

  do idim = 1, min(size(coordinates, 2), 3)
     call MPI_File_set_view(mpiFileHandle, offset, SCALAR_TYPE_MPI,                          &
          mpiDerivedTypeScalarSubarray, "native",                                            &
          MPI_INFO_NULL, ierror)
     call MPI_File_read_all(mpiFileHandle, coordinates(:,idim), size(coordinates, 1),        &
          SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
     if (.not. descriptor%isEndiannessNative) then
        do i = 1, size(coordinates, 1)
           call swapEndianness(coordinates(i, idim))
        end do
     end if
     offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
  end do

  do while (idim <= 3)
     offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
     idim = idim + 1
  end do

  if (descriptor%hasIblank) then
     call MPI_File_set_view(mpiFileHandle, offset, MPI_INTEGER,                              &
          mpiDerivedTypeIntegerSubarray, "native",                                           &
          MPI_INFO_NULL, ierror)
     call MPI_File_read_all(mpiFileHandle, iblank, size(iblank),                             &
          MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
     if (.not. descriptor%isEndiannessNative) then
        do i = 1, size(iblank)
           call swapEndianness(iblank(i))
        end do
     end if
     offset = offset + product(int(globalGridSize, MPI_OFFSET_KIND)) + 2 * SIZEOF_PLOT3D_OFF
  else
     iblank = 1
  end if

  call MPI_File_close(mpiFileHandle, ierror)

  ! Successful return code.
  success = .true.

end subroutine plot3dReadSingleGrid

subroutine plot3dReadSingleAuxiliarySolutionData(comm, filename,                             &
     offset, auxiliarySolutionData, success)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use PLOT3DDescriptor_type

  ! <<< Private members >>>
  use PLOT3DHelperImpl, only : swapEndianness

  ! <<< Public members >>>
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: filename
  integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
  real(SCALAR_KIND), intent(out) :: auxiliarySolutionData(4)
  logical, intent(out) :: success

  ! <<< Local variables >>>
  integer :: i, procRank, mpiFileHandle, ierror
  type(t_PLOT3DDescriptor) :: descriptor

#ifdef DEBUG
  if (offset < 0) then
     write(plot3dErrorMessage, '(3A,2(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value of argument offset: ", offset,                                    &
          ". Expected a non-negative value!"
     success = .false.
     return
  end if
#endif

  call MPI_Comm_rank(comm, procRank, ierror)
  call plot3dDetectFormat(comm, filename, success, descriptor)
  if (.not. success) return

  if (descriptor%fileType /= PLOT3D_SOLUTION_FILE) then
     write(plot3dErrorMessage, '(3A)') "'", trim(filename),                                  &
          "': Not a valid PLOT3D solution file."
     success = .false.
     return
  end if

  call MPI_File_open(comm, trim(filename) // char(0), MPI_MODE_RDONLY,                       &
       MPI_INFO_NULL, mpiFileHandle, ierror)

  call MPI_File_set_view(mpiFileHandle, offset, SCALAR_TYPE_MPI,                             &
       SCALAR_TYPE_MPI, "native", MPI_INFO_NULL, ierror)
  call MPI_File_read_all(mpiFileHandle, auxiliarySolutionData, 4,                            &
       SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
  if (.not. descriptor%isEndiannessNative) then
     do i = 1, 4
        call swapEndianness(auxiliarySolutionData(i))
     end do
  end if
  offset = offset + 4 * SIZEOF_SCALAR + 2 * SIZEOF_PLOT3D_OFF

  call MPI_File_close(mpiFileHandle, ierror)

  ! Successful return code.
  success = .true.

end subroutine plot3dReadSingleAuxiliarySolutionData

subroutine plot3dReadSingleSolution(comm, filename, offset, mpiDerivedTypeScalarSubarray,    &
     globalGridSize, solutionVector, success)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use PLOT3DDescriptor_type

  ! <<< Private members >>>
  use PLOT3DHelperImpl, only : swapEndianness

  ! <<< Public members >>>
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: filename
  integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
  integer, intent(in) :: mpiDerivedTypeScalarSubarray, globalGridSize(3)
  real(SCALAR_KIND), intent(out) :: solutionVector(:,:)
  logical, intent(out) :: success

  ! <<< Local variables >>>
  integer :: i, j, procRank, mpiFileHandle, ierror
  type(t_PLOT3DDescriptor) :: descriptor

#ifdef DEBUG
  if (offset < 0) then
     write(plot3dErrorMessage, '(3A,2(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value of argument offset: ", offset,                                    &
          ". Expected a non-negative value!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (size(solutionVector, 1) <= 0 .or. size(solutionVector, 2) < 3 .or.                     &
       size(solutionVector, 2) > 5) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid shape of array argument solutionVector: (",                             &
          size(solutionVector, 1), ", ",                                                     &
          size(solutionVector, 2), ")!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (any(globalGridSize <= 0)) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value of argument globalGridSize: (/ ",                                 &
          globalGridSize(1), ", ", globalGridSize(2), ", ",                                  &
          globalGridSize(3), " /)!"
     success = .false.
     return
  end if
#endif

  call MPI_Comm_rank(comm, procRank, ierror)
  call plot3dDetectFormat(comm, filename, success, descriptor)
  if (.not. success) return

  if (descriptor%fileType /= PLOT3D_SOLUTION_FILE) then
     write(plot3dErrorMessage, '(3A)') "'", trim(filename),                                  &
          "': Not a valid PLOT3D solution file."
     success = .false.
     return
  end if

  call MPI_File_open(comm, trim(filename) // char(0), MPI_MODE_RDONLY,                       &
       MPI_INFO_NULL, mpiFileHandle, ierror)

  do i = 1, min(size(solutionVector, 2) - 1, 4)
     call MPI_File_set_view(mpiFileHandle, offset, SCALAR_TYPE_MPI,                          &
          mpiDerivedTypeScalarSubarray, "native",                                            &
          MPI_INFO_NULL, ierror)
     call MPI_File_read_all(mpiFileHandle, solutionVector(:,i), size(solutionVector, 1),     &
          SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
     if (.not. descriptor%isEndiannessNative) then
        do j = 1, size(solutionVector, 1)
           call swapEndianness(solutionVector(j,i))
        end do
     end if
     offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
  end do

  do while (i <= 4)
     offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
     i = i + 1
  end do

  i = size(solutionVector, 2)
  call MPI_File_set_view(mpiFileHandle, offset, SCALAR_TYPE_MPI,                             &
       mpiDerivedTypeScalarSubarray, "native",                                               &
       MPI_INFO_NULL, ierror)
  call MPI_File_read_all(mpiFileHandle, solutionVector(:,i), size(solutionVector, 1),        &
       SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
  if (.not. descriptor%isEndiannessNative) then
     do j = 1, size(solutionVector, 1)
        call swapEndianness(solutionVector(j,i))
     end do
  end if
  offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND)) +          &
       2 * SIZEOF_PLOT3D_OFF

  call MPI_File_close(mpiFileHandle, ierror)

  ! Successful return code.
  success = .true.

end subroutine plot3dReadSingleSolution

subroutine plot3dReadSingleFunction(comm, filename, offset, mpiDerivedTypeScalarSubarray,    &
     globalGridSize, functionVector, success)

  !> Reads the function `functionVector` corresponding to a single block at offset `offset`
  !> from the beginning of a file `filename`. The offset may be obtained by calling
  !> `PLOT3D_get_offset`.

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use PLOT3DDescriptor_type

  ! <<< Private members >>>
  use PLOT3DHelperImpl, only : swapEndianness

  ! <<< Public members >>>
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: filename
  integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
  integer, intent(in) :: mpiDerivedTypeScalarSubarray, globalGridSize(3)
  real(SCALAR_KIND), intent(out) :: functionVector(:,:)
  logical, intent(out) :: success

  ! <<< Local variables >>>
  integer :: i, j, procRank, mpiFileHandle, ierror
  type(t_PLOT3DDescriptor) :: descriptor

#ifdef DEBUG
  if (offset < 0) then
     write(plot3dErrorMessage, '(3A,2(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value of argument offset: ", offset,                                    &
          ". Expected a non-negative value!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (size(functionVector, 1) <= 0 .or. size(functionVector, 2) <= 0) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid shape of array argument functionVector: (",                             &
          size(functionVector, 1), ", ",                                                     &
          size(functionVector, 2), ")!"
     success = .false.
     return
  end if
#endif

#ifdef DEBUG
  if (any(globalGridSize <= 0)) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid value of argument globalGridSize: (/ ",                                 &
          globalGridSize(1), ", ",  globalGridSize(2), ", ",                                 &
          globalGridSize(3), " /)!"
     success = .false.
     return
  end if
#endif

  call MPI_Comm_rank(comm, procRank, ierror)
  call plot3dDetectFormat(comm, filename, success, descriptor, includeFunctionFiles = .true.)
  if (.not. success) return

  if (descriptor%fileType /= PLOT3D_FUNCTION_FILE) then
     write(plot3dErrorMessage, '(3A)') "'", trim(filename),                                  &
          "': Not a valid PLOT3D function file."
     success = .true.
     return
  end if

#ifdef DEBUG
  if (size(functionVector, 2) < descriptor%nScalars) then
     write(plot3dErrorMessage, '(3A,3(I0.0,A))') "In ", __FILE__, ":", __LINE__,             &
          ": Invalid shape of array argument functionVector: (",                             &
          size(functionVector, 1), ", ", size(functionVector, 2),                            &
          "). Expected (", size(functionVector, 1), ", ",                                    &
          descriptor%nScalars, ")!"
     success = .false.
     return
  end if
#endif

  call MPI_File_open(comm, trim(filename) // char(0), MPI_MODE_RDONLY,                       &
       MPI_INFO_NULL, mpiFileHandle, ierror)

  do i = 1, min(descriptor%nScalars, size(functionVector, 2))
     call MPI_File_set_view(mpiFileHandle, offset, SCALAR_TYPE_MPI,                          &
          mpiDerivedTypeScalarSubarray, "native", MPI_INFO_NULL, ierror)
     call MPI_File_read_all(mpiFileHandle, functionVector(:,i), size(functionVector, 1),     &
          SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
     if (.not. descriptor%isEndiannessNative) then
        do j = 1, size(functionVector, 1)
           call swapEndianness(functionVector(j,i))
        end do
     end if
     offset = offset + SIZEOF_SCALAR * product(int(globalGridSize, MPI_OFFSET_KIND))
  end do
  offset = offset + 2 * SIZEOF_PLOT3D_OFF

  call MPI_File_close(mpiFileHandle, ierror)

  ! Successful return code.
  success = .true.

end subroutine plot3dReadSingleFunction
