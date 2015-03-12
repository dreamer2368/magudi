#include "config.h"

module PLOT3DDescriptor_type

  implicit none
  private

  integer, parameter, public ::                                                              &
       PLOT3D_GRID_FILE     = 0,                                                             &
       PLOT3D_SOLUTION_FILE = 1,                                                             &
       PLOT3D_FUNCTION_FILE = 2

  type, public :: t_PLOT3DDescriptor

     integer :: fileType, nGrids, nDimensions, nScalars
     logical :: hasIblank = .true., isEndiannessNative = .true.

  end type t_PLOT3DDescriptor

end module PLOT3DDescriptor_type

module PLOT3DHelper

  implicit none
  public

  interface

     subroutine plot3dDetectFormat(comm, filename, success,                                  &
          descriptor, globalGridSizes, includeFunctionFiles)

       !> Detects the format of a PLOT3D file `filename` and reads the grid sizes. This
       !> subroutine must be called collectively from all processes in the MPI communicator
       !> `comm`.

       use PLOT3DDescriptor_type

       integer, intent(in) :: comm
       character(len = *), intent(in) :: filename
       logical, intent(out) :: success

       type(t_PLOT3DDescriptor), intent(out), optional :: descriptor
       integer, allocatable, intent(out), optional :: globalGridSizes(:,:)
       logical, intent(in), optional :: includeFunctionFiles

     end subroutine plot3dDetectFormat

  end interface

  interface plot3dGetOffset

     function plot3dGetOffsetFromGridSizes_(fileType, globalGridSizes, gridIndex, hasIblank, &
          success, nScalars) result(offset)

       !> Computes the offset, in bytes, to the beginning of the data corresponding to block
       !> `gridIndex` in a 3D multi-block whole-format PLOT3D file. The offset counts past the
       !> leading record size and points directly to the beginning of the actual data.

       use MPI, only : MPI_OFFSET_KIND

       integer, intent(in) :: fileType, globalGridSizes(:,:), gridIndex
       logical, intent(in) :: hasIblank
       logical, intent(out) :: success

       integer, intent(in), optional :: nScalars

       integer(kind = MPI_OFFSET_KIND) :: offset

     end function plot3dGetOffsetFromGridSizes_

     function plot3dGetOffsetFromFile_(comm, filename, gridIndex, success) result(offset)

       !> Computes the offset, in bytes, to the beginning of the data corresponding to block
       !> `gridIndex` in a 3D multi-block whole-format PLOT3D file. The offset counts past the
       !> leading record size and points directly to the beginning of the actual data.

       use MPI, only : MPI_OFFSET_KIND

       integer, intent(in) :: comm
       character(len = *), intent(in) :: filename
       integer, intent(in) :: gridIndex
       logical, intent(out) :: success

       integer(kind = MPI_OFFSET_KIND) :: offset

     end function plot3dGetOffsetFromFile_

  end interface plot3dGetOffset

  interface

     subroutine plot3dWriteSkeleton(comm, filename, fileType,                                &
          globalGridSizes, success, nScalars)

       !> Writes the skeleton of a PLOT3D file to `filename`, including the header, and the
       !> leading and trailing record sizes. This subroutine must be called collectively from
       !> all processes in the MPI communicator `comm`.

       integer, intent(in) :: comm
       character(len = *), intent(in) :: filename
       integer, intent(in) :: fileType, globalGridSizes(:,:)
       logical, intent(out) :: success

       integer, intent(in), optional :: nScalars

     end subroutine plot3dWriteSkeleton

  end interface

  interface

     subroutine plot3dWriteSingleGrid(comm, filename, offset, mpiDerivedTypeScalarSubarray,  &
          mpiDerivedTypeIntegerSubarray, globalGridSize, coordinates, iblank, success)

       !> Writes the coordinates and IBLANK values corresponding to a single block at offset
       !> `offset` from the beginning of file `filename`.

       use MPI, only : MPI_OFFSET_KIND

       integer, intent(in) :: comm
       character(len = *), intent(in) :: filename
       integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
       integer, intent(in) :: mpiDerivedTypeScalarSubarray, mpiDerivedTypeIntegerSubarray,   &
            globalGridSize(3)
       SCALAR_TYPE, intent(in) :: coordinates(:,:)
       integer, intent(in) :: iblank(:)
       logical, intent(out) :: success

     end subroutine plot3dWriteSingleGrid

  end interface

  interface

     subroutine plot3dWriteSingleAuxiliarySolutionData(comm,                                 &
          filename, offset, auxiliarySolutionData, success)

       !> Writes the auxiliary solution data (consisting of 4 SCALAR_TYPE values)
       !> `auxiliarySolutionData` corresponding to a single block at offset `offset` from the
       !> beginning of a file `filename`.

       use MPI, only : MPI_OFFSET_KIND

       integer, intent(in) :: comm
       character(len = *), intent(in) :: filename
       integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
       SCALAR_TYPE, intent(in) :: auxiliarySolutionData(4)
       logical, intent(out) :: success

     end subroutine plot3dWriteSingleAuxiliarySolutionData

  end interface

  interface

     subroutine plot3dWriteSingleSolution(comm, filename, offset,                            &
          mpiDerivedTypeScalarSubarray, globalGridSize,                                      &
          solutionVector, success)

       !> Writes the solution `solutionVector` corresponding to a single block at offset
       !> `offset` from the beginning of a file `filename`.

       use MPI, only : MPI_OFFSET_KIND

       integer, intent(in) :: comm
       character(len = *), intent(in) :: filename
       integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
       integer, intent(in) :: mpiDerivedTypeScalarSubarray, globalGridSize(3)
       SCALAR_TYPE, intent(in) :: solutionVector(:,:)
       logical, intent(out) :: success

     end subroutine plot3dWriteSingleSolution

  end interface

  interface

     subroutine plot3dWriteSingleFunction(comm, filename, offset,                            &
          mpiDerivedTypeScalarSubarray, globalGridSize,                                      &
          functionVector, success)

       !> Writes the function `functionVector` corresponding to a single block at offset
       !> `offset` from the beginning of a file `filename`.

       use MPI, only : MPI_OFFSET_KIND

       integer, intent(in) :: comm
       character(len = *), intent(in) :: filename
       integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
       integer, intent(in) :: mpiDerivedTypeScalarSubarray, globalGridSize(3)
       SCALAR_TYPE, intent(in) :: functionVector(:,:)
       logical, intent(out) :: success

     end subroutine plot3dWriteSingleFunction

  end interface

  interface

     subroutine plot3dReadSingleGrid(comm, filename, offset, mpiDerivedTypeScalarSubarray,   &
          mpiDerivedTypeIntegerSubarray, globalGridSize, coordinates, iblank, success)

       !> Reads the coordinates and IBLANK values corresponding to a single block at offset
       !> `offset` from the beginning of a file `filename`. The offset may be obtained by
       !> calling `plot3dGetOffset`. If IBLANK values are not present in `filename`,
       !> `iblank` is set to 1.

       use MPI, only : MPI_OFFSET_KIND

       integer, intent(in) :: comm
       character(len = *), intent(in) :: filename
       integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
       integer, intent(in) :: mpiDerivedTypeScalarSubarray,                                  &
            mpiDerivedTypeIntegerSubarray, globalGridSize(3)
       SCALAR_TYPE, intent(out) :: coordinates(:,:)
       integer, intent(out) :: iblank(:)
       logical, intent(out) :: success

     end subroutine plot3dReadSingleGrid

  end interface

  interface

     subroutine plot3dReadSingleAuxiliarySolutionData(comm, filename,                        &
          offset, auxiliarySolutionData, success)

       !> Reads the auxiliary solution data (consisting of 4 SCALAR_TYPE values)
       !> `auxiliarySolutionData` corresponding to a single block at offset `offset` from the
       !> beginning of a file `filename`. The offset may be obtained by calling
       !> `plot3dGetOffset`.

       use MPI, only : MPI_OFFSET_KIND

       integer, intent(in) :: comm
       character(len = *), intent(in) :: filename
       integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
       SCALAR_TYPE, intent(out) :: auxiliarySolutionData(4)
       logical, intent(out) :: success

     end subroutine plot3dReadSingleAuxiliarySolutionData

  end interface

  interface

     subroutine plot3dReadSingleSolution(comm, filename, offset,                             &
          mpiDerivedTypeScalarSubarray, globalGridSize,                                      &
          solutionVector, success)

       !> Reads the solution `solutionVector` corresponding to a single block at offset
       !> `offset` from the beginning of a file `filename`. The offset may be obtained by
       !> calling `PLOT3D_get_offset`.

       use MPI, only : MPI_OFFSET_KIND

       integer, intent(in) :: comm
       character(len = *), intent(in) :: filename
       integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
       integer, intent(in) :: mpiDerivedTypeScalarSubarray, globalGridSize(3)
       SCALAR_TYPE, intent(out) :: solutionVector(:,:)
       logical, intent(out) :: success

     end subroutine plot3dReadSingleSolution

  end interface

  interface

     subroutine plot3dReadSingleFunction(comm, filename, offset,                             &
          mpiDerivedTypeScalarSubarray, globalGridSize,                                      &
          functionVector, success)

       use MPI, only : MPI_OFFSET_KIND

       integer, intent(in) :: comm
       character(len = *), intent(in) :: filename
       integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
       integer, intent(in) :: mpiDerivedTypeScalarSubarray, globalGridSize(3)
       SCALAR_TYPE, intent(out) :: functionVector(:,:)
       logical, intent(out) :: success

     end subroutine plot3dReadSingleFunction

  end interface

  character(len = STRING_LENGTH), public :: errorMessage

  private :: plot3dGetOffsetFromGridSizes_, plot3dGetOffsetFromFile_

end module PLOT3DHelper
