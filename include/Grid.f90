#include "config.h"

module Grid_type

  use MPI, only : MPI_COMM_NULL, MPI_DATATYPE_NULL

  use StencilOperator_type

  implicit none
  private

  integer, parameter, public ::                                                              &
       NONE    = 0,                                                                          &
       PLANE   = 1,                                                                          &
       OVERLAP = 2

  integer, parameter, public ::                                                              &
       QOI_GRID              =  0,                                                           &
       QOI_METRICS           =  1,                                                           &
       QOI_JACOBIAN          =  2,                                                           &
       QOI_TARGET_MOLLIFIER  =  3,                                                           &
       QOI_CONTROL_MOLLIFIER =  4

  type, public :: t_Grid

     type(t_StencilOperator), allocatable :: firstDerivative(:),                             &
          secondDerivative(:), dissipation(:), adjointFirstDerivative(:)

     integer, dimension(:), allocatable :: IBLANK
     SCALAR_TYPE, dimension(:,:), allocatable :: coordinates, jacobian, metrics, norm,       &
          controlMollifier, targetMollifier
#ifdef SCALAR_TYPE_IS_binary128_IEEE754
     real(SCALAR_KIND), allocatable :: mpiReduceBuffer(:)
#endif

     integer :: index, comm = MPI_COMM_NULL, globalSize(3), periodicityType(3),              &
          localSize(3), offset(3), nGridPoints
     integer :: mpiDerivedTypeScalarSubarray = MPI_DATATYPE_NULL,                            &
          mpiDerivedTypeIntegerSubarray = MPI_DATATYPE_NULL
     logical :: isCurvilinear
     real(SCALAR_KIND) :: periodicLength(3)

  end type t_Grid

end module Grid_type

module Grid_mod

  implicit none
  public

  interface

     subroutine setupGrid(this, index, globalGridSize, comm, processDistribution,            &
          periodicityType, periodicLength, simulationFlags)

       use Grid_type
       use SimulationFlags_type

       type(t_Grid) :: this
       integer, intent(in) :: index, globalGridSize(:)

       integer, intent(in), optional :: comm, processDistribution(:), periodicityType(:)
       real(SCALAR_KIND), intent(in), optional :: periodicLength(:)

       type(t_SimulationFlags), intent(in), optional :: simulationFlags

     end subroutine setupGrid

  end interface

  interface

     subroutine cleanupGrid(this)

       use Grid_type

       type(t_Grid) :: this

     end subroutine cleanupGrid

  end interface

  interface

     subroutine loadGridData(this, quantityOfInterest, filename, offsetInBytes, success)

       use MPI, only : MPI_OFFSET_KIND
       use Grid_type

       type(t_Grid) :: this
       integer, intent(in) :: quantityOfInterest
       character(len = *), intent(in) :: filename
       integer(kind = MPI_OFFSET_KIND), intent(inout) :: offsetInBytes
       logical, intent(out) :: success

     end subroutine loadGridData

  end interface

  interface

     subroutine saveGridData(this, quantityOfInterest, filename, offsetInBytes, success)

       use MPI, only : MPI_OFFSET_KIND
       use Grid_type

       type(t_Grid) :: this
       integer, intent(in) :: quantityOfInterest
       character(len = *), intent(in) :: filename
       integer(kind = MPI_OFFSET_KIND), intent(inout) :: offsetInBytes
       logical, intent(out) :: success

     end subroutine saveGridData

  end interface

  interface

     subroutine setupSpatialDiscretization(this, success, errorMessage)

       use Grid_type

       type(t_Grid) :: this

       logical, intent(out), optional :: success
       character(len = STRING_LENGTH), intent(out), optional :: errorMessage

     end subroutine setupSpatialDiscretization

  end interface

  interface

     subroutine updateGrid(this, hasNegativeJacobian, errorMessage)

       !> Updates the Jacobian, norm and normalized metrics. Generally called if the
       !> derivative schemes or the grid coordinates have changed.

       use Grid_type

       type(t_Grid) :: this

       logical, intent(out), optional :: hasNegativeJacobian
       character(len = STRING_LENGTH), intent(out), optional :: errorMessage

     end subroutine updateGrid

  end interface

  interface

     function computeInnerProduct(this, f, g, weight) result(innerProduct)

       use Grid_type

       type(t_Grid) :: this
       SCALAR_TYPE, intent(in) :: f(:,:), g(:,:)

       SCALAR_TYPE, intent(in), optional :: weight(:)

       SCALAR_TYPE :: innerProduct

     end function computeInnerProduct

  end interface

  interface computeGradient

     subroutine computeGradientOfScalar_(this, f, gradF)

       use Grid_type

       type(t_Grid) :: this
       SCALAR_TYPE, intent(in) :: f(:)
       SCALAR_TYPE, intent(out) :: gradF(:,:)

     end subroutine computeGradientOfScalar_

     subroutine computeGradientOfVector_(this, f, gradF)

       use Grid_type

       type(t_Grid) :: this
       SCALAR_TYPE, intent(in) :: f(:,:)
       SCALAR_TYPE, intent(out) :: gradF(:,:)

     end subroutine computeGradientOfVector_

  end interface computeGradient

  interface

     subroutine findMinimum(this, f, fMin, iMin, jMin, kMin)

       use Grid_type

       type(t_Grid) :: this
       SCALAR_TYPE, intent(in) :: f(:)
       SCALAR_TYPE, intent(out) :: fMin
       integer, intent(out) :: iMin, jMin, kMin

     end subroutine findMinimum

  end interface

  interface

     subroutine findMaximum(this, f, fMax, iMax, jMax, kMax)

       use Grid_type

       type(t_Grid) :: this
       SCALAR_TYPE, intent(in) :: f(:)
       SCALAR_TYPE, intent(out) :: fMax
       integer, intent(out) :: iMax, jMax, kMax

     end subroutine findMaximum

  end interface

  interface

     function isVariableWithinRange(this, f, fOutsideRange,                                  &
          iOutOfRange, jOutOfRange, kOutOfRange, minValue, maxValue)

       use Grid_type

       type(t_Grid) :: this
       SCALAR_TYPE, intent(in) :: f(:)
       SCALAR_TYPE, intent(out) :: fOutsideRange
       integer, intent(out) :: iOutOfRange, jOutOfRange, kOutOfRange

       real(SCALAR_KIND), intent(in), optional :: minValue, maxValue

       logical :: isVariableWithinRange

     end function isVariableWithinRange

  end interface

  interface

     subroutine computeNormalizedCurveLengths(this, direction, indexAtUnitCurveLength,       &
          normalizedCurveLengths, reverseDirection, coordinateDerivatives)

       use Grid_type

       type(t_Grid), intent(in) :: this
       integer, intent(in) :: direction, indexAtUnitCurveLength
       real(SCALAR_KIND), intent(out) :: normalizedCurveLengths(:)

       logical, intent(in), optional :: reverseDirection
       SCALAR_TYPE, intent(in), optional :: coordinateDerivatives(:,:)

     end subroutine computeNormalizedCurveLengths

  end interface

  interface

     subroutine computeSpongeStrengths(this, patches)

       use Grid_type
       use Patch_type

       type(t_Grid) :: this
       type(t_Patch), allocatable :: patches(:)

     end subroutine computeSpongeStrengths

  end interface

end module Grid_mod
