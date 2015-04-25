#include "config.h"

module Grid_enum

  implicit none
  public

  integer, parameter ::                                                                      &
       QOI_GRID              =  0,                                                           &
       QOI_METRICS           =  1,                                                           &
       QOI_JACOBIAN          =  2,                                                           &
       QOI_TARGET_MOLLIFIER  =  3,                                                           &
       QOI_CONTROL_MOLLIFIER =  4

end module Grid_enum

module Grid_mod

  use MPI, only : MPI_COMM_NULL, MPI_DATATYPE_NULL

  use StencilOperator_mod, only : t_StencilOperator

  implicit none
  private

  type, public :: t_Grid

     type(t_StencilOperator), allocatable :: firstDerivative(:), secondDerivative(:),        &
          dissipation(:), dissipationTranspose(:), adjointFirstDerivative(:)

     integer, dimension(:), allocatable :: iblank
     SCALAR_TYPE, dimension(:,:), allocatable :: coordinates, jacobian, metrics, norm,       &
          controlMollifier, targetMollifier
#ifdef SCALAR_TYPE_IS_binary128_IEEE754
     real(SCALAR_KIND), allocatable :: mpiReduceBuffer(:)
#endif

     integer :: index, comm = MPI_COMM_NULL, nDimensions, globalSize(3), localSize(3),       &
          offset(3), nGridPoints = 0, periodicityType(3)
     integer :: mpiDerivedTypeScalarSubarray = MPI_DATATYPE_NULL,                            &
          mpiDerivedTypeIntegerSubarray = MPI_DATATYPE_NULL
     logical :: isCurvilinear
     real(SCALAR_KIND) :: periodicLength(3)

   contains

     procedure, pass :: setup => setupGrid
     procedure, pass :: cleanup => cleanupGrid
     procedure, pass :: loadData => loadGridData
     procedure, pass :: saveData => saveGridData
     procedure, pass :: setupSpatialDiscretization
     procedure, pass :: computeCoordinateDerivatives
     procedure, pass :: update => updateGrid
     generic :: computeInnerProduct => computeScalarInnerProduct, computeVectorInnerProduct
     generic :: computeGradient => computeGradientOfScalar, computeGradientOfVector
     procedure, pass :: findMinimum
     procedure, pass :: findMaximum
     procedure, pass :: isVariableWithinRange

     procedure, private, pass :: computeScalarInnerProduct
     procedure, private, pass :: computeVectorInnerProduct
     procedure, private, pass :: computeGradientOfScalar
     procedure, private, pass :: computeGradientOfVector

  end type t_Grid

  interface

     subroutine setupGrid(this, index, globalGridSize, comm, processDistribution,            &
          periodicityType, periodicLength, simulationFlags)

       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Grid

       class(t_Grid) :: this
       integer, intent(in) :: index, globalGridSize(:)

       integer, intent(in), optional :: comm, processDistribution(:), periodicityType(:)
       real(SCALAR_KIND), intent(in), optional :: periodicLength(:)

       type(t_SimulationFlags), intent(in), optional :: simulationFlags

     end subroutine setupGrid

  end interface

  interface

     subroutine cleanupGrid(this)

       import :: t_Grid

       class(t_Grid) :: this

     end subroutine cleanupGrid

  end interface

  interface

     subroutine loadGridData(this, quantityOfInterest, filename, offsetInBytes, success)

       use MPI, only : MPI_OFFSET_KIND

       import :: t_Grid

       class(t_Grid) :: this
       integer, intent(in) :: quantityOfInterest
       character(len = *), intent(in) :: filename
       integer(kind = MPI_OFFSET_KIND), intent(inout) :: offsetInBytes
       logical, intent(out) :: success

     end subroutine loadGridData

  end interface

  interface

     subroutine saveGridData(this, quantityOfInterest, filename, offsetInBytes, success)

       use MPI, only : MPI_OFFSET_KIND

       import :: t_Grid

       class(t_Grid) :: this
       integer, intent(in) :: quantityOfInterest
       character(len = *), intent(in) :: filename
       integer(kind = MPI_OFFSET_KIND), intent(inout) :: offsetInBytes
       logical, intent(out) :: success

     end subroutine saveGridData

  end interface

  interface

     subroutine setupSpatialDiscretization(this, simulationFlags, solverOptions)

       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Grid

       class(t_Grid) :: this

       type(t_SimulationFlags), intent(in), optional :: simulationFlags
       type(t_SolverOptions), intent(in), optional :: solverOptions

     end subroutine setupSpatialDiscretization

  end interface

  interface

     subroutine computeCoordinateDerivatives(this, direction, coordinateDerivatives)

       import :: t_Grid

       class(t_Grid) :: this
       integer, intent(in) :: direction
       SCALAR_TYPE, intent(out) :: coordinateDerivatives(:,:)

     end subroutine computeCoordinateDerivatives

  end interface

  interface

     subroutine updateGrid(this, hasNegativeJacobian, errorMessage)

       !> Updates the Jacobian, norm and normalized metrics. Generally called if the
       !> derivative schemes or the grid coordinates have changed.

       import :: t_Grid

       class(t_Grid) :: this

       logical, intent(out), optional :: hasNegativeJacobian
       character(len = STRING_LENGTH), intent(out), optional :: errorMessage

     end subroutine updateGrid

  end interface

  interface computeInnerProduct

     function computeScalarInnerProduct(this, f, g, weight) result(innerProduct)

       import :: t_Grid

       class(t_Grid) :: this
       SCALAR_TYPE, intent(in) :: f(:), g(:)

       SCALAR_TYPE, intent(in), optional :: weight(:)

       SCALAR_TYPE :: innerProduct

     end function computeScalarInnerProduct

     function computeVectorInnerProduct(this, f, g, weight) result(innerProduct)

       import :: t_Grid

       class(t_Grid) :: this
       SCALAR_TYPE, intent(in) :: f(:,:), g(:,:)

       SCALAR_TYPE, intent(in), optional :: weight(:)

       SCALAR_TYPE :: innerProduct

     end function computeVectorInnerProduct

  end interface computeInnerProduct

  interface computeGradient

     subroutine computeGradientOfScalar(this, f, gradF)

       import :: t_Grid

       class(t_Grid) :: this
       SCALAR_TYPE, intent(in) :: f(:)
       SCALAR_TYPE, intent(out) :: gradF(:,:)

     end subroutine computeGradientOfScalar

     subroutine computeGradientOfVector(this, f, gradF)

       import :: t_Grid

       class(t_Grid) :: this
       SCALAR_TYPE, intent(in) :: f(:,:)
       SCALAR_TYPE, intent(out) :: gradF(:,:)

     end subroutine computeGradientOfVector

  end interface computeGradient

  interface

     subroutine findMinimum(this, f, fMin, iMin, jMin, kMin)

       import :: t_Grid

       class(t_Grid) :: this
       SCALAR_TYPE, intent(in) :: f(:)
       SCALAR_TYPE, intent(out) :: fMin
       integer, intent(out), optional :: iMin, jMin, kMin

     end subroutine findMinimum

  end interface

  interface

     subroutine findMaximum(this, f, fMax, iMax, jMax, kMax)

       import :: t_Grid

       class(t_Grid) :: this
       SCALAR_TYPE, intent(in) :: f(:)
       SCALAR_TYPE, intent(out) :: fMax
       integer, intent(out), optional :: iMax, jMax, kMax

     end subroutine findMaximum

  end interface

  interface

     function isVariableWithinRange(this, f, fOutsideRange,                                  &
          iOutOfRange, jOutOfRange, kOutOfRange, minValue, maxValue)

       import :: t_Grid

       class(t_Grid) :: this
       SCALAR_TYPE, intent(in) :: f(:)
       SCALAR_TYPE, intent(out) :: fOutsideRange
       integer, intent(out) :: iOutOfRange, jOutOfRange, kOutOfRange

       real(SCALAR_KIND), intent(in), optional :: minValue, maxValue

       logical :: isVariableWithinRange

     end function isVariableWithinRange

  end interface

end module Grid_mod
