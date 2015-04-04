#include "config.h"

module State_enum

  implicit none
  public

  integer, parameter ::                                                                      &
       QOI_FORWARD_STATE        =  100,                                                      &
       QOI_ADJOINT_STATE        =  101,                                                      &
       QOI_TARGET_STATE         =  102,                                                      &
       QOI_VORTICITY_DILATATION =  103,                                                      &
       QOI_DUMMY_FUNCTION       =  104

end module State_enum

module State_mod

  use AcousticSource_mod, only : t_AcousticSource

  implicit none
  private

  integer, parameter, private :: wp = SCALAR_KIND

  interface

     function getFileType(quantityOfInterest) result(fileType)

       integer, intent(in) :: quantityOfInterest
       integer :: fileType

     end function getFileType

  end interface

  interface

     function getNumberOfScalars(quantityOfInterest, nDimensions) result(nScalars)

       integer, intent(in) :: quantityOfInterest, nDimensions
       integer :: nScalars

     end function getNumberOfScalars

  end interface

  public :: getFileType, getNumberOfScalars

  type, public :: t_State

     type(t_AcousticSource), allocatable :: acousticSources(:)

     real(wp) :: adjointForcingFactor = 1.0_wp
     SCALAR_TYPE :: plot3dAuxiliaryData(4) = 0.0_wp

     SCALAR_TYPE, dimension(:,:), allocatable :: rightHandSide, conservedVariables,          &
          specificVolume, velocity, velocityGradient, stressTensor, pressure, temperature,   &
          heatFlux, dynamicViscosity, secondCoefficientOfViscosity, thermalDiffusivity,      &
          targetState, adjointVariables

     SCALAR_TYPE, dimension(:,:), pointer :: dummyFunction => null()

   contains

     procedure, pass :: setup => setupState
     procedure, pass :: cleanup => cleanupState
     procedure, pass :: loadData => loadStateData
     procedure, pass :: saveData => saveStateData
     procedure, pass :: makeQuiescent
     procedure, pass :: update => updateState
     procedure, pass :: computeCfl => computeStateCfl
     procedure, pass :: computeTimeStepSize => computeStateTimeStepSize
     procedure, pass :: addSources

  end type t_State

  interface

     subroutine setupState(this, grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_State

       class(t_State) :: this
       class(t_Grid) :: grid

       type(t_SimulationFlags), intent(in), optional :: simulationFlags
       type(t_SolverOptions), intent(in), optional :: solverOptions

     end subroutine setupState

  end interface

  interface

     subroutine cleanupState(this)

       import :: t_State

       class(t_State) :: this

     end subroutine cleanupState

  end interface

  interface

     subroutine loadStateData(this, grid, quantityOfInterest, filename, offset, success)

       use MPI, only : MPI_OFFSET_KIND
       use Grid_mod, only : t_Grid
       
       import :: t_State

       class(t_State) :: this
       class(t_Grid) :: grid
       integer, intent(in) :: quantityOfInterest
       character(len = *), intent(in) :: filename
       integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
       logical, intent(out) :: success

     end subroutine loadStateData

  end interface

  interface

     subroutine saveStateData(this, grid, quantityOfInterest, filename, offset, success)

       use MPI, only : MPI_OFFSET_KIND
       use Grid_mod, only : t_Grid
       
       import :: t_State

       class(t_State) :: this
       class(t_Grid) :: grid
       integer, intent(in) :: quantityOfInterest
       character(len = *), intent(in) :: filename
       integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
       logical, intent(out) :: success

     end subroutine saveStateData

  end interface

  interface

     subroutine makeQuiescent(this, nDimensions, ratioOfSpecificHeats, conservedVariables)

       import :: t_State

       class(t_State) :: this
       integer, intent(in) :: nDimensions
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats

       SCALAR_TYPE, intent(out), optional :: conservedVariables(:,:)

     end subroutine makeQuiescent

  end interface

  interface

     subroutine updateState(this, grid, simulationFlags, solverOptions, conservedVariables)

       use Grid_mod, only : t_Grid      
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_State

       class(t_State) :: this
       class(t_Grid) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

       SCALAR_TYPE, intent(in), optional :: conservedVariables(:,:)

     end subroutine updateState

  end interface

  interface

     function computeStateCfl(this, grid, simulationFlags, solverOptions) result(cfl)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_State

       class(t_State) :: this
       class(t_Grid) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

       real(SCALAR_KIND) :: cfl

     end function computeStateCfl

  end interface

  interface

     function computeStateTimeStepSize(this, grid, simulationFlags,                          &
          solverOptions) result(timeStepSize)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_State

       class(t_State) :: this
       class(t_Grid) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

       real(SCALAR_KIND) :: timeStepSize

     end function computeStateTimeStepSize

  end interface

  interface

     subroutine addSources(this, mode, time, grid)

       use Grid_mod, only : t_Grid

       import :: t_State

       class(t_State) :: this
       integer, intent(in) :: mode
       real(SCALAR_KIND), intent(in) :: time
       class(t_Grid) :: grid

     end subroutine addSources

  end interface

end module State_mod
