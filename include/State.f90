#include "config.h"

module State_type

  use AcousticSource_mod, only : t_AcousticSource

  implicit none
  private

  integer, parameter, private :: wp = SCALAR_KIND

  integer, parameter, public ::                                                              &
       QOI_FORWARD_STATE        =  100,                                                      &
       QOI_ADJOINT_STATE        =  101,                                                      &
       QOI_TARGET_STATE         =  102,                                                      &
       QOI_VORTICITY_DILATATION =  103,                                                      &
       QOI_MEAN_PRESSURE        =  104

  type, public :: t_State

     type(t_AcousticSource), allocatable :: acousticSources(:)

     integer :: nUnknowns = 0
     real(wp) :: adjointForcingFactor = 1.0_wp
     SCALAR_TYPE :: plot3dAuxiliaryData(4) = 0.0_wp

     SCALAR_TYPE, dimension(:,:), allocatable :: rightHandSide, conservedVariables,          &
          specificVolume, velocity, velocityGradient, stressTensor, pressure, temperature,   &
          heatFlux, dynamicViscosity, secondCoefficientOfViscosity, thermalDiffusivity,      &
          meanPressure, targetState, adjointVariables

  end type t_State

end module State_type

module State_mod

  implicit none
  public

  interface

     subroutine setupState(this, grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use State_type, only : t_State
       use SolverOptions_type, only : t_SolverOptions
       use SimulationFlags_type, only : t_SimulationFlags

       type(t_State) :: this
       class(t_Grid) :: grid

       type(t_SimulationFlags), intent(in), optional :: simulationFlags
       type(t_SolverOptions), intent(in), optional :: solverOptions

     end subroutine setupState

  end interface

  interface

     subroutine cleanupState(this)

       use State_type, only : t_State

       type(t_State) :: this

     end subroutine cleanupState

  end interface

  interface

     subroutine loadStateData(this, grid, quantityOfInterest, filename, offset, success)

       use MPI, only : MPI_OFFSET_KIND
       use Grid_mod, only : t_Grid
       use State_type, only : t_State

       type(t_State) :: this
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
       use State_type, only : t_State

       type(t_State) :: this
       class(t_Grid) :: grid
       integer, intent(in) :: quantityOfInterest
       character(len = *), intent(in) :: filename
       integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
       logical, intent(out) :: success

     end subroutine saveStateData

  end interface

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

  interface

     subroutine makeQuiescent(this, nDimensions, ratioOfSpecificHeats, conservedVariables)

       use State_type, only : t_State

       type(t_State) :: this
       integer, intent(in) :: nDimensions
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats

       SCALAR_TYPE, intent(out), optional :: conservedVariables(:,:)

     end subroutine makeQuiescent

  end interface

  interface

     subroutine updateState(this, grid, simulationFlags, solverOptions, conservedVariables)

       use Grid_mod, only : t_Grid
       use State_type, only : t_State
       use SolverOptions_type, only : t_SolverOptions
       use SimulationFlags_type, only : t_SimulationFlags

       type(t_State) :: this
       class(t_Grid) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

       SCALAR_TYPE, intent(in), optional :: conservedVariables(:,:)

     end subroutine updateState

  end interface

  interface

     function cfl(this, grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use State_type, only : t_State
       use SolverOptions_type, only : t_SolverOptions
       use SimulationFlags_type, only : t_SimulationFlags

       type(t_State) :: this
       class(t_Grid) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

       real(SCALAR_KIND) :: cfl

     end function cfl

  end interface

  interface

     function timeStepSize(this, grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use State_type, only : t_State
       use SolverOptions_type, only : t_SolverOptions
       use SimulationFlags_type, only : t_SimulationFlags

       type(t_State) :: this
       class(t_Grid) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

       real(SCALAR_KIND) :: timeStepSize

     end function timeStepSize

  end interface

  interface

     subroutine computeRhsForward(this, grid, patches, time, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use Patch_type, only : t_Patch
       use State_type, only : t_State
       use SolverOptions_type, only : t_SolverOptions
       use SimulationFlags_type, only : t_SimulationFlags

       type(t_State) :: this
       class(t_Grid) :: grid
       type(t_Patch), allocatable :: patches(:)
       real(SCALAR_KIND), intent(in) :: time
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine computeRhsForward

  end interface

  interface

     subroutine computeRhsAdjoint(this, grid, patches, time, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use Patch_type, only : t_Patch
       use State_type, only : t_State
       use SolverOptions_type, only : t_SolverOptions
       use SimulationFlags_type, only : t_SimulationFlags

       type(t_State) :: this
       class(t_Grid) :: grid
       type(t_Patch), allocatable :: patches(:)
       real(SCALAR_KIND), intent(in) :: time
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine computeRhsAdjoint

  end interface

  interface

     subroutine addPenaltiesForward(this, grid, patches, time, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use Patch_type, only : t_Patch
       use State_type, only : t_State
       use SolverOptions_type, only : t_SolverOptions
       use SimulationFlags_type, only : t_SimulationFlags

       type(t_State) :: this
       class(t_Grid) :: grid
       type(t_Patch), allocatable :: patches(:)
       real(SCALAR_KIND), intent(in) :: time
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine addPenaltiesForward

  end interface

  interface

     subroutine addPenaltiesAdjoint(this, grid, patches, time, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use Patch_type, only : t_Patch
       use State_type, only : t_State
       use SolverOptions_type, only : t_SolverOptions
       use SimulationFlags_type, only : t_SimulationFlags

       type(t_State) :: this
       class(t_Grid) :: grid
       type(t_Patch), allocatable :: patches(:)
       real(SCALAR_KIND), intent(in) :: time
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine addPenaltiesAdjoint

  end interface

  interface

     subroutine addSourcesForward(this, grid, patches, time)

       use Grid_mod, only : t_Grid
       use Patch_type, only : t_Patch
       use State_type, only : t_State

       type(t_State) :: this
       class(t_Grid) :: grid
       type(t_Patch), allocatable, intent(in) :: patches(:)
       real(SCALAR_KIND), intent(in) :: time

     end subroutine addSourcesForward

  end interface

  interface

     subroutine addSourcesAdjoint(this, grid, patches, time)

       use Grid_mod, only : t_Grid
       use Patch_type, only : t_Patch
       use State_type, only : t_State

       type(t_State) :: this
       class(t_Grid) :: grid
       type(t_Patch), allocatable, intent(in) :: patches(:)
       real(SCALAR_KIND), intent(in) :: time

     end subroutine addSourcesAdjoint

  end interface

  interface

     subroutine updatePatches(this, grid, patches, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use Patch_type, only : t_Patch
       use State_type, only : t_State
       use SolverOptions_type, only : t_SolverOptions
       use SimulationFlags_type, only : t_SimulationFlags

       type(t_State) :: this
       class(t_Grid) :: grid
       type(t_Patch), allocatable :: patches(:)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine updatePatches

  end interface

  interface

     subroutine addAdjointForcing(this, grid, patch, solverOptions)

       use Grid_mod, only : t_Grid
       use Patch_type, only : t_Patch
       use State_type, only : t_State
       use SolverOptions_type, only : t_SolverOptions

       type(t_State) :: this
       class(t_Grid) :: grid
       type(t_Patch) :: patch
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine addAdjointForcing

  end interface

end module State_mod
