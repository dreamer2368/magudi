#include "config.h"

module Patch_mod

  use MPI, only : MPI_COMM_NULL, MPI_DATATYPE_NULL

  implicit none
  private

  type, abstract, public :: t_Patch

     integer :: index, comm = MPI_COMM_NULL, nDimensions, globalSize(3), localSize(3),       &
          offset(3), nPatchPoints = 0, gridIndex, normalDirection, extent(6),                &
          gridLocalSize(3), gridOffset(3)
     integer, allocatable :: mpiAllScalarSubarrayTypes(:), mpiAllIntegerSubarrayTypes(:)
     logical :: isCurvilinear, penaltyInPhysicalCoordinates = .false.
#ifdef SCALAR_TYPE_IS_binary128_IEEE754
     real(SCALAR_KIND), allocatable :: mpiReduceBuffer(:)
#endif

   contains

     procedure, non_overridable, pass :: setupBase => setupPatch
     procedure, non_overridable, pass :: cleanupBase => cleanupPatch
     generic :: collect => collectScalarAtPatch,                                             &
          collectVectorAtPatch, collectTensorAtPatch
     generic :: gatherData => gatherScalarOnPatch, gatherVectorOnPatch, gatherTensorOnPatch
     generic :: scatterData => scatterScalarOnPatch, scatterVectorOnPatch,                   &
          scatterTensorOnPatch

     procedure(setup), pass, deferred :: setup
     procedure(cleanup), pass, deferred :: cleanup
     procedure(update), pass, deferred :: update
     procedure(verifyUsage), pass, deferred :: verifyUsage
     procedure(updateRhs), pass, deferred :: updateRhs

     procedure, private, pass :: collectScalarAtPatch
     procedure, private, pass :: collectVectorAtPatch
     procedure, private, pass :: collectTensorAtPatch

     procedure, private, pass :: gatherScalarOnPatch
     procedure, private, pass :: gatherVectorOnPatch
     procedure, private, pass :: gatherTensorOnPatch

     procedure, private, pass :: scatterScalarOnPatch
     procedure, private, pass :: scatterVectorOnPatch
     procedure, private, pass :: scatterTensorOnPatch

  end type t_Patch

  abstract interface

     subroutine setup(this, index, comm, patchDescriptor,                                    &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Patch

       class(t_Patch) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setup

  end interface

  abstract interface

     subroutine cleanup(this)

       import :: t_Patch

       class(t_Patch) :: this

     end subroutine cleanup

  end interface

  abstract interface

     subroutine update(this, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Patch

       class(t_Patch) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state

     end subroutine update

  end interface

  abstract interface

     function verifyUsage(this, patchDescriptor, gridSize, normalDirection, extent,          &
          simulationFlags, success, message) result(isPatchUsed)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Patch

       class(t_Patch) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       logical, intent(out) :: success
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchUsed

     end function verifyUsage

  end interface

  abstract interface

     subroutine updateRhs(this, mode, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Patch

       class(t_Patch) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine updateRhs

  end interface

  interface

     subroutine setupPatch(this, index, comm, patchDescriptor,                               &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Patch

       class(t_Patch) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupPatch

  end interface

  interface

     subroutine cleanupPatch(this)

       import :: t_Patch

       class(t_Patch) :: this

     end subroutine cleanupPatch

  end interface

  interface collectAtPatch

     subroutine collectScalarAtPatch(this, gridArray, patchArray)

       import :: t_Patch

       class(t_Patch) :: this

       SCALAR_TYPE, intent(in) :: gridArray(:)
       SCALAR_TYPE, intent(out) :: patchArray(:)

     end subroutine collectScalarAtPatch

     subroutine collectVectorAtPatch(this, gridArray, patchArray)

       import :: t_Patch

       class(t_Patch) :: this

       SCALAR_TYPE, intent(in) :: gridArray(:,:)
       SCALAR_TYPE, intent(out) :: patchArray(:,:)

     end subroutine collectVectorAtPatch

     subroutine collectTensorAtPatch(this, gridArray, patchArray)

       import :: t_Patch

       class(t_Patch) :: this

       SCALAR_TYPE, intent(in) :: gridArray(:,:,:)
       SCALAR_TYPE, intent(out) :: patchArray(:,:,:)

     end subroutine collectTensorAtPatch

  end interface collectAtPatch

  interface gatherDataOnPatch

     subroutine gatherScalarOnPatch(this, patchLocalArray, patchGlobalArray)

       import :: t_Patch

       class(t_Patch) :: this
       SCALAR_TYPE, intent(in) :: patchLocalArray(:)
       SCALAR_TYPE, allocatable :: patchGlobalArray(:)

     end subroutine gatherScalarOnPatch

     subroutine gatherVectorOnPatch(this, patchLocalArray, patchGlobalArray)

       import :: t_Patch

       class(t_Patch) :: this
       SCALAR_TYPE, intent(in) :: patchLocalArray(:,:)
       SCALAR_TYPE, allocatable :: patchGlobalArray(:,:)

     end subroutine gatherVectorOnPatch

     subroutine gatherTensorOnPatch(this, patchLocalArray, patchGlobalArray)

       import :: t_Patch

       class(t_Patch) :: this
       SCALAR_TYPE, intent(in) :: patchLocalArray(:,:,:)
       SCALAR_TYPE, allocatable :: patchGlobalArray(:,:,:)

     end subroutine gatherTensorOnPatch

  end interface gatherDataOnPatch

  interface scatterDataOnPatch

     subroutine scatterScalarOnPatch(this, patchGlobalArray, patchLocalArray)

       import :: t_Patch

       class(t_Patch) :: this
       SCALAR_TYPE, intent(in), allocatable :: patchGlobalArray(:)
       SCALAR_TYPE, intent(out) :: patchLocalArray(:)

     end subroutine scatterScalarOnPatch

     subroutine scatterVectorOnPatch(this, patchGlobalArray, patchLocalArray)

       import :: t_Patch

       class(t_Patch) :: this
       SCALAR_TYPE, intent(in), allocatable :: patchGlobalArray(:,:)
       SCALAR_TYPE, intent(out) :: patchLocalArray(:,:)

     end subroutine scatterVectorOnPatch

     subroutine scatterTensorOnPatch(this, patchGlobalArray, patchLocalArray)

       import :: t_Patch

       class(t_Patch) :: this
       SCALAR_TYPE, intent(in), allocatable :: patchGlobalArray(:,:,:)
       SCALAR_TYPE, intent(out) :: patchLocalArray(:,:,:)

     end subroutine scatterTensorOnPatch

  end interface scatterDataOnPatch

end module Patch_mod
