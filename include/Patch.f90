#include "config.h"

module Patch_mod

  use MPI, only : MPI_COMM_NULL, MPI_DATATYPE_NULL

  implicit none
  private

  type, abstract, public :: t_Patch

     integer :: index, normalDirection, gridIndex, extent(6), nDimensions, patchSize(3),     &
          offset(3), gridLocalSize(3), gridOffset(3), nPatchPoints, comm = MPI_COMM_NULL
     logical :: isCurvilinear

   contains

     procedure, non_overridable, pass :: setupBase => setupPatch
     procedure, non_overridable, pass :: cleanupBase => cleanupPatch

     procedure(setup), pass, deferred :: setup
     procedure(cleanup), pass, deferred :: cleanup
     procedure(update), pass, deferred :: update
     procedure(verifyUsage), pass, deferred :: verifyUsage
     generic :: collect => collectScalarAtPatch,                                             &
          collectVectorAtPatch, collectTensorAtPatch
     procedure(updateRhs), pass, deferred :: updateRhs

     procedure, private, pass :: collectScalarAtPatch
     procedure, private, pass :: collectVectorAtPatch
     procedure, private, pass :: collectTensorAtPatch

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

end module Patch_mod
