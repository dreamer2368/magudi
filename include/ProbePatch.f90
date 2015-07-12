#include "config.h"

module ProbePatch_mod

  use MPI, only : MPI_OFFSET_KIND

  use Patch_mod, only : t_Patch

  implicit none
  private

  type, extends(t_Patch), public :: t_ProbePatch

     integer :: iProbeBuffer = 0
     integer(kind = MPI_OFFSET_KIND) :: probeFileOffset = int(0, MPI_OFFSET_KIND)
     character(len = STRING_LENGTH) :: probeFilename
     SCALAR_TYPE, allocatable :: probeBuffer(:,:,:)

   contains

     procedure, pass :: setup => setupProbePatch
     procedure, pass :: cleanup => cleanupProbePatch
     procedure, pass :: verifyUsage => verifyProbePatchUsage
     procedure, pass :: updateRhs => updateProbePatch
     procedure, pass :: saveData => saveSolutionOnProbe

  end type t_ProbePatch

  interface

     subroutine setupProbePatch(this, index, comm, patchDescriptor,                         &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ProbePatch

       class(t_ProbePatch) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupProbePatch

  end interface

  interface

     subroutine cleanupProbePatch(this)

       import :: t_ProbePatch

       class(t_ProbePatch) :: this

     end subroutine cleanupProbePatch

  end interface

  interface

     function verifyProbePatchUsage(this, patchDescriptor, gridSize, normalDirection,       &
          extent, simulationFlags, success, message) result(isPatchUsed)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ProbePatch

       class(t_ProbePatch) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       logical, intent(out) :: success
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchUsed

     end function verifyProbePatchUsage

  end interface

  interface

     subroutine updateProbePatch(this, mode, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ProbePatch

       class(t_ProbePatch) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine updateProbePatch

  end interface

  interface

     subroutine saveSolutionOnProbe(this)

       import :: t_ProbePatch

       class(t_ProbePatch) :: this

     end subroutine saveSolutionOnProbe

  end interface

end module ProbePatch_mod
