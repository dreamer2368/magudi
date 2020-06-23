#include "config.h"

module KolmogorovForcingPatch_mod

  use Patch_mod, only : t_Patch

  implicit none

  integer, parameter, private :: real64 = selected_real_kind(15)

  type, extends(t_Patch), public :: t_KolmogorovForcingPatch

     integer :: n
     real(SCALAR_KIND) :: amplitude

     SCALAR_TYPE, allocatable :: forcePerUnitMass(:)

   contains

     procedure, pass :: setup => setupKolmogorovForcingPatch
     procedure, pass :: cleanup => cleanupKolmogorovForcingPatch
     procedure, pass :: verifyUsage => verifyKolmogorovForcingPatchUsage
     procedure, pass :: updateRhs => addKolmogorovForcing

  end type t_KolmogorovForcingPatch

  interface

     subroutine setupKolmogorovForcingPatch(this, index, comm, patchDescriptor,           &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_KolmogorovForcingPatch

       class(t_KolmogorovForcingPatch) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupKolmogorovForcingPatch

  end interface

  interface

     subroutine cleanupKolmogorovForcingPatch(this)

       import :: t_KolmogorovForcingPatch

       class(t_KolmogorovForcingPatch) :: this

     end subroutine cleanupKolmogorovForcingPatch

  end interface

  interface

     function verifyKolmogorovForcingPatchUsage(this, patchDescriptor, gridSize,          &
          normalDirection, extent, simulationFlags,                                          &
          success, message) result(isPatchUsed)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_KolmogorovForcingPatch

       class(t_KolmogorovForcingPatch) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       logical, intent(out) :: success
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchUsed

     end function verifyKolmogorovForcingPatchUsage

  end interface

  interface

     subroutine addKolmogorovForcing(this, mode, simulationFlags,                         &
          solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_KolmogorovForcingPatch

       class(t_KolmogorovForcingPatch) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine addKolmogorovForcing

  end interface

end module KolmogorovForcingPatch_mod
