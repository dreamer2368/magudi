#include "config.h"

module PatchDescriptor_mod

  implicit none
  private

  type, public :: t_PatchDescriptor

     character(len = 20) :: name
     character(len = STRING_LENGTH) :: patchType
     integer :: gridIndex, normalDirection, iMin, iMax, jMin, jMax, kMin, kMax

   contains

     procedure, pass :: validate => validatePatchDescriptor
     procedure, pass :: validateInterface => validateInterfacePatchDescriptor

  end type t_PatchDescriptor

  interface

     subroutine validatePatchDescriptor(this, globalGridSizes,                               &
          simulationFlags, solverOptions, errorCode, message)

       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_PatchDescriptor

       class(t_PatchDescriptor) :: this
       integer, intent(in) :: globalGridSizes(:,:)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       integer, intent(out) :: errorCode
       character(len = STRING_LENGTH), intent(out) :: message

     end subroutine validatePatchDescriptor

  end interface

  interface

     subroutine validateInterfacePatchDescriptor(this, globalGridSizes, simulationFlags,     &
          solverOptions, interfaceIndexReordering, interfaceDescriptor, errorCode, message)

       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_PatchDescriptor

       class(t_PatchDescriptor) :: this
       integer, intent(in) :: globalGridSizes(:,:)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_PatchDescriptor) :: interfaceDescriptor
       integer, intent(in) :: interfaceIndexReordering(3)
       integer, intent(out) :: errorCode
       character(len = STRING_LENGTH), intent(out) :: message

     end subroutine validateInterfacePatchDescriptor

  end interface

end module PatchDescriptor_mod
