#include "config.h"

module PatchDescriptor_mod

  implicit none
  private

  type, public :: t_PatchDescriptor

     character(len = 20) :: name, patchType
     integer :: gridIndex, normalDirection, iMin, iMax, jMin, jMax, kMin, kMax

   contains

     procedure, pass :: validate => validatePatchDescriptor

  end type t_PatchDescriptor

  interface

     subroutine validatePatchDescriptor(this, globalGridSizes,                               &
          simulationFlags, errorCode, message)

       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_PatchDescriptor

       class(t_PatchDescriptor) :: this
       integer, intent(in) :: globalGridSizes(:,:)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       integer, intent(out) :: errorCode
       character(len = STRING_LENGTH), intent(out) :: message

     end subroutine validatePatchDescriptor

  end interface

end module PatchDescriptor_mod
