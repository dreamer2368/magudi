#include "config.h"

module PatchDescriptor_type

  implicit none
  private

  integer, parameter, public ::                                                              &
       SPONGE                = 1,                                                            &
       ACTUATOR              = 2,                                                            &
       CONTROL_TARGET        = 3,                                                            &
       SOLENOIDAL_EXCITATION = 4,                                                            &
       SAT_FAR_FIELD         = 5,                                                            &
       SAT_SLIP_WALL         = 6,                                                            &
       SAT_ISOTHERMAL_WALL   = 7,                                                            &
       SAT_ADIABATIC_WALL    = 8,                                                            &
       SAT_BLOCK_INTERFACE   = 9

  type, public :: t_PatchDescriptor

     character(len = 20) :: name
     integer :: patchType, normalDirection, gridIndex, iMin, iMax, jMin, jMax, kMin, kMax

  end type t_PatchDescriptor

end module PatchDescriptor_type

module PatchDescriptor_mod

  implicit none
  public

  interface

     subroutine parsePatchType(patchTypeString, patchType)

       !> If `patchTypeString` matches one of the enumerated patch types defined in this
       !> module, the corresponding numeric value is returned in `patchType`. Otherwise,
       !> `patchType` is set to -1.

       character(len = *), intent(in) :: patchTypeString
       integer, intent(out) :: patchType

     end subroutine parsePatchType

  end interface

  interface

     subroutine validatePatchDescriptor(this, globalGridSizes,                               &
          simulationFlags, errorCode, message)

       use PatchDescriptor_type
       use SimulationFlags_mod, only : t_SimulationFlags

       type(t_PatchDescriptor) :: this
       integer, intent(in) :: globalGridSizes(:,:)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       integer, intent(out) :: errorCode
       character(len = STRING_LENGTH), intent(out) :: message

     end subroutine validatePatchDescriptor

  end interface

  interface

     subroutine validatePatchesConnectivity(patchDescriptors, errorCode, message)

       use PatchDescriptor_type

       type(t_PatchDescriptor), intent(in) :: patchDescriptors(:)
       integer, intent(out) :: errorCode
       character(len = STRING_LENGTH), intent(out) :: message

     end subroutine validatePatchesConnectivity

  end interface

end module PatchDescriptor_mod
