#include "config.h"

module PatchDescriptor_type

  implicit none
  private

  integer, parameter, public ::                                                              &
       SPONGE                        = 1,                                                    &
       ACTUATOR                      = 2,                                                    &
       SOLENOIDAL_EXCITATION_SUPPORT = 3,                                                    &
       SAT_FAR_FIELD                 = 4,                                                    &
       SAT_SLIP_WALL                 = 5,                                                    &
       SAT_ISOTHERMAL_WALL           = 6,                                                    &
       SAT_BLOCK_INTERFACE           = 7

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

       character(len = *), intent(in) :: patchTypeString
       integer, intent(out) :: patchType

     end subroutine parsePatchType

  end interface

  interface

     subroutine validatePatchDescriptor(this, globalGridSizes,                               &
          simulationFlags, errorCode, message)

       use PatchDescriptor_type
       use SimulationFlags_type

       type(t_PatchDescriptor) :: this
       integer, intent(in) :: globalGridSizes(:,:)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       integer, intent(out) :: errorCode
       character(len = STRING_LENGTH), intent(out), optional :: message

     end subroutine validatePatchDescriptor

  end interface

end module PatchDescriptor_mod
