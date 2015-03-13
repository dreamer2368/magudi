#include "config.h"

module Patch_type

  use MPI, only : MPI_COMM_NULL, MPI_DATATYPE_NULL

  use PatchDescriptor_type

  implicit none
  private

  type, public :: t_Patch

     integer :: index, normalDirection, gridIndex, patchType, extent(6), globalPatchSize(3), &
          patchSize(3), offset(3), gridLocalSize(3), gridOffset(3), nPatchPoints,            &
          comm = MPI_COMM_NULL
     logical :: isCurvilinear

     ! Far-field patch variables.
     SCALAR_TYPE, dimension(:,:,:), allocatable :: viscousFluxes, targetViscousFluxes

     ! Variables common to far-field and wall patches.
     real(SCALAR_KIND) :: inviscidPenaltyAmount, viscousPenaltyAmount
     SCALAR_TYPE, allocatable :: metrics(:,:)

     ! Sponge patch variables.
     real(SCALAR_KIND) :: spongeAmount
     integer :: spongeExponent
     real(SCALAR_KIND), allocatable :: spongeStrength(:)

     ! Actuator patch variables.
     SCALAR_TYPE, allocatable :: gradient(:,:)

  end type t_Patch

end module Patch_type

module Patch_mod

  implicit none
  public

  interface

     subroutine setupPatch(this, index, nDimensions, patchDescriptor, comm,                  &
          gridOffset, gridLocalSize, simulationFlags)

       use Patch_type
       use PatchDescriptor_type
       use SimulationFlags_type

       type(t_Patch) :: this
       integer, intent(in) :: index, nDimensions
       type(t_PatchDescriptor) :: patchDescriptor
       integer, intent(in) :: comm, gridOffset(3), gridLocalSize(3)
       type(t_SimulationFlags), intent(in) :: simulationFlags

     end subroutine setupPatch

  end interface

  interface

     subroutine cleanupPatch(this)

       use Patch_type

       type(t_Patch) :: this

     end subroutine cleanupPatch

  end interface

  interface collectAtPatch

     subroutine collectScalarAtPatch_(this, gridArray, patchArray)

       use Patch_type

       type(t_Patch) :: this

       SCALAR_TYPE, intent(in) :: gridArray(:)
       SCALAR_TYPE, intent(out) :: patchArray(:)

     end subroutine collectScalarAtPatch_

     subroutine collectVectorAtPatch_(this, gridArray, patchArray)

       use Patch_type

       type(t_Patch) :: this

       SCALAR_TYPE, intent(in) :: gridArray(:,:)
       SCALAR_TYPE, intent(out) :: patchArray(:,:)

     end subroutine collectVectorAtPatch_

     subroutine collectTensorAtPatch_(this, gridArray, patchArray)

       use Patch_type

       type(t_Patch) :: this

       SCALAR_TYPE, intent(in) :: gridArray(:,:,:)
       SCALAR_TYPE, intent(out) :: patchArray(:,:,:)

     end subroutine collectTensorAtPatch_

  end interface collectAtPatch

  interface

     subroutine addDamping(this, mode, rightHandSide, iblank,                                &
          solvedVariables, targetVariables)

       use Patch_type

       type(t_Patch) :: this
       integer, intent(in) :: mode
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)
       integer, intent(in) :: iblank(:)
       SCALAR_TYPE, intent(in) :: solvedVariables(:,:)

       SCALAR_TYPE, intent(in), optional :: targetVariables(:,:)

     end subroutine addDamping

  end interface

  interface

     subroutine addFarFieldPenalty(this, mode, rightHandSide, iblank, nDimensions,           &
          ratioOfSpecificHeats, conservedVariables, targetState)

       use Patch_type

       type(t_Patch) :: this
       integer, intent(in) :: mode
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)
       integer, intent(in) :: iblank(:), nDimensions
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(in) :: conservedVariables(:,:), targetState(:,:)

     end subroutine addFarFieldPenalty

  end interface

  interface

     subroutine addWallPenalty(this, mode, rightHandSide, iblank, nDimensions,               &
          ratioOfSpecificHeats, conservedVariables)

       use Patch_type

       type(t_Patch) :: this
       integer, intent(in) :: mode
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)
       integer, intent(in) :: iblank(:), nDimensions
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(in) :: conservedVariables(:,:)

     end subroutine addWallPenalty

  end interface

end module Patch_mod
