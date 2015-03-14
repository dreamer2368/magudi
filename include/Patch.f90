#include "config.h"

module Patch_type

  use MPI, only : MPI_COMM_NULL, MPI_DATATYPE_NULL

  use PatchDescriptor_type
  use SolenoidalExcitation_type

  implicit none
  private

  type, public :: t_Patch

     integer :: index, normalDirection, gridIndex, patchType, extent(6), patchSize(3),       &
          offset(3), gridLocalSize(3), gridOffset(3), nPatchPoints, comm = MPI_COMM_NULL
     logical :: isCurvilinear

     ! Common to far-field and walls.
     real(SCALAR_KIND) :: inviscidPenaltyAmount, viscousPenaltyAmount
     SCALAR_TYPE, allocatable :: metrics(:,:)

     ! Far-field variables.
     SCALAR_TYPE, dimension(:,:,:), allocatable :: viscousFluxes, targetViscousFluxes

     ! Sponge variables.
     real(SCALAR_KIND) :: spongeAmount
     integer :: spongeExponent
     real(SCALAR_KIND), allocatable :: spongeStrength(:)

     ! Actuator variables.
     SCALAR_TYPE, allocatable :: gradient(:,:)

     ! Solenoidal excitation variables.
     type(t_SolenoidalExcitation) :: solenoidalExcitation
     real(SCALAR_KIND), allocatable :: solenoidalExcitationStrength(:)

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
          ratioOfSpecificHeats, conservedVariables, targetState, adjointVariables)

       use Patch_type

       type(t_Patch) :: this
       integer, intent(in) :: mode
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)
       integer, intent(in) :: iblank(:), nDimensions
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(in) :: conservedVariables(:,:), targetState(:,:)
       SCALAR_TYPE, intent(in), optional :: adjointVariables(:,:)

     end subroutine addFarFieldPenalty

  end interface

  interface

     subroutine addWallPenalty(this, mode, rightHandSide, iblank, nDimensions,               &
          ratioOfSpecificHeats, conservedVariables, adjointVariables)

       use Patch_type

       type(t_Patch) :: this
       integer, intent(in) :: mode
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)
       integer, intent(in) :: iblank(:), nDimensions
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(in) :: conservedVariables(:,:)
       SCALAR_TYPE, intent(in), optional :: adjointVariables(:,:)

     end subroutine addWallPenalty

  end interface

  interface

     subroutine updateSolenoidalExcitationStrength(this, coordinates, iblank)

       use Patch_type

       type(t_Patch) :: this
       SCALAR_TYPE, intent(in) :: coordinates(:,:)
       integer, intent(in) :: iblank(:)

     end subroutine updateSolenoidalExcitationStrength

  end interface

  interface

     subroutine addSolenoidalExcitation(this, coordinates, iblank, time, rightHandSide)

       use Patch_type

       type(t_Patch) :: this
       SCALAR_TYPE, intent(in) :: coordinates(:,:)
       integer, intent(in) :: iblank(:)
       real(SCALAR_KIND), intent(in) :: time
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine addSolenoidalExcitation

  end interface

end module Patch_mod
