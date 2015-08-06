#include "config.h"

module Patch_factory

  use Patch_mod, only : t_Patch

  implicit none
  private

  type, public :: t_PatchFactory

     class(t_Patch), pointer :: patch => null()
     character(len = STRING_LENGTH) :: patchType = ""

   contains

     procedure, pass :: connect
     procedure, pass :: cleanup

  end type t_PatchFactory

contains

  subroutine connect(this, patchTarget, patchType, createNew)

    ! <<< Derived types >>>
    use ProbePatch_mod, only : t_ProbePatch
    use SpongePatch_mod, only : t_SpongePatch
    use ActuatorPatch_mod, only : t_ActuatorPatch
    use IsothermalWall_mod, only : t_IsothermalWall
    use CostTargetPatch_mod, only : t_CostTargetPatch
    use ImpenetrableWall_mod, only : t_ImpenetrableWall
    use InflowOutflowPatch_mod, only : t_InflowOutflowPatch
    use BlockInterfacePatch_mod, only : t_BlockInterfacePatch
    use SolenoidalExcitationPatch_mod, only : t_SolenoidalExcitationPatch

    implicit none

    ! <<< Arguments >>>
    class(t_PatchFactory) :: this
    class(t_Patch), pointer, intent(out) :: patchTarget
    character(len = *), intent(in), optional :: patchType
    logical, intent(in), optional :: createNew

    ! <<< Local variables >>>
    logical :: createNew_

    createNew_ = .false.
    if (present(createNew)) createNew_ = createNew

    if (present(patchType) .and. .not. (associated(this%patch) .and. .not. createNew_)) then

       if (associated(this%patch)) then
          call this%patch%cleanup()
          deallocate(this%patch)
       end if

       nullify(this%patch)

       this%patchType = patchType

       select case (trim(patchType))

       case ('PROBE')
          allocate(t_ProbePatch :: this%patch)

       case ('SPONGE')
          allocate(t_SpongePatch :: this%patch)

       case ('ACTUATOR')
          allocate(t_ActuatorPatch :: this%patch)

       case ('ISOTHERMAL_WALL')
          allocate(t_IsothermalWall :: this%patch)

       case ('COST_TARGET')
          allocate(t_CostTargetPatch :: this%patch)

       case ('IMPENETRABLE_WALL')
          allocate(t_ImpenetrableWall :: this%patch)

       case ('INFLOW_OUTFLOW')
          allocate(t_InflowOutflowPatch :: this%patch)

       case ('BLOCK_INTERFACE')
          allocate(t_BlockInterfacePatch :: this%patch)

       case ('SOLENOIDAL_EXCITATION')
          allocate(t_SolenoidalExcitationPatch :: this%patch)

       case default
          this%patchType = ""

       end select

    end if

    nullify(patchTarget)
    if (.not. associated(this%patch)) return
    patchTarget => this%patch

  end subroutine connect

  subroutine cleanup(this)

    class(t_PatchFactory) :: this

    if (associated(this%patch)) then
       call this%patch%cleanup()
       deallocate(this%patch)
    end if
    nullify(this%patch)

  end subroutine cleanup

end module Patch_factory
