#include "config.h"

subroutine connectFunctional(this, functionalTarget, functionalType, createNew)

  ! <<< Derived types >>>
  use DragForce_mod, only : t_DragForce
  use Functional_mod, only : t_Functional
  use PressureDrag_mod, only : t_PressureDrag
  use AcousticNoise_mod, only : t_AcousticNoise
  use ReynoldsStress_mod, only : t_ReynoldsStress
  use FlameTemperature_mod, only : t_FlameTemperature
  use ReactantDepletion_mod, only : t_ReactantDepletion
  use Functional_factory, only : t_FunctionalFactory

  implicit none

  ! <<< Arguments >>>
  class(t_FunctionalFactory) :: this
  class(t_Functional), pointer, intent(out) :: functionalTarget
  character(len = *), intent(in), optional :: functionalType
  logical, intent(in), optional :: createNew

  ! <<< Local variables >>>
  logical :: createNew_

  createNew_ = .false.
  if (present(createNew)) createNew_ = createNew

  if (present(functionalType) .and. .not. (associated(this%functional) .and.                 &
       .not. createNew_)) then

     if (associated(this%functional)) deallocate(this%functional)
     nullify(this%functional)

     this%functionalType = functionalType

     select case (trim(functionalType))

     case ('SOUND')
        allocate(t_AcousticNoise :: this%functional)

     case ('PRESSURE_DRAG')
        allocate(t_PressureDrag :: this%functional)

     case ('DRAG')
        allocate(t_DragForce :: this%functional)

     case ('REYNOLDS_STRESS')
        allocate(t_ReynoldsStress :: this%functional)

     case ('TEMPERATURE')
        allocate(t_FlameTemperature :: this%functional)

     case ('REACTANT')
        allocate(t_ReactantDepletion :: this%functional)

     case default
        this%functionalType = ""

     end select

  end if

  nullify(functionalTarget)
  if (.not. associated(this%functional)) return
  functionalTarget => this%functional

end subroutine connectFunctional

subroutine cleanupFunctionalFactory(this)

  ! <<< Derived types >>>
  use Functional_factory, only : t_FunctionalFactory

  implicit none

  ! <<< Arguments >>>
  class(t_FunctionalFactory) :: this

  if (associated(this%functional)) then
     call this%functional%cleanup()
     deallocate(this%functional)
  end if
  nullify(this%functional)

end subroutine cleanupFunctionalFactory
