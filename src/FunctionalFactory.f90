#include "config.h"

module Functional_factory

  use Functional_mod, only : t_Functional

  implicit none
  private

  type, public :: t_FunctionalFactory

     class(t_Functional), pointer :: functional => null()
     character(len = STRING_LENGTH) :: functionalType = ""

   contains

     procedure, pass :: connect
     procedure, pass :: cleanup

  end type t_FunctionalFactory

contains

  subroutine connect(this, functionalTarget, functionalType, createNew)

    ! <<< Derived types >>>
    use AcousticNoise_mod, only : t_AcousticNoise

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

    if (present(functionalType) .and. .not.                                                  &
         (associated(this%functional) .and. .not. createNew_)) then

       if (associated(this%functional)) then
          call this%functional%cleanup()
          deallocate(this%functional)
       end if

       nullify(this%functional)

       this%functionalType = functionalType

       select case (trim(functionalType))

       case ('ACOUSTIC_NOISE')
          allocate(t_AcousticNoise :: this%functional)

       case default
          this%functionalType = ""

       end select

    end if

    nullify(functionalTarget)
    if (.not. associated(this%functional)) return
    functionalTarget => this%functional

  end subroutine connect

  subroutine cleanup(this)

    class(t_FunctionalFactory) :: this

    if (associated(this%functional)) then
       call this%functional%cleanup()
       deallocate(this%functional)
    end if
    nullify(this%functional)

  end subroutine cleanup

end module Functional_factory
