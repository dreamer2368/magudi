#include "config.h"

module ReverseMigrator_factory

  use ReverseMigrator_mod, only : t_ReverseMigrator

  implicit none
  private

  type, public :: t_ReverseMigratorFactory

     class(t_ReverseMigrator), pointer :: reverseMigrator => null()
     character(len = STRING_LENGTH) :: reverseMigratorType = ""

   contains

     procedure, pass :: connect
     procedure, pass :: cleanup

  end type t_ReverseMigratorFactory

contains

  subroutine connect(this, reverseMigratorTarget, reverseMigratorType, createNew)

    ! <<< Derived types >>>
    use UniformCheckpointing_mod, only : t_UniformCheckpointing

    implicit none

    ! <<< Arguments >>>
    class(t_ReverseMigratorFactory) :: this
    class(t_ReverseMigrator), pointer, intent(out) :: reverseMigratorTarget
    character(len = *), intent(in), optional :: reverseMigratorType
    logical, intent(in), optional :: createNew

    ! <<< Local variables >>>
    logical :: createNew_

    createNew_ = .false.
    if (present(createNew)) createNew_ = createNew

    if (present(reverseMigratorType) .and. .not.                                             &
         (associated(this%reverseMigrator) .and. .not. createNew_)) then

       if (associated(this%reverseMigrator)) then
          call this%reverseMigrator%cleanup()
          deallocate(this%reverseMigrator)
       end if

       nullify(this%reverseMigrator)

       this%reverseMigratorType = reverseMigratorType

       select case (trim(reverseMigratorType))

       case ('UNIFORM_CHECKPOINTING')
          allocate(t_UniformCheckpointing :: this%reverseMigrator)

       case default
          this%reverseMigratorType = ""

       end select

    end if

    nullify(reverseMigratorTarget)
    if (.not. associated(this%reverseMigrator)) return
    reverseMigratorTarget => this%reverseMigrator

  end subroutine connect

  subroutine cleanup(this)

    class(t_ReverseMigratorFactory) :: this

    if (associated(this%reverseMigrator)) then
       call this%reverseMigrator%cleanup()
       deallocate(this%reverseMigrator)
    end if
    nullify(this%reverseMigrator)

  end subroutine cleanup

end module ReverseMigrator_factory
