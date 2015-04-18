#include "config.h"

subroutine connectReverseMigrator(this, reverseMigratorTarget, reverseMigratorType, createNew)

  ! <<< Derived types >>>
  use ReverseMigrator_mod, only : t_ReverseMigrator
  use ReverseMigrator_factory, only : t_ReverseMigratorFactory
  use UniformCheckpointer_mod, only : t_UniformCheckpointer

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

  if (present(reverseMigratorType) .and. .not. (associated(this%reverseMigrator) .and.       &
       .not. createNew_)) then

     if (associated(this%reverseMigrator)) deallocate(this%reverseMigrator)
     nullify(this%reverseMigrator)

     select case (trim(reverseMigratorType))

     case ('uniform checkpointing')
        allocate(t_UniformCheckpointer :: this%reverseMigrator)

     end select

  end if

  nullify(reverseMigratorTarget)
  if (.not. associated(this%reverseMigrator)) return
  reverseMigratorTarget => this%reverseMigrator

end subroutine connectReverseMigrator

subroutine cleanupReverseMigratorFactory(this)

  ! <<< Derived types >>>
  use ReverseMigrator_factory, only : t_ReverseMigratorFactory

  implicit none

  ! <<< Arguments >>>
  class(t_ReverseMigratorFactory) :: this

  if (associated(this%reverseMigrator)) then
     call this%reverseMigrator%cleanup()
     deallocate(this%reverseMigrator)
  end if
  nullify(this%reverseMigrator)

end subroutine cleanupReverseMigratorFactory
