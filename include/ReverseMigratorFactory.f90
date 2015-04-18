#include "config.h"

module ReverseMigrator_factory

  use ReverseMigrator_mod, only : t_ReverseMigrator

  implicit none
  private

  type, public :: t_ReverseMigratorFactory

     class(t_ReverseMigrator), pointer :: reverseMigrator => null()

   contains

     procedure, pass :: connect => connectReverseMigrator
     procedure, pass :: cleanup => cleanupReverseMigratorFactory

  end type t_ReverseMigratorFactory

  interface

     subroutine connectReverseMigrator(this, reverseMigratorTarget,                          &
          reverseMigratorType, createNew)

       use ReverseMigrator_mod, only : t_ReverseMigrator

       import :: t_ReverseMigratorFactory

       class(t_ReverseMigratorFactory) :: this
       class(t_ReverseMigrator), pointer, intent(out) :: reverseMigratorTarget
       character(len = *), intent(in), optional :: reverseMigratorType
       logical, intent(in), optional :: createNew

     end subroutine connectReverseMigrator

  end interface

  interface

     subroutine cleanupReverseMigratorFactory(this)

       import :: t_ReverseMigratorFactory

       class(t_ReverseMigratorFactory) :: this

     end subroutine cleanupReverseMigratorFactory

  end interface

end module ReverseMigrator_factory
