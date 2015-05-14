#include "config.h"

module MappingFunction_factory

  use MappingFunction_mod, only : t_MappingFunction

  implicit none
  private

  type, public :: t_MappingFunctionFactory

     class(t_MappingFunction), pointer :: mappingFunction => null()

   contains

     procedure, pass :: connect => connectMappingFunction
     procedure, pass :: cleanup => cleanupMappingFunctionFactory

  end type t_MappingFunctionFactory

  interface

     subroutine connectMappingFunction(this, mappingFunctionTarget,                          &
          mappingFunctionType, createNew)

       use MappingFunction_mod, only : t_MappingFunction

       import :: t_MappingFunctionFactory

       class(t_MappingFunctionFactory) :: this
       class(t_MappingFunction), pointer, intent(out) :: mappingFunctionTarget
       character(len = *), intent(in), optional :: mappingFunctionType
       logical, intent(in), optional :: createNew

     end subroutine connectMappingFunction

  end interface

  interface

     subroutine cleanupMappingFunctionFactory(this)

       import :: t_MappingFunctionFactory

       class(t_MappingFunctionFactory) :: this

     end subroutine cleanupMappingFunctionFactory

  end interface

end module MappingFunction_factory
