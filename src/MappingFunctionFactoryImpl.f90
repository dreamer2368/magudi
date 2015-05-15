#include "config.h"

subroutine connectMappingFunction(this, mappingFunctionTarget, mappingFunctionType, createNew)

  ! <<< Derived types >>>
  use UniformMap_mod, only : t_UniformMap
  use MappingFunction_mod, only : t_MappingFunction
  use MappingFunction_factory, only : t_MappingFunctionFactory

  implicit none

  ! <<< Arguments >>>
  class(t_MappingFunctionFactory) :: this
  class(t_MappingFunction), pointer, intent(out) :: mappingFunctionTarget
  character(len = *), intent(in), optional :: mappingFunctionType
  logical, intent(in), optional :: createNew

  ! <<< Local variables >>>
  logical :: createNew_

  createNew_ = .false.
  if (present(createNew)) createNew_ = createNew

  if (present(mappingFunctionType) .and. .not. (associated(this%mappingFunction) .and.       &
       .not. createNew_)) then

     if (associated(this%mappingFunction)) deallocate(this%mappingFunction)
     nullify(this%mappingFunction)

     select case (trim(mappingFunctionType))

     case ('UNIFORM')
        allocate(t_UniformMap :: this%mappingFunction)

     end select

  end if

  nullify(mappingFunctionTarget)
  if (.not. associated(this%mappingFunction)) return
  mappingFunctionTarget => this%mappingFunction

end subroutine connectMappingFunction

subroutine cleanupMappingFunctionFactory(this)

  ! <<< Derived types >>>
  use MappingFunction_factory, only : t_MappingFunctionFactory

  implicit none

  ! <<< Arguments >>>
  class(t_MappingFunctionFactory) :: this

  if (associated(this%mappingFunction)) then
     deallocate(this%mappingFunction)
  end if
  nullify(this%mappingFunction)

end subroutine cleanupMappingFunctionFactory
