#include "config.h"

subroutine setupDragCoefficient(this, region)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use DragCoefficient_mod, only : t_DragCoefficient

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_DragCoefficient) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: nDimensions
  character(len = STRING_LENGTH) :: message

  nDimensions = region%grids(1)%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  call this%cleanup()

  call this%setupBase(region%simulationFlags, region%solverOptions)

  this%freeStreamPressure = getOption("free_stream_pressure",                                &
       1.0_wp / region%solverOptions%ratioOfSpecificHeats)
  if (this%freeStreamPressure <= 0.0_wp) then
     write(message, '(A,E11.2,A)') "Free stream pressure is invalid: ", &
          this%freeStreamPressure, " <= 0!"
     call gracefulExit(region%comm, message)
  end if

  this%direction = 0.0_wp

  this%direction(1) = getOption("free_stream_velocity_x", 0.0_wp)
  if (nDimensions >= 2) this%direction(2) = getOption("free_stream_velocity_y", 0.0_wp)
  if (nDimensions == 3) this%direction(3) = getOption("free_stream_velocity_z", 0.0_wp)

  if (sum(this%direction ** 2) <= epsilon(0.0_wp)) then
     write(message, '(A)') "Unable to determine a unit vector along the free stream flow &
          &direction required for computing the drag coefficient!"
     call gracefulExit(region%comm, message)
  end if

end subroutine setupDragCoefficient

subroutine cleanupDragCoefficient(this)

  ! <<< Derived types >>>
  use DragCoefficient_mod, only : t_DragCoefficient

  implicit none

  ! <<< Arguments >>>
  class(t_DragCoefficient) :: this

  ! <<< Local variables >>>
  integer :: i

  call this%cleanupBase()

end subroutine cleanupDragCoefficient

function computeDragCoefficient(this, region) result(instantaneousFunctional)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use DragCoefficient_mod, only : t_DragCoefficient

  ! <<< Arguments >>>
  class(t_DragCoefficient) :: this
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousFunctional

end function computeDragCoefficient

subroutine addDragCoefficientAdjointForcing(this, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use DragCoefficient_mod, only : t_DragCoefficient
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Arguments >>>
  class(t_DragCoefficient) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

end subroutine addDragCoefficientAdjointForcing
