#include "config.h"

subroutine setupFunctional(this, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Functional_mod, only : t_Functional
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_Functional) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

end subroutine setupFunctional

subroutine cleanupFunctional(this)

  ! <<< Derived types >>>
  use Functional_mod, only : t_Functional

  implicit none

  ! <<< Arguments >>>
  class(t_Functional) :: this

end subroutine cleanupFunctional
