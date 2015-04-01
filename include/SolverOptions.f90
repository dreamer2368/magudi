#include "config.h"

module SolverOptions_enum

  implicit none
  public

  integer, parameter ::                                                                      &
       SOUND = 1,                                                                            &
       LIFT  = 2,                                                                            &
       DRAG  = 3
  
end module SolverOptions_enum

module SolverOptions_mod

  implicit none
  private

  type, public :: t_SolverOptions

     real(SCALAR_KIND) :: reynoldsNumberInverse, prandtlNumberInverse, ratioOfSpecificHeats, &
          powerLawExponent, bulkViscosityRatio, dissipationAmount, densityRange(2),          &
          temperatureRange(2), cfl, timeStepSize, convergenceTolerance(3)
     integer :: costFunctionalType

   contains

     procedure, pass :: initialize => initializeSolverOptions

  end type t_SolverOptions

  interface

     subroutine initializeSolverOptions(this, simulationFlags, comm)

       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_SolverOptions

       class(t_SolverOptions), intent(out) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       integer, intent(in), optional :: comm

     end subroutine initializeSolverOptions

  end interface

end module SolverOptions_mod
