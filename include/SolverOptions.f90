#include "config.h"

module SolverOptions_type

  implicit none
  private

  type, public :: t_SolverOptions

     real(SCALAR_KIND) :: reynoldsNumberInverse, prandtlNumberInverse, ratioOfSpecificHeats, &
          powerLawExponent, bulkViscosityRatio, dissipationAmount, densityRange(2),          &
          temperatureRange(2), cfl, timeStepSize, convergenceTolerance(3)

  end type t_SolverOptions

end module SolverOptions_type

module SolverOptions_mod

  implicit none
  public

  interface

     subroutine initializeSolverOptions(this, simulationFlags)

       use SolverOptions_type
       use SimulationFlags_type

       type(t_SolverOptions), intent(out) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags

     end subroutine initializeSolverOptions

  end interface

end module SolverOptions_mod
