#include "config.h"

module SolverOptions_mod

  implicit none
  private

  type, public :: t_SolverOptions

     real(SCALAR_KIND) :: reynoldsNumberInverse, prandtlNumberInverse, ratioOfSpecificHeats, &
          powerLawExponent, bulkViscosityRatio, dissipationAmount, densityRange(2),          &
          temperatureRange(2), cfl, timeStepSize
     integer :: nUnknowns
     character(len = STRING_LENGTH) :: discretizationType, timeintegratorType,               &
          costFunctionalType, controllerType, checkpointingScheme

   contains

     procedure, pass :: initialize => initializeSolverOptions

  end type t_SolverOptions

  interface

     subroutine initializeSolverOptions(this, nDimensions, simulationFlags, comm)

       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_SolverOptions

       class(t_SolverOptions), intent(out) :: this
       integer, intent(in) :: nDimensions
       type(t_SimulationFlags), intent(in) :: simulationFlags
       integer, intent(in), optional :: comm

     end subroutine initializeSolverOptions

  end interface

end module SolverOptions_mod
