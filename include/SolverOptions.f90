#include "config.h"

module SolverOptions_enum

  implicit none
  public

  integer, parameter, public ::                                                              &
       IDEAL_GAS         = 1,                                                                &
       IDEAL_GAS_MIXTURE = 2

end module SolverOptions_enum

module SolverOptions_mod

  implicit none

  type, public :: t_SolverOptions

     real(SCALAR_KIND) :: reynoldsNumberInverse, prandtlNumberInverse,                       &
          ratioOfSpecificHeats, powerLawExponent, bulkViscosityRatio, dissipationAmount,     &
          densityRange(2),temperatureRange(2), massFractionRange(2), cfl, timeStepSize
     SCALAR_TYPE, dimension(:), allocatable :: froudeNumberInverse, schmidtNumberInverse,    &
          molecularWeightInverse
     integer :: nSpecies, nUnknowns, equationOfState, nParticles
     character(len = STRING_LENGTH) :: discretizationType, timeintegratorType,               &
          costFunctionalType, controllerType, checkpointingScheme
     character(len = STRING_LENGTH), dimension(:), allocatable :: speciesName

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
