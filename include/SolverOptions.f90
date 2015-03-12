#include "config.h"

module SolverOptions_type

  implicit none
  private

  type, public :: t_SolverOptions

     real(SCALAR_KIND) :: reynoldsNumber, prandtlNumber, ratioOfSpecificHeats,               &
          powerLawExponent, dissipationAmount, densityRange(2), temperatureRange(2), cfl,    &
          timeStepSize, spongeAmount, farFieldInviscidPenaltyAmount,                         &
          farFieldViscousPenaltyAmount
     integer :: spongeExponent

  end type t_SolverOptions

end module SolverOptions_type

module SolverOptions_mod

  implicit none

  interface

     subroutine initializeSolverOptions(this, simulationFlags)

       use SolverOptions_type
       use SimulationFlags_type

       type(t_SolverOptions), intent(out) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags

     end subroutine initializeSolverOptions

  end interface

  interface

     subroutine updateSolverOptions(this, simulationFlags, patchData)

       use SolverOptions_type
       use PatchDescriptor_type
       use SimulationFlags_type

       type(t_SolverOptions), intent(out) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_PatchDescriptor), intent(in) :: patchData(:)

     end subroutine updateSolverOptions

  end interface

end module SolverOptions_mod
