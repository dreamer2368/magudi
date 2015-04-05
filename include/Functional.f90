#include "config.h"

module Functional_mod

  implicit none
  private

  type, abstract, public :: t_Functional

     SCALAR_TYPE :: cachedValue = real(0.0, SCALAR_KIND)

   contains

     procedure, non_overridable, pass :: setupBase => setupFunctional
     procedure, non_overridable, pass :: cleanupBase => cleanupFunctional
     procedure, pass :: writeToFile => writeFunctionalToFile

     procedure(setup), pass, deferred :: setup
     procedure(cleanup), pass, deferred :: cleanup
     procedure(compute), pass, deferred :: compute
     procedure(computeAdjointForcing), pass, deferred :: computeAdjointForcing
     procedure(isPatchValid), pass, deferred :: isPatchValid

  end type t_Functional

  abstract interface

     subroutine setup(this, region)

       use Region_mod, only : t_Region

       import :: t_Functional

       class(t_Functional) :: this
       class(t_Region) :: region

     end subroutine setup

  end interface

  abstract interface

     subroutine cleanup(this)

       import :: t_Functional

       class(t_Functional) :: this

     end subroutine cleanup

  end interface

  abstract interface

     function compute(this, time, region) result(instantaneousFunctional)

       use Region_mod, only : t_Region

       import :: t_Functional

       class(t_Functional) :: this
       real(SCALAR_KIND), intent(in) :: time
       class(t_Region), intent(in) :: region

       SCALAR_TYPE :: instantaneousFunctional

     end function compute

  end interface

  abstract interface

     subroutine computeAdjointForcing(this, region)

       use Region_mod, only : t_Region

       import :: t_Functional

       class(t_Functional) :: this
       class(t_Region) :: region

     end subroutine computeAdjointForcing

  end interface

  abstract interface

     function isPatchValid(this, patchDescriptor, gridSize, normalDirection,                 &
          extent, simulationFlags, message)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Functional

       class(t_Functional) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchValid

     end function isPatchValid

  end interface

  interface

     subroutine setupFunctional(this, simulationFlags, solverOptions)

       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Functional

       class(t_Functional) :: this
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupFunctional

  end interface

  interface

     subroutine cleanupFunctional(this)

       import :: t_Functional

       class(t_Functional) :: this

     end subroutine cleanupFunctional

  end interface

  interface

     subroutine writeFunctionalToFile(this, comm, filename, timestep, time, append)

       import :: t_Functional

       class(t_Functional) :: this
       integer, intent(in) :: comm
       character(len = *), intent(in) :: filename
       integer, intent(in) :: timestep
       real(SCALAR_KIND), intent(in) :: time
       logical, intent(in), optional :: append

     end subroutine writeFunctionalToFile

  end interface

end module Functional_mod
