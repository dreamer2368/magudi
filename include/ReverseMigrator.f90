#include "config.h"

module ReverseMigrator_mod

  implicit none

  type, abstract, public :: t_ReverseMigrator

     integer :: nStages
     integer :: startTimestep, endTimestep, saveInterval,                                    &
          numIntermediateStates, loadedTimestep = -1
     character(len = STRING_LENGTH) :: outputPrefix

   contains

     procedure, non_overridable, pass :: setupBase => setupReverseMigrator
     procedure, non_overridable, pass :: cleanupBase => cleanupReverseMigrator

     procedure(setup), pass, deferred :: setup
     procedure(cleanup), pass, deferred :: cleanup
     procedure(migrateTo), pass, deferred :: migrateTo

  end type t_ReverseMigrator

  interface

     subroutine setupReverseMigrator(this, region, timeIntegrator, outputPrefix,             &
          startTimestep, endTimestep, saveInterval,                                          &
          numIntermediateStates)

       use Region_mod, only : t_Region
       use TimeIntegrator_mod, only : t_TimeIntegrator

       import :: t_ReverseMigrator

       class(t_ReverseMigrator) :: this
       class(t_Region), intent(in) :: region
       class(t_TimeIntegrator), intent(in) :: timeIntegrator
       character(len = *), intent(in) :: outputPrefix
       integer, intent(in) :: startTimestep, endTimestep, saveInterval, numIntermediateStates

     end subroutine setupReverseMigrator

  end interface

  interface

     subroutine cleanupReverseMigrator(this)

       import :: t_ReverseMigrator

       class(t_ReverseMigrator) :: this

     end subroutine cleanupReverseMigrator

  end interface

  abstract interface

     subroutine setup(this, region, timeIntegrator, outputPrefix, startTimestep,             &
          endTimestep, saveInterval, numIntermediateStates)

       use Region_mod, only : t_Region
       use TimeIntegrator_mod, only : t_TimeIntegrator

       import :: t_ReverseMigrator

       class(t_ReverseMigrator) :: this
       class(t_Region), intent(in) :: region
       class(t_TimeIntegrator), intent(in) :: timeIntegrator
       character(len = *), intent(in) :: outputPrefix
       integer, intent(in) :: startTimestep, endTimestep, saveInterval, numIntermediateStates

     end subroutine setup

  end interface

  abstract interface

     subroutine cleanup(this)

       import :: t_ReverseMigrator

       class(t_ReverseMigrator) :: this

     end subroutine cleanup

  end interface

  abstract interface

     subroutine migrateTo(this, region, controller, timeIntegrator, timestep, stage)

       use Region_mod, only : t_Region
       use Controller_mod, only : t_Controller
       use TimeIntegrator_mod, only : t_TimeIntegrator

       import :: t_ReverseMigrator

       class(t_ReverseMigrator) :: this
       class(t_Region) :: region
       class(t_Controller) :: controller
       class(t_TimeIntegrator) :: timeIntegrator
       integer, intent(in) :: timestep, stage

     end subroutine migrateTo

  end interface

end module ReverseMigrator_mod
