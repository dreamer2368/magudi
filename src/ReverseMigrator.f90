#include "config.h"

module ReverseMigrator_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  implicit none
  private

  type, abstract, public :: t_ReverseMigrator

     integer :: nStages
     integer :: startTimestep, endTimestep, saveInterval,                                    &
          numIntermediateStates, loadedTimestep = -1
     character(len = STRING_LENGTH) :: outputPrefix

   contains

     procedure, non_overridable, pass :: setupBase
     procedure, non_overridable, pass :: cleanupBase

     procedure(setup), pass, deferred :: setup
     procedure(cleanup), pass, deferred :: cleanup
     procedure(migrateTo), pass, deferred :: migrateTo

  end type t_ReverseMigrator

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

     subroutine migrateTo(this, region, timeIntegrator, timestep, stage)

       use Region_mod, only : t_Region
       use TimeIntegrator_mod, only : t_TimeIntegrator

       import :: t_ReverseMigrator

       class(t_ReverseMigrator) :: this
       class(t_Region) :: region
       class(t_TimeIntegrator) :: timeIntegrator
       integer, intent(in) :: timestep, stage

     end subroutine migrateTo

  end interface

contains

  subroutine setupBase(this, region, timeIntegrator, outputPrefix, startTimestep,            &
       endTimestep, saveInterval, numIntermediateStates)

    ! <<< Derived types >>>
    use Region_mod, only : t_Region
    use TimeIntegrator_mod, only : t_TimeIntegrator

    implicit none

    ! <<< Arguments >>>
    class(t_ReverseMigrator) :: this
    class(t_Region), intent(in) :: region
    class(t_TimeIntegrator), intent(in) :: timeIntegrator
    character(len = *), intent(in) :: outputPrefix
    integer, intent(in) :: startTimestep, endTimestep, saveInterval, numIntermediateStates

    assert(len_trim(outputPrefix) > 0)
    assert(startTimestep >= 0)
    assert(endTimestep >= startTimestep)
    assert(saveInterval > 0)
    assert(mod(endTimestep - startTimestep, saveInterval) == 0)
    assert(numIntermediateStates >= 0)
    assert(timeIntegrator%nStages > 0)
    assert(mod(numIntermediateStates, timeIntegrator%nStages) == 0)
    assert(allocated(region%states))
    assert(size(region%states) > 0)

    call this%cleanupBase()

    this%outputPrefix = outputPrefix

    this%startTimestep = startTimestep
    this%endTimestep = endTimestep
    this%saveInterval = saveInterval
    this%numIntermediateStates = numIntermediateStates

  end subroutine setupBase

  subroutine cleanupBase(this)

    implicit none

    ! <<< Arguments >>>
    class(t_ReverseMigrator) :: this

    this%loadedTimestep = -1

  end subroutine cleanupBase

end module ReverseMigrator_mod
