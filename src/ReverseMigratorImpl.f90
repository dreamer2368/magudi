#include "config.h"

subroutine setupReverseMigrator(this, region, timeIntegrator, outputPrefix, startTimestep,   &
     endTimestep, saveInterval, numIntermediateStates)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use TimeIntegrator_mod, only : t_TimeIntegrator
  use ReverseMigrator_mod, only : t_ReverseMigrator

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

end subroutine setupReverseMigrator

subroutine cleanupReverseMigrator(this)

  ! <<< Derived types >>>
  use ReverseMigrator_mod, only : t_ReverseMigrator

  implicit none

  ! <<< Arguments >>>
  class(t_ReverseMigrator) :: this

  this%loadedTimestep = -1

end subroutine cleanupReverseMigrator
