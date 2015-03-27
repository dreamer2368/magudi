#include "config.h"

module ReverseMigrator_type

  implicit none
  private

  integer, parameter, public :: &
       UNIFORM_CHECKPOINTING = 1

  type, private :: t_IntermediateStorage

     SCALAR_TYPE, allocatable :: buffer(:,:,:)

  end type t_IntermediateStorage

  type, public :: t_ReverseMigrator

     integer :: nStages = 4
     integer :: algorithm, startTimestep, endTimestep, saveInterval,                         &
          numIntermediateStates, loadedTimestep = -1
     character(len = STRING_LENGTH) :: outputPrefix
     type(t_IntermediateStorage), allocatable :: temp_(:)

  end type t_ReverseMigrator

end module ReverseMigrator_type

module ReverseMigrator_mod

  implicit none
  public

  interface

     subroutine setupReverseMigrator(this, region, outputPrefix, algorithm,                  &
          startTimestep, endTimestep, saveInterval, numIntermediateStates)

       use Region_type, only : t_Region
       use ReverseMigrator_type, only : t_ReverseMigrator

       type(t_ReverseMigrator) :: this
       type(t_Region), intent(in) :: region
       character(len = *), intent(in) :: outputPrefix, algorithm
       integer, intent(in) :: startTimestep, endTimestep, saveInterval, numIntermediateStates

     end subroutine setupReverseMigrator

  end interface

  interface

     subroutine cleanupReverseMigrator(this)

       use ReverseMigrator_type, only : t_ReverseMigrator

       type(t_ReverseMigrator) :: this

     end subroutine cleanupReverseMigrator

  end interface

  interface

     subroutine migrateToSubstep(this, region, integrator, timestep, stage)
       
       use Region_type, only : t_Region
       use RK4Integrator_type, only : t_RK4Integrator
       use ReverseMigrator_type, only : t_ReverseMigrator
       
       type(t_ReverseMigrator) :: this
       type(t_Region) :: region
       type(t_RK4Integrator) :: integrator
       integer, intent(in) :: timestep, stage
       
     end subroutine migrateToSubstep

  end interface

end module ReverseMigrator_mod
