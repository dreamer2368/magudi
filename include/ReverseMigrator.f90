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

     integer :: algorithm, numIntermediateStates, nStages
     type(t_IntermediateStorage), allocatable :: temp_(:)

  end type t_ReverseMigrator

end module ReverseMigrator_type

module ReverseMigrator_mod

  implicit none
  public

  interface

     subroutine setupReverseMigrator(this, algorithm, region, nStages, numIntermediateStates)

       use Region_type, only : t_Region
       use ReverseMigrator_type, only : t_ReverseMigrator

       type(t_ReverseMigrator) :: this
       integer, intent(in) :: algorithm
       type(t_Region), intent(in) :: region
       integer, intent(in) :: nStages, numIntermediateStates

     end subroutine setupReverseMigrator

  end interface

  interface

     subroutine cleanupReverseMigrator(this)

       use ReverseMigrator_type, only : t_ReverseMigrator

       type(t_ReverseMigrator) :: this

     end subroutine cleanupReverseMigrator

  end interface

  interface

     subroutine migrateToTimestep(this, region, timestep, outputPrefix)

       use Region_type, only : t_Region
       use ReverseMigrator_type, only : t_ReverseMigrator

       type(t_ReverseMigrator) :: this
       type(t_Region) :: region
       integer, intent(in) :: timestep
       character(len = *), intent(in) :: outputPrefix

     end subroutine migrateToTimestep

  end interface

end module ReverseMigrator_mod
