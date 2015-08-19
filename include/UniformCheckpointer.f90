#include "config.h"

module UniformCheckpointer_mod

  use ReverseMigrator_mod, only : t_ReverseMigrator

  implicit none

  type, private :: t_IntermediateStorage
     SCALAR_TYPE, allocatable :: buffer(:,:,:)
     real(SCALAR_KIND), allocatable :: times(:)
  end type t_IntermediateStorage

  type, extends(t_ReverseMigrator), public :: t_UniformCheckpointer

     type(t_IntermediateStorage), allocatable :: data_(:)

   contains

     procedure, pass :: setup => setupUniformCheckpointer
     procedure, pass :: cleanup => cleanupUniformCheckpointer
     procedure, pass :: migrateTo => uniformCheckpointingMigrateTo

  end type t_UniformCheckpointer

  interface

     subroutine setupUniformCheckpointer(this, region, timeIntegrator, outputPrefix,         &
          startTimestep, endTimestep, saveInterval, numIntermediateStates)

       use Region_mod, only : t_Region
       use TimeIntegrator_mod, only : t_TimeIntegrator

       import :: t_UniformCheckpointer

       class(t_UniformCheckpointer) :: this
       class(t_Region), intent(in) :: region
       class(t_TimeIntegrator), intent(in) :: timeIntegrator
       character(len = *), intent(in) :: outputPrefix
       integer, intent(in) :: startTimestep, endTimestep, saveInterval, numIntermediateStates

     end subroutine setupUniformCheckpointer

  end interface

  interface

     subroutine cleanupUniformCheckpointer(this)

       import :: t_UniformCheckpointer

       class(t_UniformCheckpointer) :: this

     end subroutine cleanupUniformCheckpointer

  end interface

  interface

     subroutine uniformCheckpointingMigrateTo(this, region, timeIntegrator, timestep, stage, &
          controller)

       use Region_mod, only : t_Region
       use Controller_mod, only : t_Controller
       use TimeIntegrator_mod, only : t_TimeIntegrator

       import :: t_UniformCheckpointer

       class(t_UniformCheckpointer) :: this
       class(t_Region) :: region
       class(t_TimeIntegrator) :: timeIntegrator
       class(t_Controller) :: controller
       integer, intent(in) :: timestep, stage

     end subroutine uniformCheckpointingMigrateTo

  end interface

end module UniformCheckpointer_mod
