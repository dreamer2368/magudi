#include "config.h"

module Region_enum

  implicit none
  public

  integer, parameter ::                                                                      &
       FORWARD = +1,                                                                         &
       ADJOINT = -1,                                                                         &
       !SeungWhan: Optimization flag
       ONESTEP = 0

end module Region_enum

module Region_mod

  use MPI, only : MPI_COMM_NULL

  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use Patch_factory, only : t_PatchFactory
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  type, public :: t_Region

     type(t_Grid), allocatable :: grids(:)
     type(t_State), allocatable :: states(:)
     type(t_PatchFactory), allocatable :: patchFactories(:)
     type(t_SolverOptions) :: solverOptions
     type(t_SimulationFlags) :: simulationFlags
     type(t_PatchDescriptor), allocatable :: patchData(:)
     integer :: comm = MPI_COMM_NULL, commGridMasters = MPI_COMM_NULL, timestep = 0
     integer, allocatable :: globalGridSizes(:,:), processDistributions(:,:),                &
          gridCommunicators(:), patchCommunicators(:), patchInterfaces(:),                   &
          interfaceIndexReorderings(:,:), patchMasterRanks(:)
     logical :: outputOn = .true.

   contains

     procedure, pass :: setup => setupRegion
     procedure, pass :: cleanup => cleanupRegion
     procedure, pass :: setupBoundaryConditions
     procedure, pass :: loadData => loadRegionData
     procedure, pass :: saveData => saveRegionData
     procedure, pass :: getCfl
     procedure, pass :: getTimeStepSize
     procedure, pass :: reportGridDiagnostics
     procedure, pass :: computeRhs
     procedure, pass :: saveSpongeStrength
     procedure, pass :: resetProbes
     procedure, pass :: saveProbeData

  end type t_Region

  interface

     subroutine setupRegion(this, comm, globalGridSizes, simulationFlags,                   &
           solverOptions, verbose)

       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Region

       class(t_Region) :: this
       integer, intent(in) :: comm, globalGridSizes(:,:)

       type(t_SimulationFlags), intent(in), optional :: simulationFlags
       type(t_SolverOptions), intent(in), optional :: solverOptions
       logical, intent(in), optional :: verbose

     end subroutine setupRegion

  end interface

  interface

     subroutine cleanupRegion(this)

       import :: t_Region

       class(t_Region) :: this

     end subroutine cleanupRegion

  end interface

  interface

     subroutine setupBoundaryConditions(this, boundaryConditionFilename)

       import :: t_Region

       class(t_Region) :: this
       character(len = *), intent(in) :: boundaryConditionFilename

     end subroutine setupBoundaryConditions

  end interface

  interface

     subroutine loadRegionData(this, quantityOfInterest, filename, speciesFilename)

       import :: t_Region

       class(t_Region) :: this
       integer, intent(in) :: quantityOfInterest
       character(len = *), intent(in) :: filename

       character(len = *), intent(in), optional :: speciesFilename

     end subroutine loadRegionData

  end interface

  interface

     subroutine saveRegionData(this, quantityOfInterest, filename)

       import :: t_Region

       class(t_Region) :: this
       integer, intent(in) :: quantityOfInterest
       character(len = *), intent(in) :: filename

     end subroutine saveRegionData

  end interface

  interface

     function getCfl(this) result(cfl)

       import :: t_Region

       class(t_Region) :: this

       real(SCALAR_KIND) :: cfl

     end function getCfl

  end interface

  interface

     function getTimeStepSize(this) result(timeStepSize)

       import :: t_Region

       class(t_Region) :: this

       real(SCALAR_KIND) :: timeStepSize

     end function getTimeStepSize

  end interface

  interface

     subroutine reportGridDiagnostics(this)

       import :: t_Region

       class(t_Region) :: this

     end subroutine reportGridDiagnostics

  end interface

  interface

     subroutine computeRhs(this, mode)

       import :: t_Region

       class(t_Region) :: this
       integer, intent(in) :: mode

     end subroutine computeRhs

  end interface

  interface

     subroutine saveSpongeStrength(this, filename)

       import :: t_Region

       class(t_Region) :: this
       character(len = *), intent(in) :: filename

     end subroutine saveSpongeStrength

  end interface

  interface

     subroutine resetProbes(this)

       import :: t_Region

       class(t_Region) :: this

     end subroutine resetProbes

  end interface

  interface

     subroutine saveProbeData(this, mode, finish)

       import :: t_Region

       class(t_Region) :: this
       integer, intent(in) :: mode
       logical, intent(in), optional :: finish

     end subroutine saveProbeData

  end interface

end module Region_mod
