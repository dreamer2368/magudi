#include "config.h"

module Region_enum

  implicit none
  public

  integer, parameter ::                                                                      &
       FORWARD = +1,                                                                         &
       ADJOINT = -1

end module Region_enum

module Region_mod

  use MPI, only : MPI_COMM_NULL

  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use Patch_type, only : t_Patch
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_type, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none
  private

  type, public :: t_Region

     type(t_Grid), allocatable :: grids(:)
     type(t_State), allocatable :: states(:)
     type(t_Patch), allocatable :: patches(:)
     type(t_SolverOptions) :: solverOptions
     type(t_SimulationFlags) :: simulationFlags
     type(t_PatchDescriptor), allocatable :: patchData(:)
     integer :: comm = MPI_COMM_NULL, commGridMasters = MPI_COMM_NULL
     integer, allocatable :: globalGridSizes(:,:), processDistributions(:,:),                &
          gridCommunicators(:), patchCommunicators(:)

   contains
     
     procedure, pass :: setup => setupRegion
     procedure, pass :: cleanup => cleanupRegion
     procedure, pass :: loadData => loadRegionData
     procedure, pass :: saveData => saveRegionData
     procedure, pass :: getCfl
     procedure, pass :: getTimeStepSize
     procedure, pass :: reportGridDiagnostics
     procedure, pass :: computeRhs
     procedure, pass :: computeResiduals

  end type t_Region

  interface

     subroutine setupRegion(this, comm, globalGridSizes, boundaryConditionFilename)

       import :: t_Region

       class(t_Region) :: this
       integer, intent(in) :: comm, globalGridSizes(:,:)

       character(len = *), intent(in), optional :: boundaryConditionFilename

     end subroutine setupRegion

  end interface

  interface

     subroutine cleanupRegion(this)

       import :: t_Region

       class(t_Region) :: this

     end subroutine cleanupRegion

  end interface

  interface

     subroutine loadRegionData(this, quantityOfInterest, filename)

       import :: t_Region

       class(t_Region) :: this
       integer, intent(in) :: quantityOfInterest
       character(len = *), intent(in) :: filename

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

     subroutine computeRhs(this, mode, time)

       import :: t_Region

       class(t_Region) :: this
       integer, intent(in) :: mode
       real(SCALAR_KIND), intent(in) :: time

     end subroutine computeRhs

  end interface

  interface

     subroutine computeResiduals(this, residuals)

       import :: t_Region

       class(t_Region) :: this
       real(SCALAR_KIND), intent(out) :: residuals(3)

     end subroutine computeResiduals

  end interface

end module Region_mod
