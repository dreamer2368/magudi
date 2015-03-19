#include "config.h"

module Region_type

  use MPI, only : MPI_COMM_NULL

  use Grid_type
  use State_type
  use Patch_type
  use SolverOptions_type
  use PatchDescriptor_type
  use SimulationFlags_type

  implicit none
  private

  integer, parameter, public ::                                                              &
       FORWARD = +1,                                                                         &
       ADJOINT = -1

  type, public :: t_Region

     type(t_Grid), allocatable :: grids(:)
     type(t_State), allocatable :: states(:)
     type(t_Patch), allocatable :: patches(:)
     type(t_SolverOptions) :: solverOptions
     type(t_SimulationFlags) :: simulationFlags
     type(t_PatchDescriptor), allocatable :: patchData(:)
     integer :: comm = MPI_COMM_NULL
     integer, allocatable :: globalGridSizes(:,:), processDistributions(:,:),                &
          gridCommunicators(:), patchCommunicators(:)

  end type t_Region

end module Region_type

module Region_mod

  implicit none
  public

  interface

     subroutine setupRegion(this, comm, globalGridSizes, boundaryConditionFilename)

       use Region_type

       type(t_Region) :: this
       integer, intent(in) :: comm, globalGridSizes(:,:)

       character(len = *), intent(in), optional :: boundaryConditionFilename

     end subroutine setupRegion

  end interface

  interface

     subroutine cleanupRegion(this)

       use Region_type

       type(t_Region) :: this

     end subroutine cleanupRegion

  end interface

  interface

     subroutine loadRegionData(this, quantityOfInterest, filename)

       use Region_type

       type(t_Region) :: this
       integer, intent(in) :: quantityOfInterest
       character(len = *), intent(in) :: filename

     end subroutine loadRegionData

  end interface

  interface

     subroutine saveRegionData(this, quantityOfInterest, filename)

       use Region_type

       type(t_Region) :: this
       integer, intent(in) :: quantityOfInterest
       character(len = *), intent(in) :: filename

     end subroutine saveRegionData

  end interface

  interface

     subroutine reportGridDiagnostics(this)

       use Region_type

       type(t_Region) :: this

     end subroutine reportGridDiagnostics

  end interface

  interface

     subroutine computeRhs(this, mode, time)

       use Region_type

       type(t_Region) :: this
       integer, intent(in) :: mode
       real(SCALAR_KIND), intent(in) :: time

     end subroutine computeRhs

  end interface

  interface

     subroutine subStepHooks(this, mode, timestep, stage)

       use Region_type

       type(t_Region) :: this
       integer, intent(in) :: mode, timestep, stage

     end subroutine subStepHooks

  end interface

  interface

     subroutine reportResiduals(this)

       use Region_type

       type(t_Region) :: this

     end subroutine reportResiduals

  end interface

end module Region_mod
