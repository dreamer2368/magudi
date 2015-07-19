#include "config.h"

subroutine setupFunctional(this, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Functional_mod, only : t_Functional
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_Functional) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

end subroutine setupFunctional

subroutine cleanupFunctional(this)

  ! <<< Derived types >>>
  use Functional_mod, only : t_Functional

  implicit none

  ! <<< Arguments >>>
  class(t_Functional) :: this

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  this%cachedValue = 0.0_wp
  this%runningTimeQuadrature = 0.0_wp

end subroutine cleanupFunctional

subroutine updateAdjointForcing(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use Functional_mod, only : t_Functional
  use CostTargetPatch_mod, only : t_CostTargetPatch

  implicit none

  ! <<< Arguments >>>
  class(t_Functional) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer :: i, j, nDimensions
  class(t_Patch), pointer :: patch => null()
  class(t_CostTargetPatch), pointer :: costTargetPatch => null()

  if (.not. allocated(region%patchFactories)) return

  do i = 1, size(region%grids)

     nDimensions = region%grids(i)%nDimensions

     do j = 1, size(region%patchFactories)

        call region%patchFactories(j)%connect(patch)
        if (.not. associated(patch)) cycle

        nullify(costTargetPatch)
        select type (patch)
        type is (t_CostTargetPatch)
           costTargetPatch => patch
        end select

        if (.not. associated(costTargetPatch)) cycle

        if (costTargetPatch%gridIndex /= region%grids(i)%index .or.                          &
             costTargetPatch%nPatchPoints <= 0) cycle

        call this%computeAdjointForcing(region%simulationFlags, region%solverOptions,        &
             region%grids(i), region%states(i), costTargetPatch)

     end do !... j = 1, size(region%patchFactories)

  end do !... i = 1, size(region%grids)

end subroutine updateAdjointForcing

subroutine writeFunctionalToFile(this, comm, filename, timestep, time, append)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Functional_mod, only : t_Functional

  ! <<< Internal modules >>>
  use InputHelper, only : getFreeUnit
  use ErrorHandler, only : gracefulExit

  ! <<< Arguments >>>
  class(t_Functional) :: this
  integer, intent(in) :: comm
  character(len = *), intent(in) :: filename
  integer, intent(in) :: timestep
  real(SCALAR_KIND), intent(in) :: time
  logical, intent(in), optional :: append

  ! <<< Local variables >>>
  logical :: append_
  integer :: fileUnit, ostat, procRank, ierror
  character(len = STRING_LENGTH) :: message

  append_ = .false.
  if (present(append)) append_ = append

  call MPI_Comm_rank(comm, procRank, ierror)

  if (procRank == 0) then
     if (.not. append_) then
        open(unit = getFreeUnit(fileUnit), file = trim(filename), action = 'write',          &
             status = 'unknown', iostat = ostat)
     else
        open(unit = getFreeUnit(fileUnit), file = trim(filename), action = 'write',          &
             status = 'old', position = 'append', iostat = ostat)
     end if
  end if

  call MPI_Bcast(ostat, 1, MPI_INTEGER, 0, comm, ierror)
  if (ostat /= 0) then
     write(message, "(2A)") trim(filename), ": Failed to open file for writing!"
     call gracefulExit(comm, message)
  end if

  if (procRank == 0) then
     write(fileUnit, '(I8,1X,E13.6,2(1X,SP,' // SCALAR_FORMAT // '))')                       &
          timestep, time, this%cachedValue, this%runningTimeQuadrature
  end if

  call MPI_Bcast(ostat, 1, MPI_INTEGER, 0, comm, ierror)
  if (ostat /= 0) then
     write(message, "(2A)") trim(filename), ": Error writing to file!"
     call gracefulExit(comm, message)
  end if

  if (procRank == 0) then
     flush(fileUnit)
     close(fileUnit)
  end if

  call MPI_Barrier(comm, ierror)

end subroutine writeFunctionalToFile
