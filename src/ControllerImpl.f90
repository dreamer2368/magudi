#include "config.h"

subroutine setupController(this, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Controller_mod, only : t_Controller
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_Controller) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

end subroutine setupController

subroutine cleanupController(this)

  ! <<< Derived types >>>
  use Controller_mod, only : t_Controller

  implicit none

  ! <<< Arguments >>>
  class(t_Controller) :: this

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  this%cachedValue = 0.0_wp
  this%runningTimeQuadrature = 0.0_wp
  this%onsetTime = 0.0_wp
  this%duration = 0.0_wp

end subroutine cleanupController

subroutine writeSensitivityToFile(this, comm, filename, timestep, time, append)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Controller_mod, only : t_Controller

  ! <<< Internal modules >>>
  use InputHelper, only : getFreeUnit
  use ErrorHandler, only : gracefulExit

  ! <<< Arguments >>>
  class(t_Controller) :: this
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

end subroutine writeSensitivityToFile
