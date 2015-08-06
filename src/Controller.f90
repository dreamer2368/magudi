#include "config.h"

module Controller_mod

  implicit none

  type, abstract, public :: t_Controller

     real(SCALAR_KIND) :: cachedValue = real(0.0, SCALAR_KIND),                              &
          runningTimeQuadrature = real(0.0, SCALAR_KIND),                                    &
          onsetTime = real(0.0, SCALAR_KIND), duration = real(0.0, SCALAR_KIND)

   contains

     procedure, non_overridable, pass :: cleanupBase
     procedure, pass :: writeToFile

     procedure(setup), pass, deferred :: setup
     procedure(cleanup), pass, deferred :: cleanup
     procedure(computeSensitivity), pass, deferred :: computeSensitivity
     procedure(update), pass, deferred :: update
     procedure(updateGradient), pass, deferred :: updateGradient
     procedure(hookBeforeTimemarch), pass, deferred :: hookBeforeTimemarch
     procedure(hookAfterTimemarch), pass, deferred :: hookAfterTimemarch

  end type t_Controller

  abstract interface

     subroutine setup(this, region)

       use Region_mod, only : t_Region

       import :: t_Controller

       class(t_Controller) :: this
       class(t_Region) :: region

     end subroutine setup

  end interface

  abstract interface

     subroutine cleanup(this)

       import :: t_Controller

       class(t_Controller) :: this

     end subroutine cleanup

  end interface

  abstract interface

     function computeSensitivity(this, region) result(instantaneousSensitivity)

       use Region_mod, only : t_Region

       import :: t_Controller

       class(t_Controller) :: this
       class(t_Region), intent(in) :: region

       real(SCALAR_KIND) :: instantaneousSensitivity

     end function computeSensitivity

  end interface

  abstract interface

     subroutine update(this, region)

       use Region_mod, only : t_Region

       import :: t_Controller

       class(t_Controller) :: this
       class(t_Region), intent(in) :: region

     end subroutine update

  end interface

  abstract interface

     subroutine updateGradient(this, region)

       use Region_mod, only : t_Region

       import :: t_Controller

       class(t_Controller) :: this
       class(t_Region), intent(in) :: region

     end subroutine updateGradient

  end interface

  abstract interface

     subroutine hookBeforeTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_Controller

       class(t_Controller) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookBeforeTimemarch

  end interface

  abstract interface

     subroutine hookAfterTimemarch(this, region, mode)

       use Region_mod, only : t_Region

       import :: t_Controller

       class(t_Controller) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine hookAfterTimemarch

  end interface

contains

  subroutine cleanupBase(this)

    implicit none

    ! <<< Arguments >>>
    class(t_Controller) :: this

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND

    this%cachedValue = 0.0_wp
    this%runningTimeQuadrature = 0.0_wp
    this%onsetTime = 0.0_wp
    this%duration = 0.0_wp

  end subroutine cleanupBase

  subroutine writeToFile(this, comm, filename, timestep, time, append)

    ! <<< External modules >>>
    use MPI

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
          open(unit = getFreeUnit(fileUnit), file = trim(filename), action = 'write',        &
               status = 'unknown', iostat = ostat)
       else
          open(unit = getFreeUnit(fileUnit), file = trim(filename), action = 'write',        &
               status = 'old', position = 'append', iostat = ostat)
       end if
    end if

    call MPI_Bcast(ostat, 1, MPI_INTEGER, 0, comm, ierror)
    if (ostat /= 0) then
       write(message, "(2A)") trim(filename), ": Failed to open file for writing!"
       call gracefulExit(comm, message)
    end if

    if (procRank == 0) then
       write(fileUnit, '(I8,1X,E13.6,2(1X,SP,' // SCALAR_FORMAT // '))')                     &
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

  end subroutine writeToFile

end module Controller_mod
