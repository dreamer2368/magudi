#include "config.h"

subroutine setupResidualManager(this, prefix, region)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use ResidualManager_mod, only : t_ResidualManager

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_ResidualManager) :: this
  character(len = *), intent(in) :: prefix
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = 3), parameter :: directions = "xyz"
  integer :: i, nDimensions, nUnknowns
  character(len = STRING_LENGTH) :: message

  call this%cleanup()

  this%reportInterval = getOption("report_interval", 0)
  this%reportInterval = getOption("check_residuals_interval", this%reportInterval)

  if (this%reportInterval <= 0) return

  nDimensions = size(region%globalGridSizes, 1)
  assert_key(nDimensions, (1, 2, 3))

  nUnknowns = region%solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2)

  allocate(this%residuals(nUnknowns))
  allocate(this%tolerances(nUnknowns))

  this%residuals(:) = huge(0.0_wp)

  this%tolerances(:) = getOption("residuals/convergence_limit", sqrt(epsilon(0.0_wp)))

  do i = 1, size(this%tolerances)
     this%tolerances(i) = getOption(trim(prefix) // "/convergence_limit", this%tolerances(i))
  end do

  this%tolerances(1) = getOption("residuals/density_limit", this%tolerances(1))
  if (len_trim(prefix) > 0)                                                                  &
       this%tolerances(1) = getOption(trim(prefix) // "/density_limit", this%tolerances(1))

  do i = 1, nDimensions
     this%tolerances(i+1) = getOption("residuals/momentum_limit", this%tolerances(i+1))
     this%tolerances(i+1) = getOption("residuals/" // directions(i:i) //                     &
          "_momentum_limit", this%tolerances(i+1))
     if (len_trim(prefix) > 0)                                                               &
          this%tolerances(i+1) = getOption(trim(prefix) // "/" // directions(i:i) //         &
          "_momentum_limit", this%tolerances(i+1))
  end do

  this%tolerances(nDimensions + 2) =                                                         &
       getOption("residuals/energy_limit", this%tolerances(nDimensions + 2))
  if (len_trim(prefix) > 0) this%tolerances(nDimensions + 2) =                               &
       getOption(trim(prefix) // "/energy_limit", this%tolerances(nDimensions + 2))

  if (any(this%tolerances <= 0.0_wp)) then
     write(message, '(2A)') "Residual convergence limits must be positive.",                 &
          " Please check your input file."
     call gracefulExit(region%comm, message)
  end if

end subroutine setupResidualManager

subroutine cleanupResidualManager(this)

  ! <<< Derived types >>>
  use ResidualManager_mod, only : t_ResidualManager

  implicit none

  ! <<< Arguments >>>
  class(t_ResidualManager) :: this

  SAFE_DEALLOCATE(this%residuals)
  SAFE_DEALLOCATE(this%tolerances)

  this%hasSimulationConverged = .false.
  this%reportInterval = 0

end subroutine cleanupResidualManager

subroutine computeResiduals(this, region)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use ResidualManager_mod, only : t_ResidualManager

  implicit none

  ! <<< Arguments >>>
  class(t_ResidualManager) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions, ierror
  SCALAR_TYPE, allocatable :: F(:)
  SCALAR_TYPE :: maximumF

  if (.not. allocated(this%residuals) .or. this%reportInterval <= 0) return

  do i = 1, size(region%states)

     nDimensions = region%grids(i)%nDimensions
     assert_key(nDimensions, (1, 2, 3))

     allocate(F(region%grids(i)%nGridPoints))

     do j = 1, region%solverOptions%nUnknowns
        F = abs(region%states(i)%rightHandSide(:,j))
        call region%grids(i)%findMaximum(F, maximumF)
        this%residuals(j) = real(maximumF, wp)
     end do

     SAFE_DEALLOCATE(F)

  end do

  call MPI_Allreduce(MPI_IN_PLACE, this%residuals, size(this%residuals),                     &
       REAL_TYPE_MPI, MPI_MAX, region%comm, ierror)

  if (allocated(this%tolerances))                                                            &
       this%hasSimulationConverged = (all(this%residuals < this%tolerances))

end subroutine computeResiduals

subroutine writeResidualsToFile(this, comm, filename, timestep, time, append)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use ResidualManager_mod, only : t_ResidualManager

  ! <<< Internal modules >>>
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_ResidualManager) :: this
  integer, intent(in) :: comm
  character(len = *), intent(in) :: filename
  integer, intent(in) :: timestep
  real(SCALAR_KIND), intent(in) :: time
  logical, intent(in), optional :: append

  ! <<< Local variables >>>
  logical :: append_
  integer :: fileUnit, ostat, procRank, ierror
  character(len = STRING_LENGTH) :: formatString, message

  if (.not. allocated(this%residuals)) return

  append_ = .false.
  if (present(append)) append_ = append

  call MPI_Comm_rank(comm, procRank, ierror)

  if (procRank == 0) then
     if (.not. append_) then
        open(newunit = fileUnit, file = trim(filename), action = 'write',                    &
             status = 'unknown', iostat = ostat)
     else
        open(newunit = fileUnit, file = trim(filename), action = 'write',                    &
             status = 'old', position = 'append', iostat = ostat)
     end if
  end if

  call MPI_Bcast(ostat, 1, MPI_INTEGER, 0, comm, ierror)
  if (ostat /= 0) then
     write(message, "(2A)") trim(filename), ": Failed to open file for writing!"
     call gracefulExit(comm, message)
  end if

  if (procRank == 0) then
     write(formatString, '(A,I0.0,A)') "(I8,1X,E13.6,", size(this%residuals),                &
          "(1X,SP," // SCALAR_FORMAT // "))"
     write(fileUnit, trim(formatString)) timestep, time, this%residuals
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

end subroutine writeResidualsToFile
