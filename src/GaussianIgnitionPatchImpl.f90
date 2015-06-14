#include "config.h"

subroutine setupGaussianIgnitionPatch(this, index, comm, patchDescriptor,                &
     grid, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI
  use iso_fortran_env, only : real64

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags
  use GaussianIgnitionPatch_mod, only : t_GaussianIgnitionPatch

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_GaussianIgnitionPatch) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(real64), parameter :: pi = 4.0_real64 * atan(1.0_real64)
  real(SCALAR_KIND) :: power, heatRelease, referenceTemperature, flameTemperature
  character(len = 3), parameter :: directions = "xyz"
  character(len = STRING_LENGTH) :: key
  integer :: i, j, k, gridIndex, patchIndex, nDimensions

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (2, 3))

  call this%cleanup()
  call this%setupBase(index, comm, patchDescriptor, grid, simulationFlags, solverOptions)

  if (this%nPatchPoints > 0) then

  end if

  write(key, '(A)') "patches/" // trim(patchDescriptor%name) // "/"

  call getRequiredOption(trim(key) // "amplitude", this%amplitude, this%comm)
  call getRequiredOption(trim(key) // "radius", this%radius, this%comm)
  call getRequiredOption(trim(key) // "time_start", this%timeStart, this%comm)
  call getRequiredOption(trim(key) // "time_duration", this%timeDuration, this%comm)
  call getRequiredOption("heat_release", heatRelease, this%comm)

  this%origin = 0.0_wp

  do i = 1, 2
     this%origin(i) = getOption(trim(key) // "origin_" // directions(i:i), 0.0_wp)
  end do

  if (this%nPatchPoints > 0) then

     referenceTemperature = 1.0_wp / (solverOptions%ratioOfSpecificHeats - 1.0_wp)
     flameTemperature = referenceTemperature / (1.0_wp - heatRelease)
     power = 0.5_wp * this%amplitude * heatRelease * flameTemperature /                      &
          this%timeDuration/sqrt(2.0_wp*pi)

     allocate(this%strength(this%nPatchPoints))
     this%strength = 0.0_wp

     do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
        do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
           do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
              gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                   &
                   (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                     &
                   (k - 1 - this%gridOffset(3)))
              if (grid%iblank(gridIndex) == 0) cycle
              patchIndex = i - this%offset(1) + this%localSize(1) *                          &
                   (j - 1 - this%offset(2) + this%localSize(2) *                             &
                   (k - 1 - this%offset(3)))
              this%strength(patchIndex) = power * exp( - 0.5_wp *                            &
                   (((grid%coordinates(gridIndex,1) - this%origin(1))/this%radius) ** 2 +    &
                   ((grid%coordinates(gridIndex,2) - this%origin(2))/this%radius) ** 2) )
           end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
        end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
     end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  end if

end subroutine setupGaussianIgnitionPatch

subroutine cleanupGaussianIgnitionPatch(this)

  ! <<< Derived types >>>
  use GaussianIgnitionPatch_mod, only : t_GaussianIgnitionPatch

  implicit none

  ! <<< Arguments >>>
  class(t_GaussianIgnitionPatch) :: this

  call this%cleanupBase()

  SAFE_DEALLOCATE(this%strength)

end subroutine cleanupGaussianIgnitionPatch

subroutine addGaussianIgnition(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use GaussianIgnitionPatch_mod, only : t_GaussianIgnitionPatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_GaussianIgnitionPatch) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, gridIndex, patchIndex
  real(SCALAR_KIND) :: timePortion

  assert_key(mode, (FORWARD, ADJOINT))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  if (mode == ADJOINT) return

  call startTiming("addGaussianIgnition")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (2, 3))
  
  timePortion = exp( -0.5_wp * (state%time - this%timeStart) **2 / this%timeDuration **2 )

  do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
     do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
        do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
           gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                      &
                (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                        &
                (k - 1 - this%gridOffset(3)))
           if (grid%iblank(gridIndex) == 0) cycle
           patchIndex = i - this%offset(1) + this%localSize(1) *                             &
                (j - 1 - this%offset(2) + this%localSize(2) *                                &
                (k - 1 - this%offset(3)))

           state%rightHandSide(gridIndex, nDimensions + 2) =                                 &
                state%rightHandSide(gridIndex, ndimensions + 2) +                            &
                this%strength(patchIndex) * timePortion

        end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  call endTiming("addGaussianIgnition")

end subroutine addGaussianIgnition

function verifyGaussianIgnitionPatchUsage(this, patchDescriptor, gridSize,               &
     normalDirection, extent, simulationFlags,                                               &
     success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags
  use GaussianIgnitionPatch_mod, only : t_GaussianIgnitionPatch

  implicit none

  ! <<< Arguments >>>
  class(t_GaussianIgnitionPatch) :: this
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  integer, intent(in) :: gridSize(:), normalDirection, extent(6)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  logical, intent(out) :: success
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Result >>>
  logical :: isPatchUsed

  ! <<< Local variables >>>
  integer :: i, n

  isPatchUsed = .false.

  success = .false.

  n = size(gridSize)

  assert_key(n, (1, 2, 3))

  if (n == 1) then
     write(message, '(A)') "Can't be used with a 1D grid!"
     return
  end if

  do i = 1, size(gridSize)
     if (extent((i-1)*2+1) < 0 .or. extent((i-1)*2+2) > gridSize(i) .or.                     &
          extent((i-1)*2+1) > extent((i-1)*2+2)) then
        write(message, '(A)') "Invalid extent!"
        return
     end if
     if (extent((i-1)*2+1) == extent((i-1)*2+2)) n = n - 1
  end do

  if (n /= size(gridSize)) then
     write(message, '(2(A,I0.0),A)') "Expected a ", size(gridSize),                          &
          "D patch, but extent represents a ", n, "D patch!"
     return
  end if

  success = .true.

  isPatchUsed = .true.

end function verifyGaussianIgnitionPatchUsage
