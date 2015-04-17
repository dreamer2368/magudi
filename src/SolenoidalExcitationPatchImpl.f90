#include "config.h"

subroutine setupSolenoidalExcitationPatch(this, index, comm, patchDescriptor,                &
     grid, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI
  use iso_fortran_env, only : real64

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags
  use SolenoidalExcitationPatch_mod, only : t_SolenoidalExcitationPatch

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use RandomNumber, only : initializeRandomNumberGenerator

  implicit none

  ! <<< Arguments >>>
  class(t_SolenoidalExcitationPatch) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(real64), parameter :: pi = 4.0_real64 * atan(1.0_real64)
  character(len = 3), parameter :: directions = "xyz"
  character(len = STRING_LENGTH) :: key
  integer :: i, j, k, gridIndex, patchIndex, nDimensions, n, seed, ierror
  integer, allocatable :: seed_(:)

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (2, 3))

  call this%cleanup()
  call this%setupBase(index, comm, patchDescriptor, grid, simulationFlags, solverOptions)

  if (this%nPatchPoints > 0) then

  end if

  write(key, '(A)') "patches/" // trim(patchDescriptor%name) // "/"

  call getRequiredOption(trim(key) // "number_of_modes", this%nModes)
  this%nModes = max(0, this%nModes)

  this%origin = 0.0_wp
  this%speed = 0.0_wp

  do i = 1, 2
     this%origin(i) = getOption(trim(key) // "origin_" // directions(i:i), 0.0_wp)
     this%speed(i) = getOption(trim(key) // "velocity_" // directions(i:i), 0.0_wp)
  end do

  call getRequiredOption(trim(key) // "amplitude", this%amplitude)
  call getRequiredOption(trim(key) // "frequency", this%frequency)
  call getRequiredOption(trim(key) // "radius", this%gaussianFactor)
  this%gaussianFactor = 9.0_wp / (2.0_wp * this%gaussianFactor ** 2)

  if (this%nModes > 0) then

     seed = getOption(trim(key) // "random_seed", -1)
     call random_seed(size = n)
     allocate(seed_(n))
     seed_ = seed
     call random_seed(put = seed_)
     SAFE_DEALLOCATE(seed_)

     allocate(this%angularFrequencies(this%nModes))
     allocate(this%phases(this%nModes, 3))
     call random_number(this%angularFrequencies)
     call random_number(this%phases)
     this%phases = 2.0_real64 * pi * this%phases

     do i = 1, this%nModes
        this%angularFrequencies(i) =                                                         &
             2.0_real64 * pi * real(this%frequency, real64) *                                &
             (real(i, real64) + (this%angularFrequencies(i) - 0.5_real64)) /                 &
             (0.5_real64 * real(this%nModes, real64))
     end do
     call MPI_Bcast(this%angularFrequencies, size(this%angularFrequencies),                  &
          MPI_REAL8, 0, this%comm, ierror)
     call MPI_Bcast(this%phases, size(this%phases), MPI_REAL8, 0, this%comm, ierror)

  end if

  if (this%nPatchPoints > 0) then

     allocate(this%strength(this%nPatchPoints), source = 0.0_wp)

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
              this%strength(patchIndex) = this%amplitude * exp(- this%gaussianFactor *       &
                   sum((grid%coordinates(gridIndex,1:2) - this%origin) ** 2))
           end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
        end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
     end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  end if

end subroutine setupSolenoidalExcitationPatch

subroutine cleanupSolenoidalExcitationPatch(this)

  ! <<< Derived types >>>
  use SolenoidalExcitationPatch_mod, only : t_SolenoidalExcitationPatch

  implicit none

  ! <<< Arguments >>>
  class(t_SolenoidalExcitationPatch) :: this

  call this%cleanupBase()

  SAFE_DEALLOCATE(this%angularFrequencies)
  SAFE_DEALLOCATE(this%phases)
  SAFE_DEALLOCATE(this%strength)

end subroutine cleanupSolenoidalExcitationPatch

subroutine addSolenoidalExcitation(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use SolenoidalExcitationPatch_mod, only : t_SolenoidalExcitationPatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_SolenoidalExcitationPatch) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, nDimensions, gridIndex, patchIndex
  real(SCALAR_KIND) :: excitationAmount, phaseAngles(2), temp(4)

  assert_key(mode, (FORWARD, ADJOINT))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  if (mode == ADJOINT) return

  call startTiming("addSolenoidalExcitation")

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (2, 3))

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
           
           excitationAmount = this%strength(patchIndex)
           phaseAngles = grid%coordinates(gridIndex,1:2) -                                   &
                this%origin - this%speed * state%time

           do l = 1, this%nModes
              temp(1) = sin(this%angularFrequencies(l) * phaseAngles(1) + this%phases(l,1))
              temp(2) = cos(this%angularFrequencies(l) * phaseAngles(1) + this%phases(l,1))
              temp(3) = sin(this%angularFrequencies(l) * phaseAngles(2) + this%phases(l,2))
              temp(4) = cos(this%angularFrequencies(l) * phaseAngles(2) + this%phases(l,2))
              state%rightHandSide(gridIndex, 2) = state%rightHandSide(gridIndex, 2) +        &
                   excitationAmount * temp(1) * (this%angularFrequencies(l) * temp(4)        &
                   - 2.0_wp * this%gaussianFactor * (grid%coordinates(gridIndex, 2) -        &
                   this%origin(2)) * temp(3))
              state%rightHandSide(gridIndex, 3) = state%rightHandSide(gridIndex, 3) +        &
                   excitationAmount * temp(3) * (this%angularFrequencies(l) * temp(2) -      &
                   2.0_wp * this%gaussianFactor * (grid%coordinates(gridIndex, 1) -          &
                   this%origin(1)) * temp(1))
           end do

        end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  call endTiming("addSolenoidalExcitation")

end subroutine addSolenoidalExcitation

function verifySolenoidalExcitationPatchUsage(this, patchDescriptor, gridSize,               &
     normalDirection, extent, simulationFlags,                                               &
     success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags
  use SolenoidalExcitationPatch_mod, only : t_SolenoidalExcitationPatch

  implicit none

  ! <<< Arguments >>>
  class(t_SolenoidalExcitationPatch) :: this
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

end function verifySolenoidalExcitationPatchUsage

subroutine updateSolenoidalExcitationPatch(this, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use SolenoidalExcitationPatch_mod, only : t_SolenoidalExcitationPatch

  implicit none

  ! <<< Arguments >>>
  class(t_SolenoidalExcitationPatch) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state

  ! Nothing to be done for this patch type...

end subroutine updateSolenoidalExcitationPatch
