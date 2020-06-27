#include "config.h"

subroutine setupKolmogorovForcingPatch(this, index, comm, patchDescriptor,                &
     grid, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI
  use iso_fortran_env, only : real64

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags
  use KolmogorovForcingPatch_mod, only : t_KolmogorovForcingPatch

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_KolmogorovForcingPatch) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(real64), parameter :: pi = 4.0_real64 * atan(1.0_real64)
  character(len = STRING_LENGTH) :: key
  integer :: i, j, k, gridIndex, patchIndex, nDimensions, ierror

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (2, 3))

  call this%cleanup()
  call this%setupBase(index, comm, patchDescriptor, grid, simulationFlags, solverOptions)

  write(key, '(A)') "patches/" // trim(patchDescriptor%name) // "/"

  call getRequiredOption(trim(key) // "wavenumber", this%n, this%comm)
  this%n = max(0, this%n)

  this%amplitude = 0.0_wp

  call getRequiredOption(trim(key) // "amplitude", this%amplitude, this%comm)

  if (this%nPatchPoints > 0) then

     allocate(this%forcePerUnitMass(this%nPatchPoints))
     this%forcePerUnitMass = 0.0_wp

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

              this%forcePerUnitMass(patchIndex) = this%amplitude *                           &
                       sin(2.0_wp * pi * this%n * grid%coordinates(gridIndex,2))

           end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
        end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
     end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  end if

end subroutine setupKolmogorovForcingPatch

subroutine cleanupKolmogorovForcingPatch(this)

  ! <<< Derived types >>>
  use KolmogorovForcingPatch_mod, only : t_KolmogorovForcingPatch

  implicit none

  ! <<< Arguments >>>
  class(t_KolmogorovForcingPatch) :: this

  call this%cleanupBase()

  SAFE_DEALLOCATE(this%forcePerUnitMass)

end subroutine cleanupKolmogorovForcingPatch

subroutine addKolmogorovForcing(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use KolmogorovForcingPatch_mod, only : t_KolmogorovForcingPatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT, LINEARIZED

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_KolmogorovForcingPatch) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, gridIndex, patchIndex
  real(SCALAR_KIND) :: time, excitationAmount
  real(SCALAR_KIND), allocatable :: temporalFunctions(:,:)
  SCALAR_TYPE, allocatable :: temp(:,:,:)
  SCALAR_TYPE :: localForcing

  assert_key(mode, (FORWARD, ADJOINT, LINEARIZED))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  call startTiming("addKolmogorovForcing")

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

           localForcing = this%forcePerUnitMass(patchIndex)

           select case (mode)
           case (FORWARD)
              state%rightHandSide(gridIndex,2) = state%rightHandSide(gridIndex,2) +          &
                   state%conservedVariables(gridIndex,1) * localForcing
              ! state%rightHandSide(gridIndex,nDimensions+2) =                                 &
              !      state%rightHandSide(gridIndex,nDimensions+2) +                            &
              !      state%conservedVariables(gridIndex,2) * localForcing
           case (ADJOINT)
              state%rightHandSide(gridIndex,1) = state%rightHandSide(gridIndex,1) -          &
                   state%adjointVariables(gridIndex,2) * localForcing
              ! state%rightHandSide(gridIndex,2) = state%rightHandSide(gridIndex,2) -          &
              !      state%adjointVariables(gridIndex,nDimensions+2) * localForcing
           case (LINEARIZED)
              state%rightHandSide(gridIndex,2) = state%rightHandSide(gridIndex,2) +          &
                   state%adjointVariables(gridIndex,1) * localForcing
              ! state%rightHandSide(gridIndex,nDimensions+2) =                                 &
              !      state%rightHandSide(gridIndex,nDimensions+2) +                            &
              !      state%adjointVariables(gridIndex,2) * localForcing
           end select

        end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

  SAFE_DEALLOCATE(temp)
  SAFE_DEALLOCATE(temporalFunctions)

  call endTiming("addKolmogorovForcing")

end subroutine addKolmogorovForcing

function verifyKolmogorovForcingPatchUsage(this, patchDescriptor, gridSize,               &
     normalDirection, extent, simulationFlags,                                               &
     success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags
  use KolmogorovForcingPatch_mod, only : t_KolmogorovForcingPatch

  implicit none

  ! <<< Arguments >>>
  class(t_KolmogorovForcingPatch) :: this
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

end function verifyKolmogorovForcingPatchUsage
