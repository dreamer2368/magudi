#include "config.h"

module SolenoidalExcitationPatch_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  use Patch_mod, only : t_Patch

  implicit none
  private

  integer, parameter, private :: real64 = selected_real_kind(15)

  type, extends(t_Patch), public :: t_SolenoidalExcitationPatch

     integer :: nModes
     real(SCALAR_KIND) :: origin(2), speed(2), amplitude,                                    &
          gaussianFactor, frequency

     real(real64), allocatable :: angularFrequencies(:), phases(:,:)

     real(SCALAR_KIND), allocatable :: strength(:), spatialFunctionsCache(:,:,:)

   contains

     procedure, pass :: setup
     procedure, pass :: cleanup
     procedure, pass :: updateRhs

  end type t_SolenoidalExcitationPatch

contains

  subroutine setup(this, name, comm, grid, state, extent,                                    &
       normalDirection, simulationFlags, solverOptions)

    ! <<< External modules >>>
    use MPI, only : MPI_REAL8, MPI_COMM_NULL

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod, only : t_State
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Internal modules >>>
    use InputHelper, only : getOption, getRequiredOption
    use ErrorHandler, only : gracefulExit, issueWarning

    implicit none

    ! <<< Arguments >>>
    class(t_SolenoidalExcitationPatch) :: this
    character(len = *), intent(in) :: name
    integer, intent(in) :: comm
    class(t_Grid), intent(in) :: grid
    class(t_State), intent(in) :: state
    integer, intent(in) :: extent(6), normalDirection
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(real64), parameter :: pi = 4.0_real64 * atan(1.0_real64)
    character(len = STRING_LENGTH) :: key, message
    integer :: i, j, k, l, gridIndex, patchIndex, n, seed, ierror,                           &
         nDimensions, nUnknowns, direction
    character(len = 3), parameter :: directions = "xyz"
    integer, allocatable :: seed_(:)

    call this%cleanup()
    call this%setupBase(name, comm, grid, extent, normalDirection)

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns >= nDimensions + 2)

    direction = abs(this%normalDirection)

    if (direction /= 0) then
       write(message, '(3A)') "Normal direction for patch '", trim(name), "' will be ignored!"
       call issueWarning(grid%comm, message)
    end if

    i = 0
    do while (extent(i*2+2) > extent(i*2+1))
       i = i + 1
    end do
    if (i /= nDimensions) then
       write(message, '(A,I0.0,3A,I0.0,A)') "Expected a ", nDimensions, "D patch, but '",    &
            trim(name), "' is a ", i, "D patch!"
       call gracefulExit(grid%comm, message)
    end if

    if (nDimensions == 1) then
       write(message, '(3A)') "Patch '", trim(name), "' can't be used with a 1D grid!"
       call gracefulExit(grid%comm, message)
    end if

    write(key, '(A)') "patches/" // trim(name) // "/"

    call getRequiredOption(trim(key) // "number_of_modes", this%nModes, this%comm)
    this%nModes = max(0, this%nModes)

    this%origin = 0.0_wp
    this%speed = 0.0_wp

    do i = 1, 2
       this%origin(i) = getOption(trim(key) // "origin_" // directions(i:i), 0.0_wp)
       this%speed(i) = getOption(trim(key) // "velocity_" // directions(i:i), 0.0_wp)
    end do

    call getRequiredOption(trim(key) // "amplitude", this%amplitude, this%comm)
    call getRequiredOption(trim(key) // "frequency", this%frequency, this%comm)
    call getRequiredOption(trim(key) // "radius", this%gaussianFactor, this%comm)
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
       do i = 1, 3
          call random_number(this%phases(:,i))
       end do
       call random_number(this%angularFrequencies)
       this%phases = 2.0_real64 * pi * this%phases

       do i = 1, this%nModes
          this%angularFrequencies(i) =                                                       &
               2.0_real64 * pi * real(this%frequency, real64) *                              &
               (real(i, real64) + (this%angularFrequencies(i) - 0.5_real64)) /               &
               (0.5_real64 * real(this%nModes, real64))
       end do
       call MPI_Bcast(this%angularFrequencies, size(this%angularFrequencies),                &
            MPI_REAL8, 0, grid%comm, ierror)
       call MPI_Bcast(this%phases, size(this%phases), MPI_REAL8, 0, grid%comm, ierror)

    end if

    if (this%nPatchPoints > 0) then

       allocate(this%strength(this%nPatchPoints))
       this%strength = 0.0_wp

       if (this%nModes > 0)                                                                  &
            allocate(this%spatialFunctionsCache(this%nPatchPoints, this%nModes, 4))

       do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
          do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
             do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                 &
                     (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                   &
                     (k - 1 - this%gridOffset(3)))
                if (grid%iblank(gridIndex) == 0) cycle
                patchIndex = i - this%offset(1) + this%localSize(1) *                        &
                     (j - 1 - this%offset(2) + this%localSize(2) *                           &
                     (k - 1 - this%offset(3)))

                this%strength(patchIndex) = this%amplitude * exp(- this%gaussianFactor *     &
                     sum((grid%coordinates(gridIndex,1:2) - this%origin) ** 2))

                do l = 1, this%nModes
                   this%spatialFunctionsCache(patchIndex,l,1) =                              &
                        sin(this%angularFrequencies(l) * (grid%coordinates(gridIndex,1) -    &
                        this%origin(1)) + this%phases(l,1))
                   this%spatialFunctionsCache(patchIndex,l,2) =                              &
                        cos(this%angularFrequencies(l) * (grid%coordinates(gridIndex,1) -    &
                        this%origin(1)) + this%phases(l,1))
                   this%spatialFunctionsCache(patchIndex,l,3) =                              &
                        sin(this%angularFrequencies(l) * (grid%coordinates(gridIndex,2) -    &
                        this%origin(2)) + this%phases(l,2))
                   this%spatialFunctionsCache(patchIndex,l,4) =                              &
                        cos(this%angularFrequencies(l) * (grid%coordinates(gridIndex,2) -    &
                        this%origin(2)) + this%phases(l,2))
                end do

             end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
          end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
       end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

    end if

  end subroutine setup

  subroutine cleanup(this)

    implicit none

    ! <<< Arguments >>>
    class(t_SolenoidalExcitationPatch) :: this

    call this%cleanupBase()

    SAFE_DEALLOCATE(this%angularFrequencies)
    SAFE_DEALLOCATE(this%phases)
    SAFE_DEALLOCATE(this%strength)
    SAFE_DEALLOCATE(this%spatialFunctionsCache)

  end subroutine cleanup

  subroutine updateRhs(this, mode, simulationFlags, solverOptions, grid, state)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod, only : t_State
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : FORWARD, ADJOINT

    ! <<< Internal modules >>>
    use CNSHelper
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
    integer :: i, j, k, gridIndex, patchIndex, nDimensions
    real(SCALAR_KIND) :: excitationAmount
    real(SCALAR_KIND), allocatable :: temporalFunctions(:,:), temp(:,:,:)

    assert_key(mode, (FORWARD, ADJOINT))

    if (mode == ADJOINT) return

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (2, 3))

    call startTiming("addSolenoidalExcitation")

    if (this%nModes > 0) then
       allocate(temporalFunctions(this%nModes, 4))
       do i = 1, this%nModes
          temporalFunctions(i,1) = sin(this%angularFrequencies(i) *                          &
               this%speed(1) * state%time)
          temporalFunctions(i,2) = cos(this%angularFrequencies(i) *                          &
               this%speed(1) * state%time)
          temporalFunctions(i,3) = sin(this%angularFrequencies(i) *                          &
               this%speed(2) * state%time)
          temporalFunctions(i,4) = cos(this%angularFrequencies(i) *                          &
               this%speed(2) * state%time)
       end do
    end if

    if (this%nPatchPoints > 0 .and. this%nModes > 0) then
       allocate(temp(this%nPatchPoints, this%nModes, 4))
       do i = 1, this%nModes
          temp(:,i,1) = this%spatialFunctionsCache(:,i,1) * temporalFunctions(i,2) -         &
               this%spatialFunctionsCache(:,i,2) * temporalFunctions(i,1)
          temp(:,i,2) = this%spatialFunctionsCache(:,i,2) * temporalFunctions(i,2) +         &
               this%spatialFunctionsCache(:,i,1) * temporalFunctions(i,1)
          temp(:,i,3) = this%spatialFunctionsCache(:,i,3) * temporalFunctions(i,4) -         &
               this%spatialFunctionsCache(:,i,4) * temporalFunctions(i,3)
          temp(:,i,4) = this%spatialFunctionsCache(:,i,4) * temporalFunctions(i,4) +         &
               this%spatialFunctionsCache(:,i,3) * temporalFunctions(i,3)
       end do
    end if

    do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
       do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
          do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
             gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                    &
                  (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                      &
                  (k - 1 - this%gridOffset(3)))
             if (grid%iblank(gridIndex) == 0) cycle
             patchIndex = i - this%offset(1) + this%localSize(1) *                           &
                  (j - 1 - this%offset(2) + this%localSize(2) *                              &
                  (k - 1 - this%offset(3)))

             excitationAmount = this%strength(patchIndex)

             state%rightHandSide(gridIndex, 2) = state%rightHandSide(gridIndex, 2) +         &
                  excitationAmount * sum(temp(patchIndex,:,1) * (this%angularFrequencies *   &
                  temp(patchIndex,:,4) - 2.0_wp * this%gaussianFactor *                      &
                  (grid%coordinates(gridIndex, 2) - this%origin(2)) * temp(patchIndex,:,3)))
             state%rightHandSide(gridIndex, 3) = state%rightHandSide(gridIndex, 3) +         &
                  excitationAmount * sum(temp(patchIndex,:,3) * (this%angularFrequencies *   &
                  temp(patchIndex,:,2) - 2.0_wp * this%gaussianFactor *                      &
                  (grid%coordinates(gridIndex, 1) - this%origin(1)) * temp(patchIndex,:,1)))

          end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
       end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
    end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

    SAFE_DEALLOCATE(temp)
    SAFE_DEALLOCATE(temporalFunctions)

    call endTiming("addSolenoidalExcitation")

  end subroutine updateRhs

end module SolenoidalExcitationPatch_mod
