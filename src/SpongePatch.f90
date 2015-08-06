#include "config.h"

module SpongePatch_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  use Patch_mod, only : t_Patch

  implicit none
  private

  type, extends(t_Patch), public :: t_SpongePatch

     real(SCALAR_KIND) :: spongeAmount
     integer :: spongeExponent
     real(SCALAR_KIND), allocatable :: spongeStrength(:)

   contains

     procedure, pass :: setup
     procedure, pass :: cleanup
     procedure, pass :: updateRhs
     procedure, pass :: computeStrength

  end type t_SpongePatch

contains

  subroutine setup(this, name, comm, grid, state, extent,                                    &
       normalDirection, simulationFlags, solverOptions)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod, only : t_State
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Internal modules >>>
    use InputHelper, only : getOption
    use ErrorHandler, only : gracefulExit

    implicit none

    ! <<< Arguments >>>
    class(t_SpongePatch) :: this
    character(len = *), intent(in) :: name
    integer, intent(in) :: comm
    class(t_Grid), intent(in) :: grid
    class(t_State), intent(in) :: state
    integer, intent(in) :: extent(6), normalDirection
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    character(len = STRING_LENGTH) :: key, message
    integer :: nDimensions, nUnknowns, direction

    call this%cleanup()
    call this%setupBase(name, comm, grid, extent, normalDirection)

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns >= nDimensions + 2)

    direction = abs(this%normalDirection)

    if (direction <= 0 .or. direction > nDimensions) then
       write(message, '(3A)') "Invalid normal direction for patch '", trim(name), "'!"
       call gracefulExit(grid%comm, message)
    end if

    if (.not. simulationFlags%useTargetState) then
       write(message, '(3A)') "Patch '", trim(name), "' requires a target state for damping!"
       call gracefulExit(grid%comm, message)
    end if

    if (extent((direction-1)*2+1) == extent((direction-1)*2+2)) then
       write(message, '(A)')                                                                 &
            "Sponge must extend at least 2 grid points along the normal direction!"
       call gracefulExit(grid%comm, message)
    end if

    if ((normalDirection > 0 .and. extent((direction-1)*2+1) /= 1) .or.                      &
         (normalDirection < 0 .and.                                                          &
         extent((direction-1)*2+2) /= grid%globalSize(direction))) then
       write(message, '(3A)') "Sponge zone on patch '", trim(name),                          &
            "' is not aligned with a computational boundary!"
       call gracefulExit(grid%comm, message)
    end if

    write(key, '(A)') "patches/" // trim(name) // "/"

    ! Sponge amount.
    this%spongeAmount = getOption("defaults/sponge_amount", 1.0_wp)
    this%spongeAmount = getOption(trim(key) // "sponge_amount", this%spongeAmount)

    ! Sponge exponent.
    this%spongeExponent = getOption("defaults/sponge_exponent", 2)
    this%spongeExponent = getOption(trim(key) // "sponge_exponent", this%spongeExponent)

    if (this%nPatchPoints > 0) allocate(this%spongeStrength(this%nPatchPoints))

  end subroutine setup

  subroutine cleanup(this)

    implicit none

    ! <<< Arguments >>>
    class(t_SpongePatch) :: this

    call this%cleanupBase()

    SAFE_DEALLOCATE(this%spongeStrength)

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
    class(t_SpongePatch) :: this
    integer, intent(in) :: mode
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid), intent(in) :: grid
    class(t_State) :: state

    ! <<< Local variables >>>
    integer :: i, j, k, l, gridIndex, patchIndex

    assert_key(mode, (FORWARD, ADJOINT))

    call startTiming("addDamping")

    select case (mode)

    case (FORWARD)
       do l = 1, solverOptions%nUnknowns
          do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
             do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
                do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                   gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *              &
                        (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                &
                        (k - 1 - this%gridOffset(3)))
                   if (grid%iblank(gridIndex) == 0) cycle
                   patchIndex = i - this%offset(1) + this%localSize(1) *                     &
                        (j - 1 - this%offset(2) + this%localSize(2) *                        &
                        (k - 1 - this%offset(3)))
                   state%rightHandSide(gridIndex, l) = state%rightHandSide(gridIndex, l) -   &
                        this%spongeStrength(patchIndex) *                                    &
                        (state%conservedVariables(gridIndex, l) -                            &
                        state%targetState(gridIndex, l))
                end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
             end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
          end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
       end do !... l = 1, solverOptions%nUnknowns

    case (ADJOINT)
       do l = 1, solverOptions%nUnknowns
          do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
             do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
                do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                   gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *              &
                        (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                &
                        (k - 1 - this%gridOffset(3)))
                   if (grid%iblank(gridIndex) == 0) cycle
                   patchIndex = i - this%offset(1) + this%localSize(1) *                     &
                        (j - 1 - this%offset(2) + this%localSize(2) *                        &
                        (k - 1 - this%offset(3)))
                   state%rightHandSide(gridIndex, l) = state%rightHandSide(gridIndex, l) +   &
                        this%spongeStrength(patchIndex) * state%adjointVariables(gridIndex, l)
                end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
             end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
          end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
       end do !... l = 1, solverOptions%nUnknowns

    end select

    call endTiming("addDamping")

  end subroutine updateRhs

  subroutine computeStrength(this, grid, globalArcLengthsAlongDirection)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid

    implicit none

    ! <<< Arguments >>>
    class(t_SpongePatch) :: this
    class(t_Grid), intent(in) :: grid
    real(SCALAR_KIND), intent(in) :: globalArcLengthsAlongDirection(:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, direction, nDimensions
    real(wp), allocatable :: curveLengthIntegrand(:)

    if (this%nPatchPoints <= 0) return

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    direction = abs(this%normalDirection)
    assert(direction >= 1 .and. direction <= nDimensions)

#ifndef NDEBUG
    i = grid%nGridPoints / grid%localSize(direction) * grid%globalSize(direction)
    assert(size(globalArcLengthsAlongDirection) == i)
#endif

    allocate(curveLengthIntegrand(grid%globalSize(direction)))

    select case (direction)

    case (1)

       do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
          do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)

             do i = 1, grid%globalSize(1)
                curveLengthIntegrand(i) = real(globalArcLengthsAlongDirection(i +            &
                     grid%globalSize(1) * (j - 1 - grid%offset(2) + grid%localSize(2) *      &
                     (k - 1 - grid%offset(3)))), wp)
             end do

             if (this%normalDirection > 0) then
                do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                   this%spongeStrength(i - this%offset(1) + this%localSize(1) * (j - 1 -     &
                        this%offset(2) + this%localSize(2) * (k - 1 - this%offset(3)))) =    &
                        sum(curveLengthIntegrand(this%iMin : i - 1)) /                       &
                        sum(curveLengthIntegrand(this%iMin : this%iMax - 1))
                end do
             else
                do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                   this%spongeStrength(i - this%offset(1) + this%localSize(1) * (j - 1 -     &
                        this%offset(2) + this%localSize(2) * (k - 1 - this%offset(3)))) =    &
                        sum(curveLengthIntegrand(i + 1 : this%iMax)) /                       &
                        sum(curveLengthIntegrand(this%iMin + 1 : this%iMax))
                end do
             end if

          end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
       end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

    case (2)

       do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
          do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)

             do j = 1, grid%globalSize(2)
                curveLengthIntegrand(j) = real(globalArcLengthsAlongDirection(i -            &
                     grid%offset(1) + grid%localSize(1) * (j - 1 + grid%globalSize(2) *      &
                     (k - 1 - grid%offset(3)))), wp)
             end do

             if (this%normalDirection > 0) then
                do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
                   this%spongeStrength(i - this%offset(1) + this%localSize(1) * (j - 1 -     &
                        this%offset(2) + this%localSize(2) * (k - 1 - this%offset(3)))) =    &
                        sum(curveLengthIntegrand(this%jMin : j - 1)) /                       &
                        sum(curveLengthIntegrand(this%jMin : this%jMax - 1))
                end do
             else
                do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
                   this%spongeStrength(i - this%offset(1) + this%localSize(1) * (j - 1 -     &
                        this%offset(2) + this%localSize(2) * (k - 1 - this%offset(3)))) =    &
                        sum(curveLengthIntegrand(j + 1 : this%jMax)) /                       &
                        sum(curveLengthIntegrand(this%jMin + 1 : this%jMax))
                end do
             end if

          end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
       end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

    case (3)

       do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
          do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)

             do k = 1, grid%globalSize(3)
                curveLengthIntegrand(k) = real(globalArcLengthsAlongDirection(i -            &
                     grid%offset(1) + grid%localSize(1) * (j - 1 - grid%offset(2) +          &
                     grid%localSize(2) * (k - 1))), wp)
             end do

             if (this%normalDirection > 0) then
                do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
                   this%spongeStrength(i - this%offset(1) + this%localSize(1) * (j - 1 -     &
                        this%offset(2) + this%localSize(2) * (k - 1 - this%offset(3)))) =    &
                        sum(curveLengthIntegrand(this%kMin : k - 1)) /                       &
                        sum(curveLengthIntegrand(this%kMin : this%kMax - 1))
                end do
             else
                do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
                   this%spongeStrength(i - this%offset(1) + this%localSize(1) * (j - 1 -     &
                        this%offset(2) + this%localSize(2) * (k - 1 - this%offset(3)))) =    &
                        sum(curveLengthIntegrand(k + 1 : this%kMax)) /                       &
                        sum(curveLengthIntegrand(this%kMin + 1 : this%kMax))
                end do
             end if

          end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
       end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)

    end select !... select case (direction)

    this%spongeStrength = this%spongeAmount * (1.0_wp - this%spongeStrength) **              &
         this%spongeExponent

    SAFE_DEALLOCATE(curveLengthIntegrand)

  end subroutine computeStrength

end module SpongePatch_mod
