#include "config.h"

module CostTargetPatch_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  use Patch_mod, only : t_Patch

  implicit none
  private

  type, extends(t_Patch), public :: t_CostTargetPatch

     real(SCALAR_KIND), allocatable :: adjointForcing(:,:)

   contains

     procedure, pass :: setup
     procedure, pass :: cleanup
     procedure, pass :: updateRhs

  end type t_CostTargetPatch

contains

  subroutine setup(this, name, comm, grid, state, extent,                                    &
       normalDirection, simulationFlags, solverOptions)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod, only : t_State
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Internal modules >>>
    use ErrorHandler, only : issueWarning

    implicit none

    ! <<< Arguments >>>
    class(t_CostTargetPatch) :: this
    character(len = *), intent(in) :: name
    integer, intent(in) :: comm
    class(t_Grid), intent(in) :: grid
    class(t_State), intent(in) :: state
    integer, intent(in) :: extent(6), normalDirection
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    character(len = STRING_LENGTH) :: message
    integer :: nDimensions, nUnknowns, direction

    call this%cleanup()
    call this%setupBase(name, comm, grid, extent, normalDirection)

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns >= nDimensions + 2)

    direction = abs(this%normalDirection)

    if (simulationFlags%predictionOnly) then
       write(message, '(3A)') "Patch '", trim(name), "' will not be used!"
       call issueWarning(grid%comm, message)
    end if

    if (.not. simulationFlags%predictionOnly .and. this%nPatchPoints >= 0)                   &
         allocate(this%adjointForcing(this%nPatchPoints, nUnknowns))

  end subroutine setup

  subroutine cleanup(this)

    implicit none

    ! <<< Arguments >>>
    class(t_CostTargetPatch) :: this

    call this%cleanupBase()

    SAFE_DEALLOCATE(this%adjointForcing)

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
    class(t_CostTargetPatch) :: this
    integer, intent(in) :: mode
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid), intent(in) :: grid
    class(t_State) :: state

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, l, nDimensions, nUnknowns, gridIndex, patchIndex
    real(SCALAR_KIND) :: forcingFactor

    assert_key(mode, (FORWARD, ADJOINT))

    if (mode == FORWARD) return

    call startTiming("addCostTargetPenalty")

    if (simulationFlags%useContinuousAdjoint .or. simulationFlags%steadyStateSimulation) then
       forcingFactor = 1.0_wp
    else
       forcingFactor = state%adjointForcingFactor
    end if

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns >= nDimensions + 2)

    do l = 1, nUnknowns
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

                state%rightHandSide(gridIndex,l) = state%rightHandSide(gridIndex,l) +        &
                     forcingFactor * this%adjointForcing(patchIndex,l)

             end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
          end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
       end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
    end do !... l = 1, nUnknowns

    call endTiming("addCostTargetPenalty")

  end subroutine updateRhs

end module CostTargetPatch_mod
