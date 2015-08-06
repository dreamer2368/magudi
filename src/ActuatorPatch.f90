#include "config.h"

module ActuatorPatch_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  use MPI
  use ProbePatch_mod, only : t_ProbePatch

  implicit none
  private

  type, extends(t_ProbePatch), public :: t_ActuatorPatch

     real(SCALAR_KIND), allocatable :: gradient(:,:), controlForcing(:,:)

   contains

     procedure, pass :: setup
     procedure, pass :: cleanup
     procedure, pass :: updateRhs
     procedure, pass :: update

  end type t_ActuatorPatch

contains

  subroutine setup(this, name, comm, grid, state, extent,                                    &
       normalDirection, simulationFlags, solverOptions)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod, only : t_State
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : FORWARD, ADJOINT

    ! <<< Internal modules >>>
    use InputHelper, only : getOption
    use ErrorHandler, only : issueWarning

    implicit none

    ! <<< Arguments >>>
    class(t_ActuatorPatch) :: this
    character(len = *), intent(in) :: name
    integer, intent(in) :: comm
    class(t_Grid), intent(in) :: grid
    class(t_State), intent(in) :: state
    integer, intent(in) :: extent(6), normalDirection
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    character(len = STRING_LENGTH) :: message, outputPrefix
    integer :: nDimensions, nUnknowns

    call this%cleanup()
    call this%t_ProbePatch%setup(name, comm, grid, state, extent,                            &
         normalDirection, simulationFlags, solverOptions)

    this%saveMode = ADJOINT
    this%loadMode = FORWARD

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns >= nDimensions + 2)

    if (simulationFlags%predictionOnly) then
       write(message, '(A)') "Patch '", trim(name), "' will not be used!"
       call issueWarning(grid%comm, message)
    end if

    outputPrefix = getOption("output_prefix", PROJECT_NAME)
    write(this%probeFilename, '(4A)') trim(outputPrefix), ".gradient_", trim(name), ".dat"

    if (.not. simulationFlags%predictionOnly .and. this%nPatchPoints > 0)                    &
         allocate(this%controlForcing(this%nPatchPoints, nUnknowns), source = 0.0_wp)

  end subroutine setup

  subroutine cleanup(this)

    implicit none

    ! <<< Arguments >>>
    class(t_ActuatorPatch) :: this

    call this%t_ProbePatch%cleanup()

    SAFE_DEALLOCATE(this%gradient)
    SAFE_DEALLOCATE(this%controlForcing)

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
    class(t_ActuatorPatch) :: this
    integer, intent(in) :: mode
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid), intent(in) :: grid
    class(t_State) :: state

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, l, nDimensions, nUnknowns, gridIndex, patchIndex

    assert_key(mode, (FORWARD, ADJOINT))

    if (mode == ADJOINT .or. (mode == FORWARD .and.                                          &
         abs(state%actuationAmount) <= 0.0_wp)) return

    call startTiming("updateActuatorPatch")

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

                state%rightHandSide(gridIndex, l) = state%rightHandSide(gridIndex, l) +      &
                     grid%controlMollifier(gridIndex, 1) * this%controlForcing(patchIndex, l)

             end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
          end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
       end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
    end do !... l = 1, nUnknowns

    call endTiming("updateActuatorPatch")

  end subroutine updateRhs

  subroutine update(this, mode, simulationFlags, solverOptions, grid, state)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod, only : t_State
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : FORWARD, ADJOINT

    implicit none

    ! <<< Arguments >>>
    class(t_ActuatorPatch) :: this
    integer, intent(in) :: mode
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid), intent(in) :: grid
    class(t_State) :: state

    ! <<< Local variables >>>
    integer :: nDimensions, nUnknowns

    if (this%comm == MPI_COMM_NULL) return

    assert_key(mode, (FORWARD, ADJOINT))
    if (mode == FORWARD) return

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns >= nDimensions + 2)

    assert(allocated(this%probeBuffer))
    assert(this%iProbeBuffer >= 1 .and. this%iProbeBuffer <= size(this%probeBuffer, 3))
    assert(allocated(this%gradient))
    assert(size(this%probeBuffer, 2) == size(this%gradient, 2))

    this%probeBuffer(:,:,this%iProbeBuffer) = this%gradient

  end subroutine update

end module ActuatorPatch_mod
