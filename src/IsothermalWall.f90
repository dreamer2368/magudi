#include "config.h"

module IsothermalWall_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  use ImpenetrableWall_mod, only : t_ImpenetrableWall

  implicit none
  private

  type, extends(t_ImpenetrableWall), public :: t_IsothermalWall

     real(SCALAR_KIND) :: viscousPenaltyAmount
     real(SCALAR_KIND), allocatable :: temperature(:), dynamicViscosity(:),                  &
          secondCoefficientOfViscosity(:), thermalDiffusivity(:)

   contains

     procedure, pass :: setup
     procedure, pass :: cleanup
     procedure, pass :: updateRhs

  end type t_IsothermalWall

contains

  subroutine setup(this, name, comm, grid, state, extent,                                    &
       normalDirection, simulationFlags, solverOptions)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod, only : t_State
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Internal modules >>>
    use CNSHelper, only : computeDependentVariables, computeTransportVariables
    use InputHelper, only : getOption
    use ErrorHandler, only : gracefulExit

    implicit none

    ! <<< Arguments >>>
    class(t_IsothermalWall) :: this
    character(len = *), intent(in) :: name
    integer, intent(in) :: comm
    class(t_Grid), intent(in) :: grid
    class(t_State), intent(in) :: state
    integer, intent(in) :: extent(6), normalDirection
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    character(len = STRING_LENGTH) :: key
    integer :: nDimensions, nUnknowns, direction
    real(wp), allocatable :: targetState(:,:)
    real(wp) :: wallTemperature

    call this%cleanup()
    call this%t_ImpenetrableWall%setup(name, comm, grid, state, extent,                      &
         normalDirection, simulationFlags, solverOptions)

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns >= nDimensions + 2)

    direction = abs(this%normalDirection)

    write(key, '(A)') "patches/" // trim(name) // "/"

    ! Viscous penalty amount.
    this%viscousPenaltyAmount = getOption("defaults/viscous_penalty_amount", 1.0_wp)
    this%viscousPenaltyAmount = getOption(trim(key) // "viscous_penalty_amount",             &
         this%viscousPenaltyAmount)
    this%viscousPenaltyAmount = this%viscousPenaltyAmount *                                  &
         solverOptions%reynoldsNumberInverse
    this%viscousPenaltyAmount = this%viscousPenaltyAmount /                                  &
         grid%firstDerivative(direction)%normBoundary(1)

    if (simulationFlags%viscosityOn .and. this%nPatchPoints > 0) then

       allocate(this%temperature(this%nPatchPoints))
       allocate(this%dynamicViscosity(this%nPatchPoints))
       allocate(this%secondCoefficientOfViscosity(this%nPatchPoints))
       allocate(this%thermalDiffusivity(this%nPatchPoints))

       if (simulationFlags%useTargetState .and. allocated(state%targetState)) then

          allocate(targetState(this%nPatchPoints, nUnknowns))
          call this%collect(state%targetState, targetState)
          call computeDependentVariables(nDimensions, targetState,                           &
               solverOptions%ratioOfSpecificHeats, temperature = this%temperature)
          SAFE_DEALLOCATE(targetState)

       else

          wallTemperature = 1.0_wp / (solverOptions%ratioOfSpecificHeats - 1.0_wp)
          wallTemperature = getOption(trim(key) // "temperature", wallTemperature)
          this%temperature = wallTemperature

       end if

       call computeTransportVariables(this%temperature, solverOptions%powerLawExponent,      &
            solverOptions%bulkViscosityRatio, solverOptions%ratioOfSpecificHeats,            &
            solverOptions%reynoldsNumberInverse, solverOptions%prandtlNumberInverse,         &
            this%dynamicViscosity, this%secondCoefficientOfViscosity, this%thermalDiffusivity)

    end if

  end subroutine setup

  subroutine cleanup(this)

    implicit none

    ! <<< Arguments >>>
    class(t_IsothermalWall) :: this

    call this%t_ImpenetrableWall%cleanup()

    SAFE_DEALLOCATE(this%temperature)
    SAFE_DEALLOCATE(this%dynamicViscosity)
    SAFE_DEALLOCATE(this%secondCoefficientOfViscosity)
    SAFE_DEALLOCATE(this%thermalDiffusivity)

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
    class(t_IsothermalWall) :: this
    integer, intent(in) :: mode
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid), intent(in) :: grid
    class(t_State) :: state

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, nDimensions, nUnknowns, direction, gridIndex, patchIndex
    real(wp), allocatable :: localPenalty(:)

    assert_key(mode, (FORWARD, ADJOINT))

    call this%t_ImpenetrableWall%updateRhs(mode, simulationFlags, solverOptions, grid, state)

    call startTiming("addIsothermalWallPenalty")

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    direction = abs(this%normalDirection)
    assert(direction >= 1 .and. direction <= nDimensions)

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns >= nDimensions + 2)

    allocate(localPenalty(nUnknowns - 1))

    select case (mode)

    case (FORWARD)

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

                localPenalty(1:nDimensions) =                                                &
                     state%conservedVariables(gridIndex,2:nDimensions+1)
                localPenalty(nDimensions+1) =                                                &
                     state%conservedVariables(gridIndex,nDimensions+2) -                     &
                     state%conservedVariables(gridIndex,1) * this%temperature(patchIndex) /  &
                     solverOptions%ratioOfSpecificHeats

                state%rightHandSide(gridIndex,2:nUnknowns) =                                 &
                     state%rightHandSide(gridIndex,2:nUnknowns) -                            &
                     this%viscousPenaltyAmount * grid%jacobian(gridIndex, 1) * localPenalty

             end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
          end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
       end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

    case (ADJOINT)

    end select

    SAFE_DEALLOCATE(localPenalty)

    call endTiming("addIsothermalWallPenalty")

  end subroutine updateRhs

end module IsothermalWall_mod
