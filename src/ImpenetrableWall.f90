#include "config.h"

module ImpenetrableWall_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  use Patch_mod, only : t_Patch

  implicit none
  private

  type, extends(t_Patch), public :: t_ImpenetrableWall

     real(SCALAR_KIND) :: inviscidPenaltyAmount

   contains

     procedure, pass :: setup
     procedure, pass :: cleanup
     procedure, pass :: updateRhs

  end type t_ImpenetrableWall

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
    class(t_ImpenetrableWall) :: this
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
    integer :: nDimensions, direction

    call this%cleanup()
    call this%setupBase(name, comm, grid, extent, normalDirection)

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    direction = abs(this%normalDirection)

    if (direction <= 0 .or. direction > nDimensions) then
       write(message, '(3A)') "Invalid normal direction for patch '", trim(name), "'!"
       call gracefulExit(grid%comm, message)
    end if

    if (extent((direction-1)*2+1) /= extent((direction-1)*2+2)) then
       write(message, '(A)') "Patch '", trim(name),                                          &
            "' is not allowed to extend more than 1 grid point along normal direction!"
       call gracefulExit(grid%comm, message)
    end if

    write(key, '(A)') "patches/" // trim(name) // "/"

    ! Inviscid penalty amount.
    this%inviscidPenaltyAmount = getOption("defaults/inviscid_penalty_amount", 1.0_wp)
    this%inviscidPenaltyAmount = getOption(trim(key) // "inviscid_penalty_amount",           &
         this%inviscidPenaltyAmount)
    this%inviscidPenaltyAmount = sign(this%inviscidPenaltyAmount,                            &
         real(this%normalDirection, wp))
    this%inviscidPenaltyAmount = this%inviscidPenaltyAmount /                                &
         grid%firstDerivative(direction)%normBoundary(1)

  end subroutine setup

  subroutine cleanup(this)

    implicit none

    ! <<< Arguments >>>
    class(t_ImpenetrableWall) :: this

    call this%cleanupBase()

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
    class(t_ImpenetrableWall) :: this
    integer, intent(in) :: mode
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    class(t_Grid), intent(in) :: grid
    class(t_State) :: state

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, l, nDimensions, nUnknowns, direction, gridIndex, patchIndex
    real(wp), allocatable :: localPenalty(:), deltaPressure(:), deltaInviscidPenalty(:,:),   &
         localMetricsAlongNormalDirection(:)
    real(wp) :: localNormalMomentum

    assert_key(mode, (FORWARD, ADJOINT))
    assert(allocated(state%targetState))

    call startTiming("addImpenetrableWallPenalty")

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    direction = abs(this%normalDirection)
    assert(direction >= 1 .and. direction <= nDimensions)

    nUnknowns = solverOptions%nUnknowns
    assert(nUnknowns >= nDimensions + 2)

    allocate(localPenalty(nUnknowns))
    allocate(localMetricsAlongNormalDirection(nDimensions))

    if (mode == ADJOINT) then
       allocate(deltaPressure(nUnknowns))
       if (.not. simulationFlags%useContinuousAdjoint)                                       &
            allocate(deltaInviscidPenalty(nUnknowns, nUnknowns))
    end if

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

                localMetricsAlongNormalDirection =                                           &
                     grid%metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

                localNormalMomentum = dot_product(state%conservedVariables(gridIndex,        &
                     2:nDimensions+1), localMetricsAlongNormalDirection)

                localPenalty(1) = localNormalMomentum
                localPenalty(2:nDimensions+1) = localNormalMomentum *                        &
                     state%velocity(gridIndex,:)
                localPenalty(nDimensions+2) = localNormalMomentum *                          &
                     state%specificVolume(gridIndex,1) *                                     &
                     (state%conservedVariables(gridIndex,nDimensions+2) +                    &
                     state%pressure(gridIndex,1))

                state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) -        &
                     this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) * localPenalty

             end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
          end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
       end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

    case (ADJOINT)

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

                localMetricsAlongNormalDirection =                                           &
                     grid%metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

                deltaPressure(1) = 0.5_wp * sum(state%velocity(gridIndex,:) ** 2)
                deltaPressure(2:nDimensions+1) = - state%velocity(gridIndex,:)
                deltaPressure(nDimensions+2) = 1.0_wp
                deltaPressure = (solverOptions%ratioOfSpecificHeats - 1.0_wp) *              &
                     deltaPressure

                if (simulationFlags%useContinuousAdjoint) then

                   localNormalMomentum = dot_product(state%adjointVariables(gridIndex,       &
                        2:nDimensions+1), localMetricsAlongNormalDirection)

                   state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +     &
                        this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *           &
                        localNormalMomentum * deltaPressure

                else

                   call computeInviscidJacobian(nDimensions, gridIndex,                      &
                        state%conservedVariables,                                            &
                        grid%metrics(:,1+nDimensions*(direction-1):nDimensions*direction),   &
                        solverOptions%ratioOfSpecificHeats, deltaInviscidPenalty,            &
                        specificVolume = state%specificVolume(:,1),                          &
                        velocity = state%velocity, temperature = state%temperature(:,1))

                   do l = 1, nDimensions
                      deltaInviscidPenalty(l+1,:) = deltaInviscidPenalty(l+1,:) -            &
                           localMetricsAlongNormalDirection(l) * deltaPressure
                   end do

                   state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:) +     &
                        this%inviscidPenaltyAmount * grid%jacobian(gridIndex, 1) *           &
                        matmul(transpose(deltaInviscidPenalty),                              &
                        state%adjointVariables(gridIndex,:))

                end if

             end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
          end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
       end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

    end select

    SAFE_DEALLOCATE(deltaInviscidPenalty)
    SAFE_DEALLOCATE(deltaPressure)
    SAFE_DEALLOCATE(localMetricsAlongNormalDirection)
    SAFE_DEALLOCATE(localPenalty)

    call endTiming("addImpenetrableWallPenalty")

  end subroutine updateRhs

end module ImpenetrableWall_mod
