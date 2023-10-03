#include "config.h"

module ImmersedBoundaryImpl

  implicit none
  public

contains

  ! Regularized a Heaviside function
  ! phi: signed-distance (levelset)
  ! eps: regularization parameters (length)
  ! ---------------------------------------
  function regularizeHeaviside(phi, eps) result(H)

    implicit none

    ! Arguments
    integer, parameter :: wp = SCALAR_KIND
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    real(wp), intent(in) :: phi, eps
    real(wp)             :: H

    if (phi .le. -eps) then
       H = 0.0_wp
    else if (abs(phi) .lt. eps) then
       H = 0.5_wp * (1.0_wp + phi / eps + 1.0_wp / pi * sin(phi * pi / eps))
    else
       H = 1.0_wp
    end if

    return
  end function regularizeHeaviside

end module ImmersedBoundaryImpl

subroutine setupImmersedBoundaryPatch(this, index, comm, patchDescriptor,                            &
                                      grid, simulationFlags, solverOptions)
  ! <<< External modules >>>
  use MPI
  use iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use ImmersedBoundaryPatch_mod, only : t_ImmersedBoundaryPatch
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush, gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_ImmersedBoundaryPatch) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: message

  call this%cleanup()
  call this%setupBase(index, comm, patchDescriptor, grid, simulationFlags, solverOptions)

  if (simulationFlags%enableAdjoint) then
    write(message, '(A)') "ImmersedBoundaryPatch does not support adjoint mode!"
    call gracefulExit(this%comm, message)
  end if

  ! Set the max speed based on the reference sound speed
  this%maxSpeed = 1.0_wp

  ! Store inverse timestep size
  this%dti = 1.0_wp / solverOptions%timeStepSize

  ! Set the diffusion amount based on the stability limit
  assert(grid%minGridSpacing > 0.0_wp)
  this%dissipationAmount = getOption("immersed_boundary/dissipation_amount", 0.05_wp)
  this%dissipationAmount = this%dissipationAmount * grid%minGridSpacing**2 * this%dti / real(grid%nDimensions, wp)

  this%ibmEpsilon = getOption("immersed_boundary/regularization_parameter", 1.0_wp)

  ! Import specific heat ratio
  this%ratioOfSpecificHeats = solverOptions%ratioOfSpecificHeats

  ! Set wall temperature
  this%ibmTemperature = getOption("immersed_boundary/wall_temperature",         &
                                  1.0_wp / (solverOptions%ratioOfSpecificHeats - 1.0_wp))

  ! ! Store the unmodified grid norm
  ! allocate(this%primitiveGridNorm(this%nPatchPoints, 1))
  ! call this%collect(grid%norm, this%primitiveGridNorm)
end subroutine setupImmersedBoundaryPatch

subroutine cleanupImmersedBoundaryPatch(this)
  use ImmersedBoundaryPatch_mod, only : t_ImmersedBoundaryPatch

  implicit none

  ! <<< Arguments >>>
  class(t_ImmersedBoundaryPatch) :: this

  call this%cleanupBase()
end subroutine cleanupImmersedBoundaryPatch

function verifyImmersedBoundaryPatchUsage(this, patchDescriptor, gridSize, normalDirection,        &
                                   extent, simulationFlags, success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use ImmersedBoundaryPatch_mod, only : t_ImmersedBoundaryPatch
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_ImmersedBoundaryPatch) :: this
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  integer, intent(in) :: gridSize(:), normalDirection, extent(6)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  logical, intent(out) :: success
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Result >>>
  logical :: isPatchUsed

  ! <<< Local variables >>>
  integer :: i

  isPatchUsed = .false.

  success = .false.

  do i = 1, size(gridSize)
     if (extent((i-1)*2+1) < 0 .or. extent((i-1)*2+2) > gridSize(i) .or.                     &
          extent((i-1)*2+1) > extent((i-1)*2+2)) then
        write(message, '(A)') "Invalid extent!"
        return
     end if
  end do

  success = .true.
  isPatchUsed = .true.

end function verifyImmersedBoundaryPatchUsage

subroutine addImmersedBoundaryPenalty(this, mode, simulationFlags, solverOptions, grid, state)
  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use ImmersedBoundaryPatch_mod, only : t_ImmersedBoundaryPatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT, LINEARIZED
  use IBM_enum

  ! <<< Internal module >>>
  use ImmersedBoundaryImpl, only : regularizeHeaviside

  implicit none

  ! <<< Arguments >>>
  class(t_ImmersedBoundaryPatch) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: nUnknowns, nDimensions
  integer :: i, j, k, n, gridIndex, patchIndex
  real(wp) :: localDensity, localVelocity(grid%nDimensions),                    &
              localDensityPenalty, localTemperaturePenalty, localEnergy,        &
              localIbmDissipation(size(state%conservedVariables, 2)),           &
              localLevelsetNormal(grid%nDimensions)
  real(wp) :: objectVelocity(grid%nDimensions), source_(size(state%conservedVariables, 2)), velocityPenalty(grid%nDimensions)
  real(wp) :: buf, weight

  !TODO: implement adjoint.
  assert_key(mode, (FORWARD))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  nUnknowns = size(state%conservedVariables, 2)
  nDimensions = grid%nDimensions

  objectVelocity = 0.0_wp

  ! Update the source terms and IBM forcing
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

        localDensity = state%conservedVariables(gridIndex, 1)
        localVelocity = state%velocity(gridIndex, :)
        localEnergy = state%conservedVariables(gridIndex, grid%nDimensions + 2) /         &
                      state%conservedVariables(gridIndex, 1)
        localDensityPenalty = state%densityPenalty(gridIndex, 1)

        ! localuDotGradRho = state%uDotGradRho(gridIndex, 1)
        localIbmDissipation = state%ibmDissipation(gridIndex, :)
        localLevelsetNormal = state%levelsetNormal(gridIndex, :)

        ! Get velocity of associated object
        objectVelocity = state%objectVelocity(gridIndex, :)

        ! Compute the velocity penalty
        velocityPenalty = objectVelocity - localVelocity

        ! Zero-out the local source terms
        source_ = 0.0_WP

        ! Density treatment
        source_(1) = this%maxSpeed * localDensityPenalty +            &
                    ! - localuDotGradRho            &
                    this%dissipationAmount * localIbmDissipation(1)

        ! Momentum treatment
        do n = 1, nDimensions
          source_(n+1) = localDensity * velocityPenalty(n) * this%dti +               &
               localVelocity(n) * this%maxSpeed * localDensityPenalty +               &
               this%dissipationAmount * localIbmDissipation(n+1)
        end do

        ! Energy treatment
        source_(nDimensions+2) = localEnergy * this%maxSpeed * localDensityPenalty +                                &
            sum(state%conservedVariables(gridIndex, 2:nDimensions+1) * velocityPenalty(1:nDimensions)) * this%dti + &
            this%dissipationAmount * localIbmDissipation(nDimensions+2)

        select case (solverOptions%ibmWallType)

        case (IBM_ADIABATIC)

          localTemperaturePenalty = state%temperaturePenalty(gridIndex, 1)
          source_(nDimensions+2) = source_(nDimensions+2) +                                    &
               localDensity / this%ratioOfSpecificHeats * this%maxSpeed * localTemperaturePenalty

        case (IBM_ISOTHERMAL)

          localTemperaturePenalty = this%ibmTemperature - state%temperature(gridIndex,1)
          source_(nDimensions+2) = source_(nDimensions+2) +                                    &
               localDensity / this%ratioOfSpecificHeats * this%dti * localTemperaturePenalty

        end select

        ! Add the IBM contribution
        buf = this%ibmEpsilon * sqrt(sum((localLevelsetNormal * grid%gridSpacing(gridIndex, :))**2))
        weight = 1.0_wp - regularizeHeaviside(state%levelset(gridIndex, 1), buf)
        if (state%levelset(gridIndex, 1) .le. 0.0_wp) then
          state%rightHandSide(gridIndex,:) = source_ * weight
        else
          state%rightHandSide(gridIndex,:) = state%rightHandSide(gridIndex,:)   &
                                              + source_ * weight
        end if
      end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
    end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
end subroutine addImmersedBoundaryPenalty

subroutine computePenaltyWeight(this, grid, state, weight)
  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use ImmersedBoundaryPatch_mod, only : t_ImmersedBoundaryPatch

  ! <<< Internal module >>>
  use ImmersedBoundaryImpl, only : regularizeHeaviside

  implicit none

  ! <<< Arguments >>>
  class(t_ImmersedBoundaryPatch) :: this
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state
  real(SCALAR_KIND), dimension(:), intent(out) :: weight

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, gridIndex, patchIndex
  real(wp) :: localLevelsetNormal(grid%nDimensions)
  real(wp) :: buf

  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))
  assert(size(weight) == grid%nGridPoints)

  ! Update the source terms and IBM forcing
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

        localLevelsetNormal = state%levelsetNormal(gridIndex, :)

        ! Add the IBM contribution
        buf = this%ibmEpsilon * sqrt(sum((localLevelsetNormal * grid%gridSpacing(gridIndex, :))**2))
        weight(gridIndex) = 1.0_wp - regularizeHeaviside(state%levelset(gridIndex, 1), buf)

      end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
    end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
end subroutine computePenaltyWeight

! t_State method
subroutine updateIBMVariables(this, mode, grid, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SimulationFlags_mod, only : t_SimulationFlags
  use SolverOptions_mod, only : t_SolverOptions

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD
  use IBM_enum

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  integer, intent(in) :: mode
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  integer :: i, nDimensions
  real(wp), dimension(:, :), allocatable :: densityGradient,      &
                                            temperatureGradient

  call startTiming("updateIBMVariables")

  if (mode .ne. FORWARD) then
    return
  end if

  !if (.not. this%ibmPatchExists) return

  assert(size(this%levelset, 1) == size(this%conservedVariables, 1))
  assert(size(this%levelset, 2) == 1)

  nDimensions = grid%nDimensions

  allocate(densityGradient(grid%nGridPoints, nDimensions))
  allocate(temperatureGradient(grid%nGridPoints, nDimensions))

  !NOTE: remnant from jcode.
  ! do i = 1, grid%nGridPoints
     ! ! Compute indicator function
     ! if (levelset(i, 1) .le. 0.0_wp) then
     !    indicatorFunction(i, 1) = 0.0_wp
     ! else
     !    indicatorFunction(i, 1) = 1.0_wp
     ! end if
     !
     ! ! Update the grid norm
     ! gridNorm(i, 1) = primitiveGridNorm(i, 1) * indicatorFunction(i, 1)
  ! end do

  this%densityPenalty = 0.0_wp
  this%ibmDissipation = 0.0_wp
  call grid%computeGradient(this%conservedVariables(:,1), densityGradient)
  call grid%computeGradient(this%temperature(:,1), temperatureGradient)

  do i = 1, grid%nGridPoints
    this%densityPenalty(i, 1) = this%densityPenalty(i, 1) + sum(this%levelsetNormal(i, :) * densityGradient(i, :))
    this%densityPenalty(i, 1) = this%densityPenalty(i, 1) +                                   &
      solverOptions%ratioOfSpecificHeats / (solverOptions%ratioOfSpecificHeats - 1.0_wp) *    &
      this%conservedVariables(i, 1) / this%temperature(i, 1) *                                &
      sum(this%levelsetNormal(i, :) * this%objectAcceleration(i, :))
    ! this%uDotGradRho(i, 1) = this%uDotGradRho(i, 1) + sum(this%objectVelocity(i, :) * densityGradient(i, :))
  end do

  select case (solverOptions%ibmWallType)

  case (IBM_ADIABATIC)
    assert(allocated(this%temperaturePenalty))
    this%temperaturePenalty = 0.0_wp
    do i = 1, grid%nGridPoints
      this%temperaturePenalty(i, 1) = this%temperaturePenalty(i, 1) + sum(this%levelsetNormal(i, :) * temperatureGradient(i, :))
    end do

  case (IBM_ISOTHERMAL)
    do i = 1, grid%nGridPoints
      this%densityPenalty(i, 1) = this%densityPenalty(i, 1) +                                   &
        this%conservedVariables(i, 1) / this%temperature(i, 1) *                                &
        sum(this%levelsetNormal(i, :) * temperatureGradient(i, :))
    end do

  end select

  ! Laplacian term.
  call grid%computeLaplacian(this%conservedVariables, this%ibmDissipation)

  SAFE_DEALLOCATE(densityGradient)
  SAFE_DEALLOCATE(temperatureGradient)

  call endTiming("updateIBMVariables")

end subroutine updateIBMVariables
