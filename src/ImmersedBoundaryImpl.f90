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
  this%dissipationAmount = getOption("immersed_boundary/dissipation_amount", 0.05_wp)
  ! this%dissipationAmount = this%dissipationAmount * grid%minGridSpacing * this%dti / real(grid%nDimensions, wp)
  ! this%dissipationAmount = this%dissipationAmount * grid%minGridSpacing**2 * this%dti / real(grid%nDimensions, wp)

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
  real(wp) :: localDensity, localVelocity(grid%nDimensions), localVelocitySquared,   &
              localnDotGradRho, localuDotGradRho, localIbmDissipation(size(state%conservedVariables, 2)), &
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
        localVelocitySquared = sum(localVelocity ** 2)
        localnDotGradRho = state%nDotGradRho(gridIndex, 1)
        localuDotGradRho = state%uDotGradRho(gridIndex, 1)
        localIbmDissipation = state%ibmDissipation(gridIndex, :)
        localLevelsetNormal = state%levelsetNormal(gridIndex, :)

        ! Get velocity of associated object
        objectVelocity = state%objectVelocity(gridIndex, :)
     ! if (ibm_move) then
     !    n = objectIndex(i)
     !    objectVelocity = object(n)%velocity(1:nDimensions)
     !    ibmVelocity(i,:) = objectVelocity
     ! else
     !    objectVelocity = 0.0_WP
     ! end if

        ! Compute the velocity penalty
        velocityPenalty = objectVelocity - localVelocity

        ! Zero-out the local source terms
        source_ = 0.0_WP

        ! Density treatment
        source_(1) = this%maxSpeed * localnDotGradRho - localuDotGradRho            &
                    + this%dissipationAmount * localIbmDissipation(1)

        ! Momentum treatment
        do n = 1, nDimensions
          source_(n+1) = localDensity * velocityPenalty(n) * this%dti +                                     &
               localVelocity(n) * (this%maxSpeed * localnDotGradRho - localuDotGradRho) +                       &
               this%dissipationAmount * localIbmDissipation(n+1)
        end do

        ! Energy treatment
        source_(nDimensions+2) = 0.5_wp * localVelocitySquared * (this%maxSpeed * localnDotGradRho -           &
            localuDotGradRho) + this%dissipationAmount * localIbmDissipation(nDimensions+2) +   &
            sum(state%conservedVariables(gridIndex, 2:nDimensions+1) * velocityPenalty(1:nDimensions)) * this%dti

        ! Isothermal
        source_(nDimensions+2) = source_(nDimensions+2) +                                    &
             localDensity * (this%ibmTemperature - state%temperature(gridIndex,1)) / this%ratioOfSpecificHeats * this%dti

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

! t_State method
subroutine updateIBMVariables(this, mode, grid, simulationFlags)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  integer, intent(in) :: mode
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  integer :: i, j, n, nUnknowns
  real(wp), dimension(:, :), allocatable :: densityGradient, dissipationTerm

  call startTiming("updateIBMVariables")

  if (mode .ne. FORWARD) then
    return
  end if

  if (.not. this%ibmPatchExists) return

  assert(size(this%levelset, 1) == size(this%conservedVariables, 1))
  assert(size(this%levelset, 2) == 1)

  nUnknowns = size(this%conservedVariables, 2)

  allocate(densityGradient(grid%nGridPoints, grid%nDimensions))
  allocate(dissipationTerm(grid%nGridPoints, nUnknowns))

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

  this%nDotGradRho = 0.0_wp
  this%uDotGradRho = 0.0_wp
  this%ibmDissipation = 0.0_wp
  call grid%computeGradient(this%conservedVariables(:,1), densityGradient)

  do i = 1, grid%nGridPoints
    this%nDotGradRho(i, 1) = this%nDotGradRho(i, 1) + sum(this%levelsetNormal(i, :) * densityGradient(i, :))
    this%uDotGradRho(i, 1) = this%uDotGradRho(i, 1) + sum(this%objectVelocity(i, :) * densityGradient(i, :))
  end do

  do n = 1, grid%nDimensions
    ! Dissipation term.
    dissipationTerm = this%conservedVariables
    call grid%dissipation(n)%apply(dissipationTerm, grid%localSize)
    if (.not. simulationFlags%compositeDissipation) then
       do j = 1, nUnknowns
          dissipationTerm(:,j) = grid%arcLengths(:,n) * dissipationTerm(:,j)
       end do
       call grid%dissipationTranspose(n)%apply(dissipationTerm, grid%localSize)
       call grid%firstDerivative(n)%applyNormInverse(dissipationTerm, grid%localSize)
    end if
    this%ibmDissipation = this%ibmDissipation + dissipationTerm
  end do

  SAFE_DEALLOCATE(densityGradient)
  SAFE_DEALLOCATE(dissipationTerm)

  call endTiming("updateIBMVariables")

end subroutine updateIBMVariables
