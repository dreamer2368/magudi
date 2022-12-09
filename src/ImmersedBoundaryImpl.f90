#include "config.h"

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

end function verifyImmersedBoundaryPatchUsage

subroutine addImmersedBoundaryPenalty(this)
  use ImmersedBoundaryPatch_mod, only : t_ImmersedBoundaryPatch
  class(t_ImmersedBoundaryPatch) :: this
end subroutine addImmersedBoundaryPenalty

subroutine addImmersedBoundaryLevelset(this)
  use ImmersedBoundaryPatch_mod, only : t_ImmersedBoundaryPatch
  class(t_ImmersedBoundaryPatch) :: this
end subroutine addImmersedBoundaryLevelset

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
  integer :: i, j, n, nUnknowns
  real(wp) :: buf
  real(wp), dimension(:, :), allocatable :: densityGradient, dissipationTerm
  real(wp), dimension(:), allocatable :: objectVelocity

  call startTiming("updateIBMVariables")

  if (mode .ne. FORWARD) then
    return
  end if

  assert(size(this%levelset, 1) == size(this%conservedVariables, 1))
  assert(size(this%levelset, 2) == 1)

  nUnknowns = size(this%conservedVariables, 2)

  allocate(densityGradient(grid%nGridPoints, grid%nDimensions))
  allocate(dissipationTerm(grid%nGridPoints, nUnknowns))
  allocate(objectVelocity(grid%nDimensions))

  !TODO: compute levelset.
  this%levelset(:, 1) = grid%coordinates(:, 2) - 0.0_wp

  ! compute levelset normal.
  call grid%computeGradient(this%levelset(:, 1), this%levelsetNormal)

  do i = 1, grid%nGridPoints
     ! Make the levelset normal a unit norm
     buf = sqrt(sum(this%levelsetNormal(i,:) ** 2))
     if (buf .gt. 0.0_wp) this%levelsetNormal(i,:) = this%levelsetNormal(i,:) / buf

     ! ! Compute indicator function
     ! if (levelset(i, 1) .le. 0.0_wp) then
     !    indicatorFunction(i, 1) = 0.0_wp
     ! else
     !    indicatorFunction(i, 1) = 1.0_wp
     ! end if
     !
     ! ! Update the grid norm
     ! gridNorm(i, 1) = primitiveGridNorm(i, 1) * indicatorFunction(i, 1)
  end do

  !TODO: object velocity.
  objectVelocity = 0.0_wp

  this%nDotGradRho = 0.0_wp
  this%uDotGradRho = 0.0_wp
  this%ibmDissipation = 0.0_wp
  call grid%computeGradient(this%conservedVariables(:,1), densityGradient)
  do n = 1, grid%nDimensions
    this%nDotGradRho(:, 1) = this%nDotGradRho(:, 1) + this%levelsetNormal(:, n) * densityGradient(:, n)
    this%uDotGradRho(:, 1) = this%uDotGradRho(:, 1) + objectVelocity(n) * densityGradient(:, n)

    ! Dissipation term.
    dissipationTerm = this%conservedVariables
    call grid%dissipation(n)%apply(dissipationTerm, grid%localSize)
    if (.not. simulationFlags%compositeDissipation) then
       do j = 1, nUnknowns
          dissipationTerm(:,j) = - grid%arcLengths(:,n) * dissipationTerm(:,j)
       end do
       call grid%dissipationTranspose(n)%apply(dissipationTerm, grid%localSize)
       call grid%firstDerivative(n)%applyNormInverse(dissipationTerm, grid%localSize)
    end if
    this%ibmDissipation = this%ibmDissipation + dissipationTerm
  end do

  call endTiming("updateIBMVariables")

end subroutine updateIBMVariables
