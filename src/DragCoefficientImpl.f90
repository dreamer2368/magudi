#include "config.h"

subroutine setupDragCoefficient(this, region)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use DragCoefficient_mod, only : t_DragCoefficient

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_DragCoefficient) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: nDimensions
  character(len = STRING_LENGTH) :: message

  nDimensions = region%grids(1)%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  call this%cleanup()

  call this%setupBase(region%simulationFlags, region%solverOptions)

  this%direction = 0.0_wp

  this%direction(1) = getOption("drag_direction_x", 1.0_wp)
  if (nDimensions >= 2) this%direction(2) = getOption("drag_direction_y", 0.0_wp)
  if (nDimensions == 3) this%direction(3) = getOption("drag_direction_z", 0.0_wp)

  if (sum(this%direction ** 2) <= epsilon(0.0_wp)) then
     write(message, '(A)')                                                                   &
          "Unable to determine a unit vector for computing the drag coefficient!"
     call gracefulExit(region%comm, message)
  end if

  this%direction = this%direction / sqrt(sum(this%direction ** 2))

end subroutine setupDragCoefficient

subroutine cleanupDragCoefficient(this)

  ! <<< Derived types >>>
  use DragCoefficient_mod, only : t_DragCoefficient

  implicit none

  ! <<< Arguments >>>
  class(t_DragCoefficient) :: this

  call this%cleanupBase()

end subroutine cleanupDragCoefficient

function computeDragCoefficient(this, time, region) result(instantaneousFunctional)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use DragCoefficient_mod, only : t_DragCoefficient

  ! <<< Arguments >>>
  class(t_DragCoefficient) :: this
  real(SCALAR_KIND), intent(in) :: time
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, ierror
  class(t_Patch), pointer :: patch => null()
  real(SCALAR_KIND) :: normBoundaryFactor
  SCALAR_TYPE, allocatable :: F(:,:)

  instantaneousFunctional = 0.0_wp

  do i = 1, size(region%grids)

     nDimensions = region%grids(i)%nDimensions
     assert_key(nDimensions, (1, 2, 3))

     if (.not. allocated(region%patchFactories)) cycle

     do j = 1, size(region%patchFactories)

        call region%patchFactories(j)%connect(patch)
        if (.not. associated(patch)) return
        if (patch%gridIndex /= region%grids(i)%index) return

        select type (patch)
        class is (t_CostTargetPatch)

           k = abs(patch%normalDirection)
           normBoundaryFactor = 1.0_wp / region%grids(i)%firstDerivative(k)%normBoundary(1)

           allocate(F(region%grids(i)%nGridPoints, 2))
           F(:,1) = (region%states(i)%pressure(:,1) -                                        &
                1.0_wp / region%solverOptions%ratioOfSpecificHeats)
           F(:,2) = matmul(region%grids(i)%metrics(:,1+nDimensions*(k-1):nDimensions*k),     &
                this%direction(1:nDimensions)) * normBoundaryFactor
           instantaneousFunctional = instantaneousFunctional +                               &
                patch%computeInnerProduct(region%grids(i), F(:,1), F(:,2))
           SAFE_DEALLOCATE(F)

        end select

     end do
  end do

  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE, instantaneousFunctional, 1,                          &
       SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(instantaneousFunctional, 1, SCALAR_TYPE_MPI,                             &
          0, region%grids(i)%comm, ierror)
  end do

  this%cachedValue = instantaneousFunctional

end function computeDragCoefficient

subroutine computeDragCoefficientAdjointForcing(this, simulationFlags, solverOptions,        &
     grid, state, patch)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use DragCoefficient_mod, only : t_DragCoefficient
  use SolverOptions_mod, only : t_SolverOptions
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use CNSHelper

  implicit none

  ! <<< Arguments >>>
  class(t_DragCoefficient) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  class(t_CostTargetPatch) :: patch

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, direction, nDimensions, nUnknowns, gridIndex, patchIndex
  real(SCALAR_KIND) :: normBoundaryFactor
  SCALAR_TYPE, allocatable :: localConservedVariables(:), metricsAlongNormalDirection(:),    &
       unitNormal(:), incomingJacobianOfInviscidFlux(:,:)
  SCALAR_TYPE :: F

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = abs(patch%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2)

  normBoundaryFactor = sign(1.0_wp / grid%firstDerivative(direction)%normBoundary(1),        &
       real(patch%normalDirection, wp))

  allocate(localConservedVariables(nUnknowns))
  allocate(unitNormal(nDimensions))
  allocate(metricsAlongNormalDirection(nDimensions))
  allocate(incomingJacobianOfInviscidFlux(nUnknowns, nUnknowns))

  do k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)
     do j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
        do i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
           gridIndex = i - patch%gridOffset(1) + patch%gridLocalSize(1) *                    &
                (j - 1 - patch%gridOffset(2) + patch%gridLocalSize(2) *                      &
                (k - 1 - patch%gridOffset(3)))
           if (grid%iblank(gridIndex) == 0) cycle
           patchIndex = i - patch%offset(1) + patch%localSize(1) *                           &
                (j - 1 - patch%offset(2) + patch%localSize(2) *                              &
                (k - 1 - patch%offset(3)))

           localConservedVariables = state%conservedVariables(gridIndex,:)
           metricsAlongNormalDirection =                                                     &
                grid%metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)
           unitNormal = metricsAlongNormalDirection /                                        &
                sqrt(sum(metricsAlongNormalDirection ** 2))

           if (simulationFlags%useContinuousAdjoint) then

              F = grid%jacobian(gridIndex, 1) *                                              &
                   dot_product(state%adjointVariables(gridIndex,2:nDimensions+1) -           &
                   this%direction(1:nDimensions), unitNormal)

              select case (nDimensions)
              case (1)
                 call computeIncomingJacobianOfInviscidFlux1D(localConservedVariables,       &
                      metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,       &
                      - patch%normalDirection, incomingJacobianOfInviscidFlux,               &
                      specificVolume = state%specificVolume(gridIndex, 1),                   &
                      temperature = state%temperature(gridIndex, 1))
              case (2)
                 call computeIncomingJacobianOfInviscidFlux2D(localConservedVariables,       &
                      metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,       &
                      - patch%normalDirection, incomingJacobianOfInviscidFlux,               &
                      specificVolume = state%specificVolume(gridIndex, 1),                   &
                      temperature = state%temperature(gridIndex, 1))
              case (3)
                 call computeIncomingJacobianOfInviscidFlux3D(localConservedVariables,       &
                      metricsAlongNormalDirection, solverOptions%ratioOfSpecificHeats,       &
                      - patch%normalDirection, incomingJacobianOfInviscidFlux,               &
                      specificVolume = state%specificVolume(gridIndex, 1),                   &
                      temperature = state%temperature(gridIndex, 1))
              end select !... nDimensions

              patch%adjointForcing(patchIndex,:) = - patch%inviscidPenaltyAmount * F *       &
                   matmul(transpose(incomingJacobianOfInviscidFlux(2:nDimensions+1,:)),      &
                   unitNormal)

           else

              F = grid%jacobian(gridIndex, 1) * normBoundaryFactor *                         &
                   (solverOptions%ratioOfSpecificHeats - 1.0_wp) *                           &
                   dot_product(metricsAlongNormalDirection, this%direction(1:nDimensions))

              patch%adjointForcing(patchIndex,1) =                                           &
                   0.5_wp * sum(state%velocity(gridIndex,:) ** 2) * F
              patch%adjointForcing(patchIndex,2:nDimensions+1) =                             &
                   - state%velocity(gridIndex,:) * F
              patch%adjointForcing(patchIndex,nDimensions+2) = F

           end if

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
     end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

  SAFE_DEALLOCATE(unitNormal)
  SAFE_DEALLOCATE(metricsAlongNormalDirection)

end subroutine computeDragCoefficientAdjointForcing

function isDragCoefficientPatchValid(this, patchDescriptor, gridSize, normalDirection,       &
     extent, simulationFlags, message) result(isPatchValid)

  ! <<< Derived types >>>
  use DragCoefficient_mod, only : t_DragCoefficient
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_DragCoefficient) :: this
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  integer, intent(in) :: gridSize(:), normalDirection, extent(6)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Result >>>
  logical :: isPatchValid

  ! <<< Local variables >>>
  integer :: i

  isPatchValid = .false.
  if (normalDirection > size(gridSize) .or. normalDirection == 0) then
     write(message, '(A)') "Normal direction is invalid!"
     return
  end if

  i = abs(normalDirection)
  if (extent((i-1)*2+1) /= extent((i-1)*2+2)) then
     write(message, '(A)') "Extends more than 1 grid point along normal direction!"
     return
  end if

  if ((normalDirection > 0 .and. extent((i-1)*2+1) /= 1) .or.                                &
       (normalDirection < 0 .and. extent((i-1)*2+2) /= gridSize(i))) then
     write(message, '(2(A,I0.0),A)') "Not aligned with a computational boundary!"
     return
  end if

  isPatchValid = .true.

end function isDragCoefficientPatchValid
