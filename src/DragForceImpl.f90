#include "config.h"

subroutine setupDragForce(this, region)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use DragForce_mod, only : t_DragForce

  ! <<< Internal modules >>>
  use InputHelper, only : getOption
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_DragForce) :: this
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
          "Unable to determine a unit vector for computing drag force!"
     call gracefulExit(region%comm, message)
  end if

  this%direction = this%direction / sqrt(sum(this%direction ** 2))

end subroutine setupDragForce

subroutine cleanupDragForce(this)

  ! <<< Derived types >>>
  use DragForce_mod, only : t_DragForce

  implicit none

  ! <<< Arguments >>>
  class(t_DragForce) :: this

  call this%cleanupBase()

end subroutine cleanupDragForce

function computeDragForce(this, region) result(instantaneousFunctional)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use DragForce_mod, only : t_DragForce

  ! <<< Arguments >>>
  class(t_DragForce) :: this
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, nPatches, nDimensions, ierror
  class(t_Patch), pointer :: patch => null()
  real(SCALAR_KIND) :: normBoundaryFactor
  SCALAR_TYPE, allocatable :: F(:,:)

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(size(region%grids) == size(region%states))

  instantaneousFunctional = 0.0_wp

  do i = 1, size(region%grids)

     nDimensions = region%grids(i)%nDimensions
     assert_key(nDimensions, (1, 2, 3))

     assert(region%grids(i)%nGridPoints > 0)
     assert(allocated(region%grids(i)%firstDerivative))
     assert(size(region%grids(i)%firstDerivative) == nDimensions)
     assert(allocated(region%grids(i)%targetMollifier))
     assert(size(region%grids(i)%targetMollifier, 1) == region%grids(i)%nGridPoints)
     assert(size(region%grids(i)%targetMollifier, 2) == 1)
     assert(allocated(region%grids(i)%metrics))
     assert(size(region%grids(i)%metrics, 1) == region%grids(i)%nGridPoints)
     assert(size(region%grids(i)%metrics, 2) == nDimensions ** 2)
     assert(allocated(region%states(i)%pressure))
     assert(size(region%states(i)%pressure, 1) == region%grids(i)%nGridPoints)
     assert(size(region%states(i)%pressure, 2) == 1)

#ifdef DEBUG
     if (region%simulationFlags%viscosityOn) then
        assert(allocated(region%states(i)%stressTensor))
        assert(size(region%states(i)%stressTensor, 1) == region%grids(i)%nGridPoints)
        assert(size(region%states(i)%stressTensor, 2) == nDimensions ** 2)
     end if
#endif

     nPatches = 0
     if (allocated(region%patchFactories)) nPatches = size(region%patchFactories)

     do j = 1, nPatches

        call region%patchFactories(j)%connect(patch)
        if (.not. associated(patch)) cycle
        if (patch%gridIndex /= region%grids(i)%index) cycle

        select type (patch)
        class is (t_CostTargetPatch)

           k = abs(patch%normalDirection)
           assert(allocated(region%grids(i)%firstDerivative(k)%normBoundary))
           assert(size(region%grids(i)%firstDerivative(k)%normBoundary) > 0)
           normBoundaryFactor = 1.0_wp / region%grids(i)%firstDerivative(k)%normBoundary(1)

           allocate(F(region%grids(i)%nGridPoints, 1))

           F(:,1) = 0.0_wp
           do l = 1, nDimensions
              if (region%simulationFlags%viscosityOn) then
                 F(:,1) = F(:,1) + this%direction(l) *                                       &
                      sum(region%grids(i)%metrics(:,1+nDimensions*(k-1):nDimensions*k) *     &
                      region%states(i)%stressTensor(:,1+nDimensions*(l-1):nDimensions*l),    &
                      dim = 2)
              end if
           end do
           F(:,1) = normBoundaryFactor * F(:,1)

           instantaneousFunctional = instantaneousFunctional +                               &
                patch%computeInnerProduct(region%grids(i), F(:,1),                           &
                region%grids(i)%targetMollifier(:,1))

           SAFE_DEALLOCATE(F)

        end select

     end do

     call MPI_Allreduce(MPI_IN_PLACE, instantaneousFunctional, 1,                            &
          SCALAR_TYPE_MPI, MPI_SUM, region%grids(i)%comm, ierror)

  end do

  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE, instantaneousFunctional, 1,                          &
       SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(instantaneousFunctional, 1, SCALAR_TYPE_MPI,                             &
          0, region%grids(i)%comm, ierror)
  end do

  this%cachedValue = instantaneousFunctional

end function computeDragForce

subroutine computeDragForceSpatialDistribution(this, grid, state, F)

  use Grid_mod, only : t_Grid
  use State_mod, only : t_State

  use DragForce_mod, only : t_DragForce

  class(t_DragForce) :: this
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  SCALAR_TYPE, intent(out) :: F(:,:)

end subroutine computeDragForceSpatialDistribution

subroutine computeDragForceAdjointForcing(this, simulationFlags, solverOptions,              &
     grid, state, patch)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use DragForce_mod, only : t_DragForce
  use SolverOptions_mod, only : t_SolverOptions
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use CNSHelper

  implicit none

  ! <<< Arguments >>>
  class(t_DragForce) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  class(t_CostTargetPatch) :: patch

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: direction, nDimensions, nUnknowns
  real(SCALAR_KIND) :: normBoundaryFactor
  SCALAR_TYPE, allocatable :: temp1(:,:), temp2(:,:)

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  direction = abs(patch%normalDirection)
  assert(direction >= 1 .and. direction <= nDimensions)

  nUnknowns = solverOptions%nUnknowns
  assert(nUnknowns >= nDimensions + 2)

  normBoundaryFactor = 1.0_wp / grid%firstDerivative(direction)%normBoundary(1)

  if (patch%nPatchPoints > 0) then

     allocate(temp1(grid%nGridPoints, nUnknowns))
     allocate(temp2(grid%nGridPoints, 1))

     ! Hack for TBL:

     temp1 = 0.0_wp

     temp2(:,1) = grid%jacobian(:,1) * grid%metrics(:,1) * state%dynamicViscosity(:,1)
     call grid%adjointFirstDerivative(1)%projectOnBoundaryAndApply(temp2, grid%localSize,    &
          patch%normalDirection)
     call grid%firstDerivative(1)%applyNorm(temp2, grid%localSize)
     temp1(:,3) = grid%jacobian(:,1) * state%specificVolume(:,1) * temp2(:,1)

     temp2(:,1) = grid%jacobian(:,1) * grid%metrics(:,5) * state%dynamicViscosity(:,1)
     call grid%adjointFirstDerivative(2)%projectOnBoundaryAndApply(temp2, grid%localSize,    &
          patch%normalDirection)
     call grid%firstDerivative(1)%applyNorm(temp2, grid%localSize)
     temp1(:,2) = grid%jacobian(:,1) * state%specificVolume(:,1) * temp2(:,1)

     call patch%collect(temp1, patch%adjointForcing)
     patch%adjointForcing = patch%adjointForcing * normBoundaryFactor

     SAFE_DEALLOCATE(temp2)
     SAFE_DEALLOCATE(temp1)

  end if

end subroutine computeDragForceAdjointForcing

function isDragForcePatchValid(this, patchDescriptor, gridSize, normalDirection,             &
     extent, simulationFlags, message) result(isPatchValid)

  ! <<< Derived types >>>
  use DragForce_mod, only : t_DragForce
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_DragForce) :: this
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

  isPatchValid = .true.

end function isDragForcePatchValid
