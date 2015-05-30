#include "config.h"

subroutine setupFlameTemperature(this, region)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use FlameTemperature_mod, only : t_FlameTemperature

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables
  use InputHelper, only : getOption, getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_FlameTemperature) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, ierror

  assert(allocated(region%states))
  assert(size(region%states) > 0)

  call this%cleanup()

  call this%setupBase(region%simulationFlags, region%solverOptions)

  allocate(this%data_(size(region%globalGridSizes, 2)))

  do i = 1, size(this%data_)
     do j = 1, size(region%grids)
        if (region%grids(j)%index == i) then
           allocate(this%data_(i)%meanTemperature(1, 1))
           region%states(j)%dummyFunction => this%data_(i)%meanTemperature
           exit
        end if
     end do
     call MPI_Barrier(region%comm, ierror)
  end do

  do i = 1, size(this%data_)
     this%data_(i)%meanTemperature = 1.0_wp / ( (region%solverOptions%ratioOfSpecificHeats  -&
          1.0_wp) * (1.0_wp - region%states(1)%combustion%heatRelease) )
  end do

end subroutine setupFlameTemperature

subroutine cleanupFlameTemperature(this)

  ! <<< Derived types >>>
  use FlameTemperature_mod, only : t_FlameTemperature

  implicit none

  ! <<< Arguments >>>
  class(t_FlameTemperature) :: this

  ! <<< Local variables >>>
  integer :: i

  call this%cleanupBase()

  if (allocated(this%data_)) then
     do i = 1, size(this%data_)
        if (associated(this%data_(i)%meanTemperature)) deallocate(this%data_(i)%meanTemperature)
        nullify(this%data_(i)%meanTemperature)
     end do
  end if
  SAFE_DEALLOCATE(this%data_)

end subroutine cleanupFlameTemperature

function computeFlameTemperature(this, region) result(instantaneousFunctional)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use FlameTemperature_mod, only : t_FlameTemperature

  ! <<< Arguments >>>
  class(t_FlameTemperature) :: this
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, ierror
  SCALAR_TYPE, allocatable :: F(:,:)

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(size(region%grids) == size(region%states))

  instantaneousFunctional = 0.0_wp

  do i = 1, size(region%grids)

     assert(region%grids(i)%nGridPoints > 0)
     assert(allocated(region%grids(i)%targetMollifier))
     assert(size(region%grids(i)%targetMollifier, 1) == region%grids(i)%nGridPoints)
     assert(size(region%grids(i)%targetMollifier, 2) == 1)
     assert(allocated(region%states(i)%pressure))
     assert(size(region%states(i)%pressure, 1) == region%grids(i)%nGridPoints)
     assert(size(region%states(i)%pressure, 2) == 1)

     j = region%grids(i)%index

     assert(associated(this%data_(j)%meanTemperature))
     assert(size(this%data_(j)%meanTemperature, 1) == 1)
     assert(size(this%data_(j)%meanTemperature, 2) == 1)

     allocate(F(region%grids(i)%nGridPoints, 1))
     F = region%states(i)%pressure - this%data_(j)%meanTemperature
     instantaneousFunctional = instantaneousFunctional +                                     &
          region%grids(i)%computeInnerProduct(F, F, region%grids(i)%targetMollifier(:,1))
     SAFE_DEALLOCATE(F)

  end do

  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE, instantaneousFunctional, 1,                          &
       SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(instantaneousFunctional, 1, SCALAR_TYPE_MPI,                             &
          0, region%grids(i)%comm, ierror)
  end do

  this%cachedValue = instantaneousFunctional

end function computeFlameTemperature

subroutine computeFlameTemperatureAdjointForcing(this, simulationFlags, solverOptions,          &
     grid, state, patch)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use FlameTemperature_mod, only : t_FlameTemperature
  use SolverOptions_mod, only : t_SolverOptions
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_FlameTemperature) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  class(t_CostTargetPatch) :: patch

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: meanTemperature(:)
  SCALAR_TYPE :: F

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  allocate(meanTemperature(patch%nPatchPoints))
  i = grid%index
  call patch%collect(this%data_(i)%meanTemperature(:,1), meanTemperature)

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

           F = - 2.0_wp * grid%targetMollifier(gridIndex, 1) *                               &
                (solverOptions%ratioOfSpecificHeats - 1.0_wp) *                              &
                (state%pressure(gridIndex, 1) - meanTemperature(patchIndex))

           patch%adjointForcing(patchIndex,nDimensions+2) = F
           patch%adjointForcing(patchIndex,2:nDimensions+1) =                                &
                - state%velocity(gridIndex,:) * F
           patch%adjointForcing(patchIndex,1) =                                              &
                0.5_wp * sum(state%velocity(gridIndex,:) ** 2) * F

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
     end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

  SAFE_DEALLOCATE(meanTemperature)

end subroutine computeFlameTemperatureAdjointForcing

function isFlameTemperaturePatchValid(this, patchDescriptor, gridSize, normalDirection,         &
     extent, simulationFlags, message) result(isPatchValid)

  ! <<< Derived types >>>
  use FlameTemperature_mod, only : t_FlameTemperature
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_FlameTemperature) :: this
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  integer, intent(in) :: gridSize(:), normalDirection, extent(6)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Result >>>
  logical :: isPatchValid

  ! <<< Local variables >>>
  integer :: i, n

  isPatchValid = .false.

  n = size(gridSize)

  do i = 1, size(gridSize)
     if (extent((i-1)*2+1) == extent((i-1)*2+2)) n = n - 1
  end do

  if (n /= size(gridSize)) then
     write(message, '(2(A,I0.0),A)') "Expected a ", size(gridSize),                          &
          "D patch, but extent represents a ", n, "D patch!"
     return
  end if

  isPatchValid = .true.

end function isFlameTemperaturePatchValid
