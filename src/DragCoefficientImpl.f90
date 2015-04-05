#include "config.h"

module DragCoefficientImpl

  implicit none
  public

contains

  subroutine computeAdjointForcingOnPatch(patch, grid, velocity, pressure,                   &
       targetPressure, ratioOfSpecificHeats, forceDirection)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use Patch_mod, only : t_Patch
    use CostTargetPatch_mod, only : t_CostTargetPatch

    ! <<< Arguments >>>
    class(t_Patch), pointer, intent(in) :: patch
    class(t_Grid), intent(in) :: grid
    SCALAR_TYPE, intent(in) :: velocity(:,:), pressure(:), targetPressure
    real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats, forceDirection(3)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, direction, nDimensions, gridIndex, patchIndex
    real(SCALAR_KIND) :: normBoundaryFactor
    SCALAR_TYPE, allocatable :: localMetricsAlongDirection(:)

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    allocate(localMetricsAlongDirection(nDimensions))

    select type (patch)
    class is (t_CostTargetPatch)

       direction = abs(patch%normalDirection)
       normBoundaryFactor = 1.0_wp / grid%firstDerivative(direction)%normBoundary(1)

       do k = patch%offset(3) + 1, patch%offset(3) + patch%patchSize(3)
          do j = patch%offset(2) + 1, patch%offset(2) + patch%patchSize(2)
             do i = patch%offset(1) + 1, patch%offset(1) + patch%patchSize(1)
                gridIndex = i - patch%gridOffset(1) + patch%gridLocalSize(1) *               &
                     (j - 1 - patch%gridOffset(2) + patch%gridLocalSize(2) *                 &
                     (k - 1 - patch%gridOffset(3)))
                if (grid%iblank(gridIndex) == 0) cycle
                patchIndex = i - patch%offset(1) + patch%patchSize(1) *                      &
                     (j - 1 - patch%offset(2) + patch%patchSize(2) *                         &
                     (k - 1 - patch%offset(3)))

                localMetricsAlongDirection =                                                 &
                     grid%metrics(gridIndex,1+nDimensions*(direction-1):nDimensions*direction)

                patch%adjointForcing(patchIndex,nDimensions+2) =                             &
                     - 2.0_wp * (ratioOfSpecificHeats - 1.0_wp) *                            &
                     (pressure(gridIndex) - targetPressure) *                                &
                     normBoundaryFactor * dot_product(forceDirection(1:nDimensions),         &
                     localMetricsAlongDirection) ** 2 /                                      &
                     sqrt(sum(localMetricsAlongDirection ** 2))
                patch%adjointForcing(patchIndex,2:nDimensions+1) = - velocity(gridIndex,:) * &
                     patch%adjointForcing(patchIndex,nDimensions+2)
                patch%adjointForcing(patchIndex,1) =                                         &
                     0.5_wp * sum(velocity(gridIndex,:) ** 2) *                              &
                     patch%adjointForcing(patchIndex,nDimensions+2)

             end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%patchSize(1)
          end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%patchSize(2)
       end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%patchSize(3)

    end select

    SAFE_DEALLOCATE(localMetricsAlongDirection)

  end subroutine computeAdjointForcingOnPatch

end module DragCoefficientImpl

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

  this%freeStreamPressure = getOption("free_stream_pressure",                                &
       1.0_wp / region%solverOptions%ratioOfSpecificHeats)
  if (this%freeStreamPressure <= 0.0_wp) then
     write(message, '(A,E11.2,A)') "Free stream pressure is invalid: ",                      &
          this%freeStreamPressure, " <= 0!"
     call gracefulExit(region%comm, message)
  end if

  this%direction = 0.0_wp

  this%direction(1) = getOption("free_stream_velocity_x", 0.0_wp)
  if (nDimensions >= 2) this%direction(2) = getOption("free_stream_velocity_y", 0.0_wp)
  if (nDimensions == 3) this%direction(3) = getOption("free_stream_velocity_z", 0.0_wp)

  if (sum(this%direction ** 2) <= epsilon(0.0_wp)) then
     write(message, '(A)') "Unable to determine a unit vector along the free stream flow &
          &direction required for computing the drag coefficient!"
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
  SCALAR_TYPE, allocatable :: F(:,:)

  instantaneousFunctional = 0.0_wp

  if (.not. allocated(region%patchFactories)) return

  do i = 1, size(region%grids)

     nDimensions = region%grids(i)%nDimensions
     assert_key(nDimensions, (1, 2, 3))

     do j = 1, size(region%patchFactories)

        call region%patchFactories(j)%connect(patch)
        if (.not. associated(patch)) return
        if (patch%gridIndex /= region%grids(i)%index) return

        select type (patch)
        class is (t_CostTargetPatch)

           k = abs(patch%normalDirection)

           allocate(F(region%grids(i)%nGridPoints, 2))
           F(:,1) = (region%states(i)%pressure(:,1) - this%freeStreamPressure) *             &
                matmul(region%grids(i)%metrics(:,1+nDimensions*(k-1):nDimensions*k),         &
                this%direction(1:nDimensions))
           F(:,2) = 1.0_wp /                                                                 &
                sqrt(sum(region%grids(i)%metrics(:,1+nDimensions*(k-1):nDimensions*k) ** 2))
           instantaneousFunctional = instantaneousFunctional +                               &
                patch%computeInnerProduct(region%grids(i), F(:,1), F(:,1), F(:,2))
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

subroutine computeDragCoefficientAdjointForcing(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use DragCoefficient_mod, only : t_DragCoefficient
  use CostTargetPatch_mod, only : t_CostTargetPatch

  ! <<< Enumerations >>>
  use State_enum, only : QOI_DUMMY_FUNCTION

  ! <<< Private members >>>
  use DragCoefficientImpl, only : computeAdjointForcingOnPatch

  implicit none

  ! <<< Arguments >>>
  class(t_DragCoefficient) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer :: i, j, nDimensions
  class(t_Patch), pointer :: patch => null()

  if (.not. allocated(region%patchFactories)) return

  do i = 1, size(region%grids)

     nDimensions = region%grids(i)%nDimensions

     do j = 1, size(region%patchFactories)
        call region%patchFactories(j)%connect(patch)
        if (.not. associated(patch)) cycle
        if (patch%gridIndex /= region%grids(i)%index .or. patch%nPatchPoints <= 0) cycle

        call computeAdjointForcingOnPatch(patch, region%grids(i), region%states(i)%velocity, &
             region%states(i)%pressure(:,1), this%freeStreamPressure,                        &
             region%solverOptions%ratioOfSpecificHeats, this%direction)

     end do

  end do

end subroutine computeDragCoefficientAdjointForcing

function isDragCoefficientPatchValid(this, patchDescriptor, gridSize, normalDirection,       &
     extent, simulationFlags, message) result(isPatchValid)

  ! <<< Derived types >>>
  use AcousticNoise_mod, only : t_AcousticNoise
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_AcousticNoise) :: this
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
