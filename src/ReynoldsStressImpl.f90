#include "config.h"

subroutine setupReynoldsStress(this, region)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use ReynoldsStress_mod, only : t_ReynoldsStress

  ! <<< Enumerations >>>
  use State_enum, only : QOI_DUMMY_FUNCTION

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_ReynoldsStress) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions, ierror
  character(len = STRING_LENGTH) :: message, filename

  assert(allocated(region%states))
  assert(size(region%states) > 0)

  nDimensions = region%grids(1)%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  call this%cleanup()

  call this%setupBase(region%simulationFlags, region%solverOptions)

  allocate(this%data_(size(region%globalGridSizes, 2)))

  do i = 1, size(this%data_)
     do j = 1, size(region%grids)
        if (region%grids(j)%index == i) then
           allocate(this%data_(i)%meanVelocity(region%grids(j)%nGridPoints, nDimensions))
           region%states(j)%dummyFunction => this%data_(i)%meanVelocity
           exit
        end if
     end do
     call MPI_Barrier(region%comm, ierror)
  end do

  if (.not. region%simulationFlags%useTargetState) then
     call getRequiredOption("mean_velocity_file", filename, region%comm)
     call region%loadData(QOI_DUMMY_FUNCTION, filename)
  else
     filename = getOption("mean_velocity_file", "")
     if (len_trim(filename) > 0) then
        call region%loadData(QOI_DUMMY_FUNCTION, filename)
     else
        do i = 1, size(region%states)
           assert(allocated(region%states(i)%targetState))
           j = region%grids(i)%index
           call computeDependentVariables(region%grids(i)%nDimensions,                       &
                region%states(i)%targetState, region%solverOptions%ratioOfSpecificHeats,     &
                velocity = this%data_(j)%meanVelocity)
        end do
     end if
  end if

  this%firstDirection  = 0.0_wp
  this%secondDirection = 0.0_wp

  this%firstDirection(1) = getOption("reynolds_stress_direction1_x", 1.0_wp)
  if (nDimensions >= 2)                                                                      &
       this%firstDirection(2) = getOption("reynolds_stress_direction1_y", 0.0_wp)
  if (nDimensions == 3)                                                                      &
       this%firstDirection(3) = getOption("reynolds_stress_direction1_z", 0.0_wp)

  this%secondDirection(1) = getOption("reynolds_stress_direction2_x", 1.0_wp)
  if (nDimensions >= 2)                                                                      &
       this%secondDirection(2) = getOption("reynolds_stress_direction2_y", 0.0_wp)
  if (nDimensions == 3)                                                                      &
       this%secondDirection(3) = getOption("reynolds_stress_direction2_z", 0.0_wp)

  if (sum(this%firstDirection ** 2) <= epsilon(0.0_wp) .or.                                  &
       sum(this%secondDirection ** 2) <= epsilon(0.0_wp)) then
     write(message, '(A)')                                                                   &
          "Unable to determine unit vectors for computing Reynolds stress!"
     call gracefulExit(region%comm, message)
  end if

  this%firstDirection = this%firstDirection / sqrt(sum(this%firstDirection ** 2))
  this%secondDirection = this%secondDirection / sqrt(sum(this%secondDirection ** 2))

end subroutine setupReynoldsStress

subroutine cleanupReynoldsStress(this)

  ! <<< Derived types >>>
  use ReynoldsStress_mod, only : t_ReynoldsStress

  implicit none

  ! <<< Arguments >>>
  class(t_ReynoldsStress) :: this

  ! <<< Local variables >>>
  integer :: i

  call this%cleanupBase()

  if (allocated(this%data_)) then
     do i = 1, size(this%data_)
        if (associated(this%data_(i)%meanVelocity)) deallocate(this%data_(i)%meanVelocity)
        nullify(this%data_(i)%meanVelocity)
     end do
  end if
  SAFE_DEALLOCATE(this%data_)

end subroutine cleanupReynoldsStress

function computeReynoldsStress(this, region) result(instantaneousFunctional)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use ReynoldsStress_mod, only : t_ReynoldsStress

  ! <<< Arguments >>>
  class(t_ReynoldsStress) :: this
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, ierror
  SCALAR_TYPE, allocatable :: F(:,:)

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(size(region%grids) == size(region%states))

  instantaneousFunctional = 0.0_wp

  do i = 1, size(region%grids)

     nDimensions = region%grids(i)%nDimensions
     assert_key(nDimensions, (1, 2, 3))

     assert(region%grids(i)%nGridPoints > 0)
     assert(allocated(region%grids(i)%targetMollifier))
     assert(size(region%grids(i)%targetMollifier, 1) == region%grids(i)%nGridPoints)
     assert(size(region%grids(i)%targetMollifier, 2) == 1)
     assert(allocated(region%states(i)%pressure))
     assert(size(region%states(i)%pressure, 1) == region%grids(i)%nGridPoints)
     assert(size(region%states(i)%pressure, 2) == 1)

     j = region%grids(i)%index

     assert(associated(this%data_(j)%meanVelocity))
     assert(size(this%data_(j)%meanVelocity, 1) == region%grids(i)%nGridPoints)
     assert(size(this%data_(j)%meanVelocity, 2) == nDimensions)

     allocate(F(region%grids(i)%nGridPoints, 1))
     do k = 1, region%grids(i)%nGridPoints
        F(k,1) = 0.5_wp * dot_product(region%states(i)%velocity(k,:) -                       &
             this%data_(j)%meanVelocity(k,:), this%firstDirection(1:nDimensions)) *          &
             dot_product(region%states(i)%velocity(k,:) -                                    &
             this%data_(j)%meanVelocity(k,:), this%secondDirection(1:nDimensions))
     end do
     instantaneousFunctional = instantaneousFunctional +                                     &
          region%grids(i)%computeInnerProduct(F(:,1), region%grids(i)%targetMollifier(:,1))
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

end function computeReynoldsStress

subroutine computeReynoldsStressAdjointForcing(this, simulationFlags, solverOptions,          &
     grid, state, patch)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use ReynoldsStress_mod, only : t_ReynoldsStress
  use SolverOptions_mod, only : t_SolverOptions
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_ReynoldsStress) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  class(t_CostTargetPatch) :: patch

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: meanVelocity(:,:)
  SCALAR_TYPE :: F

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  allocate(meanVelocity(patch%nPatchPoints, nDimensions))
  i = grid%index
  call patch%collect(this%data_(i)%meanVelocity, meanVelocity)

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

           F = - 0.5_wp * grid%targetMollifier(gridIndex, 1) *                               &
                state%specificVolume(gridIndex,1) *                                          &
                dot_product(this%firstDirection(1:nDimensions),                              &
                state%velocity(gridIndex,:) - meanVelocity(patchIndex,:))

           patch%adjointForcing(patchIndex,2:nDimensions+1) =                                &
                this%secondDirection(1:nDimensions) * F
           patch%adjointForcing(patchIndex,1) =                                              &
                - dot_product(this%secondDirection(1:nDimensions),                           &
                state%velocity(gridIndex,:)) * F

           F = - 0.5_wp * grid%targetMollifier(gridIndex, 1) *                               &
                state%specificVolume(gridIndex,1) *                                          &
                dot_product(this%secondDirection(1:nDimensions),                             &
                state%velocity(gridIndex,:) - meanVelocity(patchIndex,:))

           patch%adjointForcing(patchIndex,2:nDimensions+1) =                                &
                this%firstDirection(1:nDimensions) * F
           patch%adjointForcing(patchIndex,1) =                                              &
                - dot_product(this%firstDirection(1:nDimensions),                            &
                state%velocity(gridIndex,:)) * F

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
     end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

  SAFE_DEALLOCATE(meanVelocity)

end subroutine computeReynoldsStressAdjointForcing

function isReynoldsStressPatchValid(this, patchDescriptor, gridSize, normalDirection,        &
     extent, simulationFlags, message) result(isPatchValid)

  ! <<< Derived types >>>
  use ReynoldsStress_mod, only : t_ReynoldsStress
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_ReynoldsStress) :: this
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

end function isReynoldsStressPatchValid
