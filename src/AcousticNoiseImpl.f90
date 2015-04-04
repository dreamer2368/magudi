#include "config.h"

module AcousticNoiseImpl

  implicit none
  public

contains

  subroutine computeAdjointForcingOnPatch(patch, grid, velocity,                             &
       pressure, meanPressure, ratioOfSpecificHeats)

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use Patch_mod, only : t_Patch
    use CostTargetPatch_mod, only : t_CostTargetPatch

    ! <<< Arguments >>>
    class(t_Patch), pointer, intent(in) :: patch
    class(t_Grid), intent(in) :: grid
    SCALAR_TYPE, intent(in) :: velocity(:,:), pressure(:), meanPressure(:)
    real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, nDimensions, gridIndex, patchIndex

    nDimensions = grid%nDimensions
    assert_key(nDimensions, (1, 2, 3))

    select type (patch)
    type is (t_CostTargetPatch)

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

                patch%adjointForcing(patchIndex,nDimensions+2) =                             &   
                     - 2.0_wp * (ratioOfSpecificHeats - 1.0_wp) *                            &   
                     (pressure(gridIndex) - meanPressure(gridIndex))
                patch%adjointForcing(patchIndex,2:nDimensions+1) = - velocity(gridIndex,:) * &   
                     patch%adjointForcing(patchIndex,nDimensions+2)
                patch%adjointForcing(patchIndex,1) =                                         &
                     0.5_wp * sum(velocity(gridIndex,:) ** 2) *                              &
                     patch%adjointForcing(patchIndex,nDimensions+2)

             end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%patchSize(1)
          end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%patchSize(2)
       end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%patchSize(3)

    end select
    
  end subroutine computeAdjointForcingOnPatch

end module AcousticNoiseImpl

subroutine setupAcousticNoise(this, region)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use AcousticNoise_mod, only : t_AcousticNoise

  ! <<< Enumerations >>>
  use State_enum, only : QOI_DUMMY_FUNCTION

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables
  use InputHelper, only : getOption, getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_AcousticNoise) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer :: i
  character(len = STRING_LENGTH) :: filename

  assert(allocated(region%states))
  assert(size(region%states) > 0)

  call this%cleanup()

  call this%setupBase(region%simulationFlags, region%solverOptions)

  allocate(this%data_(size(region%states)))
  do i = 1, size(this%data_)
     allocate(this%data_(i)%meanPressure(region%grids(i)%nGridPoints, 1))
     region%states(i)%dummyFunction => this%data_(i)%meanPressure
  end do
  
  if (.not. region%simulationFlags%useTargetState) then
     call getRequiredOption("mean_pressure_file", filename)
     call loadRegionData(QOI_DUMMY_FUNCTION, filename)
  else
     filename = getOption("mean_pressure_file", "")
     if (len_trim(filename) > 0) then
        call loadRegionData(QOI_DUMMY_FUNCTION, filename)
     else
        do i = 1, size(this%data_)
           assert(allocated(region%states(i)%targetState))
           call computeDependentVariables(region%grids(i)%nDimensions,                       &
                region%states(i)%targetState, region%solverOptions%ratioOfSpecificHeats,     &
                pressure = this%data_(i)%meanPressure(:,1))
        end do
     end if
  end if

end subroutine setupAcousticNoise

subroutine cleanupAcousticNoise(this)

  ! <<< Derived types >>>
  use AcousticNoise_mod, only : t_AcousticNoise

  implicit none

  ! <<< Arguments >>>
  class(t_AcousticNoise) :: this

  ! <<< Local variables >>>
  integer :: i

  call this%cleanupBase()

  if (allocated(this%data_)) then
     do i = 1, size(this%data_)
        if (associated(this%data_(i)%meanPressure)) deallocate(this%data_(i)%meanPressure)
        nullify(this%data_(i)%meanPressure)
     end do
  end if
  SAFE_DEALLOCATE(this%data_)

end subroutine cleanupAcousticNoise

function computeAcousticNoise(this, time, region) result(instantaneousFunctional)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use AcousticNoise_mod, only : t_AcousticNoise

  ! <<< Arguments >>>
  class(t_AcousticNoise) :: this
  real(SCALAR_KIND), intent(in) :: time
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, ierror
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
     assert(associated(this%data_(i)%meanPressure))
     assert(size(this%data_(i)%meanPressure, 1) == region%grids(i)%nGridPoints)
     assert(size(this%data_(i)%meanPressure, 2) == 1)

     allocate(F(region%grids(i)%nGridPoints, 1))
     F = region%states(i)%pressure - this%data_(i)%meanPressure
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

end function computeAcousticNoise

subroutine computeAcousticNoiseAdjointForcing(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use AcousticNoise_mod, only : t_AcousticNoise
  use CostTargetPatch_mod, only : t_CostTargetPatch

  ! <<< Private members >>>
  use AcousticNoiseImpl, only : computeAdjointForcingOnPatch

  implicit none

  ! <<< Arguments >>>
  class(t_AcousticNoise) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer :: i, j
  class(t_Patch), pointer :: patch => null()
  
  if (.not. allocated(region%patchFactories)) return

  do i = 1, size(region%grids)
     do j = 1, size(region%patchFactories)
        call region%patchFactories(j)%connect(patch)
        if (.not. associated(patch)) cycle
        if (patch%gridIndex /= region%grids(i)%index .or. patch%nPatchPoints <= 0) cycle
        call computeAdjointForcingOnPatch(patch, region%grids(i),                            &
             region%states(i)%velocity, region%states(i)%pressure(:,1),                      &
             this%data_(i)%meanPressure(:,1), region%solverOptions%ratioOfSpecificHeats)
     end do
  end do

end subroutine computeAcousticNoiseAdjointForcing

function isAcousticNoisePatchValid(this, patchDescriptor, gridSize, normalDirection,         &             
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

end function isAcousticNoisePatchValid
