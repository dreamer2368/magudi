#include "config.h"

subroutine setupIgnitionActuator(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use IgnitionActuator_mod, only : t_IgnitionActuator

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionActuator) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, gradientBufferSize
  class(t_Patch), pointer :: patch => null()

  call this%cleanup()

  if (region%simulationFlags%predictionOnly) return

  call this%setupBase(region%simulationFlags, region%solverOptions)

  gradientBufferSize = getOption("gradient_buffer_size", 1)

  if (allocated(region%patchFactories)) then
     do i = 1, size(region%patchFactories)
        call region%patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        select type (patch)
           class is (t_ActuatorPatch)
           if (patch%nPatchPoints <= 0) cycle

           SAFE_DEALLOCATE(patch%gradientBuffer)
           allocate(patch%gradientBuffer(patch%nPatchPoints, 1, gradientBufferSize))
           patch%gradientBuffer = 0.0_wp

           ! Control forcing not yet implemented.
           patch%controlForcing = 0.0_wp

        end select
     end do
  end if

  this%partialSensitivity = getOption("partial_sensitivity",.false.)
  if (this%partialSensitivity) then
     call getRequiredOption("sensitivity_dependence", this%sensitivityDependence,            &
          region%comm)
  end if

end subroutine setupIgnitionActuator

subroutine cleanupIgnitionActuator(this)

  ! <<< Derived types >>>
  use IgnitionActuator_mod, only : t_IgnitionActuator

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionActuator) :: this

  call this%cleanupBase()

end subroutine cleanupIgnitionActuator

function computeIgnitionActuatorSensitivity(this, region) result(instantaneousSensitivity)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use IgnitionActuator_mod, only : t_IgnitionActuator
  use IgnitionSource_mod, only : t_IgnitionSource

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionActuator) :: this
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousSensitivity

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions, ierror
  SCALAR_TYPE, allocatable :: F(:,:), ignitionForce(:,:)

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(size(region%grids) == size(region%states))

  instantaneousSensitivity = 0.0_wp

  do i = 1, size(region%grids)

     nDimensions = region%grids(i)%nDimensions
     assert_key(nDimensions, (1, 2, 3))

     assert(region%grids(i)%nGridPoints > 0)
     assert(allocated(region%grids(i)%controlMollifier))
     assert(size(region%grids(i)%controlMollifier, 1) == region%grids(i)%nGridPoints)
     assert(size(region%grids(i)%controlMollifier, 2) == 1)
     assert(allocated(region%states(i)%adjointVariables))
     assert(size(region%states(i)%adjointVariables, 1) == region%grids(i)%nGridPoints)
     assert(size(region%states(i)%adjointVariables, 2) >= nDimensions + 2)

     allocate(F(region%grids(i)%nGridPoints, 2))

     if (this%partialSensitivity) then

        assert(size(region%states(i)%ignitionSources) > 0)

        allocate(ignitionForce(region%grids(i)%nGridPoints, nDimensions+2))
        ignitionForce = 0.0_wp

        call region%states(i)%ignitionSources(1)%add(region%states(i)%time,                  &
             region%grids(i)%coordinates, region%grids(i)%iblank,                            &
             region%solverOptions%ratioOfSpecificHeats,                                      &
             region%states(i)%combustion%heatRelease, ignitionForce)

        select case (this%sensitivityDependence)

        case ('AMPLITUDE')
           F(:,2) = ignitionForce(:,nDimensions+2) /                                         &
                region%states(i)%ignitionSources(1)%amplitude

        case ('VERTICAL_POSITION')

           F(:,2) = ignitionForce(:,nDimensions+2) *                                         &
                (region%grids(i)%coordinates(:,2) -                                          &
                region%states(i)%ignitionSources(1)%location(2)) /                           &
                region%states(i)%ignitionSources(1)%radius**2

        case ('INITIAL_TIME')

           F(:,2) = ignitionForce(:,nDimensions+2) *                                         &
                (region%states(i)%time - region%states(i)%ignitionSources(1)%timeStart) /    &
                region%states(i)%ignitionSources(1)%timeDuration**2

        end select

        F(:,1) = region%states(i)%adjointVariables(:,nDimensions+2)
        instantaneousSensitivity = instantaneousSensitivity +                                &
             region%grids(i)%computeInnerProduct(F(:,1), F(:,2),                             &
             region%grids(i)%controlMollifier(:,1))

        SAFE_DEALLOCATE(ignitionForce)

     else

        F(:,1) = region%states(i)%adjointVariables(:,nDimensions+2)
        F(:,2) = region%grids(i)%controlMollifier(:,1)
        instantaneousSensitivity = instantaneousSensitivity +                                &
             region%grids(i)%computeInnerProduct(F(:,1), F(:,2))

     end if

     SAFE_DEALLOCATE(F)

  end do

  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE, instantaneousSensitivity, 1,                         &
       SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(instantaneousSensitivity, 1, SCALAR_TYPE_MPI,                            &
          0, region%grids(i)%comm, ierror)
  end do

  this%cachedValue = instantaneousSensitivity

end function computeIgnitionActuatorSensitivity

subroutine updateIgnitionActuatorForcing(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use IgnitionActuator_mod, only : t_IgnitionActuator

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionActuator) :: this
  class(t_Region), intent(in) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions
  class(t_Patch), pointer :: patch => null()

  if (.not. allocated(region%patchFactories)) return

  nDimensions = size(region%globalGridSizes, 1)
  assert_key(nDimensions, (1, 2, 3))

  do i = 1, size(region%patchFactories)
     call region%patchFactories(i)%connect(patch)
     if (.not. associated(patch)) cycle
     do j = 1, size(region%states)
        if (patch%gridIndex /= region%grids(j)%index) cycle
        select type (patch)
        class is (t_ActuatorPatch)

           patch%iGradientBuffer = patch%iGradientBuffer - 1

           assert(patch%iGradientBuffer >= 1)
           assert(patch%iGradientBuffer <= size(patch%gradientBuffer, 3))

           if (patch%iGradientBuffer == size(patch%gradientBuffer, 3))                       &
                call patch%loadGradient()

           ! No forcing for now.
           patch%controlForcing = 0.0_wp

           if (patch%iGradientBuffer == 1)                                                   &
                patch%iGradientBuffer = size(patch%gradientBuffer, 3) + 1

        end select
     end do
  end do

end subroutine updateIgnitionActuatorForcing

subroutine updateIgnitionActuatorGradient(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use IgnitionActuator_mod, only : t_IgnitionActuator

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionActuator) :: this
  class(t_Region), intent(in) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions
  class(t_Patch), pointer :: patch => null()
  SCALAR_TYPE, allocatable :: F(:,:)

  if (.not. allocated(region%patchFactories)) return

  ! Nothing to do here.

end subroutine updateIgnitionActuatorGradient

function isIgnitionActuatorPatchValid(this, patchDescriptor, gridSize,                        &
     normalDirection, extent, simulationFlags, message) result(isPatchValid)

  ! <<< Derived types >>>
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags
  use IgnitionActuator_mod, only : t_IgnitionActuator

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionActuator) :: this
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

end function isIgnitionActuatorPatchValid

subroutine hookIgnitionActuatorBeforeTimemarch(this, region, mode)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use Controller_mod, only : t_Controller
  use ActuatorPatch_mod, only : t_ActuatorPatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use InputHelper, only : getFreeUnit
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_Controller) :: this
  class(t_Region) :: region
  integer, intent(in) :: mode

  ! <<< Local variables >>>
  integer :: i, fileUnit, mpiFileHandle, procRank, ierror
  class(t_Patch), pointer :: patch => null()
  logical :: fileExists
  character(len = STRING_LENGTH) :: message

  if (.not. allocated(region%patchFactories)) return

  ! Nothing to do here.

end subroutine hookIgnitionActuatorBeforeTimemarch

subroutine hookIgnitionActuatorAfterTimemarch(this, region, mode)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use Controller_mod, only : t_Controller
  use ActuatorPatch_mod, only : t_ActuatorPatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use InputHelper, only : getFreeUnit
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_Controller) :: this
  class(t_Region) :: region
  integer, intent(in) :: mode

  ! <<< Local variables >>>
  integer :: i, procRank, ierror
  class(t_Patch), pointer :: patch => null()

  if (.not. allocated(region%patchFactories)) return

  ! Nothing to do here.

end subroutine hookIgnitionActuatorAfterTimemarch
