#include "config.h"

module IgnitionActuatorImpl

  implicit none
  public

contains

  subroutine computeSource(time, coordinates, iblank, timeStart, timeDuration, amplitude,    &
       radius, location, ratioOfSpecificHeats, heatRelease, ignitionSource)

    implicit none

    ! <<< Arguments >>>
    real(SCALAR_KIND), intent(in) :: time, timeStart, timeDuration, amplitude, radius,       &
         location(:), ratioOfSpecificHeats, heatRelease
    SCALAR_TYPE, intent(in) :: coordinates(:,:)
    integer, intent(in) :: iblank(:)
    SCALAR_TYPE, intent(out) :: ignitionSource(:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, nDimensions
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
    real(wp) :: r, power, timePortion, referenceTemperature, flameTemperature, gaussianFactor

    nDimensions = size(coordinates,2)
    assert_key(nDimensions, (1, 2, 3))
    assert(size(location) >= nDimensions)
    assert(size(coordinates,1) == size(ignitionSource))

    timePortion = exp( -0.5_wp * (time - timeStart)**2 / timeDuration**2)

    referenceTemperature = 1.0_wp / (ratioOfSpecificHeats - 1.0_wp)

    flameTemperature = referenceTemperature / (1.0_wp - heatRelease)

    power = 0.5_wp * amplitude * heatRelease * flameTemperature / timeDuration /             &
         sqrt(2.0_wp * pi)

    gaussianFactor = 0.5_wp / radius**2

    do i = 1, size(ignitionSource)

       if (iblank(i) == 0) cycle

       r = real(sum((coordinates(i,:) - location(1:nDimensions)) ** 2), wp)

       ignitionSource(i) = power * timePortion! * exp(- gaussianFactor * r)

    end do

  end subroutine computeSource

end module IgnitionActuatorImpl


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
  character(len = STRING_LENGTH) :: key

  call this%cleanup()

  if (region%simulationFlags%predictionOnly) return

  call this%setupBase(region%simulationFlags, region%solverOptions)

  write(key, '(A)') "ignition_actuator/"

  call getRequiredOption(trim(key) // "sensitivity_dependence",                              &
       this%sensitivityDependence, region%comm)

  this%location(1) =  getOption(trim(key) // "x",0.0_wp)
  this%location(2) =  getOption(trim(key) // "y",0.0_wp)
  this%location(3) =  getOption(trim(key) // "z",0.0_wp)
  call getRequiredOption(trim(key) // "amplitude", this%amplitude, region%comm)
  call getRequiredOption(trim(key) // "radius", this%radius, region%comm)
  call getRequiredOption(trim(key) // "time_start", this%timeStart, region%comm)
  call getRequiredOption(trim(key) // "time_duration", this%timeDuration, region%comm)

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

  ! <<< Private members >>>
  use IgnitionActuatorImpl, only : computeSource

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionActuator) :: this
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousSensitivity

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions, ierror
  real(SCALAR_KIND) ::  timeStart, timeDuration, amplitude, radius, location(3)
  SCALAR_TYPE, allocatable :: F(:,:), ignitionSource(:)

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(size(region%grids) == size(region%states))

  instantaneousSensitivity = 0.0_wp

  timeStart = this%timeStart
  timeDuration = this%timeDuration
  amplitude = this%amplitude
  radius = this%radius
  location = this%location

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
     allocate(ignitionSource(region%grids(i)%nGridPoints))

     F(:,1) = region%states(i)%adjointVariables(:,nDimensions+2)

     call computeSource(region%states(i)%time, region%grids(i)%coordinates,                  &
          region%grids(i)%iblank, timeStart, timeDuration, amplitude, radius, location,      &
          region%solverOptions%ratioOfSpecificHeats,                                         &
          region%states(i)%combustion%heatRelease, ignitionSource)

     select case (this%sensitivityDependence)

     case ('AMPLITUDE')
        F(:,2) = ignitionSource / this%amplitude

     case ('VERTICAL_POSITION')

        F(:,2) = ignitionSource *                                                            &
             (region%grids(i)%coordinates(:,2) - this%location(2)) /  this%radius**2

     case ('INITIAL_TIME')

        F(:,2) = ignitionSource *                                                            &
             (region%states(i)%time - this%timeStart) / this%timeDuration**2

     end select

     instantaneousSensitivity = instantaneousSensitivity +                                   &
          region%grids(i)%computeInnerProduct(F(:,1), F(:,2),                                &
          region%grids(i)%controlMollifier(:,1))

     SAFE_DEALLOCATE(ignitionSource)
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

  ! <<< Private members >>>
  use IgnitionActuatorImpl, only : computeSource

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionActuator) :: this
  class(t_Region), intent(in) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions
  real(SCALAR_KIND) ::  timeStart, timeDuration, amplitude, radius, location(3)
  SCALAR_TYPE, allocatable :: ignitionSource(:)
  class(t_Patch), pointer :: patch => null()

  if (.not. allocated(region%patchFactories)) return

  nDimensions = size(region%globalGridSizes, 1)
  assert_key(nDimensions, (1, 2, 3))

  timeStart = this%timeStart
  timeDuration = this%timeDuration
  amplitude = this%amplitude
  radius = this%radius
  location = this%location

  do i = 1, size(region%patchFactories)
     call region%patchFactories(i)%connect(patch)
     if (.not. associated(patch)) cycle
     do j = 1, size(region%states)
        if (patch%gridIndex /= region%grids(j)%index) cycle
        select type (patch)
        class is (t_ActuatorPatch)

           allocate(ignitionSource(region%grids(j)%nGridPoints))

           select case (this%sensitivityDependence)

           case ('AMPLITUDE')
              amplitude = - region%states(j)%actuationAmount *                               &
                   region%states(j)%costSensitivity

           case ('VERTICAL_POSITION')

              location(2) = - region%states(j)%actuationAmount *                             &
                   region%states(j)%costSensitivity

           case ('INITIAL_TIME')

              timeStart = - region%states(j)%actuationAmount *                               &
                   region%states(j)%costSensitivity

           end select

           call computeSource(region%states(j)%time, region%grids(j)%coordinates,            &
                region%grids(j)%iblank, timeStart, timeDuration, amplitude, radius, location,&
                region%solverOptions%ratioOfSpecificHeats,                                   &
                region%states(j)%combustion%heatRelease, ignitionSource)

           patch%controlForcing(:,:) = 0.0_wp
           patch%controlForcing(:,nDimensions+2) = ignitionSource

           SAFE_DEALLOCATE(ignitionSource)

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

  if (.not. allocated(region%patchFactories)) return

  ! Nothing to do here.

end subroutine hookIgnitionActuatorAfterTimemarch
