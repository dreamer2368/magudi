#include "config.h"

subroutine setupIgnitionActuator(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use IgnitionActuator_mod, only : t_IgnitionActuator

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionActuator) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j
  character(len = STRING_LENGTH) :: key, message

  call this%cleanup()

  if (region%simulationFlags%predictionOnly) return

  call this%setupBase(region%simulationFlags, region%solverOptions)

  this%nParameters = getOption("number_of_parameters", 1)
  if (this%nParameters <= 0) then
     write(message, '(A)') "Number of parameters must be > 0!"
     call gracefulExit(region%comm, message)
  end if

  allocate(this%cachedValue(this%nParameters))
  allocate(this%runningTimeQuadrature(this%nParameters))
  this%cachedValue = 0.0_wp
  this%runningTimeQuadrature = 0.0_wp

  allocate(this%sensitivityParameter(this%nParameters))
  do i = 1, this%nParameters
     write(key, '(A,I2.2)') "sensitivity_parameter", i
     this%sensitivityParameter(i) = getOption(trim(key),trim(key))
  end do

  allocate(this%baselineValue(this%nParameters))
  this%baselineValue = 0.0_wp

  do i = 1, size(region%grids)

     if (.not.allocated(region%states(i)%ignitionSources)) then
        write(message, '(A)') "Ignition actuator requires ignition source!"
        call gracefulExit(region%comm, message)
     end if

     do j = 1, this%nParameters

        select case (trim(this%sensitivityParameter(j)))

        case ('AMPLITUDE')
           this%baselineValue(j) = region%states(i)%ignitionSources(1)%amplitude

        case ('POSITION_X')
           this%baselineValue(j) = region%states(i)%ignitionSources(1)%location(1)

        case ('POSITION_Y')
           this%baselineValue(j) = region%states(i)%ignitionSources(1)%location(2)

        case ('POSITION_Z')
           this%baselineValue(j) = region%states(i)%ignitionSources(1)%location(3)

        case ('RADIUS_X')
           this%baselineValue(j) = region%states(i)%ignitionSources(1)%radius(1)

        case ('RADIUS_Y')
           this%baselineValue(j) = region%states(i)%ignitionSources(1)%radius(2)

        case ('RADIUS_Z')
           this%baselineValue(j) = region%states(i)%ignitionSources(1)%radius(3)

        case ('INITIAL_TIME')
           this%baselineValue(j) = region%states(i)%ignitionSources(1)%timeStart

        case ('DURATION')
           this%baselineValue(j) = region%states(i)%ignitionSources(1)%timeDuration

        case ('DAMKOHLER')
           this%baselineValue(j) = region%states(i)%combustion%Damkohler

        case ('ZEL_DOVICH')
           this%baselineValue(j) = region%states(i)%combustion%zelDovich

        case default
           write(message, '(A)') "Unknown sensitivity parameter!"
           call gracefulExit(region%comm, message)

        end select

     end do

  end do

end subroutine setupIgnitionActuator

subroutine cleanupIgnitionActuator(this)

  ! <<< Derived types >>>
  use IgnitionActuator_mod, only : t_IgnitionActuator

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionActuator) :: this

  call this%cleanupBase()

  SAFE_DEALLOCATE(this%baselineValue)
  SAFE_DEALLOCATE(this%sensitivityParameter)

end subroutine cleanupIgnitionActuator

subroutine computeIgnitionActuatorSensitivity(this, region)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use IgnitionActuator_mod, only : t_IgnitionActuator

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionActuator) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, nSpecies, nUnknowns, ierror
  real(SCALAR_KIND) ::  time, timeStart, timeDuration, amplitude, radius(3), location(3)
  SCALAR_TYPE, allocatable :: F(:,:), ignitionSource(:,:), combustionSource(:,:),            &
       instantaneousSensitivity(:)

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(size(region%grids) == size(region%states))
  assert(this%nParameters > 0)
  assert(size(this%cachedValue) == this%nParameters)

  allocate(instantaneousSensitivity(this%nParameters))
  instantaneousSensitivity = 0.0_wp

  do i = 1, size(region%states)

     nDimensions = region%grids(i)%nDimensions
     nspecies = region%solverOptions%nSpecies
     nUnknowns = region%solverOptions%nUnknowns
     time = region%states(i)%adjointCoefficientTime
     timeStart = region%states(i)%ignitionSources(1)%timeStart
     timeDuration = region%states(i)%ignitionSources(1)%timeDuration
     amplitude = region%states(i)%ignitionSources(1)%amplitude
     radius = region%states(i)%ignitionSources(1)%radius
     location = region%states(i)%ignitionSources(1)%location

     assert_key(nDimensions, (1, 2, 3))
     assert(nSpecies >= 0)
     assert(nUnknowns == nDimensions + 2 + nSpecies)
     assert(region%grids(i)%nGridPoints > 0)
     assert(size(region%states(i)%ignitionSources) > 0)
     assert(allocated(region%states(i)%adjointVariables))
     assert(size(region%states(i)%adjointVariables, 1) == region%grids(i)%nGridPoints)
     assert(size(region%states(i)%adjointVariables, 2) >= nDimensions + 2)
     assert(region%states(i)%combustion%Damkohler > 0.0_wp)

     allocate(F(region%grids(i)%nGridPoints, 2))
     allocate(ignitionSource(region%grids(i)%nGridPoints, nUnknowns))
     allocate(combustionSource(region%grids(i)%nGridPoints, nUnknowns))

     ignitionSource = 0.0_wp
     combustionSource = 0.0_wp

     call region%states(i)%ignitionSources(1)%add(time, region%grids(i)%coordinates,         &
          region%grids(i)%iblank, region%solverOptions%ratioOfSpecificHeats,                 &
          region%states(i)%combustion%heatRelease, ignitionSource)

     call region%states(i)%combustion%addForward(nDimensions, nSpecies,                      &
          region%solverOptions%ratioOfSpecificHeats, region%states(i)%conservedVariables,    &
          region%states(i)%temperature(:,1), region%states(i)%massFraction,                  &
          region%grids(i)%iblank, combustionSource)

     ! Compute sensitivities of ignition parameters.
     F(:,1) = region%states(i)%adjointVariables(:,nDimensions+2)

     do j = 1, this%nParameters
        select case (trim(this%sensitivityParameter(j)))

        case ('AMPLITUDE')
           F(:,2) = ignitionSource(:,nDimensions+2) / amplitude

           instantaneousSensitivity(j) = instantaneousSensitivity(j) +                       &
                region%grids(i)%computeInnerProduct(F(:,1), F(:,2))

        case ('POSITION_X')
           F(:,2) = ignitionSource(:,nDimensions+2) *                                        &
                (region%grids(i)%coordinates(:,1) - location(1)) / radius(1)**2

           instantaneousSensitivity(j) = instantaneousSensitivity(j) +                       &
                region%grids(i)%computeInnerProduct(F(:,1), F(:,2))

        case ('POSITION_Y')
           F(:,2) = ignitionSource(:,nDimensions+2) *                                        &
                (region%grids(i)%coordinates(:,2) - location(2)) / radius(2)**2

           instantaneousSensitivity(j) = instantaneousSensitivity(j) +                       &
                region%grids(i)%computeInnerProduct(F(:,1), F(:,2))

        case ('RADIUS_X')
           F(:,2) = ignitionSource(:,nDimensions+2) *                                        &
                (region%grids(i)%coordinates(:,1) - location(1))**2 / radius(1)**3

           instantaneousSensitivity(j) = instantaneousSensitivity(j) +                       &
                region%grids(i)%computeInnerProduct(F(:,1), F(:,2))

        case ('RADIUS_Y')
           F(:,2) = ignitionSource(:,nDimensions+2) *                                        &
                (region%grids(i)%coordinates(:,2) - location(2))**2 / radius(2)**3

           instantaneousSensitivity(j) = instantaneousSensitivity(j) +                       &
                region%grids(i)%computeInnerProduct(F(:,1), F(:,2))

        case ('INITIAL_TIME')
           F(:,2) = ignitionSource(:,nDimensions+2) * (time - timeStart) / timeDuration**2

           instantaneousSensitivity(j) = instantaneousSensitivity(j) +                       &
                region%grids(i)%computeInnerProduct(F(:,1), F(:,2))

        case ('DURATION')
           F(:,2) = - ignitionSource(:,nDimensions+2) * (timeDuration + time - timeStart) *  &
                (timeDuration - time + timeStart) / timeDuration**3

           instantaneousSensitivity(j) = instantaneousSensitivity(j) +                       &
                region%grids(i)%computeInnerProduct(F(:,1), F(:,2))

        end select

     end do

     ! Compute sensitivities of combustion parameters.
     do j = 1, this%nParameters
        select case (trim(this%sensitivityParameter(j)))

        case ('DAMKOHLER')
           F(:,1) = region%states(i)%adjointVariables(:,nDimensions+2)
           F(:,2) = combustionSource(:,nDimensions+2) /                                      &
                region%states(i)%combustion%Damkohler

           instantaneousSensitivity(j) = instantaneousSensitivity(j) +                       &
                region%grids(i)%computeInnerProduct(F(:,1), F(:,2))

           do k = 1, nSpecies
              F(:,1) = region%states(i)%adjointVariables(:,nDimensions+2+k)
              F(:,2) = combustionSource(:,nDimensions+2+k) /                                 &
                   region%states(i)%combustion%Damkohler

              instantaneousSensitivity(j) = instantaneousSensitivity(j) +                    &
                   region%grids(i)%computeInnerProduct(F(:,1), F(:,2))
           end do

        case ('ZEL_DOVICH')
           F(:,1) = region%states(i)%adjointVariables(:,nDimensions+2)
           F(:,2) = - combustionSource(:,nDimensions+2) /                                    &
                ( region%states(i)%combustion%heatRelease *                                  &
                region%states(i)%temperature(:,1) *                                          &
                (region%solverOptions%ratioOfSpecificHeats - 1.0_wp) *                       &
                (1.0_wp - region%states(i)%combustion%heatRelease) )

           instantaneousSensitivity(j) = instantaneousSensitivity(j) +                       &
                region%grids(i)%computeInnerProduct(F(:,1), F(:,2))

           do k = 1, nSpecies
              F(:,1) = region%states(i)%adjointVariables(:,nDimensions+2+k)
              F(:,2) = - combustionSource(:,nDimensions+2+k) /                               &
                   ( region%states(i)%combustion%heatRelease *                               &
                   region%states(i)%temperature(:,1) *                                       &
                   (region%solverOptions%ratioOfSpecificHeats - 1.0_wp) *                    &
                   (1.0_wp - region%states(i)%combustion%heatRelease) )

              instantaneousSensitivity(j) = instantaneousSensitivity(j) +                    &
                   region%grids(i)%computeInnerProduct(F(:,1), F(:,2))
           end do

        end select

     end do

     ! Clean up.
     SAFE_DEALLOCATE(F)
     SAFE_DEALLOCATE(ignitionSource)
     SAFE_DEALLOCATE(combustionSource)

  end do

  if (region%commGridMasters /= MPI_COMM_NULL) then
     call MPI_Allreduce(MPI_IN_PLACE, instantaneousSensitivity, this%nParameters,            &
          SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)
  end if

  do i = 1, size(region%grids)
     call MPI_Bcast(instantaneousSensitivity, this%nParameters, SCALAR_TYPE_MPI,             &
          0, region%grids(i)%comm, ierror)
  end do

  this%cachedValue = instantaneousSensitivity
  SAFE_DEALLOCATE(instantaneousSensitivity)

end subroutine computeIgnitionActuatorSensitivity

subroutine updateIgnitionActuatorForcing(this, region)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use IgnitionActuator_mod, only : t_IgnitionActuator

  implicit none

  ! <<< Arguments >>>
  class(t_IgnitionActuator) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions, nUnknowns, nSpecies

  if (.not. allocated(region%patchFactories)) return

  nDimensions = size(region%globalGridSizes, 1)
  assert_key(nDimensions, (1, 2, 3))

  nSpecies = region%solverOptions%nSpecies
  assert(nSpecies >= 0)

  nUnknowns = region%solverOptions%nUnknowns
  assert(nUnknowns == nDimensions + 2 + nSpecies)

  do i = 1, size(region%states)
     assert(abs(region%states(i)%gradientDirection) == 1)
     do j = 1, this%nParameters
        select case (trim(this%sensitivityParameter(j)))

        case ('AMPLITUDE')
           region%states(i)%ignitionSources(1)%amplitude = this%baselineValue(j) +           &
                real(region%states(i)%gradientDirection, wp) *                               &
                region%states(i)%actuationAmount * region%states(i)%controlGradient(j)

        case ('POSITION_X')
           region%states(i)%ignitionSources(1)%location(1) = this%baselineValue(j) +         &
                real(region%states(i)%gradientDirection, wp) *                               &
                region%states(i)%actuationAmount * region%states(i)%controlGradient(j)

        case ('POSITION_Y')
           region%states(i)%ignitionSources(1)%location(2) = this%baselineValue(j) +         &
                real(region%states(i)%gradientDirection, wp) *                               &
                region%states(i)%actuationAmount * region%states(i)%controlGradient(j)

        case ('POSITION_Z')
           region%states(i)%ignitionSources(1)%location(3) = this%baselineValue(j) +         &
                real(region%states(i)%gradientDirection, wp) *                               &
                region%states(i)%actuationAmount * region%states(i)%controlGradient(j)

        case ('RADIUS_X')
           region%states(i)%ignitionSources(1)%radius(1) = this%baselineValue(j) +           &
                real(region%states(i)%gradientDirection, wp) *                               &
                region%states(i)%actuationAmount * region%states(i)%controlGradient(j)

        case ('RADIUS_Y')
           region%states(i)%ignitionSources(1)%radius(2) = this%baselineValue(j) +           &
                real(region%states(i)%gradientDirection, wp) *                               &
                region%states(i)%actuationAmount * region%states(i)%controlGradient(j)

        case ('RADIUS_Z')
           region%states(i)%ignitionSources(1)%radius(3) = this%baselineValue(j) +           &
                real(region%states(i)%gradientDirection, wp) *                               &
                region%states(i)%actuationAmount * region%states(i)%controlGradient(j)

        case ('INITIAL_TIME')
           region%states(i)%ignitionSources(1)%timeStart = this%baselineValue(j) +           &
                real(region%states(i)%gradientDirection, wp) *                               &
                region%states(i)%actuationAmount * region%states(i)%controlGradient(j)

        case ('DURATION')
           region%states(i)%ignitionSources(1)%timeDuration = this%baselineValue(j) +        &
                real(region%states(i)%gradientDirection, wp) *                               &
                region%states(i)%actuationAmount * region%states(i)%controlGradient(j)

        case ('DAMKOHLER')
           region%states(i)%combustion%Damkohler = this%baselineValue(j) +                   &
                real(region%states(i)%gradientDirection, wp) *                               &
                region%states(i)%actuationAmount * region%states(i)%controlGradient(j)

        case ('ZEL_DOVICH')
           region%states(i)%combustion%zelDovich = this%baselineValue(j) +                   &
                real(region%states(i)%gradientDirection, wp) *                               &
                region%states(i)%actuationAmount * region%states(i)%controlGradient(j)

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
