#include "config.h"

subroutine setupFuelActuator(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use FuelActuator_mod, only : t_FuelActuator

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_FuelActuator) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, gradientBufferSize
  character(len = STRING_LENGTH) :: key, message, fuel
  class(t_Patch), pointer :: patch => null()

  call this%cleanup()

  if (region%simulationFlags%predictionOnly) return

  write(key, '(A)') "fuel_actuator/"
  fuel =  getOption(trim(key) // "reactant","H2")
  select case (trim(fuel))
  case ("H2")
     this%fuelIndex = region%states(1)%combustion%H2
  case ("O2")
     this%fuelIndex = region%states(1)%combustion%O2
  case default
     write(message, '(A)') "WARNING, unknown fuel!"
     call gracefulExit(region%comm, message)
  end select

  assert(this%fuelIndex > 0)
  assert(this%fuelIndex <= region%solverOptions%nSpecies)

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

        end select
     end do
  end if

end subroutine setupFuelActuator

subroutine cleanupFuelActuator(this)

  ! <<< Derived types >>>
  use FuelActuator_mod, only : t_FuelActuator

  implicit none

  ! <<< Arguments >>>
  class(t_FuelActuator) :: this

  call this%cleanupBase()

end subroutine cleanupFuelActuator

function computeFuelActuatorSensitivity(this, region) result(instantaneousSensitivity)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use FuelActuator_mod, only : t_FuelActuator

  implicit none

  ! <<< Arguments >>>
  class(t_FuelActuator) :: this
  class(t_Region) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousSensitivity

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions, ierror
  SCALAR_TYPE, allocatable :: F(:,:)

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(size(region%grids) == size(region%states))
  assert(region%solverOptions%nSpecies > 0)

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

     allocate(F(region%grids(i)%nGridPoints, 1))
     F(:,1) = region%states(i)%adjointVariables(:,nDimensions+2+this%fuelIndex) *            &
          region%grids(i)%controlMollifier(:,1)
     instantaneousSensitivity = instantaneousSensitivity +                                   &
          region%grids(i)%computeInnerProduct(F, F)
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

end function computeFuelActuatorSensitivity

subroutine updateFuelActuatorForcing(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use FuelActuator_mod, only : t_FuelActuator

  implicit none

  ! <<< Arguments >>>
  class(t_FuelActuator) :: this
  class(t_Region) :: region

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

           patch%controlForcing(:,:) = 0.0_wp
           patch%controlForcing(:,nDimensions+2+this%fuelIndex) =                            &
                - region%states(j)%actuationAmount *                                         &
                patch%gradientBuffer(:,1,patch%iGradientBuffer)

           if (patch%iGradientBuffer == 1)                                                   &
                patch%iGradientBuffer = size(patch%gradientBuffer, 3) + 1

        end select
     end do
  end do

end subroutine updateFuelActuatorForcing

subroutine updateFuelActuatorGradient(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use FuelActuator_mod, only : t_FuelActuator

  implicit none

  ! <<< Arguments >>>
  class(t_FuelActuator) :: this
  class(t_Region), intent(in) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions
  class(t_Patch), pointer :: patch => null()
  SCALAR_TYPE, allocatable :: F(:,:)

  if (.not. allocated(region%patchFactories)) return

  nDimensions = size(region%globalGridSizes, 1)
  assert_key(nDimensions, (1, 2, 3))

  do i = 1, size(region%patchFactories)
     call region%patchFactories(i)%connect(patch)
     if (.not. associated(patch)) cycle
     do j = 1, size(region%states)
        if (patch%gridIndex /= region%grids(j)%index .or. patch%nPatchPoints <= 0) cycle
        select type (patch)
        class is (t_ActuatorPatch)

           patch%iGradientBuffer = patch%iGradientBuffer + 1
           assert(patch%iGradientBuffer >= 1)
           assert(patch%iGradientBuffer <= size(patch%gradientBuffer, 3))

           allocate(F(patch%nPatchPoints, 2))
           call patch%collect(region%states(j)%adjointVariables                              &
                (:,nDimensions+2+this%fuelIndex), F(:,1))
           call patch%collect(region%grids(j)%controlMollifier(:,1), F(:,2))
           patch%gradientBuffer(:,1,patch%iGradientBuffer) = product(F, dim = 2)
           SAFE_DEALLOCATE(F)

           if (patch%iGradientBuffer == size(patch%gradientBuffer, 3)) then
              call patch%saveGradient()
              patch%iGradientBuffer = 0
           end if

        end select
     end do
  end do

end subroutine updateFuelActuatorGradient

function isFuelActuatorPatchValid(this, patchDescriptor, gridSize,                        &
     normalDirection, extent, simulationFlags, message) result(isPatchValid)

  ! <<< Derived types >>>
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags
  use FuelActuator_mod, only : t_FuelActuator

  implicit none

  ! <<< Arguments >>>
  class(t_FuelActuator) :: this
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

end function isFuelActuatorPatchValid

subroutine hookFuelActuatorBeforeTimemarch(this, region, mode)

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
  integer :: i, stat, fileUnit, mpiFileHandle, procRank, ierror
  class(t_Patch), pointer :: patch => null()
  logical :: fileExists
  character(len = STRING_LENGTH) :: message

  if (.not. allocated(region%patchFactories)) return

  do i = 1, size(region%patchFactories)
     call region%patchFactories(i)%connect(patch)
     if (.not. associated(patch)) cycle
     select type (patch)
     class is (t_ActuatorPatch)
        if (patch%comm == MPI_COMM_NULL) cycle

        call MPI_Comm_rank(patch%comm, procRank, ierror)

        select case (mode)

        case (FORWARD)
           if (procRank == 0) inquire(file = trim(patch%gradientFilename), exist = fileExists)
           call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, patch%comm, ierror)
           if (.not. fileExists) then
              write(message, '(3A,I0.0,A)') "No gradient information available for patch '", &
                   trim(patch%name), "' on grid ", patch%gridIndex, "!"
              call gracefulExit(patch%comm, message)
           end if
           patch%iGradientBuffer = size(patch%gradientBuffer, 3) + 1
           call MPI_File_open(patch%comm, trim(patch%gradientFilename) // char(0),           &
                MPI_MODE_WRONLY, MPI_INFO_NULL, mpiFileHandle, ierror)
           call MPI_File_get_size(mpiFileHandle, patch%gradientFileOffset, ierror)
           call MPI_File_close(mpiFileHandle, ierror)

        case (ADJOINT)
           if (procRank == 0) then
              open(unit = getFreeUnit(fileUnit), file = trim(patch%gradientFilename),        &
                   iostat = stat, status = 'old')
              if (stat == 0) close(fileUnit, status = 'delete')
              open(unit = getFreeUnit(fileUnit), file = trim(patch%gradientFilename),        &
                   action = 'write', status = 'unknown')
              close(fileUnit)
           end if
           call MPI_Barrier(patch%comm, ierror)
           patch%iGradientBuffer = 0
           patch%gradientFileOffset = int(0, MPI_OFFSET_KIND)

        end select

     end select
  end do

end subroutine hookFuelActuatorBeforeTimemarch

subroutine hookFuelActuatorAfterTimemarch(this, region, mode)

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

  do i = 1, size(region%patchFactories)
     call region%patchFactories(i)%connect(patch)
     if (.not. associated(patch)) cycle
     select type (patch)
     class is (t_ActuatorPatch)
        if (patch%comm == MPI_COMM_NULL) cycle

        call MPI_Comm_rank(patch%comm, procRank, ierror)

        select case (mode)

        case (ADJOINT)
           call patch%saveGradient()

        end select

     end select
  end do

end subroutine hookFuelActuatorAfterTimemarch
