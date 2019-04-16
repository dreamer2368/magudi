#include "config.h"

subroutine setupThermalActuator(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use ThermalActuator_mod, only : t_ThermalActuator

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_ThermalActuator) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: key
  integer :: i
  class(t_Patch), pointer :: patch => null()

  call this%cleanup()

  if (.not. region%simulationFlags%enableController) return

  call this%setupBase(region%simulationFlags, region%solverOptions)

  write(key, '(A)') "thermal_actuator/"

  this%useTimeRamp = getOption(trim(key) // "use_time_ramp", .false.)

  if (this%useTimeRamp) then
     call getRequiredOption(trim(key) // "ramp_width", this%rampWidthInverse)
     this%rampWidthInverse = 1.0_wp / this%rampWidthInverse
     call getRequiredOption(trim(key) // "ramp_offset", this%rampOffset)
     this%rampOffset = 1.0_wp - 0.5_wp * this%rampOffset
  end if

  this%controllerBufferSize = getOption("controller_buffer_size", 1)
  this%controllerSwitch = getOption("controller_switch", .false.)

  if (allocated(region%patchFactories)) then
     do i = 1, size(region%patchFactories)
        call region%patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle
        select type (patch)
           class is (t_ActuatorPatch)
           if (patch%nPatchPoints <= 0) cycle

           if (region%simulationFlags%enableAdjoint) then
             SAFE_DEALLOCATE(patch%gradientBuffer)
             allocate(patch%gradientBuffer(patch%nPatchPoints, 1, this%controllerBufferSize))
             patch%gradientBuffer = 0.0_wp
           end if
           SAFE_DEALLOCATE(patch%controlForcingBuffer)
           allocate(patch%controlForcingBuffer(patch%nPatchPoints, 1, this%controllerBufferSize))
           patch%controlForcingBuffer = 0.0_wp

        end select
     end do
  end if

end subroutine setupThermalActuator

subroutine cleanupThermalActuator(this)

  ! <<< Derived types >>>
  use ThermalActuator_mod, only : t_ThermalActuator

  implicit none

  ! <<< Arguments >>>
  class(t_ThermalActuator) :: this

  call this%cleanupBase()

end subroutine cleanupThermalActuator

function computeThermalActuatorSensitivity(this, region) result(instantaneousSensitivity)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use ThermalActuator_mod, only : t_ThermalActuator

  ! <<< SeungWhan: debug:time_ramp printing >>>
  use ErrorHandler, only : writeAndFlush
  use, intrinsic :: iso_fortran_env, only : output_unit

  implicit none

  ! <<< Arguments >>>
  class(t_ThermalActuator) :: this
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousSensitivity

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions, ierror
  SCALAR_TYPE, allocatable :: F(:,:)
  real(SCALAR_KIND) :: timeRampFactor

  ! <<< SeungWhan: variable for time_ramp printing >>>
  character(len = STRING_LENGTH) :: message

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(size(region%grids) == size(region%states))

  instantaneousSensitivity = 0.0_wp

  timeRampFactor = 1.0_wp
  if (this%useTimeRamp)                                                                      &
       timeRampFactor = this%rampFunction(2.0_wp * (region%states(1)%time -                  &
       this%onsetTime) / this%duration - 1.0_wp, this%rampWidthInverse, this%rampOffset)

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
     F(:,1) = region%states(i)%adjointVariables(:,nDimensions+2) *                           &
          region%grids(i)%controlMollifier(:,1) * timeRampFactor
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

end function computeThermalActuatorSensitivity

subroutine updateThermalActuatorForcing(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use ThermalActuator_mod, only : t_ThermalActuator

  ! <<< SeungWhan: debug:time_ramp printing >>>
  use ErrorHandler, only : writeAndFlush
  use, intrinsic :: iso_fortran_env, only : output_unit

  implicit none

  ! <<< Arguments >>>
  class(t_ThermalActuator) :: this
  class(t_Region), intent(in) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions
  real(SCALAR_KIND) :: timeRampFactor
  class(t_Patch), pointer :: patch => null()
  SCALAR_TYPE, allocatable :: temp(:,:)

  ! <<< SeungWhan: variable for time_ramp printing >>>
  character(len = STRING_LENGTH) :: message

  if (.not. allocated(region%patchFactories)) return

  nDimensions = size(region%globalGridSizes, 1)
  assert_key(nDimensions, (1, 2, 3))

  timeRampFactor = 0.0_wp
  if (region%states(1)%time>=this%onsetTime .and.                                            &
      region%states(1)%time<=this%onsetTime+this%duration) timeRampFactor = 1.0_wp
  if (this%useTimeRamp)                                                                      &
       timeRampFactor = this%rampFunction(2.0_wp * (region%states(1)%time -                  &
       this%onsetTime) / this%duration - 1.0_wp, this%rampWidthInverse, this%rampOffset)

  do i = 1, size(region%patchFactories)
     call region%patchFactories(i)%connect(patch)
     if (.not. associated(patch)) cycle
     do j = 1, size(region%states)
        if (patch%gridIndex /= region%grids(j)%index) cycle
        select type (patch)
        class is (t_ActuatorPatch)

           !SeungWhan: allocate temp = controlForcing
           allocate(temp,MOLD=patch%controlForcing)

           patch%iControlForcingBuffer = patch%iControlForcingBuffer - 1

           assert(patch%iControlForcingBuffer >= 1)
           assert(patch%iControlForcingBuffer <= size(patch%controlForcingBuffer, 3))

           if (patch%iControlForcingBuffer == size(patch%controlForcingBuffer, 3))                       &
                call patch%loadForcing()

           temp(:,1:nDimensions+1) = 0.0_wp
           temp(:,nDimensions+2) = patch%controlForcingBuffer(:,1,patch%iControlForcingBuffer)
           patch%controlForcing = patch%controlForcing + timeRampFactor * temp

           deallocate(temp)

           if (patch%iControlForcingBuffer == 1)                                                   &
                patch%iControlForcingBuffer = size(patch%controlForcingBuffer, 3) + 1

        end select
     end do
  end do

end subroutine updateThermalActuatorForcing

subroutine migrateToThermalActuatorForcing(this, region, startTimeStep, endTimeStep, nStages, iTimeStep, jSubStep)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use ThermalActuator_mod, only : t_ThermalActuator

  ! <<< SeungWhan: debug:time_ramp printing >>>
  use ErrorHandler, only : writeAndFlush
  use, intrinsic :: iso_fortran_env, only : output_unit

  implicit none

  ! <<< Arguments >>>
  class(t_ThermalActuator) :: this
  class(t_Region), intent(in) :: region
  integer, intent(in) :: startTimeStep, endTimeStep, nStages, iTimeStep, jSubStep

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions
  integer :: bufferIndex, bufferOffsetIndex, iControlForcingBuffer
  real(SCALAR_KIND) :: timeRampFactor
  class(t_Patch), pointer :: patch => null()
  SCALAR_TYPE, allocatable :: temp(:,:)

  ! <<< SeungWhan: variable for time_ramp printing >>>
  character(len = STRING_LENGTH) :: message

  if (.not. allocated(region%patchFactories)) return

  nDimensions = size(region%globalGridSizes, 1)
  assert_key(nDimensions, (1, 2, 3))

  timeRampFactor = 1.0_wp
  if (this%useTimeRamp)                                                                      &
       timeRampFactor = this%rampFunction(2.0_wp * (region%states(1)%time -                  &
       this%onsetTime) / this%duration - 1.0_wp, this%rampWidthInverse, this%rampOffset)

  do i = 1, size(region%patchFactories)
     call region%patchFactories(i)%connect(patch)
     if (.not. associated(patch)) cycle
     do j = 1, size(region%states)
        if (patch%gridIndex /= region%grids(j)%index) cycle
        select type (patch)
        class is (t_ActuatorPatch)

          ! Buffer index starts from end time, and starts from the beginning of the control forcing file.
          bufferIndex = (endTimeStep - iTimeStep)*nStages + (nStages-jSubStep)
          if (patch%adjointReferenceTimestep>0) bufferIndex = bufferIndex                     &
                                        + patch%adjointReferenceTimestep*nStages
          bufferOffsetIndex = bufferIndex/this%controllerBufferSize
          iControlForcingBuffer = MODULO( bufferIndex,this%controllerBufferSize ) + 1

           !SeungWhan: allocate temp = controlForcing
           allocate(temp,MOLD=patch%controlForcing)

           patch%iControlForcingBuffer = iControlForcingBuffer

           assert(patch%iControlForcingBuffer >= 1)
           assert(patch%iControlForcingBuffer <= size(patch%controlForcingBuffer, 3))

           if (patch%bufferOffsetIndex .ne. bufferOffsetIndex)                               &
                call patch%pinpointForcing(bufferOffsetIndex)

           temp(:,1:nDimensions+1) = 0.0_wp
           temp(:,nDimensions+2) = patch%controlForcingBuffer(:,1,patch%iControlForcingBuffer)
           patch%controlForcing = patch%controlForcing + timeRampFactor * temp

           deallocate(temp)

        end select
     end do
  end do

end subroutine migrateToThermalActuatorForcing

subroutine updateThermalActuatorGradient(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use ThermalActuator_mod, only : t_ThermalActuator

  implicit none

  ! <<< Arguments >>>
  class(t_ThermalActuator) :: this
  class(t_Region), intent(in) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions
  real(SCALAR_KIND) :: timeRampFactor
  class(t_Patch), pointer :: patch => null()
  SCALAR_TYPE, allocatable :: F(:,:)

  if (.not. allocated(region%patchFactories)) return

  nDimensions = size(region%globalGridSizes, 1)
  assert_key(nDimensions, (1, 2, 3))

  timeRampFactor = 1.0_wp
  if (this%useTimeRamp)                                                                      &
       timeRampFactor = this%rampFunction(2.0_wp * (region%states(1)%time -                  &
       this%onsetTime) / this%duration - 1.0_wp, this%rampWidthInverse, this%rampOffset)

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
           call patch%collect(region%states(j)%adjointVariables(:,nDimensions+2), F(:,1))
           call patch%collect(region%grids(j)%controlMollifier(:,1), F(:,2))
           F(:,2) = F(:,2) * timeRampFactor

           patch%gradientBuffer(:,1,patch%iGradientBuffer) = product(F, dim = 2)
           SAFE_DEALLOCATE(F)

           if (patch%iGradientBuffer == size(patch%gradientBuffer, 3)) then
              call patch%saveGradient()
              patch%iGradientBuffer = 0
           end if

        end select
     end do
  end do

end subroutine updateThermalActuatorGradient

function isThermalActuatorPatchValid(this, patchDescriptor, gridSize,                        &
     normalDirection, extent, simulationFlags, message) result(isPatchValid)

  ! <<< Derived types >>>
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags
  use ThermalActuator_mod, only : t_ThermalActuator

  implicit none

  ! <<< Arguments >>>
  class(t_ThermalActuator) :: this
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

end function isThermalActuatorPatchValid

subroutine hookThermalActuatorBeforeTimemarch(this, region, mode, referenceTimestep)

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
  integer, intent(in), optional :: referenceTimestep

  ! <<< Local variables >>>
  integer :: i, stat, fileUnit, mpiFileHandle, procRank, ierror
  integer :: referenceTimestep_
  integer(kind = MPI_OFFSET_KIND) :: referenceOffset
  class(t_Patch), pointer :: patch => null()
  logical :: gradientFileExists, controlForcingFileExists
  character(len = STRING_LENGTH) :: message

  referenceTimestep_ = -1
  if (PRESENT(referenceTimestep)) referenceTimestep_ = referenceTimestep

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
           if (procRank == 0) then
              inquire(file = trim(patch%controlForcingFilename), exist = controlForcingFileExists)
           end if
           call MPI_Bcast(controlForcingFileExists, 1, MPI_LOGICAL, 0, patch%comm, ierror)
           if (.not.controlForcingFileExists) then
              write(message, '(3A,I0.0,A)') "No control forcing information available for patch '", &
                   trim(patch%name), "' on grid ", patch%gridIndex, "!"
              call gracefulExit(patch%comm, message)
           end if
           if (controlForcingFileExists) then
              patch%iControlForcingBuffer = size(patch%controlForcingBuffer, 3) + 1
              call MPI_File_open(patch%comm, trim(patch%controlForcingFilename) // char(0),           &
                   MPI_MODE_WRONLY, MPI_INFO_NULL, mpiFileHandle, ierror)
              call MPI_File_get_size(mpiFileHandle, patch%controlForcingFileOffset, ierror)
              call MPI_File_close(mpiFileHandle, ierror)

              patch%controlForcingFileSize = patch%controlForcingFileOffset
              if (PRESENT(referenceTimestep)) then
                patch%forwardReferenceTimestep = referenceTimestep
                assert(patch%forwardReferenceTimestep.ge.0)
                referenceOffset = SIZEOF_SCALAR * product(int(patch%globalSize, MPI_OFFSET_KIND)) *            &
                                     size(patch%controlForcingBuffer, 2) * (4*referenceTimestep)
                patch%controlForcingFileOffset = patch%controlForcingFileOffset - referenceOffset
              end if
           end if

        case (ADJOINT)
          if (referenceTimestep_>0) then
            if (procRank == 0) then
               inquire(file = trim(patch%gradientFilename), exist = gradientFileExists)
            end if
            call MPI_Bcast(gradientFileExists, 1, MPI_LOGICAL, 0, patch%comm, ierror)
            if (.not.gradientFileExists) then
               write(message, '(3A,I0.0,A)') "No gradient information available for patch '", &
                    trim(patch%name), "' on grid ", patch%gridIndex, " for continued adjoint run!"
               call gracefulExit(patch%comm, message)
            end if
          else
            if (procRank == 0) then
              open(unit = getFreeUnit(fileUnit), file = trim(patch%gradientFilename),        &
                   iostat = stat, status = 'old')
              if (stat == 0) close(fileUnit, status = 'delete')
              open(unit = getFreeUnit(fileUnit), file = trim(patch%gradientFilename),        &
                   action = 'write', status = 'unknown')
              close(fileUnit)
            end if
          end if
           patch%iGradientBuffer = 0
           patch%gradientFileOffset = int(0, MPI_OFFSET_KIND)
           if (PRESENT(referenceTimestep)) then
             patch%adjointReferenceTimestep = referenceTimestep
             assert(patch%adjointReferenceTimestep.ge.0)
             referenceOffset = SIZEOF_SCALAR * product(int(patch%globalSize, MPI_OFFSET_KIND)) *            &
                                  size(patch%gradientBuffer, 2) * (4*referenceTimestep)
             patch%gradientFileOffset = patch%gradientFileOffset + referenceOffset
           end if
        end select

     end select
  end do

end subroutine hookThermalActuatorBeforeTimemarch

subroutine hookThermalActuatorAfterTimemarch(this, region, mode)

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

end subroutine hookThermalActuatorAfterTimemarch

PURE_FUNCTION thermalActuatorRampFunction(t, sigma, xi) result(timeRampFactor)

  implicit none

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(in) :: t, sigma, xi

  ! <<< Result >>>
  real(SCALAR_KIND) :: timeRampFactor

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  if (t <= -1.0_wp .or. t >= 1.0_wp) then
     timeRampFactor = 0.0_wp
  else
     timeRampFactor = 0.5_wp * (tanh(sigma * (t + xi)) - tanh(sigma * (t - xi)))
  end if

end function thermalActuatorRampFunction

subroutine collectThermalActuatorNorm(this, region, timeIntegrationNorm)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ThermalActuator_mod, only : t_ThermalActuator
  use ActuatorPatch_mod, only : t_ActuatorPatch

  implicit none

  ! <<< Arguments >>>
  class(t_ThermalActuator) :: this
  class(t_Region), intent(in) :: region
  real(SCALAR_KIND), intent(in) :: timeIntegrationNorm

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions
  class(t_Patch), pointer :: patch => null()
  SCALAR_TYPE, allocatable :: F(:)

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

           allocate(F(patch%nPatchPoints))
           call patch%collect(region%grids(j)%norm(:,1), F)

           patch%gradientBuffer(:,1,patch%iGradientBuffer) = F * timeIntegrationNorm
           SAFE_DEALLOCATE(F)

           if (patch%iGradientBuffer == size(patch%gradientBuffer, 3)) then
              call patch%saveGradient()
              patch%iGradientBuffer = 0
           end if

        end select
     end do
  end do

end subroutine collectThermalActuatorNorm
