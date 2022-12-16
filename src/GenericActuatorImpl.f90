#include "config.h"

subroutine setupGenericActuator(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use GenericActuator_mod, only : t_GenericActuator

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_GenericActuator) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions, nActuatorComponents
  class(t_Patch), pointer :: patch => null()

  call this%cleanup()

  if (.not. region%simulationFlags%enableController) return

  call this%setupBase(region%simulationFlags, region%solverOptions)

  nDimensions = size(region%globalGridSizes, 1)
  nActuatorComponents = nDimensions + 2
  assert_key(nDimensions, (1, 2, 3))

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
             allocate(patch%gradientBuffer(patch%nPatchPoints, nActuatorComponents,           &
                                           this%controllerBufferSize))
             patch%gradientBuffer = 0.0_wp
           end if
           SAFE_DEALLOCATE(patch%controlForcingBuffer)
           allocate(patch%controlForcingBuffer(patch%nPatchPoints, nActuatorComponents,       &
                                               this%controllerBufferSize))
           patch%controlForcingBuffer = 0.0_wp

        end select
     end do
  end if

end subroutine setupGenericActuator

subroutine cleanupGenericActuator(this)

  ! <<< Derived types >>>
  use GenericActuator_mod, only : t_GenericActuator

  implicit none

  ! <<< Arguments >>>
  class(t_GenericActuator) :: this

  call this%cleanupBase()

end subroutine cleanupGenericActuator

function computeGenericActuatorSensitivity(this, region) result(instantaneousSensitivity)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use GenericActuator_mod, only : t_GenericActuator

  ! <<< Internal modules >>>
  use Patch_factory, only : computeQuadratureOnPatches

  implicit none

  ! <<< Arguments >>>
  class(t_GenericActuator) :: this
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousSensitivity

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions, nActuatorComponents, ierror
  SCALAR_TYPE, allocatable :: F(:,:)

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(size(region%grids) == size(region%states))

  instantaneousSensitivity = 0.0_wp

  do i = 1, size(region%grids)

     nDimensions = region%grids(i)%nDimensions
     nActuatorComponents = nDimensions + 2
     assert_key(nDimensions, (1, 2, 3))

     assert(region%grids(i)%nGridPoints > 0)
     assert(allocated(region%grids(i)%controlMollifier))
     assert(size(region%grids(i)%controlMollifier, 1) == region%grids(i)%nGridPoints)
     assert(size(region%grids(i)%controlMollifier, 2) == 1)
     assert(allocated(region%states(i)%adjointVariables))
     assert(size(region%states(i)%adjointVariables, 1) == region%grids(i)%nGridPoints)
     assert(size(region%states(i)%adjointVariables, 2) >= nDimensions + 2)

     allocate(F(region%grids(i)%nGridPoints, nActuatorComponents))

     do j = 1, nActuatorComponents
       F(:,j) = region%states(i)%adjointVariables(:,j) *                                     &
                region%grids(i)%controlMollifier(:,1)
     end do

     instantaneousSensitivity = instantaneousSensitivity +                                   &
                          computeQuadratureOnPatches(region%patchFactories,'ACTUATOR',       &
                                                      region%grids(i), sum(F**2,2))

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

end function computeGenericActuatorSensitivity

subroutine updateGenericActuatorForcing(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use GenericActuator_mod, only : t_GenericActuator

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT, LINEARIZED

  implicit none

  ! <<< Arguments >>>
  class(t_GenericActuator) :: this
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

           patch%iControlForcingBuffer = patch%iControlForcingBuffer - 1

           assert(patch%iControlForcingBuffer >= 1)
           assert(patch%iControlForcingBuffer <= size(patch%controlForcingBuffer, 3))

           if (patch%iControlForcingBuffer == size(patch%controlForcingBuffer, 3))                       &
                call patch%loadForcing()

           patch%controlForcing =                                                         &
                   patch%controlForcingBuffer(:,:,patch%iControlForcingBuffer)

           if (patch%iControlForcingBuffer == 1)                                                   &
                patch%iControlForcingBuffer = size(patch%controlForcingBuffer, 3) + 1

        end select
     end do
  end do

end subroutine updateGenericActuatorForcing

subroutine updateGenericActuatorDeltaForcing(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use GenericActuator_mod, only : t_GenericActuator

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT, LINEARIZED

  implicit none

  ! <<< Arguments >>>
  class(t_GenericActuator) :: this
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
               call patch%loadDeltaForcing()

          patch%deltaControlForcing =                                                       &
                                 patch%gradientBuffer(:,:,patch%iGradientBuffer)

          if (patch%iGradientBuffer == 1)                                                   &
              patch%iGradientBuffer = size(patch%gradientBuffer, 3) + 1

        end select
     end do
  end do

end subroutine updateGenericActuatorDeltaForcing

subroutine migrateToGenericActuatorForcing(this, region, startTimeStep, endTimeStep, nStages, iTimeStep, jSubStep)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use GenericActuator_mod, only : t_GenericActuator

  implicit none

  ! <<< Arguments >>>
  class(t_GenericActuator) :: this
  class(t_Region), intent(in) :: region
  integer, intent(in) :: startTimeStep, endTimeStep, nStages, iTimeStep, jSubStep

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, nDimensions
  integer :: bufferIndex, bufferOffsetIndex, iControlForcingBuffer
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

          bufferIndex = (endTimeStep - iTimeStep)*nStages + (nStages-jSubStep)
          if (patch%adjointReferenceTimestep>0) bufferIndex = bufferIndex                     &
                                        + patch%adjointReferenceTimestep*nStages
          bufferOffsetIndex = bufferIndex/this%controllerBufferSize
          iControlForcingBuffer = MODULO( bufferIndex,this%controllerBufferSize ) + 1

          patch%iControlForcingBuffer = iControlForcingBuffer

           assert(patch%iControlForcingBuffer >= 1)
           assert(patch%iControlForcingBuffer <= size(patch%controlForcingBuffer, 3))

           if (patch%bufferOffsetIndex .ne. bufferOffsetIndex)                               &
                call patch%pinpointForcing(bufferOffsetIndex)

           patch%controlForcing =                                                             &
                     patch%controlForcingBuffer(:,:,patch%iControlForcingBuffer)

        end select
     end do
  end do

end subroutine migrateToGenericActuatorForcing

subroutine updateGenericActuatorGradient(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use GenericActuator_mod, only : t_GenericActuator

  implicit none

  ! <<< Arguments >>>
  class(t_GenericActuator) :: this
  class(t_Region), intent(in) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, nActuatorComponents
  class(t_Patch), pointer :: patch => null()
  SCALAR_TYPE, allocatable :: F(:,:)

  if (.not. allocated(region%patchFactories)) return

  nDimensions = size(region%globalGridSizes, 1)
  nActuatorComponents = nDimensions + 2
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

           call patch%collect(region%grids(j)%controlMollifier(:,1), F(:,1))
           do k = 1, nActuatorComponents
              call patch%collect(region%states(j)%adjointVariables(:,k), F(:,2))
              patch%gradientBuffer(:,k,patch%iGradientBuffer) = product(F, dim = 2)
           end do

           SAFE_DEALLOCATE(F)

           if (patch%iGradientBuffer == size(patch%gradientBuffer, 3)) then
              call patch%saveGradient()
              patch%iGradientBuffer = 0
           end if

        end select
     end do
  end do

end subroutine updateGenericActuatorGradient

function isGenericActuatorPatchValid(this, patchDescriptor, gridSize,                        &
     normalDirection, extent, simulationFlags, message) result(isPatchValid)

  ! <<< Derived types >>>
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags
  use GenericActuator_mod, only : t_GenericActuator

  implicit none

  ! <<< Arguments >>>
  class(t_GenericActuator) :: this
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

end function isGenericActuatorPatchValid

subroutine hookGenericActuatorBeforeTimemarch(this, region, mode, referenceTimestep)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use Controller_mod, only : t_Controller
  use ActuatorPatch_mod, only : t_ActuatorPatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT, LINEARIZED

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
  logical :: fileExists
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
           if (procRank == 0) inquire(file = trim(patch%controlForcingFilename), exist = fileExists)
           call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, patch%comm, ierror)
           if (.not. fileExists) then
              write(message, '(3A,I0.0,A)') "No control forcing information available for patch '", &
                   trim(patch%name), "' on grid ", patch%gridIndex, "!"
              call gracefulExit(patch%comm, message)
           end if
           patch%iControlForcingBuffer = size(patch%controlForcingBuffer, 3) + 1
           call MPI_File_open(patch%comm, trim(patch%controlForcingFilename) // char(0),           &
                MPI_MODE_WRONLY, MPI_INFO_NULL, mpiFileHandle, ierror)
           call MPI_File_get_size(mpiFileHandle, patch%controlForcingFileOffset, ierror)
           call MPI_File_close(mpiFileHandle, ierror)

           patch%controlForcingFileSize = patch%controlForcingFileOffset
           if (PRESENT(referenceTimestep)) then
             patch%forwardReferenceTimestep = referenceTimestep
             assert(patch%forwardReferenceTimestep>0)
             referenceOffset = SIZEOF_SCALAR * product(int(patch%globalSize, MPI_OFFSET_KIND)) *            &
                                  size(patch%controlForcingBuffer, 2) * (4*referenceTimestep)
             patch%controlForcingFileOffset = patch%controlForcingFileOffset - referenceOffset
           end if
        case (ADJOINT)
          if (referenceTimestep_>0) then
            if (procRank == 0) then
               inquire(file = trim(patch%gradientFilename), exist = fileExists)
            end if
            call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, patch%comm, ierror)
            if (.not.fileExists) then
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
           call MPI_Barrier(patch%comm, ierror)
           patch%iGradientBuffer = 0
           patch%gradientFileOffset = int(0, MPI_OFFSET_KIND)

           if (PRESENT(referenceTimestep)) then
             patch%adjointReferenceTimestep = referenceTimestep
             assert(patch%adjointReferenceTimestep.ge.0)
             referenceOffset = SIZEOF_SCALAR * product(int(patch%globalSize, MPI_OFFSET_KIND)) *            &
                                  size(patch%gradientBuffer, 2) * (4*referenceTimestep)
             patch%gradientFileOffset = patch%gradientFileOffset + referenceOffset
           end if

         case (LINEARIZED) ! use gradient, adjoint variables for linearized run
            if (procRank == 0) then
               inquire(file = trim(patch%gradientFilename), exist = fileExists)
            end if
            call MPI_Bcast(fileExists, 1, MPI_LOGICAL, 0, patch%comm, ierror)
            if (.not.fileExists) then
               write(message, '(3A,I0.0,A)') "No gradient information available for patch '", &
                    trim(patch%name), "' on grid ", patch%gridIndex, "!"
               call gracefulExit(patch%comm, message)
            end if
            if (fileExists) then
               patch%iGradientBuffer = size(patch%gradientBuffer, 3) + 1
               call MPI_File_open(patch%comm, trim(patch%gradientFilename) // char(0),           &
                    MPI_MODE_WRONLY, MPI_INFO_NULL, mpiFileHandle, ierror)
               call MPI_File_get_size(mpiFileHandle, patch%gradientFileOffset, ierror)
               call MPI_File_close(mpiFileHandle, ierror)

               patch%gradientFileSize = patch%gradientFileOffset
               if (PRESENT(referenceTimestep)) then
                 patch%forwardReferenceTimestep = referenceTimestep
                 assert(patch%forwardReferenceTimestep.ge.0)
                 referenceOffset = SIZEOF_SCALAR * product(int(patch%globalSize, MPI_OFFSET_KIND)) *            &
                                      size(patch%gradientBuffer, 2) * (4*referenceTimestep)
                 patch%gradientFileOffset = patch%gradientFileOffset - referenceOffset
               end if
            end if

        end select

     end select
  end do

end subroutine hookGenericActuatorBeforeTimemarch

subroutine hookGenericActuatorAfterTimemarch(this, region, mode)

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

end subroutine hookGenericActuatorAfterTimemarch

subroutine collectGenericActuatorNorm(this, region, timeIntegrationNorm)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use GenericActuator_mod, only : t_GenericActuator

  implicit none

  ! <<< Arguments >>>
  class(t_GenericActuator) :: this
  class(t_Region), intent(in) :: region
  real(SCALAR_KIND), intent(in) :: timeIntegrationNorm

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, nActuatorComponents
  class(t_Patch), pointer :: patch => null()
  SCALAR_TYPE, allocatable :: F(:)

  if (.not. allocated(region%patchFactories)) return

  nDimensions = size(region%globalGridSizes, 1)
  nActuatorComponents = nDimensions + 2
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

           do k = 1, nActuatorComponents
             patch%gradientBuffer(:,k,patch%iGradientBuffer) = F * timeIntegrationNorm
           end do

           SAFE_DEALLOCATE(F)

           if (patch%iGradientBuffer == size(patch%gradientBuffer, 3)) then
              call patch%saveGradient()
              patch%iGradientBuffer = 0
           end if

        end select
     end do
  end do

end subroutine collectGenericActuatorNorm
