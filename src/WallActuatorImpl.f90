#include "config.h"

subroutine setupWallActuator(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  Use Grid_mod, only :t_Grid
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use WallActuator_mod, only : t_WallActuator

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_WallActuator) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, gradientBufferSize
  class(t_Patch), pointer :: patch => null()
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)

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

        end select
     end do
  end if

     !region%states(j)%actuationAmount
     this%numP=1
     SAFE_DEALLOCATE(this%p)
     allocate(this%p(this%numP))
     this%p(1)=3._wp*pi/8._wp

end subroutine setupWallActuator

subroutine cleanupWallActuator(this)

  ! <<< Derived types >>>
  use WallActuator_mod, only : t_WallActuator

  implicit none

  ! <<< Arguments >>>
  class(t_WallActuator) :: this

  call this%cleanupBase()

end subroutine cleanupWallActuator

function computeWallActuatorSensitivity(this, region) result(instantaneousSensitivity)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use WallActuator_mod, only : t_WallActuator
  use CNSHelper
  use InputHelper, only : getOption
  use WavywallHelperImpl !computes specific quantities for Wall
  implicit none

  ! <<< Arguments >>>
  class(t_WallActuator) :: this
  class(t_Region), intent(in) :: region 


  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousSensitivity

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i,v,p,ii,jj,nDimensions,ierror
  SCALAR_TYPE, allocatable :: Jac(:)
  SCALAR_TYPE, allocatable :: F(:,:)
  SCALAR_TYPE, allocatable :: dQdxi(:,:,:),viscousFluxes(:,:,:),inviscidFluxes(:,:,:)
  SCALAR_TYPE, allocatable :: dJdp(:,:)
  SCALAR_TYPE, allocatable :: dMijdp(:,:,:,:)

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

     allocate(dMijdp(region%grids(i)%nGridPoints,nDimensions,nDimensions,this%numP))
     call compute_dMijdp(this,region%grids(i),dMijdp)
  
     allocate(F(region%grids(i)%nGridPoints,this%numP))
     
     allocate(dJdp(region%grids(i)%nGridPoints,this%numP))
     call compute_dJdp(this,region%grids(i),dJdp)
     
     allocate(Jac(region%grids(i)%nGridPoints))
     call compute_Jacobian(this,region%grids(i),Jac)    
 
     allocate(dQdxi(region%grids(i)%nGridPoints,nDimensions+2,nDimensions))
     allocate(inviscidFluxes(region%grids(i)%nGridPoints,nDimensions+2,nDimensions))
     allocate(viscousFluxes(region%grids(i)%nGridPoints,nDimensions+2,nDimensions))

     do ii = 1, nDimensions
          call region%grids(i)%firstDerivative(ii)%apply(dQdxi(:,:,ii),&
               region%grids(i)%localSize)
     end do

     call computeCartesianInvsicidFluxes(nDimensions,&
     region%states(i)%conservedVariables,&
     region%states(i)%velocity, region%states(i)%pressure(:,1),&
     inviscidFluxes)

     if (getOption("include_viscous_terms",.false.).and.&
          getOption("repeat_first_derivative", .true.)) then
          call computeCartesianViscousFluxes(nDimensions,&
          region%states(i)%velocity,&
          region%states(i)%stressTensor, region%states(i)%heatFlux,&
          viscousFluxes)
     end if

     F(:,:)=0._wp
     do p=1,this%numP
     do v=1,nDimensions+1
          F(:,p)=F(:,p)+region%states(i)%rightHandSide(:,v)*dJdp(:,p)*&
               region%states(i)%adjointVariables(:,v)

     do ii=1,nDimensions
     do jj=1,nDimensions
          F(:,p)=F(:,p)-Jac(:)*(inviscidFluxes(:,v,jj)-viscousFluxes(:,v,jj))*&
               dMijdp(:,ii,jj,p)*dQdxi(:,v,ii)
     end do !ii
     end do !jj

     !this is the surface contribution whose weight will be set the the
     !controller mollifier 
     do jj=1,nDimensions
          F(:,p)=F(:,p)-region%states(i)%adjointVariables(:,v)*Jac(:)*&
               (inviscidFluxes(:,v,jj)-viscousFluxes(:,v,jj))*dMijdp(:,2,jj,p)*&
               region%grids(i)%controlMollifier(:,1)
     end do

     end do !v
     end do !p

    instantaneousSensitivity = instantaneousSensitivity +&
          region%grids(i)%computeInnerProduct(F,F)

     SAFE_DEALLOCATE(viscousFluxes)
     SAFE_DEALLOCATE(inviscidFluxes)
     SAFE_DEALLOCATE(dMijdp)
     SAFE_DEALLOCATE(dQdxi)
     SAFE_DEALLOCATE(Jac)
     SAFE_DEALLOCATE(dJdp)
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

end function computeWallActuatorSensitivity

subroutine updateWallActuatorForcing(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use WallActuator_mod, only : t_WallActuator

  implicit none

  ! <<< Arguments >>>
  class(t_WallActuator) :: this
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

           !if (patch%iGradientBuffer == size(patch%gradientBuffer, 3))                       &
           !     call patch%loadGradient()

           patch%controlForcing(:,1:nDimensions+2) = 0.0_wp
           
           if (patch%iGradientBuffer == 1)                                                   &
                patch%iGradientBuffer = size(patch%gradientBuffer, 3) + 1

        end select
     end do
  end do

end subroutine updateWallActuatorForcing

subroutine updateWallActuatorGradient(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use WallActuator_mod, only : t_WallActuator

  use CNSHelper
 
  use InputHelper, only : getOption 
  implicit none

  ! <<< Arguments >>>
  class(t_WallActuator) :: this
  class(t_Region), intent(in) :: region

  ! <<< Local variables >>>
  logical::viscousOn
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, d, v, nDimensions
  class(t_Patch), pointer :: patch => null()

  nDimensions = size(region%globalGridSizes, 1)
  assert_key(nDimensions, (2))

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
           
           if (patch%iGradientBuffer == size(patch%gradientBuffer, 3)) then
              call patch%saveGradient()
              patch%iGradientBuffer = 0
           end if

        end select
     end do
  end do

end subroutine updateWallActuatorGradient

function isWallActuatorPatchValid(this, patchDescriptor, gridSize,                        &
     normalDirection, extent, simulationFlags, message) result(isPatchValid)

  ! <<< Derived types >>>
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags
  use WallActuator_mod, only : t_WallActuator

  implicit none

  ! <<< Arguments >>>
  class(t_WallActuator) :: this
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

end function isWallActuatorPatchValid

subroutine hookWallActuatorBeforeTimemarch(this, region, mode)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use WallActuator_mod, only : t_WallActuator
  use ActuatorPatch_mod, only : t_ActuatorPatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use InputHelper, only : getFreeUnit
  use ErrorHandler, only : gracefulExit
  use WavywallHelperImpl
  use Grid_enum
  use InputHelper, only : getOption, getRequiredOption  

  implicit none

  ! <<< Arguments >>>
  class(t_WallActuator) :: this
  class(t_Region) :: region
  integer, intent(in) :: mode

  ! <<< Local variables >>>
  integer :: i, g, fileUnit, mpiFileHandle, procRank, ierror
  class(t_Patch), pointer :: patch => null()
  logical :: fileExists
  logical:: hasNegativeJacobian
  character(len = STRING_LENGTH) :: message
  character(len = STRING_LENGTH) :: filename,outputPrefix

  select case (mode)
     case (FORWARD)
     do g = 1, size(region%grids)
        call updateWallCoordinates(this,region%grids(g))
        call region%grids(g)%update()
     end do
     call MPI_Barrier(region%comm, ierror)
     call getRequiredOption("grid_file", filename)
     filename=trim(filename)//'.p1'
     region%outputOn = .true.
     call region%saveData(QOI_GRID, filename)
     call MPI_Barrier(region%comm, ierror)
     outputPrefix = getOption("output_prefix", PROJECT_NAME)
     ! Save the Jacobian and normalized metrics.
     write(filename, '(2A)') trim(outputPrefix), ".Jacobian.p1.f"
     call region%saveData(QOI_JACOBIAN, filename)
     write(filename, '(2A)') trim(outputPrefix), ".metrics.p1.f"
     call region%saveData(QOI_METRICS, filename)
     region%outputOn = .false.
  end select

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
                   action = 'write', status = 'unknown')
              close(fileUnit)
           end if
           call MPI_Barrier(patch%comm, ierror)
           patch%iGradientBuffer = 0
           patch%gradientFileOffset = int(0, MPI_OFFSET_KIND)

        end select

     end select
  end do

end subroutine hookWallActuatorBeforeTimemarch

subroutine hookWallActuatorAfterTimemarch(this, region, mode)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use WallActuator_mod, only : t_WallActuator
  use ActuatorPatch_mod, only : t_ActuatorPatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  ! <<< Internal modules >>>
  use InputHelper, only : getFreeUnit
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_WallActuator) :: this
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

end subroutine hookWallActuatorAfterTimemarch

