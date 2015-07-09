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
  use WavywallHelperImpl,&
     only:compute_dMijdp,compute_dJacobiandp,updateWallCoordinates

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

     this%controlIndex=0
     this%numP=1
     SAFE_DEALLOCATE(this%p)
     allocate(this%p(this%numP))
     this%p(1)=0._wp !3._wp*pi/8._wp

     SAFE_DEALLOCATE(this%gradient)
     allocate(this%gradient(this%numP))
     this%gradient=0._wp

     SAFE_DEALLOCATE(this%instantaneousGradient)
     allocate(this%instantaneousGradient(this%numP))
     this%instantaneousGradient=0._wp

     this%sensitivity=0._wp
 
     SAFE_DEALLOCATE(this%dJacobiandp) 
     allocate(this%dJacobiandp(region%grids(1)%nGridPoints,this%numP))
     
     SAFE_DEALLOCATE(this%dMijdp)
     allocate(this%dMijdp(region%grids(1)%nGridPoints,&
          region%grids(1)%nDimensions,region%grids(1)%nDimensions,this%numP)) 

     call updateWallCoordinates(this,region%grids(1))
     call region%grids(1)%update()
     call compute_dJacobiandp(this,region%grids(1),this%dJacobiandp)
     call compute_dMijdp(this,region%grids(1),this%dMijdp)

     !you could output mesh here so that you're not handcuffed by whatever was
     !done to create the mesh

 
end subroutine setupWallActuator

subroutine cleanupWallActuator(this)

  ! <<< Derived types >>>
  use WallActuator_mod, only : t_WallActuator

  implicit none

  ! <<< Arguments >>>
  class(t_WallActuator) :: this

  call this%cleanupBase()

  SAFE_DEALLOCATE(this%dJacobiandp)
  SAFE_DEALLOCATE(this%dMijdp)

end subroutine cleanupWallActuator

subroutine computeWallActuatorSensitivity(this,timeIntegrator, region) 

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use WallActuator_mod, only : t_WallActuator
  use InputHelper, only : getOption
  use TimeIntegrator_mod, only : t_TimeIntegrator
  implicit none

  ! <<< Arguments >>>
  class(t_WallActuator) :: this
  class(t_Region), intent(in) :: region 
  class(t_TimeIntegrator),intent(in) :: timeIntegrator

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i,v,p,ii,jj,nDimensions,ierror


  !for the passive wall shape optimization 
  !gradient already has time integrated and volume integrated
  this%sensitivity=0._wp
  do i = 1, size(region%grids)
    assert(size(this%gradient,1)==this%numP)
    this%sensitivity = this%sensitivity +&
          region%grids(i)%computeInnerProduct(this%gradient,this%gradient)
  end do

  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE,this%sensitivity, 1,                         &
       SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(this%sensitivity, 1, SCALAR_TYPE_MPI,                            &
          0, region%grids(i)%comm, ierror)
  end do

end subroutine computeWallActuatorSensitivity

subroutine computeWallActuatorGradient(this,timeIntegrator,region)

! <<< External modules >>>
use MPI

! <<< Derived types >>>
use Region_mod, only : t_Region
use WallActuator_mod, only : t_WallActuator
use CNSHelper
use InputHelper, only : getOption
use TimeIntegrator_mod, only : t_TimeIntegrator
use Region_mod, only : t_Region
use TimeIntegrator_mod, only : t_TimeIntegrator
use WallActuator_mod, only : t_WallActuator

implicit none
class(t_WallActuator) :: this
class(t_Region), intent(in) :: region
class(t_TimeIntegrator),intent(in) :: timeIntegrator

! <<< Local variables >>>
integer, parameter :: wp = SCALAR_KIND
integer :: i,v,p,ii,jj,nDimensions,ierror
SCALAR_TYPE, allocatable :: F(:,:),G(:,:)
SCALAR_TYPE, allocatable ::dQdxi(:,:,:),viscousFluxes(:,:,:),inviscidFluxes(:,:,:)

assert(allocated(region%grids))
assert(allocated(region%states))
assert(size(region%grids) == size(region%states))
assert(size(region%grids)==1)


do i = 1, size(region%grids)

nDimensions = region%grids(i)%nDimensions
assert_key(nDimensions, (2))

assert(region%grids(i)%nGridPoints > 0)
assert(allocated(region%grids(i)%controlMollifier))
assert(size(region%grids(i)%controlMollifier, 1) ==region%grids(i)%nGridPoints)
assert(size(region%grids(i)%controlMollifier, 2) == 1)
assert(allocated(region%states(i)%adjointVariables))
assert(size(region%states(i)%adjointVariables, 1) ==region%grids(i)%nGridPoints)
assert(size(region%states(i)%adjointVariables, 2) >= nDimensions + 2)
assert(allocated(this%dMijdp))
assert(allocated(this%dJacobiandp))

allocate(G(region%grids(i)%nGridPoints,nDimensions+2))
allocate(F(region%grids(i)%nGridPoints,nDimensions+2))
allocate(dQdxi(region%grids(i)%nGridPoints,nDimensions+2,nDimensions))
allocate(inviscidFluxes(region%grids(i)%nGridPoints,nDimensions+2,nDimensions))
allocate(viscousFluxes(region%grids(i)%nGridPoints,nDimensions+2,nDimensions))

G=0._wp
F=0._wp
dQdxi=0._wp
inviscidFluxes=0._wp
viscousFluxes=0._wp

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

do p=1,this%numP

dQdxi=0.
F=0.

do ii=1,nDimensions
do jj=1,nDimensions
do v=1,nDimensions+2
dQdxi(:,v,ii)=dQdxi(:,v,ii)+&
          (inviscidFluxes(:,v,jj)-viscousFluxes(:,v,jj))*&
          this%dMijdp(:,ii,jj,p)
end do
end do
call region%grids(i)%firstDerivative(ii)%apply(dQdxi(:,:,ii),&
          region%grids(i)%localSize)
end do

F=sum(dQdxi,dim=3)
do v=1,nDimensions+2
F(:,v)=F(:,v)*region%grids(i)%jacobian(:,1)
F(:,v)=F(:,v)+region%states(i)%rightHandSide(:,v)*this%dJacobiandp(:,p)
end do !var

G(:,:)=region%states(i)%adjointVariables(:,:)

this%instantaneousGradient(p)=-region%grids(i)%computeInnerProduct(G,F)

this%gradient(p)=this%gradient(p)+&
     timeIntegrator%norm(timeIntegrator%stage)*region%getTimeStepSize()*&
     this%instantaneousGradient(p)

end do !wall parameters

SAFE_DEALLOCATE(viscousFluxes)
SAFE_DEALLOCATE(inviscidFluxes)
SAFE_DEALLOCATE(dQdxi)
SAFE_DEALLOCATE(F)
SAFE_DEALLOCATE(G)

end do !grids

end subroutine computeWallActuatorGradient


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


end subroutine updateWallActuatorForcing

subroutine updateWallActuatorGradient(this, region)

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use WallActuator_mod, only : t_WallActuator
  use TimeIntegrator_mod, only : t_TimeIntegrator
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
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()

  nDimensions = size(region%globalGridSizes, 1)
  assert_key(nDimensions, (2))
 
  ! Connect to the previously allocated time integrator.
  !call region%timeIntegratorFactory%connect(timeIntegrator)
  !assert(associated(timeIntegrator))

     !region%getTimeStepSize()
     !instantaneousGradient(:) = controller%computeSensitivity(region)
     !gradient(:) =gradient(:) +&
     !timeIntegrator%norm(i)*region%getTimeStepSize()*instantaneousCostSensitivity(:)

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
  use Grid_enum
  use InputHelper, only : getOption, getRequiredOption  
  use WavywallHelperImpl,&
     only:compute_dMijdp,compute_dJacobiandp,updateWallCoordinates

  implicit none

  ! <<< Arguments >>>
  class(t_WallActuator) :: this
  class(t_Region) :: region
  integer, intent(in) :: mode

  ! <<< Local variables >>>
  integer :: i, g, fileUnit, mpiFileHandle, procRank, ierror
  class(t_Patch), pointer :: patch => null()
  integer, parameter :: wp = SCALAR_KIND
  logical :: fileExists
  logical:: hasNegativeJacobian
  character(len = STRING_LENGTH) :: message
  character(len = STRING_LENGTH) :: griditeration
  character(len = STRING_LENGTH) :: filename,outputPrefix
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)

  select case (mode)

     !case (OPTIMIZE)
     !Read previously written g if one doesn't exist then
     !assume steepest descent 

     case (FORWARD)
     
     do g = 1, size(region%grids)
        this%p(:)=-region%states(g)%actuationAmount*this%gradient(:)
        call updateWallCoordinates(this,region%grids(g))
        call region%grids(g)%update()
        call compute_dJacobiandp(this,region%grids(g),this%dJacobiandp)
        call compute_dMijdp(this,region%grids(g),this%dMijdp)
     end do

     write (griditeration, "(I4)") this%controlIndex
     call MPI_Barrier(region%comm, ierror)
     call getRequiredOption("grid_file", filename)
     filename=trim(filename)//'.'//adjustl(trim(griditeration))
     region%outputOn = .true.
     call region%saveData(QOI_GRID, filename)
     call MPI_Barrier(region%comm, ierror)
     outputPrefix = getOption("output_prefix", PROJECT_NAME)
     ! Save the Jacobian and normalized metrics.
     write(filename, '(2A)') trim(outputPrefix),&
          ".Jacobian.f"//adjustl(trim(griditeration))
     call region%saveData(QOI_JACOBIAN,trim(filename))
     write(filename, '(2A)') trim(outputPrefix),&
          ".metrics.f"//adjustl(trim(griditeration))
     call region%saveData(QOI_METRICS,trim(filename))  
     region%outputOn = .false.
     this%controlIndex=this%controlIndex+1

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


