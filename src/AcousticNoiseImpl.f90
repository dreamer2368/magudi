#include "config.h"

subroutine setupAcousticNoise(this, region)

  ! <<< External modules >>>
  use MPI

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
  integer :: i, j, ierror
  character(len = STRING_LENGTH) :: filename
  SCALAR_TYPE :: endTime
  integer, parameter :: wp = SCALAR_KIND

  assert(allocated(region%states))
  assert(size(region%states) > 0)

  call this%cleanup()

  call this%setupBase(region%simulationFlags, region%solverOptions)

  allocate(this%data_(size(region%globalGridSizes, 2)))

  do i = 1, size(this%data_)
     do j = 1, size(region%grids)
        if (region%grids(j)%index == i) then
           allocate(this%data_(i)%meanPressure(region%grids(j)%nGridPoints, 1))
           region%states(j)%dummyFunction => this%data_(i)%meanPressure
           exit
        end if
     end do
     call MPI_Barrier(region%comm, ierror)
  end do

  if (.not. region%simulationFlags%useTargetState) then
     call getRequiredOption("mean_pressure_file", filename, region%comm)
     call region%loadData(QOI_DUMMY_FUNCTION, filename)
  else
     filename = getOption("mean_pressure_file", "")
     if (len_trim(filename) > 0) then
        call region%loadData(QOI_DUMMY_FUNCTION, filename)
     else
        do i = 1, size(region%states)
           assert(allocated(region%states(i)%targetState))
           j = region%grids(i)%index
           call computeDependentVariables(region%grids(i)%nDimensions,                       &
                region%states(i)%targetState, region%solverOptions%ratioOfSpecificHeats,     &
                pressure = this%data_(j)%meanPressure(:,1))
        end do
     end if
  end if

do i = 1, size(region%states)
j = region%grids(i)%index
!load necessary time dependent acoustic source stuff here (the \Omega_t(t))
this%data_(j)%t1=getOption("cost_functional_timeHorizon/startTime",0._wp)
endTime=region%getTimeStepSize()*real(getOption("number_of_timesteps", 1000))
this%data_(j)%t2=getOption("cost_functional_timeHorizon/stopTime",endTime)
this%data_(j)%sigma=getOption("cost_functional_MollifierStrength",40._wp)
this%data_(j)%normalize=getOption("cost_functional_normalizeTimeMollifier",1._wp)
assert(getOption("use_constant_CFL_mode")==.false.)
end do

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

function computeAcousticNoise(this, region) result(instantaneousFunctional)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use AcousticNoise_mod, only : t_AcousticNoise

  ! <<< Internal modules >>>
  use InputHelper, only : getFreeUnit,getOption

  ! <<< Arguments >>>
  class(t_AcousticNoise) :: this
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, ierror
  SCALAR_TYPE, allocatable :: F(:,:)

  integer :: fileUnit, iostat, procRank
  character(len = STRING_LENGTH) :: filename, message,outputPrefix

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

     j = region%grids(i)%index

     assert(associated(this%data_(j)%meanPressure))
     assert(size(this%data_(j)%meanPressure, 1) == region%grids(i)%nGridPoints)
     assert(size(this%data_(j)%meanPressure, 2) == 1)

     allocate(F(region%grids(i)%nGridPoints, 1))
     F = region%states(i)%pressure - this%data_(j)%meanPressure

     instantaneousFunctional = instantaneousFunctional +                                     &
          region%grids(i)%computeInnerProduct(F, F, region%grids(i)%targetMollifier(:,1))*&
              0.5_wp*(tanh(this%data_(j)%sigma*(region%states(j)%time-this%data_(j)%t1))-&
              tanh(this%data_(j)%sigma*(region%states(j)%time-this%data_(j)%t2)))/&
              this%data_(j)%normalize

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


!  outputPrefix = getOption("output_prefix", PROJECT_NAME)
!  call MPI_Comm_rank(region%comm, procRank, ierror)
!
!  write(filename, '(2A)') trim(outputPrefix), ".computeAcousticNoise.txt"
!  if (procRank == 0) then
!        open(unit = getFreeUnit(fileUnit), file = trim(filename), action ='readwrite',      &
!             position = 'append', iostat = iostat)
!  end if
!
!  call MPI_Bcast(iostat, 1, MPI_INTEGER, 0, region%comm, ierror)
!  if (iostat /= 0) then
!     write(message, "(2A)") trim(filename), ": Failed to open file for writing!"
!     call gracefulExit(region%comm, message)
!  end if
!
!  if (procRank == 0)&
!          write(fileUnit, '(2(1X,SP,' // SCALAR_FORMAT // '))') region%states(1)%time,&
!               0.5_wp*(tanh(this%data_(1)%sigma*(region%states(1)%time-this%data_(1)%t1))-&
!               tanh(this%data_(1)%sigma*(region%states(1)%time-this%data_(1)%t2)))/&
!               this%data_(1)%normalize
!
!  if (procRank == 0) close(fileUnit)

end function computeAcousticNoise

subroutine computeAcousticNoiseAdjointForcing(this, simulationFlags, solverOptions,          &
     region, grid, state, patch)

  ! <<< Derived types >>>
  use Region_mod, only: t_Region
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use AcousticNoise_mod, only : t_AcousticNoise
  use SolverOptions_mod, only : t_SolverOptions
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use InputHelper, only : getFreeUnit,getOption
  use MPI

  implicit none

  ! <<< Arguments >>>
  class(t_AcousticNoise) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Region), intent(in) :: region
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  class(t_CostTargetPatch) :: patch

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, gridIndex, patchIndex
  SCALAR_TYPE, allocatable :: meanPressure(:)
  SCALAR_TYPE :: F

  integer :: fileUnit, iostat, procRank,ierror
  character(len = STRING_LENGTH) :: filename, message,outputPrefix


  !debugging
  integer::firstProcMatch

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  allocate(meanPressure(patch%nPatchPoints))
  i = grid%index
  call patch%collect(this%data_(i)%meanPressure(:,1), meanPressure)

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

           F = - 2.0_wp * grid%targetMollifier(gridIndex, 1) *                               &
                (solverOptions%ratioOfSpecificHeats - 1.0_wp) *                              &
                (state%pressure(gridIndex, 1) - meanPressure(patchIndex))*&
               0.5_wp*(tanh(this%data_(grid%index)%sigma*(state%time-this%data_(grid%index)%t1))&
               -tanh(this%data_(grid%index)%sigma*(state%time-this%data_(grid%index)%t2)))/&
               this%data_(grid%index)%normalize

           patch%adjointForcing(patchIndex,nDimensions+2) = F
           patch%adjointForcing(patchIndex,2:nDimensions+1) =                                &
                - state%velocity(gridIndex,:) * F
           patch%adjointForcing(patchIndex,1) =                                              &
                0.5_wp * sum(state%velocity(gridIndex,:) ** 2) * F

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
     end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

  SAFE_DEALLOCATE(meanPressure)


!  outputPrefix = getOption("output_prefix", PROJECT_NAME)
!  call MPI_Comm_rank(region%comm, procRank, ierror)
!
!  write(filename, '(2A)') trim(outputPrefix),".computeAcousticNoiseAdjointForcing.txt"
!  if (procRank == 8) then
!       open(unit = getFreeUnit(fileUnit), file = trim(filename), action='readwrite',      &
!             position = 'append', iostat = iostat)
!  end if
!
!  if (procRank == 8)&
!          write(fileUnit, '(2(1X,SP,' // SCALAR_FORMAT // '))')state%time,&
!               0.5_wp*(tanh(this%data_(grid%index)%sigma*(state%time-this%data_(grid%index)%t1))&
!               -tanh(this%data_(grid%index)%sigma*(state%time-this%data_(grid%index)%t2)))/&
!               this%data_(grid%index)%normalize
!
!  if (procRank == 8) close(fileUnit)


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
