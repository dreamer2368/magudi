#include "config.h"

subroutine setupReverseMigrator(this, region, outputPrefix, algorithm,                       &
     startTimestep, endTimestep, saveInterval, numIntermediateStates)

  ! <<< Derived types >>>
  use Region_type, only : t_Region
  use ReverseMigrator_type, only : t_ReverseMigrator, UNIFORM_CHECKPOINTING

  ! <<< Public members >>>
  use ReverseMigrator_mod, only : cleanupReverseMigrator

  implicit none

  ! <<< Arguments >>>
  type(t_ReverseMigrator) :: this
  type(t_Region), intent(in) :: region
  character(len = *), intent(in) :: outputPrefix, algorithm
  integer, intent(in) :: startTimestep, endTimestep, saveInterval, numIntermediateStates

  ! <<< Local variables >>>
  integer :: i

  assert(len_trim(outputPrefix) > 0)
  assert(startTimestep >= 0)
  assert(endTimestep >= startTimestep)
  assert(saveInterval > 0)
  assert(mod(endTimestep - startTimestep, saveInterval) == 0)
  assert(numIntermediateStates >= 0)
  assert(mod(numIntermediateStates, this%nStages) == 0)

  assert_key(algorithm, ( \
  'uniform checkpointing' \
  ))

  assert(allocated(region%states))
  assert(size(region%states) > 0)

  call cleanupReverseMigrator(this)

  this%outputPrefix = outputPrefix

  if (trim(algorithm) == "uniform checkpointing") then
     this%algorithm = UNIFORM_CHECKPOINTING
  end if

  this%startTimestep = startTimestep
  this%endTimestep = endTimestep
  this%saveInterval = saveInterval
  this%numIntermediateStates = numIntermediateStates

  allocate(this%temp_(size(region%states)))
  do i = 1, size(this%temp_)
     assert(region%grids(i)%nGridPoints > 0)
     assert(region%states(i)%nUnknowns > 0)
     allocate(this%temp_(i)%buffer(region%grids(i)%nGridPoints, &
          region%states(i)%nUnknowns, numIntermediateStates))
  end do

end subroutine setupReverseMigrator

subroutine cleanupReverseMigrator(this)

  ! <<< Derived types >>>
  use ReverseMigrator_type, only : t_ReverseMigrator

  implicit none

  ! <<< Arguments >>>
  type(t_ReverseMigrator) :: this

  ! <<< Local variables >>>
  integer :: i

  if (allocated(this%temp_)) then
     do i = 1, size(this%temp_)
        SAFE_DEALLOCATE(this%temp_(i)%buffer)
     end do
  end if

  SAFE_DEALLOCATE(this%temp_)

  this%loadedTimestep = -1

end subroutine cleanupReverseMigrator

subroutine migrateToSubstep(this, region, integrator, timestep, stage)

  ! <<< Derived types >>>
  use State_type, only : QOI_FORWARD_STATE
  use Region_type, only : t_Region
  use RK4Integrator_type, only : t_RK4Integrator
  use ReverseMigrator_type

  ! <<< Internal modules >>>
  use Region_mod, only : loadRegionData
  use RK4Integrator_mod, only : substepForward

  implicit none

  ! <<< Arguments >>>
  type(t_ReverseMigrator) :: this
  type(t_Region) :: region
  type(t_RK4Integrator) :: integrator
  integer, intent(in) :: timestep, stage

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, timestep_, stage_
  character(len = STRING_LENGTH) :: filename
  real(wp) :: time

  assert(timestep >= this%startTimestep .and. timestep <= this%endTimestep)
  assert(stage >= 1 .and. stage <= this%nStages)

  select case (this%algorithm)

  case (UNIFORM_CHECKPOINTING)

     if (this%loadedTimestep == -1 .or.                                                       &
          timestep < this%loadedTimestep .or.                                                &
          timestep > this%loadedTimestep + this%saveInterval .or.                            &
          (timestep == this%loadedTimestep .and.                                             &
          stage < this%nStages) .or.                                                         &
          (timestep == this%loadedTimestep + this%saveInterval .and.                         &
          stage == this%nStages)) then

        if (mod(timestep, this%saveInterval) == 0 .and. stage == this%nStages) then
           timestep_ = timestep
        else if (mod(timestep, this%saveInterval) == 0) then
           timestep_ = timestep - this%saveInterval
        else
           timestep_ = timestep - mod(timestep, this%saveInterval)
        end if

        write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", timestep_, ".q"
        call loadRegionData(region, QOI_FORWARD_STATE, filename)
        time = real(region%states(1)%plot3dAuxiliaryData(4), wp)
        this%loadedTimestep = timestep_

        i = 1
        do j = 1, size(region%states)
           this%temp_(j)%buffer(:,:,i) = region%states(j)%conservedVariables
        end do

        if (this%loadedTimestep /= this%endTimestep) then
           do timestep_ = this%loadedTimestep + 1, this%loadedTimestep + this%saveInterval
              do stage_ = 1, this%nStages
                 if (timestep_ == this%loadedTimestep + this%saveInterval .and.              &
                      stage_ == this%nStages) exit
                 call substepForward(integrator, region, time, timestep_, stage_)
                 i = i + 1
                 do j = 1, size(region%states)
                    this%temp_(j)%buffer(:,:,i) = region%states(j)%conservedVariables
                 end do
              end do
           end do
        end if

     end if

     i = (timestep - 1 - this%loadedTimestep) * this%nStages + stage + 1
     assert(i >= 1 .and. i <= this%numIntermediateStates)
     do j = 1, size(region%states)
        region%states(j)%conservedVariables = this%temp_(j)%buffer(:,:,i)
     end do

  end select

end subroutine migrateToSubstep
