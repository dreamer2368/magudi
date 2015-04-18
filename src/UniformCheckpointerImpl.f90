#include "config.h"

subroutine setupUniformCheckpointer(this, region, timeIntegrator, outputPrefix,              &
     startTimestep, endTimestep, saveInterval, numIntermediateStates)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use TimeIntegrator_mod, only : t_TimeIntegrator
  use UniformCheckpointer_mod, only : t_UniformCheckpointer

  ! <<< Internal modules >>>
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_UniformCheckpointer) :: this
  class(t_Region), intent(in) :: region
  class(t_TimeIntegrator), intent(in) :: timeIntegrator
  character(len = *), intent(in) :: outputPrefix
  integer, intent(in) :: startTimestep, endTimestep, saveInterval, numIntermediateStates

  ! <<< Local variables >>>
  integer :: i, numIntermediateStates_
  character(len = STRING_LENGTH) :: message

  call this%cleanup()

  numIntermediateStates_ = min(numIntermediateStates, saveInterval * timeIntegrator%nStages)

  if (numIntermediateStates_ /= saveInterval * timeIntegrator%nStages) then
     write(message, '(2(A,I0.0),A)') "Uniform checkpointing requires storage for ",          &
          saveInterval * timeIntegrator%nStages, " intermediate states, but only ",          &
          numIntermediateStates_, " are available!"
     call gracefulExit(region%comm, message)
  end if

  call this%setupBase(region, timeIntegrator, outputPrefix, startTimestep, endTimestep,      &
       saveInterval, numIntermediateStates_)

  allocate(this%data_(size(region%states)))
  do i = 1, size(this%data_)
     assert(region%grids(i)%nGridPoints > 0)
     assert(region%solverOptions%nUnknowns > 0)
     allocate(this%data_(i)%buffer(region%grids(i)%nGridPoints,                              &
          region%solverOptions%nUnknowns, numIntermediateStates_))
  end do

end subroutine setupUniformCheckpointer

subroutine cleanupUniformCheckpointer(this)

  ! <<< Derived types >>>
  use UniformCheckpointer_mod, only : t_UniformCheckpointer

  implicit none

  ! <<< Arguments >>>
  class(t_UniformCheckpointer) :: this

  ! <<< Local variables >>>
  integer :: i

  call this%cleanupBase()

  if (allocated(this%data_)) then
     do i = 1, size(this%data_)
        SAFE_DEALLOCATE(this%data_(i)%buffer)
     end do
  end if

  SAFE_DEALLOCATE(this%data_)

end subroutine cleanupUniformCheckpointer

subroutine uniformCheckpointingMigrateTo(this, region, timeIntegrator, timestep, stage)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use TimeIntegrator_mod, only : t_TimeIntegrator
  use UniformCheckpointer_mod, only : t_UniformCheckpointer

  ! <<< Enumerations >>>
  use State_enum, only : QOI_FORWARD_STATE

  implicit none

  ! <<< Arguments >>>
  class(t_UniformCheckpointer) :: this
  class(t_Region) :: region
  class(t_TimeIntegrator) :: timeIntegrator
  integer, intent(in) :: timestep, stage

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, timestep_, stage_
  character(len = STRING_LENGTH) :: filename
  real(wp) :: time, timeStepSize

  assert(timestep >= this%startTimestep .and. timestep <= this%endTimestep)
  assert(stage >= 1 .and. stage <= timeIntegrator%nStages)

  if (this%loadedTimestep == -1 .or.                                                         &
       timestep < this%loadedTimestep .or.                                                   &
       timestep > this%loadedTimestep + this%saveInterval .or.                               &
       (timestep == this%loadedTimestep .and.                                                &
       stage < timeIntegrator%nStages) .or.                                                  &
       (timestep == this%loadedTimestep + this%saveInterval .and.                            &
       stage == timeIntegrator%nStages)) then

     if (mod(timestep, this%saveInterval) == 0 .and. stage == timeIntegrator%nStages) then
        timestep_ = timestep
     else if (mod(timestep, this%saveInterval) == 0) then
        timestep_ = timestep - this%saveInterval
     else
        timestep_ = timestep - mod(timestep, this%saveInterval)
     end if

     write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", timestep_, ".q"
     call region%loadData(QOI_FORWARD_STATE, filename)
     time = region%states(1)%time
     this%loadedTimestep = timestep_

     i = 1
     do j = 1, size(region%states)
        this%data_(j)%buffer(:,:,i) = region%states(j)%conservedVariables
     end do

     if (this%loadedTimestep /= this%endTimestep) then

        do timestep_ = this%loadedTimestep + 1, this%loadedTimestep + this%saveInterval

           do j = 1, size(region%states) !... update state
              region%states(j)%time = time
              call region%states(j)%update(region%grids(j), region%simulationFlags,          &
                   region%solverOptions)
           end do
           timeStepSize = region%getTimeStepSize()

           do stage_ = 1, timeIntegrator%nStages

              if (timestep_ == this%loadedTimestep + this%saveInterval .and.                 &
                   stage_ == timeIntegrator%nStages) exit

              call timeIntegrator%substepForward(region, time,                               &
                   timeStepSize, timestep_, stage_)

              if (stage /= timeIntegrator%nStages) then
                 do j = 1, size(region%states) !... update state
                    region%states(j)%time = time
                    call region%states(j)%update(region%grids(j), region%simulationFlags,    &
                         region%solverOptions)
                 end do
              end if

              i = i + 1
              do j = 1, size(region%states)
                 this%data_(j)%buffer(:,:,i) = region%states(j)%conservedVariables
              end do
           end do

        end do

     end if !... this%loadedTimestep /= this%endTimestep

  end if

  i = (timestep - 1 - this%loadedTimestep) * timeIntegrator%nStages + stage + 1
  assert(i >= 1 .and. i <= this%numIntermediateStates)
  do j = 1, size(region%states)
     region%states(j)%conservedVariables = this%data_(j)%buffer(:,:,i)
  end do

end subroutine uniformCheckpointingMigrateTo
