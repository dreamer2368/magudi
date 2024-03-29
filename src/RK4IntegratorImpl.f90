#include "config.h"

subroutine setupRK4Integrator(this, region)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use RK4Integrator_mod, only : t_RK4Integrator

  implicit none

  ! <<< Arguments >>>
  class(t_RK4Integrator) :: this
  class(t_Region), intent(in) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i

  assert(allocated(region%states))
  assert(size(region%states) > 0)

  call this%cleanup()

  this%nStages = 4
  call this%setupBase()
  this%norm = (/ 1.0_wp / 6.0_wp, 1.0_wp / 3.0_wp, 1.0_wp / 3.0_wp, 1.0_wp / 6.0_wp /)

  allocate(this%data_(size(region%states)))
  do i = 1, size(this%data_)
     assert(region%grids(i)%nGridPoints > 0)
     assert(region%solverOptions%nUnknowns > 0)
     allocate(this%data_(i)%buffer1(region%grids(i)%nGridPoints,                             &
          region%solverOptions%nUnknowns))
     allocate(this%data_(i)%buffer2(region%grids(i)%nGridPoints,                             &
          region%solverOptions%nUnknowns))
  end do

end subroutine setupRK4Integrator

subroutine cleanupRK4Integrator(this)

  ! <<< Derived types >>>
  use RK4Integrator_mod, only : t_RK4Integrator

  implicit none

  ! <<< Arguments >>>
  class(t_RK4Integrator) :: this

  ! <<< Local variables >>>
  integer :: i

  call this%cleanupBase()

  if (allocated(this%data_)) then
     do i = 1, size(this%data_)
        SAFE_DEALLOCATE(this%data_(i)%buffer1)
        SAFE_DEALLOCATE(this%data_(i)%buffer2)
     end do
  end if
  SAFE_DEALLOCATE(this%data_)

end subroutine cleanupRK4Integrator

subroutine substepForwardRK4(this, region, time, timeStepSize, timestep, stage)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use RK4Integrator_mod, only : t_RK4Integrator

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_RK4Integrator) :: this
  class(t_Region) :: region
  real(SCALAR_KIND), intent(inout) :: time
  real(SCALAR_KIND), intent(in) :: timeStepSize
  integer, intent(in) :: timestep, stage

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer, save :: stageLastCall = 0
  integer :: i

  assert(timestep >= 0)
  assert(stage >= 1 .and. stage <= 4)

#ifdef DEBUG
  if (stageLastCall /= 0) then
     assert(stage == stageLastCall + 1 .or. (stageLastCall == 4 .and. stage == 1))
  end if
#endif

  call startTiming("substepForward")

  stageLastCall = stage

  select case (stage)

  case (1)

     do i = 1, size(region%states)
        this%data_(i)%buffer1 = region%states(i)%conservedVariables
     end do

     region%states(:)%timeProgressive = time + timeStepSize / 2.0_wp
     call region%computeRhs(FORWARD, timestep, stage)

     do i = 1, size(region%states)
        this%data_(i)%buffer2 = region%states(i)%conservedVariables +                        &
             timeStepSize * region%states(i)%rightHandSide / 6.0_wp
        region%states(i)%conservedVariables = this%data_(i)%buffer1 +                        &
             timeStepSize * region%states(i)%rightHandSide / 2.0_wp
     end do

  case (2)

     time = time + timeStepSize / 2.0_wp
     region%states(:)%time = time
     call region%computeRhs(FORWARD, timestep, stage)

     do i = 1, size(region%states)
        this%data_(i)%buffer2 = this%data_(i)%buffer2 +                                      &
             timeStepSize * region%states(i)%rightHandSide / 3.0_wp
        region%states(i)%conservedVariables = this%data_(i)%buffer1 +                        &
             timeStepSize * region%states(i)%rightHandSide / 2.0_wp
     end do

  case (3)

     region%states(:)%timeProgressive = time + timeStepSize / 2.0_wp
     call region%computeRhs(FORWARD, timestep, stage)

     do i = 1, size(region%states)
        this%data_(i)%buffer2 = this%data_(i)%buffer2 +                                      &
             timeStepSize * region%states(i)%rightHandSide / 3.0_wp
        region%states(i)%conservedVariables = this%data_(i)%buffer1 +                        &
             timeStepSize * region%states(i)%rightHandSide
     end do

  case (4)

     time = time + timeStepSize / 2.0_wp
     region%states(:)%time = time
     call region%computeRhs(FORWARD, timestep, stage)

     do i = 1, size(region%states)
        region%states(i)%conservedVariables = this%data_(i)%buffer2 +                        &
             timeStepSize * region%states(i)%rightHandSide / 6.0_wp
     end do

  end select

  call endTiming("substepForward")

end subroutine substepForwardRK4

subroutine substepAdjointRK4(this, region, time, timeStepSize, timestep, stage)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use RK4Integrator_mod, only : t_RK4Integrator

  ! <<< Enumerations >>>
  use Region_enum, only : ADJOINT

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_RK4Integrator) :: this
  class(t_Region) :: region
  real(SCALAR_KIND), intent(inout) :: time
  real(SCALAR_KIND), intent(in) :: timeStepSize
  integer, intent(in) :: timestep, stage

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer, save :: stageLastCall = 0
  integer :: i

  assert(timestep >= 0)
  assert(stage >= 1 .and. stage <= 4)

#ifdef DEBUG
  if (stageLastCall /= 0) then
     assert(stage == stageLastCall - 1 .or. (stageLastCall == 1 .and. stage == 4))
  end if
#endif

  call startTiming("substepAdjoint")

  stageLastCall = stage

  select case (stage)

  case (4)

     do i = 1, size(region%states)
        this%data_(i)%buffer1 = region%states(i)%adjointVariables
     end do

     region%states(:)%adjointForcingFactor = 2.0_wp
     region%states(:)%timeProgressive = time - timeStepSize / 2.0_wp
     call region%computeRhs(ADJOINT, timestep, stage)

     do i = 1, size(region%states)
        this%data_(i)%buffer2 = region%states(i)%adjointVariables -                          &
             timeStepSize * region%states(i)%rightHandSide / 6.0_wp
        region%states(i)%adjointVariables = this%data_(i)%buffer1 -                          &
             timeStepSize * region%states(i)%rightHandSide / 2.0_wp
     end do

     region%states(:)%timeProgressive = time

  case (3)

     region%states(:)%adjointForcingFactor = 1.0_wp
     call region%computeRhs(ADJOINT, timestep, stage)

     do i = 1, size(region%states)
        this%data_(i)%buffer2 = this%data_(i)%buffer2 -                                      &
             timeStepSize * region%states(i)%rightHandSide / 3.0_wp
        region%states(i)%adjointVariables = this%data_(i)%buffer1 -                          &
             timeStepSize * region%states(i)%rightHandSide / 2.0_wp
     end do

     time = time - timeStepSize / 2.0_wp
     region%states(:)%time = time

  case (2)

     region%states(:)%adjointForcingFactor = 0.5_wp
     call region%computeRhs(ADJOINT, timestep, stage)

     do i = 1, size(region%states)
        this%data_(i)%buffer2 = this%data_(i)%buffer2 -                                      &
             timeStepSize * region%states(i)%rightHandSide / 3.0_wp
        region%states(i)%adjointVariables = this%data_(i)%buffer1 -                          &
             timeStepSize * region%states(i)%rightHandSide
     end do

     region%states(:)%timeProgressive = time

  case (1)

     region%states(:)%adjointForcingFactor = 1.0_wp
     call region%computeRhs(ADJOINT, timestep, stage)

     do i = 1, size(region%states)
        region%states(i)%adjointVariables = this%data_(i)%buffer2 -                          &
             timeStepSize * region%states(i)%rightHandSide / 6.0_wp
     end do

     time = time - timeStepSize / 2.0_wp
     region%states(:)%time = time

  end select

  call endTiming("substepAdjoint")

end subroutine substepAdjointRK4

subroutine substepLinearizedRK4(this, region, time, timeStepSize, timestep, stage)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use RK4Integrator_mod, only : t_RK4Integrator

  ! <<< Enumerations >>>
  use Region_enum, only : LINEARIZED

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_RK4Integrator) :: this
  class(t_Region) :: region
  real(SCALAR_KIND), intent(inout) :: time
  real(SCALAR_KIND), intent(in) :: timeStepSize
  integer, intent(in) :: timestep, stage

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer, save :: stageLastCall = 0
  integer :: i

  assert(timestep >= 0)
  assert(stage >= 1 .and. stage <= 4)

#ifdef DEBUG
  if (stageLastCall /= 0) then
     assert(stage == stageLastCall + 1 .or. (stageLastCall == 4 .and. stage == 1))
  end if
#endif

  call startTiming("substepLinearized")

  stageLastCall = stage

  select case (stage)

  case (1)

     do i = 1, size(region%states)
        this%data_(i)%buffer1 = region%states(i)%adjointVariables
     end do

     region%states(:)%timeProgressive = time + timeStepSize / 2.0_wp
     call region%computeRhs(LINEARIZED, timestep, stage)

     do i = 1, size(region%states)
        this%data_(i)%buffer2 = region%states(i)%adjointVariables +                        &
             timeStepSize * region%states(i)%rightHandSide / 6.0_wp
        region%states(i)%adjointVariables = this%data_(i)%buffer1 +                        &
             timeStepSize * region%states(i)%rightHandSide / 2.0_wp
     end do

  case (2)

     time = time + timeStepSize / 2.0_wp
     region%states(:)%time = time
     call region%computeRhs(LINEARIZED, timestep, stage)

     do i = 1, size(region%states)
        this%data_(i)%buffer2 = this%data_(i)%buffer2 +                                      &
             timeStepSize * region%states(i)%rightHandSide / 3.0_wp
        region%states(i)%adjointVariables = this%data_(i)%buffer1 +                        &
             timeStepSize * region%states(i)%rightHandSide / 2.0_wp
     end do

  case (3)

     region%states(:)%timeProgressive = time + timeStepSize / 2.0_wp
     call region%computeRhs(LINEARIZED, timestep, stage)

     do i = 1, size(region%states)
        this%data_(i)%buffer2 = this%data_(i)%buffer2 +                                      &
             timeStepSize * region%states(i)%rightHandSide / 3.0_wp
        region%states(i)%adjointVariables = this%data_(i)%buffer1 +                        &
             timeStepSize * region%states(i)%rightHandSide
     end do

  case (4)

     time = time + timeStepSize / 2.0_wp
     region%states(:)%time = time
     call region%computeRhs(LINEARIZED, timestep, stage)

     do i = 1, size(region%states)
        region%states(i)%adjointVariables = this%data_(i)%buffer2 +                        &
             timeStepSize * region%states(i)%rightHandSide / 6.0_wp
     end do

  end select

  call endTiming("substepLinearized")

end subroutine substepLinearizedRK4
