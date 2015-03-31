#include "config.h"

subroutine setupJamesonRK3Integrator(this, region)

  ! <<< Derived types >>>
  use Region_type, only : t_Region
  use JamesonRK3Integrator_mod, only : t_JamesonRK3Integrator

  implicit none

  ! <<< Arguments >>>
  class(t_JamesonRK3Integrator) :: this
  class(t_Region), intent(in) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i

  assert(allocated(region%states))
  assert(size(region%states) > 0)

  call this%cleanup()

  this%nStages = 3
  call this%setupBase()
  this%norm = (/ 0.0_wp, 0.0_wp, 1.0_wp /)

  allocate(this%temp_(size(region%states)))
  do i = 1, size(this%temp_)
     assert(region%grids(i)%nGridPoints > 0)
     assert(region%states(i)%nUnknowns > 0)
     allocate(this%temp_(i)%buffer1(region%grids(i)%nGridPoints, region%states(i)%nUnknowns))
     allocate(this%temp_(i)%buffer2(region%grids(i)%nGridPoints, region%states(i)%nUnknowns))
  end do

end subroutine setupJamesonRK3Integrator

subroutine cleanupJamesonRK3Integrator(this)

  ! <<< Derived types >>>
  use JamesonRK3Integrator_mod, only : t_JamesonRK3Integrator

  implicit none

  ! <<< Arguments >>>
  class(t_JamesonRK3Integrator) :: this

  ! <<< Local variables >>>
  integer :: i

  call this%cleanupBase()

  if (allocated(this%temp_)) then
     do i = 1, size(this%temp_)
        SAFE_DEALLOCATE(this%temp_(i)%buffer1)
        SAFE_DEALLOCATE(this%temp_(i)%buffer2)
     end do
  end if
  SAFE_DEALLOCATE(this%temp_)

end subroutine cleanupJamesonRK3Integrator

subroutine substepForwardJamesonRK3(this, region, time, timeStepSize, timestep, stage)

  ! <<< Derived types >>>
  use Region_type, only : t_Region, FORWARD
  use JamesonRK3Integrator_mod, only : t_JamesonRK3Integrator

  ! <<< Internal modules >>>
  use Region_mod, only : computeRhs
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_JamesonRK3Integrator) :: this
  class(t_Region) :: region
  real(SCALAR_KIND), intent(inout) :: time
  real(SCALAR_KIND), intent(in) :: timeStepSize
  integer, intent(in) :: timestep, stage

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer, save :: stageLastCall = 0
  integer :: i

  assert(timestep >= 0)
  assert(stage >= 1 .and. stage <= 3)

#ifdef DEBUG
  if (stageLastCall /= 0) then
     assert(stage == stageLastCall + 1 .or. (stageLastCall == 3 .and. stage == 1))
  end if
#endif

  call startTiming("substepForward")

  stageLastCall = stage

  select case (stage)

  case (1)

     do i = 1, size(region%states)
        this%temp_(i)%buffer1 = region%states(i)%conservedVariables        
     end do

     call computeRhs(region, FORWARD, time)

     do i = 1, size(region%states)
        region%states(i)%conservedVariables = this%temp_(i)%buffer1 +                        &
             timeStepSize * region%states(i)%rightHandSide
        this%temp_(i)%buffer2 = region%states(i)%conservedVariables
     end do

  case (2)

     time = time + timeStepSize / 2.0_wp
     call computeRhs(region, FORWARD, time)

     do i = 1, size(region%states)
        region%states(i)%conservedVariables = (this%temp_(i)%buffer1 +                       &
             region%states(i)%conservedVariables) / 2.0_wp +                                 &
             timeStepSize * region%states(i)%rightHandSide / 2.0_wp
     end do

  case (3)

     time = time + timeStepSize / 2.0_wp
     call computeRhs(region, FORWARD, time)

     do i = 1, size(region%states)
        region%states(i)%conservedVariables =                                                &
             (this%temp_(i)%buffer1 + this%temp_(i)%buffer2) / 2.0_wp +                      &
             timeStepSize * region%states(i)%rightHandSide / 2.0_wp
     end do

  end select

  call endTiming("substepForward")

end subroutine substepForwardJamesonRK3

subroutine substepAdjointJamesonRK3(this, region, time, timeStepSize, timestep, stage)

  ! <<< Derived types >>>
  use Region_type, only : t_Region, ADJOINT
  use JamesonRK3Integrator_mod, only : t_JamesonRK3Integrator

  implicit none

  ! <<< Arguments >>>
  class(t_JamesonRK3Integrator) :: this
  class(t_Region) :: region
  real(SCALAR_KIND), intent(inout) :: time
  real(SCALAR_KIND), intent(in) :: timeStepSize
  integer, intent(in) :: timestep, stage

end subroutine substepAdjointJamesonRK3
