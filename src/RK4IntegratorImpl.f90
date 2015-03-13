#include "config.h"

subroutine setupRK4Integrator(this, region)

  ! <<< Derived types >>>
  use RK4Integrator_type
  use Region_type

  ! <<< Public members >>>
  use RK4Integrator_mod, only : cleanupRK4Integrator

  implicit none

  ! <<< Arguments >>>
  type(t_RK4Integrator) :: this
  type(t_Region), intent(in) :: region

  ! <<< Local variables >>>
  integer :: i

  call cleanupRK4Integrator(this)

  allocate(this%temp(size(region%states)))
  do i = 1, size(this%temp)
     allocate(this%temp(i)%buffer1(region%grids(i)%nGridPoints, region%states(i)%nUnknowns))
     allocate(this%temp(i)%buffer2(region%grids(i)%nGridPoints, region%states(i)%nUnknowns))
  end do

end subroutine setupRK4Integrator

subroutine cleanupRK4Integrator(this)

  ! <<< Derived types >>>
  use RK4Integrator_type

  implicit none

  ! <<< Arguments >>>
  type(t_RK4Integrator) :: this

  ! <<< Local variables >>>
  integer :: i

  if (allocated(this%temp)) then
     do i = 1, size(this%temp)
        SAFE_DEALLOCATE(this%temp(i)%buffer1)
        SAFE_DEALLOCATE(this%temp(i)%buffer2)
     end do
  end if
  SAFE_DEALLOCATE(this%temp)

end subroutine cleanupRK4Integrator

subroutine stepForward(this, region, time, timestep, verbose)

  ! <<< External modules >>>
  use, intrinsic :: iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use RK4Integrator_type
  use Region_type

  ! <<< Internal modules >>>
  use MPIHelper, only : writeAndFlush
  use Region_mod, only : computeRhs, subStepHooks

  implicit none

  ! <<< Arguments >>>
  type(t_RK4Integrator) :: this
  type(t_Region) :: region
  real(SCALAR_KIND), intent(inout) :: time
  integer, intent(in) :: timestep
  logical, intent(in), optional :: verbose

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i
  real(wp) :: timeStepSize
  character(len = STRING_LENGTH) :: str

  do i = 1, size(region%states)
     this%temp(i)%buffer1 = region%states(i)%conservedVariables
  end do

  ! Stage 1:

  call computeRhs(region, FORWARD, time)
  call subStepHooks(region, timestep, 0)
  timeStepSize = region%states(1)%timeStepSize

  do i = 1, size(region%states)
     this%temp(i)%buffer2 = region%states(i)%conservedVariables +                            &
          timeStepSize * region%states(i)%rightHandSide / 6.0_wp
     region%states(i)%conservedVariables = this%temp(i)%buffer1 +                            &
          timeStepSize * region%states(i)%rightHandSide / 2.0_wp
  end do

  ! Stage 2:

  time = time + timeStepSize / 2.0_wp
  call computeRhs(region, FORWARD, time)
  call subStepHooks(region, timestep, 1)

  do i = 1, size(region%states)
     this%temp(i)%buffer2 = this%temp(i)%buffer2 +                                           &
          timeStepSize * region%states(i)%rightHandSide / 3.0_wp
     region%states(i)%conservedVariables = this%temp(i)%buffer1 +                            &
          timeStepSize * region%states(i)%rightHandSide / 2.0_wp
  end do

  ! Stage 3:

  call computeRhs(region, FORWARD, time)
  call subStepHooks(region, timestep, 2)

  do i = 1, size(region%states)
     this%temp(i)%buffer2 = this%temp(i)%buffer2 +                                           &
          timeStepSize * region%states(i)%rightHandSide / 3.0_wp
     region%states(i)%conservedVariables = this%temp(i)%buffer1 +                            &
          timeStepSize * region%states(i)%rightHandSide
  end do

  ! Stage 4:

  time = time + timeStepSize / 2.0_wp
  call computeRhs(region, FORWARD, time)
  call subStepHooks(region, timestep, 3)

  do i = 1, size(region%states)
     region%states(i)%conservedVariables = this%temp(i)%buffer2 +                            &
          timeStepSize * region%states(i)%rightHandSide / 6.0_wp
     region%states(i)%plot3dAuxiliaryData(1) = real(timestep, wp)
     region%states(i)%plot3dAuxiliaryData(4) = time
  end do

  call subStepHooks(region, timestep, 4)

  if (present(verbose)) then
     if (verbose) then
        if (region%simulationFlags%useConstantCfl) then
           write(str, '(2A,I8,2(A,D13.6))') PROJECT_NAME, ": timestep = ", timestep,         &
                ", dt = ", timeStepSize, ", time = ", time
        else
           write(str, '(2A,I8,2(A,D13.6))') PROJECT_NAME, ": timestep = ", timestep,         &
                ", CFL = ", region%states(1)%cfl, ", time = ", time
        end if
        call writeAndFlush(region%comm, output_unit, str)
     end if
  end if

end subroutine stepForward
