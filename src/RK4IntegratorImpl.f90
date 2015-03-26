#include "config.h"

subroutine setupRK4Integrator(this, region)

  ! <<< Derived types >>>
  use Region_type, only : t_Region
  use RK4Integrator_type, only : t_RK4Integrator

  ! <<< Public members >>>
  use RK4Integrator_mod, only : cleanupRK4Integrator

  implicit none

  ! <<< Arguments >>>
  type(t_RK4Integrator) :: this
  type(t_Region), intent(in) :: region

  ! <<< Local variables >>>
  integer :: i

  call cleanupRK4Integrator(this)

  assert(allocated(region%states))
  assert(size(region%states) > 0)

  allocate(this%temp_(size(region%states)))
  do i = 1, size(this%temp_)
     assert(region%grids(i)%nGridPoints > 0)
     assert(region%states(i)%nUnknowns > 0)
     allocate(this%temp_(i)%buffer1(region%grids(i)%nGridPoints, region%states(i)%nUnknowns))
     allocate(this%temp_(i)%buffer2(region%grids(i)%nGridPoints, region%states(i)%nUnknowns))
  end do

end subroutine setupRK4Integrator

subroutine cleanupRK4Integrator(this)

  ! <<< Derived types >>>
  use RK4Integrator_type, only : t_RK4Integrator

  implicit none

  ! <<< Arguments >>>
  type(t_RK4Integrator) :: this

  ! <<< Local variables >>>
  integer :: i

  if (allocated(this%temp_)) then
     do i = 1, size(this%temp_)
        SAFE_DEALLOCATE(this%temp_(i)%buffer1)
        SAFE_DEALLOCATE(this%temp_(i)%buffer2)
     end do
  end if
  SAFE_DEALLOCATE(this%temp_)

end subroutine cleanupRK4Integrator

subroutine subStepForward(this, region, time, timestep, stage)

  ! <<< Derived types >>>
  use Region_type, only : t_Region, FORWARD
  use RK4Integrator_type, only : t_RK4Integrator

  ! <<< Internal modules >>>
  use Region_mod, only : computeRhs

  implicit none

  ! <<< Arguments >>>
  type(t_RK4Integrator) :: this
  type(t_Region) :: region
  real(SCALAR_KIND), intent(inout) :: time
  integer, intent(in) :: timestep, stage

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i
  real(wp), save :: timeStepSize

  assert(timestep >= 0)
  assert(stage >= 1 .and. stage <= 4)

  select case (stage)

  case (1)

     do i = 1, size(region%states)
        this%temp_(i)%buffer1 = region%states(i)%conservedVariables
     end do

     call computeRhs(region, FORWARD, time)
     timeStepSize = region%states(1)%timeStepSize

     do i = 1, size(region%states)
        this%temp_(i)%buffer2 = region%states(i)%conservedVariables +                        &
             timeStepSize * region%states(i)%rightHandSide / 6.0_wp
        region%states(i)%conservedVariables = this%temp_(i)%buffer1 +                        &
             timeStepSize * region%states(i)%rightHandSide / 2.0_wp
     end do

  case (2)

     time = time + timeStepSize / 2.0_wp
     call computeRhs(region, FORWARD, time)

     do i = 1, size(region%states)
        this%temp_(i)%buffer2 = this%temp_(i)%buffer2 +                                      &
             timeStepSize * region%states(i)%rightHandSide / 3.0_wp
        region%states(i)%conservedVariables = this%temp_(i)%buffer1 +                        &
             timeStepSize * region%states(i)%rightHandSide / 2.0_wp
     end do

  case (3)

     call computeRhs(region, FORWARD, time)

     do i = 1, size(region%states)
        this%temp_(i)%buffer2 = this%temp_(i)%buffer2 +                                      &
             timeStepSize * region%states(i)%rightHandSide / 3.0_wp
        region%states(i)%conservedVariables = this%temp_(i)%buffer1 +                        &
             timeStepSize * region%states(i)%rightHandSide
     end do

  case (4)

     time = time + timeStepSize / 2.0_wp
     call computeRhs(region, FORWARD, time)

     do i = 1, size(region%states)
        region%states(i)%conservedVariables = this%temp_(i)%buffer2 +                        &
             timeStepSize * region%states(i)%rightHandSide / 6.0_wp
     end do

  end select

end subroutine subStepForward

subroutine subStepAdjoint(this, region, time, timestep, stage)

  ! <<< Derived types >>>
  use Region_type, only : t_Region, ADJOINT
  use RK4Integrator_type, only : t_RK4Integrator

  ! <<< Internal modules >>>
  use Region_mod, only : computeRhs

  implicit none

  ! <<< Arguments >>>
  type(t_RK4Integrator) :: this
  type(t_Region) :: region
  real(SCALAR_KIND), intent(inout) :: time
  integer, intent(in) :: timestep, stage

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i
  real(wp), save :: timeStepSize

  assert(timestep >= 0)
  assert(stage >= 1 .and. stage <= 4)

  select case (stage)

  case (4)

     do i = 1, size(region%states)
        this%temp_(i)%buffer1 = region%states(i)%adjointVariables
     end do

     call computeRhs(region, ADJOINT, time)
     timeStepSize = region%states(1)%timeStepSize

     do i = 1, size(region%states)
        this%temp_(i)%buffer2 = region%states(i)%adjointVariables -                          &
             timeStepSize * region%states(i)%rightHandSide / 6.0_wp
        region%states(i)%adjointVariables = this%temp_(i)%buffer1 -                          &
             timeStepSize * region%states(i)%rightHandSide / 2.0_wp
     end do

  case (3)

     time = time - timeStepSize / 2.0_wp
     call computeRhs(region, ADJOINT, time)

     do i = 1, size(region%states)
        this%temp_(i)%buffer2 = this%temp_(i)%buffer2 -                                      &
             timeStepSize * region%states(i)%rightHandSide / 3.0_wp
        region%states(i)%adjointVariables = this%temp_(i)%buffer1 -                          &
             timeStepSize * region%states(i)%rightHandSide / 2.0_wp
     end do

  case (2)

     call computeRhs(region, ADJOINT, time)

     do i = 1, size(region%states)
        this%temp_(i)%buffer2 = this%temp_(i)%buffer2 -                                      &
             timeStepSize * region%states(i)%rightHandSide / 3.0_wp
        region%states(i)%adjointVariables = this%temp_(i)%buffer1 -                          &
             timeStepSize * region%states(i)%rightHandSide
     end do

  case (1)

     time = time - timeStepSize / 2.0_wp
     call computeRhs(region, ADJOINT, time)

     do i = 1, size(region%states)
        region%states(i)%adjointVariables = this%temp_(i)%buffer2 -                          &
             timeStepSize * region%states(i)%rightHandSide / 6.0_wp
     end do

  end select

end subroutine subStepAdjoint

subroutine stepForward(this, region, time, timestep, verbose)

  ! <<< External modules >>>
  use, intrinsic :: iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Region_type, only : t_Region
  use RK4Integrator_type, only : t_RK4Integrator

  ! <<< Public members >>>
  use RK4Integrator_mod, only : subStepForward

  ! <<< Internal modules >>>
  use Region_mod, only : reportResiduals
  use ErrorHandler, only : writeAndFlush

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
  character(len = STRING_LENGTH) :: str

  assert(timestep >= 0)

  do i = 1, 4
     call subStepForward(this, region, time, timestep, i)
  end do

  do i = 1, size(region%states)
     region%states(i)%plot3dAuxiliaryData(1) = real(timestep, wp)
     region%states(i)%plot3dAuxiliaryData(4) = time
  end do

  if (present(verbose)) then
     if (verbose) then
        if (region%simulationFlags%useConstantCfl) then
           write(str, '(2A,I8,2(A,D13.6))') PROJECT_NAME, ": timestep = ", timestep,         &
                ", dt = ", region%states(1)%timeStepSize, ", time = ", time
        else
           write(str, '(2A,I8,2(A,D13.6))') PROJECT_NAME, ": timestep = ", timestep,         &
                ", CFL = ", region%states(1)%cfl, ", time = ", time
        end if
        call writeAndFlush(region%comm, output_unit, str)
        if (region%simulationFlags%steadyStateSimulation) call reportResiduals(region)
     end if
  end if

end subroutine stepForward

subroutine stepAdjoint(this, region, time, timestep, verbose)

  ! <<< External modules >>>
  use, intrinsic :: iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Region_type, only : t_Region
  use RK4Integrator_type, only : t_RK4Integrator

  ! <<< Internal modules >>>
  use Region_mod, only : computeRhs, reportResiduals
  use ErrorHandler, only : writeAndFlush

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
  character(len = STRING_LENGTH) :: str

  assert(timestep > 0)

  do i = 4, 1, -1
     call subStepAdjoint(this, region, time, timestep, i)
  end do

  do i = 1, size(region%states)
     region%states(i)%plot3dAuxiliaryData(1) = real(timestep, wp)
     region%states(i)%plot3dAuxiliaryData(4) = time
  end do

  if (present(verbose)) then
     if (verbose) then
        if (region%simulationFlags%useConstantCfl) then
           write(str, '(2A,I8,2(A,D13.6))') PROJECT_NAME, ": timestep = ", timestep,         &
                ", dt = ", region%states(1)%timeStepSize, ", time = ", time
        else
           write(str, '(2A,I8,2(A,D13.6))') PROJECT_NAME, ": timestep = ", timestep,         &
                ", CFL = ", region%states(1)%cfl, ", time = ", time
        end if
        call writeAndFlush(region%comm, output_unit, str)
        if (region%simulationFlags%steadyStateSimulation) call reportResiduals(region)
     end if
  end if

end subroutine stepAdjoint
