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

  allocate(this%temp_(size(region%states)))
  do i = 1, size(this%temp_)
     allocate(this%temp_(i)%buffer1(region%grids(i)%nGridPoints, region%states(i)%nUnknowns))
     allocate(this%temp_(i)%buffer2(region%grids(i)%nGridPoints, region%states(i)%nUnknowns))
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

  if (allocated(this%temp_)) then
     do i = 1, size(this%temp_)
        SAFE_DEALLOCATE(this%temp_(i)%buffer1)
        SAFE_DEALLOCATE(this%temp_(i)%buffer2)
     end do
  end if
  SAFE_DEALLOCATE(this%temp_)

end subroutine cleanupRK4Integrator

subroutine stepForward(this, region, time, timestep, verbose)

  ! <<< External modules >>>
  use, intrinsic :: iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use RK4Integrator_type
  use Region_type

  ! <<< Internal modules >>>
  use MPIHelper, only : writeAndFlush
  use Region_mod, only : computeRhs, subStepHooks, reportResiduals

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
     this%temp_(i)%buffer1 = region%states(i)%conservedVariables
  end do

  ! Stage 1:

  call computeRhs(region, FORWARD, time)
  call subStepHooks(region, FORWARD, timestep, 0)
  timeStepSize = region%states(1)%timeStepSize

  do i = 1, size(region%states)
     this%temp_(i)%buffer2 = region%states(i)%conservedVariables +                           &
          timeStepSize * region%states(i)%rightHandSide / 6.0_wp
     region%states(i)%conservedVariables = this%temp_(i)%buffer1 +                           &
          timeStepSize * region%states(i)%rightHandSide / 2.0_wp
  end do

  ! Stage 2:

  time = time + timeStepSize / 2.0_wp
  call computeRhs(region, FORWARD, time)
  call subStepHooks(region, FORWARD, timestep, 1)

  do i = 1, size(region%states)
     this%temp_(i)%buffer2 = this%temp_(i)%buffer2 +                                         &
          timeStepSize * region%states(i)%rightHandSide / 3.0_wp
     region%states(i)%conservedVariables = this%temp_(i)%buffer1 +                           &
          timeStepSize * region%states(i)%rightHandSide / 2.0_wp
  end do

  ! Stage 3:

  call computeRhs(region, FORWARD, time)
  call subStepHooks(region, FORWARD, timestep, 2)

  do i = 1, size(region%states)
     this%temp_(i)%buffer2 = this%temp_(i)%buffer2 +                                         &
          timeStepSize * region%states(i)%rightHandSide / 3.0_wp
     region%states(i)%conservedVariables = this%temp_(i)%buffer1 +                           &
          timeStepSize * region%states(i)%rightHandSide
  end do

  ! Stage 4:

  time = time + timeStepSize / 2.0_wp
  call computeRhs(region, FORWARD, time)
  call subStepHooks(region, FORWARD, timestep, 3)

  do i = 1, size(region%states)
     region%states(i)%conservedVariables = this%temp_(i)%buffer2 +                           &
          timeStepSize * region%states(i)%rightHandSide / 6.0_wp
     region%states(i)%plot3dAuxiliaryData(1) = real(timestep, wp)
     region%states(i)%plot3dAuxiliaryData(4) = time
  end do

  call subStepHooks(region, FORWARD, timestep, 4)

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
        if (region%simulationFlags%steadyStateSimulation) call reportResiduals(region)
     end if
  end if

end subroutine stepForward

subroutine stepAdjoint(this, region, time, timestep, verbose)

  ! <<< External modules >>>
  use, intrinsic :: iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use RK4Integrator_type
  use Region_type

  ! <<< Internal modules >>>
  use MPIHelper, only : writeAndFlush
  use Region_mod, only : computeRhs, subStepHooks, reportResiduals

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
     this%temp_(i)%buffer1 = region%states(i)%adjointVariables
  end do

  ! Stage 4:

  call subStepHooks(region, ADJOINT, timestep, 4)
  call computeRhs(region, ADJOINT, time)
  timeStepSize = region%states(1)%timeStepSize

  do i = 1, size(region%states)
     this%temp_(i)%buffer2 = region%states(i)%adjointVariables -                             &
          timeStepSize * region%states(i)%rightHandSide / 6.0_wp
     region%states(i)%adjointVariables = this%temp_(i)%buffer1 -                             &
          timeStepSize * region%states(i)%rightHandSide / 2.0_wp
  end do

  ! Stage 3:

  time = time - timeStepSize / 2.0_wp
  call subStepHooks(region, ADJOINT, timestep, 3)
  call computeRhs(region, ADJOINT, time)

  do i = 1, size(region%states)
     this%temp_(i)%buffer2 = this%temp_(i)%buffer2 -                                         &
          timeStepSize * region%states(i)%rightHandSide / 3.0_wp
     region%states(i)%adjointVariables = this%temp_(i)%buffer1 -                             &
          timeStepSize * region%states(i)%rightHandSide / 2.0_wp
  end do

  ! Stage 2:

  call subStepHooks(region, ADJOINT, timestep, 2)
  call computeRhs(region, ADJOINT, time)

  do i = 1, size(region%states)
     this%temp_(i)%buffer2 = this%temp_(i)%buffer2 -                                         &
          timeStepSize * region%states(i)%rightHandSide / 3.0_wp
     region%states(i)%adjointVariables = this%temp_(i)%buffer1 -                             &
          timeStepSize * region%states(i)%rightHandSide
  end do

  ! Stage 1:

  time = time - timeStepSize / 2.0_wp
  call subStepHooks(region, ADJOINT, timestep, 1)
  call computeRhs(region, ADJOINT, time)

  do i = 1, size(region%states)
     region%states(i)%adjointVariables = this%temp_(i)%buffer2 -                             &
          timeStepSize * region%states(i)%rightHandSide / 6.0_wp
     region%states(i)%plot3dAuxiliaryData(1) = real(timestep, wp)
     region%states(i)%plot3dAuxiliaryData(4) = time
  end do

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
        if (region%simulationFlags%steadyStateSimulation) call reportResiduals(region)
     end if
  end if

end subroutine stepAdjoint
