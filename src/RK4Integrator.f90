#include "config.h"

module RK4Integrator_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  use TimeIntegrator_mod, only : t_TimeIntegrator

  implicit none
  private

  type :: t_IntermediateStorage
     real(SCALAR_KIND), allocatable :: buffer1(:,:), buffer2(:,:)
  end type t_IntermediateStorage

  type, extends(t_TimeIntegrator), public :: t_RK4Integrator

     type(t_IntermediateStorage), allocatable :: data_(:)

   contains

     procedure, pass :: setup
     procedure, pass :: cleanup
     procedure, pass :: substepForward
     procedure, pass :: substepAdjoint

  end type t_RK4Integrator

contains

  subroutine setup(this, region)

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

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
       allocate(this%data_(i)%buffer1(region%grids(i)%nGridPoints,                           &
            region%solverOptions%nUnknowns))
       allocate(this%data_(i)%buffer2(region%grids(i)%nGridPoints,                           &
            region%solverOptions%nUnknowns))
    end do

  end subroutine setup

  subroutine cleanup(this)

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

  end subroutine cleanup

  subroutine substepForward(this, region, time, timeStepSize, timestep, stage)

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : FORWARD

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

#ifndef NDEBUG
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

       call region%computeRhs(FORWARD)

       do i = 1, size(region%states)
          this%data_(i)%buffer2 = region%states(i)%conservedVariables +                      &
               timeStepSize * region%states(i)%rightHandSide / 6.0_wp
          region%states(i)%conservedVariables = this%data_(i)%buffer1 +                      &
               timeStepSize * region%states(i)%rightHandSide / 2.0_wp
       end do

    case (2)

       time = time + timeStepSize / 2.0_wp
       region%states(:)%time = time
       call region%computeRhs(FORWARD)

       do i = 1, size(region%states)
          this%data_(i)%buffer2 = this%data_(i)%buffer2 +                                    &
               timeStepSize * region%states(i)%rightHandSide / 3.0_wp
          region%states(i)%conservedVariables = this%data_(i)%buffer1 +                      &
               timeStepSize * region%states(i)%rightHandSide / 2.0_wp
       end do

    case (3)

       call region%computeRhs(FORWARD)

       do i = 1, size(region%states)
          this%data_(i)%buffer2 = this%data_(i)%buffer2 +                                    &
               timeStepSize * region%states(i)%rightHandSide / 3.0_wp
          region%states(i)%conservedVariables = this%data_(i)%buffer1 +                      &
               timeStepSize * region%states(i)%rightHandSide
       end do

    case (4)

       time = time + timeStepSize / 2.0_wp
       region%states(:)%time = time
       call region%computeRhs(FORWARD)

       do i = 1, size(region%states)
          region%states(i)%conservedVariables = this%data_(i)%buffer2 +                      &
               timeStepSize * region%states(i)%rightHandSide / 6.0_wp
       end do

    end select

    call endTiming("substepForward")

  end subroutine substepForward

  subroutine substepAdjoint(this, region, time, timeStepSize, timestep, stage)

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : ADJOINT

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

#ifndef NDEBUG
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
       call region%computeRhs(ADJOINT)

       do i = 1, size(region%states)
          this%data_(i)%buffer2 = region%states(i)%adjointVariables -                        &
               timeStepSize * region%states(i)%rightHandSide / 6.0_wp
          region%states(i)%adjointVariables = this%data_(i)%buffer1 -                        &
               timeStepSize * region%states(i)%rightHandSide / 2.0_wp
       end do

    case (3)

       time = time - timeStepSize / 2.0_wp
       region%states(:)%time = time
       region%states(:)%adjointForcingFactor = 1.0_wp
       call region%computeRhs(ADJOINT)

       do i = 1, size(region%states)
          this%data_(i)%buffer2 = this%data_(i)%buffer2 -                                    &
               timeStepSize * region%states(i)%rightHandSide / 3.0_wp
          region%states(i)%adjointVariables = this%data_(i)%buffer1 -                        &
               timeStepSize * region%states(i)%rightHandSide / 2.0_wp
       end do

    case (2)

       region%states(:)%adjointForcingFactor = 0.5_wp
       call region%computeRhs(ADJOINT)

       do i = 1, size(region%states)
          this%data_(i)%buffer2 = this%data_(i)%buffer2 -                                    &
               timeStepSize * region%states(i)%rightHandSide / 3.0_wp
          region%states(i)%adjointVariables = this%data_(i)%buffer1 -                        &
               timeStepSize * region%states(i)%rightHandSide
       end do

    case (1)

       time = time - timeStepSize / 2.0_wp
       region%states(:)%time = time
       region%states(:)%adjointForcingFactor = 1.0_wp
       call region%computeRhs(ADJOINT)

       do i = 1, size(region%states)
          region%states(i)%adjointVariables = this%data_(i)%buffer2 -                        &
               timeStepSize * region%states(i)%rightHandSide / 6.0_wp
       end do

    end select

    call endTiming("substepAdjoint")

  end subroutine substepAdjoint

end module RK4Integrator_mod
