#include "config.h"

module ThermalActuator_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  use Controller_mod, only : t_Controller

  implicit none
  private

  type, extends(t_Controller), public :: t_ThermalActuator

   contains

     procedure, pass :: setup
     procedure, pass :: cleanup
     procedure, pass :: computeSensitivity
     procedure, pass :: update
     procedure, pass :: updateGradient
     procedure, pass :: hookBeforeTimemarch
     procedure, pass :: hookAfterTimemarch

  end type t_ThermalActuator

contains

  subroutine setup(this, region)

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use Region_mod, only : t_Region
    use ActuatorPatch_mod, only : t_ActuatorPatch

    ! <<< Internal modules >>>
    use InputHelper, only : getOption, getRequiredOption

    implicit none

    ! <<< Arguments >>>
    class(t_ThermalActuator) :: this
    class(t_Region) :: region

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, gradientBufferSize
    class(t_Patch), pointer :: patch => null()

    call this%cleanup()

    if (region%simulationFlags%predictionOnly) return

    gradientBufferSize = getOption("gradient_buffer_size", 1)

    if (allocated(region%patchFactories)) then
       do i = 1, size(region%patchFactories)
          call region%patchFactories(i)%connect(patch)
          if (.not. associated(patch)) cycle
          select type (patch)
          class is (t_ActuatorPatch)
             if (patch%nPatchPoints <= 0) cycle
             SAFE_DEALLOCATE(patch%gradient)
             SAFE_DEALLOCATE(patch%probeBuffer)
             allocate(patch%gradient(patch%nPatchPoints, 1))
             allocate(patch%probeBuffer(patch%nPatchPoints, 1, gradientBufferSize),          &
                  source = 0.0_wp)
          end select
       end do
    end if

  end subroutine setup

  subroutine cleanup(this)

    implicit none

    ! <<< Arguments >>>
    class(t_ThermalActuator) :: this

    call this%cleanupBase()

  end subroutine cleanup

  function computeSensitivity(this, region) result(instantaneousSensitivity)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    implicit none

    ! <<< Arguments >>>
    class(t_ThermalActuator) :: this
    class(t_Region), intent(in) :: region

    ! <<< Result >>>
    real(SCALAR_KIND) :: instantaneousSensitivity

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, nDimensions, ierror
    real(SCALAR_KIND), allocatable :: gradient(:,:)

    instantaneousSensitivity = 0.0_wp

    do i = 1, size(region%grids)

       nDimensions = region%grids(i)%nDimensions
       assert_key(nDimensions, (1, 2, 3))

       assert(allocated(region%grids(i)%controlMollifier))
       assert(allocated(region%states(i)%adjointVariables))

       allocate(gradient(region%grids(i)%nGridPoints, 1))

       gradient(:,1) = region%states(i)%adjointVariables(:,nDimensions+2) *                  &
            region%grids(i)%controlMollifier(:,1)
       instantaneousSensitivity = instantaneousSensitivity +                                 &
            region%grids(i)%computeInnerProduct(gradient, gradient)

       SAFE_DEALLOCATE(gradient)

    end do

    if (region%commGridMasters /= MPI_COMM_NULL)                                             &
         call MPI_Allreduce(MPI_IN_PLACE, instantaneousSensitivity, 1,                       &
         SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

    do i = 1, size(region%grids)
       call MPI_Bcast(instantaneousSensitivity, 1, SCALAR_TYPE_MPI,                          &
            0, region%grids(i)%comm, ierror)
    end do

    this%cachedValue = instantaneousSensitivity

  end function computeSensitivity

  subroutine update(this, region)

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use Region_mod, only : t_Region
    use ActuatorPatch_mod, only : t_ActuatorPatch

    implicit none

    ! <<< Arguments >>>
    class(t_ThermalActuator) :: this
    class(t_Region), intent(in) :: region

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, nDimensions
    class(t_Patch), pointer :: patch => null()

    if (.not. allocated(region%patchFactories)) return

    nDimensions = size(region%globalGridSizes, 1)
    assert_key(nDimensions, (1, 2, 3))

    do i = 1, size(region%patchFactories)
       call region%patchFactories(i)%connect(patch)
       if (.not. associated(patch)) cycle
       do j = 1, size(region%states)
          if (patch%gridIndex /= region%grids(j)%index) cycle
          select type (patch)
             class is (t_ActuatorPatch)

             patch%iProbeBuffer = patch%iProbeBuffer - 1

             assert(patch%iProbeBuffer >= 1)
             assert(patch%iProbeBuffer <= size(patch%probeBuffer, 3))

             if (patch%iProbeBuffer == size(patch%probeBuffer, 3)) call patch%loadData()

             patch%controlForcing(:,1:nDimensions+1) = 0.0_wp
             patch%controlForcing(:,nDimensions+2) = - region%states(j)%actuationAmount *    &
                  patch%probeBuffer(:,1,patch%iProbeBuffer)

             if (patch%iProbeBuffer == 1)                                                    &
                  patch%iProbeBuffer = size(patch%probeBuffer, 3) + 1

          end select
       end do
    end do

  end subroutine update

  subroutine updateGradient(this, region)

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use Region_mod, only : t_Region
    use ActuatorPatch_mod, only : t_ActuatorPatch

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : ADJOINT

    implicit none

    ! <<< Arguments >>>
    class(t_ThermalActuator) :: this
    class(t_Region), intent(in) :: region

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, nDimensions
    class(t_Patch), pointer :: patch => null()

    if (.not. allocated(region%patchFactories)) return

    nDimensions = size(region%globalGridSizes, 1)
    assert_key(nDimensions, (1, 2, 3))

    do i = 1, size(region%patchFactories)
       call region%patchFactories(i)%connect(patch)
       if (.not. associated(patch)) cycle
       do j = 1, size(region%states)
          if (patch%gridIndex /= region%grids(j)%index .or. patch%nPatchPoints <= 0) cycle
          select type (patch)
          class is (t_ActuatorPatch)

             call patch%collect(region%states(j)%adjointVariables(:,nDimensions+2),          &
                  patch%gradient(:,1))
             call patch%collectMultiply(region%grids(j)%controlMollifier(:,1),               &
                  patch%gradient(:,1))

             patch%iProbeBuffer = patch%iProbeBuffer + 1
             assert(patch%iProbeBuffer >= 1)
             assert(patch%iProbeBuffer <= size(patch%probeBuffer, 3))

             call patch%update(ADJOINT, region%simulationFlags, region%solverOptions,        &
                  region%grids(j), region%states(j))

             if (patch%iProbeBuffer == size(patch%probeBuffer, 3)) then
                call patch%saveData()
                patch%iProbeBuffer = 0
             end if

          end select
       end do
    end do

  end subroutine updateGradient

  subroutine hookBeforeTimemarch(this, region, mode)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use Region_mod, only : t_Region
    use ActuatorPatch_mod, only : t_ActuatorPatch

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : FORWARD, ADJOINT

    implicit none

    ! <<< Arguments >>>
    class(t_ThermalActuator) :: this
    class(t_Region) :: region
    integer, intent(in) :: mode

    ! <<< Local variables >>>
    integer :: i
    class(t_Patch), pointer :: patch => null()

    if (.not. allocated(region%patchFactories)) return

    do i = 1, size(region%patchFactories)
       call region%patchFactories(i)%connect(patch)
       if (.not. associated(patch)) cycle
       select type (patch)
       class is (t_ActuatorPatch)

          select case (mode)
          case (FORWARD)
             call patch%seekToEOF()
          case (ADJOINT)
             call patch%reset()
          end select

       end select
    end do

  end subroutine hookBeforeTimemarch

  subroutine hookAfterTimemarch(this, region, mode)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Patch_mod, only : t_Patch
    use Region_mod, only : t_Region
    use ActuatorPatch_mod, only : t_ActuatorPatch

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : ADJOINT

    implicit none

    ! <<< Arguments >>>
    class(t_ThermalActuator) :: this
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
             call patch%saveData()

          end select

       end select
    end do

  end subroutine hookAfterTimemarch

end module ThermalActuator_mod
