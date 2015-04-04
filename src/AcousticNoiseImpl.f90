#include "config.h"

module AcousticNoiseImpl

  implicit none
  public

contains

  subroutine loadMeanPressure(this, filename, region)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env, only : output_unit

    ! <<< Derived types >>>
    use Region_mod, only : t_Region
    use AcousticNoise_mod, only : t_AcousticNoise

    ! <<< Internal modules >>>
    use ErrorHandler, only : gracefulExit, writeAndFlush
    use PLOT3DHelper, only : plot3dErrorMessage, plot3dGetOffset, plot3dReadSingleFunction

    ! <<< Arguments >>>
    class(t_AcousticNoise) :: this
    character(len = *), intent(in) :: filename
    class(t_Region) :: region

    ! <<< Local variables >>>
    character(len = STRING_LENGTH) :: message
    logical :: success
    integer :: i, j, errorRank, procRank, ierror
    integer(kind = MPI_OFFSET_KIND) :: offset

    assert(len_trim(filename) > 0)

    write(message, '(3A)') "Reading '", trim(filename), "'..."
    call writeAndFlush(region%comm, output_unit, message, advance = 'no')

    do i = 1, size(region%gridCommunicators)

       success = .true.

       do j = 1, size(region%grids)
          if (region%grids(j)%index == i) then !... read one grid at a time
             offset = plot3dGetOffset(region%gridCommunicators(i), filename, i, success)
             if (.not. success) exit
             call plot3dReadSingleFunction(region%grids(j)%comm, trim(filename), offset,     &
                  region%grids(j)%mpiDerivedTypeScalarSubarray, region%grids(j)%globalSize,  &
                  this%data_(j)%meanPressure, success)
             exit
          end if
       end do

       call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL,                             &
            MPI_LAND, region%comm, ierror)
       if (.not. success) exit
       call MPI_Barrier(region%comm, ierror)

    end do

    if (success) then
       write(message, '(A)') " done!"
    else
       write(message, '(A)') " failed!"
    end if
    call writeAndFlush(region%comm, output_unit, message)

    if (.not. success) then
       call MPI_Comm_rank(region%comm, procRank, ierror)
       errorRank = 0
       if (len_trim(plot3dErrorMessage) > 0) errorRank = procRank
       call MPI_Allreduce(MPI_IN_PLACE, errorRank, 1, MPI_INTEGER,                           &
            MPI_MAX, region%comm, ierror)
       call MPI_Bcast(plot3dErrorMessage, STRING_LENGTH, MPI_CHARACTER,                      &  
            errorRank, region%comm, ierror)
       call gracefulExit(region%comm, plot3dErrorMessage)
    end if

  end subroutine loadMeanPressure

end module AcousticNoiseImpl

subroutine setupAcousticNoise(this, region)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use AcousticNoise_mod, only : t_AcousticNoise

  ! <<< Enumerations >>>
  use State_enum, only : QOI_DUMMY_FUNCTION

  ! <<< Private members >>>
  use AcousticNoiseImpl, only : loadMeanPressure

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables
  use InputHelper, only : getOption, getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_AcousticNoise) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer :: i
  character(len = STRING_LENGTH) :: filename

  assert(allocated(region%states))
  assert(size(region%states) > 0)

  call this%cleanup()

  call this%setupBase(region%simulationFlags, region%solverOptions)

  allocate(this%data_(size(region%states)))
  do i = 1, size(this%data_)
     allocate(this%data_(i)%meanPressure(region%grids(i)%nGridPoints, 1))
     region%states(i)%dummyFunction => this%data_(i)%meanPressure
  end do
  
  if (.not. region%simulationFlags%useTargetState) then
     call getRequiredOption("mean_pressure_file", filename)
     call loadRegionData(QOI_DUMMY_FUNCTION, filename)
  else
     filename = getOption("mean_pressure_file", "")
     if (len_trim(filename) > 0) then
        call loadRegionData(QOI_DUMMY_FUNCTION, filename)
     else
        do i = 1, size(this%data_)
           assert(allocated(region%states(i)%targetState))
           call computeDependentVariables(region%grids(i)%nDimensions,                       &
                region%states(i)%targetState, region%solverOptions%ratioOfSpecificHeats,     &
                pressure = this%data_(i)%meanPressure(:,1))
        end do
     end if
  end if

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

  ! <<< Arguments >>>
  class(t_AcousticNoise) :: this
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, ierror
  SCALAR_TYPE, allocatable :: F(:,:)

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
     assert(associated(this%data_(i)%meanPressure))
     assert(size(this%data_(i)%meanPressure, 1) == region%grids(i)%nGridPoints)
     assert(size(this%data_(i)%meanPressure, 2) == 1)

     allocate(F(region%grids(i)%nGridPoints, 1))
     F = region%states(i)%pressure - this%data_(i)%meanPressure
     instantaneousFunctional = instantaneousFunctional +                                     &
          region%grids(i)%computeInnerProduct(F, F, region%grids(i)%targetMollifier(:,1))
     SAFE_DEALLOCATE(F)

  end do

  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE, instantaneousFunctional, 1,                          &
       SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(instantaneousFunctional, 1, SCALAR_TYPE_MPI,                             &
          0, region%grids(i)%comm, ierror)
  end do

end function computeAcousticNoise

subroutine addAcousticNoiseAdjointForcing(this, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use AcousticNoise_mod, only : t_AcousticNoise
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Arguments >>>
  class(t_AcousticNoise) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

end subroutine addAcousticNoiseAdjointForcing
