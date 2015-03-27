#include "config.h"

module SolverImpl

  implicit none
  public

contains

  subroutine normalizeControlMollifier(region)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_type, only : t_Region

    ! <<< Internal modules >>>
    use Grid_mod, only : computeInnerProduct, isVariableWithinRange
    use ErrorHandler, only : gracefulExit

    ! <<< Arguments >>>
    type(t_Region) :: region

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    SCALAR_TYPE :: temp
    real(wp) :: mollifierNorm
    integer :: i, j, k, procRank, ierror
    logical :: hasNegativeMollifier
    character(len = STRING_LENGTH) :: str
    SCALAR_TYPE, allocatable :: F(:,:)

    mollifierNorm = 0.0_wp

    do i = 1, size(region%grids)

       hasNegativeMollifier = .not. isVariableWithinRange(region%grids(i),                   &
            region%grids(i)%controlMollifier(:,1), temp, i, j, k,                            &
            minValue = - 2.0_wp * epsilon(0.0_wp))
       if (hasNegativeMollifier) then
          write(str, '(A,4(I0.0,A),3(ES11.4,A))') "Control mollifier on grid ",              &
               region%grids(i)%index, " at (", i, ", ", j, ", ", k, "): ",                   &
               temp, " is negative!"
          call gracefulExit(region%grids(i)%comm, str)
       end if

       allocate(F(region%grids(i)%nGridPoints, 1), source = 1.0_wp)
       temp = computeInnerProduct(region%grids(i), F, F,                                    &
            region%grids(i)%controlMollifier(:,1))

       call MPI_Comm_rank(region%grids(i)%comm, procRank, ierror)
       if (procRank /= 0) temp = 0.0_wp

       mollifierNorm = mollifierNorm + real(temp, wp)

       SAFE_DEALLOCATE(F)

    end do

#ifdef SCALAR_IS_COMPLEX
    call MPI_Allreduce(MPI_IN_PLACE, mollifierNorm, 1, REAL_TYPE_MPI,                       &
         MPI_SUM, region%comm, ierror)
#else
    call MPI_Allreduce(MPI_IN_PLACE, mollifierNorm, 1, SCALAR_TYPE_MPI,                     &
         MPI_SUM, region%comm, ierror)
#endif

    do i = 1, size(region%grids)
       region%grids(i)%controlMollifier = region%grids(i)%controlMollifier / mollifierNorm
    end do

  end subroutine normalizeControlMollifier

  subroutine normalizeTargetMollifier(region)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_type, only : t_Region

    ! <<< Internal modules >>>
    use Grid_mod, only : computeInnerProduct, isVariableWithinRange
    use ErrorHandler, only : gracefulExit

    ! <<< Arguments >>>
    type(t_Region) :: region

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    SCALAR_TYPE :: temp
    real(wp) :: mollifierNorm
    integer :: i, j, k, procRank, ierror
    logical :: hasNegativeMollifier
    character(len = STRING_LENGTH) :: str
    SCALAR_TYPE, allocatable :: F(:,:)

    mollifierNorm = 0.0_wp

    do i = 1, size(region%grids)

       hasNegativeMollifier = .not. isVariableWithinRange(region%grids(i),                   &
            region%grids(i)%targetMollifier(:,1), temp, i, j, k,                            &
            minValue = - 2.0_wp * epsilon(0.0_wp))
       if (hasNegativeMollifier) then
          write(str, '(A,4(I0.0,A),3(ES11.4,A))') "Target mollifier on grid ",              &
               region%grids(i)%index, " at (", i, ", ", j, ", ", k, "): ",                   &
               temp, " is negative!"
          call gracefulExit(region%grids(i)%comm, str)
       end if

       allocate(F(region%grids(i)%nGridPoints, 1), source = 1.0_wp)
       temp = computeInnerProduct(region%grids(i), F, F,                                    &
            region%grids(i)%targetMollifier(:,1))

       call MPI_Comm_rank(region%grids(i)%comm, procRank, ierror)
       if (procRank /= 0) temp = 0.0_wp

       mollifierNorm = mollifierNorm + real(temp, wp)

       SAFE_DEALLOCATE(F)

    end do

#ifdef SCALAR_IS_COMPLEX
    call MPI_Allreduce(MPI_IN_PLACE, mollifierNorm, 1, REAL_TYPE_MPI,                       &
         MPI_SUM, region%comm, ierror)
#else
    call MPI_Allreduce(MPI_IN_PLACE, mollifierNorm, 1, SCALAR_TYPE_MPI,                     &
         MPI_SUM, region%comm, ierror)
#endif

    do i = 1, size(region%grids)
       region%grids(i)%targetMollifier = region%grids(i)%targetMollifier / mollifierNorm
    end do

  end subroutine normalizeTargetMollifier

end module SolverImpl

subroutine initializeSolver(region, restartFilename)

  ! <<< Derived types >>>
  use Grid_type, only : QOI_CONTROL_MOLLIFIER, QOI_TARGET_MOLLIFIER
  use State_type, only : QOI_FORWARD_STATE, QOI_TARGET_STATE, QOI_MEAN_PRESSURE
  use Region_type, only : t_Region
  use SolverOptions_type, only : SOUND_FUNCTIONAL

  ! <<< Private members >>>
  use SolverImpl, only : normalizeControlMollifier, normalizeTargetMollifier

  ! <<< Internal modules >>>
  use Grid_mod, only : computeInnerProduct
  use CNSHelper, only : computeDependentVariables
  use State_mod, only : makeQuiescent
  use Region_mod, only : loadRegionData
  use InputHelper, only : getOption, getRequiredOption

  implicit none

  ! <<< Arguments >>>
  type(t_Region) :: region
  character(len = *), intent(in), optional :: restartFilename

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: filename
  integer :: i

  ! If a target state file was specified, read the target state. Otherwise, initialize the
  ! target state to a quiescent state by default.
  if (region%simulationFlags%useTargetState) then
     filename = getOption("target_state_file", "")
     if (len_trim(filename) == 0) then
        do i = 1, size(region%states)
           call makeQuiescent(region%states(i), size(region%globalGridSizes, 1),             &
                region%solverOptions%ratioOfSpecificHeats, region%states(i)%targetState)
        end do
     else
        call loadRegionData(region, QOI_TARGET_STATE, filename)
     end if
  end if

  ! Initialize conserved variables.
  if (present(restartFilename) .and. region%simulationFlags%predictionOnly) then
     call loadRegionData(region, QOI_FORWARD_STATE, restartFilename)
  else if (region%simulationFlags%predictionOnly .or. .not. &
       region%simulationFlags%isBaselineAvailable) then
     if (region%simulationFlags%useTargetState) then
        filename = getOption("initial_condition_file", "")
        if (len_trim(filename) == 0) then
           do i = 1, size(region%states) !... initialize from target state.
              region%states(i)%conservedVariables = region%states(i)%targetState
           end do
        else
           call loadRegionData(region, QOI_FORWARD_STATE, filename) !... initialize from file.
        end if
     else
        call getRequiredOption("initial_condition_file", filename)
        call loadRegionData(region, QOI_FORWARD_STATE, filename) !... initialize from file.
     end if
  end if

  if (.not. region%simulationFlags%predictionOnly) then

     ! Initialize adjoint variables.
     do i = 1, size(region%states)
        region%states(i)%adjointVariables = 0.0_wp
     end do

     ! Initialize control mollifier.
     filename = getOption("control_mollifier_file", "")
     if (len_trim(filename) == 0) then
        do i = 1, size(region%grids)
           region%grids(i)%controlMollifier = 1.0_wp
        end do
     else
        call loadRegionData(region, QOI_CONTROL_MOLLIFIER, filename)
     end if
     call normalizeControlMollifier(region)

     ! Target mollifier.
     filename = getOption("target_mollifier_file", "")
     if (len_trim(filename) == 0) then
        do i = 1, size(region%grids)
           region%grids(i)%targetMollifier = 1.0_wp
        end do
     else
        call loadRegionData(region, QOI_TARGET_MOLLIFIER, filename)
     end if
     call normalizeTargetMollifier(region)

     select case (region%solverOptions%costFunctionalType)
     case (SOUND_FUNCTIONAL)

        ! Mean pressure.
        if (.not. region%simulationFlags%useTargetState) then
           call getRequiredOption("mean_pressure_file", filename)
           call loadRegionData(region, QOI_MEAN_PRESSURE, filename)
        else
           filename = getOption("mean_pressure_file", "")
           if (len_trim(filename) == 0) then
              do i = 1, size(region%grids)
                 call computeDependentVariables(size(region%globalGridSizes, 1),             &
                      region%states(i)%targetState,                                          &
                      region%solverOptions%ratioOfSpecificHeats,                             &
                      pressure = region%states(i)%meanPressure(:,1))
              end do
           else
              call loadRegionData(region, QOI_MEAN_PRESSURE, filename)
           end if
        end if

     end select

  end if

end subroutine initializeSolver

subroutine solveForward(region, integrator, time, timestep, nTimesteps,                      &
     saveInterval, reportInterval, outputPrefix)

  ! <<< External modules >>>
  use iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use State_type, only : t_State, QOI_FORWARD_STATE
  use Region_type, only : t_Region
  use RK4Integrator_type, only : t_RK4Integrator

  ! <<< Internal modules >>>
  use Region_mod, only : saveRegionData, reportResiduals
  use ErrorHandler, only : writeAndFlush
  use RK4Integrator_mod, only : substepForward

  implicit none

  ! <<< Arguments >>>
  type(t_Region) :: region
  type(t_RK4Integrator) :: integrator
  real(SCALAR_KIND), intent(inout) :: time
  integer, intent(inout) :: timestep
  integer, intent(in) :: nTimesteps
  integer, intent(in), optional :: saveInterval, reportInterval
  character(len = STRING_LENGTH), intent(in), optional :: outputPrefix

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: outputPrefix_, filename, str
  integer :: i, timestep_
  logical :: verbose

  assert(timestep >= 0)

  outputPrefix_ = PROJECT_NAME
  if (present(outputPrefix)) outputPrefix_ = outputPrefix

  write(filename, '(2A,I8.8,A)') trim(outputPrefix_), "-", timestep, ".q"
  call saveRegionData(region, QOI_FORWARD_STATE, filename)

  verbose = .false.

  do timestep_ = timestep + 1, timestep + nTimesteps

     if (present(reportInterval)) verbose = (reportInterval > 0 .and.                        &
          mod(timestep_, reportInterval) == 0)

     do i = 1, 4
        call substepForward(integrator, region, time, timestep_, i)
     end do

     if (verbose) then
        if (region%simulationFlags%useConstantCfl) then
           write(str, '(2A,I8,2(A,D13.6))') PROJECT_NAME, ": timestep = ", timestep_,        &
                ", dt = ", region%states(1)%timeStepSize, ", time = ", time
        else
           write(str, '(2A,I8,2(A,D13.6))') PROJECT_NAME, ": timestep = ", timestep_,        &
                ", CFL = ", region%states(1)%cfl, ", time = ", time
        end if
        call writeAndFlush(region%comm, output_unit, str)
        if (region%simulationFlags%steadyStateSimulation) call reportResiduals(region)
     end if

     if (present(saveInterval)) then
        if (saveInterval > 0) then
           if (mod(timestep_, saveInterval) == 0) then
              do i = 1, size(region%states)
                 region%states(i)%plot3dAuxiliaryData(1) = real(timestep_, wp)
                 region%states(i)%plot3dAuxiliaryData(4) = time
              end do
              write(filename, '(2A,I8.8,A)') trim(outputPrefix_), "-", timestep_, ".q"
              call saveRegionData(region, QOI_FORWARD_STATE, filename)
           end if
        end if
     end if

  end do

  timestep = timestep + nTimesteps

end subroutine solveForward

subroutine solveAdjoint(region, integrator, time, timestep, nTimesteps,                      &
     saveInterval, reportInterval, outputPrefix)

  ! <<< External modules >>>
  use iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use State_type, only : t_State, QOI_ADJOINT_STATE
  use Region_type, only : t_Region
  use RK4Integrator_type, only : t_RK4Integrator
  use ReverseMigrator_type, only : t_ReverseMigrator

  ! <<< Internal modules >>>
  use Region_mod, only : saveRegionData, reportResiduals
  use InputHelper, only : getOption
  use ErrorHandler, only : writeAndFlush
  use RK4Integrator_mod, only : substepAdjoint
  use ReverseMigrator_mod

  implicit none

  ! <<< Arguments >>>
  type(t_Region) :: region
  type(t_RK4Integrator) :: integrator
  real(SCALAR_KIND), intent(inout) :: time
  integer, intent(inout) :: timestep
  integer, intent(in) :: nTimesteps, saveInterval
  integer, intent(in), optional :: reportInterval
  character(len = STRING_LENGTH), intent(in), optional :: outputPrefix

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: outputPrefix_, filename, str
  type(t_ReverseMigrator) :: reverseMigrator
  integer :: i, timestep_
  logical :: verbose

  assert(timestep >= nTimesteps)

  outputPrefix_ = PROJECT_NAME
  if (present(outputPrefix)) outputPrefix_ = outputPrefix

  call setupReverseMigrator(reverseMigrator, region, outputPrefix_,                          &
       getOption("checkpointing_scheme", "uniform checkpointing"),                           &
       timestep - nTimesteps, timestep,                                                      &
       saveInterval, saveInterval * 4)

  write(filename, '(2A,I8.8,A)') trim(outputPrefix_), "-", timestep, ".adjoint.q"
  call saveRegionData(region, QOI_ADJOINT_STATE, filename)

  verbose = .false.

  do timestep_ = timestep - 1, timestep - nTimesteps, -1

     if (present(reportInterval)) verbose = (reportInterval > 0 .and.                        &
          mod(timestep_, reportInterval) == 0)

     do i = 4, 1, -1
        if (.not. region%simulationFlags%steadyStateSimulation) then
           if (i == 1) then
              call migrateToSubstep(reverseMigrator, region,                                 &
                   integrator, timestep_, 4)
           else
              call migrateToSubstep(reverseMigrator, region,                                 &
                   integrator, timestep_ + 1, i - 1)
           end if
        end if
        call substepAdjoint(integrator, region, time, timestep_, i)
     end do

     if (verbose) then
        if (region%simulationFlags%useConstantCfl) then
           write(str, '(2A,I8,2(A,D13.6))') PROJECT_NAME, ": timestep = ", timestep_,        &
                ", dt = ", region%states(1)%timeStepSize, ", time = ", time
        else
           write(str, '(2A,I8,2(A,D13.6))') PROJECT_NAME, ": timestep = ", timestep_,        &
                ", CFL = ", region%states(1)%cfl, ", time = ", time
        end if
        call writeAndFlush(region%comm, output_unit, str)
        if (region%simulationFlags%steadyStateSimulation) call reportResiduals(region)
     end if

     if (saveInterval > 0) then
        if (mod(timestep_, saveInterval) == 0) then
           do i = 1, size(region%states)
              region%states(i)%plot3dAuxiliaryData(1) = real(timestep_, wp)
              region%states(i)%plot3dAuxiliaryData(4) = time
           end do
           write(filename, '(2A,I8.8,A)')                                                    &
                trim(outputPrefix_), "-", timestep_, ".adjoint.q"
           call saveRegionData(region, QOI_ADJOINT_STATE, filename)
        end if
     end if

  end do

  call cleanupReverseMigrator(reverseMigrator)

  timestep = timestep - nTimesteps

end subroutine solveAdjoint
