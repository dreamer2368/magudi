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
    use PatchDescriptor_type, only : ACTUATOR

    ! <<< Internal modules >>>
    use ErrorHandler, only : gracefulExit

    ! <<< Arguments >>>
    type(t_Region) :: region

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: mollifierNorm
    integer :: i, ierror
    logical :: hasNegativeMollifier
    character(len = STRING_LENGTH) :: str

    assert(allocated(region%grids))

    mollifierNorm = 0.0_wp

    do i = 1, size(region%grids)

       assert(allocated(region%grids(i)%controlMollifier))

       hasNegativeMollifier = any(real(region%grids(i)%controlMollifier(:,1), wp) < 0.0_wp)
       call MPI_Allreduce(MPI_IN_PLACE, hasNegativeMollifier, 1,                             &
            MPI_LOGICAL, MPI_LOR, region%grids(i)%comm, ierror)
       if (hasNegativeMollifier) then
          write(str, '(A,I0.0,A)') "Control mollifying support function on grid ",           &
               region%grids(i)%index, " is not non-negative everywhere!"
          call gracefulExit(region%grids(i)%comm, str)
       end if
       mollifierNorm = mollifierNorm +                                                       &
            real(region%grids(i)%computeQuadratureOnPatches(                                 &
       region%grids(i)%controlMollifier(:,1), region%patches, ACTUATOR), wp)
    end do

    if (region%commGridMasters /= MPI_COMM_NULL)                                             &
         call MPI_Allreduce(MPI_IN_PLACE, mollifierNorm, 1, REAL_TYPE_MPI,                   &
         MPI_SUM, region%commGridMasters, ierror)

    do i = 1, size(region%grids)
       call MPI_Bcast(mollifierNorm, 1, REAL_TYPE_MPI, 0, region%grids(i)%comm, ierror)
    end do
    if (mollifierNorm <= 0.0_wp)                                                             &
         call gracefulExit(region%comm, "Control mollifying support is trivial!")

    do i = 1, size(region%grids)
       region%grids(i)%controlMollifier = region%grids(i)%controlMollifier / mollifierNorm
    end do

  end subroutine normalizeControlMollifier

  subroutine normalizeTargetMollifier(region)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_type, only : t_Region
    use PatchDescriptor_type, only : CONTROL_TARGET

    ! <<< Internal modules >>>
    use ErrorHandler, only : gracefulExit

    ! <<< Arguments >>>
    type(t_Region) :: region

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: mollifierNorm
    integer :: i, ierror
    logical :: hasNegativeMollifier
    character(len = STRING_LENGTH) :: str

    assert(allocated(region%grids))

    mollifierNorm = 0.0_wp

    do i = 1, size(region%grids)

       assert(allocated(region%grids(i)%targetMollifier))

       hasNegativeMollifier = any(real(region%grids(i)%targetMollifier(:,1), wp) < 0.0_wp)
       call MPI_Allreduce(MPI_IN_PLACE, hasNegativeMollifier, 1,                             &
            MPI_LOGICAL, MPI_LOR, region%grids(i)%comm, ierror)
       if (hasNegativeMollifier) then
          write(str, '(A,I0.0,A)') "Target mollifying support function on grid ",           &
               region%grids(i)%index, " is not non-negative everywhere!"
          call gracefulExit(region%grids(i)%comm, str)
       end if
       mollifierNorm = mollifierNorm +                                                       &
            real(region%grids(i)%computeQuadratureOnPatches(                                 &
       region%grids(i)%targetMollifier(:,1), region%patches, CONTROL_TARGET), wp)
    end do

    if (region%commGridMasters /= MPI_COMM_NULL)                                             &
         call MPI_Allreduce(MPI_IN_PLACE, mollifierNorm, 1, REAL_TYPE_MPI,                   &
         MPI_SUM, region%commGridMasters, ierror)

    do i = 1, size(region%grids)
       call MPI_Bcast(mollifierNorm, 1, REAL_TYPE_MPI, 0, region%grids(i)%comm, ierror)
    end do
    if (mollifierNorm <= 0.0_wp)                                                             &
         call gracefulExit(region%comm, "Target mollifying support is trivial!")

    do i = 1, size(region%grids)
       region%grids(i)%targetMollifier = region%grids(i)%targetMollifier / mollifierNorm
    end do

  end subroutine normalizeTargetMollifier

  function instantaneousCostFunctional(region)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_type, only : t_Region
    use SolverOptions_type

    ! <<< Arguments >>>
    type(t_Region) :: region

    ! <<< Result >>>
    SCALAR_TYPE :: instantaneousCostFunctional

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, ierror
    SCALAR_TYPE, allocatable :: F(:,:)

    assert(allocated(region%grids))
    assert(allocated(region%states))
    assert(size(region%grids) == size(region%states))

    assert_key(region%solverOptions%costFunctionalType, ( \
    SOUND, \
    LIFT,  \
    DRAG))

    instantaneousCostFunctional = 0.0_wp

    do i = 1, size(region%grids)

       assert(region%grids(i)%nGridPoints > 0)
       assert(allocated(region%grids(i)%targetMollifier))
       assert(size(region%grids(i)%targetMollifier, 1) == region%grids(i)%nGridPoints)
       assert(size(region%grids(i)%targetMollifier, 2) == 1)

       select case (region%solverOptions%costFunctionalType)

       case (SOUND)

          assert(allocated(region%states(i)%pressure))
          assert(size(region%states(i)%pressure, 1) == region%grids(i)%nGridPoints)
          assert(size(region%states(i)%pressure, 2) == 1)

          assert(allocated(region%states(i)%meanPressure))
          assert(size(region%states(i)%meanPressure, 1) == region%grids(i)%nGridPoints)
          assert(size(region%states(i)%meanPressure, 2) == 1)

          allocate(F(region%grids(i)%nGridPoints, 1))
          F = region%states(i)%pressure - region%states(i)%meanPressure
          instantaneousCostFunctional = instantaneousCostFunctional +                        &
               region%grids(i)%computeInnerProduct(F, F, region%grids(i)%targetMollifier(:,1))
          SAFE_DEALLOCATE(F)

       end select

    end do

    if (region%commGridMasters /= MPI_COMM_NULL)                                             &
         call MPI_Allreduce(MPI_IN_PLACE, instantaneousCostFunctional, 1,                    &
         SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

    do i = 1, size(region%grids)
       call MPI_Bcast(instantaneousCostFunctional, 1, SCALAR_TYPE_MPI,                       &
            0, region%grids(i)%comm, ierror)
    end do

  end function instantaneousCostFunctional

  subroutine writeLine(comm, filename, line)

    ! <<< External modules >>>
    use MPI

    ! <<< Internal modules >>>
    use ErrorHandler, only : gracefulExit

    ! <<< Arguments >>>
    integer, intent(in) :: comm
    character(len = *), intent(in) :: filename, line

    ! <<< Local variables >>>
    integer :: fileUnit, ostat, procRank, ierror
    logical, save :: firstCall = .true.
    character(len = STRING_LENGTH) :: message

    call MPI_Comm_rank(comm, procRank, ierror)

    if (procRank == 0) then

       if (firstCall) then
          open(newunit = fileUnit, file = trim(filename), action = 'write',                  &
               status = 'unknown', iostat = ostat)
       else
          open(newunit = fileUnit, file = trim(filename), action = 'write',                  &
               status = 'old', position = 'append', iostat = ostat)
       end if

    end if

    call MPI_Bcast(ostat, 1, MPI_INTEGER, 0, comm, ierror)
    if (ostat /= 0) then
       write(message, "(2A)") trim(filename), ": Failed to open file for writing!"
       call gracefulExit(comm, message)
    end if

    firstCall = .false.

    if (procRank == 0) then
       write(fileUnit, '(A)', iostat = ostat) trim(line)
    end if

    call MPI_Bcast(ostat, 1, MPI_INTEGER, 0, comm, ierror)
    if (ostat /= 0) then
       write(message, "(2A)") trim(filename), ": Error writing to file!"
       call gracefulExit(comm, message)
    end if

    if (procRank == 0) then
       flush(fileUnit)
       close(fileUnit)
    end if

  end subroutine writeLine

end module SolverImpl

subroutine initializeSolver(region, restartFilename)

  ! <<< Derived types >>>
  use State_type, only : QOI_FORWARD_STATE, QOI_TARGET_STATE, QOI_MEAN_PRESSURE
  use Region_type, only : t_Region
  use SolverOptions_type, only : SOUND

  ! <<< Enumerations >>>
  use Grid_enum

  ! <<< Private members >>>
  use SolverImpl, only : normalizeControlMollifier, normalizeTargetMollifier

  ! <<< Internal modules >>>
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
              region%states(i)%plot3dAuxiliaryData = 0.0_wp
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
     case (SOUND)

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
     saveInterval, reportInterval, outputPrefix, costFunctional)

  ! <<< External modules >>>
  use iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use State_type, only : t_State, QOI_FORWARD_STATE
  use Region_type, only : t_Region
  use TimeIntegrator_mod, only : t_TimeIntegrator

  ! <<< Private members >>>
  use SolverImpl, only : instantaneousCostFunctional, writeLine

  ! <<< Internal modules >>>
  use State_mod, only : updateState
  use Region_mod, only : saveRegionData, getTimeStepSize, getCfl, reportResiduals
  use ErrorHandler, only : writeAndFlush
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  type(t_Region) :: region
  class(t_TimeIntegrator) :: integrator
  real(SCALAR_KIND), intent(inout) :: time
  integer, intent(inout) :: timestep
  integer, intent(in) :: nTimesteps
  integer, intent(in), optional :: saveInterval, reportInterval
  character(len = STRING_LENGTH), intent(in), optional :: outputPrefix
  SCALAR_TYPE, intent(out), optional :: costFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: outputPrefix_, filename, str
  integer :: i, j, timestep_
  logical :: verbose
  real(wp) :: timeStepSize, cfl
  SCALAR_TYPE :: instantaneousCostFunctional_

  assert(timestep >= 0)

  call startTiming("solveForward")

  outputPrefix_ = PROJECT_NAME
  if (present(outputPrefix)) outputPrefix_ = outputPrefix

  if (present(costFunctional)) costFunctional = 0.0_wp

  write(filename, '(2A,I8.8,A)') trim(outputPrefix_), "-", timestep, ".q"
  call saveRegionData(region, QOI_FORWARD_STATE, filename)

  verbose = .false.

  do timestep_ = timestep + 1, timestep + nTimesteps

     if (present(reportInterval)) verbose = (reportInterval > 0 .and.                        &
          mod(timestep_, reportInterval) == 0)

     do j = 1, size(region%states) !... update state
        call updateState(region%states(j), region%grids(j),                                  &
             region%simulationFlags, region%solverOptions)
     end do
     timeStepSize = getTimeStepSize(region)

     do i = 1, integrator%nStages

        call integrator%substepForward(region, time, timeStepSize, timestep_, i)

        if (i /= integrator%nStages) then
           do j = 1, size(region%states) !... update state
              call updateState(region%states(j), region%grids(j),                            &
                   region%simulationFlags, region%solverOptions)
           end do
        end if

        if (.not. region%simulationFlags%predictionOnly .and. present(costFunctional)) then
           instantaneousCostFunctional_ = instantaneousCostFunctional(region)
           costFunctional = costFunctional +                                                 &
                integrator%norm(i) * timeStepSize * instantaneousCostFunctional_
        end if

     end do

     if (verbose) then

        if (region%simulationFlags%useConstantCfl) then
           write(str, '(2A,I8,2(A,E13.6))') PROJECT_NAME, ": timestep = ", timestep_,        &
                ", dt = ", timeStepSize, ", time = ", time
        else
           cfl = getCfl(region)
           write(str, '(2A,I8,2(A,E13.6))') PROJECT_NAME, ": timestep = ", timestep_,        &
                ", CFL = ", cfl, ", time = ", time
        end if

        call writeAndFlush(region%comm, output_unit, str)

        if (region%simulationFlags%steadyStateSimulation) call reportResiduals(region)

        if (.not. region%simulationFlags%predictionOnly .and. present(costFunctional)) then
           write(filename, '(2A)') trim(outputPrefix_), ".cost_functional.txt"
           write(str, '(I8,1X,E13.6,1X,SP,' // SCALAR_FORMAT // ')')                      &
                timestep_, time, instantaneousCostFunctional_
           call writeLine(region%comm, filename, str)
        end if

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

  call endTiming("solveForward")

end subroutine solveForward

subroutine solveAdjoint(region, integrator, time, timestep, nTimesteps,                      &
     saveInterval, reportInterval, outputPrefix)

  ! <<< External modules >>>
  use iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use State_type, only : t_State, QOI_ADJOINT_STATE
  use Region_type, only : t_Region
  use TimeIntegrator_mod, only : t_TimeIntegrator
  use ReverseMigrator_type, only : t_ReverseMigrator

  ! <<< Internal modules >>>
  use State_mod, only : updateState
  use Region_mod, only : saveRegionData, getTimeStepSize, getCfl, reportResiduals
  use InputHelper, only : getOption
  use ErrorHandler, only : writeAndFlush
  use MPITimingsHelper, only : startTiming, endTiming
  use ReverseMigrator_mod

  implicit none

  ! <<< Arguments >>>
  type(t_Region) :: region
  class(t_TimeIntegrator) :: integrator
  real(SCALAR_KIND), intent(inout) :: time
  integer, intent(inout) :: timestep
  integer, intent(in) :: nTimesteps, saveInterval
  integer, intent(in), optional :: reportInterval
  character(len = STRING_LENGTH), intent(in), optional :: outputPrefix

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: outputPrefix_, filename, str
  type(t_ReverseMigrator) :: reverseMigrator
  integer :: i, j, timestep_
  logical :: verbose
  real(wp) :: timeStepSize, cfl

  assert(timestep >= nTimesteps)

  call startTiming("solveAdjoint")

  outputPrefix_ = PROJECT_NAME
  if (present(outputPrefix)) outputPrefix_ = outputPrefix

  call setupReverseMigrator(reverseMigrator, region, outputPrefix_,                          &
       getOption("checkpointing_scheme", "uniform checkpointing"),                           &
       timestep - nTimesteps, timestep,                                                      &
       saveInterval, saveInterval * integrator%nStages)

  write(filename, '(2A,I8.8,A)') trim(outputPrefix_), "-", timestep, ".adjoint.q"
  call saveRegionData(region, QOI_ADJOINT_STATE, filename)

  verbose = .false.

  do timestep_ = timestep - 1, timestep - nTimesteps, -1

     if (present(reportInterval)) verbose = (reportInterval > 0 .and.                        &
          mod(timestep_, reportInterval) == 0)

     if (region%simulationFlags%useConstantCfl) then
        call migrateToSubstep(reverseMigrator, region, integrator, timestep_, 1)
        do j = 1, size(region%states) !... update state
           call updateState(region%states(j), region%grids(j),                               &
                region%simulationFlags, region%solverOptions)
        end do
     end if
     timeStepSize = getTimeStepSize(region)

     do i = integrator%nStages, 1, -1

        if (.not. region%simulationFlags%steadyStateSimulation) then
           if (i == 1) then
              call migrateToSubstep(reverseMigrator, region,                                 &
                   integrator, timestep_, integrator%nStages)
           else
              call migrateToSubstep(reverseMigrator, region,                                 &
                   integrator, timestep_ + 1, i - 1)
           end if
        end if

        do j = 1, size(region%states) !... update state
           call updateState(region%states(j), region%grids(j),                               &
                region%simulationFlags, region%solverOptions)
        end do

        call integrator%substepAdjoint(region, time, timeStepSize, timestep_, i)

     end do

     if (verbose) then
        if (region%simulationFlags%useConstantCfl) then
           write(str, '(2A,I8,2(A,E13.6))') PROJECT_NAME, ": timestep = ", timestep_,        &
                ", dt = ", timeStepSize, ", time = ", time
        else
           cfl = getCfl(region)
           write(str, '(2A,I8,2(A,E13.6))') PROJECT_NAME, ": timestep = ", timestep_,        &
                ", CFL = ", cfl, ", time = ", time
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

  call endTiming("solveAdjoint")

end subroutine solveAdjoint
