#include "config.h"

module SolverImpl

  implicit none
  public

contains

  subroutine normalizeControlMollifier(region)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Internal modules >>>
    use ErrorHandler, only : gracefulExit, issueWarning
    use Patch_factory, only : computeQuadratureOnPatches

    ! <<< Arguments >>>
    class(t_Region) :: region

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
            real(computeQuadratureOnPatches(region%patchFactories,                           &
            'ACTUATOR', region%grids(i), region%grids(i)%controlMollifier(:,1)), wp)
    end do

    if (region%commGridMasters /= MPI_COMM_NULL)                                             &
         call MPI_Allreduce(MPI_IN_PLACE, mollifierNorm, 1, REAL_TYPE_MPI,                   &
         MPI_SUM, region%commGridMasters, ierror)

    do i = 1, size(region%grids)
       call MPI_Bcast(mollifierNorm, 1, REAL_TYPE_MPI, 0, region%grids(i)%comm, ierror)
    end do
    if (mollifierNorm <= 0.0_wp)                                                             &
         call issueWarning(region%comm,                                                      &
         "Control mollifying support is trivial! Is an actuator patch present?")

    do i = 1, size(region%grids)
       region%grids(i)%controlMollifier = region%grids(i)%controlMollifier / mollifierNorm
    end do

  end subroutine normalizeControlMollifier

  subroutine normalizeTargetMollifier(region)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Internal modules >>>
    use ErrorHandler, only : gracefulExit, issueWarning
    use Patch_factory, only : computeQuadratureOnPatches

    ! <<< Arguments >>>
    class(t_Region) :: region

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
            real(computeQuadratureOnPatches(region%patchFactories,                           &
            'COST_TARGET', region%grids(i), region%grids(i)%targetMollifier(:,1)), wp)
    end do

    if (region%commGridMasters /= MPI_COMM_NULL)                                             &
         call MPI_Allreduce(MPI_IN_PLACE, mollifierNorm, 1, REAL_TYPE_MPI,                   &
         MPI_SUM, region%commGridMasters, ierror)

    do i = 1, size(region%grids)
       call MPI_Bcast(mollifierNorm, 1, REAL_TYPE_MPI, 0, region%grids(i)%comm, ierror)
    end do
    if (mollifierNorm <= 0.0_wp)                                                             &
         call issueWarning(region%comm,                                                      &
         "Target mollifying support is trivial! Is a cost target patch present?")

    do i = 1, size(region%grids)
       region%grids(i)%targetMollifier = region%grids(i)%targetMollifier / mollifierNorm
    end do

  end subroutine normalizeTargetMollifier

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
  use Region_mod, only : t_Region

  ! <<< Enumerations >>>
  use Grid_enum
  use State_enum

  ! <<< Private members >>>
  use SolverImpl, only : normalizeControlMollifier, normalizeTargetMollifier

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables
  use InputHelper, only : getOption, getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region
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
           call region%states(i)%makeQuiescent(size(region%globalGridSizes, 1),              &
                region%solverOptions%ratioOfSpecificHeats, region%states(i)%targetState)
        end do
     else
        call region%loadData(QOI_TARGET_STATE, filename)
     end if
  end if

  ! Initialize conserved variables.
  if (present(restartFilename) .and. region%simulationFlags%predictionOnly) then
     call region%loadData(QOI_FORWARD_STATE, restartFilename)
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
           call region%loadData(QOI_FORWARD_STATE, filename) !... initialize from file.
        end if
     else
        call getRequiredOption("initial_condition_file", filename)
        call region%loadData(QOI_FORWARD_STATE, filename) !... initialize from file.
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
        call region%loadData(QOI_CONTROL_MOLLIFIER, filename)
     end if
     call normalizeControlMollifier(region)

     ! Target mollifier.
     filename = getOption("target_mollifier_file", "")
     if (len_trim(filename) == 0) then
        do i = 1, size(region%grids)
           region%grids(i)%targetMollifier = 1.0_wp
        end do
     else
        call region%loadData(QOI_TARGET_MOLLIFIER, filename)
     end if
     call normalizeTargetMollifier(region)

  end if

end subroutine initializeSolver

subroutine solveForward(region, time, timestep, nTimesteps,                                  &
     saveInterval, outputPrefix, costFunctional)

  ! <<< External modules >>>
  use iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use State_mod, only : t_State
  use Region_mod, only : t_Region
  use Functional_mod, only : t_Functional
  use Functional_factory, only : t_FunctionalFactory
  use TimeIntegrator_mod, only : t_TimeIntegrator
  use TimeIntegrator_factory, only : t_TimeIntegratorFactory

  ! <<< Enumerations >>>
  use State_enum, only : QOI_FORWARD_STATE

  ! <<< Private members >>>
  use SolverImpl, only : writeLine

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region
  real(SCALAR_KIND), intent(inout) :: time
  integer, intent(inout) :: timestep
  integer, intent(in) :: nTimesteps
  integer, intent(in), optional :: saveInterval
  character(len = STRING_LENGTH), intent(in), optional :: outputPrefix
  SCALAR_TYPE, intent(out), optional :: costFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: outputPrefix_, filename, str
  integer :: i, j, timestep_, reportInterval, residualInterval
  class(t_Functional), pointer :: functional => null()
  type(t_FunctionalFactory) :: functionalFactory
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()
  type(t_TimeIntegratorFactory) :: timeIntegratorFactory
  real(wp) :: timeStepSize, cfl, residuals(3)
  SCALAR_TYPE :: instantaneousCostFunctional

  assert(timestep >= 0)

  call startTiming("solveForward")

  outputPrefix_ = PROJECT_NAME
  if (present(outputPrefix)) outputPrefix_ = outputPrefix

  if (present(costFunctional)) costFunctional = 0.0_wp

  write(filename, '(2A,I8.8,A)') trim(outputPrefix_), "-", timestep, ".q"
  call region%saveData(QOI_FORWARD_STATE, filename)

  residualInterval = getOption("check_residuals_interval", -1)
  if (residualInterval == -1) then
     call getRequiredOption("report_interval", reportInterval)
     residualInterval = reportInterval
     if (residualInterval < 0) residualInterval = max(1, nint(1e-3_wp * real(nTimesteps, wp)))
  else
     reportInterval = getOption("report_interval", 1)
  end if

  if (.not. region%simulationFlags%predictionOnly) then
     call functionalFactory%connect(functional, trim(region%solverOptions%costFunctionalType))
     assert(associated(functional))
     call functional%setup(region)
  end if

  call timeIntegratorFactory%connect(timeIntegrator,                                         &
       trim(region%solverOptions%timeIntegratorType))
  assert(associated(timeIntegrator))
  call timeIntegrator%setup(region)

  do timestep_ = timestep + 1, timestep + nTimesteps

     do j = 1, size(region%states) !... update state
        call region%states(j)%update(region%grids(j), region%simulationFlags,                &
             region%solverOptions)
     end do
     timeStepSize = region%getTimeStepSize()

     do i = 1, timeIntegrator%nStages

        call timeIntegrator%substepForward(region, time, timeStepSize, timestep_, i)

        if (i /= timeIntegrator%nStages) then
           do j = 1, size(region%states) !... update state
              call region%states(j)%update(region%grids(j), region%simulationFlags,          &
                   region%solverOptions)
           end do
        end if

        if (.not. region%simulationFlags%predictionOnly .and. present(costFunctional)) then
           instantaneousCostFunctional = functional%compute(time, region)
             costFunctional = costFunctional +                                               &
                  timeIntegrator%norm(i) * timeStepSize * instantaneousCostFunctional
          else
             instantaneousCostFunctional = 0.0_wp
          end if

     end do

     if (reportInterval > 0 .and. mod(timestep_, reportInterval) == 0) then
        if (region%simulationFlags%useConstantCfl) then
           write(str, '(2A,I8,3(A,E13.6))') PROJECT_NAME, ": timestep = ", timestep_,        &
                ", dt = ", timeStepSize, ", time = ", time,                                  &
                ", cost = ", instantaneousCostFunctional
        else
           cfl = region%getCfl()
           write(str, '(2A,I8,3(A,E13.6))') PROJECT_NAME, ": timestep = ", timestep_,        &
                ", CFL = ", cfl, ", time = ", time, ", cost = ", instantaneousCostFunctional
        end if
        call writeAndFlush(region%comm, output_unit, str)
        call functional%writeToFile(region%comm, trim(outputPrefix_) //                      &
             ".cost_functional.txt", timestep_, time, timestep_ > timestep + reportInterval)
     end if

     if (region%simulationFlags%steadyStateSimulation .and.                                  &
          mod(timestep_, residualInterval) == 0) then
        call region%computeResiduals(residuals)
        write(str, '(2X,3(A,(ES11.4E2)))') "residuals: density = ", residuals(1),            &
             ", momentum = ", residuals(2), ", energy = ", residuals(3)
        call writeAndFlush(region%comm, output_unit, str)
        if (all(residuals < region%solverOptions%convergenceTolerance)) then
           call writeAndFlush(region%comm, output_unit, "Solution has converged!")
           call region%saveData(QOI_FORWARD_STATE,                                    &
                trim(outputPrefix_) // ".steady_state.q")
           timestep = timestep_
           exit
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
              call region%saveData(QOI_FORWARD_STATE, filename)

           end if
        end if
     end if

     if (timestep_ == timestep + nTimesteps) timestep = timestep_

  end do

  call timeIntegratorFactory%cleanup()
  call functionalFactory%cleanup()

  call endTiming("solveForward")

end subroutine solveForward

subroutine solveAdjoint(region, time, timestep, nTimesteps, saveInterval, outputPrefix)

  ! <<< External modules >>>
  use iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use State_mod, only : t_State
  use Region_mod, only : t_Region
  use Functional_mod, only : t_Functional
  use Functional_factory, only : t_FunctionalFactory
  use TimeIntegrator_mod, only : t_TimeIntegrator
  use ReverseMigrator_type, only : t_ReverseMigrator
  use TimeIntegrator_factory, only : t_TimeIntegratorFactory

  ! <<< Enumerations >>>
  use State_enum, only : QOI_ADJOINT_STATE

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush
  use MPITimingsHelper, only : startTiming, endTiming
  use ReverseMigrator_mod

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region
  real(SCALAR_KIND), intent(inout) :: time
  integer, intent(inout) :: timestep
  integer, intent(in) :: nTimesteps, saveInterval
  character(len = STRING_LENGTH), intent(in), optional :: outputPrefix

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: outputPrefix_, filename, str
  type(t_ReverseMigrator) :: reverseMigrator
  integer :: i, j, timestep_, timemarchDirection, reportInterval, residualInterval
  class(t_Functional), pointer :: functional => null()
  type(t_FunctionalFactory) :: functionalFactory
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()
  type(t_TimeIntegratorFactory) :: timeIntegratorFactory
  real(wp) :: timeStepSize, cfl, residuals(3)

  assert(timestep >= nTimesteps)

  call startTiming("solveAdjoint")

  outputPrefix_ = PROJECT_NAME
  if (present(outputPrefix)) outputPrefix_ = outputPrefix

  write(filename, '(2A,I8.8,A)') trim(outputPrefix_), "-", timestep, ".adjoint.q"
  call region%saveData(QOI_ADJOINT_STATE, filename)

  residualInterval = getOption("check_residuals_interval", -1)
  if (residualInterval == -1) then
     call getRequiredOption("report_interval", reportInterval)
     residualInterval = reportInterval
     if (residualInterval < 0) residualInterval = max(1, nint(1e-3_wp * real(nTimesteps, wp)))
  else
     reportInterval = getOption("report_interval", 1)
  end if

  if (.not. region%simulationFlags%predictionOnly) then
     call functionalFactory%connect(functional, trim(region%solverOptions%costFunctionalType))
     assert(associated(functional))
     call functional%setup(region)
  end if

  call timeIntegratorFactory%connect(timeIntegrator,                                         &
       trim(region%solverOptions%timeIntegratorType))
  assert(associated(timeIntegrator))
  call timeIntegrator%setup(region)

  if (region%simulationFlags%steadyStateSimulation) then
     timemarchDirection = 1
  else
     timemarchDirection = -1
     call setupReverseMigrator(reverseMigrator, region, outputPrefix_,                     &
          getOption("checkpointing_scheme", "uniform checkpointing"),                           &
          timestep - nTimesteps, timestep,                                                      &
          saveInterval, saveInterval * timeIntegrator%nStages)
  end if

  do timestep_ = timestep + sign(1, timemarchDirection),                                     &
       timestep + sign(nTimesteps, timemarchDirection), timemarchDirection

     if (region%simulationFlags%useConstantCfl) then
        if (.not. region%simulationFlags%steadyStateSimulation)                              &
             call migrateToSubstep(reverseMigrator, region, timeIntegrator, timestep_, 1)
        do j = 1, size(region%states) !... update state
           call region%states(j)%update(region%grids(j), region%simulationFlags,             &
                region%solverOptions)
        end do
     end if
     timeStepSize = region%getTimeStepSize()

     do i = timeIntegrator%nStages, 1, -1

        if (.not. region%simulationFlags%steadyStateSimulation) then
           if (i == 1) then
              call migrateToSubstep(reverseMigrator, region,                                 &
                   timeIntegrator, timestep_, timeIntegrator%nStages)
           else
              call migrateToSubstep(reverseMigrator, region,                                 &
                   timeIntegrator, timestep_ + 1, i - 1)
           end if
        end if

        do j = 1, size(region%states) !... update state
           call region%states(j)%update(region%grids(j), region%simulationFlags,             &
                region%solverOptions)
        end do

        call functional%computeAdjointForcing(region)

        call timeIntegrator%substepAdjoint(region, time, timeStepSize, timestep_, i)

     end do

     if (reportInterval > 0 .and. mod(timestep_, reportInterval) == 0) then
        if (region%simulationFlags%useConstantCfl) then
           write(str, '(2A,I8,2(A,E13.6))') PROJECT_NAME, ": timestep = ", timestep_,        &
                ", dt = ", timeStepSize, ", time = ", time
        else
           cfl = region%getCfl()
           write(str, '(2A,I8,2(A,E13.6))') PROJECT_NAME, ": timestep = ", timestep_,        &
                ", CFL = ", cfl, ", time = ", time
        end if
        call writeAndFlush(region%comm, output_unit, str)
     end if

     if (region%simulationFlags%steadyStateSimulation .and.                                  &
          mod(timestep_, residualInterval) == 0) then
        call region%computeResiduals(residuals)
        write(str, '(2X,3(A,(ES11.4E2)))') "residuals: density = ", residuals(1),            &
             ", momentum = ", residuals(2), ", energy = ", residuals(3)
        call writeAndFlush(region%comm, output_unit, str)
        if (all(residuals < region%solverOptions%convergenceTolerance)) then
           call writeAndFlush(region%comm, output_unit, "Solution has converged!")
           call region%saveData(QOI_ADJOINT_STATE,                                    &
                trim(outputPrefix_) // ".steady_state.adjoint.q")
           timestep = timestep_
           exit
        end if
     end if

     if (saveInterval > 0) then
        if (mod(timestep_, saveInterval) == 0) then
           do i = 1, size(region%states)
              region%states(i)%plot3dAuxiliaryData(1) = real(timestep_, wp)
              region%states(i)%plot3dAuxiliaryData(4) = time
           end do
           write(filename, '(2A,I8.8,A)')                                                    &
                trim(outputPrefix_), "-", timestep_, ".adjoint.q"
           call region%saveData(QOI_ADJOINT_STATE, filename)
        end if
     end if

     if (timestep_ == timestep - nTimesteps) timestep = timestep_

  end do

  call cleanupReverseMigrator(reverseMigrator)

  call timeIntegratorFactory%cleanup()
  call functionalFactory%cleanup()

  call endTiming("solveAdjoint")

end subroutine solveAdjoint
