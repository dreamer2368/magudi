#include "config.h"

module SolverImpl

  implicit none
  public

contains

  subroutine showProgress(this, region, mode, startTimestep,                                 &
       timestep, time, instantaneousFunctional)

    ! <<< External modules >>>
    use iso_fortran_env, only : output_unit

    ! <<< Derived types >>>
    use Region_mod, only : t_Region
    use Solver_mod, only : t_Solver
    use Functional_mod, only : t_Functional

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE, QOI_ADJOINT_STATE
    use Region_enum, only : FORWARD, ADJOINT

    ! <<< Internal modules >>>
    use ErrorHandler, only : writeAndFlush

    ! <<< Arguments >>>
    class(t_Solver) :: this
    class(t_Region) :: region
    integer, intent(in) :: mode, timestep, startTimestep
    real(SCALAR_KIND), intent(in) :: time, instantaneousFunctional

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, nDimensions
    real(SCALAR_KIND) :: timeStepSize, cfl, residuals(3)
    character(len = STRING_LENGTH) :: str, str_, filename
    class(t_Functional), pointer :: functional => null()

    assert_key(mode, (FORWARD, ADJOINT))
    assert(startTimestep >= 0)

    nDimensions = size(region%globalGridSizes, 1)
    assert_key(nDimensions, (1, 2, 3))

    if (region%simulationFlags%useConstantCfl) then
       timeStepSize = region%getTimeStepSize()
    else
       cfl = region%getCfl()
    end if

    if (this%reportInterval > 0 .and. mod(timestep, max(1, this%reportInterval)) == 0) then

       if (region%simulationFlags%useConstantCfl) then
          write(str, '(2A,I8,2(A,E13.6))') PROJECT_NAME, ": timestep = ", timestep,          &
               ", dt = ", timeStepSize, ", time = ", abs(time)
       else
          write(str, '(2A,I8,2(A,E13.6))') PROJECT_NAME, ": timestep = ", timestep,          &
               ", CFL = ", cfl, ", time = ", abs(time)
       end if

       if (.not. region%simulationFlags%predictionOnly) then

          select case (mode)
          case (FORWARD)
             write(str_, '(A,E13.6)') ", cost = ", instantaneousFunctional
          case (ADJOINT)
             write(str_, '(A,E13.6)') ", gradient = ", instantaneousFunctional
          end select

          str = trim(str) // trim(str_)

       end if

       call writeAndFlush(region%comm, output_unit, str)

       if (.not. region%simulationFlags%predictionOnly) then

          call this%functionalFactory%connect(functional)
          assert(associated(functional))

          select case (mode)
          case (FORWARD)
             call functional%writeToFile(region%comm, trim(this%outputPrefix) //             &
                  ".cost_functional.txt", timestep, time,                                    &
                  timestep - startTimestep > this%reportInterval)
          end select

       end if

    end if

    if (this%saveInterval > 0 .and. mod(timestep, max(1, this%saveInterval)) == 0) then

       do i = 1, size(region%states)
          region%states(i)%plot3dAuxiliaryData(1) = real(timestep, wp)
          region%states(i)%time = time
       end do

       select case (mode)
       case (FORWARD)
          write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", timestep, ".q"
          call region%saveData(QOI_FORWARD_STATE, filename)
       case (ADJOINT)
          write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", timestep, ".adjoint.q"
          call region%saveData(QOI_ADJOINT_STATE, filename)
       end select

    end if

    if (region%simulationFlags%steadyStateSimulation .and.                                   &
         this%residualManager%reportInterval > 0 .and.                                       &
         mod(timestep, max(1, this%residualManager%reportInterval)) == 0) then

       call this%residualManager%compute(region)

       select case (mode)

       case (FORWARD)
          call this%residualManager%writeToFile(region%comm, trim(this%outputPrefix) //      &
               ".residuals.txt", timestep, time,                                             &
               timestep - startTimestep > this%residualManager%reportInterval)
       case (ADJOINT)
          call this%residualManager%writeToFile(region%comm, trim(this%outputPrefix) //      &
               ".adjoint_residuals.txt", timestep, time,                                     &
               timestep - startTimestep > this%residualManager%reportInterval)
       end select

       residuals(1) = this%residualManager%residuals(1)
       residuals(2) = maxval(this%residualManager%residuals(1:nDimensions))
       residuals(3) = this%residualManager%residuals(nDimensions+2)

       write(str, '(2X,3(A,(ES11.4E2)))') "residuals: density = ", residuals(1),             &
            ", momentum = ", residuals(2), ", energy = ", residuals(3)
       call writeAndFlush(region%comm, output_unit, str)

       if (this%residualManager%hasSimulationConverged) then

          call writeAndFlush(region%comm, output_unit, "Solution has converged!")

          select case (mode)
          case (FORWARD)
             call region%saveData(QOI_FORWARD_STATE,                                         &
                  trim(this%outputPrefix) // ".steady_state.q")
          case (ADJOINT)
             call region%saveData(QOI_ADJOINT_STATE,                                         &
                  trim(this%outputPrefix) // ".steady_state.adjoint.q")
          end select

       end if

    end if

  end subroutine showProgress

  subroutine checkSolutionLimits(region, mode)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use State_mod, only : t_State
    use Region_mod, only : t_Region

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE
    use Region_enum, only : FORWARD

    ! <<< Internal modules >>>
    use ErrorHandler, only : gracefulExit

    ! <<< Arguments >>>
    class(t_Region) :: region
    integer, intent(in) :: mode

    ! <<< Local variables >>>
    integer :: i, iGlobal, jGlobal, kGlobal, rankReportingError, procRank, ierror
    character(len = STRING_LENGTH) :: message
    SCALAR_TYPE :: fOutsideRange

    rankReportingError = -1
    call MPI_Comm_rank(region%comm, procRank, ierror)

    do i = 1, size(region%states)

       if (.not. region%grids(i)%isVariableWithinRange(                                      &
            region%states(i)%conservedVariables(:,1),                                        &
            fOutsideRange, iGlobal, jGlobal, kGlobal,                                        &
            minValue = region%solverOptions%densityRange(1),                                 &
            maxValue = region%solverOptions%densityRange(2))) then
          write(message, '(4(A,I0.0),3(A,(SS,ES9.2E2)),A)') "Density on grid ",              &
               region%grids(i)%index, " at (", iGlobal, ", ", jGlobal, ", ", kGlobal, "): ", &
               fOutsideRange, " out of range (",                                             &
               region%solverOptions%densityRange(1), ", ",                                   &
               region%solverOptions%densityRange(2), ")!"
          rankReportingError = procRank
          exit
       end if

       if (.not. region%grids(i)%isVariableWithinRange(region%states(i)%temperature(:,1),    &
            fOutsideRange, iGlobal, jGlobal, kGlobal,                                        &
            minValue = region%solverOptions%temperatureRange(1),                             &
            maxValue = region%solverOptions%temperatureRange(2))) then
          write(message, '(4(A,I0.0),3(A,(SS,ES9.2E2)),A)') "Temperature on grid ",          &
               region%grids(i)%index, " at (", iGlobal, ", ", jGlobal, ", ", kGlobal, "): ", &
               fOutsideRange, " out of range (",                                             &
               region%solverOptions%temperatureRange(1), ", ",                               &
               region%solverOptions%temperatureRange(2), ")!"
          rankReportingError = procRank
          exit
       end if

    end do

    call MPI_Allreduce(MPI_IN_PLACE, rankReportingError, 1,                                  &
         MPI_INTEGER, MPI_MAX, region%comm, ierror)

    if (rankReportingError /= -1) then

       if (procRank == 0 .and. rankReportingError /= 0)                                      &
            call MPI_Recv(message, STRING_LENGTH, MPI_CHARACTER, rankReportingError,         &
            rankReportingError, region%comm, MPI_STATUS_IGNORE, ierror)
       if (procRank == rankReportingError .and. rankReportingError /= 0)                     &
            call MPI_Send(message, STRING_LENGTH, MPI_CHARACTER, 0, procRank,                &
            region%comm, ierror)

       select case (mode)
       case (FORWARD)
          call region%saveData(QOI_FORWARD_STATE, PROJECT_NAME // "-crashed.q")
       end select

       call gracefulExit(region%comm, message)

    end if

  end subroutine checkSolutionLimits

end module SolverImpl

subroutine setupSolver(this, region, restartFilename, outputPrefix)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver
  use Functional_mod, only : t_Functional
  use TimeIntegrator_mod, only : t_TimeIntegrator

  ! <<< Enumerations >>>
  use Grid_enum, only : QOI_CONTROL_MOLLIFIER, QOI_TARGET_MOLLIFIER
  use State_enum, only : QOI_FORWARD_STATE, QOI_TARGET_STATE, QOI_ADJOINT_STATE

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_Solver) :: this
  class(t_Region) :: region
  character(len = *), intent(in), optional :: restartFilename, outputPrefix

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: filename
  integer :: i
  class(t_Functional), pointer :: functional => null()
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()

  if (present(outputPrefix)) then
     this%outputPrefix = outputPrefix
  else
     this%outputPrefix = getOption("output_prefix", PROJECT_NAME)
  end if

  this%saveInterval = getOption("save_interval", -1)
  if (this%saveInterval == 0) this%saveInterval = -1

  this%reportInterval = getOption("report_interval", 1)
  if (this%reportInterval == 0) this%reportInterval = -1

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
  else if (region%simulationFlags%predictionOnly .or. .not.                                  &
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
        call getRequiredOption("initial_condition_file", filename, region%comm)
        call region%loadData(QOI_FORWARD_STATE, filename) !... initialize from file.
     end if
  end if

  if (.not. region%simulationFlags%predictionOnly) then

     ! Initialize adjoint variables.
     filename = getOption("adjoint_initial_condition_file", "")
     if (len_trim(filename) == 0) then
        do i = 1, size(region%states)
           region%states(i)%adjointVariables = 0.0_wp
        end do
     else
        call region%loadData(QOI_ADJOINT_STATE, filename)
     end if

     ! Initialize control mollifier.
     filename = getOption("control_mollifier_file", "")
     if (len_trim(filename) == 0) then
        do i = 1, size(region%grids)
           region%grids(i)%controlMollifier = 1.0_wp
        end do
     else
        call region%loadData(QOI_CONTROL_MOLLIFIER, filename)
     end if

     ! Target mollifier.
     filename = getOption("target_mollifier_file", "")
     if (len_trim(filename) == 0) then
        do i = 1, size(region%grids)
           region%grids(i)%targetMollifier = 1.0_wp
        end do
     else
        call region%loadData(QOI_TARGET_MOLLIFIER, filename)
     end if

     call this%functionalFactory%connect(functional,                                         &
          trim(region%solverOptions%costFunctionalType))
     assert(associated(functional))
     call functional%setup(region)

  end if

  call this%timeIntegratorFactory%connect(timeIntegrator,                                    &
       trim(region%solverOptions%timeIntegratorType))
  assert(associated(timeIntegrator))
  call timeIntegrator%setup(region)

end subroutine setupSolver

subroutine cleanupSolver(this)

  ! <<< Derived types >>>
  use Solver_mod, only : t_Solver

  implicit none

  ! <<< Arguments >>>
  class(t_Solver) :: this

  call this%timeIntegratorFactory%cleanup()
  call this%functionalFactory%cleanup()

end subroutine cleanupSolver

function runForward(this, region, time, timestep, nTimesteps) result(costFunctional)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver
  use Functional_mod, only : t_Functional
  use TimeIntegrator_mod, only : t_TimeIntegrator

  ! <<< Enumerations >>>
  use State_enum, only : QOI_FORWARD_STATE
  use Region_enum, only : FORWARD

  ! <<< Private members >>>
  use SolverImpl, only : showProgress, checkSolutionLimits

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_Solver) :: this
  class(t_Region) :: region
  real(SCALAR_KIND), intent(inout) :: time
  integer, intent(inout) :: timestep
  integer, intent(in) :: nTimesteps

  ! <<< Result >>>
  SCALAR_TYPE :: costFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: filename
  integer :: i, j, timestep_
  SCALAR_TYPE :: instantaneousCostFunctional
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()
  class(t_Functional), pointer :: functional => null()
  real(SCALAR_KIND) :: timeStepSize

  assert(timestep >= 0)

  call startTiming("runForward")

  costFunctional = 0.0_wp

  call this%timeIntegratorFactory%connect(timeIntegrator)
  assert(associated(timeIntegrator))

  if (.not. region%simulationFlags%predictionOnly) then
     call this%functionalFactory%connect(functional)
     assert(associated(functional))
  end if

  if (region%simulationFlags%steadyStateSimulation)                                          &
       call this%residualManager%setup("", region)

  write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", timestep, ".q"
  call region%saveData(QOI_FORWARD_STATE, filename)

  do timestep_ = timestep + 1, timestep + nTimesteps

     do j = 1, size(region%states) !... update state
        region%states(j)%time = time
        call region%states(j)%update(region%grids(j), region%simulationFlags,                &
             region%solverOptions)
     end do

     timeStepSize = region%getTimeStepSize()

     do i = 1, timeIntegrator%nStages

        if (region%simulationFlags%enableSolutionLimits)                                     &
             call checkSolutionLimits(region, FORWARD)

        call timeIntegrator%substepForward(region, time, timeStepSize, timestep_, i)

        if (i /= timeIntegrator%nStages) then
           do j = 1, size(region%states) !... update state
              region%states(j)%time = time
              call region%states(j)%update(region%grids(j), region%simulationFlags,          &
                   region%solverOptions)
           end do
        end if

        if (.not. region%simulationFlags%predictionOnly) then
           instantaneousCostFunctional = functional%compute(region)
           costFunctional = costFunctional +                                                 &
                timeIntegrator%norm(i) * timeStepSize * instantaneousCostFunctional
        end if

     end do

     call showProgress(this, region, FORWARD, timestep, timestep_,                           &
          time, instantaneousCostFunctional)

     if (this%residualManager%hasSimulationConverged) then
        timestep = timestep_
        exit
     end if

     if (timestep_ == timestep + nTimesteps) timestep = timestep_

  end do

  call this%residualManager%cleanup()

  call endTiming("runForward")

end function runForward

function runAdjoint(this, region, time, timestep, nTimesteps) result(costSensitivity)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver
  use Functional_mod, only : t_Functional
  use TimeIntegrator_mod, only : t_TimeIntegrator
  use ReverseMigrator_mod, only : t_ReverseMigrator
  use ReverseMigrator_factory, only : t_ReverseMigratorFactory

  ! <<< Enumerations >>>
  use State_enum, only : QOI_ADJOINT_STATE
  use Region_enum, only : ADJOINT

  ! <<< Private members >>>
  use SolverImpl, only : showProgress, checkSolutionLimits

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_Solver) :: this
  class(t_Region) :: region
  real(SCALAR_KIND), intent(inout) :: time
  integer, intent(inout) :: timestep
  integer, intent(in) :: nTimesteps

  ! <<< Result >>>
  SCALAR_TYPE :: costSensitivity

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: filename
  integer :: i, j, timestep_, timemarchDirection
  SCALAR_TYPE :: instantaneousCostSensitivity
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()
  class(t_Functional), pointer :: functional => null()
  type(t_ReverseMigratorFactory) :: reverseMigratorFactory
  class(t_ReverseMigrator), pointer :: reverseMigrator => null()
  real(SCALAR_KIND) :: timeStepSize

  assert(.not. region%simulationFlags%predictionOnly)

  call startTiming("runAdjoint")

  costSensitivity = 0.0_wp

  call this%timeIntegratorFactory%connect(timeIntegrator)
  assert(associated(timeIntegrator))

  call this%functionalFactory%connect(functional)
  assert(associated(functional))

  call reverseMigratorFactory%connect(reverseMigrator,                                       &
       region%solverOptions%checkpointingScheme)
  assert(associated(reverseMigrator))

  if (region%simulationFlags%steadyStateSimulation) then

     time = 0.0_wp
     timestep = 0
     timemarchDirection = 1

     call this%residualManager%setup("adjoint_residuals", region)

     do j = 1, size(region%states) !... update state
        region%states(j)%time = time
        call region%states(j)%update(region%grids(j), region%simulationFlags,                &
             region%solverOptions)
     end do

  else

     timemarchDirection = -1

     call reverseMigrator%setup(region, timeIntegrator, this%outputPrefix,                   &
          timestep - nTimesteps, timestep, this%saveInterval,                                &
          this%saveInterval * timeIntegrator%nStages)

  end if

  write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", timestep, ".adjoint.q"
  call region%saveData(QOI_ADJOINT_STATE, filename)

  do timestep_ = timestep + sign(1, timemarchDirection),                                     &
       timestep + sign(nTimesteps, timemarchDirection), timemarchDirection

     if (region%simulationFlags%useConstantCfl .and.                                         &
          .not. region%simulationFlags%steadyStateSimulation) then

        call reverseMigrator%migrateTo(region, timeIntegrator, timestep_, 1)

        do j = 1, size(region%states) !... update state
           region%states(j)%time = time
           call region%states(j)%update(region%grids(j), region%simulationFlags,             &
                region%solverOptions)
        end do

     end if

     timeStepSize = region%getTimeStepSize()

     do i = timeIntegrator%nStages, 1, -1

        if (.not. region%simulationFlags%steadyStateSimulation) then

           if (i == 1) then
              call reverseMigrator%migrateTo(region, timeIntegrator,                         &
                   timestep_, timeIntegrator%nStages)
           else
              call reverseMigrator%migrateTo(region, timeIntegrator, timestep_ + 1, i - 1)
           end if

           do j = 1, size(region%states) !... update state
              region%states(j)%time = time
              call region%states(j)%update(region%grids(j), region%simulationFlags,          &
                   region%solverOptions)
           end do

        end if

        call functional%updateAdjointForcing(region)

        call timeIntegrator%substepAdjoint(region, time, timeStepSize, timestep_, i)

        if (region%simulationFlags%enableSolutionLimits)                                     &
             call checkSolutionLimits(region, ADJOINT)

     end do

     call showProgress(this, region, ADJOINT, timestep, timestep_,                           &
          time, instantaneousCostSensitivity)

     if (this%residualManager%hasSimulationConverged) then
        timestep = timestep_
        exit
     end if

     if (timestep_ == timestep - nTimesteps) timestep = timestep_

  end do

  call this%residualManager%cleanup()
  call reverseMigratorFactory%cleanup()

  call endTiming("solveAdjoint")

end function runAdjoint
