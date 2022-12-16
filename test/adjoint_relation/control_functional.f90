#include "config.h"

program control_functional

  use MPI

  use CNSHelper
  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random

  implicit none

  logical :: success, success_, isPeriodic
  integer :: i, j, nDimensions, ierror, stat
  integer :: procRank
  character(len = STRING_LENGTH) :: costType
  character(len = STRING_LENGTH) :: controllerTypes(3)

  interface

     real(SCALAR_KIND) function meanTrimmed(a)
       real(SCALAR_KIND), intent(inout) :: a(:)
     end function meanTrimmed

     subroutine sort(a)
       real(SCALAR_KIND), intent(inout) :: a(:)
     end subroutine sort

     subroutine testAdjointRelation(costType, controllerType, nDimensions, success, isPeriodic, tolerance)

       character(len=*), intent(in) :: costType
       character(len=*), intent(in) :: controllerType
       integer, intent(in) :: nDimensions
       logical, intent(out) :: success

       logical, intent(in) :: isPeriodic
       real(SCALAR_KIND), intent(in), optional :: tolerance

     end subroutine testAdjointRelation

  end interface

  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  call get_command_argument(1, costType, stat)
  if( stat.eq.0 ) costType = "SOUND"

  controllerTypes(1) = "MOMENTUM_ACTUATOR"
  controllerTypes(2) = "THERMAL_ACTUATOR"
  controllerTypes(3) = "GENERIC_ACTUATOR"

  do nDimensions = 2, 2
    do j = 1, 1 !... for each discretizationTypes
      success = .true.
      do i = 1, 3 !... test multiple times

        isPeriodic = .false.
        call testAdjointRelation(costType, controllerTypes(i), nDimensions, success_, isPeriodic)
        success = success .and. success_
        if( .not. success_) then
          if( procRank == 0 ) then
            print *, 'Failed'
            print *, 'Controller: ', trim(controllerTypes(i))
            print *, 'periodicity: ', isPeriodic
          end if
          exit
        end if
      end do
      if( procRank == 0 .and. success ) then
        print *, 'Success'
      end if
    end do
  end do

  call cleanupErrorHandler()

  call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierror)
  call MPI_Finalize(ierror)
  if (.not. success) stop -1
  stop 0

end program control_functional

subroutine testAdjointRelation(costType, controllerType, nDimensions, success, isPeriodic, tolerance)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use StencilOperator_mod, only : t_StencilOperator
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags
  use Region_mod, only : t_Region

  use TimeIntegrator_mod, only : t_TimeIntegrator
  use TimeIntegrator_factory, only : t_TimeIntegratorFactory
  use ReverseMigrator_mod, only : t_ReverseMigrator
  use ReverseMigrator_factory, only : t_ReverseMigratorFactory

  use Controller_factory, only : t_ControllerFactory
  use Functional_factory, only : t_FunctionalFactory
  use Controller_mod, only : t_Controller
  use Functional_mod, only : t_Functional
  use Patch_mod, only : t_Patch
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use LighthillTensorComponent_mod, only : t_LighthillTensorComponent
  use LighthillSource_mod, only : t_LighthillSource

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT
  use State_enum, only : QOI_FORWARD_STATE, QOI_TIME_AVERAGED_STATE, QOI_ADJOINT_STATE

  ! <<< Internal modules >>>
  use MPIHelper, only : pigeonhole
  use RandomNumber, only : random
  use RegionImpl, only : normalizeTargetMollifier, normalizeControlMollifier
  use Patch_factory, only : updatePatchFactories
  use CNSHelper
  use SolverImpl, only : loadInitialCondition

  ! <<< Arguments >>>
  character(len=*), intent(in) :: costType, controllerType
  integer, intent(in) :: nDimensions
  logical, intent(out) :: success
  logical, intent(in) :: isPeriodic
  real(SCALAR_KIND), intent(in), optional :: tolerance

  ! <<< interface >>>
  interface
     real(SCALAR_KIND) function meanTrimmed(a)
       real(SCALAR_KIND), intent(inout) :: a(:)
     end function meanTrimmed

     subroutine sort(a)
       real(SCALAR_KIND), intent(inout) :: a(:)
     end subroutine sort
  end interface

  ! <<< Local derived type variables >>>
  type(t_SimulationFlags) :: simulationFlags
  type(t_SolverOptions) :: solverOptions
  type(t_Region) :: region
  type(t_State) :: icState, targetState
  class(t_Patch), pointer :: patch => null()

  type(t_TimeIntegratorFactory) :: timeIntegratorFactory
  class(t_TimeIntegrator), pointer :: timeIntegrator => null()
  type(t_ReverseMigratorFactory) :: reverseMigratorFactory
  class(t_ReverseMigrator), pointer :: reverseMigrator => null()

  type(t_ControllerFactory) :: controllerFactory
  type(t_FunctionalFactory) :: functionalFactory
  class(t_Controller), pointer :: controller => null()
  class(t_Functional), pointer :: functional => null()

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  logical :: IS_FINAL_STEP
  real(wp) :: scalar1, sensitivity, scalar2, tolerance_,                        &
              stepSizes(32), errorHistory(32), convergenceHistory(31)
  integer, allocatable :: gridSize(:,:)
  integer :: i, j, k, l, nTimesteps, saveInterval,                      &
              timestep, startTimestep, timemarchDirection
  real(wp) :: timeStepSize, time, startTime, velMin, velMax
  character(len = STRING_LENGTH) :: filename
  character(len = STRING_LENGTH), parameter :: discretizationTypes(4) =                      &
       (/ "SBP 1-2", "SBP 2-4", "SBP 3-6", "SBP 4-8" /)

  success = .true.
  tolerance_ = 1.0E-13
  if (present(tolerance)) tolerance_ = tolerance

  ! set up simulation flags
  call simulationFlags%initialize()
  simulationFlags%enableController = .true.
  simulationFlags%enableFunctional = .true.
  simulationFlags%enableAdjoint = .true.
  simulationFlags%useTargetState = .true.
  simulationFlags%isDomainCurvilinear = .false.
  simulationFlags%viscosityOn = .false.
  simulationFlags%useConstantCfl = .false.

  ! randomize grid size
  ! Note that too small grid size will yield null matrix for stencil operators.
  allocate(gridSize(nDimensions,1))
  gridSize = 1
  do i = 1, nDimensions
     gridSize(i,1) = 21
  end do

  ! initialize solver option and region
  simulationFlags%useConstantCfl = .true.
  call solverOptions%initialize(nDimensions, simulationFlags)
  simulationFlags%useConstantCfl = .false.
  ! solverOptions%ratioOfSpecificHeats = random(1.1_wp,2.0_wp)
  solverOptions%discretizationType = trim(discretizationTypes(random(1,4)))
  solverOptions%timeStepSize = random(0.001_wp, 0.01_wp)
  solverOptions%controllerType = trim(controllerType)
  solverOptions%costFunctionalType = trim(costType)
  call region%setup(MPI_COMM_WORLD,gridSize,simulationFlags,solverOptions)
  do i = 1, size(region%grids)
     call region%grids(i)%update()
  end do

  !initialize initial condition state
  call icState%setup(region%grids(1), region%simulationFlags, region%solverOptions)
  call targetState%setup(region%grids(1), region%simulationFlags, region%solverOptions)

  velMax = 1.0E-5
  velMin = -1.0E-5
  ! Randomize conserved variables.
  do k = 1, size(region%states)
    do i = 1, region%grids(k)%nGridPoints
       region%states(1)%conservedVariables(i,1) = random(0.99_wp, 1.01_wp)
       do j = 1, nDimensions
          region%states(1)%conservedVariables(i,j+1) =                                              &
               region%states(1)%conservedVariables(i,1) * random(velMin,velMax)
       end do
       region%states(1)%conservedVariables(i, nDimensions + 2) =                                    &
            region%states(1)%conservedVariables(i,1) *                                              &
            random(2.49_wp, 2.51_wp) / region%solverOptions%ratioOfSpecificHeats +                  &
            0.5_wp / region%states(1)%conservedVariables(i,1) *                                     &
            sum(region%states(1)%conservedVariables(i,2:nDimensions+1) ** 2)
    end do
  end do
  assert(all(region%states(1)%conservedVariables(:,1) > 0.0_wp))
  icState%conservedVariables = region%states(1)%conservedVariables

  ! Randomize initial condition: uniform over space, unknowns
  do k = 1, size(region%states)
    do i = 1, region%grids(k)%nGridPoints
       targetState%conservedVariables(i,1) = random(0.99_wp, 1.01_wp)
       ! targetState%conservedVariables(i,1) = 1.0_wp
       do j = 1, nDimensions
          targetState%conservedVariables(i,j+1) =                                              &
               targetState%conservedVariables(i,1) * random(velMin,velMax)
          ! targetState%conservedVariables(i,j+1) = 0.0_wp
       end do
       targetState%conservedVariables(i, nDimensions + 2) = 1.0/1.4/0.4
    end do
  end do
  assert(all(targetState%conservedVariables(:,1) > 0.0_wp))
  region%states(1)%targetState = targetState%conservedVariables

  ! control mollifier
  do i = 1, size(region%grids)
     region%grids(i)%controlMollifier = 1.0_wp
  end do
  ! target mollifier
  do i = 1, size(region%grids)
     region%grids(i)%targetMollifier = 1.0_wp
  end do

  ! Setup boundary conditions (controller and functional only)
  call region%setupBoundaryConditions("bc.dat")
  do j = 1, size(region%patchFactories)
     call region%patchFactories(j)%connect(patch)
     if (.not. associated(patch)) cycle
     do l = 1, size(region%states)
        if (patch%gridIndex /= region%grids(l)%index) cycle
        select type (patch)
        class is (t_ActuatorPatch)
           patch%gradientFilename = "test.gradient_controlRegion.dat"
           patch%controlForcingFilename = "test.gradient_controlRegion.dat"
           if( patch%nPatchPoints > 0 ) then
             allocate(patch%controlForcing(patch%nPatchPoints, solverOptions%nUnknowns))
             patch%controlForcing = 0.0_wp
           end if
        end select
     end do
  end do

  ! Update patches.
  do i = 1, size(region%grids)
     call updatePatchFactories(region%patchFactories, region%simulationFlags,                &
          region%solverOptions, region%grids(i), region%states(i))
  end do

  ! randomize timesteps
  nTimesteps = 6
  saveInterval = 2
  nTimesteps = nTimesteps - mod(nTimesteps,saveInterval)

  ! set up controller
  call controllerFactory%connect(controller,                                         &
       trim(region%solverOptions%controllerType))
  assert(associated(controller))
  controller%controllerBufferSize = 4*saveInterval
  call controller%setup(region)

  ! set up functional
  call functionalFactory%connect(functional,                                         &
       trim(region%solverOptions%costFunctionalType))
  assert(associated(functional))
  call functional%setup(region)
  select type (functional)
  class is (t_LighthillTensorComponent)
    call random_number(functional%firstDirection)
    call random_number(functional%secondDirection)
    functional%firstDirection = 2.0_wp * functional%firstDirection - 1.0_wp
    functional%secondDirection = 2.0_wp * functional%secondDirection - 1.0_wp
    functional%firstDirection(1:nDimensions) = functional%firstDirection(1:nDimensions)     &
                               / sqrt( sum(functional%firstDirection(1:nDimensions)**2) )
    functional%secondDirection(1:nDimensions) = functional%firstDirection(1:nDimensions)    &
                               / sqrt( sum(functional%secondDirection(1:nDimensions)**2) )
  end select
  functional%runningTimeQuadrature = 0.0_wp

  ! initialize time integrator
  call timeIntegratorFactory%connect(timeIntegrator,                                    &
       trim(region%solverOptions%timeIntegratorType))
  assert(associated(timeIntegrator))
  call timeIntegrator%setup(region)

  ! ! allocate time-dependeant variables
  ! allocate(forwardState(9,nTimesteps*timeIntegrator%nStages+1))
  ! allocate(adjointState(9,nTimesteps*timeIntegrator%nStages+1))
  ! allocate(deltaState(9,nTimesteps*timeIntegrator%nStages+1))
  ! allocate(controlForcing(9,nTimesteps*timeIntegrator%nStages+1))
  ! allocate(adjointForcing(9,nTimesteps*timeIntegrator%nStages+1))
  ! forwardState = 0.0_wp
  ! adjointState = 0.0_wp
  ! deltaState = 0.0_wp
  ! controlForcing = 0.0_wp
  ! adjointForcing = 0.0_wp
  ! forwardState(:,1) = reshape(region%states(1)%conservedVariables,(/9/))
  !
  ! allocate(forwardRhs(nTimesteps*timeIntegrator%nStages))
  ! allocate(region%tempRhs(nTimesteps*timeIntegrator%nStages))
  !
  ! ! randomize RHS
  ! call random_number(forwardRhs)
  ! forwardRhs = forwardRhs * random(10.0_wp,100.0_wp)

  ! Forward run setup
  time = 0.0_wp
  region%states(:)%time = time
  region%timestep = 0

  write(filename, '(A,I8.8,A)') "test-", region%timestep, ".q"
  call region%saveData(QOI_FORWARD_STATE, filename)

  controller%controllerSwitch = .false.
  if (controller%controllerSwitch) then
     controller%onsetTime = time
     controller%duration = nTimesteps * region%solverOptions%timeStepSize
     call controller%hookBeforeTimemarch(region, FORWARD)
  end if

  do i = 1, size(region%states) !... update state
     call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
          region%solverOptions)
  end do

  ! Forward run: Dx = f
  ! region%tempRhs = forwardRhs
  scalar1 = 0.0_wp

  do timestep = 1, nTimesteps

    region%timestep = timestep
    timeStepSize = region%getTimeStepSize()

    do i = 1, timeIntegrator%nStages
      ! Update control forcing.
      if (controller%controllerSwitch) then
        call controller%cleanupForcing(region)
        call controller%updateForcing(region)
      end if

      ! Take a single sub-step using the time integrator.
      call timeIntegrator%substepForward(region, time, timeStepSize, timestep, i)

      do j = 1, size(region%states) !... update state
         call region%states(j)%update(region%grids(j), region%simulationFlags,             &
              region%solverOptions)
      end do

      ! ! Save forward states
      ! forwardState(:,(timestep-1)*timeIntegrator%nStages+i+1) = reshape(region%states(1)%conservedVariables,(/9/))
! print *, functional%compute(region)
      ! Update the cost functional.
      scalar1 = scalar1 + timeIntegrator%norm(i) * timeStepSize * functional%compute(region)

    end do

    if (mod(timestep, max(1, saveInterval)) == 0) then
          write(filename, '(A,I8.8,A)') "test-", timestep, ".q"
          call region%saveData(QOI_FORWARD_STATE, filename)
    end if

  end do

  ! Call controller hooks after time marching ends.
  if (controller%controllerSwitch) call controller%hookAfterTimemarch(region, FORWARD)

  print *, 'baseline cost functional: ', scalar1

  ! Adjoint run = D^{\dagger}x^{\dagger} = g

  ! Load the adjoint coefficients corresponding to the end of the control time horizon.
  write(filename, '(A,I8.8,A)') "test-", 0+nTimesteps, ".q"
  call region%loadData(QOI_FORWARD_STATE, filename)
  do i = 1, size(region%states) !... update state
     call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
          region%solverOptions)
  end do
  startTimestep = region%timestep
  startTime = region%states(1)%time

  ! Connect to the previously allocated reverse migrator.
  call reverseMigratorFactory%connect(reverseMigrator,                                       &
       region%solverOptions%checkpointingScheme)
  assert(associated(reverseMigrator))
  ! Setup the reverse-time migrator
  call reverseMigrator%setup(region, timeIntegrator, "test",                       &
    startTimestep - nTimesteps, startTimestep, saveInterval,                    &
    saveInterval * timeIntegrator%nStages)

  timemarchDirection = -1

  ! adjoint initial condition
  region%states(:)%adjointForcingFactor = - timeStepSize / 6.0_wp !... RK4 only.
  call functional%updateAdjointForcing(region,.false.) !...SeungWhan:obviously not final step
  do i = 1, size(region%states)
     region%states(i)%rightHandSide = 0.0_wp
  end do
  do i = 1, size(region%patchFactories)
     call region%patchFactories(i)%connect(patch)
     if (.not. associated(patch)) cycle
     select type (patch)
        class is (t_CostTargetPatch)
        do j = 1, size(region%states)
           if (patch%gridIndex /= region%grids(j)%index) cycle
           ! adjointForcing(:,size(adjointForcing,2)) = reshape(patch%adjointForcing,(/9/))
           call patch%updateRhs(ADJOINT, region%simulationFlags,                     &
                region%solverOptions, region%grids(j), region%states(j))
        end do
     end select
  end do
  do i = 1, size(region%states)
     region%states(i)%adjointVariables = region%states(i)%rightHandSide
  end do
  region%states(:)%adjointForcingFactor = 1.0_wp !... restore
  write(filename, '(A,I8.8,A)') "test-", region%timestep, ".adjoint.q"
  call region%saveData(QOI_ADJOINT_STATE, filename)
  ! adjointState(:,size(adjointState,2)) = reshape(region%states(1)%adjointVariables,(/9/))

  ! Call controller hooks before time marching starts.
  call controller%hookBeforeTimemarch(region, ADJOINT)
  !!!SeungWhan: need additional execution with FORWARD, in case of non-zero control forcing.
  controller%controllerSwitch = .false.
  call controller%cleanupForcing(region)
  if (controller%controllerSwitch) then
    controller%duration = nTimesteps * region%solverOptions%timeStepSize
    controller%onsetTime = startTime - controller%duration
    call controller%hookBeforeTimemarch(region, FORWARD)
  end if

  time = startTime

  sensitivity = 0.0_wp

  do timestep = startTimestep + sign(1, timemarchDirection),&
      startTimestep + sign(nTimesteps, timemarchDirection), timemarchDirection

    region%timestep = timestep
    timeStepSize = region%getTimeStepSize()

    do i = timeIntegrator%nStages, 1, -1
      ! Load adjoint coefficients.
      if (i == 1) then
         call reverseMigrator%migrateTo(region, controller, timeIntegrator,             &
              timestep, timeIntegrator%nStages)
      else
         call reverseMigrator%migrateTo(region, controller, timeIntegrator,             &
              timestep + 1, i - 1)
      end if

      ! Update gradient.
      call controller%updateGradient(region)

      ! Update cost sensitivity.
      sensitivity = sensitivity + timeIntegrator%norm(i) * timeStepSize * controller%computeSensitivity(region)

      ! Update adjoint forcing on cost target patches.
      ! SeungWhan: Bug fix for final step
      if( (timestep.eq.startTimestep+sign(nTimesteps,timemarchDirection)) .and.       &
          (i.eq.1) ) then
          IS_FINAL_STEP = .true.
      else
          IS_FINAL_STEP = .false.
      end if
      call functional%updateAdjointForcing(region,IS_FINAL_STEP)
      ! record adjoint forcing
      ! do j = 1, size(region%patchFactories)
      !    call region%patchFactories(j)%connect(patch)
      !    if (.not. associated(patch)) cycle
      !    do l = 1, size(region%states)
      !       if (patch%gridIndex /= region%grids(l)%index) cycle
      !       select type (patch)
      !       class is (t_CostTargetPatch)
      !         adjointForcing(:,timestep*timeIntegrator%nStages+i) = reshape(patch%adjointForcing,(/9/))
      !       end select
      !    end do
      ! end do

      ! Take a single sub-step using the time integrator.
      call timeIntegrator%substepAdjoint(region, time, timeStepSize, timestep, i)

      ! ! Save adjoint states
      ! adjointState(:,timestep*timeIntegrator%nStages+i) = reshape(region%states(1)%adjointVariables,(/9/))
    end do

    if (mod(timestep, max(1, saveInterval)) == 0) then
          write(filename, '(A,I8.8,A)') "test-", timestep, ".adjoint.q"
          call region%saveData(QOI_ADJOINT_STATE, filename)
    end if

  end do

  ! Call controller hooks after time marching ends.
  if (controller%controllerSwitch) call controller%hookAfterTimemarch(region, FORWARD)
  call controller%hookAfterTimemarch(region, ADJOINT)

  print *, 'baseline sensitivity: ', sensitivity

  ! Do finite difference approximation
  ! stepSizes(1) = 1.0/sensitivity
  stepSizes(1) = 1.0E04
  do k = 2, size(stepSizes)
     stepSizes(k) = stepSizes(k-1) * 10.0_wp**(-0.25_wp)
  end do
  errorHistory = 0.0_wp
  do k = 1, size(stepSizes)
    ! forward initial condition
    region%states(1)%conservedVariables = icState%conservedVariables
    region%states(1)%targetState = targetState%conservedVariables

    ! Forward run setup
    time = 0.0_wp
    region%states(:)%time = time
    region%timestep = 0

    controller%controllerSwitch = .true.
    if (controller%controllerSwitch) then
       controller%onsetTime = time
       controller%duration = nTimesteps * region%solverOptions%timeStepSize
       call controller%hookBeforeTimemarch(region, FORWARD)
    end if

    do i = 1, size(region%states) !... update state
       call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
            region%solverOptions)
    end do

    ! Forward run: Dx = f
    ! region%tempRhs = forwardRhs
    ! region%tempRhsIndex = 1
    scalar2 = 0.0_wp

    do timestep = 1, nTimesteps

      region%timestep = timestep
      timeStepSize = region%getTimeStepSize()

      do i = 1, timeIntegrator%nStages
        ! Update control forcing.
        if (controller%controllerSwitch) then
          call controller%cleanupForcing(region)
          call controller%updateForcing(region)
        end if

        ! take a step size of control forcing
        do j = 1, size(region%patchFactories)
           call region%patchFactories(j)%connect(patch)
           if (.not. associated(patch)) cycle
           do l = 1, size(region%states)
              if (patch%gridIndex /= region%grids(l)%index) cycle
              select type (patch)
              class is (t_ActuatorPatch)
                 patch%controlForcing = patch%controlForcing * stepSizes(k)
                 ! controlForcing(:,(timestep-1)*timeIntegrator%nStages+i+1) = reshape(patch%controlForcing,(/9/))
              end select
           end do
        end do

        ! Take a single sub-step using the time integrator.
        call timeIntegrator%substepForward(region, time, timeStepSize, timestep, i)

        do j = 1, size(region%states) !... update state
           call region%states(j)%update(region%grids(j), region%simulationFlags,             &
                region%solverOptions)
        end do

        ! Update the cost functional.
        scalar2 = scalar2 + timeIntegrator%norm(i) * timeStepSize * functional%compute(region)

        ! deltaState(:,(timestep-1)*timeIntegrator%nStages+i+1) = reshape(region%states(1)%conservedVariables,(/9/)) - forwardState(:,(timestep-1)*timeIntegrator%nStages+i+1)
      end do

    end do

    ! Call controller hooks after time marching ends.
    if (controller%controllerSwitch) call controller%hookAfterTimemarch(region, FORWARD)

    errorHistory(k) = abs( ((scalar2- scalar1)/stepSizes(k) - sensitivity)/sensitivity )
    print *, stepSizes(k), (scalar2-scalar1)/stepSizes(k), sensitivity, errorHistory(k)

    if (k > 1) then
       convergenceHistory(k-1) = log(errorHistory(k) / errorHistory(k-1)) /              &
            log(stepSizes(k) / stepSizes(k-1))
       if (k > 5) then
         if(sum(convergenceHistory(k-3:k-1))/3.0_wp < 0.0_wp) exit
       end if
    end if
  end do

  if (maxval(errorHistory).le.tolerance_) then
    success = .true.
  elseif (k > 3) then
     call sort(convergenceHistory(:k-2))
     success = success .and. nint(meanTrimmed(convergenceHistory(:k-2))).ge.1
  else
     success = .false.
     print *, 'minimum error: ', minval(errorHistory)
  end if

  ! SAFE_DEALLOCATE(forwardState)
  ! SAFE_DEALLOCATE(forwardRhs)
  ! SAFE_DEALLOCATE(adjointState)
  ! SAFE_DEALLOCATE(deltaState)
  ! SAFE_DEALLOCATE(controlForcing)
  ! SAFE_DEALLOCATE(adjointForcing)
  ! SAFE_DEALLOCATE(region%tempRhs)

  call patch%cleanup()
  call controller%cleanup()
  call functional%cleanup()
  call controllerFactory%cleanup()
  call functionalFactory%cleanup()
  call timeIntegrator%cleanup()
  call reverseMigrator%cleanup()
  call timeIntegratorFactory%cleanup()
  call reverseMigratorFactory%cleanup()
  call icState%cleanup()
  call targetState%cleanup()
  call region%cleanup()

  SAFE_DEALLOCATE(gridSize)

end subroutine testAdjointRelation

subroutine sort(a)

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(inout) :: a(:)

  ! <<< Local variables >>>
  integer :: i, j
  real(SCALAR_KIND) :: temp

  do i = 2, size(a)
     j = i - 1
     temp = a(i)
     do while (a(j) > temp)
        a(j+1) = a(j)
        j = j - 1
        if (j < 1) exit
     end do
     a(j+1) = temp
  end do

end subroutine sort

real(SCALAR_KIND) function median(a)

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(in) :: a(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: n

  n = size(a)
  if (mod(n, 2) == 1) then
     median = a((n + 1) / 2)
  else
     median = 0.5_wp * (a(n / 2) + a(n / 2 + 1))
  end if

end function median

real(SCALAR_KIND) function meanTrimmed(a)

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(inout) :: a(:)

  ! <<< Scalar variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, n
  real(wp) :: firstQuartile, thirdQuartile

  interface
     real(SCALAR_KIND) function median(a)
       real(SCALAR_KIND), intent(in) :: a(:)
     end function median
  end interface

  n = size(a)

  if (mod(n, 2) == 0) then
     firstQuartile = median(a(1:n/2))
     thirdQuartile = median(a(n/2+1:n))
  else
     firstQuartile = median(a(1:(n-1)/2))
     thirdQuartile = median(a((n+1)/2+1:n))
  end if

  meanTrimmed = 0.0_wp
  n = 0

  do i = 1, size(a)
     if (a(i) >= firstQuartile .and. a(i) <= thirdQuartile) then
        meanTrimmed = meanTrimmed + a(i)
        n = n + 1
     end if
  end do

  if (n == 0) return
  meanTrimmed = meanTrimmed / real(n, wp)

end function meanTrimmed
