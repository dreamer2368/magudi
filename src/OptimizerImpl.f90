#include "config.h"

module OptimizerImpl

  implicit none
  public

contains

  subroutine mnbrakConjugateGradient(this, region)

    ! <<< External modules >>>
    use iso_fortran_env, only : output_unit

    ! <<< Derived types >>>
    use Optimizer_mod, only : t_Optimizer
    use Region_mod, only : t_Region
    use RegionVector_mod

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE
    use Region_enum, only : FORWARD

    ! <<< Private members >>>
    use RegionVectorImpl
    use TravelingWaveImpl

    ! <<< Internal modules >>>
    use SolverImpl, only : checkSolutionLimits
    use RegionImpl, only : computeSpatialAverages
    use ErrorHandler, only : writeAndFlush, gracefulExit
    use MPITimingsHelper, only : startTiming, endTiming

    implicit none

    ! <<< Arguments >>>
    class(t_Optimizer) :: this
    class(t_Region) :: region

    ! <<< Local variables >>>
    integer, parameter :: nIter = 1e4, wp = SCALAR_KIND
    logical :: bracketed(3)
    integer :: i, n
    SCALAR_TYPE :: step, stepEvaluation
    character(len = STRING_LENGTH) :: message

    call startTiming("mnbrak")

    bracketed = (/ .true., .false., .false. /)

    step = this%initialStep
    do n = 1, nIter
      call loadRegionVector(region,this%base - this%conjGrad * step, QOI_FORWARD_STATE)
      do i = 1, size(region%states) !... update state
         call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
              region%solverOptions)
      end do

      ! Check if physical quantities are within allowed limits.
      if (region%simulationFlags%enableSolutionLimits) then
        if ( checkSolutionLimits(region, FORWARD, this%outputPrefix) ) then
          write(message,'(A)') 'Retrying with smaller step.'
          call writeAndFlush(region%comm, output_unit, message)
          step = step / this%goldenRatio
          cycle
        end if
      end if

      call computeTravelingWaveResidual(region,this%residual)
      stepEvaluation = 0.5_wp * regionInnerProduct(this%residual,this%residual,region)

      if ( stepEvaluation<this%bracket(1,2) ) then
        this%bracket(2,:) = (/ step, stepEvaluation /)
        bracketed(2) = .true.
        exit
      else
        this%bracket(3,:) = (/ step, stepEvaluation /)
        bracketed(3) = .true.
        step = step / this%goldenRatio
      end if
    end do

    if (all(bracketed)) then
      if (this%verbose) then
        call this%printBracket()
        write(message,'(A)') "Finding minimum bracket succeeded."
        call writeAndFlush(region%comm, output_unit, message)
      end if
      call endTiming("mnbrak")
      return
    elseif (.not. bracketed(2)) then
      call this%printBracket()
      write(message,'(A)') "Finding minimum bracket failed."
      call gracefulExit(region%comm,message)
    end if

    step = this%bracket(2,1) * this%goldenRatio
    do n = 1, nIter
      call loadRegionVector(region,this%base - this%conjGrad * step, QOI_FORWARD_STATE)
      do i = 1, size(region%states) !... update state
         call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
              region%solverOptions)
      end do
      call computeTravelingWaveResidual(region,this%residual)
      stepEvaluation = 0.5_wp * regionInnerProduct(this%residual,this%residual,region)

      if ( stepEvaluation .le. this%bracket(2,2) ) then
        this%bracket(1,:) = this%bracket(2,:)
        this%bracket(2,:) = (/ step, stepEvaluation /)
        step = step * this%goldenRatio
      else
        this%bracket(3,:) = (/ step, stepEvaluation /)
        bracketed(3) = .true.
        exit
      end if
    end do

    if (all(bracketed)) then
      if (this%verbose) then
        call this%printBracket()
        write(message,'(A)') "Finding minimum bracket succeeded."
        call writeAndFlush(region%comm, output_unit, message)
      end if
    else
      call this%printBracket()
      write(message,'(A)') "Finding minimum bracket failed."
      call gracefulExit(region%comm,message)
    end if

    call endTiming("mnbrak")

  end subroutine mnbrakConjugateGradient

  subroutine linminConjugateGradient(this, region)

    ! <<< External modules >>>
    use iso_fortran_env, only : output_unit

    ! <<< Derived types >>>
    use Optimizer_mod, only : t_Optimizer
    use Region_mod, only : t_Region
    use RegionVector_mod

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE
    use Region_enum, only : FORWARD

    ! <<< Private members >>>
    use RegionVectorImpl
    use TravelingWaveImpl

    ! <<< Internal modules >>>
    use RegionImpl, only : computeSpatialAverages
    use ErrorHandler, only : writeAndFlush, gracefulExit
    use MPITimingsHelper, only : startTiming, endTiming

    implicit none

    ! <<< Arguments >>>
    class(t_Optimizer) :: this
    class(t_Region) :: region

    ! <<< Local variables >>>
    integer, parameter :: nIter = 1e4, wp = SCALAR_KIND
    SCALAR_TYPE :: Cr, eps
    integer :: i, n, index
    SCALAR_TYPE :: x(3), fx(3), newBracket(4,2), step, stepEvaluation
    character(len = STRING_LENGTH) :: message
    logical :: success

    call startTiming('linmin')

    Cr = 1.0_wp - 1.0_wp / this%goldenRatio
    eps = 1.0E-15

    success = .false.

    do n = 1, nIter
      x = this%bracket(:,1)
      fx = this%bracket(:,2)

      success = ( (x(3)-x(1)) < x(2)*this%linminTol+eps )
      if ( success ) then
        if (this%verbose) then
          call this%printBracket()
          write(message,'(A)') "Line minimization succeeded."
          call writeAndFlush(region%comm, output_unit, message)
        end if

        this%initialStep = x(2)
        call loadRegionVector(region,this%base - this%conjGrad * x(2), QOI_FORWARD_STATE)
        call endTiming('linmin')
        return
      end if

      step = x(2) - 0.5_wp * ( (x(2)-x(1))**2 * (fx(2)-fx(3)) - (x(2)-x(3))**2 * (fx(2)-fx(1)) )      &
                           / ( (x(2)-x(1)) * (fx(2)-fx(3)) - (x(2)-x(3)) * (fx(2)-fx(1)) )
      if ( (step>x(3)) .or. (step<x(1)) .or.                                    &
           (abs(log10((x(3)-x(2))/(x(2)-x(1))))>1.0_wp) ) then
        if ( x(2)>0.5_wp*(x(1)+x(3)) ) then
          step = x(2) - Cr * (x(2) - x(1))
        else
          step = x(2) + Cr * (x(3) - x(2))
        end if
      end if

      call loadRegionVector(region,this%base - this%conjGrad * step, QOI_FORWARD_STATE)
      do i = 1, size(region%states) !... update state
         call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
              region%solverOptions)
      end do
      call computeTravelingWaveResidual(region,this%residual)
      stepEvaluation = 0.5_wp * regionInnerProduct(this%residual,this%residual,region)

      newBracket(1,:) = this%bracket(1,:)
      newBracket(4,:) = this%bracket(3,:)
      if (step>x(2)) then
        newBracket(2,:) = this%bracket(2,:)
        newBracket(3,:) = (/ step, stepEvaluation /)
      else
        newBracket(2,:) = (/ step, stepEvaluation /)
        newBracket(3,:) = this%bracket(2,:)
      end if

      index = minloc(newBracket(:,2), dim=1)
      assert_key(index, (2, 3))
      this%bracket = newBracket(index-1:index+1,:)
    end do

    x = this%bracket(:,1)
    success = ( (x(3)-x(1)) < x(2)*this%linminTol+eps )
    if ( success ) then
      if (this%verbose) then
        call this%printBracket()
        write(message,'(A)') "Line minimization succeeded."
        call writeAndFlush(region%comm, output_unit, message)
      end if
    else
      call this%printBracket()
      write(message,'(A)') "Line minimization failed."
      call gracefulExit(region%comm,message)
    end if

    this%initialStep = x(2)
    call loadRegionVector(region,this%base - this%conjGrad * x(2), QOI_FORWARD_STATE)

    call endTiming('linmin')

  end subroutine linminConjugateGradient

  subroutine frprmnConjugateGradient(this, region, gg, initial)

    ! <<< External modules >>>
    use iso_fortran_env, only : output_unit

    ! <<< Derived types >>>
    use Optimizer_mod, only : t_Optimizer
    use Region_mod, only : t_Region
    use RegionVector_mod

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE
    use Region_enum, only : FORWARD

    ! <<< Private members >>>
    use RegionVectorImpl

    ! <<< Internal modules >>>
    use RegionImpl, only : computeSpatialAverages
    use ErrorHandler, only : writeAndFlush, gracefulExit
    use MPITimingsHelper, only : startTiming, endTiming

    implicit none

    ! <<< Arguments >>>
    class(t_Optimizer) :: this
    class(t_Region) :: region
    real(SCALAR_KIND), intent(in) :: gg
    logical, intent(in) :: initial

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(SCALAR_KIND) :: dgg, beta, betaFR

    call startTiming('frprmn')

    call saveRegionVector(this%base,region,QOI_FORWARD_STATE)

    if (initial) then
      this%conjGrad = this%grad

      this%prevGrad = this%grad
      this%gg1 = gg

      call endTiming('frprmn')
      return
    end if

    this%prevCG = this%conjGrad
    dgg = regionInnerProduct(this%grad-this%prevGrad,this%grad,region)
    beta = dgg / this%gg1
    betaFR = gg / this%gg1
    if (beta > betaFR) then
      beta = betaFR
    elseif (beta < -betaFR) then
      beta = - betaFR
    end if
    this%conjGrad = this%prevCG * beta + this%grad

    this%prevGrad = this%grad
    this%gg1 = gg

    call endTiming('frprmn')

  end subroutine frprmnConjugateGradient

end module OptimizerImpl

subroutine setupOptimizer(this, region)

  ! <<< Derived types >>>
  use Optimizer_mod, only : t_Optimizer
  use Region_mod, only : t_Region
  use Functional_mod, only : t_Functional

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, bufferSize
  class(t_Functional), pointer :: functional => null()

  call this%cleanup()

  this%goldenRatio = 0.5_wp * ( 1.0_wp + sqrt(5.0_wp) )
  this%initialStep = getOption("optimization/initial_step",0.1_wp)
  this%linminTol = getOption("optimization/linmin_tolerance",(0.1_wp)**3)
  this%cgTol = getOption("optimization/cg_tolerance",(0.1_wp)**5)
  this%saveInterval = getOption("save_interval",-1)
  this%reportInterval = getOption("report_interval",-1)
  this%nTimesteps = getOption("number_of_timesteps",100)
  this%nTimesteps = max(0, this%nTimesteps)
  this%outputPrefix = getOption("output_prefix", PROJECT_NAME)
  this%verbose = getOption("optimization/verbosity",.false.)

  this%numGrids = size(region%grids)
  this%numParams = 0
  if (allocated(region%params%buffer)) this%numParams = size(region%params%buffer,1)

  call this%residual%set(region)
  call this%base%set(region)
  call this%grad%set(region)
  call this%conjGrad%set(region)
  call this%prevGrad%set(region)
  call this%prevCG%set(region)

end subroutine setupOptimizer

subroutine cleanupOptimizer(this)

  ! <<< Derived types >>>
  use Optimizer_mod, only : t_Optimizer
  use RegionVector_mod, only : t_RegionVector

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this

  call this%residual%cleanup()
  call this%base%cleanup()
  call this%grad%cleanup()
  call this%conjGrad%cleanup()
  call this%prevGrad%cleanup()
  call this%prevCG%cleanup()

end subroutine cleanupOptimizer

subroutine printBracket(this)

  ! <<< External modules >>>
  use MPI
  use iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Optimizer_mod, only : t_Optimizer

  ! <<< Internal modules >>>
  use ErrorHandler, only : writeAndFlush

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this

  ! <<< Local variables >>>
  character(len = STRING_LENGTH) :: message

  write(message,'(A,3(1X,'// SCALAR_FORMAT //'))') '  x : ',                    &
                      this%bracket(1,1), this%bracket(2,1), this%bracket(3,1)
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  write(message,'(A,3(1X,'// SCALAR_FORMAT //'))') 'f(x): ',                    &
                      this%bracket(1,2), this%bracket(2,2), this%bracket(3,2)
  call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

end subroutine printBracket

subroutine showProgressOptimizer(this, region, costFunctional, costSensitivity, &
                                 outputFilename, append)

  ! <<< External modules >>>
  use MPI
  use iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Optimizer_mod, only : t_Optimizer

  ! <<< Internal modules >>>
  use TravelingWaveImpl
  use InputHelper, only : getOption, getFreeUnit
  use ErrorHandler, only : writeAndFlush, gracefulExit

  implicit none

  class(t_Optimizer) :: this
  class(t_Region) :: region
  real(SCALAR_KIND), intent(in) :: costFunctional, costSensitivity
  character(len=*), intent(in) :: outputFilename
  logical, intent(in), optional :: append

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: str, filename, message
  logical :: fileExists, append_
  integer :: procRank, ierror, fileUnit, ostat

  append_ = .false.
  if (present(append)) append_ = append

  if (this%reportInterval > 0 .and. mod(region%timestep, max(1, this%reportInterval)) == 0) then
    write(str, '(2A,I8,3(A,E13.6))') PROJECT_NAME, ": timestep = ", region%timestep,          &
         ", convection speed = ", region%params%buffer(1,1),                                  &
         ", cost functional = ", costFunctional, ", cost sensitivity = ", costSensitivity
    call writeAndFlush(region%comm, output_unit, str)

    call MPI_Comm_rank(region%comm, procRank, ierror)
    if (procRank == 0) then
       if (.not. append_) then
          open(unit = getFreeUnit(fileUnit), file = trim(outputFilename), action = 'write',          &
               status = 'unknown', iostat = ostat)
       else
          open(unit = getFreeUnit(fileUnit), file = trim(outputFilename), action = 'write',          &
               status = 'old', position = 'append', iostat = ostat)
       end if
    end if

    call MPI_Bcast(ostat, 1, MPI_INTEGER, 0, region%comm, ierror)
    if (ostat /= 0) then
      write(message, "(2A)") trim(filename), ": Failed to open file for writing!"
      call gracefulExit(region%comm, message)
    end if

    if (procRank == 0) then
       write(fileUnit, '(I8,3(1X,SP,' // SCALAR_FORMAT // '))')                       &
          region%timestep, region%params%buffer(1,1), costFunctional, costSensitivity
    end if

    call MPI_Bcast(ostat, 1, MPI_INTEGER, 0, region%comm, ierror)
    if (ostat /= 0) then
       write(message, "(2A)") trim(filename), ": Error writing to file!"
       call gracefulExit(region%comm, message)
    end if

    if (procRank == 0) then
       flush(fileUnit)
       close(fileUnit)
    end if

    call MPI_Barrier(region%comm, ierror)
  end if

  if (this%saveInterval > 0 .and. mod(region%timestep, max(1, this%saveInterval)) == 0) then
    write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", region%timestep, ".q"
    call saveTravelingWave(region, filename)
  end if

end subroutine showProgressOptimizer

subroutine verifyAdjoint(this, region)

  ! <<< External modules >>>
  use iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use Optimizer_mod, only : t_Optimizer
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use RegionVector_mod

  ! <<< Enumerations >>>
  use State_enum, only : QOI_FORWARD_STATE, QOI_ADJOINT_STATE, QOI_RIGHT_HAND_SIDE
  use Region_enum, only : FORWARD, ADJOINT, LINEARIZED

  ! <<< Private members >>>
  use SolverImpl, only : showProgress, checkSolutionLimits, loadInitialCondition
  use RegionImpl, only : computeRegionIntegral
  use RegionVectorImpl
  use TravelingWaveImpl

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming
  use ErrorHandler, only : writeAndFlush
  use InputHelper, only : getOption, getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: filename, message
  integer :: i, j, k, nDimensions, nUnknowns
  SCALAR_TYPE :: costFunctional0, costFunctional1, costSensitivity, ATxy, xAy
  SCALAR_TYPE :: stepSizes(32), errorHistory(32), convergenceHistory(31)
  logical :: solutionCrashes = .false., success_

  call startTiming("verifyAdjoint")

  ! Load the initial condition.
  call getRequiredOption("initial_condition_file", filename, region%comm)
  call region%loadData(QOI_FORWARD_STATE, filename) !... initialize from file.

  if (region%simulationFlags%enableBodyForce .or. region%simulationFlags%checkConservation) then
    region%oneOverVolume = computeRegionIntegral(region)
    region%oneOverVolume = 1.0_wp / region%oneOverVolume

    region%momentumLossPerVolume = 0.0_wp
  end if

  do i = 1, size(region%states) !... update state
     call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
          region%solverOptions)
  end do
  call saveRegionVector(this%base,region,QOI_FORWARD_STATE)

  call computeTravelingWaveResidual(region,this%residual)
  costFunctional0 = 0.5_wp * regionInnerProduct(this%residual,this%residual,region)
  write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))')                                      &
                              'Forward run: cost functional = ', costFunctional0
  call writeAndFlush(region%comm, output_unit, message)

  this%grad = computeTravelingWaveJacobian(region,this%residual,ADJOINT)
  costSensitivity = regionInnerProduct(this%grad,this%grad,region)
  write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))')                                      &
                             'Adjoint run: cost sensitivity = ', costSensitivity
  call writeAndFlush(region%comm, output_unit, message)

  call random_number(this%prevGrad%params)
  do i = 1, size(region%states)
    call random_number(this%prevGrad%states(i)%conservedVariables)
  end do

  ATxy = regionInnerProduct(this%grad,this%prevGrad,region)
  xAy = regionInnerProduct(this%residual,                                       &
           computeTravelingWaveJacobian(region,this%prevGrad,LINEARIZED),region)
  write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))') '< AT * x, y > = ', ATxy
  call writeAndFlush(region%comm, output_unit, message)
  write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))') '< x, A * y > = ', xAy
  call writeAndFlush(region%comm, output_unit, message)

  stepSizes(1) = 0.00001_wp
  do k = 2, size(stepSizes)
     stepSizes(k) = stepSizes(k-1) * 10.0_wp**(-0.25_wp)
  end do
  errorHistory = 0.0_wp
  do k = 1, size(stepSizes)
    call loadRegionVector(region,this%base + this%grad*stepSizes(k),QOI_FORWARD_STATE)
    do i = 1, size(region%states) !... update state
       call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
            region%solverOptions)
    end do
    call region%computeRhs(FORWARD,1,1)
    call computeTravelingWaveResidual(region,this%residual)
    costFunctional1 = 0.5_wp * regionInnerProduct(this%residual,this%residual,region)

    errorHistory(k) = abs( ((costFunctional1-costFunctional0)/stepSizes(k) - costSensitivity)/costSensitivity )
    write(message, '(4(' // SCALAR_FORMAT // ',1X))') stepSizes(k), costSensitivity,                   &
                  (costFunctional1-costFunctional0)/stepSizes(k), errorHistory(k)
    call writeAndFlush(region%comm, output_unit, message)

    if (k > 1) then
       convergenceHistory(k-1) = log(errorHistory(k) / errorHistory(k-1)) /              &
            log(stepSizes(k) / stepSizes(k-1))
       if (k > 5) then
           if (sum(convergenceHistory(k-3:k-1))/3.0_wp < 0.0_wp) exit
       end if
    end if
  end do

  call endTiming("verifyAdjoint")

end subroutine verifyAdjoint

subroutine runNLCG(this, region, restartFilename)

  ! <<< External modules >>>
  use MPI
  use iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Optimizer_mod, only : t_Optimizer
  use RegionVector_mod

  ! <<< Enumerations >>>
  use State_enum, only : QOI_FORWARD_STATE, QOI_ADJOINT_STATE, QOI_RIGHT_HAND_SIDE
  use Region_enum, only : FORWARD, ADJOINT, LINEARIZED

  ! <<< Internal modules >>>
  use OptimizerImpl
  use RegionVectorImpl
  use TravelingWaveImpl
  use SolverImpl, only : checkSolutionLimits
  use RegionImpl, only : computeRegionIntegral
  use MPITimingsHelper, only : startTiming, endTiming
  use ErrorHandler, only : writeAndFlush, gracefulExit
  use InputHelper, only : getOption, getRequiredOption, getFreeUnit

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this
  class(t_Region) :: region
  character(len = *), intent(in), optional :: restartFilename

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, fileUnit, ostat, procRank, ierror, startTimestep, timestep
  SCALAR_TYPE :: costFunctional0, costFunctional1, costSensitivity
  character(len = STRING_LENGTH) :: filename, outputFilename, message

  call startTiming("runNLCG")

  write(outputFilename, '(2A)') trim(this%outputPrefix), ".traveling_wave.nlcg.txt"

  ! Load the initial condition.
  if (present(restartFilename)) then
    call loadTravelingWave(region,trim(restartFilename))
  else
    call getRequiredOption("initial_condition_file", filename, region%comm)
    call region%loadData(QOI_FORWARD_STATE, filename) !... initialize from file.
  end if
  startTimestep = region%timestep

  if (region%simulationFlags%enableBodyForce .or. region%simulationFlags%checkConservation) then
    region%oneOverVolume = computeRegionIntegral(region)
    region%oneOverVolume = 1.0_wp / region%oneOverVolume

    region%momentumLossPerVolume = 0.0_wp
  end if

  do timestep = startTimestep, startTimestep + this%nTimesteps ! takes one more step

    region%timestep = timestep
    do i = 1, size(region%states) !... update state
       call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
            region%solverOptions)
    end do

    ! Check if physical quantities are within allowed limits.
    if (region%simulationFlags%enableSolutionLimits) then
      if ( checkSolutionLimits(region, FORWARD, this%outputPrefix) ) then
        write(message,'(A)') 'Optimization failed.'
        call gracefulExit(region%comm, message)
      end if
    end if

    call computeTravelingWaveResidual(region,this%residual)
    costFunctional0 = 0.5_wp * regionInnerProduct(this%residual,this%residual,region)
    if (this%verbose) then
      write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))')                                    &
                                  'Forward run: cost functional = ', costFunctional0
      call writeAndFlush(region%comm, output_unit, message)
    end if

    this%grad = computeTravelingWaveJacobian(region,this%residual,ADJOINT)
    costSensitivity = regionInnerProduct(this%grad,this%grad,region)
    if (this%verbose) then
      write(message, '(A,(1X,SP,' // SCALAR_FORMAT // '))')                                    &
                                 'Adjoint run: cost sensitivity = ', costSensitivity
      call writeAndFlush(region%comm, output_unit, message)
    end if

    call this%showProgress(region, costFunctional0, costSensitivity,                           &
                           outputFilename, timestep - startTimestep > this%reportInterval)
    if ( costSensitivity < this%cgTol ) then
      write(message, '(A)') 'Gradient is below tolerance. Finishing the optimization.'
      call writeAndFlush(region%comm, output_unit, message)
      call endTiming("runNLCG")
      return
    end if

    call frprmnConjugateGradient(this, region, costSensitivity, timestep == startTimestep)

    this%bracket(1,:) = (/ 0.0_wp, costFunctional0 /)
    call mnbrakConjugateGradient(this, region)

    call linminConjugateGradient(this, region)

  end do

  call endTiming("runNLCG")

end subroutine
