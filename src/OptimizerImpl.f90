#include "config.h"

module OptimizerImpl

  implicit none
  public

contains

  subroutine mnbrakNLCG(this, region)

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

  end subroutine mnbrakNLCG

  subroutine linminNLCG(this, region)

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

  end subroutine linminNLCG

  subroutine frprmnNLCG(this, region, gg, initial)

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

  end subroutine frprmnNLCG

  function arnoldi(this, region, substep) result(h)

    ! <<< Derived types >>>
    use Optimizer_mod, only : t_Optimizer
    use Region_mod, only : t_Region
    use RegionVector_mod

    ! <<< Enumerations >>>
    use Region_enum, only : LINEARIZED

    ! <<< Private members >>>
    use RegionVectorImpl, only : regionInnerProduct
    use TravelingWaveImpl, only : computeTravelingWaveJacobian
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    class(t_Optimizer) :: this
    class(t_Region) :: region
    integer, intent(in) :: substep
    SCALAR_TYPE :: h(substep+1)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i

    call startTiming("arnoldi")

    assert((substep>=1) .and. (substep<=this%maxSubsteps))

    this%Q(substep+1) = computeTravelingWaveJacobian(region,this%Q(substep),LINEARIZED)

    do i = 1, substep
      h(i) = regionInnerProduct(this%Q(substep+1),this%Q(i),region)
      this%Q(substep+1) = this%Q(substep+1) - this%Q(i) * h(i)
    end do
    h(substep+1) = sqrt(regionInnerProduct(this%Q(substep+1),this%Q(substep+1),region))
    this%Q(substep+1) = this%Q(substep+1) * (1.0_wp/h(substep+1))
    call endTiming("arnoldi")

  end function arnoldi

  subroutine givensRotation(h, cs, sn, substep)

    ! <<< Private members >>>
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    SCALAR_TYPE, intent(inout) :: h(:), cs(:), sn(:)
    integer, intent(in) :: substep

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i
    SCALAR_TYPE :: rotation(2,2), temp

    call startTiming("givensRotation")

    assert(size(h)==substep+1)
    assert(size(cs)==substep)
    assert(size(sn)==substep)

    do i = 1, substep - 1
      rotation(1,:) = (/ cs(i), sn(i) /)
      rotation(2,:) = (/ -rotation(1,2), rotation(1,1) /)
      h(i:i+1) = matmul(rotation,h(i:i+1))
    end do

    ! compute Givens rotation matrix.
    if (h(substep)==0.0_wp) then
      cs(substep) = 0.0_wp
      sn(substep) = 1.0_wp
    else
      temp = sqrt(sum(h(substep:substep+1)**2))
      cs(substep) = abs(h(substep)) / temp
      sn(substep) = cs(substep) * h(substep+1) / h(substep)
    end if

    ! eliminate the last h.
    h(substep) = cs(substep) * h(substep) + sn(substep) * h(substep+1)
    h(substep+1) = 0.0_wp

    call endTiming("givensRotation")

  end subroutine givensRotation

  pure function backwardSubstitution(A,b) result(x)

    ! <<< Arguments >>>
    SCALAR_TYPE, intent(in) :: A(:,:), b(:)
    SCALAR_TYPE :: x(size(b))

    ! <<< Local variables >>>
    integer :: i, N

    assert(size(A,1)==size(A,2))
    assert(size(A,1)==size(b))

    N = size(b)

    x(N) = b(N) / A(N,N)
    do i = N-1, 1, -1
      x(i) = ( b(i) - sum(A(i,i+1:N) * x(i+1:N)) ) / A(i,i)
    end do
  end function backwardSubstitution

end module OptimizerImpl

subroutine setupOptimizer(this, region, mode)

  ! <<< Derived types >>>
  use Optimizer_mod, only : t_Optimizer
  use Region_mod, only : t_Region
  use Functional_mod, only : t_Functional

  ! <<< Enumerations >>>
  use Optimizer_enum

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this
  class(t_Region) :: region
  integer, intent(in) :: mode

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, bufferSize
  class(t_Functional), pointer :: functional => null()

  call this%cleanup()

  this%goldenRatio = 0.5_wp * ( 1.0_wp + sqrt(5.0_wp) )
  this%initialStep = getOption("optimization/initial_step",0.1_wp)
  this%linminTol = getOption("optimization/linmin_tolerance",(0.1_wp)**3)
  this%cgTol = getOption("optimization/cg_tolerance",(0.1_wp)**5)
  this%maxSubsteps = getOption("optimization/max_substeps",30)
  this%maxRestart = getOption("optimization/max_gmres_restart",100)
  this%saveInterval = getOption("save_interval",-1)
  this%reportInterval = getOption("report_interval",-1)
  this%nTimesteps = getOption("number_of_timesteps",100)
  this%nTimesteps = max(0, this%nTimesteps)
  this%outputPrefix = getOption("output_prefix", PROJECT_NAME)
  this%verbose = getOption("optimization/verbosity",.false.)

  this%numGrids = size(region%grids)
  this%numParams = 0
  if (allocated(region%params%buffer)) this%numParams = size(region%params%buffer,1)

  select case(mode)
  case(VERIFY)
    call this%residual%set(region)
    call this%base%set(region)
    call this%grad%set(region)
    call this%prevGrad%set(region)
  case(NLCG)
    call this%residual%set(region)
    call this%base%set(region)
    call this%grad%set(region)
    call this%conjGrad%set(region)
    call this%prevGrad%set(region)
    call this%prevCG%set(region)
  case(CGS)
    call this%residual%set(region)
    call this%base%set(region)
    call this%x%set(region)
    call this%grad%set(region)
    call this%conjGrad%set(region)
    call this%Ap%set(region)
  case(GMRES)
    call this%residual%set(region)
    call this%base%set(region)
    call this%x%set(region)
    call this%r%set(region)

    assert(this%maxSubsteps>0)
    allocate(this%H(this%maxSubsteps+1,this%maxSubsteps))
    allocate(this%sn(this%maxSubsteps))
    allocate(this%cs(this%maxSubsteps))
    allocate(this%beta(this%maxSubsteps+1))
    allocate(this%Q(this%maxSubsteps+1))
    do i = 1, this%maxSubsteps+1
      call this%Q(i)%set(region)
    end do
  case(BICGSTAB)
    call this%residual%set(region)
    call this%base%set(region)
    call this%x%set(region)
    call this%r%set(region)
    call this%rcg%set(region)
    call this%conjGrad%set(region)
    call this%Ap%set(region)
    call this%s%set(region)
    call this%As%set(region)
  end select

end subroutine setupOptimizer

subroutine cleanupOptimizer(this)

  ! <<< Derived types >>>
  use Optimizer_mod, only : t_Optimizer
  use RegionVector_mod, only : t_RegionVector

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this

  ! <<< Local variables >>>
  integer :: i

  call this%residual%cleanup()
  call this%base%cleanup()
  call this%grad%cleanup()
  call this%conjGrad%cleanup()
  call this%prevGrad%cleanup()
  call this%prevCG%cleanup()

  call this%x%cleanup()
  call this%r%cleanup()
  call this%rcg%cleanup()
  call this%Ap%cleanup()
  call this%s%cleanup()
  call this%As%cleanup()

  SAFE_DEALLOCATE(this%H)
  SAFE_DEALLOCATE(this%sn)
  SAFE_DEALLOCATE(this%cs)
  SAFE_DEALLOCATE(this%beta)
  if (allocated(this%Q)) then
    do i = 1, this%maxSubsteps+1
      call this%Q(i)%cleanup()
    end do
    SAFE_DEALLOCATE(this%Q)
  end if

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

subroutine showProgressOptimizer(this, region, mode, step, scalars,             &
                                 outputFilename, append)

  ! <<< External modules >>>
  use MPI
  use iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Optimizer_mod, only : t_Optimizer

  ! <<< Enumerations >>>
  use Optimizer_enum

  ! <<< Internal modules >>>
  use TravelingWaveImpl
  use InputHelper, only : getOption, getFreeUnit
  use ErrorHandler, only : writeAndFlush, gracefulExit

  implicit none

  class(t_Optimizer) :: this
  class(t_Region) :: region
  integer, intent(in) :: mode, step
  real(SCALAR_KIND), intent(in) :: scalars(:)
  character(len=*), intent(in) :: outputFilename
  logical, intent(in), optional :: append

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: str, filename, message
  logical :: fileExists, append_
  integer :: procRank, ierror, fileUnit, ostat

  append_ = .false.
  if (present(append)) append_ = append

  if (this%reportInterval > 0 .and. mod(step, max(1, this%reportInterval)) == 0) then
    select case(mode)
    case (VERIFY, NLCG)
      assert(size(scalars)==2)
      write(str, '(2A,I8,3(A,E13.6))') PROJECT_NAME, ": timestep = ", step,                     &
           ", convection speed = ", region%params%buffer(1,1),                                  &
           ", cost functional = ", scalars(1), ", cost sensitivity = ", scalars(2)
    case (NEWTON, CGS, GMRES, BICGSTAB)
      assert(size(scalars)>=1)
      write(str, '(2A,I8,2(A,E13.6))') PROJECT_NAME, ": timestep = ", step,                     &
           ", convection speed = ", region%params%buffer(1,1),                                  &
           ", residual = ", scalars(1)
    case default
      write(str, '(A)') 'The mode cannot be identified!'
      call gracefulExit(region%comm, str)
    end select
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
      select case(mode)
      case (VERIFY, NLCG)
        write(fileUnit, '(I8,3(1X,SP,' // SCALAR_FORMAT // '))')                       &
              step, region%params%buffer(1,1), scalars(1), scalars(2)
      case (NEWTON, GMRES, BICGSTAB, CGS)
        write(fileUnit, '(I8,2(1X,SP,' // SCALAR_FORMAT // '))')                       &
              step, region%params%buffer(1,1), scalars(1)
      end select
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

  select case(mode)
  case (VERIFY, NLCG, NEWTON)
    if (this%saveInterval > 0 .and. mod(step, max(1, this%saveInterval)) == 0) then
      write(filename, '(2A,I8.8,A)') trim(this%outputPrefix), "-", step, ".q"
      call saveTravelingWave(region, filename)
    end if
  end select

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
  use Optimizer_enum, only : NLCG

  ! <<< Internal modules >>>
  use OptimizerImpl, only : mnbrakNLCG, linminNLCG, frprmnNLCG
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

    call this%showProgress(region, NLCG, timestep, (/ costFunctional0, costSensitivity /),     &
                           outputFilename, timestep - startTimestep > this%reportInterval)
    if ( costSensitivity < this%cgTol ) then
      write(message, '(A)') 'Gradient is below tolerance. Finishing the optimization.'
      call writeAndFlush(region%comm, output_unit, message)
      call endTiming("runNLCG")
      return
    end if

    call frprmnNLCG(this, region, costSensitivity, timestep == startTimestep)

    this%bracket(1,:) = (/ 0.0_wp, costFunctional0 /)
    call mnbrakNLCG(this, region)

    call linminNLCG(this, region)

  end do

  call endTiming("runNLCG")

end subroutine

subroutine runCGS(this, region, restartFilename)

  ! <<< External modules >>>
  use MPI
  use iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Optimizer_mod, only : t_Optimizer
  use RegionVector_mod

  ! <<< Enumerations >>>
  use State_enum, only : QOI_FORWARD_STATE
  use Region_enum, only : FORWARD, LINEARIZED, ADJOINT
  use Optimizer_enum, only : NEWTON, CGS

  ! <<< Internal modules >>>
  use RegionVectorImpl
  use TravelingWaveImpl
  use SolverImpl, only : checkSolutionLimits
  use RegionImpl, only : computeRegionIntegral
  use MPITimingsHelper, only : startTiming, endTiming
  use ErrorHandler, only : writeAndFlush, gracefulExit
  use InputHelper, only : getFreeUnit

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this
  class(t_Region) :: region
  character(len = *), intent(in) :: restartFilename

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, m, n, fileUnit, ostat, procRank, ierror,                        &
             startTimestep, timestep, substep, numSteps, restart,                        &
             totalSubsteps
  SCALAR_TYPE :: residualNorm, residualNorm1, rr, rr1, rr0, alpha, beta
  character(len = STRING_LENGTH) :: strFormat, filename, outputFilename, message

  call startTiming("runCGS")

  write(outputFilename, '(2A)') trim(this%outputPrefix), ".traveling_wave.cgs.txt"
  write(filename, '(2A)') trim(this%outputPrefix), ".traveling_wave.newton.txt"

  ! Load the initial condition.
  call loadTravelingWave(region,trim(restartFilename))
  startTimestep = region%timestep
  totalSubsteps = 0

  if (region%simulationFlags%enableBodyForce .or. region%simulationFlags%checkConservation) then
    region%oneOverVolume = computeRegionIntegral(region)
    region%oneOverVolume = 1.0_wp / region%oneOverVolume

    region%momentumLossPerVolume = 0.0_wp
  end if

  do i = 1, size(region%states) !... update state
     call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
          region%solverOptions)
  end do

  ! Check if physical quantities are within allowed limits.
  if (region%simulationFlags%enableSolutionLimits) then
    if ( checkSolutionLimits(region, FORWARD, this%outputPrefix) ) then
      write(message,'(A)') 'Newton iteration reached an infeasible solution!'
      call gracefulExit(region%comm, message)
    end if
  end if

  call computeTravelingWaveResidual(region,this%residual)
  residualNorm = sqrt(regionInnerProduct(this%residual,this%residual,region))

  do timestep = startTimestep + 1, startTimestep + this%nTimesteps

    region%timestep = timestep

    call saveRegionVector(this%base,region,QOI_FORWARD_STATE)
    this%x%params = 0.0_wp
    do i = 1, size(region%states)
      this%x%states(i)%conservedVariables = 0.0_wp
    end do
    this%grad = computeTravelingWaveJacobian(region,this%residual * (-1.0_wp),ADJOINT)
    this%conjGrad = this%grad
    rr = regionInnerProduct(this%grad,this%grad,region)
    rr0 = rr

    do substep = 1, this%maxSubsteps
      rr1 = rr

      this%Ap = computeTravelingWaveJacobian(region,this%conjGrad,LINEARIZED)
      this%Ap = computeTravelingWaveJacobian(region,this%Ap,ADJOINT)
      alpha = rr1 / regionInnerProduct(this%conjGrad,this%Ap,region)

      this%x = this%x + this%conjGrad * alpha
      this%grad = this%grad - this%Ap * alpha

      rr = regionInnerProduct(this%grad,this%grad,region)
      beta = rr / rr1
      this%conjGrad = this%grad + this%conjGrad * beta

      call this%showProgress(region, CGS, totalSubsteps + substep, (/ sqrt(rr) /), &
               outputFilename, totalSubsteps + substep - 1 > this%reportInterval)
      if (sqrt(rr)<this%linminTol*sqrt(rr0)) then
        exit
      end if
    end do
    totalSubsteps = totalSubsteps + substep

    if (sqrt(rr)>=this%linminTol*sqrt(rr0)) then
      write(message,'(A)') 'CGS iteration failed.'
      call gracefulExit(region%comm, message)
    end if

    call loadRegionVector(region,this%base + this%x,QOI_FORWARD_STATE)

    do i = 1, size(region%states) !... update state
       call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
            region%solverOptions)
    end do

    ! Check if physical quantities are within allowed limits.
    if (region%simulationFlags%enableSolutionLimits) then
      if ( checkSolutionLimits(region, FORWARD, this%outputPrefix) ) then
        write(message,'(A)') 'Newton iteration reached an infeasible solution!'
        call gracefulExit(region%comm, message)
      end if
    end if

    call computeTravelingWaveResidual(region,this%residual)
    residualNorm1 = residualNorm
    residualNorm = sqrt(regionInnerProduct(this%residual,this%residual,region))

    call this%showProgress(region, NEWTON, timestep, (/ residualNorm /),                       &
                           filename, timestep - startTimestep - 1 > this%reportInterval)

    if (residualNorm > residualNorm1) then
      write(message,'(2(A,ES13.6))') 'Newton iteration failed! Previous: ', residualNorm1,     &
                                     ', Current: ', residualNorm
      call gracefulExit(region%comm, message)
    end if

  end do

  call endTiming("runCGS")

end subroutine runCGS

subroutine runGMRES(this, region, restartFilename)

  ! <<< External modules >>>
  use MPI
  use iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Optimizer_mod, only : t_Optimizer
  use RegionVector_mod

  ! <<< Enumerations >>>
  use State_enum, only : QOI_FORWARD_STATE
  use Region_enum, only : FORWARD, LINEARIZED
  use Optimizer_enum, only : NEWTON, GMRES

  ! <<< Internal modules >>>
  use OptimizerImpl, only : arnoldi, givensRotation, backwardSubstitution
  use RegionVectorImpl
  use TravelingWaveImpl
  use SolverImpl, only : checkSolutionLimits
  use RegionImpl, only : computeRegionIntegral
  use MPITimingsHelper, only : startTiming, endTiming
  use ErrorHandler, only : writeAndFlush, gracefulExit
  use InputHelper, only : getFreeUnit

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this
  class(t_Region) :: region
  character(len = *), intent(in) :: restartFilename

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, m, n, fileUnit, ostat, procRank, ierror,                        &
             startTimestep, timestep, substep, numSteps, restart,                        &
             totalSubsteps
             ! Ns, index, nUnknowns, nGridPoints, offset, offset1
  SCALAR_TYPE :: residualNorm, residualNorm1, rNorm, error
  character(len = STRING_LENGTH) :: strFormat, filename, outputFilename, message
  SCALAR_TYPE, allocatable :: y(:), A(:,:)

  call startTiming("runGMRES")

  write(outputFilename, '(2A)') trim(this%outputPrefix), ".traveling_wave.gmres.txt"
  write(filename, '(2A)') trim(this%outputPrefix), ".traveling_wave.newton.txt"

  ! Load the initial condition.
  call loadTravelingWave(region,trim(restartFilename))
  startTimestep = region%timestep
  totalSubsteps = 0

  if (region%simulationFlags%enableBodyForce .or. region%simulationFlags%checkConservation) then
    region%oneOverVolume = computeRegionIntegral(region)
    region%oneOverVolume = 1.0_wp / region%oneOverVolume

    region%momentumLossPerVolume = 0.0_wp
  end if

  do i = 1, size(region%states) !... update state
     call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
          region%solverOptions)
  end do

  ! Check if physical quantities are within allowed limits.
  if (region%simulationFlags%enableSolutionLimits) then
    if ( checkSolutionLimits(region, FORWARD, this%outputPrefix) ) then
      write(message,'(A)') 'Newton iteration reached an infeasible solution!'
      call gracefulExit(region%comm, message)
    end if
  end if

  call computeTravelingWaveResidual(region,this%residual)
  residualNorm = sqrt(regionInnerProduct(this%residual,this%residual,region))

  do timestep = startTimestep + 1, startTimestep + this%nTimesteps

    region%timestep = timestep

    call saveRegionVector(this%base,region,QOI_FORWARD_STATE)
    error = residualNorm
    this%x%params = 0.0_wp
    do i = 1, size(region%states)
      this%x%states(i)%conservedVariables = 0.0_wp
    end do

    ! Ns = 1
    ! do i = 1, size(region%states)
    !   Ns = Ns + size(region%states(i)%conservedVariables,1)*size(region%states(i)%conservedVariables,2)
    ! end do
    ! allocate(A(Ns,Ns))
    ! offset = 0
    ! do i = 1, size(region%states)
    !   nGridPoints = size(region%states(i)%conservedVariables,1)
    !   nUnknowns = size(region%states(i)%conservedVariables,2)
    !   do j = 1, nGridPoints
    !     do k = 1, nUnknowns
    !       this%x%params = 0.0_wp
    !       do l = 1, size(region%states)
    !         this%x%states(l)%conservedVariables = 0.0_wp
    !       end do
    !       this%x%states(i)%conservedVariables(j,k) = 1.0_wp
    !       index = k + (j-1)*nUnknowns + offset
    !
    !       this%r = computeTravelingWaveJacobian(region,this%x,LINEARIZED)
    !       offset1 = 0
    !       do l = 1, size(region%states)
    !         do m = 1, size(region%states(l)%conservedVariables,1)
    !           A(offset1+1:offset1+nUnknowns,index) = this%r%states(l)%conservedVariables(m,:)
    !           offset1 = offset1 + nUnknowns
    !         end do
    !       end do
    !       A(Ns,index) = this%r%params(1)
    !     end do
    !   end do
    !   offset = offset + nGridPoints*nUnknowns
    ! end do
    ! this%x%params = 1.0_wp
    ! do l = 1, size(region%states)
    !   this%x%states(l)%conservedVariables = 0.0_wp
    ! end do
    ! this%r = computeTravelingWaveJacobian(region,this%x,LINEARIZED)
    ! offset1 = 0
    ! do l = 1, size(region%states)
    !   do m = 1, size(region%states(l)%conservedVariables,1)
    !     A(offset1+1:offset1+nUnknowns,Ns) = this%r%states(l)%conservedVariables(m,:)
    !     offset1 = offset1 + nUnknowns
    !   end do
    ! end do
    ! A(Ns,Ns) = this%r%params(1)
    ! print *, Ns
    ! print *, size(A)
    ! open(unit=getFreeUnit(fileUnit),file='Jacobian.bin',status='replace',form='unformatted',access='stream')
    ! do i = 1,Ns
    !   ! print *, i,'-th column'
    !   ! print *, A(:,i)
    !   write(fileUnit) A(:,i)
    ! end do
    ! close(fileUnit)
    ! SAFE_DEALLOCATE(A)
    ! open(unit=getFreeUnit(fileUnit),file='RightHandSide.bin',status='replace',form='unformatted',access='stream')
    ! do i = 1,size(region%states)
    !   do j = 1, size(region%states(i)%conservedVariables,1)
    !     write(fileUnit) this%residual%states(i)%conservedVariables(j,:)
    !   end do
    ! end do
    ! write(fileUnit) this%residual%params(1)
    ! close(fileUnit)
    !
    ! this%x%params = 0.0_wp
    ! do i = 1, size(region%states)
    !   this%x%states(i)%conservedVariables = 0.0_wp
    ! end do

    do restart = 1, this%maxRestart
      this%r = this%residual * (-1.0_wp)                                  &
                         - computeTravelingWaveJacobian(region,this%x,LINEARIZED)
      rNorm = sqrt(regionInnerProduct(this%r,this%r,region))

      this%H = 0.0_wp
      this%sn = 0.0_wp
      this%cs = 0.0_wp
      this%beta = 0.0_wp
      this%Q(1) = this%r * (1.0_wp/rNorm)
      this%beta(1) = rNorm
      do substep = 1, this%maxSubsteps
        this%H(1:substep+1,substep) = arnoldi(this,region,substep)

        call givensRotation(this%H(1:substep+1,substep),this%cs(1:substep),this%sn(1:substep),substep)

        ! update the residual vector
        this%beta(substep+1) = - this%sn(substep) * this%beta(substep)
        this%beta(substep) = this%cs(substep) * this%beta(substep)

        error = abs(this%beta(substep+1))
        numSteps = substep

        call this%showProgress(region, GMRES, totalSubsteps + substep, (/ error /), &
                 outputFilename, totalSubsteps + substep - 1 > this%reportInterval)

        if (error<this%linminTol*residualNorm) then
          exit
        end if

      end do

      totalSubsteps = totalSubsteps + numSteps

      allocate(y(numSteps))
      y = backwardSubstitution(this%H(1:numSteps,1:numSteps),this%beta(1:numSteps))
      do i = 1, numSteps
        ! .. newton equation solve for (-residual)
        this%x = this%x + this%Q(i) * y(i)
      end do
      SAFE_DEALLOCATE(y)

      if (error<this%linminTol*residualNorm) then
        exit
      end if
    end do

    if (error>=this%linminTol*residualNorm) then
      write(message,'(A)') 'GMRES iteration failed.'
      call gracefulExit(region%comm, message)
    end if

    call loadRegionVector(region,this%base + this%x,QOI_FORWARD_STATE)

    do i = 1, size(region%states) !... update state
       call region%states(i)%update(region%grids(i), region%simulationFlags,                   &
            region%solverOptions)
    end do

    ! Check if physical quantities are within allowed limits.
    if (region%simulationFlags%enableSolutionLimits) then
      if ( checkSolutionLimits(region, FORWARD, this%outputPrefix) ) then
        write(message,'(A)') 'Newton iteration reached an infeasible solution!'
        call gracefulExit(region%comm, message)
      end if
    end if

    call computeTravelingWaveResidual(region,this%residual)
    residualNorm1 = residualNorm
    residualNorm = sqrt(regionInnerProduct(this%residual,this%residual,region))

    call this%showProgress(region, NEWTON, timestep, (/ error /),                              &
                           filename, timestep - startTimestep - 1 > this%reportInterval)

    if (residualNorm > residualNorm1) then
      write(message,'(2(A,ES13.6))') 'Newton iteration failed! Previous: ', residualNorm1,     &
                                     ', Current: ', residualNorm
      call gracefulExit(region%comm, message)
    end if

  end do

  call endTiming("runGMRES")

end subroutine runGMRES

subroutine runBICGSTAB(this, region, restartFilename)

  ! <<< External modules >>>
  use MPI
  use iso_fortran_env, only : output_unit

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Optimizer_mod, only : t_Optimizer
  use RegionVector_mod

  ! <<< Enumerations >>>
  ! use State_enum, only : QOI_FORWARD_STATE, QOI_ADJOINT_STATE, QOI_RIGHT_HAND_SIDE
  use Region_enum, only : FORWARD, LINEARIZED

  ! <<< Internal modules >>>
  ! use OptimizerImpl, only : arnoldi, givensRotation, backwardSubstitution
  use RegionVectorImpl
  use TravelingWaveImpl
  use SolverImpl, only : checkSolutionLimits
  use RegionImpl, only : computeRegionIntegral
  use MPITimingsHelper, only : startTiming, endTiming
  use ErrorHandler, only : writeAndFlush, gracefulExit
  ! use InputHelper, only : getOption, getRequiredOption, getFreeUnit

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this
  class(t_Region) :: region
  character(len = *), intent(in) :: restartFilename

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, fileUnit, ostat, procRank, ierror,                              &
             startTimestep, timestep
  SCALAR_TYPE :: a, w, rr, rr0, rr1, b, residualNorm
  character(len = STRING_LENGTH) :: filename, outputFilename, message

  call startTiming("runBICGSTAB")

  write(outputFilename, '(2A)') trim(this%outputPrefix), ".traveling_wave.bicgstab.txt"

  ! Load the initial condition.
  call loadTravelingWave(region,trim(restartFilename))
  startTimestep = region%timestep

  if (region%simulationFlags%enableBodyForce .or. region%simulationFlags%checkConservation) then
    region%oneOverVolume = computeRegionIntegral(region)
    region%oneOverVolume = 1.0_wp / region%oneOverVolume

    region%momentumLossPerVolume = 0.0_wp
  end if

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
  ! residualNorm = sqrt(regionInnerProduct(this%residual,this%residual,region))
  ! error = 1.0_wp
  this%x%params = 0.0_wp
  do i = 1, size(region%states)
    this%x%states(i)%conservedVariables = 0.0_wp
  end do

  this%r = this%residual * (-1.0_wp)                                            &
                        - computeTravelingWaveJacobian(region,this%x,LINEARIZED)
  ! this%rcg = computeTravelingWaveJacobian(region,this%r,LINEARIZED)
  this%rcg = this%r
  this%conjGrad = this%r

  residualNorm = regionInnerProduct(this%r,this%r,region)
  rr0 = regionInnerProduct(this%r,this%rcg,region)
  if (rr0<this%cgTol*residualNorm) then
    write(message,'(A)') 'Initial vectors are chosen to be trivial.'
    call gracefulExit(region%comm, message)
  end if

  do timestep = startTimestep, startTimestep + this%nTimesteps
    region%timestep = timestep

    this%Ap = computeTravelingWaveJacobian(region,this%conjGrad,LINEARIZED)
    a = rr0 / regionInnerProduct(this%Ap,this%rcg,region)

    this%s = this%r - this%Ap * a
    this%As = computeTravelingWaveJacobian(region,this%s,LINEARIZED)
    w = regionInnerProduct(this%As,this%s,region)                              &
        / regionInnerProduct(this%As,this%As,region)

    this%x = this%x + this%conjGrad * a + this%s * w
    this%r = this%s - this%As * w

    rr = regionInnerProduct(this%r,this%r,region)
    write(message,'(I8.8,1X,A,1X,'// SCALAR_FORMAT //')') timestep, ' error: ', rr
    call writeAndFlush(region%comm, output_unit, message)

    if (rr<this%cgTol*residualNorm) then
      write(message,'(A)') 'BICGSTAB procedure finished.'
      call writeAndFlush(region%comm, output_unit, message)
      exit
    end if

    rr1 = regionInnerProduct(this%r,this%rcg,region)
    b = a / w * rr1 / rr0
    this%conjGrad = this%r + (this%conjGrad - this%Ap * w) * b

  end do

  call endTiming("runBICGSTAB")

end subroutine runBICGSTAB
