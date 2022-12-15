#include "config.h"

program full_rhs_linearized

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver

  use Grid_enum
  use State_enum

  use InputHelper, only : parseInputFile, getFreeUnit, getOption, getRequiredOption
  use InputHelperImpl, only: dict, find
  use ErrorHandler
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  integer :: i, dictIndex, procRank, numProcs, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix, message
  logical :: success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_Solver) :: solver

  interface

     subroutine testLinearizedRelation(solver,region,success,tolerance)
       use Region_mod, only : t_Region
       use Solver_mod, only : t_Solver

       class(t_Solver) :: solver
       class(t_Region) :: region
       logical, intent(out) :: success

       real(SCALAR_KIND), intent(in), optional :: tolerance

     end subroutine testLinearizedRelation

  end interface

  ! Initialize MPI.
  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)

  call initializeErrorHandler()

  if (command_argument_count() .ge. 1) then
     write(message, '(A)') "Usage: magudi [FILE]"
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)') "High-performance Fortran-based adjoint optimization tool."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)')                                                                   &
          "magudi.inp, bc.dat, and grid file are required for this test."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)') "High-performance Fortran-based adjoint optimization tool."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     call cleanupErrorHandler()
     call MPI_Finalize(ierror)
     stop -1
  end if

  call startTiming("total")

  ! Parse options from the input file.
  filename = PROJECT_NAME // ".inp"
  call parseInputFile(filename)

  ! adjoint run options
  call find("enable_adjoint_solver",dictIndex)
  dict(dictIndex)%val = "true"

  outputPrefix = getOption("output_prefix", PROJECT_NAME)

  ! Verify that the grid file is in valid PLOT3D format and fetch the grid dimensions:
  ! `globalGridSizes(i,j)` is the number of grid points on grid `j` along dimension `i`.
  call getRequiredOption("grid_file", filename)
  call plot3dDetectFormat(MPI_COMM_WORLD, filename, success,                                 &
       globalGridSizes = globalGridSizes)
  if (.not. success) call gracefulExit(MPI_COMM_WORLD, plot3dErrorMessage)

  ! Setup the region.
  call region%setup(MPI_COMM_WORLD, globalGridSizes)

  ! Read the grid file.
  call getRequiredOption("grid_file", filename)
  call region%loadData(QOI_GRID, filename)

  ! Update the grids by computing the Jacobian, metrics, and norm.
  do i = 1, size(region%grids)
     call region%grids(i)%update()
  end do
  call MPI_Barrier(region%comm, ierror)

  ! Write out some useful information.
  call region%reportGridDiagnostics()

  ! Save the Jacobian and normalized metrics.
  write(filename, '(2A)') trim(outputPrefix), ".Jacobian.f"
  call region%saveData(QOI_JACOBIAN, filename)
  write(filename, '(2A)') trim(outputPrefix), ".metrics.f"
  call region%saveData(QOI_METRICS, filename)

  ! Initialize the solver.
  call solver%setup(region, outputPrefix = outputPrefix)

  ! Main code logic.
  success = .true.
  call testLinearizedRelation(solver,region,success)

  if ( success ) then
    write(message, '(A)') "Full RHS adjoint test is passed."
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
  else
    write(message, '(A)') "Full RHS adjoint test failed."
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
    stop -1
  end if

  call solver%cleanup()
  call region%cleanup()

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  ! Finalize MPI.
  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

end program full_rhs_linearized


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

subroutine testLinearizedRelation(solver,region,success,tolerance)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver
  use State_mod, only : t_State

  use Region_enum, only : FORWARD, ADJOINT, LINEARIZED
  use State_enum

  ! <<< Internal modules >>>
  use RandomNumber, only : random
  use PLOT3DHelper
  use InputHelper, only : getOption
  use RegionImpl, only : computeRegionIntegral

  ! <<< Arguments >>>
  class(t_Solver) :: solver
  class(t_Region) :: region
  logical, intent(out) :: success
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
  ! type(t_SimulationFlags) :: simulationFlags
  ! type(t_SolverOptions) :: solverOptions
  ! type(t_Grid) :: grid
  type(t_State) :: state0(size(region%grids)), state1(size(region%grids)),            &
                    deltaState(size(region%grids))

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp) :: scalar1, scalar2,&
              stepSizes(32), errorHistory(32), convergenceHistory(31)
  integer :: i, j, k, ierror, procRank, clock
  integer :: nDimensions, stage
  character(len=STRING_LENGTH) :: filename
  real(SCALAR_KIND), allocatable :: F(:,:), deltaPrimitiveVariables(:,:)
  integer, allocatable :: seed(:)

  ! character(len = STRING_LENGTH) :: errorMessage

  call random_seed(size=i)
  allocate(seed(i))
  call system_clock(count=clock)
  seed = clock + 37*(/(j-1,j=1,i)/)
  call random_seed(put=seed)
  deallocate(seed)

  success = .true.

  call MPI_Comm_rank(region%comm, procRank, ierror)

  nDimensions = size(region%globalGridSizes,1)

  ! initialize states
  do i = 1, size(state0)
     call state0(i)%setup(region%grids(i), region%simulationFlags, region%solverOptions)
     call state1(i)%setup(region%grids(i), region%simulationFlags, region%solverOptions)
     call deltaState(i)%setup(region%grids(i), region%simulationFlags, region%solverOptions)
  end do
  call MPI_Barrier(region%comm, ierror)

  ! Randomize conserved variables.
  do k = 1, size(state0)
    do i = 1, region%grids(k)%nGridPoints
       state0(k)%conservedVariables(i,1) = random(0.01_wp, 10.0_wp)
       do j = 1, nDimensions
          state0(k)%conservedVariables(i,j+1) =                                              &
               state0(k)%conservedVariables(i,1) * random(-10.0_wp, 10.0_wp)
       end do
       state0(k)%conservedVariables(i, nDimensions + 2) = state0(k)%conservedVariables(i,1) *    &
            random(0.01_wp, 10.0_wp) / region%solverOptions%ratioOfSpecificHeats +              &
            0.5_wp / state0(k)%conservedVariables(i,1) *                                     &
            sum(state0(k)%conservedVariables(i,2:nDimensions+1) ** 2)
    end do
  end do
  do i = 1, size(state0)
    assert(all(state0(i)%conservedVariables(:,1) > 0.0_wp))
  end do

  ! Compute dependent variables.
  do i = 1, size(state0)
    call state0(i)%update(region%grids(i),region%simulationFlags,region%solverOptions)
  end do

  ! Randomize adjoint variables.
  do i = 1, size(state0)
    allocate(F(region%grids(i)%nGridPoints, region%solverOptions%nUnknowns))
    call random_number(F)
    state0(i)%adjointVariables = F
    SAFE_DEALLOCATE(F)
  end do

  ! Randomize delta conserved variables.
  do k = 1, size(state0)
    allocate(deltaPrimitiveVariables(region%grids(k)%nGridPoints, region%solverOptions%nUnknowns))
    do i = 1, region%grids(k)%nGridPoints
      do j = 1, nDimensions + 2
         deltaPrimitiveVariables(i,j) = random(-1.0_wp, 1.0_wp)
      end do
    end do
    deltaState(k)%conservedVariables(:,1) = deltaPrimitiveVariables(:,1)
    do j = 1, nDimensions
       deltaState(k)%conservedVariables(:,j+1) = state0(k)%conservedVariables(:,j+1) /                     &
            state0(k)%conservedVariables(:,1) * deltaPrimitiveVariables(:,1) +                    &
            state0(k)%conservedVariables(:,1) * deltaPrimitiveVariables(:,j+1)
    end do
    deltaState(k)%conservedVariables(:,nDimensions+2) = state0(k)%conservedVariables(:,nDimensions+2) /    &
         state0(k)%conservedVariables(:,1) * deltaPrimitiveVariables(:,1) +                       &
         sum(state0(k)%conservedVariables(:,2:nDimensions+1) *                                    &
         deltaPrimitiveVariables(:,2:nDimensions+1), dim = 2) +                          &
         state0(k)%conservedVariables(:,1) / region%solverOptions%ratioOfSpecificHeats *                               &
         deltaPrimitiveVariables(:,nDimensions+2)
    SAFE_DEALLOCATE(deltaPrimitiveVariables)
  end do

  ! Set region to state0
  do i = 1, size(state0)
    region%states(i)%conservedVariables = state0(i)%conservedVariables
    call region%states(i)%update(region%grids(i), region%simulationFlags, region%solverOptions)

    region%states(i)%adjointVariables = deltaState(i)%conservedVariables
  end do

  ! Save conserved variable as next target state.
  write(filename,'(2A)') trim(solver%outputPrefix),'.target.q'
  call region%saveData(QOI_FORWARD_STATE, filename)

  stage = random(1,4)
  if (region%simulationFlags%enableBodyForce) then
    region%oneOverVolume = computeRegionIntegral(region)
    call random_number(region%initialXmomentum)
    region%initialXmomentum = region%initialXmomentum * region%oneOverVolume
    if (stage==1) then
      region%momentumLossPerVolume = 0.0_wp
    else
      call random_number(region%momentumLossPerVolume)
    end if
    region%oneOverVolume = 1.0_wp / region%oneOverVolume

    region%adjointMomentumLossPerVolume = 0.0_wp
  end if

  ! Compute baseline rhs
  call region%computeRhs(FORWARD,1,stage)
  do i = 1, size(state0)
    state0(i)%rightHandSide = region%states(i)%rightHandSide
  end do

  ! Compute adjoint rhs for inviscid flux
  call region%computeRhs(LINEARIZED,1,stage)

  ! <R^{\dagger}u, \delta v>
  scalar1 = 0.0_wp
  do i = 1, size(state0)
    scalar1 = scalar1 + region%grids(i)%computeInnerProduct(state0(i)%adjointVariables,      &
                                                            region%states(i)%rightHandSide)
  end do
  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE, scalar1, 1,                          &
       SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(scalar1, 1, SCALAR_TYPE_MPI,                             &
          0, region%grids(i)%comm, ierror)
  end do

  ! <u, \delta R(v)>
  ! Prepare step sizes
  stepSizes(1) = 0.001_wp
  do k = 2, size(stepSizes)
     stepSizes(k) = stepSizes(k-1) * 10.0_wp**(-0.25_wp)
  end do
  errorHistory = 0.0_wp
  do k = 1, size(stepSizes)
    !(1) finite difference on conserved variables
    do i = 1, size(state0)
      region%states(i)%conservedVariables = state0(i)%conservedVariables + stepSizes(k) * deltaState(i)%conservedVariables
      assert(all(region%states(i)%conservedVariables(:,1) > 0.0_wp))

      ! Compute dependent variables.
      call region%states(i)%update(region%grids(i),region%simulationFlags,region%solverOptions)
      assert(all(region%states(i)%specificVolume(:,1) > 0.0_wp))
      assert(all(region%states(i)%temperature(:,1) > 0.0_wp))
    end do

    ! (2)Compute baseline rhs
    call region%computeRhs(FORWARD,1,stage)

    ! (3) <u, \delta R(v)>
    scalar2 = 0.0_wp
    do i = 1, size(state0)
      scalar2 = scalar2 + region%grids(i)%computeInnerProduct(state0(i)%adjointVariables,                             &
                                                    region%states(i)%rightHandSide - state0(i)%rightHandSide)
    end do
    if (region%commGridMasters /= MPI_COMM_NULL)                                               &
         call MPI_Allreduce(MPI_IN_PLACE, scalar2, 1,                          &
         SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

    do i = 1, size(region%grids)
       call MPI_Bcast(scalar2, 1, SCALAR_TYPE_MPI,                             &
            0, region%grids(i)%comm, ierror)
    end do

    errorHistory(k) = abs( (scalar2/stepSizes(k) - scalar1)/scalar1 )
    if (procRank==0)                                          &
      print *, stepSizes(k), scalar1, scalar2/stepSizes(k), errorHistory(k)

    if (k > 1) then
       convergenceHistory(k-1) = log(errorHistory(k) / errorHistory(k-1)) /              &
            log(stepSizes(k) / stepSizes(k-1))
       if (k > 5) then
           if (sum(convergenceHistory(k-3:k-1))/3.0_wp < 0.0_wp) exit
       end if
    end if
  end do

  if (k > 2) then
     call sort(convergenceHistory(:k-2))
     success = success .and. nint(meanTrimmed(convergenceHistory(:k-2))).ge.1
  else
     success = .false.
     if (procRank==0)                                          &
       print *, convergenceHistory
  end if

end subroutine testLinearizedRelation
