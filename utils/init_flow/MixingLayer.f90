#include "config.h"

module MixingLayer_mod

  implicit none
  public

  ! <<< Global parameters >>>
  integer :: nx, ny, nz
  integer, dimension(3) :: iStart, iEnd, interiorStart, interiorEnd,                         &
       targetStart, targetEnd, controlStart, controlEnd
  real(SCALAR_KIND), dimension(3) :: minLengthInterior, minLengthOuter,                      &
       maxLengthInterior, maxLengthOuter

end module MixingLayer_mod

subroutine mixingLayerGrid(region)

  ! <<< Private members >>>
  use MixingLayer_mod

  ! <<< External modules >>>
  use MPI

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  ! <<< Derived types >>>
  use Region_mod, only : t_Region

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, procRank, numProcs, ierror
  real(wp) :: gridSpacing(3), b, c, sigma, minMeshsize, maxMeshsize, y1, y2
  real(wp), allocatable, dimension(:) :: s, g
  logical :: stretchX, stretchY
  character(len = STRING_LENGTH) :: key, message

  assert(size(region%grids) == 1)

  nDimensions = region%grids(1)%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  ! Get process information.
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  ! Read in the grid size.
  minLengthInterior = 0.0_wp; minLengthOuter = 0.0_wp
  maxLengthInterior = 0.0_wp; maxLengthOuter = 0.0_wp
  do i = 1, nDimensions
     write(key, '(A,I3.3,A,I1.1,A)') "grid", 1, "/dir", i, "/min_length_outer"
     call getRequiredOption(trim(key), minLengthOuter(i))
     write(key, '(A,I3.3,A,I1.1,A)') "grid", 1, "/dir", i, "/min_length_interior"
     minLengthInterior(i) = getOption(trim(key), minLengthOuter(i))
     if (minLengthInterior(i) < minLengthOuter(i)) then
        write(message, '(A,F6.2)') trim(key) // " must be >=", minLengthOuter(i)
        call gracefulExit(region%comm, message)
     end if

     write(key, '(A,I3.3,A,I1.1,A)') "grid", 1, "/dir", i, "/max_length_outer"
     call getRequiredOption(trim(key), maxLengthOuter(i))
     write(key, '(A,I3.3,A,I1.1,A)') "grid", 1, "/dir", i, "/max_length_interior"
     maxLengthInterior(i) = getOption(trim(key), maxLengthOuter(i))
     if (maxLengthInterior(i) > maxLengthOuter(i)) then
        write(message, '(A,F6.2)') trim(key) // " must be <=", maxLengthOuter(i)
        call gracefulExit(region%comm, message)
     end if
  end do

  ! Store some grid information.
  nx = region%grids(1)%globalSize(1)
  ny = region%grids(1)%globalSize(2)
  nz = region%grids(1)%globalSize(3)
  gridSpacing = 0.0_wp
  do i = 1, nDimensions
     gridSpacing(i) = (maxLengthOuter(i) - minLengthOuter(i)) /                              &
          real(region%grids(1)%globalSize(i), wp)
  end do

  ! Get the grid partition.
  iStart = 1; iEnd = 1
  do i = 1, region%grids(1)%nDimensions
     iStart(i) = region%grids(1)%offset(i) + 1
     iEnd(i)   = region%grids(1)%offset(i) + region%grids(1)%localSize(i)
  end do

  ! Should we stretch the mesh?
  stretchX = getOption('stretch_x',.false.)
  stretchY = getOption('stretch_y',.false.)

  ! Generate the grid.
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Create X.
           if (.not. stretchX) region%grids(1)%coordinates(grid_index(i,j,k), 1) =           &
                (maxLengthOuter(1) - minLengthOuter(1)) * real(i-1, wp) / real(nx-1, wp) +   &
                minLengthOuter(1)

           ! Create Y.
           if (.not. stretchY) region%grids(1)%coordinates(grid_index(i,j,k), 2) =           &
                (maxLengthOuter(2) - minLengthOuter(2)) * real(j-1,wp) / real(ny-1, wp) + &
                minLengthOuter(2)

           ! Create Z.
           if (nz > 1) region%grids(1)%coordinates(grid_index(i,j,k), 3) =                   &
                (maxLengthOuter(3) - minLengthOuter(3) - gridSpacing(3)) * real(k-1, wp) /   &
                real(nz-1, wp) + minLengthOuter(3)
        end do
     end do
  end do

  ! Stretch the grid.
  if (stretchX) then
     ! Parameters
     sigma=0.2_wp
     b=20.0_wp
     c=0.9_wp

     ! Create uniform spacing.
     allocate(s(nx))
     do i = 1, nx
        s(i) = real(i - 1, wp) / (nx - 1)
     end do

     ! Compute mapping g(s).
     allocate(g(nx))
     call mapping_function(s, b, c, sigma, g)

     ! Find min/max spacing.
     minMeshsize =  huge(1.0_wp)
     maxMeshsize = -huge(1.0_wp)

     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)
              ! Create x.
              region%grids(1)%coordinates(grid_index(i,j,k), 1) = minLengthOuter(1) +        &
                   0.5_wp * (maxLengthOuter(1) - minLengthOuter(1)) * (1.0_wp + g(i))

              ! Find min/max spacing.
              if (i > 1) then
                 y1 = region%grids(1)%coordinates(grid_index(i,j,k), 1)
                 y2 = region%grids(1)%coordinates(grid_index(i-1,j,k), 1)
                 minMeshsize=min(minMeshsize, abs(y2 - y1))
                 maxMeshsize=max(maxMeshsize, abs(y2 - y1))
              end if
           end do
        end do
     end do
     call MPI_Allreduce(MPI_IN_PLACE, maxMeshsize, 1, REAL_TYPE_MPI, MPI_MAX,                &
          MPI_COMM_WORLD, ierror)
     call MPI_Allreduce(MPI_IN_PLACE, minMeshsize, 1, REAL_TYPE_MPI, MPI_MIN,                &
          MPI_COMM_WORLD, ierror)
     if (procRank == 0) then
        print *
        print *, 'min/max x-spacing:',minMeshsize, maxMeshsize
        print *
     end if

     deallocate(s, g)
  end if

  if (stretchY) then
     ! Grid stretching parameters.
     sigma = 0.2_wp
     b = 20.0_wp
     c = 0.62_wp

     ! Create uniform spacing.
     allocate(s(ny))
     do j = 1, ny
        s(j) = real(j - 1, wp) / (ny - 1)
     end do

     ! Compute mapping g(s).
     allocate(g(ny))
     call mapping_function(s, b, c, sigma, g)

     ! Find min/max spacing.
     minMeshsize =  huge(1.0_wp)
     maxMeshsize = -huge(1.0_wp)

     do k = iStart(3), iEnd(3)
        do j = iStart(2), iEnd(2)
           do i = iStart(1), iEnd(1)
              ! Create y.
              region%grids(1)%coordinates(grid_index(i,j,k), 2) = minLengthOuter(2) +        &
                   0.5_wp * (maxLengthOuter(2) - minLengthOuter(2)) * (1.0_wp + g(j))

              ! Find min/max spacing.
              if (j > 1) then
                 y1 = region%grids(1)%coordinates(grid_index(i,j,k), 2)
                 y2 = region%grids(1)%coordinates(grid_index(i,j-1,k), 2)
                 minMeshsize=min(minMeshsize, abs(y2 - y1))
                 maxMeshsize=max(maxMeshsize, abs(y2 - y1))
              end if
           end do
        end do
     end do
     call MPI_Allreduce(MPI_IN_PLACE, maxMeshsize, 1, REAL_TYPE_MPI, MPI_MAX,                &
          MPI_COMM_WORLD, ierror)
     call MPI_Allreduce(MPI_IN_PLACE, minMeshsize, 1, REAL_TYPE_MPI, MPI_MIN,                &
          MPI_COMM_WORLD, ierror)
     if (procRank == 0) then
        print *
        print *, 'min/max y-spacing:',minMeshsize, maxMeshsize
        print *
     end if

     deallocate(s, g)
  end if

  ! Find extents of interior region (used to define sponge extents).
  interiorStart = 1; interiorEnd = -1
  j = iStart(2); k = iStart(3)
  interiorStart(1) = nx
  do i = iEnd(1), iStart(1), -1
     if (region%grids(1)%coordinates(grid_index(i,j,k), 1) >= minLengthInterior(1))          &
          interiorStart(1) = i
  end do
  interiorEnd(1) = 1
  do i = iStart(1), iEnd(1)
     if (region%grids(1)%coordinates(grid_index(i,j,k), 1) <= maxLengthInterior(1))          &
          interiorEnd(1) = i
  end do
  i = iStart(1)
  interiorStart(2) = ny
  do j = iEnd(2), iStart(2), -1
     if (region%grids(1)%coordinates(grid_index(i,j,k), 2) >= minLengthInterior(2))          &
          interiorStart(2) = j
  end do
  interiorEnd(2) = 1
  do j = iStart(2), iEnd(2)
     if (region%grids(1)%coordinates(grid_index(i,j,k), 2) <= maxLengthInterior(2))          &
          interiorEnd(2) = j
  end do
  call MPI_Allreduce(MPI_IN_PLACE, interiorStart(1), 1, MPI_INTEGER, MPI_MIN,                &
       MPI_COMM_WORLD, ierror)
  call MPI_Allreduce(MPI_IN_PLACE, interiorEnd(1), 1, MPI_INTEGER, MPI_MAX,                  &
       MPI_COMM_WORLD, ierror)
  call MPI_Allreduce(MPI_IN_PLACE, interiorStart(2), 1, MPI_INTEGER, MPI_MIN,                &
       MPI_COMM_WORLD, ierror)
  call MPI_Allreduce(MPI_IN_PLACE, interiorEnd(2), 1, MPI_INTEGER, MPI_MAX,                  &
       MPI_COMM_WORLD, ierror)

  if (procRank == 0) then
     write (*,'(A)') ''
     write (*,'(A, i0.0, A, i0.0, A, i0.0, A, i0.0, A)') 'Interior domain extents: [',       &
          interiorStart(1), ', ', interiorEnd(1), '] x [',                                   &
          interiorStart(2), ', ', interiorEnd(2), ']'
  end if

contains

  function grid_index(i, j, k) result(l)
    implicit none

    integer, intent(in) :: i, j, k
    integer :: l
    
    l = i - region%grids(1)%offset(1) + region%grids(1)%localSize(1) *                       &
         (j - 1 - region%grids(1)%offset(2) + region%grids(1)%localSize(2) *                 &
         (k - 1 - region%grids(1)%offset(3)))
    
    return
  end function grid_index

  Subroutine mapping_function(s, b, c, sigma, g)

    use MathHelper, only : pi

    implicit none

    integer, parameter :: wp = SCALAR_KIND
    real(wp), intent(in) :: s(:), b, c, sigma
    real(wp), intent(out) :: g(size(s))

    g = ((s - 0.5_wp) * (1.0_wp + 2.0_wp * b) - b * sigma *                                  &
         (exp(- ((s - 0.5_wp + c) / sigma) ** 2) / sqrt(pi) +                                &
         ((s - 0.5_wp + c) / sigma) * erf((s - 0.5_wp + c) / sigma) -                        &
         exp(- ((s - 0.5_wp - c) / sigma) ** 2) / sqrt(pi) -                                 &
         ((s - 0.5_wp - c) / sigma) * erf((s - 0.5_wp - c) / sigma))) /                      &
         (0.5_wp + b - b * sigma * (exp(- ((0.5_wp + c) / sigma) ** 2) /                     &
         sqrt(pi) + ((0.5_wp + c) / sigma) * erf((0.5_wp + c) / sigma) -                     &
         exp(- ((0.5_wp - c) / sigma) ** 2) / sqrt(pi) - ((0.5_wp - c) / sigma) *            &
         erf((0.5_wp - c) / sigma)))

    return
  end subroutine mapping_function

end subroutine mixingLayerGrid

subroutine mixingLayerQ(region)

  ! <<< Private members >>>
  use MixingLayer_mod

  ! <<< External modules >>>
  use MPI

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit
  use ThermoChemistry, only : getMolecularWeight

  ! <<< Enumerations >>>
  use SolverOptions_enum

  ! <<< Derived types >>>
  use Region_mod, only : t_Region

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, k, nDimensions, nSpecies, H2, O2, N2, procRank, numProcs, ierror
  real(wp) :: ratioOfSpecificHeats, x, y, density, velocity, temperature, pressure,         &
       fuel, oxidizer, inert, Z, Z0, YF0, Yo0, upperVelocity, lowerVelocity, growthRate,    &
       outflowVelocity
  real(wp), dimension(:), allocatable :: Wi
  character(len = STRING_LENGTH) :: message

  ! Get process information.
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  assert(size(region%grids) == 1)

  nDimensions = region%grids(1)%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nSpecies = region%solverOptions%nSpecies
  assert(nSpecies >= 0)
  if (nspecies.gt.2) then
     write(message, '(A)') "WARNING, max of 2 species (for now)"
     call gracefulExit(region%comm, message)
  end if

  ! Get species indices.
  if (allocated(region%solverOptions%speciesName)) then
     do k = 1, nSpecies + 1
        select case (trim(region%solverOptions%speciesName(k)))
        case ('H2', 'HYDROGEN')
           H2 = k
        case ('O2', 'OXYGEN')
           O2 = k
        case ('N2', 'NITROGEN')
           N2 = k
        case default
           write(message, '(3A)') "Unknown species: ",                                       &
                trim(region%solverOptions%speciesName(k)), "!"
           call gracefulExit(region%comm, message)
        end select
     end do
  else
     H2 = 1
     O2 = 2
     N2 = nSpecies + 1
  end if

  ! Get molecular weights.
  if (region%solverOptions%equationOfState == IDEAL_GAS_MIXTURE) then
     allocate(Wi(nSpecies+1))
     Wi = region%solverOptions%molecularWeightInverse
  end if

  ! Get the ratio of specific heats.
  ratioOfSpecificHeats = region%solverOptions%ratioOfSpecificHeats

  ! Mixing layer velocities.
  upperVelocity = getOption("upper_velocity", 0.0_wp)
  lowerVelocity = getOption("lower_velocity", 0.0_wp)
  growthRate = getOption("growth_rate", 0.1_wp)
  outflowVelocity = getOption("outflow_velocity", 0.0_wp)

  ! Mixture properties.
  Z0 = getOption("initial_mixture_fraction", 1.0_wp)
  call getRequiredOption("initial_fuel_mass_fraction", Yf0)
  call getRequiredOption("initial_oxidizer_mass_fraction", Yo0)

  do i = 1, region%grids(1)%nGridPoints

     ! Local coordinates.
     x = region%grids(1)%coordinates(i,1)
     y = region%grids(1)%coordinates(i,2)

     ! Velocity.
     velocity = lowerVelocity +                                                              &
          0.5_wp * (upperVelocity - lowerVelocity) * (1.0_wp + tanh(2.0_wp * y /             &
          (1.0_wp + growthRate*max(0.0_wp, x - minLengthInterior(1)))))

     ! Accelerate the outflow.
     velocity = velocity + (1.0_wp - lowerVelocity) * outflowVelocity *                      &
          max( 0.0_wp, (x - minLengthInterior(1)) /                                          &
          (maxLengthOuter(1) - maxLengthInterior(1)) )

     ! Pressure
     pressure = 1.0_wp / ratioOfSpecificHeats

     ! Temperature.
     temperature =  1.0_wp / (ratioOfSpecificHeats - 1.0_wp)

     ! Mixture fraction.
     Z = 0.5_wp * Z0 * ( 1.0_wp - erf(y / (1.0_wp +                                          &
          growthRate * max(0.0_wp, min(x, maxLengthInterior(1)) - minLengthInterior(1)))) )

     ! Mass fractions.
     fuel = YF0 * Z
     oxidizer = YO0 * (1.0_wp - Z)

     ! Species decay at the outflow.
     fuel = fuel - fuel * max( 0.0_wp, (x - maxLengthInterior(1)) /                          &
          (maxLengthOuter(1) - maxLengthInterior(1)) )
     oxidizer = oxidizer - oxidizer * max( 0.0_wp, (x - maxLengthInterior(1)) /              &
          (maxLengthOuter(1) - maxLengthInterior(1)) )

     ! Correct species mass fractions.
     inert = 1.0_wp - fuel - oxidizer
     if (inert < 0.0_wp) oxidizer = oxidizer + inert

     ! Get density from the equation of state
     select case (region%solverOptions%equationOfState)
     case(IDEAL_GAS)
        density = ratioOfSpecificHeats * pressure /                                          &
             (temperature * (ratioOfSpecificHeats - 1.0_wp))
     case (IDEAL_GAS_MIXTURE)
        density = ratioOfSpecificHeats * pressure /                                          &
             ( temperature * (ratioOfSpecificHeats - 1.0_wp) *                               &
             (fuel * (Wi(H2) - Wi(N2)) + oxidizer * (Wi(O2) - Wi(N2)) + Wi(N2)) )
     end select

     ! State variables.
     region%states(1)%conservedVariables(i,1) = density
     region%states(1)%conservedVariables(i,2) = region%states(1)%conservedVariables(i,1) *   &
          velocity
     region%states(1)%conservedVariables(i,3:nDimensions+1) = 0.0_wp
     region%states(1)%conservedVariables(i,nDimensions+2) = pressure /                       &
          (ratioOfSpecificHeats - 1.0_wp) +                                                  &
          0.5_wp * region%states(1)%conservedVariables(i,1) * velocity ** 2
     if (nSpecies.gt.0) region%states(1)%conservedVariables(i,nDimensions+2+H2) = fuel *     &
          region%states(1)%conservedVariables(i,1)
     if (nSpecies.gt.1) region%states(1)%conservedVariables(i,nDimensions+2+O2) = oxidizer * &
          region%states(1)%conservedVariables(i,1)

  end do

  ! Cleanup.
  SAFE_DEALLOCATE(Wi)

end subroutine mixingLayerQ

subroutine mixingLayerTargetMollifier(region)

  ! <<< Private members >>>
  use MixingLayer_mod

  ! <<< External modules >>>
  use MPI

  ! <<< Internal modules >>>
  use InputHelper, only : getOption
  use ErrorHandler, only : gracefulExit

  ! <<< Derived types >>>
  use Region_mod, only : t_Region

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, procRank, ierror
  real(wp) :: x
  real(wp), parameter :: s = 10.0_wp, r = 0.2_wp, eps = 1.0E-6_wp

  ! Initialize mollifier extents.
  targetStart = 1; targetEnd = -1

  ! Return for now.
  return

  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

end subroutine mixingLayerTargetMollifier

subroutine mixingLayerControlMollifier(region)

  ! <<< Private members >>>
  use MixingLayer_mod

  ! <<< External modules >>>
  use MPI

  ! <<< Internal modules >>>
  use InputHelper, only : getOption
  use ErrorHandler, only : gracefulExit

  ! <<< Derived types >>>
  use Region_mod, only : t_Region

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, procRank, ierror
  real(wp) :: x
  real(wp), parameter :: s = 10.0_wp, r = 0.2_wp, eps = 1.0E-6_wp

  ! Initialize mollifier extents.
  controlStart = 1; controlEnd = -1

  ! Return for now.
  return

  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

end subroutine mixingLayerControlMollifier

subroutine mixingLayerBC

  ! <<< Private members >>>
  use MixingLayer_mod, only : interiorStart, interiorEnd, targetStart, targetEnd,            &
       controlStart, controlEnd

  ! <<< External modules >>>
  use MPI

  ! <<< Internal modules >>>
  use InputHelper, only : getFreeUnit

  implicit none

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, bc, nbc, iunit, procRank, ierror
  integer, allocatable, dimension(:) :: normDir, imin, imax, jmin, jmax, kmin, kmax
  character(len = 22), allocatable, dimension(:) :: name, type

  ! Only root process writes boundary conditions
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
  if (procRank /= 0) return

  ! Initialize a large number of boundary conditions.
  nbc = 99

  ! Allocate boundary conditions.
  allocate(name(nbc), type(nbc), normDir(nbc),                                               &
       imin(nbc), imax(nbc), jmin(nbc), jmax(nbc), kmin(nbc), kmax(nbc))

  ! Set the boundary conditions.

  ! Inflow BC
  bc = 1
  name   (bc) = 'inflow'
  type   (bc) = 'SAT_FAR_FIELD'
  normDir(bc) =  1
  imin   (bc) =  1
  imax   (bc) =  1
  jmin   (bc) =  1
  jmax   (bc) = -1
  kmin   (bc) =  1
  kmax   (bc) = -1

  ! Inflow sponge
  bc = bc+1
  name   (bc) = 'inflowSponge'
  type   (bc) = 'SPONGE'
  normDir(bc) =  1
  imin   (bc) =  1
  imax   (bc) =  interiorStart(1)
  jmin   (bc) =  1
  jmax   (bc) = -1
  kmin   (bc) =  1
  kmax   (bc) = -1

  ! Outflow BC
  bc = bc+1
  name   (bc) = 'outflow'
  type   (bc) = 'SAT_FAR_FIELD'
  normDir(bc) = -1
  imin   (bc) = -1
  imax   (bc) = -1
  jmin   (bc) =  1
  jmax   (bc) = -1
  kmin   (bc) =  1
  kmax   (bc) = -1

  ! Outflow sponge
  bc = bc+1
  name   (bc) = 'outflowSponge'
  type   (bc) = 'SPONGE'
  normDir(bc) = -1
  imin   (bc) =  interiorEnd(1)
  imax   (bc) = -1
  jmin   (bc) =  1
  jmax   (bc) = -1
  kmin   (bc) =  1
  kmax   (bc) = -1

  ! Bottom BC
  bc = bc+1
  name   (bc) = 'bottom'
  type   (bc) = 'SAT_FAR_FIELD'
  normDir(bc) =  2
  imin   (bc) =  1
  imax   (bc) = -1
  jmin   (bc) =  1
  jmax   (bc) =  1
  kmin   (bc) =  1
  kmax   (bc) = -1

  ! Bottom sponge
  bc = bc+1
  name   (bc) = 'bottomSponge'
  type   (bc) = 'SPONGE'
  normDir(bc) =  2
  imin   (bc) =  1
  imax   (bc) = -1
  jmin   (bc) =  1
  jmax   (bc) =  interiorStart(2)
  kmin   (bc) =  1
  kmax   (bc) = -1

  ! Top BC
  bc = bc+1
  name   (bc) = 'top'
  type   (bc) = 'SAT_FAR_FIELD'
  normDir(bc) = -2
  imin   (bc) =  1
  imax   (bc) = -1
  jmin   (bc) = -1
  jmax   (bc) = -1
  kmin   (bc) =  1
  kmax   (bc) = -1

  ! Top sponge
  bc = bc + 1
  name   (bc) = 'topSponge'
  type   (bc) = 'SPONGE'
  normDir(bc) = -2
  imin   (bc) =  1
  imax   (bc) = -1
  jmin   (bc) =  interiorEnd(2)
  jmax   (bc) = -1
  kmin   (bc) =  1
  kmax   (bc) = -1

  ! Excitation support
  bc = bc + 1
  name   (bc) = 'excitationSupport'
  type   (bc) = 'SOLENOIDAL_EXCITATION'
  normDir(bc) =  0
  imin   (bc) =  1
  imax   (bc) =  interiorStart(1)
  jmin   (bc) =  1
  jmax   (bc) = -1
  kmin   (bc) =  1
  kmax   (bc) = -1

  ! Target region
  bc = bc+1
  name   (bc) = 'targetRegion'
  type   (bc) = 'COST_TARGET'
  normDir(bc) =  0
  imin   (bc) =  targetStart(1)
  imax   (bc) =  targetEnd(1)
  jmin   (bc) =  targetStart(2)
  jmax   (bc) =  targetEnd(2)
  kmin   (bc) =  targetStart(3)
  kmax   (bc) =  targetEnd(3)

  ! Control region
  bc = bc+1
  name   (bc) = 'controlRegion'
  type   (bc) = 'ACTUATOR'
  normDir(bc) =  0
  imin   (bc) =  controlStart(1)
  imax   (bc) =  controlEnd(1)
  jmin   (bc) =  controlStart(2)
  jmax   (bc) =  controlEnd(2)
  kmin   (bc) =  controlStart(3)
  kmax   (bc) =  controlEnd(3)

  ! Overwrite number of boundary conditions.
  nbc = bc

  ! Open the file.
  open(unit = getFreeUnit(iunit), file = "bc.dat")

  write (*,'(A)') ''
  write (*,'(A)') 'Writing boundary conditions'

  ! Write the header.
  write(iunit,'(1a87)') "# Name                 Type                  Grid normDir iMin iMax jMin jMax kMin kMax"
  write(iunit,'(1a87)') "# ==================== ===================== ==== ======= ==== ==== ==== ==== ==== ===="

  ! Write the boundary conditions.
  do i = 1, nbc
     write(iunit,'(a2,2a21,8I5)') '  ', adjustl(name(i)), adjustl(type(i)),                  &
          1, normDir(i), imin(i), imax(i), jmin(i), jmax(i), kmin(i), kmax(i)
  end do

  ! Close the file.
  close(iunit)

  ! Clean up.
  deallocate(name, type, normDir, imin, imax, jmin, jmax, kmin, kmax)

end subroutine mixingLayerBC
