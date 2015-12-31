#include "config.h"

module Premixed_mod

  implicit none
  public

  ! <<< Global parameters >>>
  integer :: nx, ny, nz
  integer, dimension(3) :: iStart, iEnd, targetStart, targetEnd, controlStart, controlEnd
  real(SCALAR_KIND), dimension(3) :: minLengthInterior, minLengthOuter,                      &
       maxLengthInterior, maxLengthOuter

end module Premixed_mod

subroutine premixedGrid(region)

  ! <<< Private members >>>
  use Premixed_mod

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
  real(wp) :: gridSpacing(3)
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

  ! Generate the grid.
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Create Y.
           if (nx > 1) region%grids(1)%coordinates(grid_index(i,j,k), 1) =                   &
                (maxLengthOuter(1) - minLengthOuter(1) - gridSpacing(1)) * real(i - 1, wp) / &
                real(nx - 1, wp) + minLengthOuter(1)

           ! Create Y.
           if (ny > 1) region%grids(1)%coordinates(grid_index(i,j,k), 2) =                   &
                (maxLengthOuter(2) - minLengthOuter(2) - gridSpacing(2)) * real(j - 1, wp) / &
                real(ny - 1, wp) + minLengthOuter(2)

           ! Create Z.
           if (nz > 1) region%grids(1)%coordinates(grid_index(i,j,k), 3) =                   &
                (maxLengthOuter(3) - minLengthOuter(3) - gridSpacing(3)) * real(k - 1, wp) / &
                real(nz - 1, wp) + minLengthOuter(3)
        end do
     end do
  end do

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

end subroutine premixedGrid

subroutine premixedQ(region)

  ! <<< Private members >>>
  use Premixed_mod

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
  integer :: k, nDimensions, nSpecies, H2, O2, N2, procRank, numProcs, ierror
  real(wp) :: ratioOfSpecificHeats, density, temperature, pressure, fuel, oxidizer,          &
       inert, T0, Z0, YF0, Yo0
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

  ! Mixture properties.
  Z0 = getOption("initial_mixture_fraction", 1.0_wp)
  call getRequiredOption("initial_fuel_mass_fraction", Yf0)
  call getRequiredOption("initial_oxidizer_mass_fraction", Yo0)

  ! Components.
  fuel = YF0 * Z0
  oxidizer = YO0 * (1.0_wp - Z0)

  ! Get the ratio of specific heats.
  ratioOfSpecificHeats = region%solverOptions%ratioOfSpecificHeats

  ! Pressure.
  pressure = 1.0_wp / ratioOfSpecificHeats

  ! Temperature.
  T0 =  1.0_wp / (ratioOfSpecificHeats - 1.0_wp)
  temperature = getOption("initial_temperature", T0)

  ! Get density from the equation of state.
  select case (region%solverOptions%equationOfState)
  case(IDEAL_GAS)
     density = ratioOfSpecificHeats * pressure /                                             &
          (temperature * (ratioOfSpecificHeats - 1.0_wp))
  case (IDEAL_GAS_MIXTURE)
     density = ratioOfSpecificHeats * pressure /                                             &
          ( temperature * (ratioOfSpecificHeats - 1.0_wp) *                                  &
          (fuel * (Wi(H2) - Wi(N2)) + oxidizer * (Wi(O2) - Wi(N2)) + Wi(N2)) )
  end select
  if (procRank == 0) then
     print *
     print *, 'Mixture density = ', density
     print *
  end if

  ! State variables.
  region%states(1)%conservedVariables(:, 1) = density
  region%states(1)%conservedVariables(:, 2:nDimensions+1) = 0.0_wp
  region%states(1)%conservedVariables(:, nDimensions+2) = pressure /                        &
       (ratioOfSpecificHeats - 1.0_wp)
  if (nSpecies > 0) region%states(1)%conservedVariables(:, nDimensions+2+H2) = fuel *       &
       region%states(1)%conservedVariables(:, 1)
  if (nSpecies > 1) region%states(1)%conservedVariables(:, nDimensions+2+O2) = oxidizer *   &
       region%states(1)%conservedVariables(:, 1)

  ! Cleanup.
  SAFE_DEALLOCATE(Wi)

end subroutine premixedQ

subroutine premixedTargetMollifier(region)

  ! <<< Private members >>>
  use Premixed_mod

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

  ! Initialize mollifier extents.
  targetStart = 1; targetEnd = -1

  return

end subroutine premixedTargetMollifier

subroutine premixedControlMollifier(region)

  ! <<< Private members >>>
  use Premixed_mod

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

  ! Initialize mollifier extents.
  controlStart = 1; controlEnd = -1

  return

end subroutine premixedControlMollifier

subroutine premixedBC

  ! <<< Private members >>>
  use Premixed_mod, only : targetStart, targetEnd, controlStart, controlEnd

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

end subroutine premixedBC
