#include "config.h"

module Jet_mod

  implicit none
  public

  ! <<< Global parameters >>>
  integer :: nx, ny, nz
  integer, dimension(3) :: iStart, iEnd, targetStart, targetEnd, controlStart, controlEnd
  real(SCALAR_KIND), dimension(3) :: minLengthInterior, minLengthOuter,                      &
       maxLengthInterior, maxLengthOuter

end module Jet_mod

subroutine jetGrid(region)

  ! <<< Private members >>>
  use Jet_mod

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

end subroutine jetGrid

subroutine jetQ(region)

  ! <<< Private members >>>
  use Jet_mod

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
  character(len = STRING_LENGTH) :: errorMessage
  integer :: i, nDimensions, ierror
  real(wp) :: ratioOfSpecificHeats, machNumber, temperatureRatio,                            &
       axialCoordinateAtNozzleExit, momentumThicknessAtExit, slopeOfMomentumThickness,       &
       momentumThickness, radialCoordinate, normalizedExitVelocity,                          &
       normalizedExitDensity, speedOfSoundAtExit, nozzleLipRadius, potentialCoreLength

  if (region%grids(1)%nDimensions /= 3) then
     write(errorMessage, '(A)')                                                              & 
          "Jet initial condition generator requires a three-dimensional grid!"
  end if

  ratioOfSpecificHeats = region%solverOptions%ratioOfSpecificHeats
  machNumber = getOption("Mach_number", 1.3_wp)
  momentumThicknessAtExit = getOption("nozzle_exit_momentum_thickness", 0.04_wp)
  slopeOfMomentumThickness = getOption("slope_of_momentum_thickness", 0.03_wp)
  axialCoordinateAtNozzleExit = getOption("axial_coordinate_at_nozzle_exit", 0.0_wp)
  temperatureRatio = getOption("temperature_ratio",                                          &
       1.0_wp / (1.0_wp + 0.5_wp * (ratioOfSpecificHeats - 1.0_wp) * machNumber ** 2))
  nozzleLipRadius = getOption("nozzle_lip_radius", 0.5_wp)
  potentialCoreLength = getOption("potential_core_length",                                   &
       4.2_wp + 1.1_wp * machNumber ** 2)

  do i = 1, region%grids(1)%nGridPoints

     momentumThickness = momentumThicknessAtExit + slopeOfMomentumThickness *                &
          max(real(region%grids(1)%coordinates(i,3), wp) - axialCoordinateAtNozzleExit,      &
          0.0_wp)
     radialCoordinate = sqrt(region%grids(1)%coordinates(i,1) ** 2 +                         &
          region%grids(1)%coordinates(i,2) ** 2) / nozzleLipRadius + epsilon(0.0_wp)
     normalizedExitVelocity = 0.5_wp * (1.0_wp + tanh(0.25_wp / momentumThickness *          &
          (1.0_wp / radialCoordinate - radialCoordinate)))
     if (region%grids(1)%coordinates(i,3) > potentialCoreLength)                             &
          normalizedExitVelocity = normalizedExitVelocity * (1.0_wp -                        &
          exp(1.35_wp / (1.0_wp - region%grids(1)%coordinates(i,3) / potentialCoreLength)))
     normalizedExitDensity = 1.0_wp / (0.5_wp * (ratioOfSpecificHeats - 1.0_wp) *            &
          normalizedExitVelocity * (1.0_wp - normalizedExitVelocity) * machNumber ** 2 +     &
          normalizedExitVelocity + (1.0_wp - normalizedExitVelocity) / temperatureRatio)
     speedOfSoundAtExit = sqrt(temperatureRatio)

     region%states(1)%conservedVariables(i,1) = normalizedExitDensity / temperatureRatio
     region%states(1)%conservedVariables(i,2) = 0.0_wp
     region%states(1)%conservedVariables(i,3) = 0.0_wp
     region%states(1)%conservedVariables(i,4) = region%states(1)%conservedVariables(i,1) *   &
          machNumber * speedOfSoundAtExit * normalizedExitVelocity
     region%states(1)%conservedVariables(i,5) =                                              &
          1.0_wp / ratioOfSpecificHeats / (ratioOfSpecificHeats - 1.0_wp) +                  &
          0.5_wp / region%states(1)%conservedVariables(i,1) *                                &
          region%states(1)%conservedVariables(i,4) ** 2
     region%states(1)%conservedVariables(i,5+region%solverOptions%nSpecies) = 0.0_wp

  end do

end subroutine jetQ

subroutine jetTargetMollifier(region)

  ! <<< Private members >>>
  use Jet_mod

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

end subroutine jetTargetMollifier

subroutine jetControlMollifier(region)

  ! <<< Private members >>>
  use Jet_mod

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

end subroutine jetControlMollifier

subroutine jetBC

  ! <<< Private members >>>
  use Jet_mod, only : targetStart, targetEnd, controlStart, controlEnd

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

end subroutine jetBC
