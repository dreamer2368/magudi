#include "config.h"

program jet

  use MPI
  use, intrinsic :: iso_fortran_env

  use Grid_type
  use State_type
  use Region_type

  use Grid_mod, only : setupSpatialDiscretization, updateGrid
  use Region_mod
  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage

  implicit none

  integer :: i, ierror
  character(len = STRING_LENGTH) :: filename
  logical :: success
  type(t_Region) :: region
  integer, allocatable :: globalGridSizes(:,:)

  ! Initialize MPI.
  call MPI_Init(ierror)

  ! Parse options from the input file.
  filename = "jet.inp"
  call parseInputFile(filename)

  ! Verify that the grid file is in valid PLOT3D format and fetch the grid dimensions:
  ! `globalGridSizes(i,j)` is the number of grid points on grid `j` along dimension `i`.
  if (command_argument_count() == 1) then
     call get_command_argument(1, filename)
  else     
     call getRequiredOption("grid_file", filename)
  end if
  call plot3dDetectFormat(MPI_COMM_WORLD, filename,                                          &
       success, globalGridSizes = globalGridSizes)
  if (.not. success) call gracefulExit(MPI_COMM_WORLD, plot3dErrorMessage)

  ! Setup the region and load the grid file.
  call setupRegion(region, MPI_COMM_WORLD, globalGridSizes)
  call loadRegionData(region, QOI_GRID, filename)

  ! Setup spatial discretization.
  do i = 1, size(region%grids)
     call setupSpatialDiscretization(region%grids(i))
  end do
  call MPI_Barrier(region%comm, ierror)

  ! Compute normalized metrics, norm matrix and Jacobian.
  do i = 1, size(region%grids)
     call updateGrid(region%grids(i))
  end do
  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  ! Write out some useful information.
  call reportGridDiagnostics(region)

  ! Generate the initial condition and target state.
  do i = 1, size(region%grids)
     call jetInitialCondition(region%states(i), region%grids(i),                             &
          region%simulationFlags%useTargetState)
  end do

  ! Save initial condition.
  i = len_trim(filename)
  if (filename(i-3:i) == ".xyz") then
     filename = filename(:i-4) // ".ic.q"
  else
     filename = PROJECT_NAME // ".ic.q"
  end if
  call saveRegionData(region, QOI_FORWARD_STATE, filename)

  ! Save target state.
  if (region%simulationFlags%useTargetState) then
     i = len_trim(filename)
     if (filename(i-4:i) == ".ic.q") then
        filename = filename(:i-5) // ".target.q"
     else
        filename = PROJECT_NAME // ".target.q"
     end if
     call saveRegionData(region, QOI_TARGET_STATE, filename)
  end if

  ! Cleanup.
  call cleanupRegion(region)

  ! Finalize MPI.
  call MPI_Finalize(ierror)

contains

  subroutine jetInitialCondition(state, grid, generateTargetState)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Grid_type
    use State_type

    ! <<< Internal modules >>>
    use InputHelper, only : getOption
    use ErrorHandler, only : gracefulExit

    ! <<< Arguments >>>
    type(t_State) :: state
    type(t_Grid) :: grid
    logical, intent(in), optional :: generateTargetState

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    logical :: generateTargetState_
    character(len = STRING_LENGTH) :: errorMessage
    integer :: i, nDimensions, ierror
    real(wp) :: ratioOfSpecificHeats, machNumber, temperatureRatio,                          &
         axialCoordinateAtNozzleExit, momentumThicknessAtExit, slopeOfMomentumThickness,     &
         momentumThickness, radialCoordinate, normalizedExitVelocity,                        &
         normalizedExitDensity, speedOfSoundAtExit, nozzleLipRadius

    generateTargetState_ = .false.
    if (present(generateTargetState)) generateTargetState_ = generateTargetState

    call MPI_Cartdim_get(grid%comm, nDimensions, ierror)
    if (nDimensions /= 3) then
       write(errorMessage, '(A)') &
            "Jet initial condition generator requires a three-dimensional grid!"
    end if

    ratioOfSpecificHeats = getOption("ratio_of_specific_heats", 1.4_wp)
    machNumber = getOption("Mach_number", 1.3_wp)
    momentumThicknessAtExit = getOption("nozzle_exit_momentum_thickness", 0.04_wp)
    slopeOfMomentumThickness = getOption("slope_of_momentum_thickness", 0.03_wp)
    axialCoordinateAtNozzleExit = getOption("axial_coordinate_at_nozzle_exit", 0.0_wp)
    temperatureRatio = getOption("temperature_ratio", 0.7_wp)
    nozzleLipRadius = getOption("nozzle_lip_radius", 0.5_wp)

    do i = 1, grid%nGridPoints

       momentumThickness = momentumThicknessAtExit + slopeOfMomentumThickness *              &
            max(real(grid%coordinates(i,3), wp) - axialCoordinateAtNozzleExit, 0.0_wp)
       radialCoordinate = sqrt(grid%coordinates(i,1) ** 2 + grid%coordinates(i,2) ** 2) /    &
            nozzleLipRadius + epsilon(0.0_wp)
       normalizedExitVelocity = 0.5_wp * (1.0_wp + tanh(0.25_wp / momentumThickness *        &
            (1.0_wp / radialCoordinate - radialCoordinate)))
       normalizedExitDensity = 1.0_wp / (0.5_wp * (ratioOfSpecificHeats - 1.0_wp) *          &
            normalizedExitVelocity * (1.0_wp - normalizedExitVelocity) * machNumber ** 2 +   &
            normalizedExitVelocity + (1.0_wp - normalizedExitVelocity) / temperatureRatio)
       speedOfSoundAtExit = sqrt(temperatureRatio)

       state%conservedVariables(i,1) = normalizedExitDensity / temperatureRatio
       state%conservedVariables(i,2) = state%conservedVariables(i,1) *                       &
            machNumber * normalizedExitVelocity
       state%conservedVariables(i,3) = 0.0_wp
       state%conservedVariables(i,4) = 0.0_wp
       state%conservedVariables(i,5) =                                                       &
            1.0_wp / ratioOfSpecificHeats / (ratioOfSpecificHeats - 1.0_wp) +                &
            0.5_wp / state%conservedVariables(i,1) * state%conservedVariables(i,2) ** 2

       if (generateTargetState_) state%targetState(i,:) = state%conservedVariables(i,:)

    end do
    
  end subroutine jetInitialCondition

end program jet
