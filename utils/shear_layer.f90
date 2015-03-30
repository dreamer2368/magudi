#include "config.h"

program shear_layer

  use MPI
  use, intrinsic :: iso_fortran_env

  use Grid_enum
  use State_type
  use Region_type

  use Region_mod
  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage

  !> Generates the initial condition and target state for a shear layer.

  implicit none

  integer :: i, ierror
  character(len = STRING_LENGTH) :: filename
  logical :: success
  type(t_Region) :: region
  integer, allocatable :: globalGridSizes(:,:)

  ! Initialize MPI.
  call MPI_Init(ierror)

  ! Parse options from the input file.
  filename = "shear_layer.inp"
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
     call region%grids(i)%setupSpatialDiscretization()
  end do
  call MPI_Barrier(region%comm, ierror)

  ! Compute normalized metrics, norm matrix and Jacobian.
  do i = 1, size(region%grids)
     call region%grids(i)%update()
  end do
  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  ! Write out some useful information.
  call reportGridDiagnostics(region)

  ! Generate the initial condition and target state.
  do i = 1, size(region%grids)
     call shearLayerInitialCondition(region%states(i), region%grids(i),                      &
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

  subroutine shearLayerInitialCondition(state, grid, generateTargetState)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
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
    integer :: i, nDimensions, ierror
    real(wp) :: ratioOfSpecificHeats, convectiveVelocity, velocityDifference,                &
         temperatureRatio, slopeOfVorticityThickness,                                        &
         lowerFluidVelocity, upperFluidVelocity,                                             &
         lowerFluidTemperature, upperFluidTemperature,                                       &
         velocity, temperature

    generateTargetState_ = .false.
    if (present(generateTargetState)) generateTargetState_ = generateTargetState

    call MPI_Cartdim_get(grid%comm, nDimensions, ierror)

    ratioOfSpecificHeats = getOption("ratio_of_specific_heats", 1.4_wp)
    convectiveVelocity = getOption("convective_velocity", 0.0_wp)
    velocityDifference = getOption("velocity_difference", 0.0_wp)
    temperatureRatio = getOption("temperature_ratio", 1.0_wp)
    slopeOfVorticityThickness = getOption("slope_of_vorticity_thickness", 0.0_wp)

    lowerFluidVelocity = convectiveVelocity - 0.5_wp * velocityDifference
    upperFluidVelocity = convectiveVelocity + 0.5_wp * velocityDifference

    lowerFluidTemperature = 1.0_wp / (ratioOfSpecificHeats - 1.0_wp)
    upperFluidTemperature = lowerFluidTemperature * temperatureRatio

    do i = 1, grid%nGridPoints

       if (generateTargetState_) then

       velocity = lowerFluidVelocity +                                                       &
            0.5_wp * velocityDifference * (1.0_wp +                                          &
            tanh(2.0_wp * real(grid%coordinates(i,2), wp) / (1.0_wp +                        &
            slopeOfVorticityThickness * max(0.0_wp, real(grid%coordinates(i,1), wp)))))
       temperature = (lowerFluidTemperature * (upperFluidVelocity - velocity) +              &
            upperFluidTemperature * (velocity - lowerFluidVelocity)) /                       &
            (upperFluidVelocity - lowerFluidVelocity) +                                      &
            0.5_wp * (upperFluidVelocity - velocity) * (velocity - lowerFluidVelocity)

          state%targetState(i,1) = 1.0_wp / ((ratioOfSpecificHeats - 1.0_wp) * temperature)
          state%targetState(i,2) = state%targetState(i,1) * velocity
          state%targetState(i,3:nDimensions+1) = 0.0_wp
          state%targetState(i,nDimensions+2) =                                               &
               state%targetState(i,1) * temperature / ratioOfSpecificHeats +                 &
               0.5_wp * state%targetState(i,1) * velocity ** 2

       end if

       velocity = lowerFluidVelocity +                                                       &
            0.5_wp * velocityDifference *                                                    &
            (1.0_wp + tanh(2.0_wp * real(grid%coordinates(i,2), wp)))
       temperature = (lowerFluidTemperature * (upperFluidVelocity - velocity) +              &
            upperFluidTemperature * (velocity - lowerFluidVelocity)) /                       &
            (upperFluidVelocity - lowerFluidVelocity) +                                      &
            0.5_wp * (upperFluidVelocity - velocity) * (velocity - lowerFluidVelocity)

       state%conservedVariables(i,1) =                                                       &
            1.0_wp / ((ratioOfSpecificHeats - 1.0_wp) * temperature)
       state%conservedVariables(i,2) = state%conservedVariables(i,1) * velocity
       state%conservedVariables(i,3:nDimensions+1) = 0.0_wp
       state%conservedVariables(i,nDimensions+2) =                                           &
            state%conservedVariables(i,1) * temperature / ratioOfSpecificHeats +             &
            0.5_wp * state%conservedVariables(i,1) * velocity ** 2

    end do

  end subroutine shearLayerInitialCondition

end program shear_layer
