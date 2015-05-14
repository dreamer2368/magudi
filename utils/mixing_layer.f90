#include "config.h"

program mixing_layer

  use MPI
  use, intrinsic :: iso_fortran_env

  use Grid_enum
  use State_enum

  use Region_mod, only : t_Region

  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush, gracefulExit
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage

  !> Generates the grid, BC, initial condition, and target state for a spatial mixing layer.

  implicit none

  integer :: i, ierror
  integer :: i1, i2, j1, j2
  character(len = STRING_LENGTH) :: filename
  logical :: success
  type(t_Region) :: region
  integer, allocatable :: globalGridSizes(:,:)

  print *, '! ================================================= !'
  print *, '!                                                   !'
  print *, '!    2D MIXING LAYER GENERATOR                      !'
  print *, '!    Creates:                                       !'
  print *, '!               - Cartesian grid                    !'
  print *, '!               - Target/initial solution           !'
  print *, '!               - Boundary conditions               !'
  print *, '!    for a spatially-evolving mixing layer          !'
  print *, '!                                                   !'
  print *, '! ================================================= !'

  ! Initialize MPI.
  call MPI_Init(ierror)

  ! Parse options from the input file.
  filename = "mixing_layer.inp"
  call parseInputFile(filename)

  ! Generate the grid
  call mixingLayerGrid(i1,i2,j1,j2,filename)

  ! Generate the BC
  call mixingLayerBC(i1,i2,j1,j2)

  ! Generate the initial condition and target state.
  do i = 1, size(region%grids)
     call mixingLayerInitialCondition(region%states(i), region%grids(i),                      &
          region%simulationFlags%useTargetState)
  end do

  ! Save initial condition.
  i = len_trim(filename)
  if (filename(i-3:i) == ".xyz") then
     filename = filename(:i-4) // ".ic.q"
  else
     filename = PROJECT_NAME // ".ic.q"
  end if
  call region%saveData(QOI_FORWARD_STATE, filename)

  ! Save target state.
  if (region%simulationFlags%useTargetState) then
     i = len_trim(filename)
     if (filename(i-4:i) == ".ic.q") then
        filename = filename(:i-5) // ".target.q"
     else
        filename = PROJECT_NAME // ".target.q"
     end if
     call region%saveData(QOI_TARGET_STATE, filename)
  end if

  ! Cleanup.
  call region%cleanup()

  ! Finalize MPI.
  call MPI_Finalize(ierror)

contains

  ! =============== !
  ! Grid generation !
  ! =============== !
  subroutine mixingLayerGrid(i1,i2,j1,j2,filename)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod

    ! <<< Internal modules >>>
    use InputHelper, only : getOption
    use ErrorHandler, only : gracefulExit

    use InputHelper, only : parseInputFile, getOption, getRequiredOption
    use ErrorHandler, only : writeAndFlush, gracefulExit
    use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage

    ! <<< Arguments >>>
    type(t_State) :: state
    type(t_Grid) :: grid
    integer, intent(out) :: i1, i2, j1, j2
    character(len = STRING_LENGTH), intent(inout) :: filename

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    logical :: generateTargetState_
    integer :: i, j, k, n, nDimensions, ierror
    integer :: nx, ny, nz
    real(wp) :: xmini, xmaxi, ymini, ymaxi
    real(wp) :: xmino, xmaxo, ymino, ymaxo
    real(wp) :: g0, b, c, sig, dy_min, dy_max
    real(wp), allocatable, dimension(:) :: s, g
    logical :: stretch_y

    ! Read in grid size and dimensions.
    nx = getOption("nx", 1025)
    ny = getOption("ny", 513)
    nz = 1
    xmini = getOption("xmin_interior", 0.0_wp)
    xmaxi = getOption("xmax_interior", 120.0_wp)
    xmini = getOption("ymin_interior", -30.0_wp)
    xmaxi = getOption("ymax_interior", 30.0_wp)
    xmino = getOption("xmin_outer", -40.0_wp)
    xmaxo = getOption("xmax_outer", 160.0_wp)
    xmino = getOption("ymin_outer", -50.0_wp)
    xmaxo = getOption("ymax_outer", 50.0_wp)

    ! Allocate the global grid size and assign values.
    ! `globalGridSizes(i,j)` is the number of grid points on grid `j` along dimension `i`.
    allocate(globalGridSizes(2,1))
    globalGridSizes(1,1) = nx
    globalGridSizes(2,1) = ny

    ! Setup the region.
    call region%setup(MPI_COMM_WORLD, globalGridSizes)

    ! Should we stretch the mesh?
    stretch_y = getOption('stretch_y',.false.)

    ! Generate the grid.
    do k = 1, grid%localSize(3)
       do j = 1, grid%localSize(2)
          do i = 1, grid%localSize(1)
             ! Create X
             grid%coordinates(i + grid%localSize(1) * (j - 1 +                               &
                  grid%localSize(2) * (k - 1)), 1) = (xmaxo-xmino)*real(i-1,wp)/real(grid%globalSize(1)-1) + xmino

             ! Create Y
             if (.not.stretch_y) then
                grid%coordinates(i + grid%localSize(1) * (j - 1 +                               &
                  grid%localSize(2) * (k - 1)), 2) = (ymaxo-ymino)*real(j-1,wp)/real(grid%globalSize(2)-1) + ymino
             end if

             ! Create Z
             grid%coordinates(i + grid%localSize(1) * (j - 1 +                               &
                  grid%localSize(2) * (k - 1)), 3) = 0.0_wp
          end do
       end do
    end do

    ! Stretch the grid.
    if (stretch_y) then
       ! Parameters
       sig=0.2_wp
       b=20.0_wp
       c=0.62_wp
       n=100

       ! Create uniform spacing.
       allocate(s(grid%localSize(2)))
       do j = 1, grid%localSize(2)
          s(j) = real(j-1,wp)/real(grid%localSize(2)-1,wp)
       end do

       ! Compute g(s).
       allocate(g(grid%localSize(2)))
       call g_of_s(1.0_wp,b,c,sig,n,g0)
       do j=1,grid%localSize(2)
          call g_of_s(s(j),b,c,sig,n,g(j))
       end do
       ! Compute y.
       do k = 1, grid%localSize(3)
          do j = 1, grid%localSize(2)
             do i = 1, grid%localSize(1)
                ! Create X
                grid%coordinates(i + grid%localSize(1) * (j - 1 +                               &
                     grid%localSize(2) * (k - 1)), 2) = 0.5_wp*(ymaxo-ymino)*g(j)/g0
             end do
          end do
       end do
    end if

    ! Save the grid
    call region%saveData(QOI_GRID, filename)

    ! Compute normalized metrics, norm matrix and Jacobian.
    do i = 1, size(region%grids)
       call region%grids(i)%update()
    end do
    call MPI_Barrier(MPI_COMM_WORLD, ierror)

    ! Write out some useful information.
    call region%reportGridDiagnostics()

!!$    ! Find min/max spacing.
!!$    dy_min=huge(1.0_wp)
!!$    dy_max=-huge(1.0_wp)
!!$    do j=1,
!!$       dy_min=min(dy_min,X(1,1,J+1,1,2)-X(1,1,J,1,2))
!!$       dy_max=max(dy_max,X(1,1,J+1,1,2)-X(1,1,J,1,2))
!!$    end do
!!$    print *
!!$    print *, 'min/max y-spacing:',dy_min,dy_max
!!$    print *
!!$
!!$    ! Find extents of outer region.
!!$    i1=1
!!$    do i = 1, nx
!!$       if (X(1,I,1,1,1) <= xmini) i1=i
!!$    End Do
!!$    i2=nx
!!$    Do i = nx,1,-1
!!$       if (X(1,I,1,1,1) >= xmaxi) i2=i
!!$    end do
!!$    j1=1
!!$    do j = 1, ny
!!$       if (X(1,1,J,1,2) <= ymini) j1=j
!!$    end do
!!$    j2=ny
!!$    do j = ny,1,-1
!!$       if (X(1,1,J,1,2) >= ymaxi) j2=j
!!$    end do
!!$
!!$    print *, 'Extents:',i1,i2,j1,j2

    return
  end subroutine mixingLayerGrid


  ! ======================== !
  ! Grid stretching function !
  ! ======================== !
  Subroutine g_of_s(s,b,c,sig,n,g)
    implicit none

    integer, parameter :: wp = SCALAR_KIND
    integer, intent(in) :: n
    real(wp), intent(in) :: s,b,c,sig
    real(wp), intent(out) :: g

    integer :: i
    real(wp) :: dx1,dx2,int1,int2
    real(wp), dimension(n) :: x1,x2

    ! Compute stretched grid for mixing layer.
    do i = 1, N
       x1(i) = ((s-0.5_8)/sig)*real(i-1,wp)/real(n-1,wp) - c/sig
       x2(i) = ((s-0.5_8)/sig)*real(i-1,wp)/real(n-1,wp) + c/sig
    end do
    dx1=x1(2)-x1(1)
    dx2=x2(2)-x2(1)

    int1=0.0_8
    int2=0.0_8
    do i=1,n
       int1 = int1 + erf(x1(i))*dx1
       int2 = int2 + erf(x2(i))*dx2
    end do

    g = (s-0.5_8)*(1.0_8+2.0_8*b)+b*sig*(int1-int2)

    return
  end subroutine g_of_s


  ! =================== !
  ! Boundary conditions !
  ! =================== !
  subroutine mixingLayerBC(i1,i2,j1,j2)
    implicit none

    integer, intent(in) :: i1, i2, j1, j2
    integer :: iunit

    ! Open the file
    iunit=11
    open(iunit,file="bc.dat")
    
    write (*,'(A)') ''
    write (*,'(A)') 'Writing boundary conditions'

    ! Write the header
    write(iunit,'(1a87)') "# Name                 Type                  Grid normDir iMin iMax jMin jMax kMin kMax"
    write(iunit,'(1a87)') "# ==================== ===================== ==== ======= ==== ==== ==== ==== ==== ===="

    ! Input the boundary conditions
    write(iunit,'(2a22,8I6)') 'inflow',            'SAT_FAR_FIELD',       1,    1,     1,   1,   1,  -1,   1,  -1
    write(iunit,'(2a22,8I6)') 'inflowSponge',      'SPONGE',              1,    1,     1,  i1,   1,  -1,   1,  -1
    write(iunit,'(2a22,8I6)') 'outflow',           'SAT_FAR_FIELD',       1,   -1,    -1,  -1,   1,  -1,   1,  -1
    write(iunit,'(2a22,8I6)') 'outflowSponge',     'SPONGE',              1,   -1,    12,  -1,   1,  -1,   1,  -1
    write(iunit,'(2a22,8I6)') 'bottom',            'SAT_FAR_FIELD',       1,    2,     1,  -1,   1,   1,   1,  -1
    write(iunit,'(2a22,8I6)') 'bottomSponge',      'SPONGE',              1,    2,     1,  -1,   1,  j1,   1,  -1
    write(iunit,'(2a22,8I6)') 'top',               'SAT_FAR_FIELD',       1,   -2,     1,  -1,  -1,  -1,   1,  -1
    write(iunit,'(2a22,8I6)') 'topSponge',         'SPONGE',              1,   -2,     1,  -1,  j2,  -1,   1,  -1

    ! Close the file
    close(iunit)

    return
  end subroutine mixingLayerBC


  ! ================== !
  ! Initial conditions !
  ! ================== !
  subroutine mixingLayerInitialCondition(state, grid, generateTargetState)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod

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

    return
  end subroutine mixingLayerInitialCondition

end program mixing_layer
