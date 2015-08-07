#include "config.h"

program mixing_layer

  use MPI
  use, intrinsic :: iso_fortran_env

  use Grid_enum
  use State_enum

  use Region_mod, only : t_Region

  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush, gracefulExit
  use PLOT3DHelper, only : plot3dDetectFormat

  !> Generates the grid, BC, initial condition, target state, and adjoint mollifiers
  !> for a reactive spatial mixing layer.

  implicit none

  integer :: i, numProcs, ierror
  integer :: imin_sponge, imax_sponge, jmin_sponge, jmax_sponge
  character(len = STRING_LENGTH) :: inputname, filename
  type(t_Region) :: region
  integer, allocatable :: globalGridSizes(:,:)

  print *
  print *, '! ================================================= !'
  print *, '!                                                   !'
  print *, '!    2D MIXING LAYER GENERATOR                      !'
  print *, '!    Creates:                                       !'
  print *, '!               - Cartesian grid                    !'
  print *, '!               - Target/initial solution           !'
  print *, '!               - Boundary conditions               !'
  print *, '!    for a spatially evolving mixing layer          !'
  print *, '!                                                   !'
  print *, '! ================================================= !'
  print *

  ! Initialize MPI.
  call MPI_Init(ierror)

  ! Exit if executed in parallel
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)
  if (numProcs.gt.1) then
     print *, 'mixing_layer.f90 only implemented in serial for now...'
     stop
  end if

  ! Parse options from the input file.
  inputname = "magudi.inp"
  call parseInputFile(inputname)

  ! Generate the grid
  call mixingLayerGrid(imin_sponge,imax_sponge,jmin_sponge,jmax_sponge)

  ! Save the grid
  call getRequiredOption("grid_file", filename)
  call region%saveData(QOI_GRID, filename)

  ! Compute normalized metrics, norm matrix and Jacobian.
  do i = 1, size(region%grids)
     call region%grids(i)%update()
  end do
  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  ! Write out some useful information.
  call region%reportGridDiagnostics()

  ! Generate the BC
  call mixingLayerBC(imin_sponge,imax_sponge,jmin_sponge,jmax_sponge)

  ! Generate the initial condition and target state.
  do i = 1, size(region%grids)
     call mixingLayerInitialCondition(region%states(i), region%grids(i),                     &
          region%simulationFlags%useTargetState)
  end do

  ! Save initial condition.
  call getRequiredOption("grid_file", filename)
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
  subroutine mixingLayerGrid(i1,i2,j1,j2)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod

    ! <<< Internal modules >>>
    use InputHelper, only : getOption
    use ErrorHandler, only : gracefulExit

    ! <<< Arguments >>>
    integer, intent(out) :: i1, i2, j1, j2

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k
    integer :: nx, ny, nz, nx_, ny_, nz_
    real(wp) :: xmini, xmaxi, ymini, ymaxi
    real(wp) :: xmino, xmaxo, ymino, ymaxo
    real(wp) :: zmin, zmax
    real(wp) :: b, c, sigma, min_meshsize, max_meshsize, y1, y2
    real(wp), allocatable, dimension(:) :: s, g
    logical :: stretch_x, stretch_y

    ! Read in grid size and dimensions.
    nx = getOption("nx", 1025)
    ny = getOption("ny", 513)
    nz = getOption("nz", 1)
    xmini = getOption("xmin_interior", 0.0_wp)
    xmaxi = getOption("xmax_interior", 120.0_wp)
    ymini = getOption("ymin_interior", -30.0_wp)
    ymaxi = getOption("ymax_interior", 30.0_wp)
    xmino = getOption("xmin_outer", -40.0_wp)
    xmaxo = getOption("xmax_outer", 160.0_wp)
    ymino = getOption("ymin_outer", -50.0_wp)
    ymaxo = getOption("ymax_outer", 50.0_wp)
    zmin = getOption("zmin", 0.0_wp)
    zmax = getOption("zmax", 0.0_wp)

    ! Allocate the global grid size and assign values.
    ! `globalGridSizes(i,j)` is the number of grid points on grid `j` along dimension `i`.
    if (nz.eq.1) then
       allocate(globalGridSizes(2,1))
    else
       allocate(globalGridSizes(3,1))
       globalGridSizes(3,1) = nz
    end if
    globalGridSizes(1,1) = nx
    globalGridSizes(2,1) = ny

    ! Setup the region.
    call region%setup(MPI_COMM_WORLD, globalGridSizes)

    ! Local mesh
    nx_ = region%grids(1)%localSize(1)
    ny_ = region%grids(1)%localSize(2)
    nz_ = region%grids(1)%localSize(3)

    ! Should we stretch the mesh?
    stretch_x = getOption('stretch_x',.false.)
    stretch_y = getOption('stretch_y',.false.)

    ! Generate the grid.
    do k = 1, nz_
       do j = 1, ny_
          do i = 1, nx_
             ! Create X
             if (.not.stretch_x) region%grids(1)%coordinates(i+nx_*(j-1+ny_*(k-1)),1) =      &
                  (xmaxo - xmino)*real(i-1,wp)/real(nx-1) + xmino

             ! Create Y
             if (.not.stretch_y) region%grids(1)%coordinates(i+nx_*(j-1+ny_*(k-1)),2) =      &
                     (ymaxo - ymino)*real(j-1,wp)/real(ny-1) + ymino

             ! Create Z
             if (nz.ne.1) region%grids(1)%coordinates(i+nx_*(j-1+ny_*(k-1)),3) =             &
                  (zmax - zmin)*real(k-1,wp)/real(nz-1) + zmin
          end do
       end do
    end do

    ! Stretch the grid.
    if (stretch_x) then
       ! Parameters
       sigma=0.07_wp
       b=80.0_wp
       c=0.58_wp

       ! Create uniform spacing.
       allocate(s(nx_))
       do i = 1, nx_
          s(i) = real(i-1,wp) / (nx_ - 1)
       end do

       ! Compute g(s).
       allocate(g(nx_))
       call mapping_function(s, b, c, sigma, g)

       ! Find min/max spacing.
       min_meshsize=huge(1.0_wp)
       max_meshsize=-huge(1.0_wp)

       ! Compute x.
       do k = 1, nz_
          do j = 1, ny_
             do i = 1, nx_
                ! Create x.
                region%grids(1)%coordinates(i + nx_ * (j - 1 + ny_ * (k-1)), 1)              &
                     = xmino + 0.5_wp * (xmaxo - xmino) * (1.0_wp + g(i))

                ! Find min/max spacing.
                if (i.gt.1) then
                   y1 = region%grids(1)%coordinates(i + nx_ * (j - 1 + ny_ * (k - 1)), 1)
                   y2 = region%grids(1)%coordinates(i - 1 + nx_ * (j - 1 + ny_ * (k - 1)), 1)
                   min_meshsize=min(min_meshsize, abs(y2 - y1))
                   max_meshsize=max(max_meshsize, abs(y2 - y1))
                end if
             end do
          end do
       end do
       print *
       print *, 'min/max x-spacing:',min_meshsize,max_meshsize
       print *

       deallocate(s, g)
    end if

    if (stretch_y) then
       ! Parameters
       sigma=0.2_wp
       b=20.0_wp
       c=0.62_wp

       ! Create uniform spacing.
       allocate(s(ny_))
       do j = 1, ny_
          s(j) = real(j-1,wp) / (ny_ - 1)
       end do

       ! Compute g(s).
       allocate(g(ny_))
       call mapping_function(s, b, c, sigma, g)

       ! Find min/max spacing.
       min_meshsize=huge(1.0_wp)
       max_meshsize=-huge(1.0_wp)

       ! Compute y.
       do k = 1, nz_
          do j = 1, ny_
             do i = 1, nx_
                ! Create y.
                region%grids(1)%coordinates(i + nx_ * (j - 1 + ny_ * (k-1)), 2)              &
                     = ymino + 0.5_wp * (ymaxo - ymino) * (1.0_wp + g(j))

                ! Find min/max spacing.
                if (j.gt.1) then
                   y1 = region%grids(1)%coordinates(i + nx_ * (j - 1 + ny_ * (k - 1)), 2)
                   y2 = region%grids(1)%coordinates(i + nx_ * (j - 2 + ny_ * (k - 1)), 2)
                   min_meshsize=min(min_meshsize, abs(y2 - y1))
                   max_meshsize=max(max_meshsize, abs(y2 - y1))
                end if
             end do
          end do
       end do
       print *
       print *, 'min/max y-spacing:',min_meshsize,max_meshsize
       print *

       deallocate(s, g)
    end if

    ! Find extents of outer region.
    j=1; k=1;
    i1=1
    do i = 1, nx
       if (region%grids(1)%coordinates(i+nx_*(j-1+ny_*(k-1)),1) <= xmini) i1=i+1
    end do
    i2=nx
    do i = nx,1,-1
       if (region%grids(1)%coordinates(i+nx_*(j-1+ny_*(k-1)),1) >= xmaxi) i2=i
    end do
    i=1; k=1;
    j1=1
    do j = 1, ny
       if (region%grids(1)%coordinates(i+nx_*(j-1+ny_*(k-1)),2) <= ymini) j1=j+1
    end do
    j2=ny
    do j = ny,1,-1
       if (region%grids(1)%coordinates(i+nx_*(j-1+ny_*(k-1)),2) >= ymaxi) j2=j
    end do

    print *
    print *, 'Interior domain extents: [',i1,',',i2,'] x [',j1,',',j2,']'
    print*

    return
  end subroutine mixingLayerGrid


  ! ======================== !
  ! Grid stretching function !
  ! ======================== !
  Subroutine mapping_function(s, b, c, sigma, g)
    implicit none

    integer, parameter :: wp = SCALAR_KIND
    real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
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


  ! =================== !
  ! Boundary conditions !
  ! =================== !
  subroutine mixingLayerBC(imin_sponge,imax_sponge,jmin_sponge,jmax_sponge)
    implicit none

    integer, intent(in), optional :: imin_sponge,imax_sponge,jmin_sponge,jmax_sponge
    integer :: i, bc, nbc, iunit
    integer, allocatable, dimension(:) :: grid,normDir,imin,imax,jmin,jmax,kmin,kmax
    character(len = 22), allocatable, dimension(:) :: name,type

    ! Number of BC
    nbc = 11

    ! Allocate BC
    allocate(name(nbc),type(nbc),grid(nbc),normDir(nbc),&
         imin(nbc),imax(nbc),jmin(nbc),jmax(nbc),kmin(nbc),kmax(nbc))

    ! Set the BC
    ! GRID 1
    grid(:)    = 1

    ! BC 1
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

    ! BC 2
    bc = 2
    name   (bc) = 'inflowSponge'
    type   (bc) = 'SPONGE'
    normDir(bc) =  1
    imin   (bc) =  1
    imax   (bc) =  imin_sponge
    jmin   (bc) =  1
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! BC 3
    bc = 3
    name   (bc) = 'outflow'
    type   (bc) = 'SAT_FAR_FIELD'
    normDir(bc) = -1
    imin   (bc) = -1
    imax   (bc) = -1
    jmin   (bc) =  1
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! BC 4
    bc = 4
    name   (bc) = 'outflowSponge'
    type   (bc) = 'SPONGE'
    normDir(bc) = -1
    imin   (bc) =  imax_sponge
    imax   (bc) = -1
    jmin   (bc) =  1
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! BC 5
    bc = 5
    name   (bc) = 'bottom'
    type   (bc) = 'SAT_FAR_FIELD'
    normDir(bc) =  2
    imin   (bc) =  1
    imax   (bc) = -1
    jmin   (bc) =  1
    jmax   (bc) =  1
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! BC 6
    bc = 6
    name   (bc) = 'bottomSponge'
    type   (bc) = 'SPONGE'
    normDir(bc) =  2
    imin   (bc) =  1
    imax   (bc) = -1
    jmin   (bc) =  1
    jmax   (bc) =  jmin_sponge
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! BC 7
    bc = 7
    name   (bc) = 'top'
    type   (bc) = 'SAT_FAR_FIELD'
    normDir(bc) = -2
    imin   (bc) =  1
    imax   (bc) = -1
    jmin   (bc) = -1
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! BC 8
    bc = 8
    name   (bc) = 'topSponge'
    type   (bc) = 'SPONGE'
    normDir(bc) = -2
    imin   (bc) =  1
    imax   (bc) = -1
    jmin   (bc) =  jmax_sponge
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! BC 9
    bc = 9
    name   (bc) = 'excitationSupport'
    type   (bc) = 'SOLENOIDAL_EXCITATION'
    normDir(bc) =  0
    imin   (bc) =  1
    imax   (bc) =  imin_sponge
    jmin   (bc) =  1
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! BC 10
    bc = 10
    name   (bc) = 'targetRegion'
    type   (bc) = 'COST_TARGET'
    normDir(bc) =  0
    imin   (bc) =  1
    imax   (bc) = -1
    jmin   (bc) =  1
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! BC 11
    bc = 11
    name   (bc) = 'controlRegion'
    type   (bc) = 'ACTUATOR'
    normDir(bc) =  0
    imin   (bc) =  1
    imax   (bc) = -1
    jmin   (bc) =  1
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

   ! Open the file
    iunit=11
    open(iunit,file="bc.dat")
    
    write (*,'(A)') ''
    write (*,'(A)') 'Writing boundary conditions'

    ! Write the header
    write(iunit,'(1a87)') "# Name                 Type                  Grid normDir iMin iMax jMin jMax kMin kMax"
    write(iunit,'(1a87)') "# ==================== ===================== ==== ======= ==== ==== ==== ==== ==== ===="

    ! Write the BC
    do i=1,nbc
       write(iunit,'(a2,2a21,8I5)') '  ',adjustl(name(i)),adjustl(type(i)),&
            grid(i),normDir(i),imin(i),imax(i),jmin(i),jmax(i),kmin(i),kmax(i)
    end do

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
    integer :: i, nSpecies, H2, O2, nDimensions, ierror
    real(wp) :: ratioOfSpecificHeats, upperVelocity, lowerVelocity,                          &
         density, velocity, temperature, Z, fuel, oxidizer, YF0, YO0, Z0,                    &
         x, y, xmini, xmaxi, ymini, ymaxi, xmaxo, xmino, ymaxo, ymino
    real(wp), parameter :: growthRate = 0.1_wp

    generateTargetState_ = .false.
    if (present(generateTargetState)) generateTargetState_ = generateTargetState

    call MPI_Cartdim_get(grid%comm, nDimensions, ierror)

    ! Species
    nSpecies = getOption("number_of_species", 0)

    ! Only implemented for 2 species
    if (nspecies.gt.2) then
       print *, 'WARNING, max of 2 species (for now)'
       stop
    end if

    ! Mixing layer velocities.
    upperVelocity = getOption("upper_velocity", 0.0_wp)
    lowerVelocity = getOption("lower_velocity", 0.0_wp)

    ! Domain extents
    xmini = getOption("xmin_interior", 0.0_wp)
    xmaxi = getOption("xmax_interior", 120.0_wp)
    ymini = getOption("ymin_interior", -30.0_wp)
    ymaxi = getOption("ymax_interior", 30.0_wp)
    xmino = getOption("xmin_outer", -40.0_wp)
    xmaxo = getOption("xmax_outer", 160.0_wp)
    ymino = getOption("ymin_outer", -50.0_wp)
    ymaxo = getOption("ymax_outer", 50.0_wp)

    ! Species indeces.
    H2 = nDimensions+2+1
    O2 = H2 + 1

    ! Mixture properties.
    Z0 = getOption("initial_mixture_fraction", 1.0_wp)
    call getRequiredOption("initial_fuel_mass_fraction", Yf0)
    call getRequiredOption("initial_oxidizer_mass_fraction", Yo0)

    ! Specific heat ration (gamma).
    ratioOfSpecificHeats = getOption("ratio_of_specific_heats", 1.4_wp)

    do i = 1, grid%nGridPoints

       ! Local coordinates.
       x = grid%coordinates(i,1)
       y = grid%coordinates(i,2)

       ! Mixture fraction.
       Z = 0.5_wp * Z0 * ( 1.0_wp - erf(y / (1.0_wp +                                        &
            growthRate * max(0.0_wp, min(x, xmaxi) - xmini))) )

       ! Mixture components.
       fuel = YF0 * Z
       oxidizer = YO0 * (1.0_wp - Z)

       ! Species decay at the outflow.
       fuel = fuel - fuel * max( 0.0_wp, (x - xmaxi) / (xmaxo - xmaxi) )
       oxidizer = oxidizer - oxidizer * max( 0.0_wp, (x - xmaxi) / (xmaxo - xmaxi) )

       ! Density.
       if (nspecies.gt.0) then
          density = 1.0_wp
       else
          density = 1.0_wp
       end if

       ! Velocity.
       velocity = lowerVelocity +                                                            &
            0.5_wp * (upperVelocity - lowerVelocity) * (1.0_wp + tanh(2.0_wp *               &
            grid%coordinates(i,2) /                                                          &
            (1.0_wp + growthRate*max(0.0_wp, grid%coordinates(i,1)-xmini))))

       ! Temperature.
       temperature =  1.0_wp / (ratioOfSpecificHeats - 1.0_wp)

       ! State variables.
       state%conservedVariables(i,1) = density
       state%conservedVariables(i,2) = state%conservedVariables(i,1) * velocity
       state%conservedVariables(i,3:nDimensions+1) = 0.0_wp
       state%conservedVariables(i,nDimensions+2) =                                           &
            state%conservedVariables(i,1) * temperature / ratioOfSpecificHeats +             &
            0.5_wp * state%conservedVariables(i,1) * velocity ** 2
       if (nSpecies.gt.0) state%conservedVariables(i,H2) = fuel *                            &
            state%conservedVariables(i,1)
       if (nSpecies.gt.1) state%conservedVariables(i,O2) = oxidizer *                        &
            state%conservedVariables(i,1)

       ! Target solution.
       if (generateTargetState_) state%targetState(i,:) = state%conservedVariables(i,:)

    end do

    return
  end subroutine mixingLayerInitialCondition


  ! ================ !
  ! Target mollifier !
  ! ================ !
  subroutine mixingLayerTargetMollifier(state, grid, imin, imax, jmin, jmax)

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

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, gridIndex
    integer, intent(out) :: imin, imax, jmin, jmax
    integer :: nx, ny, nz
    real(wp) :: xmin, xmax, ymin, ymax
    real(wp) :: x
    real(wp), dimension(:,:,:,:), allocatable :: mollifier
    real(wp), parameter :: s = 10.0_wp, r = 0.2_wp, eps = 1.0E-6_wp

    ! Make sure target mollifier is allocated.
    if (region%simulationFlags%predictionOnly) then
       print *, 'WARNING: target mollifier requires disable_adjoint_solver = false'
       stop
    end if

    ! Read the mollifier extents.
    call getRequiredOption("targer_mollifier_xmin", xmin)
    call getRequiredOption("targer_mollifier_xmax", xmax)
    call getRequiredOption("targer_mollifier_ymin", ymin)
    call getRequiredOption("targer_mollifier_ymax", ymax)

    ! Get domain size
    nx = grid%globalSize(1)
    ny = grid%globalSize(2)
    if (size(grid%globalSize)>2) nz = grid%globalSize(3)

    ! Initialize the target mollifier.
    allocate(mollifier(nx,ny,nz,2))
    mollifier = 1.0_wp

    ! Hyperbolic tangent profile in y.
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             gridIndex = i+nx*(j-1+ny*(k-1))
             x = 2.0_wp * (grid%coordinates(gridIndex,2) - ymin) / (ymax - ymin) - 1.0_wp
             mollifier(i,j,k,1) = tanh(s * (x + 1.0_wp - 0.5_wp*r)) -                 &
                  tanh(s * (x - 1.0_wp + 0.5_wp*r))
          end do
       end do
    end do
    mollifier(:,:,:,1) = 0.5_wp * (mollifier(:,:,:,1) - minval(mollifier(:,:,:,1)))

    ! Hyperbolic tangent profile in x.
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             gridIndex = i+nx*(j-1+ny*(k-1))
             x = 2.0_wp * (grid%coordinates(gridIndex,1) - xmin) / (xmax - xmin) - 1.0_wp
             mollifier(i,j,k,2) = tanh(s * (x + 1.0_wp - 0.5_wp*r)) -                         &
                  tanh(s * (x - 1.0_wp + 0.5_wp*r))
          end do
       end do
    end do
    mollifier(:,:,:,2) = 0.5_wp * (mollifier(:,:,:,2) - minval(mollifier(:,:,:,2)))

    ! Transfer to Magudi.
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             gridIndex = i+nx*(j-1+ny*(k-1))
             grid%targetMollifier(gridIndex,1) = mollifier(i,j,k,1) * mollifier(i,j,k,2)
          end do
       end do
    end do

    ! Find initial extents in y.
    jmin = 1; jmax = 1
    i = 1; k = 1
    do j = 1, ny
       gridIndex = i+nx*(j-1+ny*(k-1))
       if (grid%coordinates(gridIndex,2) <= ymin) jmin = j
       if (grid%coordinates(gridIndex,2) < ymax) jmax = j+1
    end do

    ! Find new extents in x.
    j = int(0.5_wp*(jmin+jmax)); k = 1
    do i = 1, nx
       if (mollifier(i,j,k,1) * mollifier(i,j,k,2) < eps) then
          imin = i
       else
          exit
       end if
    end do
    do i = imin+1, nx
       if (mollifier(i,j,k,1) * mollifier(i,j,k,2) < eps) then
          imax = i
          exit
       end if
    end do

    ! Find new extents in y.
    i = int(0.5_wp*(imin+imax)); k = 1
    do j = 1, ny
       if (mollifier(i,j,k,1) * mollifier(i,j,k,2) < eps) then
          jmin = j
       else
          exit
       end if
    end do
    do j = jmin+1, ny
       if (mollifier(i,j,k,1) * mollifier(i,j,k,2) < eps) then
          jmax = j
          exit
       end if
    end do

    print *
    print *, 'Target mollifier extents: [',imin,',',imax,'] x [',jmin,',',jmax,']'
    print *

    ! Clean up.
    deallocate(mollifier)

    return
  end subroutine mixingLayerTargetMollifier


  ! ================= !
  ! Control mollifier !
  ! ================= !
  subroutine mixingLayerControlMollifier(state, grid, imin, imax, jmin, jmax)

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

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, gridIndex
    integer, intent(out) :: imin, imax, jmin, jmax
    integer :: nx, ny, nz
    real(wp) :: xmin, xmax, ymin, ymax
    real(wp) :: x
    real(wp), dimension(:,:,:,:), allocatable :: mollifier
    real(wp), parameter :: s = 10.0_wp, r = 0.2_wp, eps = 1.0e-6_wp

    ! Make sure control mollifier is allocated.
    if (region%simulationFlags%predictionOnly) then
       print *, 'WARNING: control mollifier requires disable_adjoint_solver = false'
       stop
    end if

    ! Read the mollifier extents.
    call getRequiredOption("control_mollifier_xmin", xmin)
    call getRequiredOption("control_mollifier_xmax", xmax)
    call getRequiredOption("control_mollifier_ymin", ymin)
    call getRequiredOption("control_mollifier_ymax", ymax)

    ! Get domain size
    nx = grid%globalSize(1)
    ny = grid%globalSize(2)
    if (size(grid%globalSize)>2) nz = grid%globalSize(3)

    ! Initialize the control mollifier.
    allocate(mollifier(nx,ny,nz,2))
    mollifier = 1.0_wp

    ! Hyperbolic tangent profile in y.
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             gridIndex = i+nx*(j-1+ny*(k-1))
             x = 2.0_wp * (grid%coordinates(gridIndex,2) - ymin) / (ymax - ymin) - 1.0_wp
             mollifier(i,j,k,1) = tanh(s * (x + 1.0_wp - 0.5_wp*r)) -                 &
                  tanh(s * (x - 1.0_wp + 0.5_wp*r))
          end do
       end do
    end do
    mollifier(:,:,:,1) = 0.5_wp * (mollifier(:,:,:,1) - minval(mollifier(:,:,:,1)))

    ! Hyperbolic tangent profile in x.
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             gridIndex = i+nx*(j-1+ny*(k-1))
             x = 2.0_wp * (grid%coordinates(gridIndex,1) - xmin) / (xmax - xmin) - 1.0_wp
             mollifier(i,j,k,2) = tanh(s * (x + 1.0_wp - 0.5_wp*r)) -                         &
                  tanh(s * (x - 1.0_wp + 0.5_wp*r))
          end do
       end do
    end do
    mollifier(:,:,:,2) = 0.5_wp * (mollifier(:,:,:,2) - minval(mollifier(:,:,:,2)))

    ! Transfer to Magudi.
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             gridIndex = i+nx*(j-1+ny*(k-1))
             grid%controlMollifier(gridIndex,1) = mollifier(i,j,k,1) * mollifier(i,j,k,2)
          end do
       end do
    end do

    ! Find initial extents in y.
    jmin = 1; jmax = 1
    i = 1; k = 1
    do j = 1, ny
       gridIndex = i+nx*(j-1+ny*(k-1))
       if (grid%coordinates(gridIndex,2) <= ymin) jmin = j
       if (grid%coordinates(gridIndex,2) < ymax) jmax = j+1
    end do

    ! Find new extents in x.
    j = int(0.5_wp*(jmin+jmax)); k = 1
    do i = 1, nx
       if (mollifier(i,j,k,1) * mollifier(i,j,k,2) < eps) then
          imin = i
       else
          exit
       end if
    end do
    do i = imin+1, nx
       if (mollifier(i,j,k,1) * mollifier(i,j,k,2) < eps) then
          imax = i
          exit
       end if
    end do

    ! Find new extents in y.
    i = int(0.5_wp*(imin+imax)); k = 1
    do j = 1, ny
       if (mollifier(i,j,k,1) * mollifier(i,j,k,2) < eps) then
          jmin = j
       else
          exit
       end if
    end do
    do j = jmin+1, ny
       if (mollifier(i,j,k,1) * mollifier(i,j,k,2) < eps) then
          jmax = j
          exit
       end if
    end do

    print *
    print *, 'Control mollifier extents: [',imin,',',imax,'] x [',jmin,',',jmax,']'
    print *

    ! Clean up.
    deallocate(mollifier)

    return
  end subroutine mixingLayerControlMollifier

end program mixing_layer
