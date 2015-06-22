#include "config.h"

program counterflow_diffusion

  use MPI
  use, intrinsic :: iso_fortran_env

  use Grid_enum
  use State_enum

  use Region_mod, only : t_Region

  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush, gracefulExit
  use PLOT3DHelper, only : plot3dDetectFormat

  !> Generates the grid, BC, initial condition, and target state
  !> for a reactive spatial mixing layer.

  implicit none

  integer :: i, numProcs, ierror
  integer :: imin_sponge, imax_sponge, jmin_sponge, jmax_sponge
  character(len = STRING_LENGTH) :: inputname, filename
  type(t_Region) :: region
  integer, allocatable :: globalGridSizes(:,:)

  print *
  print *, '! ======================================= !'
  print *, '!                                         !'
  print *, '!    2D COUNTERFLOW DUFFUSION FLAME       !'
  print *, '!    Creates:                             !'
  print *, '!             - Cartesian grid            !'
  print *, '!             - Target/initial solution   !'
  print *, '!             - Boundary conditions       !'
  print *, '!    for a counterflow diffusion flame    !'
  print *, '!                                         !'
  print *, '! ======================================= !'
  print *

  ! Initialize MPI.
  call MPI_Init(ierror)

  ! Exit if executed in parallel
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)
  if (numProcs.gt.1) then
     print *, 'counterflow_diffusion.f90 only implemented in serial for now...'
     stop
  end if

  ! Parse options from the input file.
  inputname = "magudi.inp"
  call parseInputFile(inputname)

  ! Generate the grid
  call counterflowDiffusionGrid(imin_sponge, imax_sponge, jmin_sponge, jmax_sponge)

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
  call counterflowDiffusionBC(imin_sponge, imax_sponge, jmin_sponge, jmax_sponge)

  ! Generate the initial condition and target state.
  do i = 1, size(region%grids)
     call counterflowDiffusionInitialCondition(region%states(i), region%grids(i),                        &
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
  subroutine counterflowDiffusionGrid(i1, i2, j1, j2)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod

    ! <<< Internal modules >>>
    use InputHelper, only : getOption
    use ErrorHandler, only : gracefulExit

    implicit none

    ! <<< Arguments >>>
    integer, intent(out) :: i1, i2, j1, j2

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, n, gridIndex
    integer :: nx, ny, nz, nx_, ny_, nz_
    real(wp) :: xmini, xmaxi, ymini, ymaxi
    real(wp) :: xmino, xmaxo, ymino, ymaxo
    real(wp) :: zmin, zmax
    real(wp) :: dx, dy, dz
    real(wp) :: g0, b, c, sig, dy_min, dy_max, y1, y2
    real(wp), allocatable, dimension(:) :: s, g
    logical :: stretch_y

    ! Read in grid size and dimensions.
    call getRequiredOption("nx", nx)
    call getRequiredOption("ny", ny)
    nz = getOption("nz", 1)
    call getRequiredOption("xmin_interior", xmini)
    call getRequiredOption("xmax_interior", xmaxi)
    call getRequiredOption("ymin_interior", ymini)
    call getRequiredOption("ymax_interior", ymaxi)
    call getRequiredOption("xmin_outer", xmino)
    call getRequiredOption("xmax_outer", xmaxo)
    call getRequiredOption("ymin_outer", ymino)
    call getRequiredOption("ymax_outer", ymaxo)
    zmin = getOption("zmin", 0.0_wp)
    zmax = getOption("zmax", 0.0_wp)

    ! Should we stretch the mesh?
    stretch_y = getOption('stretch_y',.false.)

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

    ! Grid spacing
    dx = (xmaxo - xmino) / real(nx-1,wp)
    dy = (ymaxo - ymino) / real(ny-1,wp)
    dz = (zmax  - zmin ) / real(nz-1,wp)

    ! Generate the grid.
    do k = 1, nz_
       do j = 1, ny_
          do i = 1, nx_

             gridIndex = i+nx_*(j-1+ny_*(k-1))

             ! Create X
             region%grids(1)%coordinates(gridIndex,1) =                                      &
                  (xmaxo - xmino) * real(i-1,wp) / real(nx-1) + xmino


             ! Create Y
             if (.not.stretch_y) region%grids(1)%coordinates(gridIndex,2) =                  &
                  (ymaxo - ymino) * real(j-1,wp) / real(ny-1) + ymino

             ! Create Z
             if (nz.ne.1) region%grids(1)%coordinates(gridIndex,3) =                         &
                  (zmax - zmin - dz) * real(k-1,wp) / real(nz-1) + zmin

          end do
       end do
    end do

    ! Stretch the grid.
    if (stretch_y) then
       ! Parameters
       sig=0.18_wp
       b=20.0_wp
       c=0.6_wp
       n=100

       ! Create uniform spacing.
       allocate(s(ny_))
       do j = 1, ny_
          s(j) = real(j-1,wp)/real(ny-1,wp)
       end do

       ! Compute g(s).
       allocate(g(ny_))
       call g_of_s(1.0_wp,b,c,sig,n,g0)
       do j=1,ny_
          call g_of_s(s(j),b,c,sig,n,g(j))
       end do

       ! Find min/max spacing.
       dy_min=huge(1.0_wp)
       dy_max=-huge(1.0_wp)

       ! Compute y.
       do k = 1, nz_
          do j = 1, ny_
             do i = 1, nx_

                gridIndex = i+nx_*(j-1+ny_*(k-1))

                ! Create y.
                region%grids(1)%coordinates(gridIndex,2) = 0.5_wp*(ymaxo-ymino)*g(j)/g0

                ! Find min/max spacing.
                if (j.gt.1) then
                   y1 = region%grids(1)%coordinates(i+nx_*(j-1+ny_*(k-1)),2)
                   y2 = region%grids(1)%coordinates(i+nx_*(j-2+ny_*(k-1)),2)
                   dy_min=min(dy_min,abs(y2-y1))
                   dy_max=max(dy_max,abs(y2-y1))
                end if
             end do
          end do
       end do
       print *
       print *, 'min/max y-spacing:',dy_min,dy_max
       print *
    end if

    ! Find extents of outer region.
    j=1; k=1;
    i1=1
    do i = 1, nx
       if (region%grids(1)%coordinates(i+nx_*(j-1+ny_*(k-1)),1) <= xmini) i1=i
    end do
    i2=nx
    do i = nx,1,-1
       if (region%grids(1)%coordinates(i+nx_*(j-1+ny_*(k-1)),1) >= xmaxi) i2=i
    end do
    i=1; k=1;
    j1=1
    do j = 1, ny
       if (region%grids(1)%coordinates(i+nx_*(j-1+ny_*(k-1)),2) <= ymini) j1=j
    end do
    j2=ny
    do j = ny,1,-1
       if (region%grids(1)%coordinates(i+nx_*(j-1+ny_*(k-1)),2) >= ymaxi) j2=j
    end do

    print *
    print *, 'Interior domain extents: [',i1,',',i2,'] x [',j1,',',j2,']'
    print*

    return
  end subroutine counterflowDiffusionGrid

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
    do i=1,n
       x1(i) = ((s-0.5_wp)/sig)*real(i-1,wp)/real(n-1,wp) - c/sig
       x2(i) = ((s-0.5_wp)/sig)*real(i-1,wp)/real(n-1,wp) + c/sig
    end do
    dx1=x1(2)-x1(1)
    dx2=x2(2)-x2(1)

    int1=0.0_wp
    int2=0.0_wp
    do i=1,n
       int1 = int1 + erf(x1(i))*dx1
       int2 = int2 + erf(x2(i))*dx2
    end do

    g = (s-0.5_wp)*(1.0_wp+2.0_wp*b)+b*sig*(int1-int2)

    return
  end subroutine g_of_s

  ! =================== !
  ! Boundary conditions !
  ! =================== !
  subroutine counterflowDiffusionBC(imin_sponge, imax_sponge, jmin_sponge, jmax_sponge)

    ! <<< Internal modules >>>
    use InputHelper, only : getOption

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: imin_sponge, imax_sponge, jmin_sponge, jmax_sponge

    ! <<< Local variables >>>
    integer :: i, bc, nbc, iunit
    integer, allocatable, dimension(:) :: grid, normDir, imin, imax, jmin, jmax, kmin, kmax
    character(len = 22), allocatable, dimension(:) :: name, type
    character(len = STRING_LENGTH) :: entry

    ! Initialize BC
    nbc = 99
    allocate( name(nbc), type(nbc), grid(nbc), normDir(nbc),                                 &
         imin(nbc), imax(nbc), jmin(nbc), jmax(nbc), kmin(nbc), kmax(nbc) )

    ! Set the BC
    ! GRID 1
    grid(:) = 1
    bc = 0

    bc = bc + 1
    name   (bc) = 'farField.E'
    type   (bc) = 'SAT_FAR_FIELD'
    normDir(bc) =  1
    imin   (bc) =  1
    imax   (bc) =  1
    jmin   (bc) =  1
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    bc = bc + 1
    name   (bc) = 'farField.W'
    type   (bc) = 'SAT_FAR_FIELD'
    normDir(bc) = -1
    imin   (bc) = -1
    imax   (bc) = -1
    jmin   (bc) =  1
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    bc = bc + 1
    name   (bc) = 'sponge.E'
    type   (bc) = 'SPONGE'
    normDir(bc) =  1
    imin   (bc) =  1
    imax   (bc) =  imin_sponge
    jmin   (bc) =  1
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    bc = bc + 1
    name   (bc) = 'sponge.W'
    type   (bc) = 'SPONGE'
    normDir(bc) =  -1
    imin   (bc) =  imax_sponge
    imax   (bc) = -1
    jmin   (bc) =  1
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    bc = bc + 1
    name   (bc) = 'farField.S'
    type   (bc) = 'SAT_FAR_FIELD'
    normDir(bc) =  2
    imin   (bc) =  1
    imax   (bc) = -1
    jmin   (bc) =  1
    jmax   (bc) =  1
    kmin   (bc) =  1
    kmax   (bc) = -1

    bc = bc + 1
    name   (bc) = 'farField.N'
    type   (bc) = 'SAT_FAR_FIELD'
    normDir(bc) =  -2
    imin   (bc) =  1
    imax   (bc) = -1
    jmin   (bc) = -1
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    bc = bc + 1
    name   (bc) = 'sponge.S'
    type   (bc) = 'SPONGE'
    normDir(bc) =  2
    imin   (bc) =  1
    imax   (bc) = -1
    jmin   (bc) =  1
    jmax   (bc) = jmin_sponge
    kmin   (bc) =  1
    kmax   (bc) = -1

    bc = bc + 1
    name   (bc) = 'sponge.N'
    type   (bc) = 'SPONGE'
    normDir(bc) = -2
    imin   (bc) =  1
    imax   (bc) = -1
    jmin   (bc) =  jmax_sponge
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    entry = getOption("target_mollifier_file", "")
    if (len_trim(entry) > 0) then
       bc = bc + 1
       name   (bc) = 'targetRegion'
       type   (bc) = 'COST_TARGET'
       normDir(bc) =  0
       imin   (bc) =  43
       imax   (bc) =  471
       jmin   (bc) =  37
       jmax   (bc) =  221
       kmin   (bc) =  1
       kmax   (bc) = -1
    end if

    entry = getOption("control_mollifier_file", "")
    if (len_trim(entry) > 0) then
       bc = bc + 1
       name   (bc) = 'controlRegion'
       type   (bc) = 'ACTUATOR'
       normDir(bc) =  0
       imin   (bc) =  43
       imax   (bc) =  471
       jmin   (bc) =  234
       jmax   (bc) =  249
       kmin   (bc) =  1
       kmax   (bc) = -1
    end if

    nbc = bc

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
  end subroutine counterflowDiffusionBC


  ! ================== !
  ! Initial conditions !
  ! ================== !
  subroutine counterflowDiffusionInitialCondition(state, grid, generateTargetState)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use State_mod

    ! <<< Internal modules >>>
    use InputHelper, only : getOption
    use ErrorHandler, only : gracefulExit

    implicit none

    ! <<< Arguments >>>
    type(t_State) :: state
    type(t_Grid) :: grid
    logical, intent(in), optional :: generateTargetState

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    logical :: generateTargetState_
    integer :: i, nSpecies, H2, O2, nDimensions, ierror
    real(SCALAR_KIND) :: ratioOfSpecificHeats, density, temperature, fuel, oxidizer,         &
         u, v, Z, T0, flameTemperature, heatRelease, Yf0, Yo0, jetVelocity, a, eta, x, y, Re

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

    ! Species indeces.
    H2 = nDimensions+2+1
    O2 = H2 + 1

    ! Rad in mixture properties.
    call getRequiredOption("initial_fuel_mass_fraction", Yf0)
    call getRequiredOption("initial_oxidizer_mass_fraction", Yo0)
    call getRequiredOption("initial_jet_velocity", jetVelocity)
    call getRequiredOption("Reynolds_number", Re)
    ratioOfSpecificHeats = getOption("ratio_of_specific_heats", 1.4_wp)
    call getRequiredOption("heat_release", heatRelease)

    ! Temperatures.
    T0 =  1.0_wp / (ratioOfSpecificHeats - 1.0_wp)
    flameTemperature = T0 / (1.0_wp - heatRelease)

    ! Stretching parameter
    a = 2.0_wp * jetVelocity /                                                               &
         ( maxval(grid%coordinates(:,2)) - minval(grid%coordinates(:,2)) )

    do i = 1, grid%nGridPoints

       ! Get the local coordinates.
       x = grid%coordinates(i,1)
       y = grid%coordinates(i,2)

       ! Similarity transformation.
       eta = y * sqrt(Re * a)

       ! Velocities.
       u = 0.5_wp * a * x
       v = - a * y

       ! Temperature.
       temperature = T0 + 0.5_wp * ( (1.0_wp + erf(eta / sqrt(2.0_wp) + 1.0_wp)) -           &
            (1.0_wp + erf(eta / sqrt(2.0_wp) - 1.0_wp)) ) * (flameTemperature - T0)

       ! Density.
       density = T0 / temperature

       ! Mixture fraction
       Z = 0.5_wp * ( 1.0_wp + erf(eta / sqrt(2.0_wp)) )

       ! Components.
       fuel = YF0 * Z
       oxidizer = YO0 * (1.0_wp-Z)

       ! State variables
       state%conservedVariables(i,1) = density
       state%conservedVariables(i,2:nDimensions+1) = 0.0_wp
       state%conservedVariables(i,2) = density * u
       state%conservedVariables(i,3) = density * v
       state%conservedVariables(i,nDimensions+2) = density * temperature /                   &
            ratioOfSpecificHeats + 0.5_wp * density * (u**2 + v**2)
       if (nSpecies.gt.0) state%conservedVariables(i,H2) = density * fuel
       if (nSpecies.gt.1) state%conservedVariables(i,O2) = density * oxidizer

       print *, 'Pressure:',0.4_wp*(state%conservedVariables(i,nDimensions+2)-0.5_wp*density*(u**2+v**2))

       ! Target solution
       if (generateTargetState_) state%targetState(i,:) = state%conservedVariables(i,:)

    end do

    return
  end subroutine counterflowDiffusionInitialCondition

end program counterflow_diffusion
