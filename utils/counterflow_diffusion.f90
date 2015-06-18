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
    integer :: i, j, k, gridIndex
    integer :: nx, ny, nz, nx_, ny_, nz_
    real(wp) :: xmini, xmaxi, ymini, ymaxi
    real(wp) :: xmino, xmaxo, ymino, ymaxo
    real(wp) :: zmin, zmax
    real(wp) :: dx, dy, dz

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
             region%grids(1)%coordinates(gridIndex,2) =                                      &
                  (ymaxo - ymino) * real(j-1,wp) / real(ny-1) + ymino

             ! Create Z
             if (nz.ne.1) region%grids(1)%coordinates(gridIndex,3) =                         &
                  (zmax - zmin - dz) * real(k-1,wp) / real(nz-1) + zmin

          end do
       end do
    end do

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
       imin   (bc) =  29
       imax   (bc) =  173
       jmin   (bc) =  86
       jmax   (bc) =  116
       kmin   (bc) =  1
       kmax   (bc) = -1
    end if

    entry = getOption("control_mollifier_file", "")
    if (len_trim(entry) > 0) then
       bc = bc + 1
       name   (bc) = 'controlRegion'
       type   (bc) = 'ACTUATOR'
       normDir(bc) =  0
       imin   (bc) =  29
       imax   (bc) =  173
       jmin   (bc) =  143
       jmax   (bc) =  173
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
         u, v, Z, Z0, T0, Yf0, Yo0, jetVelocity, a, eta, x, y, Re

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

    ! Mixture properties
    Z0 = getOption("initial_mixture_fraction", 1.0_wp)
    call getRequiredOption("initial_fuel_mass_fraction", Yf0)
    call getRequiredOption("initial_oxidizer_mass_fraction", Yo0)
    call getRequiredOption("initial_jet_velocity", jetVelocity)
    call getRequiredOption("Reynolds_number", Re)

    ! Gamma
    ratioOfSpecificHeats = getOption("ratio_of_specific_heats", 1.4_wp)

    ! Temperature
    T0 =  1.0_wp / (ratioOfSpecificHeats - 1.0_wp)
    temperature = getOption("initial_temperature", T0)

    ! Density
    density = 1.0_wp

    ! Stretching parameter
    a = 2.0_wp * jetVelocity /                                                               &
         ( maxval(grid%coordinates(:,2)) - minval(grid%coordinates(:,2)) )

    do i = 1, grid%nGridPoints

       ! Get the local coordinates
       x = grid%coordinates(i,1)
       y = grid%coordinates(i,2)

       ! Velocities
       u = 0.5_wp * a * x
       v = - a * y

       ! Similarity transformation
       eta = y * sqrt(Re * a)

       ! Mixture fraction
       Z = 0.5_wp * ( 1.0_wp + erf(eta / sqrt(2.0_wp)) )

       ! Components
       fuel = YF0 * Z * Z0
       oxidizer = YO0 * (1.0_wp-Z) * (1.0_wp - Z0)

       ! State variables
       state%conservedVariables(i,1) = density
       state%conservedVariables(i,2:nDimensions+1) = 0.0_wp
       state%conservedVariables(i,2) = density * u
       state%conservedVariables(i,3) = density * v
       state%conservedVariables(i,nDimensions+2) = density * temperature /                   &
            ratioOfSpecificHeats
       state%conservedVariables(i,nDimensions+2) = density * temperature /                   &
            ratioOfSpecificHeats + 0.5_wp * density * (u**2 + v**2)
       if (nSpecies.gt.0) state%conservedVariables(i,H2) = density * fuel
       if (nSpecies.gt.1) state%conservedVariables(i,O2) = density * oxidizer

       ! Target solution
       if (generateTargetState_) state%targetState(i,:) = state%conservedVariables(i,:)

    end do

    return
  end subroutine counterflowDiffusionInitialCondition

end program counterflow_diffusion
