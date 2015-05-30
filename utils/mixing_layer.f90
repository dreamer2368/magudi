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

  !> Generates the grid, BC, initial condition, and target state for a spatial mixing layer.

  implicit none

  integer :: i, numProcs, ierror
  integer :: i1, i2, j1, j2
  character(len = STRING_LENGTH) :: inputname,filename
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
  print *, '!    for a spatially-evolving mixing layer          !'
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
  call mixingLayerGrid(i1,i2,j1,j2)

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
  call mixingLayerBC(i1,i2,j1,j2)

  ! Generate the initial condition and target state.
  do i = 1, size(region%grids)
     call mixingLayerInitialCondition(region%states(i), region%grids(i),                     &
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
    integer :: i, j, k, n
    integer :: nx, ny, nz, nx_, ny_, nz_
    real(wp) :: xmini, xmaxi, ymini, ymaxi
    real(wp) :: xmino, xmaxo, ymino, ymaxo
    real(wp) :: g0, b, c, sig, dy_min, dy_max, y1, y2
    real(wp), allocatable, dimension(:) :: s, g
    logical :: stretch_y

    ! Read in grid size and dimensions.
    nx = getOption("nx", 1025)
    ny = getOption("ny", 513)
    nz = 1
    xmini = getOption("xmin_interior", 0.0_wp)
    xmaxi = getOption("xmax_interior", 120.0_wp)
    ymini = getOption("ymin_interior", -30.0_wp)
    ymaxi = getOption("ymax_interior", 30.0_wp)
    xmino = getOption("xmin_outer", -40.0_wp)
    xmaxo = getOption("xmax_outer", 160.0_wp)
    ymino = getOption("ymin_outer", -50.0_wp)
    ymaxo = getOption("ymax_outer", 50.0_wp)

    ! Allocate the global grid size and assign values.
    ! `globalGridSizes(i,j)` is the number of grid points on grid `j` along dimension `i`.
    allocate(globalGridSizes(2,1))
    globalGridSizes(1,1) = nx
    globalGridSizes(2,1) = ny

    ! Setup the region.
    call region%setup(MPI_COMM_WORLD, globalGridSizes)

    ! Local mesh
    nx_ = region%grids(1)%localSize(1)
    ny_ = region%grids(1)%localSize(2)
    nz_ = region%grids(1)%localSize(3)

    ! Should we stretch the mesh?
    stretch_y = getOption('stretch_y',.false.)

    ! Generate the grid.
    do k = 1, nz_
       do j = 1, ny_
          do i = 1, nx_
             ! Create X
             region%grids(1)%coordinates(i+nx_*(j-1+ny_*(k-1)),1)                            &
                  = (xmaxo-xmino)*real(i-1,wp)/real(nx-1) + xmino

             ! Create Y
             if (.not.stretch_y) then
                region%grids(1)%coordinates(i+nx_*(j-1+ny_*(k-1)),2)                         &
                     = (ymaxo-ymino)*real(j-1,wp)/real(ny-1) + ymino
             end if

             ! Create Z
             !region%grids(1)%coordinates(i+nx_*(j-1+ny_*(k-1)),3)                            &
             !     = 0.0_wp
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
                ! Create y.
                region%grids(1)%coordinates(i+nx_*(j-1+ny_*(k-1)),2)                         &
                     = 0.5_wp*(ymaxo-ymino)*g(j)/g0

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
    print *, 'Extents:',i1,i2,j1,j2
    print*

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
  subroutine mixingLayerBC(i1,i2,j1,j2)
    implicit none

    integer, intent(in) :: i1, i2, j1, j2
    integer :: i, bc, nbc, iunit
    integer, allocatable, dimension(:) :: grid,normDir,imin,imax,jmin,jmax,kmin,kmax
    character(len = 22), allocatable, dimension(:) :: name,type

    ! Number of BC
    nbc = 12

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
    imax   (bc) =  i1
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
    imin   (bc) =  i2
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
    jmax   (bc) =  j1
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
    jmin   (bc) =  j2
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! BC 9
    bc = 9
    name   (bc) = 'excitationSupport'
    type   (bc) = 'SOLENOIDAL_EXCITATION'
    normDir(bc) =  0
    imin   (bc) =  1
    imax   (bc) =  i1
    jmin   (bc) =  1
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! BC 10
    bc = 10
    name   (bc) = 'localizedIgnition'
    type   (bc) = 'GAUSSIAN_IGNITION'
    normDir(bc) =  0
    imin   (bc) =  i1
    imax   (bc) =  i2
    jmin   (bc) =  j1
    jmax   (bc) =  j2
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! BC 11
    bc = 11
    name   (bc) = 'targetRegion'
    type   (bc) = 'COST_TARGET'
    normDir(bc) =  0
    imin   (bc) =  i1
    imax   (bc) =  i2
    jmin   (bc) =  j1
    jmax   (bc) =  j2
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! BC 12
    bc = 12
    name   (bc) = 'controlRegion'
    type   (bc) = 'ACTUATOR'
    normDir(bc) =  0
    imin   (bc) =  i1
    imax   (bc) =  i2
    jmin   (bc) =  j1
    jmax   (bc) =  j2
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
         density, velocity, temperature, Z, fuel, oxidizer, YF0, YO0

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

    ! Species parameters.
    H2 = nDimensions+2+1
    O2 = H2 + 1
    YF0 = getOption("initial_fuel_mass_fraction", 0.0_wp)
    YO0 = getOption("initial_oxidizer_mass_fraction", 0.0_wp)

    ! Gamma
    ratioOfSpecificHeats = getOption("ratio_of_specific_heats", 1.4_wp)

    do i = 1, grid%nGridPoints

       ! Mixture fraction
       Z = 0.5_wp*(1.0_wp-erf(grid%coordinates(i,2)))

       ! Components
       fuel = YF0*Z
       oxidizer = YO0*(1.0_wp-Z)

       ! Density
       if (nspecies.gt.0) then
          density = 1.0_wp
       else
          density = 1.0_wp
       end if

       ! Velocity
       velocity = lowerVelocity +                                                          &
            0.5_wp * (upperVelocity - lowerVelocity) * (1.0_wp + tanh(2.0_wp *             &
            grid%coordinates(i,2)))

       ! Temperature
       temperature =  1.0_wp / (ratioOfSpecificHeats - 1.0_wp)

       ! State variables
       state%conservedVariables(i,1) = density
       state%conservedVariables(i,2) = state%conservedVariables(i,1) * velocity
       state%conservedVariables(i,3:nDimensions+1) = 0.0_wp
       state%conservedVariables(i,nDimensions+2) =                                           &
            state%conservedVariables(i,1) * temperature / ratioOfSpecificHeats +             &
            0.5_wp * state%conservedVariables(i,1) * velocity ** 2
       if (nSpecies.gt.0) state%conservedVariables(i,H2) = fuel*state%conservedVariables(i,1)
       if (nSpecies.gt.1) state%conservedVariables(i,O2) = oxidizer*state%conservedVariables(i,1)

       ! Target solution
       if (generateTargetState_) state%targetState(i,:) = state%conservedVariables(i,:)

    end do

    return
  end subroutine mixingLayerInitialCondition

end program mixing_layer
