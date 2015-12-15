#include "config.h"

program jet_crossflow

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env

 ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Grid_enum
  use State_enum
  use SolverOptions_enum

  ! <<< Internal modules >>>
  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush, gracefulExit
  use PLOT3DHelper, only : plot3dDetectFormat

  implicit none

  ! <<< Global Parameters >>>
  type(t_Region) :: region
  type(t_Grid) :: grid
  type(t_State) :: state
  type(t_SolverOptions) :: solverOptions
  type(t_SimulationFlags) :: simulationFlags
  integer :: i, numProcs, ierror, nx, ny, nz
  integer :: iJet1, iJet2, kJet1, kJet2
  integer :: imin_sponge, imax_sponge, jmax_sponge
  real(KIND=8) :: xmini, xmaxi, ymaxi
  real(KIND=8) :: xmino, xmaxo, ymino, ymaxo
  real(KIND=8) :: zmin, zmax, jetDiameter, xJet
  character(len = STRING_LENGTH) :: inputname, filename
  integer, allocatable :: globalGridSizes(:,:)
  logical :: includeSandpaper

  print *
  print *, '! ============================================= !'
  print *, '!                                               !'
  print *, '!    JET IN CROSSFLOW GENERATOR                 !'
  print *, '!    Creates:                                   !'
  print *, '!               - Grid with sandpaper           !'
  print *, '!               - Target/initial solution       !'
  print *, '!               - Boundary conditions           !'
  print *, '!    for a jet in cross flow of a               !'
  print *, '!    spatially-evolving boundary layer          !'
  print *, '!                                               !'
  print *, '! ============================================= !'
  print *

  ! Initialize MPI.
  call MPI_Init(ierror)

  ! Exit if executed in parallel
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)
  if (numProcs > 1) then
     print *, 'jet_crossflow.f90 only implemented in serial for now...'
     stop
  end if

  ! Parse options from the input file.
  inputname = "magudi.inp"
  call parseInputFile(inputname)

  ! Generate the grid.
  call jetCrossFlowGrid(imin_sponge, imax_sponge, jmax_sponge)

  ! Save the grid.
  region%grids(1) = grid
  call getRequiredOption("grid_file", filename)
  call region%saveData(QOI_GRID, filename)

  ! Compute normalized metrics, norm matrix, and Jacobian.
  do i = 1, size(region%grids)
     call region%grids(i)%update()
  end do
  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  ! Write out some useful information.
  call region%reportGridDiagnostics()

  ! Generate the boundary conditions.
  call jetCrossFlowBC(imin_sponge, imax_sponge, jmax_sponge)

  ! Generate the initial condition and target state.
  call jetCrossFlowInitialCondition

  ! Save initial condition.
  region%states(1) = state
  call getRequiredOption("grid_file", filename)
  i = len_trim(filename)
  if (filename(i-3:i) == ".xyz") then
     filename = filename(:i-4) // ".ic.q"
  else
     filename = PROJECT_NAME // ".ic.q"
  end if
  call region%saveData(QOI_FORWARD_STATE, filename)

  ! Save target state.
  if (simulationFlags%useTargetState) then
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
  subroutine jetCrossFlowGrid(i1, i2, j2)

    ! Jet in crossflow schematic and the corresponding indeces 
    !
    !
    ! y = tripHeight >     _~~~~~~~~~~~~~~~~~~~~~~~~~~_
    !                     |                            |
    ! y = 0  > ___________|                            |_________________|___:___|__________
    !
    !          ^          ^                                              ^       ^         ^
    !         x=0        x=tripLocation                               i=iJet1   i=iJet2  i=nx

    ! <<< External modules >>>
    use MPI

    ! <<< Arguments >>>
    integer, intent(out) :: i1, i2, j2

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, n, nGrit, iTrip1, iTrip2
    integer, allocatable :: jetBox(:,:)
    real(wp) :: x, y, z, dx, dz, ytilde, r, delta
    real(wp) :: tripLocation, tripWidth, tripHeight
    real(wp) :: gritSize, rnd, gauss, amp, sig, x0, y0, z0, z12, alpha, theta
    logical :: includeGrit, stretchY, conformToJet

    ! Read in grid size and dimensions.
    call getRequiredOption("nx", nx)
    call getRequiredOption("ny", ny)
    call getRequiredOption("nz", nz)
    call getRequiredOption("xmin_interior", xmini)
    call getRequiredOption("xmax_interior", xmaxi)
    call getRequiredOption("ymax_interior", ymaxi)
    call getRequiredOption("xmin_outer", xmino)
    call getRequiredOption("xmax_outer", xmaxo)
    call getRequiredOption("ymin_outer", ymino)
    call getRequiredOption("ymax_outer", ymaxo)
    zmin = getOption("zmin", 0.0_wp)
    zmax = getOption("zmax", 0.0_wp)

    ! Grid spacing
    dx = (xmaxo - xmino) / real(nx, wp)
    dz = (zmax - zmin) / real(nz, wp)

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

    ! Setup the region and simplify.
    call region%setup(MPI_COMM_WORLD, globalGridSizes)
    grid = region%grids(1)
    simulationFlags = region%simulationFlags
    solverOptions = region%solverOptions
    state = region%states(1)

    ! Should we stretch the mesh?
    stretchY = getOption('stretch_y', .false.)
    if (stretchY) r = 2.0_wp

    ! Generate the grid.
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             ! Create X
             grid%coordinates(i+nx*(j-1+ny*(k-1)),1) =                                       &
                  (xmaxo - xmino) * real(i-1, wp) / real(nx-1, wp) + xmino

             ! Create Y
             if (stretchY) then
                ytilde = real(ny-j, wp)/real(ny-1, wp)
                grid%coordinates(i+nx*(j-1+ny*(k-1)),2) =                                    &
                     (ymaxo - ymino) * (1.0_wp - tanh(r * ytilde) / tanh(r)) + ymino
             else
                grid%coordinates(i+nx*(j-1+ny*(k-1)),2) =                                    &
                     (ymaxo - ymino) * real(j-1,wp) / real(ny-1, wp) + ymino
             end if

             ! Create Z
             if (nz.ne.1) grid%coordinates(i+nx*(j-1+ny*(k-1)),3) =                          &
                  (zmax - zmin - dz) * real(k-1, wp) / real(nz-1, wp) + zmin
          end do
       end do
    end do

    ! Generate the sandpaper.
    includeSandpaper = getOption("include_sandpaper", .false.)
    if (includeSandpaper) Then
       call getRequiredOption("trip_location", tripLocation)
       call getRequiredOption("trip_width", tripWidth)
       call getRequiredOption("trip_height", tripHeight)
       includeGrit = getOption("include_grit", .false.)
       if (includeGrit) then
          call getRequiredOption("grit_diameter", gritSize)
          call getRequiredOption("number_of_particles", nGrit)
          if (tripHeight > gritSize) tripHeight = tripHeight - gritSize
       end if

       ! Find the sandpaper extents.
       j = 1; k = 1
       do i = 1, nx - 1
          if (grid%coordinates(i+nx*(j-1+ny*(k-1)),1) < tripLocation)                        &
               iTrip1 = i + 1
       end do
       do i = nx, 2, -1
          if (grid%coordinates(i+nx*(j-1+ny*(k-1)),1) > tripLocation +                       &
               tripWidth) iTrip2 = i - 1
       end do
       print *
       print *, 'Sandpaper extents: [', iTrip1, ',', iTrip2, ']'


       ! Deform the mesh to the sandpaper height.
       do k = 1, nz
          do i = 1, nx
             ! Create a smooth step.
             j = 1
             sig = 20.0_wp
             x = grid%coordinates(i+nx*(j-1+ny*(k-1)),1)
             grid%coordinates(i+nx*(j-1+ny*(k-1)),2) = 0.5_wp * tripHeight *                 &
                  (tanh(sig*(x - tripLocation)) - tanh(sig*(x - tripLocation-tripWidth)))

             ! Shift grid points above (smoothly).
             do j = 2, ny
                y = grid%coordinates(i+nx*(j-1+ny*(k-1)),2)
                delta = grid%coordinates(1+nx*(j-1+ny*(1-1)),2) -                            &
                     grid%coordinates(1+nx*(j-2+ny*(1-1)),2)
                ytilde = grid%coordinates(i+nx*(j-2+ny*(k-1)),2) +                           &
                     delta
                alpha = tanh(0.1_wp * (j - 2) / (ny - 2))
                grid%coordinates(i+nx*(j-1+ny*(k-1)),2) =                                    &
                     ytilde * (1.0_wp - alpha) + y * alpha
             end do
          end do
       end do
  
       ! Embed the particles.
       If (includeGrit) then

          ! Standard deviation
          sig = 1.5_wp * gritSize

          ! Loop through number of particles.
          do n = 1, nGrit

             ! Compute amplitude.
             call random_number(rnd)
             amp = gritSize * (1.0_wp + 0.1_wp * (rnd - 0.5_wp))
             call random_number(rnd)
             if (rnd < 0.5_wp) amp = -amp

             ! Get a random location.
             call random_number(rnd)
             x0 = tripLocation + 5.0_wp * sig + (tripWidth - 10.0_wp * sig) * rnd
             call random_number(rnd)
             z0 = (zmax - zmin - dz) * rnd + zmin
             if (nz == 1) z0 = 0.0_wp

             ! Modify the grid.
             do k = 1, nz
                do i = iTrip1, iTrip2

                   ! Get the coordinates.
                   j = 1
                   x = grid%coordinates(i+nx*(j-1+ny*(k-1)),1)
                   y0 = grid%coordinates(i+nx*(j-1+ny*(k-1)),2)
                   z = grid%coordinates(i+nx*(j-1+ny*(k-1)),3)
                      
                   ! Represent sandpaper particles as Gaussian.
                   gauss = amp * exp(-((x - x0)**2 / (2.0_wp * sig**2) +                     &
                        (z - z0)**2 / (2.0_wp * sig**2)))

                   ! Account for periodicity in z.
                   z12 = (zmax - zmin) - abs(z - z0)
                   If (nz > 1) gauss = gauss + amp *                                         &
                        exp(-((x-x0)**2 / (2.0_wp * sig**2) + z12**2 / (2.0_wp * sig**2)))

                   ! Update the vertical coordinate.
                   grid%coordinates(i+nx*(j-1+ny*(k-1)),2) =                                 &
                        grid%coordinates(i+nx*(j-1+ny*(k-1)),2) + gauss
                   
                   ! Shift grid points above (smoothly).
                   do j = 2, ny
                      y = grid%coordinates(i+nx*(j-1+ny*(k-1)),2)
                      delta = grid%coordinates(1+nx*(j-1+ny*(1-1)),2) -                      &
                           grid%coordinates(1+nx*(j-2+ny*(1-1)),2)
                      ytilde = grid%coordinates(i+nx*(j-2+ny*(k-1)),2) +                     &
                           delta
                      alpha = tanh(2.0_wp * (j - 2) / (ny - 2))
                      grid%coordinates(i+nx*(j-1+ny*(k-1)),2) =                              &
                           ytilde * (1.0_wp - alpha) + y * alpha
                   end do

                end do
             end do
          end do
       end if !... includeGrit

    end if !...  includeSandpaper

    ! Find the original jet extents.
    call getRequiredOption("jet_diameter", jetDiameter)
    call getRequiredOption("jet_position", xJet)
    iJet1 = nx; iJet2 = 1
    kJet1 = nz; kJet2 = 1
    j = 1
    do k = 1, nz
       do i = 1, nx
          x = grid%coordinates(i+nx*(j-1+ny*(k-1)),1)
          z = grid%coordinates(i+nx*(j-1+ny*(k-1)),3)
          r = sqrt((x-xJet)**2 + (z - 0.5_wp * (zmax + zmin))**2)
          if (r <= 0.5_wp * jetDiameter) then
             iJet1 = min(iJet1, i)
             iJet2 = max(iJet2, i)
             kJet1 = min(kJet1, k)
             kJet2 = max(kJet2, k)
          end if
       end do
    end do

    ! Conform the mesh to the jet perimeter.
    conformToJet = getOption("conform_grid_to_jet", .false.)
    if (conformToJet) then
       j = 1
       r = 0.495_wp * jetDiameter
       allocate(jetBox(iJet1-1:iJet2+1, kJet1-1:kJet2+1))
       jetBox = 0

       ! Find a point at `i` by shifting in `k`.
       do i = iJet1, iJet2
          k = 1
          x = grid%coordinates(i+nx*(j-1+ny*(k-1)),1)

          ! First quadrant.
          delta = x - xJet
          theta = asin(delta / r)
          delta = r * cos(theta)
          z0 = 0.5_wp * (zmax + zmin) + delta
          delta = zmax - zmin
          do k = kJet1 - 1, kJet2 + 1
             ! Find the closest point to the jet perimeter.
             z = grid%coordinates(i+nx*(j-1+ny*(k-1)),3)
             if (abs(z - z0) < delta) then
                delta = abs(z - z0)
                n = k
             end if
          end do
          ! Shift the corresponding grid point.
          k = n
          if (jetBox(i, k) == 0) grid%coordinates(i+nx*(j-1+ny*(k-1)),3) = z0
          jetBox(i, k) = 1 ! Tag it.

          ! Second quadrant.
          delta = x - xJet
          theta = acos(delta / r)
          delta = r * sin(theta)
          z0 = 0.5_wp * (zmax + zmin) - delta
          delta = zmax - zmin
          do k = kJet1 - 1, kJet2 + 1
             ! Find the closest point to the jet perimeter.
             z = grid%coordinates(i+nx*(j-1+ny*(k-1)),3)
             if (abs(z - z0) < delta) then
                delta = abs(z - z0)
                n = k
             end if
          end do
          ! Shift the corresponding grid point.
          k = n
          if (jetBox(i, k) == 0) grid%coordinates(i+nx*(j-1+ny*(k-1)),3) = z0
          jetBox(i, k) = 1 ! Tag it.
       end do

       ! Find a point at `k` by shifting in `i`.
       do k = kJet1, kJet2
          i = 1
          z = grid%coordinates(i+nx*(j-1+ny*(k-1)),3)

          ! Third quadrant.
          delta = z - 0.5_wp * (zmax + zmin)
          theta = asin(delta / r)
          delta = r * cos(theta)
          x0 = xJet - delta
          delta = xmaxo - xmino
          do i = iJet1 - 1, iJet2 + 1
             ! Find the closest point to the jet perimeter.
             x = grid%coordinates(i+nx*(j-1+ny*(k-1)),1)
             if (abs(x - x0) < delta) then
                delta = abs(x - x0)
                n = i
             end if
          end do
          ! Shift the corresponding grid point.
          i = n
          if (jetBox(i, k) == 0) grid%coordinates(i+nx*(j-1+ny*(k-1)),1) = x0
          jetBox(i, k) = 1 ! Tag it.

          ! Fourth quadrant.
          delta = z - 0.5_wp * (zmax + zmin)
          theta = acos(delta / r)
          delta = r * sin(theta)
          x0 = xJet + delta
          delta = xmaxo - xmino
          do i = iJet1 - 1, iJet2 + 1
             ! Find the closest point to the jet perimeter.
             x = grid%coordinates(i+nx*(j-1+ny*(k-1)),1)
             if (abs(x - x0) < delta) then
                delta = abs(x - x0)
                n = i
             end if
          end do
          ! Shift the corresponding grid point.
          i = n
          if (jetBox(i, k) == 0) grid%coordinates(i+nx*(j-1+ny*(k-1)),1) = x0
          jetBox(i, k) = 1 ! Tag it.
       end do

       ! Smooth surrounding grid points.
       j = 1
       do k = kJet1 - 1, kJet2 + 1
          do i = iJet1 - 1, iJet2 + 1
             if (jetBox(i, k) == 1) then
                do n = 1, 41
                   ! Stretching parameter.
                   alpha = tanh(1.0_wp * (n - 1) / 40)

                   ! Shift left.
                   x = grid%coordinates(i-n+nx*(j-1+ny*(k-1)),1)
                   delta = grid%coordinates(i-n+nx*(j+ny*(k-1)),1) -                         &
                        grid%coordinates(i-n+1+nx*(j+ny*(k-1)),1)
                   ytilde = grid%coordinates(i-n+1+nx*(j-1+ny*(k-1)),1) +                    &
                        delta
                   grid%coordinates(i-n+nx*(j-1+ny*(k-1)),1) =                               &
                        ytilde * (1.0_wp - alpha) + x * alpha

                   ! Shift right.
                   x = grid%coordinates(i+n+nx*(j-1+ny*(k-1)),1)
                   delta = grid%coordinates(i+n+nx*(j+ny*(k-1)),1) -                         &
                        grid%coordinates(i+n-1+nx*(j+ny*(k-1)),1)
                   ytilde = grid%coordinates(i+n-1+nx*(j-1+ny*(k-1)),1) +                    &
                        delta
                   grid%coordinates(i+n+nx*(j-1+ny*(k-1)),1) =                               &
                        ytilde * (1.0_wp - alpha) + x * alpha

                   ! Shift up.
                   x = grid%coordinates(i+nx*(j-1+ny*(k-1-n)),3)
                   delta = grid%coordinates(i+nx*(j+ny*(k-1-n)),3) -                         &
                        grid%coordinates(i+nx*(j+ny*(k-n)),3)
                   ytilde = grid%coordinates(i+nx*(j-1+ny*(k-n)),3) +                        &
                        delta
                   grid%coordinates(i+nx*(j-1+ny*(k-1-n)),3) =                               &
                        ytilde * (1.0_wp - alpha) + x * alpha

                   ! Shift down.
                   x = grid%coordinates(i+nx*(j-1+ny*(k-1+n)),3)
                   delta = grid%coordinates(i+nx*(j+ny*(k-1+n)),3) -                         &
                        grid%coordinates(i+nx*(j+ny*(k-2+n)),3)
                   ytilde = grid%coordinates(i+nx*(j-1+ny*(k-2+n)),3) +                      &
                        delta
                   grid%coordinates(i+nx*(j-1+ny*(k-1+n)),3) =                               &
                        ytilde * (1.0_wp - alpha) + x * alpha

                end do
             end if
          end do
       end do
       deallocate(jetBox)

       ! Find the new jet extents.
       iJet1 = nx; iJet2 = 1
       kJet1 = nz; kJet2 = 1
       j = 1
       do k = 1, nz
          do i = 1, nx
             x = grid%coordinates(i+nx*(j-1+ny*(k-1)),1)
             z = grid%coordinates(i+nx*(j-1+ny*(k-1)),3)
             r = sqrt((x-xJet)**2 + (z - 0.5_wp * (zmax + zmin))**2)
             if (r <= 1.0_wp * jetDiameter) then
                iJet1 = min(iJet1, i)
                iJet2 = max(iJet2, i)
                kJet1 = min(kJet1, k)
                kJet2 = max(kJet2, k)
             end if
          end do
       end do
    end if !... conform to jet.

    print *
    print *, 'Jet extents: [', iJet1, ',', iJet2, '] x [',kJet1, ',', kJet2,']'

    ! Find extents of outer region.
    j=1; k=1;
    i1=1
    do i = 1, nx
       if (grid%coordinates(i+nx*(j-1+ny*(k-1)),1) <= xmini) i1=i+1
    end do
    i2=nx
    do i = nx,1,-1
       if (grid%coordinates(i+nx*(j-1+ny*(k-1)),1) >= xmaxi) i2=i
    end do
    i=1; k=1;
    j2=ny
    do j = ny,1,-1
       if (grid%coordinates(i+nx*(j-1+ny*(k-1)),2) >= ymaxi) j2=j
    end do

    print *
    print *, 'Interior domain extents: [',i1,',',i2,'] x [',1,',',j2,']'
    print*

    ! Clean up.
    deallocate(globalGridSizes)

    return
  end subroutine jetCrossFlowGrid


  ! ================== !
  ! Initial conditions !
  ! ================== !
  subroutine jetCrossFlowInitialCondition

    ! <<< External modules >>>
    use MPI

    ! <<< Internal modules >>>
    use InputHelper, only : getOption
    use ErrorHandler, only : gracefulExit
    use ThermoChemistry, only : getMolecularWeight

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, l, kk, nDimensions, nSpecies, H2, O2, N2, ierror
    real(wp) :: ratioOfSpecificHeats, crossflowVelocity, jetVelocity, x, y, z, y0, r, sig,   &
         density, velocity(3), temperature, pressure, fuel, oxidizer, inert, YF0, Yo0, yDecay
    real(wp), dimension(:), allocatable :: Wi
    real(SCALAR_KIND), parameter :: spreadingRate = 0.094_wp

    ! Solution to Blasius boundary layer
    real(wp) :: blasius0, blasius1, delta, eta, xx, x0, Re_c
    real(wp) :: f2l, f2r, f0l, f0r
    real(wp), dimension(0:9) :: by0 = (/                                                     &
         0.0_8, 0.165571818583440_8, 0.650024518764203_8, 1.39680822972500_8,                &
         2.30574664618049_8, 3.28327391871370_8, 4.27962110517696_8,                         &
         5.27923901129384_8, 6.27921363832835_8, 7.27921257797747_8 /)
    real(wp), dimension(0:9) :: by1 = (/                                                     &
         0.0_8, 0.329780063306651_8, 0.629765721178679_8, 0.84604458266019_8,                &
         0.95551827831671_8, 0.99154183259084_8, 0.99897290050990_8,                         &
         0.9999216098795_8, 0.99999627301467_8, 0.99999989265063_8 /)
    real(wp), dimension(0:9) :: by2 = (/                                                     &
         0.332057384255589_8, 0.323007152241930_8, 0.266751564401387_8, 0.161360240845588_8, &
         0.06423404047594_8, 0.01590689966410_8, 0.00240199722109_8,                         &
         0.00022016340923_8, 0.00001224984692_8, 0.00000041090325_8 /)

    call MPI_Cartdim_get(grid%comm, nDimensions, ierror)

    ! Species
    nSpecies = solverOptions%nSpecies

    ! Only implemented for 2 species
    if (nspecies.gt.2) then
       print *, 'WARNING, max of 2 species (for now)'
       stop
    end if

    ! Get species indices.
    if (allocated(solverOptions%speciesName)) then
       do k = 1, nSpecies + 1
          select case (trim(solverOptions%speciesName(k)))
          case ('H2', 'HYDROGEN')
             H2 = k
          case ('O2', 'OXYGEN')
             O2 = k
          case ('N2', 'NITROGEN')
             N2 = k
          case default
             print *, "Unknown species: ", trim(solverOptions%speciesName(k)), "!"
             stop
          end select
       end do
    else
       H2 = 1
       O2 = 2
       N2 = nSpecies + 1
    end if

    ! Get molecular weights.
    if (solverOptions%equationOfState == IDEAL_GAS_MIXTURE) then
       allocate(Wi(nSpecies+1))
       Wi = solverOptions%molecularWeightInverse
    end if

    ! Read in the parameters needed for generating the crossflow & jet velocities.
    call getRequiredOption("crossflow_velocity", crossflowVelocity)
    call getRequiredOption("Reynolds_number", Re_c)
    call getRequiredOption("Blasius_virtual_origin", x0)
    call getRequiredOption("jet_velocity", jetVelocity)

    ! Mixture properties.
    call getRequiredOption("initial_fuel_mass_fraction", Yf0)
    call getRequiredOption("initial_oxidizer_mass_fraction", Yo0)

    ! Get the ratio of specific heats.
    ratioOfSpecificHeats = solverOptions%ratioOfSpecificHeats

    do k = 1, nz
       do j = 1, ny
          do i = 1, nx

             ! Get the local coordinates.
             l = i+nx*(j-1+ny*(k-1))
             x = grid%coordinates(l,1)
             y = grid%coordinates(l,2)
             z = grid%coordinates(l,3)
             y0 = grid%coordinates(i+nx*(1-1+ny*(k-1)),2)

             ! Initialize the mixture.
             velocity = 0.0_wp
             pressure = 1.0_wp / ratioOfSpecificHeats
             temperature =  1.0_wp / (ratioOfSpecificHeats - 1.0_wp)
             fuel = 0.0_wp
             oxidizer = Yo0
             density = 1.0_wp

             ! Virtual origin of the Blasius profile.
             xx = x + x0

             ! Return the first derivative of the Blasius function.
             delta = sqrt(xx / Re_c / crossflowVelocity)
             eta = (y - y0) / delta
             if (eta <= 0.0_WP) then
                blasius0 = 0.0_wp
                blasius1 = 0.0_wp
             else if (eta >= 9.0_wp) then
                blasius0 = by0(9) + (eta-9.0_WP)
                blasius1 = 1.0_WP
             else
                kk = int(eta)

                f0l = by0(kk)
                f0r = by0(kk + 1)
                f2l = by2(kk)
                f2r = by2(kk + 1)
                blasius0 = &
                     1.0_WP/6.0_wp*f2l*(real(kk+1, wp) - eta)**3 +                           &
                     1.0_WP/6.0_wp*f2r*(eta-real(kk, wp))**3 +                               &
                     (f0l-1.0_WP/6.0_WP*f2l)*(real(kk+1, wp) - eta) +                        &
                     (f0r-1.0_WP/6.0_WP*f2r)*(eta-real(kk, wp))

                f0l = by1(kk)
                f0r = by1(kk+1)
                f2l = -0.5_wp * by0(kk)*by2(kk)
                f2r = -0.5_wp * by0(kk+1)*by2(kk+1)
                blasius1 = &
                     1.0_wp/6.0_wp*f2l*(real(kk+1, wp)-eta)**3 +                             &
                     1.0_wp/6.0_wp*f2r*(eta-real(kk, wp))**3 +                               &
                     (f0l-1.0_wp/6.0_WP*f2l)*(real(kk+1, wp) - eta) +                        &
                     (f0r-1.0_wp/6.0_WP*f2r)*(eta-real(kk, wp))

             end if
    
             ! Update the velocity.
             velocity(1) = crossflowVelocity * blasius1
             velocity(2) = 0.5_WP * sqrt(crossflowVelocity / xx / Re_c) *                    &
                  (eta * blasius1 - blasius0)

             ! Bottom wall.
             r = sqrt((x - xJet)**2 + (z - 0.5_wp * (zmax + zmin))**2)
             if (j == 1 .and. r <= jetDiameter) then
                ! Start with zero velocity at the bottom.
                velocity = 0.0_wp

                ! Jet conditions.
                sig = 10.0_wp
                yDecay = 1.0_wp
                if (y > 0.5_wp * jetDiameter) yDecay = max(1.0_wp - 2.0_wp *                 &
                     (y - 0.5_wp * jetDiameter) / jetDiameter, 0.0_wp)

                ! Fuel stream.
                fuel = Yf0 * yDecay * 0.5_wp * (tanh(sig * (r + 0.5_wp * jetDiameter)) -     &
                     tanh(sig * (r - 0.5_wp * jetDiameter)))
                oxidizer = (1.0_wp - fuel) * Yo0

                ! Poiseuille velocity profile.
                !velocity(2) = max(0.0_wp,                                                 &
                !     2.0_wp * jetVelocity * (1.0_wp - (r / (0.5_wp*jetDiameter))**2) *    &
                !     yDecay)
                !velocity(2) = jetVelocity * tanh(4.0_wp * (1.0_wp - 2.0_wp * r / jetDiameter)) * yDecay

                ! Tanh velocity profile.
                velocity(2) = jetVelocity * yDecay * 0.5_wp * (tanh(sig *                    &
                     (r + 0.5_wp * jetDiameter)) - tanh(sig * (r - 0.5_wp * jetDiameter)))

                ! Self-similar velocity profile.
                !eta = (x - xJet)
                !a = (sqrt(2.0_wp) - 1.0_wp) / (0.5_wp * jetDiameter)**2
                !velocity(2) = jetVelocity * (1.0_wp + a * eta**2) ** (-2.0_wp)
                !velocity(1) = 0.5_wp * jetVelocity * (eta - a * eta**3) / (1.0_wp + a * eta**2)**2

             end if

             ! Correct species mass fractions.
             inert = 1.0_wp - fuel - oxidizer
             if (inert < 0.0_wp) oxidizer = oxidizer + inert

             ! Get density from the equation of state
             select case (solverOptions%equationOfState)
             case(IDEAL_GAS)
                density = ratioOfSpecificHeats * pressure /                                  &
                     (temperature * (ratioOfSpecificHeats - 1.0_wp))
             case (IDEAL_GAS_MIXTURE)
                density = ratioOfSpecificHeats * pressure /                                  &
                     ( temperature * (ratioOfSpecificHeats - 1.0_wp) *                       &
                     (fuel * (Wi(H2) - Wi(N2)) + oxidizer * (Wi(O2) - Wi(N2)) + Wi(N2)) )
             end select

             ! State variables.
             state%conservedVariables(l,1) = density
             state%conservedVariables(l,2:nDimensions+1) =                                   &
                  state%conservedVariables(l,1) * velocity(1:nDimensions)
             state%conservedVariables(l,nDimensions+2) = pressure /                          &
                  (ratioOfSpecificHeats - 1.0_wp) +                                          &
                  0.5_wp * state%conservedVariables(l,1) * sum(velocity ** 2)
             if (nSpecies.gt.0) state%conservedVariables(l, nDimensions+2+H2) = fuel *       &
                  state%conservedVariables(l,1)
             if (nSpecies.gt.1) state%conservedVariables(l, nDimensions+2+O2) = oxidizer *   &
                  state%conservedVariables(l,1)

             ! Target solution.
             if (simulationFlags%useTargetState) state%targetState(l,:) =                    &
                  state%conservedVariables(l,:)

          end do
       end do
    end do

    return
  end subroutine jetCrossFlowInitialCondition


  ! =================== !
  ! Boundary conditions !
  ! =================== !
  subroutine jetCrossFlowBC(imin_sponge, imax_sponge, jmax_sponge)

    implicit none

    ! <<< Arguments >>>
    integer, intent(in), optional :: imin_sponge, imax_sponge, jmax_sponge

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, bc, nbc, iunit
    integer, allocatable, dimension(:) :: normDir, imin, imax, jmin, jmax, kmin, kmax
    character(len = STRING_LENGTH) :: str
    character(len = 22), allocatable, dimension(:) :: name, type

    ! Initialize a large number of boundary conditions.
    nbc = 99

    ! Allocate boundary conditions.
    allocate(name(nbc), type(nbc), normDir(nbc),                                             &
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
    imax   (bc) =  imin_sponge
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
    imin   (bc) =  imax_sponge
    imax   (bc) = -1
    jmin   (bc) =  1
    jmax   (bc) = -1
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
    bc = bc+1
    name   (bc) = 'topSponge'
    type   (bc) = 'SPONGE'
    normDir(bc) = -2
    imin   (bc) =  1
    imax   (bc) = -1
    jmin   (bc) =  jmax_sponge
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! Bottom BCs (no-slip wall)
    i = 0
    bc = bc+1
    i = i + 1
    write(str,'(A10,i1)') 'bottomWall', i
    name   (bc) = trim(str)
    type   (bc) = 'SAT_ISOTHERMAL_WALL'
    normDir(bc) =  2
    imin   (bc) =  1
    imax   (bc) =  iJet1 - 1
    jmin   (bc) =  1
    jmax   (bc) =  1
    kmin   (bc) =  1
    kmax   (bc) = -1

    if (nz > 1) then
       bc = bc+1
       i = i + 1
       write(str,'(A10,i1)') 'bottomWall', i
       name   (bc) = trim(str)
       type   (bc) = 'SAT_ISOTHERMAL_WALL'
       normDir(bc) =  2
       imin   (bc) =  iJet1
       imax   (bc) =  iJet2
       jmin   (bc) =  1
       jmax   (bc) =  1
       kmin   (bc) =  1
       kmax   (bc) =  kJet1 - 1

       bc = bc+1
       i = i + 1
       write(str,'(A10,i1)') 'bottomWall', i
       name   (bc) = trim(str)
       type   (bc) = 'SAT_ISOTHERMAL_WALL'
       normDir(bc) =  2
       imin   (bc) =  iJet1
       imax   (bc) =  iJet2
       jmin   (bc) =  1
       jmax   (bc) =  1
       kmin   (bc) =  kJet2 + 1
       kmax   (bc) = -1
    end if

    bc = bc+1
       i = i + 1
       write(str,'(A10,i1)') 'bottomWall', i
    name   (bc) = trim(str)
    type   (bc) = 'SAT_ISOTHERMAL_WALL'
    normDir(bc) =  2
    imin   (bc) =  iJet2 + 1
    imax   (bc) = -1
    jmin   (bc) =  1
    jmax   (bc) =  1
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! Bottom wall (inflow)
    bc = bc+1
    name   (bc) = 'jetInflow'
    type   (bc) = 'SAT_FAR_FIELD'
    normDir(bc) =  2
    imin   (bc) =  iJet1
    imax   (bc) =  iJet2
    jmin   (bc) =  1
    jmax   (bc) =  1
    kmin   (bc) =  kJet1
    kmax   (bc) =  kJet2

    ! Target region
    bc = bc+1
    name   (bc) = 'targetRegion'
    type   (bc) = 'COST_TARGET'
    normDir(bc) =  0
    imin   (bc) =  1
    imax   (bc) = -1
    jmin   (bc) =  1
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! Control region
    bc = bc+1
    name   (bc) = 'controlRegion'
    type   (bc) = 'ACTUATOR'
    normDir(bc) =  0
    imin   (bc) =  1
    imax   (bc) = -1
    jmin   (bc) =  1
    jmax   (bc) = -1
    kmin   (bc) =  1
    kmax   (bc) = -1

    ! Overwrite number of boundary conditions.
    nbc = bc

    ! Open the file.
    iunit=11
    open(iunit,file="bc.dat")
    
    write (*,'(A)') ''
    write (*,'(A)') 'Writing boundary conditions'

    ! Write the header.
    write(iunit,'(1a87)') "# Name                 Type                  Grid normDir iMin iMax jMin jMax kMin kMax"
    write(iunit,'(1a87)') "# ==================== ===================== ==== ======= ==== ==== ==== ==== ==== ===="

    ! Write the boundary conditions.
    do i=1,nbc
       write(iunit,'(a2,2a21,8I5)') '  ',adjustl(name(i)),adjustl(type(i)),                  &
            1, normDir(i), imin(i), imax(i), jmin(i), jmax(i), kmin(i), kmax(i)
    end do

    ! Close the file.
    close(iunit)

    ! Clean up.
    deallocate(name, type, normDir, imin, imax, jmin, jmax, kmin, kmax)

    return
  end subroutine jetCrossFlowBC


  ! ================ !
  ! Target mollifier !
  ! ================ !
  subroutine jetCrossFlowTargetMollifier(imin, imax, jmin, jmax)

    ! <<< External modules >>>
    use MPI

    ! <<< Internal modules >>>
    use InputHelper, only : getOption
    use ErrorHandler, only : gracefulExit

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
    if (simulationFlags%predictionOnly) then
       print *, 'WARNING: target mollifier requires disable_adjoint_solver = false'
       stop
    end if

    ! Read the mollifier extents.
    call getRequiredOption("targer_mollifier_xmin", xmin)
    call getRequiredOption("targer_mollifier_xmax", xmax)
    call getRequiredOption("targer_mollifier_ymin", ymin)
    call getRequiredOption("targer_mollifier_ymax", ymax)

    ! Initialize the target mollifier.
    allocate(mollifier(nx,ny,nz,2))
    mollifier = 1.0_wp

    ! Hyperbolic tangent profile in y.
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             gridIndex = i+nx*(j-1+ny*(k-1))
             x = 2.0_wp * (grid%coordinates(gridIndex,2) - ymin) / (ymax - ymin) - 1.0_wp
             mollifier(i,j,k,1) = tanh(s * (x + 1.0_wp - 0.5_wp*r)) -                        &
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
             mollifier(i,j,k,2) = tanh(s * (x + 1.0_wp - 0.5_wp*r)) -                        &
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
  end subroutine jetCrossFlowTargetMollifier


  ! ================= !
  ! Control mollifier !
  ! ================= !
  subroutine jetCrossFlowControlMollifier(imin, imax, jmin, jmax)

    ! <<< External modules >>>
    use MPI

    ! <<< Internal modules >>>
    use InputHelper, only : getOption
    use ErrorHandler, only : gracefulExit

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, gridIndex
    integer, intent(out) :: imin, imax, jmin, jmax
    real(wp) :: xmin, xmax, ymin, ymax
    real(wp) :: x
    real(wp), dimension(:,:,:,:), allocatable :: mollifier
    real(wp), parameter :: s = 10.0_wp, r = 0.2_wp, eps = 1.0e-6_wp

    ! Make sure control mollifier is allocated.
    if (simulationFlags%predictionOnly) then
       print *, 'WARNING: control mollifier requires disable_adjoint_solver = false'
       stop
    end if

    ! Read the mollifier extents.
    call getRequiredOption("control_mollifier_xmin", xmin)
    call getRequiredOption("control_mollifier_xmax", xmax)
    call getRequiredOption("control_mollifier_ymin", ymin)
    call getRequiredOption("control_mollifier_ymax", ymax)

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
  end subroutine jetCrossFlowControlMollifier

end program jet_crossflow
