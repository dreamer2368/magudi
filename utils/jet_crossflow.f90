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
  type(t_SolverOptions) :: solverOptions
  type(t_SimulationFlags) :: simulationFlags
  integer :: i, numProcs, procRank, ierror, nx, ny, nz, nx_, ny_, nz_
  integer :: imin, imax, jmin, jmax, kmin, kmax
  integer :: iJet1, iJet2, kJet1, kJet2
  integer :: imin_sponge, imax_sponge, jmax_sponge
  real(KIND=8) :: xmini, xmaxi, ymaxi
  real(KIND=8) :: xmino, xmaxo, ymino, ymaxo
  real(KIND=8) :: zmin, zmax, jetDiameter, xJet
  character(len = STRING_LENGTH) :: filename, message
  integer, allocatable :: globalGridSizes(:,:)
  logical :: includeSandpaper, conformToJet

  ! Initialize MPI.
  call MPI_Init(ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  if (procRank == 0) then
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
     print *, '!    ** Parallel decomposition in z only **     !' 
     print *, '!                                               !'
     print *, '! ============================================= !'
     print *
  end if

  if (command_argument_count() > 1) then
     write(message, '(A)') "Usage: magudi [INPUT]"
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)') "High-performance Fortran-based adjoint optimization tool."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)')                                                                   &
          "Maximum of 1 INPUT file allowed."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     call cleanupErrorHandler()
     call MPI_Finalize(ierror)
     stop -1
  end if

  ! Get the input file name.
  if (command_argument_count() == 1) then
     call get_command_argument(1, filename)
     if (filename(1:1) == '-' .or. len_trim(filename) == 0) then
        write(message, '(A)') "No input file name was detected, using 'input'."
        call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
        filename = 'input'
     end if
  else
     write(message, '(A)') "No input file name was detected, using 'input'."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     filename = 'input'
  end if

  ! Parse options from the input file.
  call parseInputFile(filename)

  ! Generate and partition the grid.
  call jetCrossFlowGrid(imin_sponge, imax_sponge, jmax_sponge)

  ! Save the grid.
  call getRequiredOption("grid_file", filename)
  call region%saveData(QOI_GRID, filename)

  ! Compute normalized metrics, norm matrix, and Jacobian.
!!$  do i = 1, size(region%grids)
!!$     call region%grids(i)%update()
!!$  end do
!!$  call MPI_Barrier(MPI_COMM_WORLD, ierror)
!!$
!!$  ! Write out some useful information.
!!$  call region%reportGridDiagnostics()

  ! Generate the boundary conditions.
  call jetCrossFlowBC(imin_sponge, imax_sponge, jmax_sponge)

  ! Generate the initial condition and target state.
  call jetCrossFlowInitialCondition

  ! Save initial condition.
  call getRequiredOption("grid_file", filename)
  i = len_trim(filename)
  if (filename(i-3:i) == ".xyz") then
     filename = filename(:i-4) // ".ic.q"
  else
     filename = PROJECT_NAME // ".ic.q"
  end if
  call region%saveData(QOI_FORWARD_STATE, filename)

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
    !                      _~~~~~~~~~~~~~~~~~~~~~~~~~~_
    !                     |                            |
    ! j = 1  > ___________|                            |_________________|___:___|__________
    !
    !          ^          ^                            ^                 ^       ^         ^
    !         i=1        i=iTrip1                    i=iTrip2          i=iJet1  i=iJet2  i=nx

    ! <<< External modules >>>
    use MPI

    ! <<< Arguments >>>
    integer, intent(out) :: i1, i2, j2

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, n, nGrit, iTrip1, iTrip2, jetExtent
    integer, allocatable :: jetBox(:,:)
    real(wp) :: x, y, z, dx, dz, ytilde, r, delta, theta
    real(wp) :: tripLocation, tripWidth, tripHeight, totalHeight, peakHeight, valleyHeight
    real(wp) :: gritHeight, gritWidth, rnd, gauss, amp, sig, x0, y0, z0, z12, alpha
    logical :: includeGrit, stretchY
    character(len = STRING_LENGTH) :: key

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
    simulationFlags = region%simulationFlags
    solverOptions = region%solverOptions

    ! Store local grid size.
    nx_ = region%grids(1)%localSize(1)
    ny_ = region%grids(1)%localSize(2)
    nz_ = region%grids(1)%localSize(3)

    ! Get the grid partition
    imin = region%grids(1)%offset(1) + 1
    imax = region%grids(1)%offset(1) + nx_
    jmin = region%grids(1)%offset(2) + 1
    jmax = region%grids(1)%offset(2) + ny_
    kmin = region%grids(1)%offset(3) + 1
    kmax = region%grids(1)%offset(3) + nz_
    print *, 'Proc', procRank + 1, 'local size:',nx_, ny_, nz_

    ! Exit if processors are decomposed in x or y.
    if (nx_ /= nx .or. ny_ /= ny) then
       print *
       print *, 'jet_crossflow.f90 only implemented for parallel decomposition in z!'
       stop
    end if

    ! Should we stretch the mesh?
    stretchY = getOption('stretch_y', .false.)
    if (stretchY) r = 2.0_wp

    ! Generate the grid.
    do k = kmin, kmax
       do j = jmin, jmax
          do i = imin, imax

             ! Create X
             region%grids(1)%coordinates(grid_index(i,j,k), 1) =                             &
                  (xmaxo - xmino) * real(i-1, wp) / real(nx-1, wp) + xmino

             ! Create Y
             if (stretchY) then
                ytilde = real(ny-j, wp)/real(ny-1, wp)
                region%grids(1)%coordinates(grid_index(i,j,k), 2) =                          &
                     (ymaxo - ymino) * (1.0_wp - tanh(r * ytilde) / tanh(r)) + ymino
             else
                region%grids(1)%coordinates(grid_index(i,j,k), 2) =                          &
                     (ymaxo - ymino) * real(j-1,wp) / real(ny-1, wp) + ymino
             end if

             ! Create Z
             if (nz > 1) region%grids(1)%coordinates(grid_index(i,j,k), 3) =                &
                  (zmax - zmin - dz) * real(k-1, wp) / real(nz-1, wp) + zmin
          end do
       end do
    end do

    ! Generate the sandpaper.
    includeSandpaper = getOption("include_sandpaper", .false.)
    if (includeSandpaper) then
       write(key, '(A)') "sandpaper/"
       call getRequiredOption(trim(key)//"location", tripLocation)
       call getRequiredOption(trim(key)//"width", tripWidth)
       call getRequiredOption(trim(key)//"mean_plane_height", tripHeight)
       includeGrit = getOption(trim(key)//"include_grit", .false.)
       if (includeGrit) then
          call getRequiredOption(trim(key)//"grit_height", gritHeight)
          call getRequiredOption(trim(key)//"grit_width", gritWidth)
          call getRequiredOption(trim(key)//"number_of_particles", nGrit)
       end if

       ! Find the sandpaper extents.
       iTrip1 = 1; iTrip2 = nx
       j = jmin; k = kmin
       do i = imin, imax
          if (region%grids(1)%coordinates(grid_index(i,j,k), 1) < tripLocation)              &
               iTrip1 = i + 1
       end do
       do i = imax, imin, -1
          if (region%grids(1)%coordinates(grid_index(i,j,k),1) > tripLocation + tripWidth)   &
               iTrip2 = i - 1
       end do
       call MPI_Allreduce(MPI_IN_PLACE, iTrip1, 1, MPI_INTEGER, MPI_MAX,                     &
            MPI_COMM_WORLD, ierror)
       call MPI_Allreduce(MPI_IN_PLACE, iTrip2, 1, MPI_INTEGER, MPI_MIN,                     &
            MPI_COMM_WORLD, ierror)
       if (procRank == 0) then
          print *
          print *, 'Sandpaper extents: [', iTrip1, ',', iTrip2, ']'
       end if

       ! Deform the mesh to the sandpaper height.
       do k = kmin, kmax
          do i = imin, imax
             j = 1

             ! Create a smooth step.
             sig = 20.0_wp
             x = region%grids(1)%coordinates(grid_index(i,j,k), 1)
             region%grids(1)%coordinates(grid_index(i,j,k), 2) = 0.5_wp * tripHeight *       &
                  (tanh(sig * (x - tripLocation)) - tanh(sig*(x - tripLocation - tripWidth)))

             ! Shift grid points above (smoothly).
             do j = 2, ny
                ! Get unperturbed height variation.
                y0 = region%grids(1)%coordinates(grid_index(1,j,k), 2)
                delta = y0 - region%grids(1)%coordinates(grid_index(1,j-1,k), 2)

                ! Get the current height.
                y = region%grids(1)%coordinates(grid_index(i,j,k), 2)

                ! Adjust the current height.
                ytilde = region%grids(1)%coordinates(grid_index(i,j-1,k),2) + delta
                alpha = tanh(0.1_wp * (j - 2) / (ny - 2))
                region%grids(1)%coordinates(grid_index(i,j,k), 2) =                          &
                     ytilde * (1.0_wp - alpha) + y * alpha
             end do
          end do
       end do
  
       ! Embed the particles.
       If (includeGrit) then

          ! Standard deviation
          sig = gritWidth

          ! Loop through number of particles.
          do n = 1, nGrit

             ! Compute amplitude.
             amp = gritHeight
             call random_number(rnd)
             if (rnd < 0.5_wp) amp = -amp

             ! Get a random location.
             call random_number(rnd)
             x0 = tripLocation + 5.0_wp * sig + (tripWidth - 10.0_wp * sig) * rnd
             call random_number(rnd)
             z0 = (zmax - zmin - dz) * rnd + zmin
             if (nz == 1) z0 = 0.0_wp

             ! Modify the grid.
             do k = kmin, kmax
                do i = iTrip1, iTrip2
                   j = 1

                   ! Get the coordinates.
                   x = region%grids(1)%coordinates(grid_index(i,j,k), 1)
                   y = region%grids(1)%coordinates(grid_index(i,j,k), 2)
                   z = 0.0_wp
                   if (nz > 1) z = region%grids(1)%coordinates(grid_index(i,j,k), 3)
                      
                   ! Represent sandpaper particles as Gaussian.
                   gauss = amp * exp(-((x - x0)**2 / (2.0_wp * sig**2) +                     &
                        (z - z0)**2 / (2.0_wp * sig**2)))

                   ! Account for periodicity in z.
                   z12 = (zmax - zmin) - abs(z - z0)
                   If (nz > 1) gauss = gauss + amp *                                         &
                        exp(-((x-x0)**2 / (2.0_wp * sig**2) + z12**2 / (2.0_wp * sig**2)))

                   ! Update the vertical coordinate.
                   region%grids(1)%coordinates(grid_index(i,j,k), 2) = y + gauss
                   
                   ! Shift grid points above (smoothly).
                   do j = 2, ny
                      ! Get unperturbed height variation.
                      y0 = region%grids(1)%coordinates(grid_index(1,j,k), 2)
                      delta = y0 - region%grids(1)%coordinates(grid_index(1,j-1,k), 2)

                      ! Get the current height.
                      y = region%grids(1)%coordinates(grid_index(i,j,k), 2)
                      
                      ! Adjust the current height.
                      ytilde = region%grids(1)%coordinates(grid_index(i,j-1,k),2) + delta
                      alpha = tanh(4.0_wp * (j - 2) / (ny - 2))
                      region%grids(1)%coordinates(grid_index(i,j,k), 2) =                    &
                           ytilde * (1.0_wp - alpha) + y * alpha
                   end do

                end do
             end do
          end do
       end if !... includeGrit

       ! Output surface roughness dimentions.
       totalHeight = 0.0_wp
       peakHeight = 0.0_wp
       valleyHeight = huge(1.0_wp)
       j = 1
       do k = kmin, kmax
          do i = iTrip1, iTrip2
             y = region%grids(1)%coordinates(grid_index(i,j,k), 2)
             totalHeight = max(totalHeight, y)
             peakHeight = max(peakHeight, y - tripHeight)
             valleyHeight = min(valleyHeight, y - tripHeight)
          end do
       end do
       call MPI_Allreduce(MPI_IN_PLACE, totalHeight, 1, MPI_REAL8, MPI_MAX,                  &
            MPI_COMM_WORLD, ierror)
       call MPI_Allreduce(MPI_IN_PLACE, peakHeight, 1, MPI_REAL8, MPI_MAX,                   &
            MPI_COMM_WORLD, ierror)
       call MPI_Allreduce(MPI_IN_PLACE, valleyHeight, 1, MPI_REAL8, MPI_MIN,                 &
            MPI_COMM_WORLD, ierror)
       if (procRank == 0) then
          print *
          print *, 'Max surface height:', real(totalHeight, 4)
          print *, 'Peak height:', real(peakHeight, 4)
          print *, 'Valley height:', real(valleyHeight, 4)
       end if

    end if !...  includeSandpaper

    ! Find the jet extents.
    call getRequiredOption("jet_diameter", jetDiameter)
    call getRequiredOption("jet_position", xJet)
    iJet1 = nx; iJet2 = 1
    kJet1 = nz; kJet2 = 1
    j = 1
    do k = kmin, kmax
       do i = 1, nx
          x = region%grids(1)%coordinates(grid_index(i,j,k), 1)
          z = 0.0_wp
          if (nz > 1) z = region%grids(1)%coordinates(grid_index(i,j,k), 3)
          r = sqrt((x - xJet)**2 + (z - 0.5_wp * (zmax + zmin))**2)
          if (r <= 0.5_wp * jetDiameter) then
             iJet1 = min(iJet1, i)
             iJet2 = max(iJet2, i)
             kJet1 = min(kJet1, k)
             kJet2 = max(kJet2, k)
          end if
       end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE, iJet1, 1, MPI_INTEGER, MPI_MIN,                         &
         MPI_COMM_WORLD, ierror)
    call MPI_Allreduce(MPI_IN_PLACE, kJet1, 1, MPI_INTEGER, MPI_MIN,                         &
         MPI_COMM_WORLD, ierror)
    call MPI_Allreduce(MPI_IN_PLACE, iJet2, 1, MPI_INTEGER, MPI_MAX,                         &
         MPI_COMM_WORLD, ierror)
    call MPI_Allreduce(MPI_IN_PLACE, kJet2, 1, MPI_INTEGER, MPI_MAX,                         &
         MPI_COMM_WORLD, ierror)

    ! Conform the mesh to the jet perimeter.
    conformToJet = getOption("conform_grid_to_jet", .false.)

    ! Do not conform 2D grids.
    if (nz == 1 .and. conformToJet) then
       print *
       print *, 'Warning, conform_to_jet only enabled for 3D geometries.'
       conformToJet = .false.
    end if

    ! Stop if running in parallel.
    if (conformToJet .and. numProcs > 1) then
       print *
       print *, 'Warning, conform_to_jet currently implemented in serial!'
       stop
    end if
    if (conformToJet) then
       j = 1
       r = 0.495_wp * jetDiameter
       allocate(jetBox(iJet1-1:iJet2+1, kJet1-1:kJet2+1))
       jetBox = 0

       ! Find a point at `i` by shifting in `k`.
       do i = iJet1, iJet2
          k = kmin
          x = region%grids(1)%coordinates(grid_index(i,j,k), 1)

          ! First quadrant.
          delta = x - xJet
          theta = asin(delta / r)
          delta = r * cos(theta)
          z0 = 0.5_wp * (zmax + zmin) + delta
          delta = zmax - zmin
          do k = kJet1 - 1, kJet2 + 1
             ! Find the closest point to the jet perimeter.
             z = region%grids(1)%coordinates(grid_index(i,j,k), 3)
             if (abs(z - z0) < delta) then
                delta = abs(z - z0)
                n = k
             end if
          end do
          ! Shift the corresponding grid point.
          k = n
          if (jetBox(i, k) == 0) region%grids(1)%coordinates(grid_index(i,j,k), 3) = z0
          jetBox(i, k) = 1 ! Tag it.

          ! Second quadrant.
          delta = x - xJet
          theta = acos(delta / r)
          delta = r * sin(theta)
          z0 = 0.5_wp * (zmax + zmin) - delta
          delta = zmax - zmin
          do k = kJet1 - 1, kJet2 + 1
             ! Find the closest point to the jet perimeter.
             z = region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),3)
             if (abs(z - z0) < delta) then
                delta = abs(z - z0)
                n = k
             end if
          end do
          ! Shift the corresponding grid point.
          k = n
          if (jetBox(i, k) == 0) region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),3) = z0
          jetBox(i, k) = 1 ! Tag it.
       end do

       ! Find a point at `k` by shifting in `i`.
       do k = kJet1, kJet2
          i = 1
          z = region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),3)

          ! Third quadrant.
          delta = z - 0.5_wp * (zmax + zmin)
          theta = asin(delta / r)
          delta = r * cos(theta)
          x0 = xJet - delta
          delta = xmaxo - xmino
          do i = iJet1 - 1, iJet2 + 1
             ! Find the closest point to the jet perimeter.
             x = region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),1)
             if (abs(x - x0) < delta) then
                delta = abs(x - x0)
                n = i
             end if
          end do
          ! Shift the corresponding grid point.
          i = n
          if (jetBox(i, k) == 0) region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),1) = x0
          jetBox(i, k) = 1 ! Tag it.

          ! Fourth quadrant.
          delta = z - 0.5_wp * (zmax + zmin)
          theta = acos(delta / r)
          delta = r * sin(theta)
          x0 = xJet + delta
          delta = xmaxo - xmino
          do i = iJet1 - 1, iJet2 + 1
             ! Find the closest point to the jet perimeter.
             x = region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),1)
             if (abs(x - x0) < delta) then
                delta = abs(x - x0)
                n = i
             end if
          end do
          ! Shift the corresponding grid point.
          i = n
          if (jetBox(i, k) == 0) region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),1) = x0
          jetBox(i, k) = 1 ! Tag it.
       end do

       ! Smooth surrounding grid points.
       jetExtent = floor(0.5_wp * real(nz, wp)) - int(jetDiameter / dz)
       sig = 1.0_wp / 20.0_wp
       do k = kJet1 - 1, kJet2 + 1
          do i = iJet1 - 1, iJet2 + 1
             if (jetBox(i, k) == 1) then
                do n = 1, jetExtent

                   ! Shift left (x-direction).
                   j = 1
                   alpha = tanh(sig * real((n - 1), wp))
                   x = region%grids(1)%coordinates(grid_index(i-n,j,k), 1)
                   delta = region%grids(1)%coordinates(grid_index(i-n,j+1,k), 1) -           &
                        region%grids(1)%coordinates(grid_index(i-n+1,j+1,k), 1)
                   ytilde = region%grids(1)%coordinates(grid_index(i-n+1,j,k), 1) +          &
                        delta
                   region%grids(1)%coordinates(grid_index(i-n,j,k), 1) =                     &
                        ytilde * (1.0_wp - alpha) + x * alpha
                   do j = 2, jetExtent
                      alpha = tanh(sig * real((j - 2), wp))
                      x = region%grids(1)%coordinates(grid_index(i-n,j,k), 1)
                      region%grids(1)%coordinates(grid_index(i-n,j,k), 1) =                  &
                           region%grids(1)%coordinates(grid_index(i-n,j-1,k), 1) *           &
                           (1.0_wp - alpha) + x * alpha
                   end do

                   ! Shift right (x-direction).
                   j = 1
                   alpha = tanh(sig * real((n - 1), wp))
                   x = region%grids(1)%coordinates(grid_index(i+n,j,k), 1)
                   delta = region%grids(1)%coordinates(grid_index(i+n,j+1,k), 1) -           &
                        region%grids(1)%coordinates(grid_index(i+n-1,j+1,k),1)
                   ytilde = region%grids(1)%coordinates(grid_index(i+n-1,j,k), 1) +          &
                        delta
                   region%grids(1)%coordinates(grid_index(i+n,j,k), 1) =                     &
                        ytilde * (1.0_wp - alpha) + x * alpha
                   do j = 2, jetExtent
                      alpha =  tanh(sig * real((j - 1), wp))
                      x = region%grids(1)%coordinates(grid_index(i+n,j,k), 1)
                      region%grids(1)%coordinates(grid_index(i+n,j,k), 1) =                  &
                           region%grids(1)%coordinates(grid_index(i+n,j-1,k), 1) *           &
                           (1.0_wp - alpha) + x * alpha
                   end do

                   ! Shift up (z-direction).
                   j = 1
                   alpha = tanh(sig * real((n - 1), wp))
                   z = region%grids(1)%coordinates(grid_index(i,j,k-n), 3)
                   delta = region%grids(1)%coordinates(grid_index(i,j+1,k-n), 3) -           &
                        region%grids(1)%coordinates(grid_index(i,j+1,k-n+1), 3)
                   ytilde = region%grids(1)%coordinates(grid_index(i,j,k-n+1), 3) +          &
                        delta
                   region%grids(1)%coordinates(grid_index(i,j,k-n), 3) =                     &
                        ytilde * (1.0_wp - alpha) + z * alpha
                   do j = 2, jetExtent
                      alpha =  tanh(sig * real((j - 2), wp))
                      z = region%grids(1)%coordinates(grid_index(i,j,k-n), 3)
                      region%grids(1)%coordinates(grid_index(i,j,k-n), 3) =                  &
                           region%grids(1)%coordinates(grid_index(i,j-1,k-n), 3) *           &
                           (1.0_wp - alpha) + z * alpha
                   end do

                   ! Shift down (z-direction).
                   j = 1
                   alpha = tanh(sig * real((n - 1), wp))
                   z = region%grids(1)%coordinates(grid_index(i,j,k+n), 3)
                   delta = region%grids(1)%coordinates(grid_index(i,j+1,k+n), 3) -           &
                        region%grids(1)%coordinates(grid_index(i,j+1,k+n-1), 3)
                   ytilde = region%grids(1)%coordinates(grid_index(i,j,k+n-1), 3) +          &
                        delta
                   region%grids(1)%coordinates(grid_index(i,j,k+n), 3) =                     &
                        ytilde * (1.0_wp - alpha) + z * alpha
                   do j = 2, jetExtent
                      alpha =  tanh(sig * real((j - 2), wp))
                      z = region%grids(1)%coordinates(grid_index(i,j,k+n), 3)
                      region%grids(1)%coordinates(grid_index(i,j,k+n), 3) =                  &
                           region%grids(1)%coordinates(grid_index(i,j-1,k+n), 3) *           &
                           (1.0_wp - alpha) + z * alpha
                   end do

                end do
             end if
          end do
       end do
       deallocate(jetBox)

    end if !... conform to jet.

    ! Find the new jet extents.
    iJet1 = nx; iJet2 = 1
    kJet1 = nz; kJet2 = 1
    j = 1
    do k = kmin, kmax
       do i = 1, nx
          x = region%grids(1)%coordinates(grid_index(i,j,k), 1)
          z = 0.0_wp
          if (nz > 1) z = region%grids(1)%coordinates(grid_index(i,j,k), 3)
          r = sqrt((x - xJet)**2 + (z - 0.5_wp * (zmax + zmin))**2)
          if (r <= 0.5_wp * jetDiameter) then
             iJet1 = min(iJet1, i)
             iJet2 = max(iJet2, i)
             kJet1 = min(kJet1, k)
             kJet2 = max(kJet2, k)
          end if
       end do
    end do
    call MPI_Allreduce(MPI_IN_PLACE, iJet1, 1, MPI_INTEGER, MPI_MIN,                         &
         MPI_COMM_WORLD, ierror)
    call MPI_Allreduce(MPI_IN_PLACE, kJet1, 1, MPI_INTEGER, MPI_MIN,                         &
         MPI_COMM_WORLD, ierror)
    call MPI_Allreduce(MPI_IN_PLACE, iJet2, 1, MPI_INTEGER, MPI_MAX,                         &
         MPI_COMM_WORLD, ierror)
    call MPI_Allreduce(MPI_IN_PLACE, kJet2, 1, MPI_INTEGER, MPI_MAX,                         &
         MPI_COMM_WORLD, ierror)

    if (procRank == 0) then
       print *
       print *, 'Jet extents: [', iJet1, ',', iJet2, '] x [',kJet1, ',', kJet2,']'
    end if

    ! Find extents of outer region.
    j = 1; k = kmin
    i1 = 1
    do i = 1, nx
       if (region%grids(1)%coordinates(grid_index(i,j,k), 1) <= xmini) i1 = i + 1
    end do
    i2 = nx
    do i = nx,1,-1
       if (region%grids(1)%coordinates(grid_index(i,j,k), 1) >= xmaxi) i2 = i
    end do
    i = 1
    j2 = ny
    do j = ny, 1, -1
       if (region%grids(1)%coordinates(grid_index(i,j,k), 2) >= ymaxi) j2=j
    end do
    call MPI_Allreduce(MPI_IN_PLACE, i1, 1, MPI_INTEGER, MPI_MAX,                            &
         MPI_COMM_WORLD, ierror)
    call MPI_Allreduce(MPI_IN_PLACE, i2, 1, MPI_INTEGER, MPI_MIN,                            &
         MPI_COMM_WORLD, ierror)
    call MPI_Allreduce(MPI_IN_PLACE, j2, 1, MPI_INTEGER, MPI_MIN,                            &
         MPI_COMM_WORLD, ierror)

    if (procRank == 0) then
       print *
       print *, 'Interior domain extents: [',i1,',',i2,'] x [',1,',',j2,']'
       print*
    end if

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
    integer :: i, j, k, kk, nDimensions, nSpecies, H2, O2, N2, ierror
    real(wp) :: ratioOfSpecificHeats, crossflowVelocity, jetVelocity, x, y, z, y0, r, sig, a,&
         density, velocity(3), temperature, pressure, fuel, oxidizer, inert, YF0, Yo0, yDecay
    real(wp), dimension(:), allocatable :: Wi
    real(SCALAR_KIND), parameter :: spreadingRate = 0.094_wp
    character(len = STRING_LENGTH) :: velocityProfile, jetShape
    logical :: insideJet

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

    call MPI_Cartdim_get(region%grids(1)%comm, nDimensions, ierror)

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
    call getRequiredOption("jet_velocity_profile", velocityProfile)
    jetShape = getOption("jet_shape", "SQUARE")

    ! Mixture properties.
    call getRequiredOption("initial_fuel_mass_fraction", Yf0)
    call getRequiredOption("initial_oxidizer_mass_fraction", Yo0)

    ! Get the ratio of specific heats.
    ratioOfSpecificHeats = solverOptions%ratioOfSpecificHeats

    do k = kmin, kmax
       do j = jmin, jmax
          do i = imin, imax

             ! Get the local coordinates.
             x = region%grids(1)%coordinates(grid_index(i,j,k), 1)
             y = region%grids(1)%coordinates(grid_index(i,j,k), 2)
             y0= region%grids(1)%coordinates(grid_index(i,1,k), 2)
             z = 0.0_wp
             if (nz > 1) z = region%grids(1)%coordinates(grid_index(i,j,k), 3)

             ! Initialize the mixture.
             velocity = 0.0_wp
             pressure = 1.0_wp / ratioOfSpecificHeats
             temperature =  1.0_wp / (ratioOfSpecificHeats - 1.0_wp)
             fuel = 0.0_wp
             oxidizer = Yo0

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

             ! Make sure wall velocities are zero.
             if (j == 1) velocity = 0.0_wp

             ! Jet conditions.
             r = sqrt((x - xJet)**2 + (z - 0.5_wp * (zmax + zmin))**2)
             insideJet = .false.
             select case(trim(jetShape))
             case('ROUND')
                if (r <= 1.0_wp * jetDiameter) insideJet = .true.
             case('SQUARE')
                if (i >= iJet1 .and. i <= iJet2 .and. k >= kJet1 .and. k <= kJet2)           &
                     insideJet = .true.
             case ('NONE')
                ! Nothing to do.
             case default
                print *, 'Error: Unknown jet shape!'
                stop
             end select

             !if (insideJet) then
                yDecay = max(0.0_wp, 0.5_wp * (1.0_wp + tanh((4.0_wp - y / jetDiameter))))

                ! Fuel stream.
                sig = 6.0_wp
                fuel = Yf0 * yDecay
                if (trim(jetShape) == 'ROUND') then
                   fuel = fuel * 0.5_wp * (                                                  &
                        tanh(sig * (r + 0.5_wp * jetDiameter) / jetDiameter) -               &
                        tanh(sig * (r - 0.5_wp * jetDiameter) / jetDiameter))
                else if (trim(jetShape) == 'SQUARE') then
                   fuel = fuel * 0.5_wp * (                                                  &
                        tanh(sig * (x - xJet + 0.5_wp * jetDiameter) / jetDiameter) -        &
                        tanh(sig * (x - xJet - 0.5_wp * jetDiameter) / jetDiameter))  *      &
                        0.5_wp * (                                                           &
                        tanh(sig * (z - 0.5_wp * (zmax + zmin) + 0.5_wp * jetDiameter) /     &
                        jetDiameter) - tanh(sig * (z - 0.5_wp * (zmax + zmin) - 0.5_wp *     &
                        jetDiameter) / jetDiameter))
                end if
                oxidizer = (1.0_wp - fuel) * Yo0

                ! Jet velocity.
                sig = 6.0_wp
                velocity(2) = jetVelocity * yDecay
                select case(trim(velocityProfile))
                case('TANH')

                   if (trim(jetShape) == 'ROUND') then
                      velocity(2) = velocity(2) * 0.5_wp * (                                 &
                           tanh(sig * (r + 0.5_wp * jetDiameter) / jetDiameter) -            &
                           tanh(sig * (r - 0.5_wp * jetDiameter) / jetDiameter))
                   else if (trim(jetShape) == 'SQUARE') then
                      velocity(2) = velocity(2) * 0.5_wp * (                                 &
                           tanh(sig * (x - xJet + 0.5_wp * jetDiameter) / jetDiameter) -     &
                           tanh(sig * (x - xJet - 0.5_wp * jetDiameter) / jetDiameter))  *   &
                           0.5_wp * (                                                        &
                           tanh(sig * (z - 0.5_wp * (zmax + zmin) + 0.5_wp * jetDiameter) /  &
                           jetDiameter) - tanh(sig * (z - 0.5_wp * (zmax + zmin) - 0.5_wp *  &
                           jetDiameter) / jetDiameter))
                end if

                case ('POISEUILLE')

                   if (trim(jetShape) == 'ROUND') then
                      velocity(2) = max(0.0_wp,                                              &
                           velocity(2) * 2.0_wp * (1.0_wp - (r / (0.5_wp * jetDiameter))**2))
                   else if (trim(jetShape) == 'SQUARE') then
                      velocity(2) = 2.0_wp * velocity(2) * &
                           max(0.0_wp, (1.0_wp - ((x - xJet) / (0.5_wp * jetDiameter))**2)) *&
                           max(0.0_wp, (1.0_wp - ((z - 0.5_wp * (zmax + zmin)) /             &
                           (0.5_wp * jetDiameter))**2))
                   end if

                case ('SELF_SIMILAR')

                   eta = (x - xJet)
                   a = (sqrt(2.0_wp) - 1.0_wp) / (0.5_wp * jetDiameter)**2
                   velocity(2) = velocity(2) * (1.0_wp + a * eta**2) ** (-2.0_wp)
                   velocity(1) = 0.5_wp * jetVelocity * (eta - a * eta**3) /                 &
                        (1.0_wp + a * eta**2)**2 * yDecay

                case default

                   velocity = 0.0_wp

                end select

                !end if !... if (insideJet)

             ! Correct species mass fractions.
             inert = 1.0_wp - fuel - oxidizer
             if (inert < 0.0_wp) then
                if (procRank == 0) print *, 'Something is wrong with the species mass fraction!'
                oxidizer = oxidizer + inert
             end if

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
             region%states(1)%conservedVariables(grid_index(i,j,k), 1) = density
             region%states(1)%conservedVariables(grid_index(i,j,k), 2:nDimensions+1) =       &
                  density * velocity(1:nDimensions)
             region%states(1)%conservedVariables(grid_index(i,j,k),nDimensions+2) =          &
                  pressure / (ratioOfSpecificHeats - 1.0_wp) + 0.5_wp * density *            &
                  sum(velocity ** 2)
             if (nSpecies.gt.0) region%states(1)%conservedVariables(grid_index(i,j,k),       &
                  nDimensions+2+H2) = density * fuel
             if (nSpecies.gt.1) region%states(1)%conservedVariables(grid_index(i,j,k),       &
                  nDimensions+2+O2) = density * oxidizer

          end do !... do i = imin, imax
       end do !... do j = jmin, jmax
    end do !... do k = kmin, kmax

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
    character(len = 22), allocatable, dimension(:) :: name, type

    ! Only root process writes boundary conditions
    if (procRank /= 0) return

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

    ! Bottom wall
    bc = bc+1
    name   (bc) = 'bottomWall'
    type   (bc) = 'SAT_ISOTHERMAL_WALL'
    normDir(bc) =  2
    imin   (bc) =  1
    imax   (bc) = -1
    jmin   (bc) =  1
    jmax   (bc) =  1
    kmin   (bc) =  1
    kmax   (bc) = -1


    ! Jet inflow at the bottom wall
    bc = bc+1
    name   (bc) = 'jetInflow'
    type   (bc) = 'SAT_FAR_FIELD'
    normDir(bc) =  2
    imin   (bc) =  1
    imax   (bc) = -1
    jmin   (bc) =  1
    jmax   (bc) =  1
    kmin   (bc) =  1
    kmax   (bc) = -1

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
    do i = 1, nbc
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
    integer :: i, j, k
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
             x = 2.0_wp * (region%grids(1)%coordinates(grid_index(i,j,k), 2) - ymin) /       &
                  (ymax - ymin) - 1.0_wp
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
             x = 2.0_wp * (region%grids(1)%coordinates(grid_index(i,j,k), 1) - xmin) /       &
                  (xmax - xmin) - 1.0_wp
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
             region%grids(1)%targetMollifier(grid_index(i,j,k),1) = mollifier(i,j,k,1) *     &
                  mollifier(i,j,k,2)
          end do
       end do
    end do

    ! Find initial extents in y.
    jmin = 1; jmax = 1
    i = 1; k = 1
    do j = 1, ny
       if (region%grids(1)%coordinates(grid_index(i,j,k), 2) <= ymin) jmin = j
       if (region%grids(1)%coordinates(grid_index(i,j,k), 2) < ymax) jmax = j + 1
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

    if (procRank == 0) then
       print *
       print *, 'Target mollifier extents: [',imin,',',imax,'] x [',jmin,',',jmax,']'
       print *
    end if

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
    integer :: i, j, k
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
             x = 2.0_wp * (region%grids(1)%coordinates(grid_index(i,j,k), 2) - ymin) /       &
                  (ymax - ymin) - 1.0_wp
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
             x = 2.0_wp * (region%grids(1)%coordinates(grid_index(i,j,k), 1) - xmin) /       &
                  (xmax - xmin) - 1.0_wp
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
             region%grids(1)%controlMollifier(grid_index(i,j,k), 1) = mollifier(i,j,k,1) *   &
                  mollifier(i,j,k,2)
          end do
       end do
    end do

    ! Find initial extents in y.
    jmin = 1; jmax = 1
    i = 1; k = 1
    do j = 1, ny
       if (region%grids(1)%coordinates(grid_index(i,j,k), 2) <= ymin) jmin = j
       if (region%grids(1)%coordinates(grid_index(i,j,k), 2) < ymax) jmax = j+1
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

    if (procRank == 0) then
       print *
       print *, 'Control mollifier extents: [',imin,',',imax,'] x [',jmin,',',jmax,']'
       print *
    end if

    ! Clean up.
    deallocate(mollifier)

    return
  end subroutine jetCrossFlowControlMollifier


  ! Returns the lexicographic grid index
  function grid_index(i, j, k) result(l)
    implicit none

    integer, intent(in) :: i, j, k
    integer :: l
    
    l = i - region%grids(1)%offset(1) + region%grids(1)%localSize(1) *                       &
         (j - 1 - region%grids(1)%offset(2) + region%grids(1)%localSize(2) *                 &
         (k - 1 - region%grids(1)%offset(3)))
    
    return
  end function grid_index

end program jet_crossflow
