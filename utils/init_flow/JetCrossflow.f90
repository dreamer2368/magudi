#include "config.h"

module JetCrossflow_mod

  implicit none
  public

  ! <<< Global parameters >>>
  integer :: nx, ny, nz
  integer, dimension(3) :: iStart, iEnd, interiorStart, interiorEnd,                         &
       targetStart, targetEnd, controlStart, controlEnd
  integer :: iJet1, iJet2, kJet1, kJet2
  real(SCALAR_KIND) :: jetDiameter, jetPosition
  logical :: includeSandpaper, conformToJet

end module JetCrossflow_mod

subroutine jetCrossFlowGrid(region)

  ! <<< Private members >>>
  use JetCrossflow_mod

  ! <<< External modules >>>
  use MPI

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  ! <<< Derived types >>>
  use Region_mod, only : t_Region

  implicit none

  ! Jet in crossflow schematic and the corresponding indeces 
  !
  !
  !                      _~~~~~~~~~~~~~~~~~~~~~~~~~~_
  ! j = 1  > ___________|                            |_________________|___:___|__________
  !
  !          ^          ^                            ^                 ^       ^         ^
  !         i=1        i=iTrip1                    i=iTrip2          i=iJet1  i=iJet2  i=nx

  ! <<< Arguments >>>
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, n, nGrit, iTrip1, iTrip2, jetExtent, nDimensions, procRank, numProcs,  &
       ierror
  integer, allocatable :: jetBox(:,:)
  real(wp), dimension(3) :: minLengthInterior, minLengthOuter,                               &
       maxLengthInterior, maxLengthOuter
  real(wp) :: x, y, z, gridSpacing(3), ytilde, r, delta, theta
  real(wp) :: tripLocation, tripWidth, tripHeight, totalHeight, peakHeight, valleyHeight
  real(wp) :: gritHeight, gritWidth, rnd, gauss, amp, sig, x0, y0, z0, z12, alpha
  logical :: includeGrit, stretchY
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

  ! Exit if processors are decomposed in x or y.
  if (region%grids(1)%localSize(1) /= nx .or. region%grids(1)%localSize(2) /= ny) then
     write(message, '(A)')                                                                   &
          "JetCrossflow.f90 only implemented for parallel decomposition in z!"
     call gracefulExit(region%comm, message)
  end if

  ! Should we stretch the mesh?
  stretchY = getOption('stretch_y', .false.)
  if (stretchY) r = 2.0_wp

  ! Generate the grid.
  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

           ! Create X.
           region%grids(1)%coordinates(grid_index(i,j,k), 1) =                               &
                (maxLengthOuter(1) - minLengthOuter(1)) * real(i-1, wp) / real(nx-1, wp) +   &
                minLengthOuter(1)

           ! Create Y.
           if (stretchY) then
              ytilde = real(ny-j, wp)/real(ny-1, wp)
              region%grids(1)%coordinates(grid_index(i,j,k), 2) =                            &
                   (maxLengthOuter(2) - minLengthOuter(2)) * (1.0_wp - tanh(r * ytilde) /    &
                   tanh(r)) + minLengthOuter(2)
           else
              region%grids(1)%coordinates(grid_index(i,j,k), 2) =                            &
                   (maxLengthOuter(2) - minLengthOuter(2)) * real(j-1,wp) / real(ny-1, wp) + &
                   minLengthOuter(2)
           end if

           ! Create Z.
           if (nz > 1) region%grids(1)%coordinates(grid_index(i,j,k), 3) =                   &
                (maxLengthOuter(3) - minLengthOuter(3) - gridSpacing(3)) * real(k-1, wp) /   &
                real(nz-1, wp) + minLengthOuter(3)
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
     j = iStart(2); k = iStart(3)
     do i = iStart(1), iEnd(1)
        if (region%grids(1)%coordinates(grid_index(i,j,k), 1) < tripLocation)                &
             iTrip1 = i + 1
     end do
     do i = iEnd(1), iStart(1), -1
        if (region%grids(1)%coordinates(grid_index(i,j,k),1) > tripLocation + tripWidth)     &
             iTrip2 = i - 1
     end do
     call MPI_Allreduce(MPI_IN_PLACE, iTrip1, 1, MPI_INTEGER, MPI_MAX,                       &
          MPI_COMM_WORLD, ierror)
     call MPI_Allreduce(MPI_IN_PLACE, iTrip2, 1, MPI_INTEGER, MPI_MIN,                       &
          MPI_COMM_WORLD, ierror)
     if (procRank == 0) then
        write (*,'(A)') ''
        write (*,'(A, i0.0, A, i0.0, A)') 'Sandpaper extents: [', iTrip1, ', ', iTrip2, ']'
     end if

     ! Deform the mesh to the sandpaper height.
     do k = iStart(3), iEnd(3)
        do i = iStart(1), iEnd(1)
           j = 1

           ! Create a smooth step.
           sig = 20.0_wp
           x = region%grids(1)%coordinates(grid_index(i,j,k), 1)
           region%grids(1)%coordinates(grid_index(i,j,k), 2) = 0.5_wp * tripHeight *         &
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
              alpha = tanh(0.1_wp * real(j - 2, wp) / real(ny - 2, wp))
              region%grids(1)%coordinates(grid_index(i,j,k), 2) =                            &
                   ytilde * (1.0_wp - alpha) + y * alpha
           end do
        end do
     end do
  
     ! Embed the particles.
     If (includeGrit) then

        ! Standard deviation
        sig = gritWidth

        ! Loop through number of particles.
        write (*,'(A)') ''
        write (*,'(A, i0.0, A)') 'Adding ', nGrit, ' particles...'
        do n = 1, nGrit

           ! Compute amplitude.
           amp = gritHeight
           call random_number(rnd)
           if (rnd < 0.5_wp) amp = -amp

           ! Get a random location.
           call random_number(rnd)
           x0 = tripLocation + 5.0_wp * sig + (tripWidth - 10.0_wp * sig) * rnd
           call random_number(rnd)
           z0 = (maxLengthOuter(3) - minLengthOuter(3) - gridSpacing(3)) * rnd +             &
                minLengthOuter(3)
           if (nz == 1) z0 = 0.0_wp

           ! Modify the grid.
           do k = iStart(3), iEnd(3)
              do i = iTrip1, iTrip2
                 j = 1

                 ! Get the coordinates.
                 x = region%grids(1)%coordinates(grid_index(i,j,k), 1)
                 y = region%grids(1)%coordinates(grid_index(i,j,k), 2)
                 z = 0.0_wp
                 if (nz > 1) z = region%grids(1)%coordinates(grid_index(i,j,k), 3)
                      
                 ! Represent sandpaper particles as Gaussian.
                 gauss = amp * exp(-((x - x0)**2 / (2.0_wp * sig**2) +                       &
                      (z - z0)**2 / (2.0_wp * sig**2)))

                 ! Account for periodicity in z.
                 z12 = (maxLengthOuter(3) - minLengthOuter(3)) - abs(z - z0)
                 If (nz > 1) gauss = gauss + amp *                                           &
                      exp(-((x - x0)**2 / (2.0_wp * sig**2) + z12**2 / (2.0_wp * sig**2)))

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
                    alpha = tanh(4.0_wp * real(j - 2, wp) / real(ny - 2, wp))
                    region%grids(1)%coordinates(grid_index(i,j,k), 2) =                      &
                         ytilde * (1.0_wp - alpha) + y * alpha
                 end do

              end do
           end do

           ! Output the progress.
           if (mod(real(n, wp), real(nGrit, wp) / 10.0_wp) == 0.0_wp .and. procRank == 0)    &
                write(*,'(f5.1, A)') real(n, wp) / real(nGrit, wp) * 100.0_wp, '% complete'

        end do
     end if !... includeGrit

     ! Output surface roughness dimentions.
     totalHeight = 0.0_wp
     peakHeight = 0.0_wp
     valleyHeight = huge(1.0_wp)
     j = 1
     do k = iStart(3), iEnd(3)
        do i = iTrip1, iTrip2
           y = region%grids(1)%coordinates(grid_index(i,j,k), 2)
           totalHeight = max(totalHeight, y)
           peakHeight = max(peakHeight, y - tripHeight)
           valleyHeight = min(valleyHeight, y - tripHeight)
        end do
     end do
     call MPI_Allreduce(MPI_IN_PLACE, totalHeight, 1, MPI_REAL8, MPI_MAX,                    &
          MPI_COMM_WORLD, ierror)
     call MPI_Allreduce(MPI_IN_PLACE, peakHeight, 1, MPI_REAL8, MPI_MAX,                     &
          MPI_COMM_WORLD, ierror)
     call MPI_Allreduce(MPI_IN_PLACE, valleyHeight, 1, MPI_REAL8, MPI_MIN,                   &
          MPI_COMM_WORLD, ierror)
     if (procRank == 0) then
        write (*,'(A)') ''
        print *, 'Max surface height:', real(totalHeight, 4)
        print *, 'Peak height:', real(peakHeight, 4)
        print *, 'Valley height:', real(valleyHeight, 4)
     end if

  end if !...  includeSandpaper

  ! Find the jet extents.
  call getRequiredOption("jet_diameter", jetDiameter)
  call getRequiredOption("jet_position", jetPosition)
  iJet1 = nx; iJet2 = 1
  kJet1 = nz; kJet2 = 1
  j = 1
  do k = iStart(3), iEnd(3)
     do i = 1, nx
        x = region%grids(1)%coordinates(grid_index(i,j,k), 1)
        z = 0.0_wp
        if (nz > 1) z = region%grids(1)%coordinates(grid_index(i,j,k), 3)
        r = sqrt((x - jetPosition)**2 + z**2)
        if (r <= 0.5_wp * jetDiameter) then
           iJet1 = min(iJet1, i)
           iJet2 = max(iJet2, i)
           kJet1 = min(kJet1, k)
           kJet2 = max(kJet2, k)
        end if
     end do
  end do
  call MPI_Allreduce(MPI_IN_PLACE, iJet1, 1, MPI_INTEGER, MPI_MIN,                           &
       MPI_COMM_WORLD, ierror)
  call MPI_Allreduce(MPI_IN_PLACE, kJet1, 1, MPI_INTEGER, MPI_MIN,                           &
       MPI_COMM_WORLD, ierror)
  call MPI_Allreduce(MPI_IN_PLACE, iJet2, 1, MPI_INTEGER, MPI_MAX,                           &
       MPI_COMM_WORLD, ierror)
  call MPI_Allreduce(MPI_IN_PLACE, kJet2, 1, MPI_INTEGER, MPI_MAX,                           &
       MPI_COMM_WORLD, ierror)

  ! Conform the mesh to the jet perimeter.
  conformToJet = getOption("conform_grid_to_jet", .false.)

  ! Do not conform 2D grids.
  if (nz == 1 .and. conformToJet) then
     write (*,'(A)') ''
     write (*,'(A)') 'Warning, conform_to_jet only enabled for 3D geometries.'
     conformToJet = .false.
  end if

  ! Stop if running in parallel.
  if (conformToJet .and. numProcs > 1) then
     write (*,'(A)') ''
     write(message, '(A)') 'Warning, conform_to_jet currently implemented in serial!'
     call gracefulExit(region%comm, message)
  end if
  if (conformToJet) then
     j = 1
     r = 0.495_wp * jetDiameter
     allocate(jetBox(iJet1-1:iJet2+1, kJet1-1:kJet2+1))
     jetBox = 0

     ! Find a point at `i` by shifting in `k`.
     do i = iJet1, iJet2
        k = iStart(3)
        x = region%grids(1)%coordinates(grid_index(i,j,k), 1)

        ! First quadrant.
        delta = x - jetPosition
        theta = asin(delta / r)
        delta = r * cos(theta)
        z0 = delta
        delta = maxLengthOuter(3) - minLengthOuter(3)
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
        jetBox(i, k) = 1 !... Tag it.

        ! Second quadrant.
        delta = x - jetPosition
        theta = acos(delta / r)
        delta = r * sin(theta)
        z0 = - delta
        delta = maxLengthOuter(3) - minLengthOuter(3)
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
        jetBox(i, k) = 1 !... Tag it.
     end do

     ! Find a point at `k` by shifting in `i`.
     do k = kJet1, kJet2
        i = 1
        z = region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),3)

        ! Third quadrant.
        delta = z
        theta = asin(delta / r)
        delta = r * cos(theta)
        x0 = jetPosition - delta
        delta = maxLengthOuter(1) - minLengthOuter(1)
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
        jetBox(i, k) = 1 !... Tag it.

        ! Fourth quadrant.
        delta = z
        theta = acos(delta / r)
        delta = r * sin(theta)
        x0 = jetPosition + delta
        delta = maxLengthOuter(1) - minLengthOuter(1)
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
        jetBox(i, k) = 1 !... Tag it.
     end do

     ! Smooth surrounding grid points.
     jetExtent = floor(0.5_wp * real(nz, wp)) - int(jetDiameter / gridSpacing(3))
     sig = 1.0_wp / 20.0_wp
     do k = kJet1 - 1, kJet2 + 1
        do i = iJet1 - 1, iJet2 + 1
           if (jetBox(i, k) == 1) then
              do n = 1, jetExtent

                 ! Shift left (x-direction).
                 j = 1
                 alpha = tanh(sig * real((n - 1), wp))
                 x = region%grids(1)%coordinates(grid_index(i-n,j,k), 1)
                 delta = region%grids(1)%coordinates(grid_index(i-n,j+1,k), 1) -             &
                      region%grids(1)%coordinates(grid_index(i-n+1,j+1,k), 1)
                 ytilde = region%grids(1)%coordinates(grid_index(i-n+1,j,k), 1) +            &
                      delta
                 region%grids(1)%coordinates(grid_index(i-n,j,k), 1) =                       &
                      ytilde * (1.0_wp - alpha) + x * alpha
                 do j = 2, jetExtent
                    alpha = tanh(sig * real((j - 2), wp))
                    x = region%grids(1)%coordinates(grid_index(i-n,j,k), 1)
                    region%grids(1)%coordinates(grid_index(i-n,j,k), 1) =                    &
                         region%grids(1)%coordinates(grid_index(i-n,j-1,k), 1) *             &
                         (1.0_wp - alpha) + x * alpha
                 end do

                 ! Shift right (x-direction).
                 j = 1
                 alpha = tanh(sig * real((n - 1), wp))
                 x = region%grids(1)%coordinates(grid_index(i+n,j,k), 1)
                 delta = region%grids(1)%coordinates(grid_index(i+n,j+1,k), 1) -             &
                      region%grids(1)%coordinates(grid_index(i+n-1,j+1,k),1)
                 ytilde = region%grids(1)%coordinates(grid_index(i+n-1,j,k), 1) +            &
                      delta
                 region%grids(1)%coordinates(grid_index(i+n,j,k), 1) =                       &
                      ytilde * (1.0_wp - alpha) + x * alpha
                 do j = 2, jetExtent
                    alpha =  tanh(sig * real((j - 1), wp))
                    x = region%grids(1)%coordinates(grid_index(i+n,j,k), 1)
                    region%grids(1)%coordinates(grid_index(i+n,j,k), 1) =                    &
                         region%grids(1)%coordinates(grid_index(i+n,j-1,k), 1) *             &
                         (1.0_wp - alpha) + x * alpha
                 end do

                 ! Shift up (z-direction).
                 j = 1
                 alpha = tanh(sig * real((n - 1), wp))
                 z = region%grids(1)%coordinates(grid_index(i,j,k-n), 3)
                 delta = region%grids(1)%coordinates(grid_index(i,j+1,k-n), 3) -             &
                      region%grids(1)%coordinates(grid_index(i,j+1,k-n+1), 3)
                 ytilde = region%grids(1)%coordinates(grid_index(i,j,k-n+1), 3) +            &
                      delta
                 region%grids(1)%coordinates(grid_index(i,j,k-n), 3) =                       &
                      ytilde * (1.0_wp - alpha) + z * alpha
                 do j = 2, jetExtent
                    alpha =  tanh(sig * real((j - 2), wp))
                    z = region%grids(1)%coordinates(grid_index(i,j,k-n), 3)
                    region%grids(1)%coordinates(grid_index(i,j,k-n), 3) =                    &
                         region%grids(1)%coordinates(grid_index(i,j-1,k-n), 3) *             &
                         (1.0_wp - alpha) + z * alpha
                 end do

                 ! Shift down (z-direction).
                 j = 1
                 alpha = tanh(sig * real((n - 1), wp))
                 z = region%grids(1)%coordinates(grid_index(i,j,k+n), 3)
                 delta = region%grids(1)%coordinates(grid_index(i,j+1,k+n), 3) -             &
                      region%grids(1)%coordinates(grid_index(i,j+1,k+n-1), 3)
                 ytilde = region%grids(1)%coordinates(grid_index(i,j,k+n-1), 3) +            &
                      delta
                 region%grids(1)%coordinates(grid_index(i,j,k+n), 3) =                       &
                      ytilde * (1.0_wp - alpha) + z * alpha
                 do j = 2, jetExtent
                    alpha =  tanh(sig * real((j - 2), wp))
                    z = region%grids(1)%coordinates(grid_index(i,j,k+n), 3)
                    region%grids(1)%coordinates(grid_index(i,j,k+n), 3) =                    &
                         region%grids(1)%coordinates(grid_index(i,j-1,k+n), 3) *             &
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
  do k = iStart(3), iEnd(3)
     do i = 1, nx
        x = region%grids(1)%coordinates(grid_index(i,j,k), 1)
        z = 0.0_wp
        if (nz > 1) z = region%grids(1)%coordinates(grid_index(i,j,k), 3)
        r = sqrt((x - jetPosition)**2 + z**2)
        if (r <= 0.5_wp * jetDiameter) then
           iJet1 = min(iJet1, i)
           iJet2 = max(iJet2, i)
           kJet1 = min(kJet1, k)
           kJet2 = max(kJet2, k)
        end if
     end do
  end do
  call MPI_Allreduce(MPI_IN_PLACE, iJet1, 1, MPI_INTEGER, MPI_MIN,                           &
       MPI_COMM_WORLD, ierror)
  call MPI_Allreduce(MPI_IN_PLACE, kJet1, 1, MPI_INTEGER, MPI_MIN,                           &
       MPI_COMM_WORLD, ierror)
  call MPI_Allreduce(MPI_IN_PLACE, iJet2, 1, MPI_INTEGER, MPI_MAX,                           &
       MPI_COMM_WORLD, ierror)
  call MPI_Allreduce(MPI_IN_PLACE, kJet2, 1, MPI_INTEGER, MPI_MAX,                           &
       MPI_COMM_WORLD, ierror)

  if (procRank == 0) then
     print *
     print *, 'Jet extents: [', iJet1, ',', iJet2, '] x [',kJet1, ',', kJet2,']'
  end if

  ! Find extents of interior region (used to define sponge extents).
  interiorStart = 1; interiorEnd = -1
  j = iStart(2); k = iStart(3)
  interiorStart(1) = nx
  do i = iEnd(1), iStart(1), -1
     if (region%grids(1)%coordinates(grid_index(i,j,k), 1) >= minLengthInterior(1))          &
          interiorStart(1) = i
  end do
  interiorEnd(1) = 1
  do i = iStart(1), iEnd(1)
     if (region%grids(1)%coordinates(grid_index(i,j,k), 1) <= maxLengthInterior(1))          &
          interiorEnd(1) = i
  end do
  i = iStart(1)
  interiorStart(2) = ny
  do j = iEnd(2), iStart(2), -1
     if (region%grids(1)%coordinates(grid_index(i,j,k), 2) >= minLengthInterior(2))          &
          interiorStart(2) = j
  end do
  interiorEnd(2) = 1
  do j = iStart(2), iEnd(2)
     if (region%grids(1)%coordinates(grid_index(i,j,k), 2) <= maxLengthInterior(2))          &
          interiorEnd(2) = j
  end do
  call MPI_Allreduce(MPI_IN_PLACE, interiorStart(1), 1, MPI_INTEGER, MPI_MIN,                &
       MPI_COMM_WORLD, ierror)
  call MPI_Allreduce(MPI_IN_PLACE, interiorEnd(1), 1, MPI_INTEGER, MPI_MAX,                  &
       MPI_COMM_WORLD, ierror)
  call MPI_Allreduce(MPI_IN_PLACE, interiorStart(2), 1, MPI_INTEGER, MPI_MIN,                &
       MPI_COMM_WORLD, ierror)
  call MPI_Allreduce(MPI_IN_PLACE, interiorEnd(2), 1, MPI_INTEGER, MPI_MAX,                  &
       MPI_COMM_WORLD, ierror)

  if (procRank == 0) then
     write (*,'(A)') ''
     write (*,'(A, i0.0, A, i0.0, A, i0.0, A, i0.0, A)') 'Interior domain extents: [',       &
          interiorStart(1), ', ', interiorEnd(1), '] x [',                                   &
          interiorStart(2), ', ', interiorEnd(2), ']'
  end if

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

end subroutine jetCrossFlowGrid

subroutine jetCrossFlowQ(region)

  ! <<< Private members >>>
  use JetCrossflow_mod

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
  integer :: i, j, k, kk, nDimensions, nSpecies, H2, O2, N2, procRank, numProcs, ierror
  real(wp) :: ratioOfSpecificHeats, crossflowVelocity, jetVelocity, x, y, z, y0, r, sig, a,  &
       density, velocity(3), temperature, pressure, fuel, oxidizer, inert, YF0, Yo0, yDecay
  real(wp), dimension(:), allocatable :: Wi
  real(SCALAR_KIND), parameter :: spreadingRate = 0.094_wp
  character(len = STRING_LENGTH) :: velocityProfile, jetShape, message
  logical :: insideJet

  ! Solution to Blasius boundary layer
  real(wp) :: blasius0, blasius1, delta, eta, xx, x0, Re_c
  real(wp) :: f2l, f2r, f0l, f0r
  real(wp), dimension(0:9) :: by0 = (/                                                       &
       0.0_8, 0.165571818583440_8, 0.650024518764203_8, 1.39680822972500_8,                  &
       2.30574664618049_8, 3.28327391871370_8, 4.27962110517696_8,                           &
       5.27923901129384_8, 6.27921363832835_8, 7.27921257797747_8 /)
  real(wp), dimension(0:9) :: by1 = (/                                                       &
       0.0_8, 0.329780063306651_8, 0.629765721178679_8, 0.84604458266019_8,                  &
       0.95551827831671_8, 0.99154183259084_8, 0.99897290050990_8,                           &
       0.9999216098795_8, 0.99999627301467_8, 0.99999989265063_8 /)
  real(wp), dimension(0:9) :: by2 = (/                                                       &
       0.332057384255589_8, 0.323007152241930_8, 0.266751564401387_8, 0.161360240845588_8,   &
       0.06423404047594_8, 0.01590689966410_8, 0.00240199722109_8,                           &
       0.00022016340923_8, 0.00001224984692_8, 0.00000041090325_8 /)

  ! Get process information.
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  assert(size(region%grids) == 1)

  nDimensions = region%grids(1)%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nSpecies = region%solverOptions%nSpecies
  assert(nSpecies >= 0)
  if (nspecies.gt.2) then
     write(message, '(A)') "WARNING, max of 2 species (for now)"
     call gracefulExit(region%comm, message)
  end if

  ! Get species indices.
  if (allocated(region%solverOptions%speciesName)) then
     do k = 1, nSpecies + 1
        select case (trim(region%solverOptions%speciesName(k)))
        case ('H2', 'HYDROGEN')
           H2 = k
        case ('O2', 'OXYGEN')
           O2 = k
        case ('N2', 'NITROGEN')
           N2 = k
        case default
           write(message, '(3A)') "Unknown species: ",                                       &
                trim(region%solverOptions%speciesName(k)), "!"
           call gracefulExit(region%comm, message)
        end select
     end do
  else
     H2 = 1
     O2 = 2
     N2 = nSpecies + 1
  end if

  ! Get molecular weights.
  if (region%solverOptions%equationOfState == IDEAL_GAS_MIXTURE) then
     allocate(Wi(nSpecies+1))
     Wi = region%solverOptions%molecularWeightInverse
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
  ratioOfSpecificHeats = region%solverOptions%ratioOfSpecificHeats

  do k = iStart(3), iEnd(3)
     do j = iStart(2), iEnd(2)
        do i = iStart(1), iEnd(1)

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

           ! Offset by the virtual origin of the Blasius profile.
           xx = x + x0

           ! Return the first derivative of the Blasius function.
           delta = sqrt(xx / Re_c / crossflowVelocity)
           eta = (y - y0) / delta !... subtract off y0 for smoother transient.
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
                   1.0_WP/6.0_wp*f2l*(real(kk+1, wp) - eta)**3 +                             &
                   1.0_WP/6.0_wp*f2r*(eta-real(kk, wp))**3 +                                 &
                   (f0l-1.0_WP/6.0_WP*f2l)*(real(kk+1, wp) - eta) +                          &
                   (f0r-1.0_WP/6.0_WP*f2r)*(eta-real(kk, wp))

              f0l = by1(kk)
              f0r = by1(kk+1)
              f2l = -0.5_wp * by0(kk)*by2(kk)
              f2r = -0.5_wp * by0(kk+1)*by2(kk+1)
              blasius1 = &
                   1.0_wp/6.0_wp*f2l*(real(kk+1, wp)-eta)**3 +                               &
                   1.0_wp/6.0_wp*f2r*(eta-real(kk, wp))**3 +                                 &
                   (f0l-1.0_wp/6.0_WP*f2l)*(real(kk+1, wp) - eta) +                          &
                   (f0r-1.0_wp/6.0_WP*f2r)*(eta-real(kk, wp))

           end if
    
           ! Update the velocity.
           velocity(1) = crossflowVelocity * blasius1
           velocity(2) = 0.5_WP * sqrt(crossflowVelocity / xx / Re_c) *                      &
                (eta * blasius1 - blasius0)

           ! Make sure wall velocities are zero.
           if (j == 1) velocity = 0.0_wp

           ! Jet conditions.
           r = sqrt((x - jetPosition)**2 + z**2)
           insideJet = .false.
           select case(trim(jetShape))
           case('ROUND')
              if (r <= 1.0_wp * jetDiameter) insideJet = .true.
           case('SQUARE')
              if (i >= iJet1 .and. i <= iJet2 .and. k >= kJet1 .and. k <= kJet2)             &
                   insideJet = .true.
           case ('NONE')
              ! Nothing to do.
           case default
              write (*,'(A)') ''
              write(message, '(A)') 'Unknown jet shape' // trim(jetShape)
              call gracefulExit(region%comm, message)
           end select

           !if (insideJet) then
           yDecay = max(0.0_wp,                                                              &
                0.5_wp * (1.0_wp + tanh(10.0_wp * (0.5_wp - y / jetDiameter))))

           ! Fuel stream.
           sig = 6.0_wp
           fuel = Yf0 * yDecay
           if (trim(jetShape) == 'ROUND') then
              fuel = fuel * 0.5_wp * (                                                       &
                   tanh(sig * (r + 0.5_wp * jetDiameter) / jetDiameter) -                    &
                   tanh(sig * (r - 0.5_wp * jetDiameter) / jetDiameter))
           else if (trim(jetShape) == 'SQUARE') then
              fuel = fuel * 0.5_wp * (                                                       &
                   tanh(sig * (x - jetPosition + 0.5_wp * jetDiameter) / jetDiameter) -      &
                   tanh(sig * (x - jetPosition - 0.5_wp * jetDiameter) / jetDiameter))  *    &
                   0.5_wp * (                                                                &
                   tanh(sig * (z + 0.5_wp * jetDiameter) / jetDiameter) -                    &
                   tanh(sig * (z - 0.5_wp * jetDiameter) / jetDiameter))
           end if
           oxidizer = (1.0_wp - fuel) * Yo0

           ! Jet velocity.
           sig = 6.0_wp
           velocity(2) = jetVelocity * yDecay
           select case(trim(velocityProfile))
           case('TANH')

              if (trim(jetShape) == 'ROUND') then
                 velocity(2) = velocity(2) * 0.5_wp * (                                      &
                      tanh(sig * (r + 0.5_wp * jetDiameter) / jetDiameter) -                 &
                      tanh(sig * (r - 0.5_wp * jetDiameter) / jetDiameter))
              else if (trim(jetShape) == 'SQUARE') then
                 velocity(2) = velocity(2) * 0.5_wp * (                                      &
                      tanh(sig * (x - jetPosition + 0.5_wp * jetDiameter) / jetDiameter) -   &
                      tanh(sig * (x - jetPosition - 0.5_wp * jetDiameter) / jetDiameter))  * &
                      0.5_wp * (                                                             &
                      tanh(sig * (z + 0.5_wp * jetDiameter) / jetDiameter) -                 &
                      tanh(sig * (z - 0.5_wp * jetDiameter) / jetDiameter))
              end if

           case ('POISEUILLE')

              if (trim(jetShape) == 'ROUND') then
                 velocity(2) = max(0.0_wp,                                                   &
                      velocity(2) * 2.0_wp * (1.0_wp - (r / (0.5_wp * jetDiameter))**2))
              else if (trim(jetShape) == 'SQUARE') then
                 velocity(2) = 2.0_wp * velocity(2) * max(0.0_wp, (1.0_wp -                  &
                      ((x - jetPosition) / (0.5_wp * jetDiameter))**2)) *                    &
                      max(0.0_wp, (1.0_wp - (z / (0.5_wp * jetDiameter))**2))
              end if

           case ('SELF_SIMILAR')

              eta = (x - jetPosition)
              a = (sqrt(2.0_wp) - 1.0_wp) / (0.5_wp * jetDiameter)**2
              velocity(2) = velocity(2) * (1.0_wp + a * eta**2) ** (-2.0_wp)
              velocity(1) = 0.5_wp * jetVelocity * (eta - a * eta**3) /                      &
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
           select case (region%solverOptions%equationOfState)
           case(IDEAL_GAS)
              density = ratioOfSpecificHeats * pressure /                                    &
                   (temperature * (ratioOfSpecificHeats - 1.0_wp))
           case (IDEAL_GAS_MIXTURE)
              density = ratioOfSpecificHeats * pressure /                                    &
                   ( temperature * (ratioOfSpecificHeats - 1.0_wp) *                         &
                   (fuel * (Wi(H2) - Wi(N2)) + oxidizer * (Wi(O2) - Wi(N2)) + Wi(N2)) )
           end select

           ! State variables.
           region%states(1)%conservedVariables(grid_index(i,j,k), 1) = density
           region%states(1)%conservedVariables(grid_index(i,j,k), 2:nDimensions+1) =         &
                density * velocity(1:nDimensions)
           region%states(1)%conservedVariables(grid_index(i,j,k),nDimensions+2) =            &
                pressure / (ratioOfSpecificHeats - 1.0_wp) + 0.5_wp * density *              &
                sum(velocity ** 2)
           if (nSpecies.gt.0) region%states(1)%conservedVariables(grid_index(i,j,k),         &
                nDimensions+2+H2) = density * fuel
           if (nSpecies.gt.1) region%states(1)%conservedVariables(grid_index(i,j,k),         &
                nDimensions+2+O2) = density * oxidizer

        end do !... do i = iStart(1), iEnd(1)
     end do !... do j = iStart(2), iEnd(2)
  end do !... do k = iStart(3), iEnd(3)

  ! Cleanup.
  SAFE_DEALLOCATE(Wi)

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

end subroutine jetCrossFlowQ

subroutine jetCrossFlowTargetMollifier(region)

  ! <<< Private members >>>
  use JetCrossflow_mod

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

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, procRank, ierror
  real(wp) :: x
  real(wp), parameter :: s = 10.0_wp, r = 0.2_wp, eps = 1.0E-6_wp

  ! Initialize mollifier extents.
  targetStart = 1; targetEnd = -1

  ! Return for now.
  return

  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

end subroutine jetCrossFlowTargetMollifier

subroutine jetCrossFlowControlMollifier(region)

  ! <<< Private members >>>
  use JetCrossflow_mod

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

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, procRank, ierror
  real(wp) :: x
  real(wp), parameter :: s = 10.0_wp, r = 0.2_wp, eps = 1.0E-6_wp

  ! Initialize mollifier extents.
  controlStart = 1; controlEnd = -1

  ! Return for now.
  return

  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

end subroutine jetCrossFlowControlMollifier

subroutine jetCrossFlowBC

  ! <<< Private members >>>
  use JetCrossflow_mod, only : interiorStart, interiorEnd, targetStart, targetEnd,           &
       controlStart, controlEnd

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
  imax   (bc) =  interiorStart(1)
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
  imin   (bc) =  interiorEnd(1)
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
  jmin   (bc) =  interiorEnd(2)
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

end subroutine jetCrossFlowBC
