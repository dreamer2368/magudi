#include "config.h"

program shear_layer

  use MPI
  use, intrinsic :: iso_fortran_env

  use Grid_enum
  use State_enum

  use Region_mod, only : t_Region

  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush, gracefulExit
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
  filename = PROJECT_NAME // ".inp"
  call parseInputFile(filename)

  ! Verify that the grid file is in valid PLOT3D format and fetch the grid dimensions:
  ! `globalGridSizes(i,j)` is the number of grid points on grid `j` along dimension `i`.
  if (command_argument_count() == 1) then
     call get_command_argument(1, filename)
  else
     call getRequiredOption("grid_file", filename)
  end if
  i = len_trim(filename)
  if (filename(i-3:i) .ne. ".xyz") then
     filename = trim(filename) // ".xyz"
  end if

  ! Setup the region and load the grid file.
  call shearLayerGrid(region)

  ! Write out some useful information.
  call region%reportGridDiagnostics()

  ! Save grid.
  call region%saveData(QOI_GRID, filename)

  ! Generate the initial condition and target state.
  do i = 1, size(region%grids)
     call shearLayerInitialCondition(region%states(i), region%grids(i))
  end do

  ! Save initial condition.
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

  subroutine shearLayerGrid(region)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env, only : output_unit

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Enumerations >>>
    use Grid_enum, only : NONE, PLANE, OVERLAP

    ! <<< Internal modules >>>
    use InputHelper, only : getOption, getRequiredOption
    use ErrorHandler, only : gracefulExit, writeAndFlush

    implicit none

    ! <<< Arguments >>>
    type(t_Region) :: region

    ! Local variables
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, nx, ny, nz, gridIndex, gridIndex2, periodicityType(3)
    integer, dimension(3,1) :: globalGridSize
    real(wp) :: Lx, Ly, Lz, dx, dy, dz
    character(len = STRING_LENGTH) :: key, val, message
    real(wp) :: sigma, b, c, minMeshsize, maxMeshsize, y1,  y2
    real(wp), allocatable, dimension(:) :: s, g
    logical :: stretchGrid
    integer :: ierror

    interface
       subroutine mapping_function(s, b, c, sigma, g)
         real(SCALAR_KIND), intent(in) :: s(:), b, c, sigma
         real(SCALAR_KIND), intent(out) :: g(size(s))
       end subroutine mapping_function
    end interface

    globalGridSize(1,1) = getOption("shear_layer/grid_size/x",100)
    globalGridSize(2,1) = getOption("shear_layer/grid_size/y",100)
    globalGridSize(3,1) = getOption("shear_layer/grid_size/z",100)

    call region%setup(MPI_COMM_WORLD,globalGridSize)

    ! Simplify
    nx = globalGridSize(1,1)
    ny = globalGridSize(2,1)
    nz = globalGridSize(3,1)

    ! Read in the grid size
    Lx = getOption("shear_layer/domain_size/x",0.0_wp)
    Ly = getOption("shear_layer/domain_size/y",0.0_wp)
    Lz = getOption("shear_layer/domain_size/z",0.0_wp)

    ! Read in periodicity option
    do i = 1, 3
       write(key, '(A,I3.3,A,I1.1,A)') "grid", region%grids(1)%index, "/dir", i, "/"
       val = getOption(trim(key) // "periodicity_type", "")
       periodicityType(i) = NONE
       if (trim(val) == "PLANE") then
          periodicityType(i) = PLANE
       else if (trim(val) == "OVERLAP") then
          periodicityType(i) = OVERLAP
       end if
    end do !... i =  1, size(globalSize)

    ! Compute the grid spacing
    if (periodicityType(1) .eq. PLANE) then
       dx = Lx / real(nx, wp)
    else
       dx = Lx / real(nx-1, wp)
    end if
    if (periodicityType(2) .eq. PLANE) then
       dy = Ly / real(ny, wp)
    else
       dy = Ly / real(ny-1, wp)
    end if
    if (periodicityType(3) .eq. PLANE) then
       dz = Lz / real(nz, wp)
    else
       dz = Lz / real(nz-1, wp)
    end if

    ! Should we stretch the mesh?
    stretchGrid = getOption("shear_layer/stretch_grid",.false.)

    ! Generate the grid
    do k = region%grids(1)%offset(3) + 1, region%grids(1)%offset(3) + region%grids(1)%localSize(3)
       do j = region%grids(1)%offset(2) + 1, region%grids(1)%offset(2) + region%grids(1)%localSize(2)
          do i = region%grids(1)%offset(1) + 1, region%grids(1)%offset(1) + region%grids(1)%localSize(1)

             gridIndex = i - region%grids(1)%offset(1) + region%grids(1)%localSize(1) *         &
                   (j - region%grids(1)%offset(2) - 1 + region%grids(1)%localSize(2) *          &
                   (k - region%grids(1)%offset(3) - 1))

             ! Create X
             if (nx .gt. 1) then
                if (periodicityType(1) .eq. PLANE) then
                   region%grids(1)%coordinates(gridIndex, 1) = (Lx - dx) *  real(i - 1, wp) /  &
                        real(nx - 1, wp)
                else
                   region%grids(1)%coordinates(gridIndex, 1) = Lx * real(i - 1, wp) /          &
                        real(nx - 1, wp)
                end if
             end if

             ! Create Y
             if (ny .gt. 1 .and. .not. stretchGrid) then
                region%grids(1)%coordinates(gridIndex, 2) = Ly * real(j - 1, wp) /             &
                     real(ny - 1, wp) - 0.5_wp * Ly
             end if

             ! Create Z
             if (nz .gt. 1) then
                if (periodicityType(3) .eq. PLANE) then
                   region%grids(1)%coordinates(gridIndex, 3) = (Lz - dz) *  real(k - 1, wp) /  &
                        real(nz - 1, wp) - 0.5_wp * Lz
                else
                   region%grids(1)%coordinates(gridIndex, 3) = Lz * real(k - 1, wp) /          &
                        real(nz - 1, wp) - 0.5_wp * Lz
                end if
             end if

          end do
       end do
    end do

    ! Grid stretching
    if (stretchGrid) then
       ! Grid stretching parameters
       sigma = 0.21_wp
       b = 12.0_wp
       c = 0.6_wp

       ! Create uniform spacing
       allocate(s(ny))
       do j = 1, ny
          s(j) = real(j - 1, wp) / real(ny - 1, wp)
       end do

       ! Compute mapping g(s)
       allocate(g(ny))
       call mapping_function(s, b, c, sigma, g)

       ! Find min/max spacing
       minMeshsize =  huge(1.0_wp)
       maxMeshsize = -huge(1.0_wp)

       do k = region%grids(1)%offset(3) + 1, region%grids(1)%offset(3) + region%grids(1)%localSize(3)
          do j = region%grids(1)%offset(2) + 1, region%grids(1)%offset(2) + region%grids(1)%localSize(2)
             do i = region%grids(1)%offset(1) + 1, region%grids(1)%offset(1) + region%grids(1)%localSize(1)
                gridIndex = i - region%grids(1)%offset(1) + region%grids(1)%localSize(1) *          &
                       (j - region%grids(1)%offset(2) - 1 + region%grids(1)%localSize(2) *          &
                       (k - region%grids(1)%offset(3) - 1))
                ! Create y
                region%grids(1)%coordinates(gridIndex, 2) = 0.5_wp * Ly * (1.0_wp + g(j)) - 0.5_wp * Ly

                ! Find min/max spacing
                if (j .gt. 2) then
                   gridIndex2 = gridIndex - region%grids(1)%localSize(1)
                   y1 = region%grids(1)%coordinates(gridIndex, 2)
                   y2 = region%grids(1)%coordinates(gridIndex2, 2)
                   minMeshsize = min(minMeshsize, abs(y2 - y1))
                   maxMeshsize = max(maxMeshsize, abs(y2 - y1))
                end if
             end do
          end do
       end do
       call MPI_Allreduce(MPI_IN_PLACE, minMeshsize, 1, REAL_TYPE_MPI,                            &
                          MPI_MIN, region%grids(1)%comm, ierror)
       call MPI_Allreduce(MPI_IN_PLACE, maxMeshsize, 1, REAL_TYPE_MPI,                            &
                          MPI_MAX, region%grids(1)%comm, ierror)
       if (region%commGridMasters /= MPI_COMM_NULL) then
          call MPI_Allreduce(MPI_IN_PLACE, minMeshsize, 1, REAL_TYPE_MPI,                         &
          MPI_MIN, region%commGridMasters, ierror)
          call MPI_Allreduce(MPI_IN_PLACE, maxMeshsize, 1, REAL_TYPE_MPI,                         &
          MPI_MAX, region%commGridMasters, ierror)
       end if

       do i = 1, size(region%grids)
         call MPI_Bcast(minMeshsize, 1, REAL_TYPE_MPI, 0, region%grids(i)%comm, ierror)
         call MPI_Bcast(maxMeshsize, 1, REAL_TYPE_MPI, 0, region%grids(i)%comm, ierror)
       end do
       write(message,'(A,2(1X,'//SCALAR_FORMAT//'))') 'min/max y-spacing:', minMeshsize, maxMeshsize
       call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

       deallocate(s)
       deallocate(g)
    end if

    ! Update the grids by computing the Jacobian, metrics, and norm.
    do i = 1, size(region%grids)
       call region%grids(i)%update()
    end do

  end subroutine shearLayerGrid

  subroutine shearLayerInitialCondition(state, grid)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env, only : output_unit

    ! <<< Derived types >>>
    use State_mod, only : t_State
    use Grid_mod, only : t_Grid

    ! <<< Internal modules >>>
    use InputHelper, only : getOption, getRequiredOption
    use ErrorHandler, only : writeAndFlush, gracefulExit
    use RandomNumber, only : initializeRandomNumberGenerator, random

    implicit none

    ! <<< Arguments >>>
    type(t_Grid) :: grid
    type(t_State) :: state

    ! Local variables
    integer :: i, j, k, ierror
    integer, parameter :: nkx = 10, nkz = 5
    integer, parameter :: wp = SCALAR_KIND
    real(SCALAR_KIND), parameter :: pi = 4.0_wp * atan(1.0_wp)
    integer:: nDimensions, nUnknowns, midplaneCount, zint, xint, numBlocks = 5
    real(wp) :: density, temperature, pressure, velocity, velocityDifference,                  &
         lowerTemperature, upperTemperature, T0, P0, x, y, z, Lpx, Lpz,                        &
         scalar, initialThickness, rand, ratioOfSpecificHeats,                                 &
         shear_layer_initial_rms_fraction, midplaneRMS, minDensity, maxDensity, Yv
    real(wp), dimension(:,:), allocatable :: gradA, phaseJitter, velocityFluctuations
    real(wp), dimension(:), allocatable :: A
    real(wp), dimension(nkx, -nkz:nkz) :: phases
    character(len=STRING_LENGTH) :: message

    interface
      subroutine broadband_real(x, y, z, Lpx, Lpz, lengthScale, phases,num_kx, num_kz, val)
        real(SCALAR_KIND), intent(in) :: x, y, z, Lpx, Lpz, lengthScale
        real(SCALAR_KIND), intent(out) :: val
        integer,intent(in):: num_kx, num_kz
        real(SCALAR_KIND), dimension(num_kx,-num_kz:num_kz),intent(in) :: phases
      end subroutine
    end interface

    ! Set the number of conserved variables and allocate the array
    nDimensions = grid%nDimensions
    nUnknowns = nDimensions + 2

    ! Initialize random number generator
    call initializeRandomNumberGenerator()

    ! Zero-out the conserved variables
    do i = 1, size(region%grids)
      region%states(i)%conservedVariables = 0.0_wp
    end do

    ! Reference temperature and pressure
    ratioOfSpecificHeats = getOption("ratio_of_specific_heats", 1.4_wp)
    T0 = 1.0_wp / (ratioOfSpecificHeats - 1.0_wp)
    P0 = 1.0_wp / ratioOfSpecificHeats

    ! Read in the shear layer parameters
    call getRequiredOption('shear_layer/thickness', initialThickness)
    call getRequiredOption('shear_layer/velocity_difference', velocityDifference)
    lowerTemperature = getOption('shear_layer/lower_temperature', T0)
    upperTemperature = getOption('shear_layer/upper_temperature', T0)
    pressure = getOption('shear_layer/pressure', P0)

    ! Mean profiles
    minDensity = huge(1.0_wp); maxDensity = -huge(1.0_wp)
    do i = 1, grid%nGridPoints

      ! Vertical coordinate
      y = grid%coordinates(i, 2)

      ! Velocity profile
      velocity = 0.5_wp * velocityDifference * tanh(0.5_wp * y / initialThickness)

      ! Temperature profile
      temperature = lowerTemperature + 0.5_wp * (upperTemperature - lowerTemperature) *       &
          (1.0_wp + tanh(0.5_wp * y / initialThickness))

      ! Assign the density
      ! Crocco-Busemann (Sandham 1990 & Vaghefi 2014)
      density = 1.0_wp + 0.5_wp * (ratioOfSpecificHeats - 1.0_wp) *                        &
               (0.25_wp * velocityDifference**2 - velocity**2)
      density = 1.0_wp / density

      ! Keep track of min/max density
      minDensity = min(minDensity, density)
      maxDensity = max(maxDensity, density)

      ! Assign the density
      state%conservedVariables(i,1) = density

      ! Assign the mean velocity
      state%conservedVariables(i,2) = velocity
    end do
    call MPI_Allreduce(MPI_IN_PLACE, minDensity, 1, REAL_TYPE_MPI,                            &
                       MPI_MIN, grid%comm, ierror)
    call MPI_Allreduce(MPI_IN_PLACE, maxDensity, 1, REAL_TYPE_MPI,                            &
                       MPI_MAX, grid%comm, ierror)
    if (region%commGridMasters /= MPI_COMM_NULL) then
       call MPI_Allreduce(MPI_IN_PLACE, minDensity, 1, REAL_TYPE_MPI,                         &
       MPI_MIN, region%commGridMasters, ierror)
       call MPI_Allreduce(MPI_IN_PLACE, maxDensity, 1, REAL_TYPE_MPI,                         &
       MPI_MAX, region%commGridMasters, ierror)
    end if

    call MPI_Bcast(minDensity, 1, REAL_TYPE_MPI, 0, grid%comm, ierror)
    call MPI_Bcast(maxDensity, 1, REAL_TYPE_MPI, 0, grid%comm, ierror)
    write(message,'(A,2(1X,'//SCALAR_FORMAT//'))') 'min/max density:', minDensity, maxDensity
    call writeAndFlush(MPI_COMM_WORLD, output_unit, message)

    ! Storage for vector potential
    allocate(A(grid%nGridPoints))
    allocate(gradA(grid%nGridPoints, nDimensions))

    ! Vector potential
    call getRequiredOption("shear_layer/domain_size/x",Lpx)
    Lpz = 0.0_wp
    if (nDimensions .eq. 3) call getRequiredOption("shear_layer/domain_size/z",Lpz)

    phases = 0.0_wp
    do j = 1, nkx
      do i = -nkz, nkz
         call random_number(rand)
         rand = 2.0_wp * rand - 1.0_wp
         phases(j,i) = rand * pi
      end do
    end do

    allocate(phaseJitter(0:numBlocks, 0:numBlocks))
    phaseJitter = 0.0_wp
    do i = 0, numBlocks
      do j = 0, numBlocks
         call random_number(rand)
         phaseJitter(i,j) = rand
      end do
    end do

    z = 0.0_wp
    zint = 1
    do i = 1, grid%nGridPoints
      ! Get the local coordinates
      x = grid%coordinates(i, 1)
      y = grid%coordinates(i, 2)
      z = 0.0_wp
      if (nDimensions .eq. 3) z = grid%coordinates(i, 3) + 0.5_wp * Lpz

      xint = int(x / (Lpx / real(numBlocks, wp)))
      if (nDimensions .eq. 3) zint = int(z / (Lpz / real(numBlocks, wp)))

      ! Broadbanded scalar with Phase Jittering similar to Kim thesis and Lui thesis,
      ! except theirs were phase jitted in time.
      if (phaseJitter(xint, zint) .lt. 0.33_wp) then
         call broadband_real(x, y, z, Lpx, Lpz, initialThickness, phases + pi, nkx, nkz,   &
              scalar)
      else if (phaseJitter(xint, zint) .lt. 0.667_wp .and.                                 &
           phaseJitter(xint, zint) .ge. 0.33_wp) then
         call broadband_real(x, y, z, Lpx, Lpz, initialThickness, phases - pi, nkx, nkz,   &
              scalar)
      else
         call broadband_real(x, y, z, Lpx, Lpz, initialThickness, phases, nkx, nkz, scalar)
      end if
      A(i) = scalar
    end do

    deallocate(phaseJitter) ! ... not needed anymore

    ! Vector cross product to get a solenoidal VELOCITY field
    call grid%computeGradient(A, gradA)

    allocate(velocityFluctuations(grid%nGridPoints, 3))

    if (nDimensions .eq. 3) then
      velocityFluctuations(:, 1) = (gradA(:, 2) - gradA(:, 3))
      velocityFluctuations(:, 2) =-(gradA(:, 1) - gradA(:, 3))
      velocityFluctuations(:, 3) = (gradA(:, 1) - gradA(:, 2))
    else if (nDimensions .eq. 2) then
      velocityFluctuations(:, 1) = gradA(:, 2)
      velocityFluctuations(:, 2) =-gradA(:, 1)
    end if

    ! Rescale the velocity perturbations to have rms of x% of delta U
    midplaneRMS = 0.0_wp; midplaneCount = 0
    do i = 1, grid%nGridPoints
      y = grid%coordinates(i, 2)
      if (y .ge. -initialThickness .and. y .le. initialThickness) then
         midplaneRMS = midplaneRMS + sum(velocityFluctuations(i,:)**2)
         midplaneCount = midplaneCount + 1
      end if
    end do
    call MPI_Allreduce(MPI_IN_PLACE, midplaneRMS, 1, REAL_TYPE_MPI,                            &
                       MPI_SUM, grid%comm, ierror)
    call MPI_Allreduce(MPI_IN_PLACE, midplaneCount, 1, MPI_INTEGER,                            &
                       MPI_SUM, grid%comm, ierror)
    if (region%commGridMasters /= MPI_COMM_NULL) then
       call MPI_Allreduce(MPI_IN_PLACE, midplaneRMS, 1, REAL_TYPE_MPI,                         &
       MPI_SUM, region%commGridMasters, ierror)
       call MPI_Allreduce(MPI_IN_PLACE, midplaneCount, 1, MPI_INTEGER,                         &
       MPI_SUM, region%commGridMasters, ierror)
    end if

    call MPI_Bcast(midplaneRMS, 1, REAL_TYPE_MPI, 0, grid%comm, ierror)
    call MPI_Bcast(midplaneCount, 1, MPI_INTEGER, 0, grid%comm, ierror)
    midplaneRMS = midplaneRMS / real(midplaneCount, wp)
    midplaneRMS = midplaneRMS**0.5

    ! Rescale the velocity
    call getRequiredOption('shear_layer/rms_fraction', shear_layer_initial_rms_fraction)

    velocityFluctuations = velocityFluctuations * shear_layer_initial_rms_fraction /        &
        midplaneRMS * velocityDifference

    ! Mollify the perturbations after making it solenoidal
    do i = 1, grid%nGridPoints
      y = grid%coordinates(i, 2)
      scalar = 0.5_wp * (tanh(5.0_wp * (y + 1.0_wp)) - (tanh(5.0_wp * (y - 1.0_wp))))
      velocityFluctuations(i,:) = velocityFluctuations(i,:) * scalar
    end do

    ! Add perturbations to the mean specified from above.
    state%conservedVariables(:, 2:nDimensions+1) = state%conservedVariables(:, 2:nDimensions+1) +       &
        velocityFluctuations(:, 1:nDimensions)

    ! Clean up
    deallocate(A)
    deallocate(gradA)
    deallocate(velocityFluctuations)

    do i = 1, grid%nGridPoints

       ! Set the momentum
       state%conservedVariables(i,2:nDimensions+1) = state%conservedVariables(i,1) *                       &
            state%conservedVariables(i,2:nDimensions+1)

       ! Assign the energy
       state%conservedVariables(i, nDimensions+2) = pressure / (ratioOfSpecificHeats - 1.0_wp) +     &
            0.5_wp * state%conservedVariables(i,2)**2 / state%conservedVariables(i,1)
       if (nDimensions .gt. 1) state%conservedVariables(i, nDimensions+2) =                          &
            state%conservedVariables(i, nDimensions+2) + 0.5_wp * state%conservedVariables(i,3)**2 /       &
            state%conservedVariables(i,1)
       if (nDimensions .eq. 3) state%conservedVariables(i, nDimensions+2) =                          &
            state%conservedVariables(i, nDimensions+2) + 0.5_wp * state%conservedVariables(i,4)**2 /       &
            state%conservedVariables(i,1)

    end do

  end subroutine shearLayerInitialCondition

end program shear_layer

! Mapping functional for grid stretching
! --------------------------------------
subroutine mapping_function(s, b, c, sigma, g)

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  real(SCALAR_KIND), parameter :: pi = 4.0_wp * atan(1.0_wp)
  real(SCALAR_KIND), intent(in) :: s(:), b, c, sigma
  real(SCALAR_KIND), intent(out) :: g(size(s))

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

! Routine for generating realistic velocity perturbations
! -------------------------------------------------------
subroutine broadband_real(x, y, z, Lpx, Lpz, lengthScale, phases,num_kx, num_kz, val)

  implicit none

  ! Arguments
  real(SCALAR_KIND), intent(in) :: x, y, z, Lpx, Lpz, lengthScale
  real(SCALAR_KIND), intent(out) :: val
  integer,intent(in):: num_kx, num_kz
  real(SCALAR_KIND), dimension(num_kx,-num_kz:num_kz),intent(in) :: phases

  ! Local variables
  integer, parameter :: wp = SCALAR_KIND
  real(SCALAR_KIND), parameter :: pi = 4.0_wp * atan(1.0_wp)
  real(wp) :: kxo, kzo, kyo
  integer:: i, k

  kxo = 0.0_wp; kzo = 0.0_wp
  if (Lpx .gt. 0.0_wp) kxo = 2.0_wp / (Lpx)
  if (Lpz .gt. 0.0_wp) kzo = 2.0_wp / (Lpz)
  kyo = 1.0_wp / (10.0_wp * lengthScale)
  val = 0.0_wp
  do i = 1, num_kx
     do k = -num_kz, num_kz
        val = val + cos(2.0_wp * pi * real(i, wp) * kxo * x +                              &
             2.0_wp * pi * real(k, wp) * kzo * z +                                         &
             2.0_wp * pi * kyo * y + phases(i,k))
     end do
  end do

  return
end subroutine broadband_real
