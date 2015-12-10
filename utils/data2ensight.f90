#include "config.h"

program data2ensight
! ============================================================= !
!                                                               !
! Converts PLOT3D data file to EnSight format                   !
!                                                               !
! ============================================================= !
!                                                               !
! Written by Jesse Capecelatro (jcaps@illinois.edu)             !
! 9 August 2014                                                 !
!                                                               !
! ============================================================= !

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env

 ! <<< Derived types >>>
  use Region_mod, only : t_Region

  ! <<< Enumerations >>>
  use Grid_enum
  use State_enum

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables
  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush, gracefulExit
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage

  implicit none

  ! I/O
  integer :: useAdjoint, ierror
  character(len = STRING_LENGTH) :: grid_name, fname, prefix
  character(LEN=80) :: directory

  ! Magudi stuff
  type(t_Region) :: region
  integer, allocatable :: globalGridSizes(:,:)
  Logical :: success

  ! EnSight stuff
  integer :: num, numFiles
  integer :: nx, ny, nz, nDimensions, nSpecies, npart, reclength, ii
  integer :: startIter, stopIter, skipIter, iter, var, nvar
  real(KIND=4), dimension(:,:,:), allocatable :: x,y,z,iblank,rbuffer
  real(KIND=4) :: xmin, xmax, ymin, ymax, zmin, zmax
  real(KIND=8), dimension(:), allocatable :: time, temperature
  character(LEN=80), dimension(:), allocatable :: names
  character(LEN=80) :: binary_form
  character(LEN=80) :: file_description1,file_description2
  character(LEN=80) :: node_id,element_id
  character(LEN=80) :: part,description_part,cblock,extents,cbuffer
  logical :: useIblank

  ! Local variables
  integer :: i, j, k
  real(KIND=8) :: gamma, heatRelease, zelDovich, referenceTemperature, flameTemperature,     &
       activationTemperature, density, Yf, Yo, Da

  ! Initialize MPI
  call MPI_Init(ierror)

  ! Parse the input file
  call parseInputFile("magudi.inp")

  ! Get file prefix
  call getRequiredOption("output_prefix", prefix)

  ! Get combustion parameters
  call getRequiredOption("ratio_of_specific_heats", gamma)
  call getRequiredOption("heat_release", heatRelease)
  call getRequiredOption("Zel_Dovich", zelDovich)
  call getRequiredOption("Damkohler_number", Da)
  referenceTemperature = 1.0_8 / (gamma - 1.0_8)
  flameTemperature = referenceTemperature / (1.0_8 - heatRelease)
  activationTemperature = zelDovich / heatRelease * flameTemperature

  ! Read information from standard input
  print*,'==========================================='
  print*,'| Magudi - data to ENSIGHT GOLD converter |'
  print*,'==========================================='
  print*
  write (*,"(a14)",advance="no")  " start iter "
  read "(i6)", startIter
  write (*,"(a13)",advance="no")  " stop iter "
  read "(i6)", stopIter
  write (*,"(a13)",advance="no")  " skip iter "
  read "(i6)", skipIter
  write (*,"(a15)",advance="no")  " read adjoint? "
  read "(i6)", useAdjoint
  directory='ensight-3D'

  ! Check for errors.
  if (useAdjoint /= 0 .and. useAdjoint /=1) then
     print *, 'WARNING: read adjoint /= 0 or 1!'
     stop
  end if

  ! Set the grid name.
  grid_name = trim(prefix)//'.xyz'

  ! Verify that the grid file is in valid PLOT3D format and fetch the grid dimensions:
  ! `globalGridSizes(i,j)` is the number of grid points on grid `j` along dimension `i`.
  call plot3dDetectFormat(MPI_COMM_WORLD, grid_name,                                         &
       success, globalGridSizes = globalGridSizes)
  if (.not. success) call gracefulExit(MPI_COMM_WORLD, plot3dErrorMessage)

  ! Setup the region and load the grid file.
  call region%setup(MPI_COMM_WORLD, globalGridSizes)
  call region%loadData(QOI_GRID, grid_name)

  ! Compute normalized metrics, norm matrix and Jacobian.
  do i = 1, size(region%grids)
     call region%grids(i)%update()
  end do
  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  ! Write out some useful information.
  call region%reportGridDiagnostics()

  ! Get number of dimensions.
  nDimensions = region%grids(1)%nDimensions

  ! Get number of species.
  nSpecies = region%solverOptions%nSpecies

  ! Get number of variables.
  write(fname,'(2A,I8.8,A)') trim(prefix),'-', startIter, '.q'
  call region%loadData(QOI_FORWARD_STATE, fname)
  nvar = size(region%states(1)%conservedVariables(1,:))

  ! Sanity check.
  if (nvar.ne.region%solverOptions%nUnknowns) then
     print *, 'Something is wrong! (nvar /= nUnknowns)'
     stop
  end if

  ! Check if iblank is used
  useIblank = .false.
  if (size(region%grids(1)%iblank) > 0) useIblank = .true.

  ! Include adjoint variables.
  if (useAdjoint == 1) nvar = nvar * 2
  
  print *, 'Number of variables:',nvar
  print *

  ! Get number of files and time sereies.
  numFiles = 0
  do iter = startIter, stopIter, skipIter
     numFiles = numFiles + 1
  end do
  allocate(time(numFiles))
  time = 0.0_8
  numFiles = 0
  do iter = startIter, stopIter, skipIter
     numFiles = numFiles + 1
     time(numFiles) = real(numFiles-1,8)
  end do

  ! Assign variable names.
  allocate(names(nvar+3))
  select case (nDimensions)
  case (2)
     names(1) = 'RHO'
     names(2) = 'U'
     names(3) = 'V'
     names(4) = 'E'
  case (3)
     names(1) = 'RHO'
     names(2) = 'U'
     names(3) = 'V'
     names(4) = 'W'
     names(5) = 'E'
  end select

  select case (nSpecies)
  case (1)
     names(nDimensions+3) = 'H2'
  case (2)
     names(nDimensions+3) = 'H2'
     names(nDimensions+4) = 'O2'
  end select

  ! Adjoint variable names.
  if (useAdjoint == 1) then
     do i = 1, nvar/2
        write(names(i+nvar/2),'(2A)') trim(names(i)), '_adjoint'
     end do
  end if

  ! Output additional data
  names(nvar+1) = 'Temperature'
  names(nvar+2) = 'Reaction'
  if (useIblank) names(nvar+3) = 'iblank'

  ! Create the EnSight directory
  !call system("mkdir -p "//trim(directory))

  ! ======================= !
  ! Write the geometry file !
  ! ======================= !

  ! Get mesh extents (single precision)
  nx = region%grids(1)%globalSize(1)
  ny = region%grids(1)%globalSize(2)
  nz = region%grids(1)%globalSize(3)

  ! Store the coordinates and iblank in 3D
  allocate(x(nx,ny,nz))
  allocate(y(nx,ny,nz))
  allocate(z(nx,ny,nz))
  if (useIblank) allocate(iblank(nx,ny,nz))
  do k=1,nz
     do j=1,ny
        do i=1,nx
           x(i,j,k) = real(region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),1),4)
           y(i,j,k) = real(region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),2),4)
           if (nz.gt.1) then
              z(i,j,k) = real(region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),3),4)
           else
              z(i,j,k) = 0.0_4
           end if
           if (useIblank) iblank(i,j,k) = real(region%grids(1)%iblank(i+nx*(j-1+ny*(k-1))),4)
        end do
     end do
  end do

  xmin = minval(x)
  xmax = maxval(x)
  ymin = minval(y)
  ymax = maxval(y)
  zmin = minval(z)
  zmax = maxval(z)

  ! Write EnSight geometry
  binary_form      ='C Binary'
  file_description1='Ensight Gold Geometry File'
  file_description2='WriteRectEnsightGeo Routine'
  node_id          ='node id is off'
  element_id       ='element id off'
  part             ='part'
  npart            =1
  write(description_part,'(A)') 'Grid'
  cblock        ='block'
  extents          ='extents'

  ! Size of file.
  reclength=80*9+4*(4+3*NX*NY*NZ)

  ! Write the file.
  write(fname,'(2A)') trim(directory), '/geometry'
  print *, 'Writing: ',fname
  open (unit=10, file=trim(fname), form='unformatted', access="direct", recl=reclength)
  write(unit=10,rec=1) binary_form, &
       file_description1, &
       file_description2, &
       node_id, &
       element_id, &
       part,npart, &
       description_part, &
       cblock, &
       NX,NY,NZ, &
       (((x(i,j,k), i=1,nx), j=1,ny), k=1,nz), &
       (((y(i,j,k), i=1,nx), j=1,ny), k=1,nz), &
       (((z(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  
  ! Close the file.
  close (10)


  ! ===================================== !
  ! Loop through and write the data files !
  ! ===================================== !

  ! Allocate single-precision array
  allocate(rbuffer(nx,ny,nz))

  ! Allocate temperature array
  allocate(temperature(size(region%states(1)%conservedVariables, 1)))

  ! Loop through files and write
  num = 1
  do iter = startIter, stopIter, skipIter
     print *, 'Writing timestep',iter
    
     write(fname,'(2A,I8.8,A)') trim(prefix),'-', iter, '.q'
     call region%loadData(QOI_FORWARD_STATE, fname)

     if (useAdjoint == 1) then
        write(fname,'(2A,I8.8,A)') trim(prefix),'-', iter, '.adjoint.q'
        call region%loadData(QOI_ADJOINT_STATE, fname)
     end if

    ! Binary file length
    reclength=80*3+4*(1+nx*ny*nz)

     do var=1, nvar

        ! Convert the data to single precision
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 ii = i+nx*(j-1+ny*(k-1))
                 ! Store solution in buffer
                 if (useAdjoint==1 .and. var>nvar/2) then
                    ! Store adjoint variable.
                    rbuffer(i,j,k) = real(region%states(1)%adjointVariables(ii,var-nvar/2),4)
                 else
                    ! Store state variable.
                    rbuffer(i,j,k) = real(region%states(1)%conservedVariables(ii,var),4)

                    ! Divide out rho
                    if (var > 1) then
                       rbuffer(i,j,k) = rbuffer(i,j,k) /                                     &
                            real(region%states(1)%conservedVariables(ii,1),4)
                    end if
                 end if
              end do
           end do
        end do

        ! Write Ensight scalar file
        cbuffer=trim(names(var))
        write(fname,'(4A,I6.6)') trim(directory),'/',trim(names(var)),'.',num
        open (unit=10, file=trim(fname), form='unformatted', access="direct", recl=reclength)
        write(unit=10, rec=1) cbuffer,part,npart,cblock,(((rbuffer(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
        close(10)

     end do ! var

     ! Get the temperature.
     call computeDependentVariables(nDimensions, nSpecies,                                   &
             region%states(1)%conservedVariables,                                            &
             region%solverOptions%equationOfState,                                           &
             region%solverOptions%ratioOfSpecificHeats,                                      &
             region%solverOptions%molecularWeightInverse,                                    &
             temperature = temperature)
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ii = i+nx*(j-1+ny*(k-1))
              rbuffer(i,j,k) = real(temperature(ii), 4)
           end do
        end do
     end do

     ! Dimenionalize
     rbuffer = rbuffer * (real(gamma, 4) - 1.0_4) * 293.15_4

     ! Write temperature to Ensight file
     cbuffer=trim(names(nvar+1))
     write(fname,'(4A,I6.6)') trim(directory),'/',trim(names(nvar+1)),'.',num
     open (unit=10, file=trim(fname), form='unformatted', access="direct", recl=reclength)
     write(unit=10, rec=1) cbuffer,part,npart,cblock,(((rbuffer(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
     close(10)

     ! Reaction rate
     rbuffer = 0.0_4
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ii = i+nx*(j-1+ny*(k-1))
              ! Get local variables
              density = region%states(1)%conservedVariables(ii, 1)
              Yf = region%states(1)%conservedVariables(ii, nDimensions+3) / density
              Yo = region%states(1)%conservedVariables(ii, nDimensions+4) / density

              ! Reaction rate
              rbuffer(i,j,k) = real(Da * density * Yf * Yo *                                 &
                   exp(- activationTemperature / temperature(ii)), 4)
           end do
        end do
     end do

     ! Write temperature to Ensight file
     cbuffer=trim(names(nvar+2))
     write(fname,'(4A,I6.6)') trim(directory),'/',trim(names(nvar+2)),'.',num
     open (unit=10, file=trim(fname), form='unformatted', access="direct", recl=reclength)
     write(unit=10, rec=1) cbuffer,part,npart,cblock,(((rbuffer(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
     close(10)

     ! Write the iblank array to EnSight file
     if (useIblank) then
        cbuffer=trim(names(nvar+3))
        write(fname,'(4A,I6.6)') trim(directory),'/',trim(names(nvar+3)),'.',num
        open (unit=10, file=trim(fname), form='unformatted', access="direct", recl=reclength)
        write(unit=10, rec=1) cbuffer,part,npart,cblock,(((iblank(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
        close(10)
     end if

     ! ... counter
     num = num + 1
  end do

  ! =================== !
  ! Write the case file !
  ! =================== !
  write(fname,'(2A)') trim(directory), '/magudi.case'
  open(10,file=trim(fname),form="formatted",iostat=ierror,status="REPLACE")

  ! Write the case
  cbuffer='FORMAT'
  write(10,'(a80)') cbuffer
  cbuffer='type: ensight gold'
  write(10,'(a80)') cbuffer
  cbuffer='GEOMETRY'
  write(10,'(a80)') cbuffer
  write(cbuffer,'(A)') 'model: geometry'
  write(10,'(a80)') cbuffer

  cbuffer='VARIABLE'
  write(10,'(a)') cbuffer
  do var = 1, nvar + 2
     write(cbuffer,'(5A)') 'scalar per node: 1 ',  trim(names(var)), ' ',                    &
          trim(names(var)), '.******'
     write(10,'(a80)') cbuffer
  end do
  if (useIblank) then
     write(cbuffer,'(5A)') 'scalar per node: 1 ',  trim(names(nvar+3)), ' ',                 &
          trim(names(nvar+3)), '.******'
     write(10,'(a80)') cbuffer
  end if
  
  cbuffer='TIME'
  write(10,'(a80)') cbuffer
  cbuffer='time set: 1'
  write(10,'(a80)') cbuffer
  write(10,'(a16,i12)') 'number of steps:',numFiles
  cbuffer='filename start number: 1'
  write(10,'(a80)') cbuffer
  cbuffer='filename increment: 1'
  write(10,'(a80)') cbuffer
  write(10,'(a12,10000000(3(ES12.5),/))') 'time values:',time

  ! Close the file
  close(10)

  ! Clean up & leave
  deallocate(x, y, z, rbuffer, time)
  if (useIblank) deallocate(iblank)

end program data2ensight
