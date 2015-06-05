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
  use MPI
  use, intrinsic :: iso_fortran_env

  use Grid_enum
  use State_enum

  use Region_mod, only : t_Region

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
  integer :: nx, ny, nz, nspec, ndim, npart, reclength, ii
  integer :: startIter, stopIter, skipIter, iter, var, nvar
  real(KIND=4), Dimension(:,:,:), Allocatable :: x,y,z,iblank,rbuffer
  real(KIND=4) :: xmin, xmax, ymin, ymax, zmin, zmax
  real(KIND=8), Dimension(:), Allocatable :: time
  character(LEN=80), dimension(:), Allocatable :: names
  character(LEN=80) :: binary_form
  character(LEN=80) :: file_description1,file_description2
  character(LEN=80) :: node_id,element_id
  character(LEN=80) :: part,description_part,cblock,extents,cbuffer
  logical :: use_iblank

  ! Local variables
  integer :: i, j, k

  ! Initialize MPI.
  call MPI_Init(ierror)

  ! Parse the input file.
  call parseInputFile("magudi.inp")

  ! Get file prefix
  call getRequiredOption("output_prefix", prefix)

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
  call plot3dDetectFormat(MPI_COMM_WORLD, grid_name,                                          &
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
  ndim = size(region%grids(1)%coordinates(1,:))

  ! Get number of variables.
  write(fname,'(2A,I8.8,A)') trim(prefix),'-', startIter, '.q'
  call region%loadData(QOI_FORWARD_STATE, fname)
  nvar = size(region%states(1)%conservedVariables(1,:))

  ! Number of species.
  nspec = nvar - (ndim+2)

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
  allocate(names(nvar+1))
  select case (ndim)
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

  select case (nspec)
  case (1)
     names(ndim+3) = 'H2'
  case (2)
     names(ndim+3) = 'H2'
     names(ndim+4) = 'O2'
  end select

  ! Adjoint variable names.
  if (useAdjoint == 1) then
     do i = 1, nvar/2
        write(names(i+nvar/2),'(2A)') trim(names(i)), '_adjoint'
     end do
  end if

  ! Include temperature
  names(nvar+1) = 'Temperature'

  ! Check if iblank is used
  use_iblank = .false.
  if (size(region%grids(1)%iblank) > 0) use_iblank = .true.

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
  if (use_iblank) allocate(iblank(nx,ny,nz))
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
           if (use_iblank) iblank(i,j,k) = real(region%grids(1)%iblank(i+nx*(j-1+ny*(k-1))),4)
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
  if (use_iblank) then
     cblock        ='block iblanked'
  else
     cblock        ='block'
  end If
  extents          ='extents'

  ! Size of file.
  if (use_iblank) then
     reclength=80*8+4*(4+4*NX*NY*NZ)
  else
     reclength=80*9+4*(4+3*NX*NY*NZ)
  end if

  ! Write the file.
  write(fname,'(2A)') trim(directory), '/geometry'
  print *, 'Writing: ',fname
  open (unit=10, file=trim(fname), form='unformatted', access="direct", recl=reclength)
  if (use_iblank) then
     write(unit=10,rec=1) binary_form, &
          file_description1, &
          file_description2, &
          node_id, &
          element_id, &
          part,npart, &
          description_part, &
          cblock, &
          nx,ny,nz, &
          (((x(i,j,k), i=1,nx), j=1,ny), k=1,nz), &
          (((y(i,j,k), i=1,nx), j=1,ny), k=1,nz), &
          (((z(i,j,k), i=1,nx), j=1,ny), k=1,nz), &
          (((iblank(i,j,k), i=1,nx), j=1,ny), k=1,nz)
  else
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
  end if
  
  ! Close the file.
  close (10)


  ! ===================================== !
  ! Loop through and write the data files !
  ! ===================================== !

  ! Allocate single-precision array
  allocate(rbuffer(nx,ny,nz))

  ! Loop through files and write
  num = 1
  do iter = startIter, stopIter, skipIter
     print *, 'Writing timestep',iter
    
     write(fname,'(2A,I8.8,A)') trim(prefix),'-', startIter, '.q'
     call region%loadData(QOI_FORWARD_STATE, fname)

     if (useAdjoint == 1) then
        write(fname,'(2A,I8.8,A)') trim(prefix),'-', startIter, '.adjoint.q'
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
                 if (var > nvar/2) then
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

     ! Temperature
     rbuffer = 0.0_4
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ii = i+nx*(j-1+ny*(k-1))
              ! Store temperature in buffer
              select case (ndim)
              case (2)
                 rbuffer(i,j,k) = real(1.4_8 * (region%states(1)%conservedVariables(ii,ndim+2) - &
                      0.5_8 * (region%states(1)%conservedVariables(ii,2)**2 +                    &
                      region%states(1)%conservedVariables(ii,3)**2) /                            &
                      region%states(1)%conservedVariables(ii,1)) /                               &
                      region%states(1)%conservedVariables(ii,1),4)

              case (3)
                 rbuffer(i,j,k) = real(1.4_8 * (region%states(1)%conservedVariables(ii,ndim+2) - &
                      0.5_8 * (region%states(1)%conservedVariables(ii,2)**2 +                    &
                      region%states(1)%conservedVariables(ii,3)**2 +                             &
                      region%states(1)%conservedVariables(ii,4)**2) /                            &
                      region%states(1)%conservedVariables(ii,1)) /                               &
                      region%states(1)%conservedVariables(ii,1),4)

              end select
           end do
        end do
     end do

     ! Dimenionalize
     rbuffer = rbuffer * 0.4_4 * 293.15_4

     ! Write temperature to Ensight file
     cbuffer=trim(names(nvar+1))
     write(fname,'(4A,I6.6)') trim(directory),'/',trim(names(nvar+1)),'.',num
     open (unit=10, file=trim(fname), form='unformatted', access="direct", recl=reclength)
     write(unit=10, rec=1) cbuffer,part,npart,cblock,(((rbuffer(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
     close(10)

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
  do var=1,nvar+1
     write(cbuffer,'(5A)') 'scalar per node: 1 ',  trim(names(var)), ' ', trim(names(var)), '.******'
     write(10,'(a80)') cbuffer
  end do
  
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
  deallocate(X, Y, Z, rbuffer, time)
  if (use_iblank) deallocate(iblank)

end program data2ensight
