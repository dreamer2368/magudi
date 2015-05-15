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
  integer :: ierr, ierror
  character(len = STRING_LENGTH) :: grid_name,fname
  character(LEN=80) :: prefix, directory

  ! Magudi stuff
  type(t_Region) :: region
  integer, allocatable :: globalGridSizes(:,:)
  Logical :: success

  ! EnSight stuff
  integer :: num, numFiles
  integer :: i, j, k, l, n, nx, ny, nz, npart, reclength
  integer :: startIter, stopIter, skipIter, iter, var, nvar, ibuffer
  real(KIND=4), Dimension(:,:,:), Allocatable :: x,y,z,rbuffer
  real(KIND=4) :: xmin, xmax, ymin, ymax, zmin, zmax
  real(KIND=8), Dimension(:), Allocatable :: time
  character(LEN=2)  :: ib
  character(LEN=80), dimension(:), Allocatable :: names
  character(LEN=80) :: binary_form
  character(LEN=80) :: file_description1,file_description2
  character(LEN=80) :: node_id,element_id
  character(LEN=80) :: part,description_part,cblock,extents,cbuffer

  ! Initialize MPI.
  call MPI_Init(ierror)

  ! Read information from standard input
  print*,'==========================================='
  print*,'| Magudi - data to ENSIGHT GOLD converter |'
  print*,'==========================================='
  print*
  write (*,"(a15)",advance="no")  " output prefix "
  read "(a)", prefix
  write (*,"(a14)",advance="no")  " start iter "
  read "(i6)", startIter
  write (*,"(a13)",advance="no")  " stop iter "
  read "(i6)", stopIter
  write (*,"(a13)",advance="no")  " skip iter "
  read "(i6)", skipIter
  write (*,"(a14)",advance="no")  " use iblank "
  read "(a)", ib
  directory='ensight-3D'

  ! Check for errors.
  If ((ib .NE. 'y') .AND. (ib .NE. 'n')) Then
     Write (*,'(A)') 'Invalid value for IBLANK: use either "y" (yes) or "n" (no)'
     Stop
  End If

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

  ! Get number of variables.
  write(fname,'(2A,I8.8,A)') trim(prefix),'-', startIter, '.q'

  ! Load the solution file.
  call region%loadData(QOI_FORWARD_STATE, fname)
  nvar = size(region%states(1)%conservedVariables(1,:))
  print *, 'Number of variables:',nvar
  print *

  ! Assign variable names.
  allocate(names(nvar))
  select case (nvar)
  case(4)
     names(1) = 'RHO'
     names(2) = 'U'
     names(3) = 'V'
     names(4) = 'E'
  case (5)
     names(1) = 'RHO'
     names(2) = 'U'
     names(3) = 'V'
     names(4) = 'W'
     names(5) = 'E'
  end select

  ! Create the EnSight directory
  !call system("mkdir -p "//trim(directory))

  ! ======================= !
  ! Write the geometry file !
  ! ======================= !

  ! Get mesh extents (single precision)
  nx = region%grids(1)%globalSize(1)
  ny = region%grids(1)%globalSize(2)
  nz = region%grids(1)%globalSize(3)

  ! Store the coordinates in 3D
  allocate(X(nx,ny,nz))
  allocate(Y(nx,ny,nz))
  allocate(Z(nx,ny,nz))
  do k=1,nz
     do j=1,ny
        do i=1,nx
           x(i,j,k) = real(region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),1),4)
           y(i,j,k) = real(region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),2),4)
           if (nz.gt.1) then
              z(i,j,k) = region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),3)
           else
              z(i,j,k) = 0.0_4
           end if
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
  write(description_part,'(A,I1)') 'Grid', ng
  if (ib .eq. 'y') then
     cblock        ='block iblanked'
  else
     cblock        ='block'
  end If
  extents          ='extents'

  ! Size of file.
  if (ib .EQ. 'y') then
     reclength=80*8+4*(4+4*NX*NY*NZ)
  else
     reclength=80*9+4*(4+3*NX*NY*NZ)
  end if

  ! Write the file.
  write(fname,'(2A)') trim(directory), '/geometry'
  print *, 'Writing: ',fname
  !fname = trim(directory) // '/geometry'
  open (unit=10, file=trim(fname), form='unformatted', access="direct", recl=reclength)
  if (ib .EQ. 'y') then
     write(unit=10,rec=1) binary_form, &
          file_description1, &
          file_description2, &
          node_id, &
          element_id, &
          part,npart, &
          description_part, &
          cblock, &
          NX,NY,NZ, &
          (((real(X(I,J,K),4), I=1,NX), J=1,NY), K=1,NZ), &
          (((real(Y(I,J,K),4), I=1,NX), J=1,NY), K=1,NZ), &
          (((real(Z(I,J,K),4), I=1,NX), J=1,NY), K=1,NZ)!, &
          !(((IBLANK1(I,J,K), I=1,NX), J=1,NY), K=1,NZ)
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
          (((real(X(I,J,K),4), I=1,NX), J=1,NY), K=1,NZ), &
          (((real(Y(I,J,K),4), I=1,NX), J=1,NY), K=1,NZ), &
          (((real(Z(I,J,K),4), I=1,NX), J=1,NY), K=1,NZ)
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
    
    write(fname,'(2A,I8.8,A)') trim(prefix),'-', iter, '.q'
     write (*,'(A)') fname
     call region%loadData(QOI_FORWARD_STATE, fname)

     ! Binary file length
     reclength=80*3+4*(1+NX*NY*NZ)

     do var=1, nvar

        ! Convert the data to single precision
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 ! Store solution in buffer
                 rbuffer(i,j,k) = real(region%states(1)%conservedVariables(i+nx*(j-1+ny*(k-1)),var),4)

                 ! Divide out rho
                 if (var > 1) then
                    rbuffer(i,j,k) = rbuffer(i,j,k)/real(region%states(1)%conservedVariables(i+nx*(j-1+ny*(k-1)),1),4)
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

     ! ... counter
     num = num + 1
  end do
  

  ! =================== !
  ! Write the case file !
  ! =================== !
  write(fname,'(2A)') trim(directory), '/magudi.case'
  open(10,file=trim(fname),form="formatted",iostat=ierr,status="REPLACE")

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
  do var=1,nvar
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

end program data2ensight
