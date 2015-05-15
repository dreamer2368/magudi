#include "config.h"

Program data2ensight
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

  !  Reals
  Real(KIND=4), Dimension(:,:,:), Allocatable :: X,Y,Z,rbuffer
  Real(KIND=4) :: XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX
  Real(KIND=8), Dimension(:), Allocatable :: time

  !  Integers
  Integer :: ng, num, numFiles
  Integer :: ngrid, I, J, K, l, N, NX, NY, NZ, npart, reclength
  Integer :: startIter, stopIter, skipIter, iter, var, nvar, ierr, ierror, ibuffer
  integer, allocatable :: globalGridSizes(:,:)

  !  Characters
  character(LEN=80), dimension(:), Allocatable :: names
  Character(LEN=2)  :: prec, gf, vf, ib, ib2
  Character(LEN=80) :: prefix,directory
  character(LEN=80) :: binary_form
  character(LEN=80) :: file_description1,file_description2
  character(LEN=80) :: node_id,element_id
  character(LEN=80) :: part,description_part,cblock,extents,cbuffer
  character(len = STRING_LENGTH) :: grid_name,fname

  ! Logicals
  Logical :: success

  ! Type
  type(t_Region) :: region

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

  ! Check for errors
  If ((ib .NE. 'y') .AND. (ib .NE. 'n')) Then
     Write (*,'(A)') 'Invalid value for IBLANK: use either "y" (yes) or "n" (no)'
     Stop
  End If

  ! Set the grid name
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

  ! Initialize
  numFiles = 0
  Do iter = startIter, stopIter, skipIter
     numFiles = numFiles + 1
  End Do
  allocate(time(numFiles))
  time = 0.0_8
  numFiles = 0
  Do iter = startIter, stopIter, skipIter
     numFiles = numFiles + 1
     time(numFiles) = real(numFiles-1,8)
  End Do

  ! Find number of variables (should be 5+)
  Write(fname,'(2A,I8.8,A)') trim(prefix),'-', startIter, '.q'

  ! Load the solution file.
  call region%loadData(QOI_FORWARD_STATE, fname)
  nvar = size(region%states(1)%conservedVariables(1,:))
  print *, 'Number of variables:',nvar
  print *

  ! Assign variable names
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
  NX = region%grids(1)%globalSize(1)
  NY = region%grids(1)%globalSize(2)
  NZ = region%grids(1)%globalSize(3)

  ! Store the coordinates in 3D
  allocate(X(nx,ny,nz))
  allocate(Y(nx,ny,nz))
  allocate(Z(nx,ny,nz))
  print *, size(region%grids(1)%coordinates(1,:))
  do k=1,nz
     do j=1,ny
        do i=1,nx
           X(i,j,k) = real(region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),1),4)
           Y(i,j,k) = real(region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),2),4)
           if (nz.gt.1) then
              Z(i,j,k) = region%grids(1)%coordinates(i+nx*(j-1+ny*(k-1)),3)
           else
              Z(i,j,k) = 0.0_4
           end if
        end do
     end do
  end do

  XMIN = minval(X)
  XMAX = maxval(X)
  YMIN = minval(Y)
  YMAX = maxval(Y)
  ZMIN = minval(Z)
  ZMAX = maxval(Z)

  ! Write EnSight geometry
  binary_form      ='C Binary'
  file_description1='Ensight Gold Geometry File'
  file_description2='WriteRectEnsightGeo Routine'
  node_id          ='node id is off'
  element_id       ='element id off'
  part             ='part'
  npart            =1
  Write(description_part,'(A,I1)') 'Grid', ng
  If (ib .EQ. 'y') then
     cblock        ='block iblanked'
  else
     cblock        ='block'
  End If
  extents          ='extents'

  ! Size of file
  If (ib .EQ. 'y') Then
     reclength=80*8+4*(4+4*NX*NY*NZ)
  Else
     reclength=80*9+4*(4+3*NX*NY*NZ)
  End If

  ! Write file
  Write(fname,'(2A)') trim(directory), '/geometry'
  print *, 'Writing: ',fname
  !fname = trim(directory) // '/geometry'
  Open (unit=10, file=trim(fname), form='unformatted', access="direct", recl=reclength)
  If (ib .EQ. 'y') then
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
  Else
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
  End If
  
  ! Close the file
  Close (10)


  ! ===================================== !
  ! Loop through and write the data files !
  ! ===================================== !

  ! Allocate single-precision array
  Allocate(rbuffer(NX,NY,NZ))

  ! Loop through files and write
  num = 1
  Do iter = startIter, stopIter, skipIter
    print *, 'Writing timestep',iter
    
    Write(fname,'(2A,I8.8,A)') trim(prefix),'-', iter, '.q'
     Write (*,'(A)') fname
     call region%loadData(QOI_FORWARD_STATE, fname)

     ! Binary file length
     reclength=80*3+4*(1+NX*NY*NZ)

     Do var=1, nvar

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
        Write(fname,'(4A,I6.6)') trim(directory),'/',trim(names(var)),'.',num
        Open (unit=10, file=trim(fname), form='unformatted', access="direct", recl=reclength)
        write(unit=10, rec=1) cbuffer,part,npart,cblock,(((rbuffer(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
        Close(10)

     End Do

     ! ... counter
     num = num + 1
  End Do
  

  ! =================== !
  ! Write the case file !
  ! =================== !
  Write(fname,'(2A)') trim(directory), '/magudi.case'
  open(10,file=trim(fname),form="formatted",iostat=ierr,status="REPLACE")

  ! Write the case
  cbuffer='FORMAT'
  write(10,'(a80)') cbuffer
  cbuffer='type: ensight gold'
  write(10,'(a80)') cbuffer
  cbuffer='GEOMETRY'
  write(10,'(a80)') cbuffer
  Write(cbuffer,'(A)') 'model: geometry'
  write(10,'(a80)') cbuffer

  cbuffer='VARIABLE'
  write(10,'(a)') cbuffer
  Do var=1,nvar
     Write(cbuffer,'(5A)') 'scalar per node: 1 ',  trim(names(var)), ' ', trim(names(var)), '.******'
     write(10,'(a80)') cbuffer
  End Do
  
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
  Deallocate(X, Y, Z, rbuffer, time)

End Program data2ensight
