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

  !  Arrays
  Real(KIND=8), Dimension(:), Allocatable :: time
  Real(KIND=4), Dimension(:,:,:), Allocatable :: rbuffer
  Real(KIND=4) :: XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX

  !  Integers
  Integer :: ng, num, numFiles
  Integer :: ngrid, I, J, K, N, NX, NY, NZ, NDIM, npart, reclength
  Integer :: startIter, stopIter, skipIter, iter, var, nvar, ierr, ierror, ibuffer
  integer, allocatable :: globalGridSizes(:,:)

  !  Characters
  character(LEN=80), dimension(:), Allocatable :: names
  Character(LEN=2)  :: prec, gf, vf, ib, ib2
  Character(LEN=80) :: prefix,fname,directory
  character(LEN=80) :: binary_form
  character(LEN=80) :: file_description1,file_description2
  character(LEN=80) :: node_id,element_id
  character(LEN=80) :: part,description_part,cblock,extents,cbuffer
  character(len = STRING_LENGTH) :: grid_name

  ! Logicals
  Logical :: success

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
  print *, grid_name

  ! Verify that the grid file is in valid PLOT3D format and fetch the grid dimensions:
  ! `globalGridSizes(i,j)` is the number of grid points on grid `j` along dimension `i`.
  call plot3dDetectFormat(MPI_COMM_WORLD, grid_name,                                          &
       success, globalGridSizes = globalGridSizes)
  if (.not. success) call gracefulExit(MPI_COMM_WORLD, plot3dErrorMessage)

  ! Setup the region and load the grid file.
  call region%setup(MPI_COMM_WORLD, globalGridSizes)
  call region%loadData(QOI_GRID, trim(grid_name))
!!$
!!$  ! Compute normalized metrics, norm matrix and Jacobian.
!!$  do i = 1, size(region%grids)
!!$     call region%grids(i)%update()
!!$  end do
!!$  call MPI_Barrier(MPI_COMM_WORLD, ierror)
!!$
!!$  ! Write out some useful information.
!!$  call region%reportGridDiagnostics()

  ! Initialize
  numFiles = 0
  Do iter = startIter, stopIter, skipIter
     numFiles = numFiles + 1
  End Do

  ! Find number of variables (should be 5+)
  Write(fname,'(2A,I8.8,A)') trim(prefix),'-', startIter, '.q'

!!$  Call Read_Soln(NDIM,ngrid,ND1,Q1,TAU,prec,gf,vf,fname,0)
!!$  nvar=SIZE(Q1(1,1,1,1,:))
!!$  If (nvar < 5) then
!!$     Write (*,'(A)') 'Number of variables less than 5'
!!$     Stop
!!$  Else
!!$     Print *, 'Number of variables:',nvar
!!$     Write (*,'(A)') ''
!!$  End If
!!$
!!$  ! Assign variable names
!!$  allocate(names(nvar))
!!$  names(1) = 'RHO'
!!$  names(2) = 'U'
!!$  names(3) = 'V'
!!$  names(4) = 'W'
!!$  names(5) = 'E'
!!$
!!$  ! Loop over all files, find time array first
!!$  allocate(time(numFiles))
!!$  time = 0.0_8
!!$  num=0
!!$  Do iter = startIter, stopIter, skipIter
!!$    Write(fname,'(A,I8.8,A)') 'RocFlo-CM.', iter, '.q'
!!$
!!$    ! ... read solution header to find time
!!$    Call Read_Soln_Header(NDIM, ngrid, ND1, fname, tau)
!!$
!!$    ! ... save time
!!$    num = num + 1
!!$    time(num) = tau(4)
!!$
!!$  End Do
!!$
!!$  ! Create the EnSight directory
!!$  call system("mkdir -p "//trim(directory))
!!$
!!$
!!$  ! ======================= !
!!$  ! Write the geometry file !
!!$  ! ======================= !
!!$
!!$  !  Step one, read in the mesh1 grid file
!!$  Write (*,'(A)') ''
!!$  Call Read_Grid(NDIM,ngrid,ND1,X1,IBLANK1, &
!!$       prec,gf,vf,ib2,grid_name,1)
!!$  Write (*,'(A)') ''
!!$
!!$  ! Loop through number of grids
!!$  Do ng=1, ngrid
!!$     ! Get mesh extents (single precision)
!!$     NX = ND1(ng,1)
!!$     NY = ND1(ng,2)
!!$     NZ = ND1(ng,3)
!!$
!!$     XMIN = minval(X1(ng,1:NX,1:NY,1:NZ,1))
!!$     XMAX = maxval(X1(ng,1:NX,1:NY,1:NZ,1))
!!$     YMIN = minval(X1(ng,1:NX,1:NY,1:NZ,2))
!!$     YMAX = maxval(X1(ng,1:NX,1:NY,1:NZ,2))
!!$     ZMIN = minval(X1(ng,1:NX,1:NY,1:NZ,3))
!!$     ZMAX = maxval(X1(ng,1:NX,1:NY,1:NZ,3))
!!$
!!$     ! Write EnSight geometry
!!$     binary_form      ='C Binary'
!!$     file_description1='Ensight Gold Geometry File'
!!$     file_description2='WriteRectEnsightGeo Routine'
!!$     node_id          ='node id is off'
!!$     element_id       ='element id off'
!!$     part             ='part'
!!$     npart            =1
!!$     Write(description_part,'(A,I1)') 'Grid', ng
!!$     If (ib .EQ. 'y') then
!!$        cblock        ='block iblanked'
!!$     else
!!$        cblock        ='block'
!!$     End If
!!$     extents          ='extents'
!!$
!!$     ! Size of file
!!$     If (ib .EQ. 'y') Then
!!$        reclength=80*8+4*(4+4*NX*NY*NZ)
!!$     Else
!!$        reclength=80*9+4*(4+3*NX*NY*NZ)
!!$     End If
!!$
!!$     ! Write file
!!$     Write(fname,'(2A,I1)') trim(directory), '/geometry', ng
!!$     print *, 'Writing: ',fname
!!$     !fname = trim(directory) // '/geometry'
!!$     Open (unit=10, file=trim(fname), form='unformatted', access="direct", recl=reclength)
!!$     If (ib .EQ. 'y') then
!!$        write(unit=10,rec=1) binary_form, &
!!$             file_description1, &
!!$             file_description2, &
!!$             node_id, &
!!$             element_id, &
!!$             part,npart, &
!!$             description_part, &
!!$             cblock, &
!!$             NX,NY,NZ, &
!!$             (((sngl(X1(ng,I,J,K,1)), I=1,NX), J=1,NY), K=1,NZ), &
!!$             (((sngl(X1(ng,I,J,K,2)), I=1,NX), J=1,NY), K=1,NZ), &
!!$             (((sngl(X1(ng,I,J,K,3)), I=1,NX), J=1,NY), K=1,NZ), &
!!$             (((IBLANK1(ng,I,J,K), I=1,NX), J=1,NY), K=1,NZ)
!!$     Else
!!$        write(unit=10,rec=1) binary_form, &
!!$             file_description1, &
!!$             file_description2, &
!!$             node_id, &
!!$             element_id, &
!!$             part,npart, &
!!$             description_part, &
!!$             cblock, &
!!$             NX,NY,NZ, &
!!$             (((sngl(X1(ng,I,J,K,1)), I=1,NX), J=1,NY), K=1,NZ), &
!!$             (((sngl(X1(ng,I,J,K,2)), I=1,NX), J=1,NY), K=1,NZ), &
!!$             (((sngl(X1(ng,I,J,K,3)), I=1,NX), J=1,NY), K=1,NZ)
!!$     End If
!!$  
!!$     ! Close the file
!!$     Close (10)
!!$
!!$  End Do
!!$
!!$
!!$  ! ===================================== !
!!$  ! Loop through and write the data files !
!!$  ! ===================================== !
!!$
!!$  ! Allocate single-precision array
!!$  Allocate(rbuffer(1:maxval(ND1(:,1)),1:maxval(ND1(:,2)),1:maxval(ND1(:,3))))
!!$
!!$  ! Loop through files and write
!!$  num = 1
!!$  Do iter = startIter, stopIter, skipIter
!!$    print *, 'Writing timestep',iter
!!$    
!!$     Write(fname,'(A,I8.8,A)') 'RocFlo-CM.', iter, '.q'
!!$     Write (*,'(A)') fname
!!$     Call Read_Soln(NDIM,ngrid,ND1,Q1,TAU,prec,gf,vf,fname,0)
!!$
!!$     ! Loop through each grid
!!$     Do ng=1, ngrid
!!$
!!$        ! Get domain size
!!$        NX = ND1(ng,1)
!!$        NY = ND1(ng,2)
!!$        NZ = ND1(ng,3)
!!$
!!$        ! Binary file length
!!$        reclength=80*3+4*(1+NX*NY*NZ)
!!$
!!$        Do var=1, nvar
!!$
!!$           ! Convert the data to single precision
!!$           If (var <= 5) Then
!!$              rbuffer(1:NX,1:NY,1:NZ) = sngl(Q1(ng,1:NX,1:NY,1:NZ,var))
!!$           Else
!!$              rbuffer(1:NX,1:NY,1:NZ) = sngl(Q2(ng,1:NX,1:NY,1:NZ,var-5)/Q1(ng,1:NX,1:NY,1:NZ,1))
!!$           End If
!!$
!!$           ! Divide by rho for rhoU, rhoV, and rhoW
!!$           If (var>1 .and. var<=5) then
!!$              rbuffer(1:NX,1:NY,1:NZ)=sngl(Q1(ng,1:NX,1:NY,1:NZ,var)/Q1(ng,1:NX,1:NY,1:NZ,1))
!!$           End If
!!$
!!$           ! Write Ensight scalar file
!!$           cbuffer=trim(names(var))
!!$           Write(fname,'(3A,I1,A,I6.6)') trim(directory),'/',trim(names(var)),ng,'.',num
!!$           Open (unit=10, file=trim(fname), form='unformatted', access="direct", recl=reclength)
!!$           write(unit=10, rec=1) cbuffer,part,npart,cblock,(((rbuffer(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
!!$           Close(10)
!!$
!!$        End Do
!!$
!!$        ! Write IBLANK
!!$        rbuffer(1:NX,1:NY,1:NZ) = Real(IBLANK1(ng,1:NX,1:NY,1:NZ),4)
!!$        Write(fname,'(2A,I1,A,I6.6)') trim(directory),'/IB',ng,'.',num
!!$        Open (unit=10, file=trim(fname), form='unformatted', access="direct", recl=reclength)
!!$        write(unit=10, rec=1) cbuffer,part,npart,cblock,(((rbuffer(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
!!$        Close(10)
!!$     End Do
!!$
!!$     ! ... counter
!!$     num = num + 1
!!$  End Do
!!$  
!!$
!!$  ! =================== !
!!$  ! Write the case file !
!!$  ! =================== !
!!$  Do ng=1, ngrid
!!$     Write(fname,'(2A,I1,A)') trim(directory), '/plascomcm', ng, '.case'
!!$     open(10,file=trim(fname),form="formatted",iostat=ierr,status="REPLACE")
!!$
!!$     ! Write the case
!!$     cbuffer='FORMAT'
!!$     write(10,'(a80)') cbuffer
!!$     cbuffer='type: ensight gold'
!!$     write(10,'(a80)') cbuffer
!!$     cbuffer='GEOMETRY'
!!$     write(10,'(a80)') cbuffer
!!$     Write(cbuffer,'(A,I1)') 'model: geometry', ng
!!$     write(10,'(a80)') cbuffer
!!$
!!$     cbuffer='VARIABLE'
!!$     write(10,'(a)') cbuffer
!!$     Do var=1,nvar
!!$        Write(cbuffer,'(2A,I1,2A,I1,A)') 'scalar per node: 1 ',  trim(names(var)), ng, ' ', trim(names(var)), ng, '.******'
!!$        write(10,'(a80)') cbuffer
!!$     End Do
!!$
!!$     ! Write IB
!!$     Write(cbuffer,'(1A,I1,1A,I1,A)') 'scalar per node: 1 IB', ng, ' IB', ng, '.******'
!!$     write(10,'(a80)') cbuffer
!!$  
!!$     cbuffer='TIME'
!!$     write(10,'(a80)') cbuffer
!!$     cbuffer='time set: 1'
!!$     write(10,'(a80)') cbuffer
!!$     write(10,'(a16,x,i12)') 'number of steps:',numFiles
!!$     cbuffer='filename start number: 1'
!!$     write(10,'(a80)') cbuffer
!!$     cbuffer='filename increment: 1'
!!$     write(10,'(a80)') cbuffer
!!$     write(10,'(a12,x,10000000(3(ES12.5,x),/))') 'time values:',time
!!$
!!$     ! Close the file
!!$     close(10)
!!$
!!$  End Do
!!$
!!$  ! Clean up & leave
!!$  Deallocate(ND1, X1, Q1, IBLANK1, time)

End Program data2ensight
