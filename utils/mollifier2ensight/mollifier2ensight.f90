#include "config.h"

program mollifier2ensight
! ============================================================= !
!                                                               !
! Converts PLOT3D mollifier file to EnSight format              !
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
  character(len = STRING_LENGTH) :: grid_name, fname, prefix
  character(LEN=80) :: directory

  ! Magudi stuff
  type(t_Region) :: region
  integer, allocatable :: globalGridSizes(:,:)
  Logical :: success

  ! EnSight stuff
  integer :: nx, ny, nz, ndim, npart, reclength, ii
  real(KIND=4), Dimension(:,:,:), Allocatable :: x,y,z,iblank,rbuffer
  real(KIND=4) :: xmin, xmax, ymin, ymax, zmin, zmax
  character(LEN=80) :: binary_form
  character(LEN=80) :: file_description1,file_description2
  character(LEN=80) :: node_id,element_id
  character(LEN=80) :: part,description_part,cblock,extents,cbuffer
  logical :: use_iblank, writeTarget, writeControl

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
  directory='ensight-3D'

  ! Decide which mollifiers to write
  writeTarget = .false.
  fname = getOption("target_mollifier_file", "")
  if (trim(fname).ne."") writeTarget = .true.
  writeControl = .false.
  fname = getOption("control_mollifier_file", "")
  if (trim(fname).ne."") writeControl = .true.

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

  ! Load the mollifier files.
  if (writeTarget) then
     fname = trim(prefix)//'.target_mollifier.f'
     call region%loadData(QOI_TARGET_MOLLIFIER, fname)
  end if

  if (writeControl) then
     fname = trim(prefix)//'.control_mollifier.f'
     call region%loadData(QOI_CONTROL_MOLLIFIER, fname)
  end if

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
  print *
  print *, 'Writing: ',fname
  print *
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


  ! ==================== !
  ! Write the mollifiers !
  ! ==================== !

  ! Allocate single-precision array
  allocate(rbuffer(nx,ny,nz))

  ! Binary file length
  reclength=80*3+4*(1+nx*ny*nz)

  if (writeTarget) then
     ! Convert the data to single precision
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! Store solution in buffer
              ii = i+nx*(j-1+ny*(k-1))
              rbuffer(i,j,k) = real(region%grids(1)%targetMollifier(ii,1),4)
           end do
         end do
     end do

     ! Write Ensight scalar file
     cbuffer='TARGET_MOLLIFIER'
     write(fname,'(4A,I6.6)') trim(directory),'/',trim(cbuffer),'.',1
     open (unit=10, file=trim(fname), form='unformatted', access="direct", recl=reclength)
     write(unit=10, rec=1) cbuffer,part,npart,cblock,(((rbuffer(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
     close(10)
  end if

  if (writeControl) then
     ! Convert the data to single precision
     do k=1,nz
        do j=1,ny
           do i=1,nx
              ! Store solution in buffer
              ii = i+nx*(j-1+ny*(k-1))
              rbuffer(i,j,k) = real(region%grids(1)%controlMollifier(ii,1),4)
           end do
        end do
     end do

     ! Write Ensight scalar file
     cbuffer='CONTROL_MOLLIFIER'
     write(fname,'(4A,I6.6)') trim(directory),'/',trim(cbuffer),'.',1
     open (unit=10, file=trim(fname), form='unformatted', access="direct", recl=reclength)
     write(unit=10, rec=1) cbuffer,part,npart,cblock,(((rbuffer(i,j,k),i=1,NX),j=1,NY),k=1,NZ)
     close(10)
  end if

  ! =================== !
  ! Write the case file !
  ! =================== !
  write(fname,'(2A)') trim(directory), '/mollifier.case'
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
  write(10,'(a80)') cbuffer

  if (writeTarget) then
     write(cbuffer,'(5A)') 'scalar per node: 1 ',  'TARGET_MOLLIFIER', ' ', 'TARGET_MOLLIFIER', '.******'
     write(10,'(a80)') cbuffer
  end if
  if (writeControl) then  
     write(cbuffer,'(5A)') 'scalar per node: 1 ',  'CONTROL_MOLLIFIER', ' ', 'CONTROL_MOLLIFIER', '.******'
     write(10,'(a80)') cbuffer
 end if
  
  cbuffer='TIME'
  write(10,'(a80)') cbuffer
  cbuffer='time set: 1'
  write(10,'(a80)') cbuffer
  write(10,'(a16,i12)') 'number of steps:',1
  cbuffer='filename start number: 1'
  write(10,'(a80)') cbuffer
  cbuffer='filename increment: 1'
  write(10,'(a80)') cbuffer
  write(10,'(a12,10000000(3(ES12.5),/))') 'time values:',0.0_8

  ! Close the file
  close(10)

  ! Clean up & leave
  deallocate(X, Y, Z, rbuffer)
  if (use_iblank) deallocate(iblank)

end program mollifier2ensight
