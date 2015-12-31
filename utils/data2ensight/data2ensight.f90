#include "config.h"

program data2ensight

  use MPI
  use, intrinsic :: iso_fortran_env

  use Grid_enum
  use State_enum
  use Region_enum, only : FORWARD, ADJOINT

  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver

  use InputHelper, only : getInputName, parseInputFile, getOption, getRequiredOption
  use ErrorHandler, only : writeAndFlush, gracefulExit
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use CNSHelper, only : computeDependentVariables

  implicit none

  ! <<< Local variables >>>
  type(t_Region) :: region
  type(t_Solver) :: solver
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, procRank, numProcs, num, numFiles, useAdjoint, ierror
  integer, allocatable :: globalGridSizes(:,:)
  logical :: success
  character(len = STRING_LENGTH) :: prefix, filename
  integer :: startIter, stopIter, skipIter, iter

  ! Initialize MPI.
  call MPI_Init(ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)

  ! Parse the input file.
  call getInputName(filename)
  call parseInputFile(filename)

  ! Read information from standard input.
  if (procRank == 0) then
     print *
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
     print *
  end if

  ! Broadcast the information.
  call MPI_Bcast(startIter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call MPI_Bcast(stopIter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call MPI_Bcast(skipIter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
  call MPI_Bcast(useAdjoint, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

  ! Check for errors.
  if (useAdjoint /= 0 .and. useAdjoint /=1) then
     print *, 'ERROR: read adjoint /= 0 or 1!'
     stop
  end if

  ! Get the grid filename.
  call getRequiredOption("grid_file", filename)

  ! Verify that the grid file is in valid PLOT3D format and fetch the grid dimensions:
  ! `globalGridSizes(i,j)` is the number of grid points on grid `j` along dimension `i`.
  call plot3dDetectFormat(MPI_COMM_WORLD, filename, success,                                 &
       globalGridSizes = globalGridSizes)
  if (.not. success) call gracefulExit(MPI_COMM_WORLD, plot3dErrorMessage)

  ! Setup the region and load the grid file.
  call region%setup(MPI_COMM_WORLD, globalGridSizes)
  print *
  call region%loadData(QOI_GRID, filename)

  ! Compute normalized metrics, norm matrix and Jacobian.
  do i = 1, size(region%grids)
     call region%grids(i)%update()
  end do
  call MPI_Barrier(MPI_COMM_WORLD, ierror)

  ! Setup EnSight directory and write geometry.
  allocate(solver%ensight(size(region%grids)))
  do i = 1, size(region%grids)
     call solver%ensight(i)%setup(region%grids(i), region%states(i)%time)
  end do

  ! Get number of files to read in.
  numFiles = 0
  do iter = startIter, stopIter, skipIter
     numFiles = numFiles + 1
  end do

  ! Get file prefix.
  call getRequiredOption("output_prefix", prefix)

  ! Loop through files and write.
  num = 1
  do iter = startIter, stopIter, skipIter
    
     ! Read in the solution file.
     write(filename,'(2A,I8.8,A)') trim(prefix),'-', iter, '.q'
     print *
     call region%loadData(QOI_FORWARD_STATE, filename)

     ! Update the state and write the EnSight files.
     do i = 1, size(region%states)
        call computeDependentVariables(region%grids(i)%nDimensions,                          &
             region%solverOptions%nSpecies, region%states(i)%conservedVariables,             &
             region%solverOptions%equationOfState, region%solverOptions%ratioOfSpecificHeats,&
             region%solverOptions%molecularWeightInverse,                                    &
             region%states(i)%specificVolume(:,1), region%states(i)%velocity,                &
             region%states(i)%pressure(:,1), region%states(i)%temperature(:,1),              &
             region%states(i)%massFraction)
        call solver%ensight(i)%output(region%states(i), region%grids(i),FORWARD,             &
             region%states(1)%time, region%solverOptions%nSpecies)
     end do

     if (useAdjoint == 1) then
        write(filename,'(2A,I8.8,A)') trim(prefix),'-', iter, '.adjoint.q'
        print *
        call region%loadData(QOI_ADJOINT_STATE, filename)
        do i = 1, size(region%states)
           call solver%ensight(i)%output(region%states(i), region%grids(i), ADJOINT,         &
                region%states(1)%time, region%solverOptions%nSpecies)
        end do
     end if

  end do

  ! Cleanup.
  call region%cleanup()

  ! Finalize MPI.
  call MPI_Finalize(ierror)

end program data2ensight
