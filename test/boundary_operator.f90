#include "config.h"

program boundary_operator

  use MPI
  use, intrinsic :: iso_fortran_env, only : real64

  use MPIHelper, only : pigeonhole
  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random
  use StencilOperator_mod, only : t_StencilOperator

  implicit none

  integer, parameter :: wp = SCALAR_KIND
  real(real64), parameter :: testDuration = 2.0_real64
  character(len = STRING_LENGTH), parameter :: discretizationTypes(4) =                      &
       (/ "SBP 1-2", "SBP 2-4", "SBP 3-6", "SBP 4-8" /)
  real(wp), parameter :: tolerance = epsilon(0.0_wp)
  integer :: i, j, k, n, localIndex, direction, nDimensions, nComponents, faceOrientation,   &
       localSize(3), offset(3), globalSize(3), cartesianCommunicator,                        &
       numProcs, procRank, ierror
  integer, dimension(:), allocatable :: processDistribution, processCoordinates
  logical :: success
  logical, allocatable :: isPeriodic(:)
  real(real64) :: startTime
  character(len = STRING_LENGTH) :: str
  type(t_StencilOperator) :: A
  real(SCALAR_KIND), allocatable :: f(:,:)
  SCALAR_TYPE, allocatable :: u(:,:), v(:,:)

  call MPI_Init(ierror)

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)

  startTime = MPI_Wtime()

  do

     str = trim(discretizationTypes(random(1, size(discretizationTypes))))
     call MPI_Bcast(str, len(str), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)

     call A%setup(trim(str) // " first derivative")

     do nDimensions = 1, 3

        allocate(processDistribution(nDimensions))
        allocate(processCoordinates(nDimensions))
        allocate(isPeriodic(nDimensions))

        direction = random(1, nDimensions)
        call MPI_Bcast(direction, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

        if (random(0, 1) > 0) then
           faceOrientation = +1
        else
           faceOrientation = -1
        end if
        call MPI_Bcast(faceOrientation, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

        isPeriodic = .false.
        processDistribution = 0
        call MPI_Dims_create(numProcs, nDimensions, processDistribution, ierror)
        call MPI_Cart_create(MPI_COMM_WORLD, nDimensions, processDistribution,               &
             isPeriodic, .true., cartesianCommunicator, ierror)
        call A%update(cartesianCommunicator, direction)
        call random_number(A%rhsBoundary1)
        call random_number(A%rhsBoundary2)
        call MPI_Bcast(A%rhsBoundary1, size(A%rhsBoundary1), REAL_TYPE_MPI,                  &
             0, MPI_COMM_WORLD, ierror)
        call MPI_Bcast(A%rhsBoundary2, size(A%rhsBoundary1), REAL_TYPE_MPI,                  &
             0, MPI_COMM_WORLD, ierror)

        globalSize = 1
        do i = 1, nDimensions
           if (i == direction) then
              n = numProcs * (2 * A%boundaryWidth + 1)
              globalSize(i) = random(n, max(2 * n, 2 ** 6))
           else
              globalSize(i) = random(n, max(2 * n, 2 ** 6))
           end if
        end do
        call MPI_Bcast(globalSize, nDimensions, MPI_INTEGER, 0, cartesianCommunicator, ierror)

        call MPI_Comm_rank(cartesianCommunicator, procRank, ierror)
        call MPI_Cart_coords(cartesianCommunicator, procRank,                                &
             nDimensions, processCoordinates, ierror)

        offset = 0
        localSize = 1

        do i = 1, nDimensions
           call pigeonhole(globalSize(i), processDistribution(i),                            &
                processCoordinates(i), offset(i), localSize(i))
        end do

        nComponents = random(1, 9)
        call MPI_Bcast(nComponents, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

        allocate(f(product(localSize), nComponents))
        allocate(u(product(localSize), nComponents))
        allocate(v(product(localSize), nComponents))

        call random_number(f)

        u = f
        call A%apply(u, localSize)

        v = f
        call A%applyAtDomainBoundary(v, localSize, faceOrientation)

        GridLoop: do k = offset(3) + 1, offset(3) + localSize(3)
           do j = offset(2) + 1, offset(2) + localSize(2)
              do i = offset(1) + 1, offset(1) + localSize(1)

                 localIndex = i - offset(1) + localSize(1) * (j - 1 - offset(2) +            &
                      localSize(2) * (k - 1 - offset(3)))

                 if (any(abs(u(localIndex,:) - v(localIndex,:)) > tolerance)) then

                    if (direction == 1 .and. faceOrientation == +1)                          &
                         success = success .and. (i /= 1)
                    if (direction == 1 .and. faceOrientation == -1)                          &
                         success = success .and. (i /= globalSize(1))
                    if (direction == 2 .and. faceOrientation == +1)                          &
                         success = success .and. (j /= 1)
                    if (direction == 2 .and. faceOrientation == -1)                          &
                         success = success .and. (j /= globalSize(2))
                    if (direction == 3 .and. faceOrientation == +1)                          &
                         success = success .and. (k /= 1)
                    if (direction == 3 .and. faceOrientation == -1)                          &
                         success = success .and. (k /= globalSize(3))
                    if (.not. success) exit GridLoop

                 end if

                 if (any(abs(v(localIndex,:)) > tolerance)) then

                    if (direction == 1 .and. faceOrientation == +1)                          &
                         success = success .and. (i == 1)
                    if (direction == 1 .and. faceOrientation == -1)                          &
                         success = success .and. (i == globalSize(1))
                    if (direction == 2 .and. faceOrientation == +1)                          &
                         success = success .and. (j == 1)
                    if (direction == 2 .and. faceOrientation == -1)                          &
                         success = success .and. (j == globalSize(2))
                    if (direction == 3 .and. faceOrientation == +1)                          &
                         success = success .and. (k == 1)
                    if (direction == 3 .and. faceOrientation == -1)                          &
                         success = success .and. (k == globalSize(3))
                    if (.not. success) exit GridLoop

                 end if

              end do !... i = offset(1) + 1, offset(1) + localSize(1)
           end do !... j = offset(2) + 1, offset(2) + localSize(2)
        end do GridLoop !... k = offset(3) + 1, offset(3) + localSize(3)

        call MPI_Barrier(MPI_COMM_WORLD, ierror)
        call MPI_Comm_free(cartesianCommunicator, ierror)

        SAFE_DEALLOCATE(v)
        SAFE_DEALLOCATE(u)
        SAFE_DEALLOCATE(f)

        SAFE_DEALLOCATE(isPeriodic)
        SAFE_DEALLOCATE(processCoordinates)
        SAFE_DEALLOCATE(processDistribution)

        call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL,                            &
             MPI_LAND, MPI_COMM_WORLD, ierror)
        if (.not. success) exit

     end do !... nDimensions = 1, 3

     if (.not. success) exit
     call MPI_Barrier(MPI_COMM_WORLD, ierror)
     if (MPI_Wtime() - startTime > testDuration) exit !... run for `testDuration` seconds

  end do

  call A%cleanup()

  call MPI_Finalize(ierror)
  if (.not. success) stop -1
  stop 0

end program boundary_operator
