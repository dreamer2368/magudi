#include "config.h"

subroutine readPatchInterfaceInformation(region)

  ! <<< Derived types >>>
  use Region_mod, only : t_Region

  ! <<< Internal modules >>>
  use InputHelper, only : getOption
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer :: i, j, k, l, nDimensions
  character(len = STRING_LENGTH) :: key, str, message

  if (.not. allocated(region%patchData)) return

  SAFE_DEALLOCATE(region%patchInterfaces)
  SAFE_DEALLOCATE(region%interfaceIndexReorderings)

  allocate(region%patchInterfaces(size(region%patchData)), source = 0)
  allocate(region%interfaceIndexReorderings(3, size(region%patchData)), source = 0)

  nDimensions = size(region%globalGridSizes, 1)

  ! Read interface information:

  do i = 1, size(region%patchData)

     write(key, '(A)') "patches/" // trim(region%patchData(i)%name) // "/conforms_with"
     str = getOption(key, "")

     if (len_trim(str) > 0) then

        do j = 1, size(region%patchData)
           if (trim(str) == trim(region%patchData(j)%name)) region%patchInterfaces(i) = j
        end do
        if (region%patchInterfaces(i) == 0) then
           write(message, '(5A)') "Invalid interface specification for patch '",             &
                trim(region%patchData(i)%name), "': no patch found matching the name '",     &
                trim(str), "'!"
           call gracefulExit(region%comm, message)
        end if

        do j = 1, 3
           region%interfaceIndexReorderings(j,i) = j
           if (j <= nDimensions) then
              write(key, '(A,I1)') "patches/" // trim(region%patchData(i)%name) //           &
                   "/interface_index", j
              region%interfaceIndexReorderings(j,i) = getOption(trim(key), j)
              if (j == 3 .and. region%interfaceIndexReorderings(j,i) /= 3) then
                 write(message, '(3A)') "Interface index reordering for patch '",            &
                      trim(region%patchData(i)%name), "' is currently not supported!"
                 call gracefulExit(region%comm, message)
              end if
           end if
        end do

        if (.not. all(region%interfaceIndexReorderings(1:nDimensions,i) /= 0 .and.           &
             abs(region%interfaceIndexReorderings(1:nDimensions,i)) <= nDimensions)) then
           write(message, '(3A)') "Invalid interface index reordering for patch '",          &
                trim(region%patchData(i)%name), "'!"
           call gracefulExit(region%comm, message)
        end if

     end if

  end do

  ! Commutativity of interfaces:

  do i = 1, size(region%patchData)
     if (region%patchInterfaces(i) == 0) cycle
     j = region%patchInterfaces(i)

     if (region%patchInterfaces(j) == 0) then

        region%patchInterfaces(j) = i

        do l = 1, 3
           do k = 1, 3
              if (abs(region%interfaceIndexReorderings(k,i)) == l) then
                 region%interfaceIndexReorderings(l,j) =                                     &
                      sign(k, region%interfaceIndexReorderings(k,i))
                 exit
              end if
           end do
        end do

     else if (region%patchInterfaces(j) /= i) then

        write(message, '(3A)') "Invalid interface specification for patch '",                &
             trim(region%patchData(j)%name), "': violates commutativity!"
        call gracefulExit(region%comm, message)

     end if

  end do

end subroutine readPatchInterfaceInformation

subroutine exchangeInterfaceData(region, mode)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region
  integer, intent(in) :: mode

  ! <<< Local variables >>>
  integer :: i, j, iRequest, mpiTag, procRank, ierror
  integer, allocatable :: mpiSendRequests(:)
  class(t_Patch), pointer :: patch => null()
  class(t_BlockInterfacePatch), pointer :: blockInterfacePatch => null()

  if (.not. allocated(region%patchInterfaces)) return

  call MPI_Comm_rank(region%comm, procRank, ierror)

  if (any(region%patchMasterRanks == procRank))                                              &
       allocate(mpiSendRequests(count(region%patchInterfaces > 0 .and.                       &
       region%patchMasterRanks == procRank)), source = MPI_REQUEST_NULL)

  iRequest = 0

  do i = 1, size(region%patchInterfaces)
     if (procRank == region%patchMasterRanks(i) .and. region%patchInterfaces(i) > 0) then

        assert(allocated(region%patchFactories))

        mpiTag = i + size(region%patchInterfaces) * (region%patchInterfaces(i) - 1)

        nullify(blockInterfacePatch)

        do j = 1, size(region%patchFactories)
           call region%patchFactories(j)%connect(patch)
           if (.not. associated(patch)) cycle
           if (patch%index /= i) cycle
           select type (patch)
           class is (t_BlockInterfacePatch)
              blockInterfacePatch => patch
              exit
           end select
        end do

        assert(associated(blockInterfacePatch))
        assert(allocated(blockInterfacePatch%sendBuffer))

        iRequest = iRequest + 1
        call MPI_Isend(blockInterfacePatch%sendBuffer, size(blockInterfacePatch%sendBuffer), &
             SCALAR_TYPE_MPI, region%patchMasterRanks(region%patchInterfaces(i)), mpiTag,    &
             region%comm, mpiSendRequests(iRequest), ierror)

     end if
  end do

  do i = 1, size(region%patchInterfaces)
     if (procRank == region%patchMasterRanks(i) .and. region%patchInterfaces(i) > 0) then

        assert(allocated(region%patchFactories))

        mpiTag = region%patchInterfaces(i) + size(region%patchInterfaces) * (i - 1)

        nullify(blockInterfacePatch)

        do j = 1, size(region%patchFactories)
           call region%patchFactories(j)%connect(patch)
           if (.not. associated(patch)) cycle
           if (patch%index /= i) cycle
           select type (patch)
           class is (t_BlockInterfacePatch)
              blockInterfacePatch => patch
              exit
           end select
        end do

        assert(associated(blockInterfacePatch))
        assert(allocated(blockInterfacePatch%receiveBuffer))

        call MPI_Recv(blockInterfacePatch%receiveBuffer,                                     &
             size(blockInterfacePatch%receiveBuffer), SCALAR_TYPE_MPI,                       &
             region%patchMasterRanks(region%patchInterfaces(i)),                             &
             mpiTag, region%comm, MPI_STATUS_IGNORE, ierror)

     end if
  end do

  if (allocated(mpiSendRequests))                                                            &
       call MPI_Waitall(size(mpiSendRequests), mpiSendRequests, MPI_STATUSES_IGNORE, ierror)

  do i = 1, size(region%patchInterfaces)
     if (procRank == region%patchMasterRanks(i) .and. region%patchInterfaces(i) > 0) then

        assert(allocated(region%patchFactories))

        nullify(blockInterfacePatch)

        do j = 1, size(region%patchFactories)
           call region%patchFactories(j)%connect(patch)
           if (.not. associated(patch)) cycle
           if (patch%index /= i) cycle
           select type (patch)
           class is (t_BlockInterfacePatch)
              blockInterfacePatch => patch
              exit
           end select
        end do

        assert(associated(blockInterfacePatch))

        call blockInterfacePatch%reshapeReceivedData(region%interfaceIndexReorderings(:,i))

     end if
     call MPI_Barrier(region%comm, ierror)
  end do

  SAFE_DEALLOCATE(mpiSendRequests)

end subroutine exchangeInterfaceData

subroutine checkFunctionContinuityAtInterfaces(region, tolerance)

  ! <<< External modules >>>
  use MPI
  use, intrinsic :: iso_fortran_env

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD

  ! <<< Public members >>>
  use InterfaceHelper, only : exchangeInterfaceData

  ! <<< Internal modules >>>
  use ErrorHandler, only : gracefulExit, writeAndFlush
  use RandomNumber, only : initializeRandomNumberGenerator, random

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region
  real(SCALAR_KIND), intent(in) :: tolerance

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(real64), parameter :: checkDuration = 2.0_real64
  character(len = STRING_LENGTH), parameter :: functionTypes(1) = (/ "polynomial" /)
  integer :: i, j, k, l, nDimensions, nCoefficients, errorRank, ierror
  real(real64) :: startTime
  character(len = STRING_LENGTH) :: functionType, message
  class(t_Patch), pointer :: patch => null()
  real(wp), allocatable :: coefficients(:,:)
  SCALAR_TYPE, allocatable :: patchCoordinates(:,:), patchFunction(:,:), gatherBuffer(:,:),  &
       scatterBuffer(:,:), patchFunctionReceived(:,:)

  errorRank = -1

  if (.not. allocated(region%patchInterfaces)) return

  write(message, "(A)") "Checking continuity at interfaces..."
  call writeAndFlush(region%comm, output_unit, message, advance = 'no')

  call initializeRandomNumberGenerator()

  nDimensions = size(region%globalGridSizes, 1)

  startTime = MPI_Wtime()

  do

     functionType = trim(functionTypes(random(1, size(functionTypes))))
     call MPI_Bcast(functionType, len(functionType), MPI_CHARACTER, 0, region%comm, ierror)

     select case (trim(functionType))

     case ("polynomial")
        nCoefficients = random(1, 4) + 1 !... up to 4th order polynomials
        call MPI_Bcast(nCoefficients, 1, MPI_INTEGER, 0, region%comm, ierror)
     end select

     allocate(coefficients(nDimensions, nCoefficients))
     call random_number(coefficients)
     call MPI_Bcast(coefficients, size(coefficients), REAL_TYPE_MPI, 0, region%comm, ierror)

     if (allocated(region%patchFactories)) then
        do i = 1, size(region%patchFactories)
           call region%patchFactories(i)%connect(patch)
           if (.not. associated(patch)) cycle
           do j = 1, size(region%states)
              if (patch%gridIndex /= region%grids(j)%index .or. patch%nPatchPoints < 0) cycle
              select type (patch)
                 class is (t_BlockInterfacePatch)

                 allocate(patchCoordinates(patch%nPatchPoints, nDimensions))
                 allocate(patchFunction(patch%nPatchPoints, 1))
                 patchFunction = 0.0_wp
                 call patch%collect(region%grids(j)%coordinates, patchCoordinates)
                 if (allocated(patch%sendBuffer))                                            &
                      allocate(gatherBuffer(patch%nPatchPoints, 1))

                 select case (trim(functionType))

                 case ("polynomial")
                    do k = 1, nCoefficients
                       do l = 1, nDimensions
                          patchFunction(:,1) = patchFunction(:,1) +                          &
                               coefficients(l,k) * patchCoordinates(:,l) ** real(k - 1, wp)
                       end do
                    end do

                 end select

                 call patch%gatherData(patchFunction, gatherBuffer)
                 if (allocated(patch%sendBuffer)) patch%sendBuffer(:,1) = gatherBuffer(:,1)
                 SAFE_DEALLOCATE(gatherBuffer)
                 SAFE_DEALLOCATE(patchFunction)
                 SAFE_DEALLOCATE(patchCoordinates)

              end select
           end do
        end do
     end if

     call exchangeInterfaceData(region, FORWARD)
     call MPI_Barrier(region%comm, ierror)

     if (allocated(region%patchFactories)) then
        do i = 1, size(region%patchFactories)
           call region%patchFactories(i)%connect(patch)
           if (.not. associated(patch)) cycle
           do j = 1, size(region%states)
              if (patch%gridIndex /= region%grids(j)%index .or. patch%nPatchPoints < 0) cycle
              select type (patch)
                 class is (t_BlockInterfacePatch)

                 allocate(patchCoordinates(patch%nPatchPoints, nDimensions))
                 allocate(patchFunction(patch%nPatchPoints, 1))
                 allocate(patchFunctionReceived(patch%nPatchPoints, 1))
                 patchFunction = 0.0_wp
                 call patch%collect(region%grids(j)%coordinates, patchCoordinates)
                 if (allocated(patch%receiveBuffer))                                         &
                      allocate(scatterBuffer(patch%nPatchPoints, 1))

                 select case (trim(functionType))

                 case ("polynomial")
                    do k = 1, nCoefficients
                       do l = 1, nDimensions
                          patchFunction(:,1) = patchFunction(:,1) +                          &
                               coefficients(l,k) * patchCoordinates(:,l) ** real(k - 1, wp)
                       end do
                    end do

                 end select

                 if (allocated(patch%receiveBuffer))                                         &
                      scatterBuffer(:,1) = patch%receiveBuffer(:,1)
                 call patch%scatterData(scatterBuffer, patchFunctionReceived)

                 if (any(abs(patchFunctionReceived - patchFunction) > tolerance)) then
                    call MPI_Comm_rank(region%comm, errorRank, ierror)
                    write(message, '(3A,2(E13.6,A),2A)')                                     &
                         "Interface continuity check failed for a random '",                 &
                         trim(functionType), "' function. Difference in value ",             &
                         maxval(abs(patchFunctionReceived - patchFunction)),                 &
                         " exceeds tolerance ", tolerance, " on patch '",                    &
                         trim(region%patchData(patch%index)%name), "'!"
                 end if

                 SAFE_DEALLOCATE(scatterBuffer)
                 SAFE_DEALLOCATE(patchFunctionReceived)
                 SAFE_DEALLOCATE(patchFunction)
                 SAFE_DEALLOCATE(patchCoordinates)

              end select
           end do
        end do
     end if

     SAFE_DEALLOCATE(coefficients)

     if (MPI_Wtime() - startTime > checkDuration) exit

     call MPI_Allreduce(MPI_IN_PLACE, errorRank, 1, MPI_INTEGER, MPI_MAX, region%comm, ierror)
     if (errorRank /= -1) exit

  end do

  call MPI_Allreduce(MPI_IN_PLACE, errorRank, 1, MPI_INTEGER, MPI_MAX, region%comm, ierror)
  if (errorRank /= -1) then
     write(message, "(A)") " failed!"
     call writeAndFlush(region%comm, output_unit, message)
     call MPI_Bcast(message, len(message), MPI_CHARACTER, errorRank, region%comm, ierror)
     call gracefulExit(region%comm, message)
  else
     write(message, "(A)") " done!"
     call writeAndFlush(region%comm, output_unit, message)
  end if

end subroutine checkFunctionContinuityAtInterfaces
