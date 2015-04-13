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

subroutine exchangeInterfaceData(region)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer :: i, j, mpiSendTag, mpiRecvTag, procRank, ierror
  integer, allocatable :: mpiSendRequests(:), mpiRecvRequests(:)
  class(t_Patch), pointer :: patch => null()
  class(t_BlockInterfacePatch), pointer :: blockInterfacePatch => null()

  if (.not. allocated(region%patchInterfaces)) return

  call MPI_Comm_rank(region%comm, procRank, ierror)

  allocate(mpiSendRequests(size(region%patchInterfaces)), source = MPI_REQUEST_NULL)
  allocate(mpiRecvRequests(size(region%patchInterfaces)), source = MPI_REQUEST_NULL)

  do i = 1, size(region%patchInterfaces)
     if (procRank == region%patchMasterRanks(i) .and. region%patchInterfaces(i) > 0 .and.    &
          allocated(region%patchFactories)) then

        nullify(blockInterfacePatch)
        mpiSendTag = i + size(region%patchInterfaces) * (region%patchInterfaces(i) - 1)
        mpiRecvTag = region%patchInterfaces(i) + size(region%patchInterfaces) * (i - 1)

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

        if (associated(blockInterfacePatch)) then
           call MPI_Isend(blockInterfacePatch%sendBuffer,                                    &
                size(blockInterfacePatch%sendBuffer), SCALAR_TYPE_MPI,                       &
                region%patchMasterRanks(region%patchInterfaces(i)), mpiSendTag,              &
                region%comm, mpiSendRequests(i), ierror)
           call MPI_Irecv(blockInterfacePatch%receiveBuffer,                                 &
                size(blockInterfacePatch%receiveBuffer), SCALAR_TYPE_MPI,                    &
                region%patchMasterRanks(region%patchInterfaces(i)), mpiRecvTag,              &
                region%comm, mpiRecvRequests(i), ierror)
        end if

     end if
  end do

  call MPI_Waitall(size(mpiSendRequests), mpiSendRequests, MPI_STATUSES_IGNORE, ierror)
  call MPI_Waitall(size(mpiRecvRequests), mpiRecvRequests, MPI_STATUSES_IGNORE, ierror)

  if (allocated(region%patchFactories)) then
     do j = 1, size(region%patchFactories)
        call region%patchFactories(j)%connect(patch)
        if (.not. associated(patch)) cycle
        select type (patch)
           class is (t_BlockInterfacePatch)
           blockInterfacePatch => patch
           if (allocated(blockInterfacePatch%receiveBuffer))                                 &
                call blockInterfacePatch%reshapeReceivedData(                                &
                region%interfaceIndexReorderings(:,blockInterfacePatch%index))
        end select
     end do
  end if

  SAFE_DEALLOCATE(mpiRecvRequests)
  SAFE_DEALLOCATE(mpiSendRequests)

  call MPI_Barrier(region%comm, ierror)

end subroutine exchangeInterfaceData

subroutine checkInterfaceContinuity(region, tolerance, success)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  ! <<< Public members >>>
  use InterfaceHelper, only : exchangeInterfaceData

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region
  real(SCALAR_KIND), intent(in) :: tolerance
  logical, intent(out) :: success

  ! <<< Local variables >>>
  integer :: i, j, k, nDimensions, nGlobalPatchPoints, ierror
  integer, allocatable :: nComponents(:)
  class(t_Patch), pointer :: patch => null()
  class(t_BlockInterfacePatch), pointer :: blockInterfacePatch => null()
  SCALAR_TYPE, allocatable :: patchCoordinates(:,:), interfaceCoordinates(:,:)

  success = .true.

  if (.not. allocated(region%patchInterfaces)) return

  nDimensions = size(region%globalGridSizes, 1)
  allocate(nComponents(size(region%patchInterfaces)), source = 0)

  do i = 1, size(region%grids)
     if (allocated(region%patchFactories)) then
        do j = 1, size(region%patchFactories)
           call region%patchFactories(j)%connect(patch)
           if (.not. associated(patch)) cycle
           if (patch%gridIndex /= region%grids(i)%index .or. patch%nPatchPoints <= 0) cycle

           select type (patch)
           class is (t_BlockInterfacePatch)

              blockInterfacePatch => patch
              k = blockInterfacePatch%index

              if (allocated(blockInterfacePatch%sendBuffer)) then
                 nComponents(k) = size(blockInterfacePatch%sendBuffer, 2) !... save for later.
                 assert(allocated(blockInterfacePatch%receiveBuffer))
                 assert(size(blockInterfacePatch%receiveBuffer, 2) == nComponents(k))
              end if

              nGlobalPatchPoints = product(blockInterfacePatch%globalSize)

              SAFE_DEALLOCATE(blockInterfacePatch%sendBuffer)
              SAFE_DEALLOCATE(blockInterfacePatch%receiveBuffer)

              allocate(blockInterfacePatch%sendBuffer(nGlobalPatchPoints, nDimensions))
              allocate(blockInterfacePatch%receiveBuffer(nGlobalPatchPoints, nDimensions))

              allocate(patchCoordinates(blockInterfacePatch%nPatchPoints, nDimensions))
              call blockInterfacePatch%collect(region%grids(i)%coordinates, patchCoordinates)
              call blockInterfacePatch%gatherData(patchCoordinates,                          &
                   blockInterfacePatch%sendBuffer)

              SAFE_DEALLOCATE(patchCoordinates)

           end select

        end do !... j = 1, size(region%patchFactories)
     end if !... allocated(region%patchFactories)
  end do !... i = 1, size(region%grids)

  call exchangeInterfaceData(region)

  do i = 1, size(region%grids)
     if (allocated(region%patchFactories)) then
        do j = 1, size(region%patchFactories)
           call region%patchFactories(j)%connect(patch)
           if (.not. associated(patch)) cycle
           if (patch%gridIndex /= region%grids(i)%index .or. patch%nPatchPoints <= 0) cycle

           select type (patch)
           class is (t_BlockInterfacePatch)

              blockInterfacePatch => patch
              k = blockInterfacePatch%index

              allocate(interfaceCoordinates(blockInterfacePatch%nPatchPoints, nDimensions))
              call blockInterfacePatch%scatterData(blockInterfacePatch%receiveBuffer,        &
                   interfaceCoordinates)

              allocate(patchCoordinates(blockInterfacePatch%nPatchPoints, nDimensions))
              call blockInterfacePatch%collect(region%grids(i)%coordinates, patchCoordinates)

              if (any(abs(patchCoordinates - interfaceCoordinates) > tolerance))             &
                   success = .false.

              SAFE_DEALLOCATE(patchCoordinates)
              SAFE_DEALLOCATE(interfaceCoordinates)

              SAFE_DEALLOCATE(blockInterfacePatch%sendBuffer)
              SAFE_DEALLOCATE(blockInterfacePatch%receiveBuffer)

              nGlobalPatchPoints = product(blockInterfacePatch%globalSize)

              if (nComponents(k) > 0) then
                 allocate(blockInterfacePatch%sendBuffer(nGlobalPatchPoints, nComponents(k)))
                 allocate(blockInterfacePatch%receiveBuffer(nGlobalPatchPoints,              &
                      nComponents(k)))
              end if

           end select

           if (.not. success) exit
        end do !... j = 1, size(region%patchFactories)
     end if !... allocated(region%patchFactories)
  end do !... i = 1, size(region%grids)

  call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierror)

end subroutine checkInterfaceContinuity
