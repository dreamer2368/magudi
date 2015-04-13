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
  class(t_Patch), pointer :: patch => null()
  class(t_BlockInterfacePatch), pointer :: blockInterfacePatch => null()

  if (.not. allocated(region%patchInterfaces)) return

  call MPI_Comm_rank(region%comm, procRank, ierror)

  do i = 1, size(region%patchInterfaces)
     if (any(region%patchMasterRanks == procRank) .and. region%patchInterfaces(i) > 0 .and.  &
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

        assert(associated(blockInterfacePatch))
        call MPI_Sendrecv_replace(blockInterfacePatch%exchangeBuffer,                        &
             size(blockInterfacePatch%exchangeBuffer), SCALAR_TYPE_MPI,                      &
             region%patchMasterRanks(region%patchInterfaces(i)), mpiSendTag,                 &
             region%patchMasterRanks(region%patchInterfaces(i)), mpiRecvTag,                 &
             region%comm, MPI_STATUS_IGNORE, ierror)

     end if
  end do

  if (allocated(region%patchFactories)) then
     do j = 1, size(region%patchFactories)
        call region%patchFactories(j)%connect(patch)
        if (.not. associated(patch)) cycle
        select type (patch)
           class is (t_BlockInterfacePatch)
           blockInterfacePatch => patch
           if (allocated(blockInterfacePatch%exchangeBuffer))                                &
                call blockInterfacePatch%reshapeExchangeBuffer(                              &
                region%interfaceIndexReorderings(:,blockInterfacePatch%index))
        end select
     end do
  end if

  call MPI_Barrier(region%comm, ierror)

end subroutine exchangeInterfaceData
