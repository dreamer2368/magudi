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

  if (allocated(region%patchFactories)) then

     do i = 1, size(region%patchFactories)
        call region%patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle

        do j = 1, size(region%states)
           if (patch%gridIndex /= region%grids(j)%index) cycle
           select type (patch)
           class is (t_BlockInterfacePatch)
              call patch%collectInterfaceData(mode, region%simulationFlags,                    &
                   region%solverOptions, region%grids(j), region%states(j))
           end select
        end do

     end do

  end if

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

subroutine checkInterfaceContinuity(region, tolerance, success)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Patch_mod, only : t_Patch
  use Region_mod, only : t_Region
  use BlockInterfacePatch_mod, only : t_BlockInterfacePatch

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD

  ! <<< Public members >>>
  use InterfaceHelper, only : exchangeInterfaceData

  implicit none

  ! <<< Arguments >>>
  class(t_Region) :: region
  real(SCALAR_KIND), intent(in) :: tolerance
  logical, intent(out) :: success

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, ierror
  class(t_Patch), pointer :: patch => null()

  success = .true.

  if (.not. allocated(region%patchInterfaces)) return

  do i = 1, size(region%states)
     region%states(i)%rightHandSide = 0.0_wp
  end do

  call exchangeInterfaceData(region, FORWARD)

  if (allocated(region%patchFactories)) then
     do i = 1, size(region%patchFactories)
        call region%patchFactories(i)%connect(patch)
        if (.not. associated(patch)) cycle

        select type (patch)
        class is (t_BlockInterfacePatch)

           do j = 1, size(region%states)
              if (patch%gridIndex /= region%grids(j)%index) cycle
              call patch%updateRhs(FORWARD, region%simulationFlags, region%solverOptions,    &
                   region%grids(j), region%states(j))
           end do

        end select
     end do
  end if

  call MPI_Barrier(region%comm, ierror)

  do i = 1, size(region%states)
     if (maxval(abs(region%states(i)%rightHandSide)) > tolerance) then
        success = .false.
        exit
     end if
  end do

  call MPI_Allreduce(MPI_IN_PLACE, success, 1, MPI_LOGICAL, MPI_LAND, region%comm, ierror)

end subroutine checkInterfaceContinuity
