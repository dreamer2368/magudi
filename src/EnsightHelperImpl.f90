#include "config.h"

module EnsightHelperImpl

  implicit none
  public

contains

  pure subroutine swapIntegerEndianness_(a)

    ! <<< Arguments >>>
    integer, intent(inout) :: a

    ! <<< Local variables >>>
    integer :: b

    call mvbits(a,  0, 8, b, 24)
    call mvbits(a,  8, 8, b, 16)
    call mvbits(a, 16, 8, b,  8)
    call mvbits(a, 24, 8, b,  0)

    a = b

  end subroutine swapIntegerEndianness_

  pure subroutine swapScalarEndianness_(a)

    ! <<< Arguments >>>
    SCALAR_TYPE, intent(inout) :: a

    ! <<< Local variables >>>
    integer :: i
    integer(kind = 1) :: b(SIZEOF_SCALAR), c(SIZEOF_SCALAR)

    b = transfer(a, b)
    do i = 1, SIZEOF_SCALAR
       c(i) = b(SIZEOF_SCALAR+1-i)
    end do
    a = transfer(c, a)

  end subroutine swapScalarEndianness_

end module EnsightHelperImpl

subroutine setupEnsight(this, comm, gridIndex, localSize, globalSize, offset, time)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use EnsightHelper, only : t_Ensight

  implicit none

  ! <<< Arguments >>>
  class(t_Ensight) :: this
  integer, intent(in) :: comm, gridIndex, localSize(3), globalSize(3), offset(3)
  real(SCALAR_KIND), intent(in) :: time

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer, parameter :: maxrecs = 10000
  integer :: i, j, procRank, iunit, iostat, ierror
  logical :: fileExists
  character(len = STRING_LENGTH) :: line
  logical :: lineFound

  ! Get my rank in `comm`.
  call MPI_Comm_rank(comm, procRank, ierror)

  ! Open the file
  this%directory = 'ensight-3D'
  write(this%filename, '(A,I2.2,A)') "magudi_grid", gridIndex, ".case"
  if (procRank == 0) then
     inquire(file = trim(this%directory)//'/'//trim(this%filename), exist = fileExists)
     if (fileExists) then
        ! Get the number of time values.
        this%nOutputTimes = 0
        lineFound = .false.
        open(unit = iunit, file = trim(this%directory)//'/'//trim(this%filename))
        do i = 1, maxrecs
           read(iunit,'(A)', iostat = iostat) line

           if ( iostat < 0 ) exit 

           if (lineFound) this%nOutputTimes = this%nOutputTimes + 1

           if (.not. lineFound .and. trim(adjustl(line)) == 'time values:') lineFound = .true.

           if (i == maxrecs) then
              write(*,*) 'Error:  Maximum number of records exceeded reading ',              &
                   trim(this%directory)//'/'//trim(this%filename)
              stop
           end if
        end do
        rewind(iunit)

        ! Read in the time values.
        allocate(this%outputTimes(this%nOutputTimes))
        lineFound = .false.
        j = 0
        do i = 1, maxrecs
           if (.not. lineFound) read(iunit,'(A)', iostat = iostat) line

           if ( iostat < 0 ) exit 

           if (lineFound) then
              j = j + 1
              read (iunit, *, iostat = iostat) this%outputTimes(j)
           end if

           if (.not. lineFound .and. trim(adjustl(line)) == 'time values:') lineFound = .true.

           if (i == maxrecs) then
              write(*,*) 'Error:  Maximum number of records exceeded reading ',              &
                   trim(this%directory)//'/'//trim(this%filename)
              stop
           end if
        end do
        close(iunit)

        ! Remove future time
        future: do i = 1, size(this%outputTimes)
           if (this%outputTimes(i) >= time * 0.99999_wp) then
              this%nOutputTimes = i - 1
              exit future
           end if
        end do future
     else
        ! Create directory
        call Execute_Command_Line('mkdir -p ' // trim(adjustl(this%directory)))
        ! Set the time
        this%nOutputTimes = 0
        allocate(this%outputTimes(1))
     end if
  end if

  ! Create the view.
  this%dataSize = product(localSize)
  this%gdataSize= product(globalSize)
  call MPI_TYPE_CREATE_SUBARRAY(3, globalSize, localSize, offset, MPI_ORDER_FORTRAN,         &
       MPI_REAL, this%fileview, ierror)
  call MPI_TYPE_COMMIT(this%fileview, ierror)

  ! Allocate single-precision buffer arrays.
  allocate(this%buffer1_sp(this%dataSize))
  allocate(this%buffer3_sp(this%dataSize, 3))

end subroutine setupEnsight

subroutine outputEnsight(this, state, gridIndex, mode, time)

  ! <<< Derived types >>>
  use EnsightHelper, only : t_Ensight

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use State_mod, only : t_State

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

 ! <<< Private members >>>
  use EnsightHelperImpl

  implicit none

  ! <<< Arguments >>>
  class(t_Ensight) :: this
  class(t_State) :: state
  integer, intent(in) :: gridIndex, mode
  real(SCALAR_KIND), intent(in) :: time

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, mpiFileHandle, ierror

  print *, 'hiii'

end subroutine outputEnsight
