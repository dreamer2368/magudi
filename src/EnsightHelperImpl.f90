#include "config.h"

module EnsightHelperImpl

  implicit none
  public

contains

  subroutine writeEnsightCase(this, time, nSpecies)

    ! <<< Internal modules >>>
    use InputHelper, only : getFreeUnit

    ! <<< Derived types >>>
    use EnsightHelper, only : t_Ensight

    implicit none

    ! <<< Arguments >>>
    class(t_Ensight) :: this
    integer, intent(in), optional :: nSpecies
    real(SCALAR_KIND), intent(in) :: time

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, nSpecies_, iunit, ierror
    real(wp), dimension(:), allocatable :: buffer
    character(len = STRING_LENGTH) :: str, name

    ! Update the time info.
    allocate(buffer(this%nOutputTimes))
    buffer(1:this%nOutputTimes) = this%outputTimes(this%nOutputTimes)
    deallocate(this%outputTimes)
    this%nOutputTimes = this%nOutputTimes + 1
    allocate(this%outputTimes(this%nOutputTimes))
    this%outputTimes(1:this%nOutputTimes - 1) = buffer(1:this%nOutputTimes - 1)
    this%outputTimes(this%nOutputTimes) = time
    deallocate(buffer)
  
    ! Open the file.
    str = trim(adjustl(this%directory))//'/'//trim(adjustl(this%filename))
    open(unit = getFreeUnit(iunit), file=trim(str), form="formatted", iostat=ierror,               &
         status="REPLACE")

    ! Write the case.
    str='FORMAT'
    write(iunit, '(a80)') str
    str='type: ensight gold'
    write(iunit, '(a80)') str
    str='GEOMETRY'
    write(iunit, '(a80)') str
    str='model: geometry'
    write(iunit, '(a80)') str
    str='VARIABLE'
    write(iunit, '(a80)') str

    ! Density.
    str='scalar per node: 1 DENSITY DENSITY/DENSITY.******'
    write(iunit,'(a80)') str

    ! Velocity.
    str='vector per node: 1 VELOCITY VELOCITY/VELOCITY.******'
    write(iunit,'(a80)') str

    ! Temperature.
    str='scalar per node: 1 TEMPERATURE TEMPERATURE/TEMPERATURE.******'
    write(iunit,'(a80)') str

    ! Mass fraction.
    nSpecies_ = 0
    if (present(nSpecies)) nSpecies_ = nSpecies
    if (nSpecies_ > 0) then
       do i = 1, nSpecies_
          write(name, "(A,I2.2)") "MASS_FRACTION_", i
          str = 'scalar per node: 1 ' // trim(adjustl(name)) // ' ' //                       &
               trim(adjustl(name)) // '/' // trim(adjustl(name))//'.******'
          write(iunit,'(a80)') str
       end do
    end if

    ! Time section
    str='TIME'
    write(iunit,'(a80)') str
    str='time set: 1'
    write(iunit,'(a80)') str
    str='number of steps:'
    write(iunit,'(a16, i12)') str, this%nOutputTimes
    str='filename start number: 1'
    write(iunit,'(a80)') str
    str='filename increment: 1'
    write(iunit,'(a80)') str
    str='time values:'
    write(iunit,'(a80)') str
    write(iunit, '(10000000(1(ES12.5)))') this%OutputTimes

    ! Close the file
    close(iunit)

  end subroutine writeEnsightCase

  subroutine writeEnsightScalar(scalar, name)

    ! <<< Arguments >>>
    real(KIND=4), dimension(:), intent(in) :: scalar
    character(len = STRING_LENGTH), intent(in) :: name

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i

  end subroutine writeEnsightScalar

  subroutine writeEnsightVector(vector, name)

    ! <<< Arguments >>>
    real(KIND=4), dimension(:,:), intent(in) :: vector
    character(len = STRING_LENGTH), intent(in) :: name

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i

  end subroutine writeEnsightVector

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

subroutine outputEnsight(this, state, comm, gridIndex, mode, time, nSpecies)

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
  integer, intent(in) :: comm, gridIndex, mode
  integer, intent(in), optional :: nSpecies
  real(SCALAR_KIND), intent(in) :: time

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, procRank, nSpecies_, ierror
  character(len = STRING_LENGTH) :: name

  nSpecies_ = 0
  if (present(nSpecies)) nSpecies_ = nSpecies

  ! Get my rank in `comm`.
  call MPI_Comm_rank(comm, procRank, ierror)

  select case(mode)
  case (FORWARD)

     ! Update the case file (root process only).
     if (procRank == 0) call writeEnsightCase(this, time, nSpecies_)

     ! Density.
     name = 'DENSITY'
     this%buffer1_sp = real(state%conservedVariables(:,1), 4)
     call writeEnsightScalar(this%buffer1_sp, name)

     ! Velocity.
     name = 'VELOCITY'
     this%buffer3_sp = 0.0_4
     do i = 1, size(state%velocity, 2)
        this%buffer3_sp(:, i) = real(state%velocity(:, i), 4)
     end do
     call writeEnsightVector(this%buffer3_sp, name)

     ! Temperature.
     name = 'TEMPERATURE'
     this%buffer1_sp = real(state%temperature(:,1), 4)
     call writeEnsightScalar(this%buffer1_sp, name)

     ! Mass fraction.
     if (nSpecies_ > 0) then
        do i = 1, size(state%massFraction, 2)
           write(name, "(A,I2.2)") "MASS_FRACTION_", i
           this%buffer1_sp = real(state%massFraction(:,i), 4)
           call writeEnsightScalar(this%buffer1_sp, name)
        end do
     end if

  case (ADJOINT)

  end select

end subroutine outputEnsight
