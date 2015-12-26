#include "config.h"

module EnsightHelperImpl

  implicit none
  public

contains

  subroutine writeEnsightCase(this, time, gridIndex, nSpecies)

    ! <<< Internal modules >>>
    use InputHelper, only : getFreeUnit

    ! <<< Derived types >>>
    use EnsightHelper, only : t_Ensight

    implicit none

    ! <<< Arguments >>>
    class(t_Ensight) :: this
    integer, intent(in) :: gridIndex
    integer, intent(in), optional :: nSpecies
    real(SCALAR_KIND), intent(in) :: time

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, nSpecies_, iunit, ierror
    real(wp), dimension(:), allocatable :: buffer
    character(len = STRING_LENGTH) :: str, name

    ! Update the time info.
    if (this%nOutputTimes > 0) then
       allocate(buffer(this%nOutputTimes))
       buffer(1:this%nOutputTimes) = this%outputTimes(1:this%nOutputTimes)
       deallocate(this%outputTimes)
       this%nOutputTimes = this%nOutputTimes + 1
       allocate(this%outputTimes(this%nOutputTimes))
       this%outputTimes(1:this%nOutputTimes - 1) = buffer(1:this%nOutputTimes - 1)
       this%outputTimes(this%nOutputTimes) = time
       deallocate(buffer)
    else
       deallocate(this%outputTimes)
       this%nOutputTimes = 1
       allocate(this%outputTimes(this%nOutputTimes))
       this%outputTimes(this%nOutputTimes) = time
    end if
  
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
    write(str, '(A,I2.2)') "model: geometry_", gridIndex
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
    do i = 1, this%nOutputTimes
       write(iunit, '(ES12.5)') this%OutputTimes(i)
    end do

    ! Close the file
    close(iunit)

  end subroutine writeEnsightCase

  subroutine writeEnsightGeometry(this, grid, gridIndex)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env, only : output_unit

    ! <<< Internal modules >>>
    use ErrorHandler, only : writeAndFlush

    ! <<< Derived types >>>
    use EnsightHelper, only : t_Ensight
    use Grid_mod, only : t_Grid

    implicit none

    ! <<< Arguments >>>
    class(t_Ensight) :: this
    class(t_Grid) :: grid
    integer, intent(in) :: gridIndex

    ! <<< Local variables >>>
    integer :: i, procRank, iunit, ierror
    integer(kind = MPI_OFFSET_KIND) :: offset
    character(len = 80) :: buffer, filename, message
    logical :: fileExists, useIblank

    ! Get my rank in `comm`.
    call MPI_Comm_rank(grid%comm, procRank, ierror)

    ! Decide if iblank should be written.
    useIblank = .false.
    i = minval(grid%iblank)
    call MPI_Allreduce(MPI_IN_PLACE, i, 1, MPI_INTEGER, MPI_MIN, grid%comm, ierror)
    if (i < 1) useIblank = .true.

    ! Open the file to write.
    filename = trim(adjustl(this%directory)) // "/geometry_"
    write(filename, '(A,I2.2)') trim(adjustl(filename)), gridIndex

    ! Output to screen.
    write(message, '(3A)') "Writing '", trim(filename), "'..."
    call writeAndFlush(grid%comm, output_unit, message, advance = 'no')

    inquire(file = filename, exist = fileExists)
    if (fileExists .and. procRank == 0) call MPI_FILE_DELETE(filename, MPI_INFO_NULL, ierror)
    call MPI_FILE_OPEN(grid%comm, trim(adjustl(filename)),                                   &
         IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, iunit, ierror)

    ! Write the geometry header.
    if (procRank == 0) then
       buffer = 'C Binary'
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       buffer = 'Ensight Gold Geometry File'
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       buffer = 'Curvilinear Geometry from Magudi'
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       buffer = 'node id off'
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       buffer = 'element id off'
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       buffer = 'part'
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       i = 1
       call MPI_FILE_WRITE(iunit, i, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
       write(buffer, '(A,I2.2)') "Grid_", gridIndex
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       if (useIblank) then
          buffer = 'block curvilinear iblanked'
       else
          buffer = 'block'
       end if
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       call MPI_FILE_WRITE(iunit, grid%globalSize(1), 1, MPI_INTEGER, MPI_STATUS_IGNORE,     &
            ierror)
       call MPI_FILE_WRITE(iunit, grid%globalSize(2), 1, MPI_INTEGER, MPI_STATUS_IGNORE,     &
            ierror)
       call MPI_FILE_WRITE(iunit, grid%globalSize(3), 1, MPI_INTEGER, MPI_STATUS_IGNORE,     &
            ierror)
    end if

    ! Store coordinates in single-precision array.
    this%buffer3_sp = 0.0_4
    do i = 1, size(grid%coordinates, 2)
       this%buffer3_sp(:, i) = real(grid%coordinates(:, i), 4)
    end do

    ! Write the grid coordinates.
    offset = 80 * 8 + 4 * 4
    do i = 1, 3
       call MPI_File_set_view(iunit, offset, MPI_REAL, this%fileview, "native",              &
            MPI_INFO_NULL, ierror)
       call MPI_File_write_all(iunit, this%buffer3_sp(:, i), this%datasize,                  &
            MPI_REAL, MPI_STATUS_IGNORE, ierror)
       offset = offset + 4 * int(this%gDataSize, MPI_OFFSET_KIND)
    end do

    ! Write the iblank.
    if (useIblank) then
       call MPI_File_set_view(iunit, offset, MPI_INTEGER, this%fileview, "native",           &
            MPI_INFO_NULL, ierror)
       call MPI_File_write_all(iunit, grid%iblank, this%datasize, MPI_INTEGER,               &
            MPI_STATUS_IGNORE, ierror)
    end if

    ! Close the file.
    call MPI_File_close(iunit, ierror)

    write(message, '(A)') " done!"
    call writeAndFlush(grid%comm, output_unit, message)

  end subroutine writeEnsightGeometry

  subroutine writeEnsightScalar(this, name, comm)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env, only : output_unit

    ! <<< Internal modules >>>
    use ErrorHandler, only : writeAndFlush

    ! <<< Derived types >>>
    use EnsightHelper, only : t_Ensight

    implicit none

    ! <<< Arguments >>>
    class(t_Ensight) :: this
    integer, intent(in) :: comm
    character(len = STRING_LENGTH), intent(in) :: name

    ! <<< Local variables >>>
    integer :: i, procRank, iunit, ierror
    integer(kind = MPI_OFFSET_KIND) :: offset
    character(len = 80) :: buffer, filename, message
    logical :: fileExists

    ! Get my rank in `comm`.
    call MPI_Comm_rank(comm, procRank, ierror)

    ! Generate the file.
    if (procRank == 0) then
       call Execute_Command_Line('mkdir -p ' // trim(adjustl(this%directory)) // '/' //      &
            trim(adjustl(name)))
       filename = trim(adjustl(this%directory)) // '/' // trim(adjustl(name)) // "/" //      &
            trim(adjustl(name)) // "."
       write(filename(len_trim(filename) + 1 : len_trim(filename) + 6), '(i6.6)')            &
            this%nOutputTimes
    end if

    ! Communicate the filename.
    call MPI_Bcast(filename, 80, MPI_CHARACTER, 0, comm, ierror)     

    ! Output to screen.
    write(message, '(3A)') "Writing '", trim(filename), "'..."
    call writeAndFlush(comm, output_unit, message, advance = 'no')

    ! Open the file.
    inquire(file = filename, exist = fileExists)
    if (fileExists .and. procRank == 0) call MPI_FILE_DELETE(filename, MPI_INFO_NULL, ierror)
    call MPI_FILE_OPEN(comm, trim(adjustl(filename)), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), &
         MPI_INFO_NULL, iunit, ierror)

    ! Write header (root process only).
    if (procRank == 0) then
       buffer = trim(adjustl(name))
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       buffer = 'part'
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       i = 1
       call MPI_FILE_WRITE(iunit, i, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
       buffer = 'block'
       call MPI_FILE_WRITE(iunit, buffer, 80, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
    end if

    ! Write the file.
    offset = 3 * 80 + 4
    call MPI_FILE_SET_VIEW(iunit, offset, MPI_REAL, this%fileview, "native",                 &
         MPI_INFO_NULL, ierror)
    call MPI_FILE_WRITE_ALL(iunit, this%buffer1_sp, this%datasize, MPI_REAL,                 &
         MPI_STATUS_IGNORE, ierror)
  
    ! Close the file.
    call MPI_FILE_CLOSE(iunit, ierror)

    write(message, '(A)') " done!"
    call writeAndFlush(comm, output_unit, message)

  end subroutine writeEnsightScalar

  subroutine writeEnsightVector(this, name, comm)

    ! <<< External modules >>>
    use MPI
    use, intrinsic :: iso_fortran_env, only : output_unit

    ! <<< Internal modules >>>
    use InputHelper, only : getFreeUnit
    use ErrorHandler, only : writeAndFlush

    ! <<< Derived types >>>
    use EnsightHelper, only : t_Ensight

    implicit none

    ! <<< Arguments >>>
    class(t_Ensight) :: this
    integer, intent(in) :: comm
    character(len = STRING_LENGTH), intent(in) :: name

    ! <<< Local variables >>>
    integer :: i, procRank, size, iunit, ierror
    integer(kind = MPI_OFFSET_KIND) :: offset
    character(len = 80) :: buffer, filename, message
    logical :: fileExists

    ! Get my rank in `comm`.
    call MPI_Comm_rank(comm, procRank, ierror)

    ! Generate the file.
    if (procRank == 0) then
       call Execute_Command_Line('mkdir -p ' // trim(adjustl(this%directory)) // '/' //      &
            trim(adjustl(name)))
       filename = trim(adjustl(this%directory)) // '/' // trim(adjustl(name)) // "/" //      &
            trim(adjustl(name)) // "."
       write(filename(len_trim(filename) + 1 : len_trim(filename) + 6), '(i6.6)')            &
            this%nOutputTimes
    end if

    ! Communicate the filename.
    call MPI_Bcast(filename, 80, MPI_CHARACTER, 0, comm, ierror)

    ! Output to screen.
    write(message, '(3A)') "Writing '", trim(filename), "'..."
    call writeAndFlush(comm, output_unit, message, advance = 'no')
  
    ! Open the file.
    inquire(file=filename, exist=fileExists)
    if (fileExists .and. procRank == 0) call MPI_FILE_DELETE(filename, MPI_INFO_NULL, ierror)
    call MPI_FILE_OPEN(comm, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, &
         iunit, ierror)

    ! Write header (root process only).
    if (procRank == 0) then
       buffer = trim(adjustl(name))
       size = 80
       call MPI_FILE_WRITE(iunit, buffer, size, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       buffer = 'part'
       size = 80
       call MPI_FILE_WRITE(iunit, buffer, size, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
       i = 1
       size = 1
       call MPI_FILE_WRITE(iunit, i, size, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
       buffer = 'block'
       size = 80
       call MPI_FILE_WRITE(iunit, buffer, size, MPI_CHARACTER, MPI_STATUS_IGNORE, ierror)
    end if

    ! Write the file.
    offset = 3 * 80 + 4
    do i = 1, 3
       call MPI_File_set_view(iunit, offset, MPI_REAL, this%fileview, "native",              &
            MPI_INFO_NULL, ierror)
       call MPI_File_write_all(iunit, this%buffer3_sp(:, i), this%datasize,                  &
            MPI_REAL, MPI_STATUS_IGNORE, ierror)
       offset = offset + 4 * int(this%gDataSize, MPI_OFFSET_KIND)
    end do
  
    ! Close the file.
    call MPI_FILE_CLOSE(iunit, ierror)

    write(message, '(A)') " done!"
    call writeAndFlush(comm, output_unit, message)

  end subroutine writeEnsightVector

end module EnsightHelperImpl

subroutine setupEnsight(this, grid, gridIndex, time)

  ! <<< External modules >>>
  use MPI

  ! <<< Internal modules >>>
  use InputHelper, only : getFreeUnit

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use EnsightHelper, only : t_Ensight

  ! <<< Private members >>>
  use EnsightHelperImpl

  implicit none

  ! <<< Arguments >>>
  class(t_Ensight) :: this
  class(t_Grid) :: grid
  integer, intent(in) :: gridIndex
  real(SCALAR_KIND), intent(in) :: time

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer, parameter :: maxrecs = 10000
  integer :: i, j, procRank, iunit, iostat, ierror
  logical :: fileExists
  character(len = STRING_LENGTH) :: line
  logical :: lineFound

  ! Get my rank in `comm`.
  call MPI_Comm_rank(grid%comm, procRank, ierror)

  ! Open the case file.
  this%directory = 'ensight-3D'
  write(this%filename, '(A,I2.2,A)') "magudi_grid", gridIndex, ".case"
  if (procRank == 0) then
     inquire(file = trim(this%directory)//'/'//trim(this%filename), exist = fileExists)
     if (fileExists) then
        ! Get the number of time values.
        this%nOutputTimes = 0
        lineFound = .false.
        open(unit = getFreeUnit(iunit), file = trim(this%directory) // '/' //                &
             trim(this%filename))
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
           read(iunit,'(A)', iostat = iostat) line

           if ( iostat < 0 ) exit 

           if (lineFound) then
              j = j + 1
              backspace(iunit)
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
  this%dataSize = product(grid%localSize)
  this%gDataSize = product(grid%globalSize)
  call MPI_TYPE_CREATE_SUBARRAY(3, grid%globalSize, grid%localSize, grid%offset,             &
       MPI_ORDER_FORTRAN, MPI_REAL, this%fileview, ierror)
  call MPI_TYPE_COMMIT(this%fileview, ierror)

  ! Allocate single-precision buffer arrays.
  allocate(this%buffer1_sp(this%dataSize))
  allocate(this%buffer3_sp(this%dataSize, 3))

  ! Write the geometry.
  call writeEnsightGeometry(this, grid, gridIndex)

end subroutine setupEnsight

subroutine outputEnsight(this, state, comm, gridIndex, mode, time, nSpecies)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use EnsightHelper, only : t_Ensight
  use State_mod, only : t_State

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables

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
  integer :: i, procRank, nDimensions, nSpecies_, ierror
  character(len = STRING_LENGTH) :: name

  ! Get my rank in `comm`.
  call MPI_Comm_rank(comm, procRank, ierror)

  nDimensions = size(state%velocity, 2)
  assert_key(nDimensions, (1, 2, 3))

  nSpecies_ = 0
  if (present(nSpecies)) nSpecies_ = nSpecies
  assert(nSpecies_ == size(state%massFraction, 2))

  select case(mode)
  case (FORWARD)

     ! Update the case file (root process only).
     if (procRank == 0) call writeEnsightCase(this, time, gridIndex, nSpecies_)

     ! Density.
     name = 'DENSITY'
     this%buffer1_sp = real(state%conservedVariables(:,1), 4)
     call writeEnsightScalar(this, name, comm)

     ! Velocity.
     name = 'VELOCITY'
     this%buffer3_sp = 0.0_4
     do i = 1, nDimensions
        this%buffer3_sp(:, i) = real(state%velocity(:, i), 4)
     end do
     call writeEnsightVector(this, name, comm)

     ! Temperature.
     name = 'TEMPERATURE'
     this%buffer1_sp = real(state%temperature(:,1), 4)
     call writeEnsightScalar(this, name, comm)

     ! Mass fractions.
     do i = 1, nSpecies_
        write(name, "(A,I2.2)") "MASS_FRACTION_", i
        this%buffer1_sp = real(state%massFraction(:,i), 4)
        call writeEnsightScalar(this, name, comm)
     end do

  case (ADJOINT)

     ! Step the counter back in time.
     if (procRank == 0) this%nOutputTimes = this%nOutputTimes - 1

     ! Write the adjoint variables.
     do i = 1, nDimensions + 2 + nSpecies_
        write(name, "(A,I2.2)") "ADJOINT_", i
        this%buffer1_sp = real(state%adjointVariables(:,i), 4)
        call writeEnsightScalar(this, name, comm)
     end do

  end select

end subroutine outputEnsight
