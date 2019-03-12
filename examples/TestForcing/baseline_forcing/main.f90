program main

	implicit none

    integer, parameter :: wp = selected_real_kind(15)
    real, parameter :: pi = 4.0_wp*ATAN(1.0_wp)

	integer :: i,j
	real(kind=wp) :: input(4,4), time=0.4_wp

	! print to screen
	print *, 'calling program main'

    print *, pi

	open(unit=301,file='baseline.dat',status='replace',form='unformatted',access='stream')

    do i=32,1,-1
        input = 0.1_wp*pi/time*COS( pi*(i-1)/32.0_wp )
        write(301) input
    end do

    close(301)

	! print to screen
	print *, 'program main...done.'

contains

end program
