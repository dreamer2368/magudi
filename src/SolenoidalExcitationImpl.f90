#include "config.h"

subroutine setupSolenoidalExcitation(this, comm, nModes, location, &
     speed, amplitude, mostUnstableFrequency, radius, seed)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use SolenoidalExcitation_type

  ! <<< Public members >>>
  use SolenoidalExcitation_mod, only : cleanupSolenoidalExcitation

  ! <<< Internal modules >>>
  use RandomNumber, only : initializeRandomNumberGenerator

  implicit none

  ! <<< Arguments >>>
  type(t_SolenoidalExcitation) :: this
  integer, intent(in) :: comm, nModes
  real(SCALAR_KIND), intent(in) :: location(:), speed(:), &
       amplitude, mostUnstableFrequency, radius
  integer, intent(in), optional :: seed

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND, real64 = selected_real_kind(15)
  real(wp), parameter :: pi = 4.0_real64 * atan(1.0_real64)
  integer :: i, n, ierror
  integer, allocatable :: seed_(:)

  call cleanupSolenoidalExcitation(this)

  this%nModes = nModes

  this%location = 0.0_wp
  this%location(1:size(location)) = location

  this%speed = 0.0_wp
  this%speed(1:size(speed)) = speed

  this%amplitude = amplitude
  this%mostUnstableFrequency = mostUnstableFrequency
  this%gaussianFactor = 9.0_wp / (2.0_wp * radius ** 2)

  if (present(seed)) then
     call random_seed(size = n)
     allocate(seed_(n))
     seed_ = seed
     call random_seed(put = seed_)
     SAFE_DEALLOCATE(seed_)
  else
     call initializeRandomNumberGenerator()
  end if

  allocate(this%angularFrequencies(this%nModes))
  allocate(this%phases(this%nModes, 3))
  call random_number(this%angularFrequencies)
  call random_number(this%phases)
  this%phases = 2.0_real64 * pi * this%phases

  do i = 1, nModes
     this%angularFrequencies(i) =                                                            &
          2.0_real64 * pi * real(this%mostUnstableFrequency, real64) *                       &
          (real(i, real64) + (this%angularFrequencies(i) - 0.5_real64)) /                    &
          (0.5_real64 * real(this%nModes, real64))
  end do
  call MPI_Bcast(this%angularFrequencies, size(this%angularFrequencies),                     &
       MPI_REAL8, 0, comm, ierror)
  call MPI_Bcast(this%phases, size(this%phases), MPI_REAL8, 0, comm, ierror)

end subroutine setupSolenoidalExcitation

subroutine cleanupSolenoidalExcitation(this)

  ! <<< Derived types >>>
  use SolenoidalExcitation_type

  implicit none

  ! <<< Arguments >>>
  type(t_SolenoidalExcitation) :: this

  SAFE_DEALLOCATE(this%angularFrequencies)
  SAFE_DEALLOCATE(this%phases)

end subroutine cleanupSolenoidalExcitation
