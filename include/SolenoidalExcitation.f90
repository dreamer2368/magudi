#include "config.h"

module SolenoidalExcitation_type

  implicit none
  private

  integer, parameter, private :: real64 = selected_real_kind(15)

  type, public :: t_SolenoidalExcitation

     integer :: nModes
     real(SCALAR_KIND) :: location(3), speed(3), amplitude, &
          gaussianFactor, mostUnstableFrequency

     real(real64), allocatable :: angularFrequencies(:), phases(:,:)

  end type t_SolenoidalExcitation

end module SolenoidalExcitation_type

module SolenoidalExcitation_mod

  implicit none
  public

  interface

     subroutine setupSolenoidalExcitation(this, comm, nModes, location, &
          speed, amplitude, mostUnstableFrequency, radius, seed)

       use SolenoidalExcitation_type

       type(t_SolenoidalExcitation) :: this
       integer, intent(in) :: comm, nModes
       real(SCALAR_KIND), intent(in) :: location(:), speed(:), &
            amplitude, mostUnstableFrequency, radius

       integer, intent(in), optional :: seed

     end subroutine setupSolenoidalExcitation

  end interface

  interface

     subroutine cleanupSolenoidalExcitation(this)

       use SolenoidalExcitation_type

       type(t_SolenoidalExcitation) :: this

     end subroutine cleanupSolenoidalExcitation

  end interface

  interface

     subroutine addSolenoidalExcitation(this, time, coordinates, iblank, rightHandSide)

       use SolenoidalExcitation_type

       type(t_SolenoidalExcitation) :: this
       real(SCALAR_KIND), intent(in) :: time
       SCALAR_TYPE, intent(in) :: coordinates(:,:)
       integer, intent(in) :: iblank(:)
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine addSolenoidalExcitation

  end interface

end module SolenoidalExcitation_mod
