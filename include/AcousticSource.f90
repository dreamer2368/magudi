#include "config.h"

module AcousticSource_type

  implicit none

  type, public :: t_AcousticSource
     real(SCALAR_KIND) :: location(3), amplitude, gaussianFactor, angularFrequency, phase
  end type t_AcousticSource

end module AcousticSource_type

module AcousticSource_mod

  implicit none

  interface

     subroutine setupAcousticSource(this, amplitude, frequency, radius, x, y, z, phase)

       use AcousticSource_type

       type(t_AcousticSource) :: this
       real(SCALAR_KIND), intent(in) :: amplitude, frequency, radius

       real(SCALAR_KIND), intent(in), optional :: x, y, z, phase

     end subroutine setupAcousticSource

  end interface

  interface

     subroutine addAcousticSource(this, time, location, iblank, rightHandSide)

       use AcousticSource_type

       type(t_AcousticSource) :: this
       real(SCALAR_KIND), intent(in) :: time
       SCALAR_TYPE, intent(in) :: location(:,:)
       integer, intent(in) :: iblank(:)
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine addAcousticSource

  end interface

end module AcousticSource_mod
