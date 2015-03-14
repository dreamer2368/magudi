#include "config.h"

module AcousticSource_type

  implicit none
  private

  type, public :: t_AcousticSource
     real(SCALAR_KIND) :: location(3), amplitude, gaussianFactor, angularFrequency, phase
  end type t_AcousticSource

end module AcousticSource_type

module AcousticSource_mod

  implicit none
  public

  interface

     subroutine setupAcousticSource(this, location, amplitude, frequency, radius, phase)

       use AcousticSource_type

       type(t_AcousticSource) :: this
       real(SCALAR_KIND), intent(in) :: location(:), amplitude, frequency, radius

       real(SCALAR_KIND), intent(in), optional :: phase

     end subroutine setupAcousticSource

  end interface

  interface

     subroutine addAcousticSource(this, time, coordinates, iblank, rightHandSide)

       use AcousticSource_type

       type(t_AcousticSource) :: this
       real(SCALAR_KIND), intent(in) :: time
       SCALAR_TYPE, intent(in) :: coordinates(:,:)
       integer, intent(in) :: iblank(:)
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine addAcousticSource

  end interface

end module AcousticSource_mod
