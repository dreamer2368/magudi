#include "config.h"

module AcousticSource_mod

  implicit none
  private

  type, public :: AcousticSource

     real(SCALAR_KIND) :: location(3), amplitude, gaussianFactor, angularFrequency, phase

   contains

     procedure, public, pass :: setup => setupAcousticSource
     procedure, public, pass :: add => addAcousticSource

  end type AcousticSource

  interface
     subroutine setupAcousticSource(this, location, amplitude, frequency, radius, phase)
       import :: AcousticSource
       class(AcousticSource) :: this
       real(SCALAR_KIND), intent(in) :: location(:), amplitude, frequency, radius
       real(SCALAR_KIND), intent(in), optional :: phase
     end subroutine setupAcousticSource
  end interface

  interface
     subroutine addAcousticSource(this, time, coordinates, iblank, rightHandSide)
       import :: AcousticSource
       class(AcousticSource) :: this
       real(SCALAR_KIND), intent(in) :: time
       SCALAR_TYPE, intent(in) :: coordinates(:,:)
       integer, intent(in) :: iblank(:)
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)
     end subroutine addAcousticSource
  end interface

end module AcousticSource_mod
