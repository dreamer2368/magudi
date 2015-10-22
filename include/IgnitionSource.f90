#include "config.h"

module IgnitionSource_mod

  implicit none

  type, public :: t_IgnitionSource

     real(SCALAR_KIND) :: location(3), radius(3), amplitude, timeStart, timeDuration,        &
          shockMach
     logical :: depositVorticity

   contains

     procedure, public, pass :: setup => setupIgnitionSource
     procedure, public, pass :: add => addIgnitionSource
     procedure, public, pass :: addAdjoint => addAdjointIgnitionSource

  end type t_IgnitionSource

  interface

     subroutine setupIgnitionSource(this, ratioOfSpecificHeats, location, radius, amplitude, &
          timeStart, timeDuration, shockMach)

       use Grid_mod, only : t_Grid

       import :: t_IgnitionSource

       class(t_IgnitionSource) :: this
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats, location(:), radius(:),        &
            amplitude, timeStart, timeDuration, shockMach

     end subroutine setupIgnitionSource

  end interface

  interface

     subroutine addIgnitionSource(this, time, coordinates, iblank, density,                  &
          ratioOfSpecificHeats, heatRelease, rightHandSide)

       import :: t_IgnitionSource

       class(t_IgnitionSource) :: this
       real(SCALAR_KIND), intent(in) :: time, ratioOfSpecificHeats, heatRelease
       SCALAR_TYPE, intent(in) :: coordinates(:,:), density(:)
       integer, intent(in) :: iblank(:)
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine addIgnitionSource

  end interface

  interface

     subroutine addAdjointIgnitionSource(this, time, coordinates, iblank, adjointVariables,  &
          ratioOfSpecificHeats, rightHandSide)

       import :: t_IgnitionSource

       class(t_IgnitionSource) :: this
       real(SCALAR_KIND), intent(in) :: time, ratioOfSpecificHeats
       SCALAR_TYPE, intent(in) :: coordinates(:,:), adjointVariables(:,:)
       integer, intent(in) :: iblank(:)
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine addAdjointIgnitionSource

  end interface

end module IgnitionSource_mod
