#include "config.h"

module IgnitionSource_mod

  implicit none

  type, public :: t_IgnitionSource

     real(SCALAR_KIND) :: location(3), radius(3), amplitude, timeStart, timeDuration

   contains

     procedure, public, pass :: setup => setupIgnitionSource
     procedure, public, pass :: add => addIgnitionSource

  end type t_IgnitionSource

  interface

     subroutine setupIgnitionSource(this, location, radius, amplitude, timeStart,            &
          timeDuration)

       use Grid_mod, only : t_Grid

       import :: t_IgnitionSource

       class(t_IgnitionSource) :: this
       real(SCALAR_KIND), intent(in) :: location(:), radius(:), amplitude, timeStart,        &
            timeDuration

     end subroutine setupIgnitionSource

  end interface

  interface

     subroutine addIgnitionSource(this, time, coordinates, iblank, ratioOfSpecificHeats,     &
          heatRelease, rightHandSide)

       import :: t_IgnitionSource

       class(t_IgnitionSource) :: this
       real(SCALAR_KIND), intent(in) :: time, ratioOfSpecificHeats, heatRelease
       SCALAR_TYPE, intent(in) :: coordinates(:,:)
       integer, intent(in) :: iblank(:)
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine addIgnitionSource

  end interface

end module IgnitionSource_mod
