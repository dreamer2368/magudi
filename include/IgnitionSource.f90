#include "config.h"

module IgnitionSource_mod

  implicit none
  private

  type, public :: t_IgnitionSource

     real(SCALAR_KIND) :: location(3), amplitude, radius(3), timeStart, timeDuration

     SCALAR_TYPE, allocatable :: strength(:)

   contains

     procedure, public, pass :: setup => setupIgnitionSource
     procedure, public, pass :: add => addIgnitionSource

  end type t_IgnitionSource

  interface

     subroutine setupIgnitionSource(this, location, amplitude, radius, timeStart,            &
          timeDuration)

       use Grid_mod, only : t_Grid

       import :: t_IgnitionSource

       class(t_IgnitionSource) :: this
       real(SCALAR_KIND), intent(in) :: location(:), amplitude, radius, timeStart,           &
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
