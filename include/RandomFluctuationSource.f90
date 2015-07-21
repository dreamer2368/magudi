#include "config.h"

module RandomFluctuationSource_mod

  implicit none
  private

  type, public :: t_RandomFluctuationSource

     real(SCALAR_KIND) :: amplitude

   contains

     procedure, public, pass :: setup => setupRandomFluctuationSource
     procedure, public, pass :: add => addRandomFluctuationSource

  end type t_RandomFluctuationSource

  interface

     subroutine setupRandomFluctuationSource(this, amplitude)

       import :: t_RandomFluctuationSource

       class(t_RandomFluctuationSource) :: this
       real(SCALAR_KIND), intent(in) :: amplitude

     end subroutine setupRandomFluctuationSource

  end interface

  interface

     subroutine addRandomFluctuationSource(this, time, iblank, rightHandSide)

       import :: t_RandomFluctuationSource

       class(t_RandomFluctuationSource) :: this
       real(SCALAR_KIND), intent(in) :: time
       integer, intent(in) :: iblank(:)
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine addRandomFluctuationSource

  end interface

end module RandomFluctuationSource_mod
