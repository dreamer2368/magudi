#include "config.h"

module Combustion_mod

  implicit none
  private

  type, public :: t_Combustion

     integer :: nReactions, H2, O2, N2
     real(SCALAR_KIND) :: Yf0, Yo0, YFs, Damkohler, stoichiometricRatio,                     &
          heatRelease, zelDovich
     SCALAR_TYPE, dimension(:), allocatable :: stoichiometricCoefficient

   contains

     procedure, public, pass :: setup => setupCombustion
     procedure, public, pass :: add => addCombustion

  end type t_Combustion

  interface

     subroutine setupCombustion(this, nSpecies, comm)

       import :: t_Combustion

       class(t_Combustion) :: this
       integer, intent(in) :: nSpecies, comm

     end subroutine setupCombustion

  end interface

  interface

     subroutine addCombustion(this, density, temperature, massFraction,                      &
          ratioOfSpecificHeats, coordinates, iblank, rightHandSide)

       import :: t_Combustion

       class(t_Combustion) :: this
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(in) :: density(:), temperature(:), massFraction(:,:), coordinates(:,:)
       integer, intent(in) :: iblank(:)
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine addCombustion

  end interface

end module Combustion_mod
