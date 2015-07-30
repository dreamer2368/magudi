#include "config.h"

module Combustion_mod

  implicit none

  type, public :: t_Combustion

     integer :: nReactions, H2, O2, N2
     real(SCALAR_KIND) :: Yf0, Yo0, YFs, Damkohler, stoichiometricRatio,                     &
          heatRelease, zelDovich
     SCALAR_TYPE, dimension(:), allocatable :: stoichiometricCoefficient

   contains

     procedure, public, pass :: setup => setupCombustion
     procedure, public, pass :: addForward => addCombustionForward
     procedure, public, pass :: addAdjoint => addCombustionAdjoint

  end type t_Combustion

  interface

     subroutine setupCombustion(this, nSpecies, comm)

       import :: t_Combustion

       class(t_Combustion) :: this
       integer, intent(in) :: nSpecies, comm

     end subroutine setupCombustion

  end interface

  interface

     subroutine addCombustionForward(this, nDimensions, density, temperature, massFraction,  &
          ratioOfSpecificHeats, iblank, rightHandSide)

       import :: t_Combustion

       class(t_Combustion) :: this
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(in) :: density(:), temperature(:), massFraction(:,:)
       integer, intent(in) :: nDimensions, iblank(:)
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine addCombustionForward

  end interface

  interface

     subroutine addCombustionAdjoint(this, nDimensions, nSpecies, nUnknowns,                 &
          ratioOfSpecificHeats, conservedVariables, adjointVariables, velocity,              &
          massFraction, specificVolume, temperature, iblank, rightHandSide)

       import :: t_Combustion

       class(t_Combustion) :: this
       integer, intent(in) :: nDimensions, nSpecies, nUnknowns, iblank(:)
       SCALAR_TYPE, dimension(:,:), intent(in) :: conservedVariables, adjointVariables,      &
            velocity, massFraction, specificVolume, temperature
       SCALAR_TYPE, intent(in) :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine addCombustionAdjoint

  end interface

end module Combustion_mod
