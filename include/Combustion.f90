#include "config.h"

module Combustion_enum

  implicit none
  public

  integer, parameter, public ::                                                              &
       NONE             = 1,                                                                 &
       ONE_STEP         = 2,                                                                 &
       DETAILED         = 3

end module Combustion_enum

module Combustion_mod

  implicit none

  type, public :: t_Combustion

     integer :: chemistryModel, nReactions, H2, O2, N2
     real(SCALAR_KIND) :: YFs, Damkohler, stoichiometricRatio, heatRelease, zelDovich,       &
          residenceTime, Tin
     SCALAR_TYPE, dimension(:), allocatable :: Y0, stoichiometricCoefficient, Yin,           &
          activationTemperature, temperatureExponent
     SCALAR_TYPE, dimension(:,:), allocatable :: stoichiometricCoefficientF,                 &
          stoichiometricCoefficientR
     logical :: wellStirredReactor

   contains

     procedure, public, pass :: setup => setupCombustion
     procedure, public, pass :: addForward => addCombustionForward
     procedure, public, pass :: addAdjoint => addCombustionAdjoint

  end type t_Combustion

  type, public :: t_Species

     real(SCALAR_KIND) :: schmidtNumberInverse, molecularWeightInverse
     SCALAR_TYPE, dimension(:), allocatable :: stoichiometricCoefficientF,                   &
          stoichiometricCoefficientR
     character(len = STRING_LENGTH) :: name

  end type t_Species

  interface

     subroutine setupCombustion(this, nSpecies, species, comm)

       import :: t_Combustion

       class(t_Combustion) :: this
       integer, intent(in) :: nSpecies, comm
       character(len = STRING_LENGTH), intent(in), optional :: species(:)

     end subroutine setupCombustion

  end interface

  interface

     subroutine addCombustionForward(this, nDimensions, nSpecies, ratioOfSpecificHeats,      &
          conservedVariables, temperature, massFraction, iblank, rightHandSide)

       import :: t_Combustion

       class(t_Combustion) :: this
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(in) :: conservedVariables(:,:), temperature(:), massFraction(:,:)
       integer, intent(in) :: nDimensions, nSpecies, iblank(:)
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine addCombustionForward

  end interface

  interface

     subroutine addCombustionAdjoint(this, nDimensions, nSpecies, nUnknowns, equationOfState,&
          ratioOfSpecificHeats, conservedVariables, adjointVariables, velocity,              &
          massFraction, specificVolume, temperature, molecularWeightInverse, iblank,         &
          rightHandSide)

       import :: t_Combustion

       class(t_Combustion) :: this
       integer, intent(in) :: nDimensions, nSpecies, nUnknowns, iblank(:), equationOfState
       SCALAR_TYPE, dimension(:,:), intent(in) :: conservedVariables, adjointVariables,      &
            velocity, massFraction, specificVolume, temperature
       SCALAR_TYPE, intent(in) :: ratioOfSpecificHeats, molecularWeightInverse(:)
       SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

     end subroutine addCombustionAdjoint

  end interface

end module Combustion_mod
