#include "config.h"

subroutine setupCombustion(this, nSpecies, comm)

  ! <<< Derived types >>>
  use Combustion_mod, only : t_Combustion

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_Combustion) :: this
  integer, intent(in) :: nSpecies, comm

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)

  if (nSpecies <= 0) return

  ! Species indices.
  this%H2 = 1
  this%O2 = 2
  this%N2 = nSpecies + 1

  ! Read combustion parameters from input.
  call getRequiredOption("number_of_reactions", this%nReactions, comm)
  call getRequiredOption("stoichiometric_ratio", this%stoichiometricRatio, comm)
  call getRequiredOption("heat_release", this%heatRelease, comm)
  call getRequiredOption("Zel_Dovich", this%zelDovich, comm)
  call getRequiredOption("Damkohler_number", this%Damkohler, comm)
  call getRequiredOption("initial_fuel_mass_fraction", this%Yf0, comm)
  call getRequiredOption("initial_oxidizer_mass_fraction", this%Yo0, comm)

  ! Stoichiometric fuel mass fraction.
  this%Yfs = 1.0_wp / (1.0_wp + this%stoichiometricRatio * this%Yf0 / this%Yo0)

  ! Stoichiometric coefficients.
  allocate(this%stoichiometricCoefficient(nSpecies))
  this%stoichiometricCoefficient = 0.0_wp
  if (this%nReactions == 1) then
     ! One-step chemistry.
     this%stoichiometricCoefficient(this%H2) = 1.0_wp
     this%stoichiometricCoefficient(this%O2) = this%stoichiometricRatio
  else if (this%nReactions == 0) then
     ! Nothing to do
  else
     print *, 'WARNING, maximum of 1 reaction for now!'
     stop
  end if

end subroutine setupCombustion

subroutine addCombustion(this, density, temperature, massFraction, ratioOfSpecificHeats,     &
     coordinates, iblank, rightHandSide)

  ! <<< Derived types >>>
  use Combustion_mod, only : t_Combustion

  implicit none

  ! <<< Arguments >>>
  class(t_Combustion) :: this
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
  SCALAR_TYPE, intent(in) :: density(:), temperature(:), massFraction(:,:), coordinates(:,:)
  integer, intent(in) :: iblank(:)
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, k, nDimensions, nSpecies
  real(SCALAR_KIND) :: chemicalSource, referenceTemperature, flameTemperature, Yfsi
  SCALAR_TYPE :: massFraction_(size(massFraction,2))

  nSpecies = size(massFraction,2)
  assert(nSpecies >= 0)

  if (nSpecies <= 0 .or. this%nReactions == 0) return

  nDimensions = size(coordinates, 2)
  assert_key(nDimensions, (1, 2, 3))
  assert(this%nReactions > 0)

  if (this%nReactions == 1) then

     ! One-step chemistry.
     ! H2 + sO2 -> (1+s)P

     yfsi = 1.0_wp / this%Yfs

     do i = 1, size(rightHandSide, 1)
        if (iblank(i) == 0) cycle

        ! Bound mass fractions between 0 and 1.
        do k = 1, nSpecies
           massFraction_(k) = max(massFraction(i,k), 0.0_wp)
           massFraction_(k) = min(massFraction_(k), 1.0_wp)
        end do

        referenceTemperature = 1.0_wp / (ratioOfSpecificHeats - 1.0_wp)
        flameTemperature = referenceTemperature / (1.0_wp - this%heatRelease)

        chemicalSource = this%Damkohler * density(i) * massFraction_(this%H2) *              &
             massFraction_(this%O2) * exp(-this%zelDovich / this%heatRelease *               &
             flameTemperature / temperature(i))

        ! Heat release due to combustion.
        rightHandSide(i,nDimensions+2) = rightHandSide(i,nDimensions+2) +                    &
             this%heatRelease * flameTemperature * chemicalSource * Yfsi

        ! Species source terms.
        do k = 1, nSpecies
           rightHandSide(i,nDimensions+2+k) = rightHandSide(i,nDimensions+2+k) -             &
                this%stoichiometricCoefficient(k) * chemicalSource
        end do
     end do

  end if

end subroutine addCombustion
