#include "config.h"

subroutine setupCombustion(this, nSpecies, comm)

  ! <<< Derived types >>>
  use Combustion_mod, only : t_Combustion

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_Combustion) :: this
  integer, intent(in) :: nSpecies, comm

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  character(len = STRING_LENGTH) :: message

  if (nSpecies == 0) return

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

     ! Nothing to do.

  else

     write(message, '(A)') "WARNING, maximum of 1 reaction for now!"
     call gracefulExit(comm, message)

  end if

end subroutine setupCombustion

subroutine addCombustionForward(this, nDimensions, density, temperature, massFraction,       &
     ratioOfSpecificHeats, iblank, rightHandSide)

  ! <<< Derived types >>>
  use Combustion_mod, only : t_Combustion

  implicit none

  ! <<< Arguments >>>
  class(t_Combustion) :: this
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
  SCALAR_TYPE, intent(in) :: density(:), temperature(:), massFraction(:,:)
  integer, intent(in) :: nDimensions, iblank(:)
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, k, nSpecies
  real(SCALAR_KIND) :: referenceTemperature, flameTemperature, activationTemperature,        &
       chemicalSource, H
  SCALAR_TYPE :: massFraction_(size(massFraction,2))

  nSpecies = size(massFraction,2)
  assert(nSpecies >= 0)

  if (nSpecies == 0 .or. this%nReactions == 0) return

  assert_key(nDimensions, (1, 2, 3))
  assert(nSpecies > 0)
  assert(this%nReactions > 0)

  if (this%nReactions == 1) then

     ! One-step chemistry.
     ! H2 + sO2 => (1+s)P

     referenceTemperature = 1.0_wp / (ratioOfSpecificHeats - 1.0_wp)
     flameTemperature = referenceTemperature / (1.0_wp - this%heatRelease)
     activationTemperature = this%zelDovich / this%heatRelease * flameTemperature
     H = this%heatRelease * flameTemperature / this%Yfs

     do i = 1, size(rightHandSide, 1)
        if (iblank(i) == 0) cycle

        ! Bound mass fractions between 0 and 1.
        do k = 1, nSpecies
           massFraction_(k) = massFraction(i,k)
           !massFraction_(k) = max(massFraction(i,k), 0.0_wp)
           !massFraction_(k) = min(massFraction_(k), 1.0_wp)
        end do

        chemicalSource = this%Damkohler * density(i) * massFraction_(this%H2) *              &
             massFraction_(this%O2) * exp(- activationTemperature / temperature(i))

        ! Heat release due to combustion.
        rightHandSide(i,nDimensions+2) = rightHandSide(i,nDimensions+2) +                    &
             H * chemicalSource

        ! Species source terms.
        do k = 1, nSpecies
           rightHandSide(i,nDimensions+2+k) = rightHandSide(i,nDimensions+2+k) -             &
                this%stoichiometricCoefficient(k) * chemicalSource
        end do
     end do

  end if

end subroutine addCombustionForward

subroutine addCombustionAdjoint(this, nDimensions, nSpecies, nUnknowns,                      &
     ratioOfSpecificHeats, conservedVariables, adjointVariables, velocity, massFraction,     &
     specificVolume, temperature, iblank, rightHandSide)

  ! <<< Derived types >>>
  use Combustion_mod, only : t_Combustion

  ! <<< Internal modules >>>
  use CNSHelper, only : computeJacobianOfSource

  implicit none

  ! <<< Arguments >>>
  class(t_Combustion) :: this
  integer, intent(in) :: nDimensions, nSpecies, nUnknowns, iblank(:)
  SCALAR_TYPE, dimension(:,:), intent(in) :: conservedVariables, adjointVariables,           &
       velocity, massFraction, specificVolume, temperature
  SCALAR_TYPE, intent(in) :: ratioOfSpecificHeats
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nGridPoints
  SCALAR_TYPE, allocatable :: localSourceJacobian(:,:), temp1(:), temp2(:),                  &
       localConservedVariables(:), localVelocity(:), localMassFraction(:)

  if (nSpecies == 0 .or. this%nReactions == 0) return

  assert_key(nDimensions, (1, 2, 3))
  assert(nSpecies > 0)
  assert(nUnknowns == nDimensions + 2 + nSpecies)
  assert(this%nReactions > 0)
  nGridPoints = size(conservedVariables,1)
  assert(nGridPoints > 0)

  allocate(localConservedVariables(nUnknowns))
  allocate(localVelocity(nDimensions))
  allocate(localMassFraction(nSpecies))
  allocate(localSourceJacobian(nUnknowns, nUnknowns))
  allocate(temp1(nUnknowns))
  allocate(temp2(nUnknowns))

  do j = 1, nGridPoints

     if (iblank(j) == 0) cycle

     localConservedVariables = conservedVariables(j,:)
     localVelocity = velocity(j,:)
     localMassFraction = massFraction(j,:)

     call computeJacobianOfSource(nDimensions, nSpecies,                                     &
          localConservedVariables, ratioOfSpecificHeats, this,                               &
          localSourceJacobian, specificVolume = specificVolume(j,1),                         &
          velocity = localVelocity, temperature = temperature(j,1),                          &
          massFraction = localMassFraction)

     temp1 = matmul(transpose(localSourceJacobian), adjointVariables(j,:))

     temp2(1) = temp1(1)
     do i = 1, nDimensions
        temp2(1) = temp2(1) - localVelocity(i) * specificVolume(j,1) *                       &
             temp1(i+1)
     end do
     temp2(1) = temp2(1) + (0.5_wp * ratioOfSpecificHeats *                                  &
          sum(localVelocity ** 2) - temperature(j,1)) * specificVolume(j,1) *                &
          temp1(nDimensions+2)
     do k = 1, nSpecies
        temp2(1) = temp2(1) - localMassFraction(k) * specificVolume(j,1) *                   &
             temp1(nDimensions+2+k)
     end do

     do i = 1, nDimensions
        temp2(i+1) = specificVolume(j,1) * temp1(i+1) - ratioOfSpecificHeats *               &
             localVelocity(i) * specificVolume(j,1) * temp1(nDimensions+2)
     end do

     temp2(nDimensions+2) = ratioOfSpecificHeats * specificVolume(j,1) *                     &
          temp1(nDimensions+2)

     do k = 1, nSpecies
        temp2(nDimensions+2+k) = specificVolume(j,1) * temp1(nDimensions+2+k)
     end do

     rightHandSide(j,:) = rightHandSide(j,:) - temp2

  end do !... j = 1, nGridPoints

  SAFE_DEALLOCATE(localConservedVariables)
  SAFE_DEALLOCATE(localVelocity)
  SAFE_DEALLOCATE(localMassFraction)
  SAFE_DEALLOCATE(localSourceJacobian)
  SAFE_DEALLOCATE(temp1)
  SAFE_DEALLOCATE(temp2)

end subroutine addCombustionAdjoint
