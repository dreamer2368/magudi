#include "config.h"

subroutine setupCombustion(this, nSpecies, species, comm)

  ! <<< Derived types >>>
  use Combustion_mod, only : t_Combustion

  ! <<< Internal modules >>>
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  ! <<< Enumerations >>>
  use Combustion_enum

  implicit none

  ! <<< Arguments >>>
  class(t_Combustion) :: this
  integer, intent(in) :: nSpecies, comm
  character(len = STRING_LENGTH), intent(in), optional :: species(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp), parameter :: pi = 4.0_wp * atan(1.0_wp)
  integer :: k
  character(len = STRING_LENGTH) :: message, val

  if (nSpecies == 0) return

  ! Get species indices.
  if (present(species)) then
     do k = 1, nSpecies + 1
        select case (trim(species(k)))
        case ('H2', 'HYDROGEN')
           this%H2 = k
        case ('O2', 'OXYGEN')
           this%O2 = k
        case ('N2', 'NITROGEN')
           this%N2 = k
        case default
           write(message, '(3A)') "Unknown species: ", trim(species(k)), "!"
           call gracefulExit(comm, message)
        end select
     end do
  else
     this%H2 = 1
     this%O2 = 2
     this%N2 = nSpecies + 1
  end if

  ! Combustion model
  val = getOption("combustion_model", "NONE")
  select case (trim(val))
  case ("NONE")
     ! No combustion.
     this%chemistryModel = NONE

  case ("ONE_STEP")
     ! One-step irreversible reaction
     this%chemistryModel = ONE_STEP

     this%nReactions = 1

     allocate(this%Y0(nSpecies))

     ! Read combustion parameters from input.
     call getRequiredOption("stoichiometric_ratio", this%stoichiometricRatio, comm)
     call getRequiredOption("heat_release", this%heatRelease, comm)
     call getRequiredOption("Zel_Dovich", this%zelDovich, comm)
     call getRequiredOption("Damkohler_number", this%Damkohler, comm)
     call getRequiredOption("initial_fuel_mass_fraction", this%Y0(this%H2), comm)
     call getRequiredOption("initial_oxidizer_mass_fraction", this%Y0(this%O2), comm)

     ! Stoichiometric fuel mass fraction.
     this%Yfs = 1.0_wp / (1.0_wp + this%stoichiometricRatio * this%Y0(this%H2)               &
          / this%Y0(this%O2))

     ! Stoichiometric coefficients.
     allocate(this%stoichiometricCoefficient(nSpecies))
     this%stoichiometricCoefficient = 0.0_wp

     this%stoichiometricCoefficient(this%H2) = 1.0_wp
     this%stoichiometricCoefficient(this%O2) = this%stoichiometricRatio

  case ("DETAILED")
     ! Detailed chemistry.
     this%chemistryModel = DETAILED

     call getRequiredOption("number_of_reactions", this%nReactions, comm)

     do k = 1, this%nSpecies
        write(message, "(A,I1.1)") "stoichiometric_coefficient_", k
        call getRequiredOption(trim(message), this%stoichiometricCoefficient(k), comm)
     end do

  case default
     write(message, '(3A)') "Unknown combustion model ",trim(val), "!"
     call gracefulExit(comm, message)

  end select

  ! Well-stirred reactor assumption
  this%wellStirredReactor = getOption("well_stirred_reactor", .false.)
  if (this%wellStirredReactor) then
     allocate(this%Yin(nSpecies))
     call getRequiredOption("residence_time", this%residenceTime, comm)
     this%residenceTime = 1.0_wp / this%residenceTime
     call getRequiredOption("inlet_temperature", this%Tin, comm)
     call getRequiredOption("inlet_fuel", this%Yin(this%H2), comm)
     call getRequiredOption("inlet_oxidizer", this%Yin(this%O2), comm)
  end if

end subroutine setupCombustion

subroutine addCombustionForward(this, nDimensions, nSpecies, ratioOfSpecificHeats,           &
     conservedVariables, temperature, massFraction, iblank, rightHandSide)

  ! <<< Derived types >>>
  use Combustion_mod, only : t_Combustion

  implicit none

  ! <<< Arguments >>>
  class(t_Combustion) :: this
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
  SCALAR_TYPE, intent(in) :: conservedVariables(:,:), temperature(:), massFraction(:,:)
  integer, intent(in) :: nDimensions, nSpecies, iblank(:)
  SCALAR_TYPE, intent(inout) :: rightHandSide(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, k
  real(SCALAR_KIND) :: referenceTemperature, flameTemperature, activationTemperature,        &
       chemicalSource, H

  assert(nSpecies >= 0)
  assert(size(massFraction,2) == nSpecies)

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

        chemicalSource = this%Damkohler * conservedVariables(i,nDimensions+2+this%H2) *      &
             conservedVariables(i,nDimensions+2+this%O2) *                                   &
             exp(- activationTemperature / temperature(i))

        ! Heat release due to combustion.
        rightHandSide(i,nDimensions+2) = rightHandSide(i,nDimensions+2) +                    &
             H * chemicalSource

        ! Species source terms.
        do k = 1, nSpecies
           rightHandSide(i,nDimensions+2+k) = rightHandSide(i,nDimensions+2+k) -             &
                this%stoichiometricCoefficient(k) * chemicalSource
        end do

        ! Well-stirred reactor
        if (this%wellStirredReactor) then
           rightHandSide(i,nDimensions+2) = rightHandSide(i,nDimensions+2) +                 &
                this%residenceTime * (this%Tin / ratioOfSpecificHeats -                      &
                conservedVariables(i,nDimensions+2))
           do k = 1, nSpecies
              rightHandSide(i,nDimensions+2+k) = rightHandSide(i,nDimensions+2+k) +          &
                   this%residenceTime * (this%Yin(k) - conservedVariables(i,nDimensions+2+k))
           end do
        end if
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
