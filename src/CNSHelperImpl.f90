#include "config.h"

PURE_SUBROUTINE computeDependentVariables(nDimensions, nSpecies, conservedVariables,         &
     ratioOfSpecificHeats, specificVolume,                                                   &
     velocity, pressure, temperature, massFraction)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: nDimensions, nSpecies
  SCALAR_TYPE, intent(in) :: conservedVariables(:,:)
  real(SCALAR_KIND), intent(in), optional :: ratioOfSpecificHeats
  SCALAR_TYPE, intent(out), optional :: specificVolume(:),                                   &
       velocity(:,:), pressure(:), temperature(:), massFraction(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  real(wp) :: ratioOfSpecificHeats_
  integer :: i

  assert(size(conservedVariables, 1) > 0)
  assert_key(nDimensions, (1, 2, 3))
  assert(size(conservedVariables, 2) >= nDimensions + 2)

  ratioOfSpecificHeats_ = 1.4_wp
  if (present(ratioOfSpecificHeats)) then
     assert(ratioOfSpecificHeats > 1.0_wp)
     ratioOfSpecificHeats_ = ratioOfSpecificHeats
  end if

  ! Specific volume.
  if (present(specificVolume)) then
     assert(size(specificVolume) == size(conservedVariables, 1))
     specificVolume = 1.0_wp / conservedVariables(:,1)
  end if

  ! Velocity.
  if (present(velocity)) then
     assert(size(velocity, 1) == size(conservedVariables, 1))
     assert(size(velocity, 2) == nDimensions)
     if (present(specificVolume)) then
        do i = 1, nDimensions
           velocity(:,i) = specificVolume * conservedVariables(:,i+1)
        end do
     else
        do i = 1, nDimensions
           velocity(:,i) = conservedVariables(:,i+1) / conservedVariables(:,1)
        end do
     end if
  end if

  ! Pressure.
  if (present(pressure)) then
     assert(size(pressure) == size(conservedVariables, 1))
     if (present(velocity)) then
        pressure = (ratioOfSpecificHeats_ - 1.0_wp) * (conservedVariables(:,nDimensions+2) - &
             0.5_wp * conservedVariables(:,1) * sum(velocity ** 2, dim = 2))
     else if (present(specificVolume)) then
        pressure = (ratioOfSpecificHeats_ - 1.0_wp) * (conservedVariables(:,nDimensions+2) - &
             0.5_wp * sum(conservedVariables(:,2:nDimensions+1) ** 2, dim = 2) *             &
             specificVolume)
     else
        pressure = (ratioOfSpecificHeats_ - 1.0_wp) * (conservedVariables(:,nDimensions+2) - &
             0.5_wp * sum(conservedVariables(:,2:nDimensions+1) ** 2, dim = 2) /             &
             conservedVariables(:,1))
     end if
  end if

  ! Temperature.
  if (present(temperature)) then
     assert(size(temperature) == size(conservedVariables, 1))
     if (present(pressure)) then
        if (present(specificVolume)) then
           temperature = ratioOfSpecificHeats_ *                                             &
                pressure / (ratioOfSpecificHeats_ - 1.0_wp) * specificVolume
        else
           temperature = ratioOfSpecificHeats_ *                                             &
                pressure / (ratioOfSpecificHeats_ - 1.0_wp) / conservedVariables(:,1)
        end if
     else
        temperature = ratioOfSpecificHeats_ * (conservedVariables(:,nDimensions+2) -         &
             0.5_wp * sum(conservedVariables(:,2:nDimensions+1) ** 2, dim = 2) /             &
             conservedVariables(:,1)) / conservedVariables(:,1)
     end if
  end if

  ! Mass fraction.
  if (present(massFraction) .and. nSpecies > 0) then
     assert(size(massFraction, 1) == size(conservedVariables, 1))
     assert(size(massFraction, 2) == nSpecies)
     if (present(specificVolume)) then
        do i = 1, nSpecies
           massFraction(:,i) = specificVolume * conservedVariables(:,nDimensions + 2 + i)
        end do
     else
        do i = 1, nSpecies
           massFraction(:,i) = conservedVariables(:,nDimensions + 2 + i) / conservedVariables(:,1)
        end do
     end if
  end if

end subroutine computeDependentVariables

PURE_SUBROUTINE computeTransportVariables(nSpecies, temperature, powerLawExponent,           &
     bulkViscosityRatio, ratioOfSpecificHeats, reynoldsNumberInverse, prandtlNumberInverse,  &
     schmidtNumberInverse, dynamicViscosity, secondCoefficientOfViscosity,                   &
     thermalDiffusivity, massDiffusivity)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: nSpecies
  SCALAR_TYPE, intent(in) :: temperature(:)
  real(SCALAR_KIND), intent(in) :: powerLawExponent, ratioOfSpecificHeats,                   &
       reynoldsNumberInverse
  real(SCALAR_KIND), intent(in), optional :: bulkViscosityRatio, prandtlNumberInverse,       &
       schmidtNumberInverse(:)
  SCALAR_TYPE, intent(out), optional :: dynamicViscosity(:),                                 &
       secondCoefficientOfViscosity(:), thermalDiffusivity(:), massDiffusivity(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: k

  assert(size(temperature) > 0)
  assert(powerLawExponent >= 0.0_wp)
  assert(reynoldsNumberInverse > 0.0_wp)

  if (powerLawExponent <= 0.0_wp) then !... handle powerLawExponent = 0 separately.

     ! Dynamic viscosity.
     if (present(dynamicViscosity)) then
        assert(size(dynamicViscosity) == size(temperature))
        dynamicViscosity = reynoldsNumberInverse
     end if

     ! Second coefficient of viscosity.
     if (present(secondCoefficientOfViscosity)) then
        assert(size(secondCoefficientOfViscosity) == size(temperature))
        assert(present(bulkViscosityRatio))
        assert(bulkViscosityRatio >= 0.0_wp)
        secondCoefficientOfViscosity =                                                       &
             (bulkViscosityRatio - 2.0_wp / 3.0_wp) * reynoldsNumberInverse
     end if

     ! Thermal diffusivity.
     if (present(thermalDiffusivity)) then
        assert(size(thermalDiffusivity) == size(temperature))
        assert(present(prandtlNumberInverse))
        assert(prandtlNumberInverse > 0.0_wp)
        thermalDiffusivity = reynoldsNumberInverse * prandtlNumberInverse
     end if

     ! Mass diffusivity.
     if (present(massDiffusivity) .and. nSpecies > 0) then
        assert(size(massDiffusivity,1) == size(temperature))
        assert(size(massDiffusivity,2) == nSpecies)
        assert(present(schmidtNumberInverse))
        assert(size(schmidtNumberInverse) == nSpecies)
        assert(schmidtNumberInverse > 0.0_wp)
        do k = 1, nSpecies
           massDiffusivity(:,k) = reynoldsNumberInverse * schmidtNumberInverse(k)
        end do
     end if

  else

     ! Dynamic viscosity.
     if (present(dynamicViscosity)) then
        assert(size(dynamicViscosity) == size(temperature))
        dynamicViscosity =                                                                   &
             ((ratioOfSpecificHeats - 1.0_wp) * temperature) ** powerLawExponent *           &
          reynoldsNumberInverse
     end if

     ! Second coefficient of viscosity.
     if (present(secondCoefficientOfViscosity)) then
        assert(size(secondCoefficientOfViscosity) == size(temperature))
        assert(present(bulkViscosityRatio))
        assert(bulkViscosityRatio > 0.0_wp)
        if (present(dynamicViscosity)) then
           secondCoefficientOfViscosity =                                                    &
                (bulkViscosityRatio - 2.0_wp / 3.0_wp) * dynamicViscosity
        else
           secondCoefficientOfViscosity = (bulkViscosityRatio - 2.0_wp / 3.0_wp) *           &
                ((ratioOfSpecificHeats - 1.0_wp) * temperature) ** powerLawExponent *        &
                reynoldsNumberInverse
        end if
     end if

     ! Thermal diffusivity.
     if (present(thermalDiffusivity)) then
        assert(size(thermalDiffusivity) == size(temperature))
        assert(present(prandtlNumberInverse))
        assert(prandtlNumberInverse > 0.0_wp)
        if (present(dynamicViscosity)) then
           thermalDiffusivity = dynamicViscosity * prandtlNumberInverse
        else
           thermalDiffusivity =                                                              &
                ((ratioOfSpecificHeats - 1.0_wp) * temperature) ** powerLawExponent *        &
                reynoldsNumberInverse * prandtlNumberInverse
        end if
     end if

     ! Mass diffusivity.
     if (present(massDiffusivity) .and. nSpecies > 0) then
        assert(size(massDiffusivity,1) == size(temperature))
        assert(size(massDiffusivity,2) == nSpecies)
        assert(present(schmidtNumberInverse))
        assert(size(schmidtNumberInverse) == nSpecies)
        assert(schmidtNumberInverse > 0.0_wp)
        if (present(dynamicViscosity)) then
           do k = 1, nSpecies
              massDiffusivity(:,k) = dynamicViscosity * schmidtNumberInverse(k)
           end do
        else
           do k = 1, nSpecies
              massDiffusivity(:,k) =                                                         &
                   ((ratioOfSpecificHeats - 1.0_wp) * temperature) ** powerLawExponent *     &
                   reynoldsNumberInverse * schmidtNumberInverse(k)
           end do
        end if
     end if

  end if

end subroutine computeTransportVariables

PURE_SUBROUTINE computeRoeAverage(nDimensions, conservedVariablesL,                          &
     conservedVariablesR, ratioOfSpecificHeats, roeAverage)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: nDimensions
  SCALAR_TYPE, intent(in) :: conservedVariablesL(:), conservedVariablesR(:)
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
  SCALAR_TYPE, intent(out) :: roeAverage(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i
  SCALAR_TYPE :: sqrtDensityL, sqrtDensityR, specificVolumeL, specificVolumeR,               &
       enthalpyL, enthalpyR

  sqrtDensityL = sqrt(conservedVariablesL(1))
  sqrtDensityR = sqrt(conservedVariablesR(1))

  specificVolumeL = 1.0_wp / conservedVariablesL(1)
  specificVolumeR = 1.0_wp / conservedVariablesR(1)

  enthalpyL = ratioOfSpecificHeats * conservedVariablesL(nDimensions+2) -                    &
       0.5_wp * (ratioOfSpecificHeats - 1.0_wp) * specificVolumeL *                          &
       sum(conservedVariablesL(2:nDimensions+1) ** 2)
  enthalpyR = ratioOfSpecificHeats * conservedVariablesR(nDimensions+2) -                    &
       0.5_wp * (ratioOfSpecificHeats - 1.0_wp) * specificVolumeR *                          &
       sum(conservedVariablesR(2:nDimensions+1) ** 2)

  roeAverage(1) = sqrtDensityL * sqrtDensityR
  do i = 1, nDimensions
     roeAverage(i+1) = (sqrtDensityR * conservedVariablesL(i+1) +                            &
          sqrtDensityL * conservedVariablesR(i+1)) / (sqrtDensityL + sqrtDensityR)
  end do

  roeAverage(nDimensions+2) = ((sqrtDensityR * enthalpyL + sqrtDensityL * enthalpyR) /       &
       (sqrtDensityL + sqrtDensityR) +                                                       &
       0.5_wp * (ratioOfSpecificHeats - 1.0_wp) / roeAverage(1) *                            &
       sum(roeAverage(2:nDimensions+1) ** 2)) / ratioOfSpecificHeats

end subroutine computeRoeAverage

PURE_SUBROUTINE computeStressTensor(nDimensions, velocityGradient, dynamicViscosity,         &
     secondCoefficientOfViscosity, stressTensor)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: nDimensions
  SCALAR_TYPE, intent(inout) :: velocityGradient(:,:)
  SCALAR_TYPE, intent(in) :: dynamicViscosity(:), secondCoefficientOfViscosity(:)
  SCALAR_TYPE, intent(out), optional :: stressTensor(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  SCALAR_TYPE, allocatable :: divergenceOfVelocity(:)

  assert(size(velocityGradient, 1) > 0)
  assert_key(nDimensions, (1, 2, 3))
  assert(size(velocityGradient, 2) == nDimensions ** 2)
  assert(size(dynamicViscosity) == size(velocityGradient, 1))
  assert(size(secondCoefficientOfViscosity) == size(velocityGradient, 1))

  if (present(stressTensor)) then !... out-of-place computation.

     assert(size(stressTensor, 1) == size(velocityGradient, 1))
     assert(size(stressTensor, 2) == nDimensions ** 2)

     select case (nDimensions)

     case (1)
        stressTensor(:,1) = (2.0_wp * dynamicViscosity +                                     &
             secondCoefficientOfViscosity) * velocityGradient(:,1)

     case (2)
        stressTensor(:,1) = 2.0_wp * dynamicViscosity * velocityGradient(:,1) +              &
             secondCoefficientOfViscosity * (velocityGradient(:,1) + velocityGradient(:,4))
        stressTensor(:,2) = dynamicViscosity * (velocityGradient(:,2) + velocityGradient(:,3))
        stressTensor(:,3) = stressTensor(:,2)
        stressTensor(:,4) = 2.0_wp * dynamicViscosity * velocityGradient(:,4) +              &
             secondCoefficientOfViscosity * (velocityGradient(:,1) + velocityGradient(:,4))

     case (3)
        stressTensor(:,1) = 2.0_wp * dynamicViscosity * velocityGradient(:,1) +              &
             secondCoefficientOfViscosity * (velocityGradient(:,1) +                         &
             velocityGradient(:,5) + velocityGradient(:,9))
        stressTensor(:,2) = dynamicViscosity * (velocityGradient(:,2) + velocityGradient(:,4))
        stressTensor(:,3) = dynamicViscosity * (velocityGradient(:,3) + velocityGradient(:,7))
        stressTensor(:,4) = stressTensor(:,2)
        stressTensor(:,5) = 2.0_wp * dynamicViscosity * velocityGradient(:,5) +              &
             secondCoefficientOfViscosity * (velocityGradient(:,1) +                         &
             velocityGradient(:,5) + velocityGradient(:,9))
        stressTensor(:,6) = dynamicViscosity * (velocityGradient(:,6) + velocityGradient(:,8))
        stressTensor(:,7) = stressTensor(:,3)
        stressTensor(:,8) = stressTensor(:,6)
        stressTensor(:,9) = 2.0_wp * dynamicViscosity * velocityGradient(:,9) +              &
             secondCoefficientOfViscosity * (velocityGradient(:,1) +                         &
             velocityGradient(:,5) + velocityGradient(:,9))

     end select !... nDimensions

  else !... in-place computation.

     select case (nDimensions)

     case (1)
        velocityGradient(:,1) = (2.0_wp * dynamicViscosity +                                 &
             secondCoefficientOfViscosity) * velocityGradient(:,1)

     case (2)
        allocate(divergenceOfVelocity(size(velocityGradient, 1)))
        divergenceOfVelocity = secondCoefficientOfViscosity *                                &
             (velocityGradient(:,1) + velocityGradient(:,4))
        velocityGradient(:,1) = 2.0_wp * dynamicViscosity * velocityGradient(:,1) +          &
             divergenceOfVelocity
        velocityGradient(:,2) = dynamicViscosity *                                           &
             (velocityGradient(:,2) + velocityGradient(:,3))
        velocityGradient(:,3) = velocityGradient(:,2)
        velocityGradient(:,4) = 2.0_wp * dynamicViscosity * velocityGradient(:,4) +          &
             divergenceOfVelocity

     case (3)
        allocate(divergenceOfVelocity(size(velocityGradient, 1)))
        divergenceOfVelocity = secondCoefficientOfViscosity *                                &
             (velocityGradient(:,1) + velocityGradient(:,5) + velocityGradient(:,9))
        velocityGradient(:,1) = 2.0_wp * dynamicViscosity * velocityGradient(:,1) +          &
             divergenceOfVelocity
        velocityGradient(:,2) = dynamicViscosity *                                           &
             (velocityGradient(:,2) + velocityGradient(:,4))
        velocityGradient(:,3) = dynamicViscosity *                                           &
             (velocityGradient(:,3) + velocityGradient(:,7))
        velocityGradient(:,4) = velocityGradient(:,2)
        velocityGradient(:,5) = 2.0_wp * dynamicViscosity * velocityGradient(:,5) +          &
             divergenceOfVelocity
        velocityGradient(:,6) = dynamicViscosity *                                           &
             (velocityGradient(:,6) + velocityGradient(:,8))
        velocityGradient(:,7) = velocityGradient(:,3)
        velocityGradient(:,8) = velocityGradient(:,6)
        velocityGradient(:,9) = 2.0_wp * dynamicViscosity * velocityGradient(:,9) +          &
             divergenceOfVelocity

     end select !... nDimensions

  end if

  SAFE_DEALLOCATE(divergenceOfVelocity)

end subroutine computeStressTensor

PURE_SUBROUTINE computeVorticityMagnitudeAndDilatation(nDimensions, velocityGradient,        &
     vorticityMagnitude, dilatation)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: nDimensions
  SCALAR_TYPE, intent(in) :: velocityGradient(:,:)
  SCALAR_TYPE, intent(out), optional :: vorticityMagnitude(:), dilatation(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  assert(size(velocityGradient, 1) > 0)
  assert_key(nDimensions, (1, 2, 3))
  assert(size(velocityGradient, 2) == nDimensions ** 2)

  select case (nDimensions)

  case (1)
     if (present(dilatation)) then
        assert(size(dilatation) == size(velocityGradient, 1))
        dilatation = velocityGradient(:,1)
     end if
     if (present(vorticityMagnitude)) then
        assert(size(vorticityMagnitude) == size(velocityGradient, 1))
        vorticityMagnitude = 0.0_wp
     end if

  case (2)
     if (present(dilatation)) then
        assert(size(dilatation) == size(velocityGradient, 1))
        dilatation = velocityGradient(:,1) + velocityGradient(:,4)
     end if
     if (present(vorticityMagnitude)) then
        assert(size(vorticityMagnitude) == size(velocityGradient, 1))
        vorticityMagnitude = abs(velocityGradient(:,3) - velocityGradient(:,2))
     end if

  case (3)
     if (present(dilatation)) then
        assert(size(dilatation) == size(velocityGradient, 1))
        dilatation = velocityGradient(:,1) + velocityGradient(:,5) + velocityGradient(:,9)
     end if
     if (present(vorticityMagnitude)) then
        assert(size(vorticityMagnitude) == size(velocityGradient, 1))
        vorticityMagnitude = sqrt(abs(velocityGradient(:,8) - velocityGradient(:,6)) ** 2 +  &
             abs(velocityGradient(:,3) - velocityGradient(:,7)) ** 2 +                       &
             abs(velocityGradient(:,4) - velocityGradient(:,2)) ** 2)
     end if

  end select !... nDimensions

end subroutine computeVorticityMagnitudeAndDilatation

PURE_SUBROUTINE computeCartesianInvsicidFluxes(nDimensions, nSpecies, conservedVariables,    &
     velocity, pressure, inviscidFluxes)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: nDimensions, nSpecies
  SCALAR_TYPE, intent(in) :: conservedVariables(:,:), velocity(:,:), pressure(:)
  SCALAR_TYPE, intent(out) :: inviscidFluxes(:,:,:)

  ! <<< Local variables >>>
  integer :: k

  assert(size(conservedVariables, 1) > 0)
  assert_key(nDimensions, (1, 2, 3))
  assert(size(conservedVariables, 2) >= nDimensions + 2)
  assert(size(velocity, 1) == size(conservedVariables, 1))
  assert(size(velocity, 2) == nDimensions)
  assert(size(pressure) == size(conservedVariables, 1))
  assert(size(inviscidFluxes, 1) == size(conservedVariables, 1))
  assert(size(inviscidFluxes, 2) >= nDimensions + 2)
  assert(size(inviscidFluxes, 3) == nDimensions)

  select case (nDimensions)

  case (1)
     inviscidFluxes(:,1,1) = conservedVariables(:,2)
     inviscidFluxes(:,2,1) = conservedVariables(:,2) * velocity(:,1) + pressure
     inviscidFluxes(:,3,1) = velocity(:,1) * (conservedVariables(:,3) + pressure)
     do k=1,nSpecies
        inviscidFluxes(:,k+3,1) = conservedVariables(:,k+3) * velocity(:,1)
     end do

  case (2)
     inviscidFluxes(:,1,1) = conservedVariables(:,2)
     inviscidFluxes(:,2,1) = conservedVariables(:,2) * velocity(:,1) + pressure
     inviscidFluxes(:,3,1) = conservedVariables(:,2) * velocity(:,2)
     inviscidFluxes(:,4,1) = velocity(:,1) * (conservedVariables(:,4) + pressure)
     inviscidFluxes(:,1,2) = conservedVariables(:,3)
     inviscidFluxes(:,2,2) = inviscidFluxes(:,3,1)
     inviscidFluxes(:,3,2) = conservedVariables(:,3) * velocity(:,2) + pressure
     inviscidFluxes(:,4,2) = velocity(:,2) * (conservedVariables(:,4) + pressure)
     do k=1,nSpecies
        inviscidFluxes(:,k+4,1) = conservedVariables(:,k+4) * velocity(:,1)
        inviscidFluxes(:,k+4,2) = conservedVariables(:,k+4) * velocity(:,2)
     end do

  case (3)
     inviscidFluxes(:,1,1) = conservedVariables(:,2)
     inviscidFluxes(:,2,1) = conservedVariables(:,2) * velocity(:,1) + pressure
     inviscidFluxes(:,3,1) = conservedVariables(:,2) * velocity(:,2)
     inviscidFluxes(:,4,1) = conservedVariables(:,2) * velocity(:,3)
     inviscidFluxes(:,5,1) = velocity(:,1) * (conservedVariables(:,5) + pressure)
     inviscidFluxes(:,1,2) = conservedVariables(:,3)
     inviscidFluxes(:,2,2) = inviscidFluxes(:,3,1)
     inviscidFluxes(:,3,2) = conservedVariables(:,3) * velocity(:,2) + pressure
     inviscidFluxes(:,4,2) = conservedVariables(:,3) * velocity(:,3)
     inviscidFluxes(:,5,2) = velocity(:,2) * (conservedVariables(:,5) + pressure)
     inviscidFluxes(:,1,3) = conservedVariables(:,4)
     inviscidFluxes(:,2,3) = inviscidFluxes(:,4,1)
     inviscidFluxes(:,3,3) = inviscidFluxes(:,4,2)
     inviscidFluxes(:,4,3) = conservedVariables(:,4) * velocity(:,3) + pressure
     inviscidFluxes(:,5,3) = velocity(:,3) * (conservedVariables(:,5) + pressure)
     do k=1,nSpecies
        inviscidFluxes(:,k+5,1) = conservedVariables(:,k+5) * velocity(:,1)
        inviscidFluxes(:,k+5,2) = conservedVariables(:,k+5) * velocity(:,2)
        inviscidFluxes(:,k+5,3) = conservedVariables(:,k+5) * velocity(:,3)
     end do

  end select !... nDimensions

end subroutine computeCartesianInvsicidFluxes

PURE_SUBROUTINE computeCartesianViscousFluxes(nDimensions, nSpecies, velocity,               &
     massFraction, stressTensor, heatFlux, speciesFlux, viscousFluxes)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: nDimensions, nSpecies
  SCALAR_TYPE, intent(in) :: velocity(:,:), stressTensor(:,:), heatFlux(:,:)
  SCALAR_TYPE, intent(in), optional :: massFraction(:,:), speciesFlux(:,:,:)
  SCALAR_TYPE, intent(out) :: viscousFluxes(:,:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: k

  assert(size(velocity, 1) > 0)
  assert_key(nDimensions, (1, 2, 3))
  assert(size(velocity, 2) == nDimensions)
  assert(size(stressTensor, 1) == size(velocity, 1))
  assert(size(stressTensor, 2) == nDimensions ** 2)
  assert(size(heatFlux, 1) == size(velocity, 1))
  assert(size(heatFlux, 2) == nDimensions)
  assert(size(speciesFlux, 3) == nDimensions)
  assert(size(viscousFluxes, 1) == size(velocity, 1))
  assert(size(viscousFluxes, 2) >= nDimensions + 2)
  assert(size(viscousFluxes, 3) == nDimensions)
  if (nSpecies > 0) then
     assert(size(massFraction, 1) == size(velocity, 1))
     assert(size(massFraction, 2) == nSpecies)
     assert(size(speciesFlux, 1) == size(velocity, 1))
     assert(size(speciesFlux, 2) == nSpecies)
  end if

  select case (nDimensions)

  case (1)
     viscousFluxes(:,1,1) = 0.0_wp
     viscousFluxes(:,2,1) = stressTensor(:,1)
     viscousFluxes(:,3,1) = velocity(:,1) * stressTensor(:,1) - heatFlux(:,1)
     do k=1,nSpecies
        viscousFluxes(:,k+3,1) = - speciesFlux(:,k,1)
     end do

  case (2)
     viscousFluxes(:,1,1) = 0.0_wp
     viscousFluxes(:,2,1) = stressTensor(:,1)
     viscousFluxes(:,3,1) = stressTensor(:,3)
     viscousFluxes(:,4,1) = velocity(:,1) * stressTensor(:,1) +                              &
          velocity(:,2) * stressTensor(:,3) - heatFlux(:,1)
     viscousFluxes(:,1,2) = 0.0_wp
     viscousFluxes(:,2,2) = stressTensor(:,2)
     viscousFluxes(:,3,2) = stressTensor(:,4)
     viscousFluxes(:,4,2) = velocity(:,1) * stressTensor(:,2) +                              &
          velocity(:,2) * stressTensor(:,4) - heatFlux(:,2)
     do k=1,nSpecies
        viscousFluxes(:,k+4,1) = - speciesFlux(:,k,1)
        viscousFluxes(:,k+4,2) = - speciesFlux(:,k,2)
     end do

  case (3)
     viscousFluxes(:,1,1) = 0.0_wp
     viscousFluxes(:,2,1) = stressTensor(:,1)
     viscousFluxes(:,3,1) = stressTensor(:,4)
     viscousFluxes(:,4,1) = stressTensor(:,7)
     viscousFluxes(:,5,1) = velocity(:,1) * stressTensor(:,1) +                              &
          velocity(:,2) * stressTensor(:,4) +                                                &
          velocity(:,3) * stressTensor(:,7) - heatFlux(:,1)
     viscousFluxes(:,1,2) = 0.0_wp
     viscousFluxes(:,2,2) = stressTensor(:,2)
     viscousFluxes(:,3,2) = stressTensor(:,5)
     viscousFluxes(:,4,2) = stressTensor(:,8)
     viscousFluxes(:,5,2) = velocity(:,1) * stressTensor(:,2) +                              &
          velocity(:,2) * stressTensor(:,5) +                                                &
          velocity(:,3) * stressTensor(:,8) - heatFlux(:,2)
     viscousFluxes(:,1,3) = 0.0_wp
     viscousFluxes(:,2,3) = stressTensor(:,3)
     viscousFluxes(:,3,3) = stressTensor(:,6)
     viscousFluxes(:,4,3) = stressTensor(:,9)
     viscousFluxes(:,5,3) = velocity(:,1) * stressTensor(:,3) +                              &
          velocity(:,2) * stressTensor(:,6) +                                                &
          velocity(:,3) * stressTensor(:,9) - heatFlux(:,3)
     do k=1,nSpecies
        viscousFluxes(:,k+5,1) = - speciesFlux(:,k,1)
        viscousFluxes(:,k+5,2) = - speciesFlux(:,k,2)
        viscousFluxes(:,k+5,3) = - speciesFlux(:,k,3)
     end do

  end select !... nDimensions

end subroutine computeCartesianViscousFluxes

PURE_SUBROUTINE computeSpectralRadius(nDimensions, ratioOfSpecificHeats,                     &
     velocity, temperature, metrics, spectralRadius, isDomainCurvilinear)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: nDimensions
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
  SCALAR_TYPE, intent(in) :: velocity(:,:), temperature(:), metrics(:,:)
  SCALAR_TYPE, intent(out) :: spectralRadius(:,:)
  logical, intent(in), optional :: isDomainCurvilinear

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  logical :: isDomainCurvilinear_

  assert_key(nDimensions, (1, 2, 3))
  assert(ratioOfSpecificHeats > 1.0_wp)

  assert(size(velocity, 1) > 0)
  assert(size(velocity, 2) == nDimensions)
  assert(size(temperature) == size(velocity, 1))
  assert(size(metrics, 1) == size(velocity, 1))
  assert(size(metrics, 2) == nDimensions ** 2)
  assert(size(spectralRadius, 1) == size(velocity, 1))
  assert(size(spectralRadius, 2) == nDimensions)

  isDomainCurvilinear_ = .true.
  if (present(isDomainCurvilinear)) isDomainCurvilinear_ = isDomainCurvilinear

  ! Temporary storage for speed of sound.
  spectralRadius(:,nDimensions) = sqrt((ratioOfSpecificHeats - 1.0_wp) * temperature)

  select case (nDimensions)

  case (1)
     spectralRadius(:,1) = abs(metrics(:,1) * velocity(:,1)) +                               &
          spectralRadius(:,1) * abs(metrics(:,1))

  case (2)
     if (isDomainCurvilinear_) then
        spectralRadius(:,1) = abs(metrics(:,1) * velocity(:,1) +                             &
             metrics(:,2) * velocity(:,2)) +                                                 &
             spectralRadius(:,2) * sqrt(metrics(:,1) ** 2 + metrics(:,2) ** 2)
        spectralRadius(:,2) = abs(metrics(:,3) * velocity(:,1) +                             &
             metrics(:,4) * velocity(:,2)) +                                                 &
             spectralRadius(:,2) * sqrt(metrics(:,3) ** 2 + metrics(:,4) ** 2)
     else
        spectralRadius(:,1) = abs(metrics(:,1) * velocity(:,1)) +                            &
             spectralRadius(:,2) * abs(metrics(:,1))
        spectralRadius(:,2) = abs(metrics(:,4) * velocity(:,2)) +                            &
             spectralRadius(:,2) * abs(metrics(:,4))
     end if

  case (3)
     if (isDomainCurvilinear_) then
        spectralRadius(:,1) = abs(metrics(:,1) * velocity(:,1) +                             &
             metrics(:,2) * velocity(:,2) + metrics(:,3) * velocity(:,3)) +                  &
             spectralRadius(:,3) * sqrt(metrics(:,1) ** 2 +                                  &
             metrics(:,2) ** 2 + metrics(:,3) ** 2)
        spectralRadius(:,2) = abs(metrics(:,4) * velocity(:,1) +                             &
             metrics(:,5) * velocity(:,2) + metrics(:,6) * velocity(:,3)) +                  &
             spectralRadius(:,3) * sqrt(metrics(:,4) ** 2 +                                  &
             metrics(:,5) ** 2 + metrics(:,6) ** 2)
        spectralRadius(:,3) = abs(metrics(:,7) * velocity(:,1) +                             &
             metrics(:,8) * velocity(:,2) + metrics(:,9) * velocity(:,3)) +                  &
             spectralRadius(:,3) * sqrt(metrics(:,7) ** 2 +                                  &
             metrics(:,8) ** 2 + metrics(:,9) ** 2)
     else
        spectralRadius(:,1) = abs(metrics(:,1) * velocity(:,1)) +                            &
             spectralRadius(:,3) * abs(metrics(:,1))
        spectralRadius(:,2) = abs(metrics(:,5) * velocity(:,2)) +                            &
             spectralRadius(:,3) * abs(metrics(:,5))
        spectralRadius(:,3) = abs(metrics(:,9) * velocity(:,3)) +                            &
             spectralRadius(:,3) * abs(metrics(:,9))
     end if

  end select !... nDimensions

end subroutine computeSpectralRadius

PURE_SUBROUTINE transformFluxes(nDimensions, fluxes, metrics,                                &
     transformedFluxes, isDomainCurvilinear)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: nDimensions
  SCALAR_TYPE, intent(in) :: fluxes(:,:,:), metrics(:,:)
  SCALAR_TYPE, intent(out) :: transformedFluxes(:,:,:)
  logical, intent(in), optional :: isDomainCurvilinear

  ! <<< Local variables >>>
  logical :: isDomainCurvilinear_
  integer :: i

  isDomainCurvilinear_ = .true.
  if (present(isDomainCurvilinear)) isDomainCurvilinear_ = isDomainCurvilinear

  assert(size(fluxes, 1) > 0)
  assert_key(nDimensions, (1, 2, 3))
  assert(size(fluxes, 2) >= nDimensions + 2)
  assert(size(fluxes, 3) == nDimensions)
  assert(size(metrics, 1) == size(fluxes, 1))
  assert(size(metrics, 2) == nDimensions ** 2)
  assert(all(shape(transformedFluxes) == shape(fluxes)))

  select case (nDimensions)

  case (1)
     do i = 1, size(fluxes, 2)
        transformedFluxes(:,i,1) = metrics(:,1) * fluxes(:,i,1)
     end do

  case (2)
     if (isDomainCurvilinear_) then
        do i = 1, size(fluxes, 2)
           transformedFluxes(:,i,1) = metrics(:,1) * fluxes(:,i,1) +                         &
                metrics(:,2) * fluxes(:,i,2)
           transformedFluxes(:,i,2) = metrics(:,3) * fluxes(:,i,1) +                         &
                metrics(:,4) * fluxes(:,i,2)
        end do
     else
        do i = 1, size(fluxes, 2)
           transformedFluxes(:,i,1) = metrics(:,1) * fluxes(:,i,1)
           transformedFluxes(:,i,2) = metrics(:,4) * fluxes(:,i,2)
        end do
     end if

  case (3)
     if (isDomainCurvilinear_) then
        do i = 1, size(fluxes, 2)
           transformedFluxes(:,i,1) = metrics(:,1) * fluxes(:,i,1) +                         &
                metrics(:,2) * fluxes(:,i,2) + metrics(:,3) * fluxes(:,i,3)
           transformedFluxes(:,i,2) = metrics(:,4) * fluxes(:,i,1) +                         &
                metrics(:,5) * fluxes(:,i,2) + metrics(:,6) * fluxes(:,i,3)
           transformedFluxes(:,i,3) = metrics(:,7) * fluxes(:,i,1) +                         &
                metrics(:,8) * fluxes(:,i,2) + metrics(:,9) * fluxes(:,i,3)
        end do
     else
        do i = 1, size(fluxes, 2)
           transformedFluxes(:,i,1) = metrics(:,1) * fluxes(:,i,1)
           transformedFluxes(:,i,2) = metrics(:,5) * fluxes(:,i,2)
           transformedFluxes(:,i,3) = metrics(:,9) * fluxes(:,i,3)
        end do
     end if

  end select !... nDimensions

end subroutine transformFluxes

PURE_FUNCTION computeCfl(nDimensions, iblank, jacobian, metrics, velocity, temperature,      &
     timeStepSize, ratioOfSpecificHeats, dynamicViscosity,                                   &
     thermalDiffusivity) result(cfl)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: nDimensions, iblank(:)
  SCALAR_TYPE, intent(in) :: jacobian(:), metrics(:,:), velocity(:,:), temperature(:)
  real(SCALAR_KIND), intent(in) :: timeStepSize, ratioOfSpecificHeats
  SCALAR_TYPE, intent(in), optional :: dynamicViscosity(:), thermalDiffusivity(:)

  ! <<< Result >>>
  real(SCALAR_KIND) :: cfl

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j
  real(wp) :: localSpeedOfSound, localWaveSpeed

  assert(size(iblank) > 0)
  assert_key(nDimensions, (1, 2, 3))
  assert(size(jacobian) == size(iblank))
  assert(size(metrics, 1) == size(iblank))
  assert(size(metrics, 2) == nDimensions ** 2)
  assert(size(velocity, 1) == size(iblank))
  assert(size(velocity, 2) == nDimensions)
  assert(size(temperature) == size(iblank))
  assert(timeStepSize > 0.0_wp)
  assert(ratioOfSpecificHeats > 1.0_wp)

  cfl = 0.0_wp

  ! Advection.
  do i = 1, size(iblank)
     if (iblank(i) == 0) cycle ! ... skip hole points.
     assert(real(temperature(i), wp) > 0.0_wp)
     localSpeedOfSound = sqrt((ratioOfSpecificHeats - 1.0_wp) * real(temperature(i), wp))
     localWaveSpeed = 0.0_wp
     do j = 1, nDimensions
        localWaveSpeed = localWaveSpeed + localSpeedOfSound *                                &
             sqrt(real(sum(metrics(i,1+nDimensions*(j-1):nDimensions*j) ** 2), wp)) +        &
             abs(real(dot_product(velocity(i,:),                                             &
             metrics(i,1+nDimensions*(j-1):nDimensions*j)), wp))
     end do
     localWaveSpeed = real(jacobian(i), wp) * localWaveSpeed
     cfl = max(cfl, localWaveSpeed * timeStepSize)
  end do

  ! Diffusion.
  if (present(dynamicViscosity)) then
     assert(present(thermalDiffusivity))
     assert(size(dynamicViscosity) == size(iblank))
     assert(size(thermalDiffusivity) == size(iblank))
     do i = 1, size(iblank)
        if (iblank(i) == 0) cycle !... skip hole points.
        localWaveSpeed = 0.0_wp
        do j = 1, nDimensions
           localWaveSpeed = localWaveSpeed +                                                 &
                sum(real(metrics(i,1+nDimensions*(j-1):nDimensions*j), wp) ** 2)
        end do
        localWaveSpeed = real(jacobian(i) ** 2 * sum(metrics(i,:) ** 2), wp) *               &
             max(2.0_wp * real(dynamicViscosity(i), wp), real(thermalDiffusivity(i), wp))
        cfl = max(cfl, localWaveSpeed * timeStepSize)
     end do
  end if

end function computeCfl

PURE_FUNCTION computeTimeStepSize(nDimensions, iblank, jacobian, metrics, velocity,          &
     temperature, cfl, ratioOfSpecificHeats, dynamicViscosity, thermalDiffusivity)           &
     result(timeStepSize)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: nDimensions, iblank(:)
  SCALAR_TYPE, intent(in) :: jacobian(:), metrics(:,:), velocity(:,:), temperature(:)
  real(SCALAR_KIND), intent(in) :: cfl, ratioOfSpecificHeats
  SCALAR_TYPE, intent(in), optional :: dynamicViscosity(:), thermalDiffusivity(:)

  ! <<< Result >>>
  real(SCALAR_KIND) :: timeStepSize

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j
  real(wp) :: localSpeedOfSound, localWaveSpeed


  assert(size(iblank) > 0)
  assert_key(nDimensions, (1, 2, 3))
  assert(size(jacobian) == size(iblank))
  assert(size(metrics, 1) == size(iblank))
  assert(size(metrics, 2) == nDimensions ** 2)
  assert(size(velocity, 1) == size(iblank))
  assert(size(velocity, 2) == nDimensions)
  assert(size(temperature) == size(iblank))
  assert(cfl > 0.0_wp)
  assert(ratioOfSpecificHeats > 1.0_wp)

  timeStepSize = huge(0.0_wp)

  ! Advection.
  do i = 1, size(iblank)
     if (iblank(i) == 0) cycle ! ... skip hole points.
     assert(real(temperature(i), wp) > 0.0_wp)
     localSpeedOfSound = sqrt((ratioOfSpecificHeats - 1.0_wp) * real(temperature(i), wp))
     localWaveSpeed = 0.0_wp
     do j = 1, nDimensions
        localWaveSpeed = localWaveSpeed + localSpeedOfSound *                                &
             sqrt(real(sum(metrics(i,1+nDimensions*(j-1):nDimensions*j) ** 2), wp)) +        &
             abs(real(dot_product(velocity(i,:),                                             &
             metrics(i,1+nDimensions*(j-1):nDimensions*j)), wp))
     end do
     localWaveSpeed = real(jacobian(i), wp) * localWaveSpeed
     timeStepSize = min(timeStepSize, cfl / localWaveSpeed)
  end do

  ! Diffusion.
  if (present(dynamicViscosity)) then
     assert(present(thermalDiffusivity))
     assert(size(dynamicViscosity) == size(iblank))
     assert(size(thermalDiffusivity) == size(iblank))
     do i = 1, size(iblank)
        if (iblank(i) == 0) cycle !... skip hole points.
        localWaveSpeed = 0.0_wp
        do j = 1, nDimensions
           localWaveSpeed = localWaveSpeed +                                                 &
                sum(real(metrics(i,1+nDimensions*(j-1):nDimensions*j), wp) ** 2)
        end do
        localWaveSpeed = real(jacobian(i) ** 2 * sum(metrics(i,:) ** 2), wp) *               &
             max(2.0_wp * real(dynamicViscosity(i), wp), real(thermalDiffusivity(i), wp))
        timeStepSize = min(timeStepSize, cfl / localWaveSpeed)
     end do
  end if

end function computeTimeStepSize

PURE_SUBROUTINE computeJacobianOfInviscidFlux(nDimensions, nSpecies,                         &
     conservedVariables, metrics, ratioOfSpecificHeats,                                      &
     jacobianOfInviscidFlux, deltaConservedVariables,                                        &
     specificVolume, velocity, temperature, massFraction, deltaJacobianOfInviscidFlux)

  ! <<< Arguments >>>
  integer, intent(in) :: nDimensions, nSpecies
  SCALAR_TYPE, intent(in) :: conservedVariables(:), metrics(:)
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
  SCALAR_TYPE, intent(out) :: jacobianOfInviscidFlux(:,:)
  SCALAR_TYPE, intent(in), optional :: deltaConservedVariables(:,:), specificVolume,         &
       velocity(:), temperature, massFraction(:)
  SCALAR_TYPE, intent(out), optional :: deltaJacobianOfInviscidFlux(:,:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: k
  SCALAR_TYPE :: specificVolume_, velocity_(nDimensions), temperature_,                      &
       massFraction_(nSpecies), contravariantVelocity, phiSquared,                           &
       deltaConservedVariables_(nDimensions + nSpecies + 2,nDimensions + nSpecies + 2),      &
       deltaSpecificVolume(nDimensions + nSpecies + 2),                                      &
       deltaVelocity(nDimensions,nDimensions + nSpecies + 2),                                &
       deltaTemperature(nDimensions + nSpecies + 2),                                         &
       deltaMassFraction(nSpecies, nDimensions + nSpecies + 2),                               &
       deltaContravariantVelocity(nDimensions + nSpecies + 2),                               &
       deltaPhiSquared(nDimensions + nSpecies + 2)

  assert_key(nDimensions, (1, 2, 3))
  assert(nSpecies >= 0)
  assert(size(conservedVariables) == nDimensions + 2 + nSpecies)
  assert(size(metrics) == size(velocity))

  ! Compute specific volume if it was not specified.
  if (present(specificVolume)) then
     specificVolume_ = specificVolume
  else
     specificVolume_ = 1.0_wp / conservedVariables(1)
  end if

  ! Compute velocity if it was not specified.
  if (present(velocity)) then
     assert(size(velocity) == nDimensions)
     velocity_ = velocity
  else
     do k = 1, nDimensions
        velocity_(k) = specificVolume_ * conservedVariables(k+1)
     end do
  end if

  ! Compute temperature if it was not specified.
  if (present(temperature)) then
     temperature_ = temperature
  else
     temperature_ = ratioOfSpecificHeats * (specificVolume_ * conservedVariables(nDimensions + 2) &
          - 0.5_wp * (sum(velocity_ ** 2)))
  end if

  ! Compute mass fraction if it was not specified.
  if (present(massFraction) .and. nSpecies > 0) then
     assert(size(massFraction) == nSpecies)
     massFraction_ = massFraction
  else
     do k = 1, nSpecies
        massFraction_(k) = conservedVariables(nDimensions+2+k) * specificVolume_
     end do
  end if

  ! Other dependent variables.
  contravariantVelocity = sum(metrics * velocity_) !... not normalized.

  phiSquared = 0.5_wp * (ratioOfSpecificHeats - 1.0_wp) * sum(velocity_**2)

  ! Zero-out Jacobian of inviscid flux.
  jacobianOfInviscidFlux = 0.0_wp

  ! Compute Jacobian of inviscid flux.
  select case (nDimensions)

  case (1)

     jacobianOfInviscidFlux(1,1) = 0.0_wp
     jacobianOfInviscidFlux(2,1) = phiSquared * metrics(1) -                                 &
          contravariantVelocity * velocity_(1)
     jacobianOfInviscidFlux(3,1) = contravariantVelocity * ((ratioOfSpecificHeats - 2.0_wp)  &
          / (ratioOfSpecificHeats - 1.0_wp) * phiSquared - temperature_)
     do k = 1, nSpecies
        jacobianOfInviscidFlux(3+k,1) = -massFraction_(k) * contravariantVelocity
     end do

     jacobianOfInviscidFlux(1,2) = metrics(1)
     jacobianOfInviscidFlux(2,2) = contravariantVelocity -                                   &
          (ratioOfSpecificHeats - 2.0_wp) * velocity_(1) * metrics(1)
     jacobianOfInviscidFlux(3,2) = (temperature_ + phiSquared /                              &
          (ratioOfSpecificHeats - 1.0_wp)) * metrics(1) - (ratioOfSpecificHeats - 1.0_wp) *  &
          contravariantVelocity * velocity_(1)
     do k = 1, nSpecies
        jacobianOfInviscidFlux(3+k,2) = massFraction_(k) * metrics(1)
     end do

     jacobianOfInviscidFlux(1,3) = 0.0_wp
     jacobianOfInviscidFlux(2,3) = (ratioOfSpecificHeats - 1.0_wp) * metrics(1)
     jacobianOfInviscidFlux(3,3) = ratioOfSpecificHeats * contravariantVelocity

     do k = 1, nSpecies
        jacobianOfInviscidFlux(3+k,3+k) = contravariantVelocity
     end do

  case (2)

     jacobianOfInviscidFlux(1,1) = 0.0_wp
     jacobianOfInviscidFlux(2,1) = phiSquared * metrics(1) -                                 &
          contravariantVelocity * velocity_(1)
     jacobianOfInviscidFlux(3,1) = phiSquared * metrics(2) -                                 &
          contravariantVelocity * velocity_(2)
     jacobianOfInviscidFlux(4,1) = contravariantVelocity * ((ratioOfSpecificHeats - 2.0_wp)  &
          / (ratioOfSpecificHeats - 1.0_wp) * phiSquared - temperature_)
     do k = 1, nSpecies
        jacobianOfInviscidFlux(4+k,1) = -massFraction_(k) * contravariantVelocity
     end do

     jacobianOfInviscidFlux(1,2) = metrics(1)
     jacobianOfInviscidFlux(2,2) = contravariantVelocity -                                   &
          (ratioOfSpecificHeats - 2.0_wp) * velocity_(1) * metrics(1)
     jacobianOfInviscidFlux(3,2) = velocity_(2) * metrics(1) -                               &
          (ratioOfSpecificHeats - 1.0_wp) * velocity_(1) * metrics(2)
     jacobianOfInviscidFlux(4,2) = (temperature_ + phiSquared /                              &
          (ratioOfSpecificHeats - 1.0_wp)) * metrics(1) - (ratioOfSpecificHeats - 1.0_wp) *  &
          contravariantVelocity * velocity_(1)
     do k = 1, nSpecies
        jacobianOfInviscidFlux(4+k,2) = massFraction_(k) * metrics(1)
     end do

     jacobianOfInviscidFlux(1,3) = metrics(2)
     jacobianOfInviscidFlux(2,3) = velocity_(1) * metrics(2) -                               &
          (ratioOfSpecificHeats - 1.0_wp) * velocity_(2) * metrics(1)
     jacobianOfInviscidFlux(3,3) = contravariantVelocity - (ratioOfSpecificHeats - 2.0_wp) * &
          velocity_(2) * metrics(2)
     jacobianOfInviscidFlux(4,3) = (temperature_ + phiSquared /                              &
          (ratioOfSpecificHeats - 1.0_wp)) * metrics(2) - (ratioOfSpecificHeats - 1.0_wp) *  &
          contravariantVelocity * velocity_(2)
     do k = 1, nSpecies
        jacobianOfInviscidFlux(4+k,3) = massFraction_(k) * metrics(2)
     end do

     jacobianOfInviscidFlux(1,4) = 0.0_wp
     jacobianOfInviscidFlux(2,4) = (ratioOfSpecificHeats - 1.0_wp) * metrics(1)
     jacobianOfInviscidFlux(3,4) = (ratioOfSpecificHeats - 1.0_wp) * metrics(2)
     jacobianOfInviscidFlux(4,4) = ratioOfSpecificHeats * contravariantVelocity

     do k = 1, nSpecies
        jacobianOfInviscidFlux(4+k,4+k) = contravariantVelocity
     end do

  case (3)

     jacobianOfInviscidFlux(1,1) = 0.0_wp
     jacobianOfInviscidFlux(2,1) = phiSquared * metrics(1) -                                 &
          contravariantVelocity * velocity_(1)
     jacobianOfInviscidFlux(3,1) = phiSquared * metrics(2) -                                 &
          contravariantVelocity * velocity_(2)
     jacobianOfInviscidFlux(4,1) = phiSquared * metrics(3) -                                 &
          contravariantVelocity * velocity_(3)
     jacobianOfInviscidFlux(5,1) = contravariantVelocity * ((ratioOfSpecificHeats - 2.0_wp)  &
          / (ratioOfSpecificHeats - 1.0_wp) * phiSquared - temperature_)
     do k = 1, nSpecies
        jacobianOfInviscidFlux(5+k,1) = -massFraction_(k) * contravariantVelocity
     end do

     jacobianOfInviscidFlux(1,2) = metrics(1)
     jacobianOfInviscidFlux(2,2) = contravariantVelocity -                                   &
          (ratioOfSpecificHeats - 2.0_wp) * velocity_(1) * metrics(1)
     jacobianOfInviscidFlux(3,2) = velocity_(2) * metrics(1) -                               &
          (ratioOfSpecificHeats - 1.0_wp) * velocity_(1) * metrics(2)
     jacobianOfInviscidFlux(4,2) = velocity_(3) * metrics(1) -                               &
          (ratioOfSpecificHeats - 1.0_wp) * velocity_(1) * metrics(3)
     jacobianOfInviscidFlux(5,2) = (temperature_ + phiSquared /                              &
          (ratioOfSpecificHeats - 1.0_wp)) * metrics(1) - (ratioOfSpecificHeats - 1.0_wp) *  &
          contravariantVelocity * velocity_(1)
     do k = 1, nSpecies
        jacobianOfInviscidFlux(5+k,2) = massFraction_(k) * metrics(1)
     end do

     jacobianOfInviscidFlux(1,3) = metrics(2)
     jacobianOfInviscidFlux(2,3) = velocity_(1) * metrics(2) -                               &
          (ratioOfSpecificHeats - 1.0_wp) * velocity_(2) * metrics(1)
     jacobianOfInviscidFlux(3,3) = contravariantVelocity - (ratioOfSpecificHeats - 2.0_wp) * &
          velocity_(2) * metrics(2)
     jacobianOfInviscidFlux(4,3) = velocity_(3) * metrics(2) -                               &
          (ratioOfSpecificHeats - 1.0_wp) * velocity_(2) * metrics(3)
     jacobianOfInviscidFlux(5,3) = (temperature_ + phiSquared /                              &
          (ratioOfSpecificHeats - 1.0_wp)) * metrics(2) - (ratioOfSpecificHeats - 1.0_wp) *  &
          contravariantVelocity * velocity_(2)
     do k = 1, nSpecies
        jacobianOfInviscidFlux(5+k,3) = massFraction_(k) * metrics(2)
     end do

     jacobianOfInviscidFlux(1,4) = metrics(3)
     jacobianOfInviscidFlux(2,4) = velocity_(1) * metrics(3) -                               &
          (ratioOfSpecificHeats - 1.0_wp) * velocity_(3) * metrics(1)
     jacobianOfInviscidFlux(3,4) = velocity_(2) * metrics(3) -                               &
          (ratioOfSpecificHeats - 1.0_wp) * velocity_(3) * metrics(2)
     jacobianOfInviscidFlux(4,4) = contravariantVelocity - (ratioOfSpecificHeats - 2.0_wp) * &
          velocity_(3) * metrics(3)
     jacobianOfInviscidFlux(5,4) = (temperature_ + phiSquared /                              &
          (ratioOfSpecificHeats - 1.0_wp)) * metrics(3) - (ratioOfSpecificHeats - 1.0_wp) *  &
          contravariantVelocity * velocity_(3)
     do k = 1, nSpecies
        jacobianOfInviscidFlux(5+k,4) = massFraction_(k) * metrics(3)
     end do

     jacobianOfInviscidFlux(1,5) = 0.0_wp
     jacobianOfInviscidFlux(2,5) = (ratioOfSpecificHeats - 1.0_wp) * metrics(1)
     jacobianOfInviscidFlux(3,5) = (ratioOfSpecificHeats - 1.0_wp) * metrics(2)
     jacobianOfInviscidFlux(4,5) = (ratioOfSpecificHeats - 1.0_wp) * metrics(3)
     jacobianOfInviscidFlux(5,5) = ratioOfSpecificHeats * contravariantVelocity

     do k = 1, nSpecies
        jacobianOfInviscidFlux(5+k,5+k) = contravariantVelocity
     end do
        
  end select !... nDimensions

  ! Compute variations in Jacobian of inviscid flux.
  if (present(deltaJacobianOfInviscidFlux)) then

     ! Zero-out the variations of the Jacobian of inviscid flux.
     deltaJacobianOfInviscidFlux = 0.0_wp

     select case (nDimensions)

     case (1)

        ! If not specified, use identity matrix for the variation of conservedVariables.
        if (present(deltaConservedVariables)) then
           deltaConservedVariables_ = deltaConservedVariables
        else
           deltaConservedVariables_ = 0.0_wp
           deltaConservedVariables_(1,1) = 1.0_wp
           deltaConservedVariables_(2,2) = 1.0_wp
           deltaConservedVariables_(3,3) = 1.0_wp
           do k = 1, nSpecies
              deltaConservedVariables_(3+k,3+k) = 1.0_wp
           end do
        end if

        ! Compute variations of specific volume, velocity, temperature and mass fraction.
        deltaSpecificVolume = -1.0_wp / conservedVariables(1) ** 2 *                         &
             deltaConservedVariables_(1,:)
        deltaVelocity(1,:) = deltaSpecificVolume * conservedVariables(2) +                   &
             specificVolume_ * deltaConservedVariables_(2,:)
        deltaTemperature = ratioOfSpecificHeats * (deltaSpecificVolume *                     &
             conservedVariables(3) + specificVolume_ * deltaConservedVariables_(2,:) -       &
             (velocity_(1) * deltaVelocity(1,:)))
        do k = 1, nSpecies
           deltaMassFraction(k,:) = deltaSpecificVolume * conservedVariables(3+k) +          &
                specificVolume_ * deltaConservedVariables_(3+k,:)
        end do

        ! Compute variations of other dependent variables:

        deltaContravariantVelocity = metrics(1) * deltaVelocity(1,:)

        deltaPhiSquared = (ratioOfSpecificHeats - 1.0_wp) *                                  &
             (velocity_(1) * deltaVelocity(1,:))

        deltaJacobianOfInviscidFlux(1,1,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(2,1,:) = deltaPhiSquared * metrics(1) -                  &
             deltaContravariantVelocity * velocity_(1) - contravariantVelocity *             &
             deltaVelocity(1,:)
        deltaJacobianOfInviscidFlux(3,1,:) = deltaContravariantVelocity *                    &
             ((ratioOfSpecificHeats - 2.0_wp) / (ratioOfSpecificHeats - 1.0_wp) *            &
             phiSquared - temperature_) + contravariantVelocity *                            &
             ((ratioOfSpecificHeats - 2.0_wp) / (ratioOfSpecificHeats - 1.0_wp) *            &
             deltaPhiSquared - deltaTemperature)
        do k = 1, nSpecies
           deltaJacobianOfInviscidFlux(3+k,1,:) = -massFraction_(k) *                        &
                deltacontravariantVelocity - deltaMassFraction(k,:) *                        &
                contravariantVelocity
        end do

        deltaJacobianOfInviscidFlux(1,2,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(2,2,:) = deltaContravariantVelocity -                    &
             (ratioOfSpecificHeats - 2.0_wp) * deltaVelocity(1,:) * metrics(1)
        deltaJacobianOfInviscidFlux(3,2,:) = (deltaTemperature + deltaPhiSquared /           &
             (ratioOfSpecificHeats - 1.0_wp)) * metrics(1) -                                 &
             (ratioOfSpecificHeats - 1.0_wp) * (deltaContravariantVelocity * velocity_(1) +  &
             contravariantVelocity * deltaVelocity(1,:))
        do k = 1, nSpecies
           deltaJacobianOfInviscidFlux(3+k,2,:) = deltaMassFraction(k,:) * metrics(1)
        end do

        deltaJacobianOfInviscidFlux(1,3,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(2,3,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(3,3,:) = ratioOfSpecificHeats *                          &
             deltaContravariantVelocity
        do k = 1, nSpecies
           deltaJacobianOfInviscidFlux(3+k,3,:) = 0.0_wp
        end do

        do k = 1, nSpecies
           deltaJacobianOfInviscidFlux(3+k,3+k,:) = deltaContravariantVelocity
        end do

     case (2)

        ! If not specified, use identity matrix for the variation of conservedVariables.
        if (present(deltaConservedVariables)) then
           deltaConservedVariables_ = deltaConservedVariables
        else
           deltaConservedVariables_ = 0.0_wp
           deltaConservedVariables_(1,1) = 1.0_wp
           deltaConservedVariables_(2,2) = 1.0_wp
           deltaConservedVariables_(3,3) = 1.0_wp
           deltaConservedVariables_(4,4) = 1.0_wp
           do k = 1, nSpecies
              deltaConservedVariables_(4+k,4+k) = 1.0_wp
           end do
        end if

        ! Compute variations of specific volume, velocity, temperature and mass fraction.
        deltaSpecificVolume = -1.0_wp / conservedVariables(1) ** 2 *                         &
             deltaConservedVariables_(1,:)
        deltaVelocity(1,:) = deltaSpecificVolume * conservedVariables(2) +                   &
             specificVolume_ * deltaConservedVariables(2,:)
        deltaVelocity(2,:) = deltaSpecificVolume * conservedVariables(3) +                   &
             specificVolume_ * deltaConservedVariables(3,:)
        deltaTemperature = ratioOfSpecificHeats * (deltaSpecificVolume *                     &
             conservedVariables(4) + specificVolume_ * deltaConservedVariables_(4,:) -       &
             (velocity_(1) * deltaVelocity(1,:) + velocity_(2) * deltaVelocity(2,:)))
        do k = 1, nSpecies
           deltaMassFraction(k,:) = deltaSpecificVolume * conservedVariables(4+k) +          &
                specificVolume_ * deltaConservedVariables_(4+k,:)
        end do

        ! Compute variations of other dependent variables:

        deltaContravariantVelocity = metrics(1) * deltaVelocity(1,:) + metrics(2) *          &
             deltaVelocity(2,:)

        deltaPhiSquared = (ratioOfSpecificHeats - 1.0_wp) *                                  &
             (velocity_(1) * deltaVelocity(1,:) + velocity_(2) * deltaVelocity(2,:))

        deltaJacobianOfInviscidFlux(1,1,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(2,1,:) = deltaPhiSquared * metrics(1) -                  &
             deltaContravariantVelocity * velocity_(1) -                                     &
             contravariantVelocity * deltaVelocity(1,:)
        deltaJacobianOfInviscidFlux(3,1,:) = deltaPhiSquared * metrics(2) -                  &
             deltaContravariantVelocity * velocity_(2) -                                     &
             contravariantVelocity * deltaVelocity(2,:)
        deltaJacobianOfInviscidFlux(4,1,:) = deltaContravariantVelocity *                    &
             ((ratioOfSpecificHeats - 2.0_wp) / (ratioOfSpecificHeats - 1.0_wp) *            &
             phiSquared - temperature_) + contravariantVelocity *                            &
             ((ratioOfSpecificHeats - 2.0_wp) / (ratioOfSpecificHeats - 1.0_wp) *            &
             deltaPhiSquared - deltaTemperature)
        do k = 1, nSpecies
           deltaJacobianOfInviscidFlux(4+k,1,:) = -massFraction_(k) *                        &
                deltacontravariantVelocity - deltaMassFraction(k,:) *                        &
                contravariantVelocity
        end do

        deltaJacobianOfInviscidFlux(1,2,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(2,2,:) = deltaContravariantVelocity -                    &
             (ratioOfSpecificHeats - 2.0_wp) * deltaVelocity(1,:) * metrics(1)
        deltaJacobianOfInviscidFlux(3,2,:) = deltaVelocity(2,:) * metrics(1) -               &
             (ratioOfSpecificHeats - 1.0_wp) * deltaVelocity(1,:) * metrics(2)
        deltaJacobianOfInviscidFlux(4,2,:) = (deltaTemperature + deltaPhiSquared /           &
             (ratioOfSpecificHeats - 1.0_wp)) * metrics(1) -                                 &
             (ratioOfSpecificHeats - 1.0_wp) * (deltaContravariantVelocity * velocity_(1) +  &
             contravariantVelocity * deltaVelocity(1,:))
        do k = 1, nSpecies
           deltaJacobianOfInviscidFlux(4+k,2,:) = deltaMassFraction(k,:) * metrics(1)
        end do

        deltaJacobianOfInviscidFlux(1,3,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(2,3,:) = deltaVelocity(1,:) * metrics(2) -               &
             (ratioOfSpecificHeats - 1.0_wp) * deltaVelocity(2,:) * metrics(1)
        deltaJacobianOfInviscidFlux(3,3,:) = deltaContravariantVelocity -                    &
             (ratioOfSpecificHeats - 2.0_wp) * deltaVelocity(2,:) * metrics(2)
        deltaJacobianOfInviscidFlux(4,3,:) = (deltaTemperature + deltaPhiSquared /           &
             (ratioOfSpecificHeats - 1.0_wp)) * metrics(2) -                                 &
             (ratioOfSpecificHeats - 1.0_wp) * (deltaContravariantVelocity * velocity_(2) +  &
             contravariantVelocity * deltaVelocity(2,:))
        do k = 1, nSpecies
           deltaJacobianOfInviscidFlux(4+k,3,:) = deltaMassFraction(k,:) * metrics(2)
        end do

        deltaJacobianOfInviscidFlux(1,4,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(2,4,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(3,4,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(4,4,:) = ratioOfSpecificHeats *                          &
             deltaContravariantVelocity
        do k = 1, nSpecies
           deltaJacobianOfInviscidFlux(4+k,4,:) = 0.0_wp
        end do

        do k = 1, nSpecies
           deltaJacobianOfInviscidFlux(4+k,4+k,:) = deltaContravariantVelocity
        end do

     case (3)

        ! If not specified, use identity matrix for the variation of conservedVariables.
        if (present(deltaConservedVariables)) then
           deltaConservedVariables_ = deltaConservedVariables
        else
           deltaConservedVariables_ = 0.0_wp
           deltaConservedVariables_(1,1) = 1.0_wp
           deltaConservedVariables_(2,2) = 1.0_wp
           deltaConservedVariables_(3,3) = 1.0_wp
           deltaConservedVariables_(4,4) = 1.0_wp
           deltaConservedVariables_(5,5) = 1.0_wp
           do k = 1, nSpecies
              deltaConservedVariables_(5+k,5+k) = 1.0_wp
           end do
        end if

        ! Compute variations of specific volume, velocity, temperature and mass fraction.
        deltaSpecificVolume = -1.0_wp / conservedVariables(1) ** 2 *                         &
             deltaConservedVariables_(1,:)
        deltaVelocity(1,:) = deltaSpecificVolume * conservedVariables(2) +                   &
             specificVolume_ * deltaConservedVariables(2,:)
        deltaVelocity(2,:) = deltaSpecificVolume * conservedVariables(3) +                   &
             specificVolume_ * deltaConservedVariables(3,:)
        deltaVelocity(3,:) = deltaSpecificVolume * conservedVariables(4) +                   &
             specificVolume_ * deltaConservedVariables(4,:)
        deltaTemperature = ratioOfSpecificHeats * (deltaSpecificVolume *                     &
             conservedVariables(5) + specificVolume_ * deltaConservedVariables_(5,:) -       &
             (velocity_(1) * deltaVelocity(1,:) + velocity_(2) * deltaVelocity(2,:) +        &
             velocity_(3) * deltaVelocity(3,:)))
        do k = 1, nSpecies
           deltaMassFraction(k,:) = deltaSpecificVolume * conservedVariables(5+k) +          &
                specificVolume_ * deltaConservedVariables_(5+k,:)
        end do

        ! Compute variations of other dependent variables:

        deltaContravariantVelocity =                                                         &
             metrics(1) * deltaVelocity(1,:) +                                               &
             metrics(2) * deltaVelocity(2,:) +                                               &
             metrics(3) * deltaVelocity(3,:)

        deltaPhiSquared = (ratioOfSpecificHeats - 1.0_wp) *                                  &
             (velocity_(1) * deltaVelocity(1,:) + velocity_(2) * deltaVelocity(2,:) +        &
             velocity_(3) * deltaVelocity(3,:))

        deltaJacobianOfInviscidFlux(1,1,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(2,1,:) = deltaPhiSquared * metrics(1) -                  &
             deltaContravariantVelocity * velocity_(1) - contravariantVelocity *             &
             deltaVelocity(1,:)
        deltaJacobianOfInviscidFlux(3,1,:) = deltaPhiSquared * metrics(2) -                  &
             deltaContravariantVelocity * velocity_(2) - contravariantVelocity *             &
             deltaVelocity(2,:)
        deltaJacobianOfInviscidFlux(4,1,:) = deltaPhiSquared * metrics(3) -                  &
             deltaContravariantVelocity * velocity_(3) - contravariantVelocity *             &
             deltaVelocity(3,:)
        deltaJacobianOfInviscidFlux(5,1,:) = deltaContravariantVelocity *                    &
             ((ratioOfSpecificHeats - 2.0_wp) / (ratioOfSpecificHeats - 1.0_wp) *            &
             phiSquared -  temperature_) + contravariantVelocity *                           &
             ((ratioOfSpecificHeats - 2.0_wp) / (ratioOfSpecificHeats - 1.0_wp) *            &
             deltaPhiSquared - deltaTemperature)
        do k = 1, nSpecies
           deltaJacobianOfInviscidFlux(5+k,1,:) = -massFraction_(k) *                        &
                deltacontravariantVelocity - deltaMassFraction(k,:) *                        &
                contravariantVelocity
        end do

        deltaJacobianOfInviscidFlux(1,2,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(2,2,:) = deltaContravariantVelocity -                    &
             (ratioOfSpecificHeats - 2.0_wp) * deltaVelocity(1,:) * metrics(1)
        deltaJacobianOfInviscidFlux(3,2,:) = deltaVelocity(2,:) * metrics(1) -               &
             (ratioOfSpecificHeats - 1.0_wp) * deltaVelocity(1,:) * metrics(2)
        deltaJacobianOfInviscidFlux(4,2,:) = deltaVelocity(3,:) * metrics(1) -               &
             (ratioOfSpecificHeats - 1.0_wp) * deltaVelocity(1,:) * metrics(3)
        deltaJacobianOfInviscidFlux(5,2,:) = (deltaTemperature + deltaPhiSquared /           &
             (ratioOfSpecificHeats - 1.0_wp)) * metrics(1) -                                 &
             (ratioOfSpecificHeats - 1.0_wp) * (deltaContravariantVelocity * velocity_(1) +  &
             contravariantVelocity * deltaVelocity(1,:))
        do k = 1, nSpecies
           deltaJacobianOfInviscidFlux(5+k,2,:) = deltaMassFraction(k,:) * metrics(1)
        end do

        deltaJacobianOfInviscidFlux(1,3,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(2,3,:) = deltaVelocity(1,:) * metrics(2) -               &
             (ratioOfSpecificHeats - 1.0_wp) * deltaVelocity(2,:) * metrics(1)
        deltaJacobianOfInviscidFlux(3,3,:) = deltaContravariantVelocity -                    &
             (ratioOfSpecificHeats - 2.0_wp) * deltaVelocity(2,:) * metrics(2)
        deltaJacobianOfInviscidFlux(4,3,:) = deltaVelocity(3,:) * metrics(2) -               &
             (ratioOfSpecificHeats - 1.0_wp) * deltaVelocity(2,:) * metrics(3)
        deltaJacobianOfInviscidFlux(5,3,:) = (deltaTemperature + deltaPhiSquared /           &
             (ratioOfSpecificHeats - 1.0_wp)) * metrics(2) -                                 &
             (ratioOfSpecificHeats - 1.0_wp) * (deltaContravariantVelocity * velocity_(2) +  &
             contravariantVelocity * deltaVelocity(2,:))
        do k = 1, nSpecies
           deltaJacobianOfInviscidFlux(5+k,3,:) = deltaMassFraction(k,:) * metrics(2)
        end do

        deltaJacobianOfInviscidFlux(1,4,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(2,4,:) = deltaVelocity(1,:) * metrics(3) -               &
             (ratioOfSpecificHeats - 1.0_wp) * deltaVelocity(3,:) * metrics(1)
        deltaJacobianOfInviscidFlux(3,4,:) = deltaVelocity(2,:) * metrics(3) -               &
             (ratioOfSpecificHeats - 1.0_wp) * deltaVelocity(3,:) * metrics(2)
        deltaJacobianOfInviscidFlux(4,4,:) = deltaContravariantVelocity -                    &
             (ratioOfSpecificHeats - 2.0_wp) * deltaVelocity(3,:) * metrics(3)
        deltaJacobianOfInviscidFlux(5,4,:) = (deltaTemperature + deltaPhiSquared /           &
             (ratioOfSpecificHeats - 1.0_wp)) * metrics(3) -                                 &
             (ratioOfSpecificHeats - 1.0_wp) * (deltaContravariantVelocity * velocity_(3) +  &
             contravariantVelocity * deltaVelocity(3,:))
        do k = 1, nSpecies
           deltaJacobianOfInviscidFlux(5+k,4,:) = deltaMassFraction(k,:) * metrics(3)
        end do

        deltaJacobianOfInviscidFlux(1,5,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(2,5,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(3,5,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(4,5,:) = 0.0_wp
        deltaJacobianOfInviscidFlux(5,5,:) = ratioOfSpecificHeats *                          &
             deltaContravariantVelocity
        do k = 1, nSpecies
           deltaJacobianOfInviscidFlux(5+k,5,:) = 0.0_wp
        end do

        do k = 1, nSpecies
           deltaJacobianOfInviscidFlux(5+k,5+k,:) = deltaContravariantVelocity
        end do

     end select !... nDimensions

  end if

end subroutine computeJacobianOfInviscidFlux

PURE_SUBROUTINE computeIncomingJacobianOfInviscidFlux(nDimensions, nSpecies,                 &
     conservedVariables, metrics, ratioOfSpecificHeats, incomingDirection,                   &
     incomingJacobianOfInviscidFlux, deltaIncomingJacobianOfInviscidFlux,                    &
     deltaConservedVariables, specificVolume, velocity, temperature, massFraction)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: nDimensions, nSpecies, incomingDirection
  SCALAR_TYPE, intent(in) :: conservedVariables(:), metrics(:)
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
  SCALAR_TYPE, intent(out) :: incomingJacobianOfInviscidFlux(:,:)
  SCALAR_TYPE, intent(out), optional :: deltaIncomingJacobianOfInviscidFlux(:,:,:)
  SCALAR_TYPE, intent(in), optional :: deltaConservedVariables(:,:), specificVolume,         &
       velocity(:), temperature, massFraction(:)

  ! <<< Local variables >>>
  integer :: i, j, k
  integer, parameter :: wp = SCALAR_KIND
  SCALAR_TYPE :: arcLength, normalizedMetrics(nDimensions), specificVolume_,                 &
       velocity_(nDimensions), temperature_, massFraction_(nSpecies),                        &
       contravariantVelocity, speedOfSound, phiSquared,                                      &
       rightEigenvectors(nDimensions + nSpecies + 2,nDimensions + nSpecies + 2),             &
       eigenvalues(nDimensions + nSpecies + 2),                                              &
       leftEigenvectors(nDimensions + nSpecies + 2,nDimensions + nSpecies + 2),              &
       deltaConservedVariables_(nDimensions + nSpecies + 2,nDimensions + nSpecies + 2),      &
       deltaSpecificVolume(nDimensions + nSpecies + 2),                                      &
       deltaVelocity(nDimensions,nDimensions + nSpecies + 2),                                &
       deltaTemperature(nDimensions + nSpecies + 2),                                         &
       deltaMassFraction(nSpecies, nDimensions + nSpecies + 2),                              &
       deltaContravariantVelocity(nDimensions + nSpecies + 2),                               &
       deltaSpeedOfSound(nDimensions + nSpecies + 2),                                        &
       deltaPhiSquared(nDimensions + nSpecies + 2), deltaRightEigenvectors                   &
       (nDimensions + nSpecies + 2,nDimensions + nSpecies + 2,nDimensions + nSpecies + 2),   &
       deltaEigenvalues(nDimensions + nSpecies + 2,nDimensions + nSpecies + 2),              &
       deltaLeftEigenvectors                                                                 &
       (nDimensions + nSpecies + 2,nDimensions + nSpecies + 2,nDimensions + nSpecies + 2),   &
       temp(nDimensions + nSpecies + 2)

  ! Normalize the metrics.
  arcLength = sqrt(sum(metrics ** 2))
  normalizedMetrics = metrics / arcLength

  ! Compute specific volume if it was not specified.
  if (present(specificVolume)) then
     specificVolume_ = specificVolume
  else
     specificVolume_ = 1.0_wp / conservedVariables(1)
  end if

  ! Compute velocity if it was not specified.
  if (present(velocity)) then
     velocity_ = velocity
  else
     do i = 1, nDimensions
        velocity_(i) = specificVolume_ * conservedVariables(i + 1)
     end do
  end if

  ! Compute temperature if it was not specified.
  if (present(temperature)) then
     temperature_ = temperature
  else
     temperature_ = ratioOfSpecificHeats * (specificVolume_ *                                &
          conservedVariables(nDimensions + 2) - 0.5_wp * sum(velocity_ ** 2))
  end if

  ! Compute mass fraction if it was not specified.
  if (present(massFraction) .and. nSpecies > 0) then
     massFraction_ = massFraction
  else
     do k = 1, nSpecies
        massFraction_(k) = specificVolume_ * conservedVariables(nDimensions+2+k)
     end do
  end if

  ! Other dependent variables.
  contravariantVelocity = sum(normalizedMetrics * velocity_)
  speedOfSound = sqrt((ratioOfSpecificHeats - 1.0_wp) * temperature_)
  phiSquared = 0.5_wp * (ratioOfSpecificHeats - 1.0_wp) * sum(velocity_ ** 2)

  select case (nDimensions)

  case (1)

     ! Eigenvalues.
     eigenvalues(1) = contravariantVelocity
     eigenvalues(2) = contravariantVelocity + speedOfSound
     eigenvalues(3) = contravariantVelocity - speedOfSound
     do k = 1, nSpecies
        eigenvalues(3+k) = contravariantVelocity
     end do
     eigenvalues = arcLength * eigenvalues

     if (present(deltaIncomingJacobianOfInviscidFlux)) then

        ! If not specified, use identity matrix for the variation of conservedVariables.
        if (present(deltaConservedVariables)) then
           deltaConservedVariables_ = deltaConservedVariables
        else
           deltaConservedVariables_ = 0.0_wp
           deltaConservedVariables_(1,1) = 1.0_wp
           deltaConservedVariables_(2,2) = 1.0_wp
           deltaConservedVariables_(3,3) = 1.0_wp
           do k = 1, nSpecies
              deltaConservedVariables_(3+k,3+k) = 1.0_wp
           end do
        end if

        ! Compute variations of specific volume, velocity, temperature and mass fraction.
        deltaSpecificVolume = - specificVolume_ ** 2 * deltaConservedVariables_(1,:)
        deltaVelocity(1,:) = deltaSpecificVolume * conservedVariables(2) +                   &
             specificVolume_ * deltaConservedVariables_(2,:)
        deltaTemperature = ratioOfSpecificHeats *                                            &
             (deltaSpecificVolume * conservedVariables(3) +                                  &
             specificVolume_ * deltaConservedVariables_(3,:) -                               &
             (velocity_(1) * deltaVelocity(1,:)))
        do k = 1, nSpecies
           deltaMassFraction(k,:) = deltaSpecificVolume * conservedVariables(3+k) +          &
             specificVolume_ * deltaConservedVariables_(3+k,:)
        end do

        ! Compute variations of other dependent variables.
        deltaContravariantVelocity = normalizedMetrics(1) * deltaVelocity(1,:)
        deltaSpeedOfSound = 0.5_wp / speedOfSound *                                          &
             (ratioOfSpecificHeats - 1.0_wp) * deltaTemperature
        deltaPhiSquared = (ratioOfSpecificHeats - 1.0_wp) *                                  &
             (velocity_(1) * deltaVelocity(1,:))

        ! Variation of matrix containing eigenvalues.
        deltaEigenvalues(1,:) = deltaContravariantVelocity
        deltaEigenvalues(2,:) = deltaContravariantVelocity + deltaSpeedOfSound
        deltaEigenvalues(3,:) = deltaContravariantVelocity - deltaSpeedOfSound
        do k = 1, nSpecies
           deltaEigenvalues(3+k,:) = deltaContravariantVelocity
        end do
        deltaEigenvalues = arcLength * deltaEigenvalues

     end if

     ! Zero-out the eigenvalues corresponding to outgoing characteristics and
     ! corresponding variations.
     do i = 1, 3 + nSpecies
        if (incomingDirection * real(eigenvalues(i), wp) < 0.0_wp) then
           eigenvalues(i) = 0.0_wp
           if (present(deltaIncomingJacobianOfInviscidFlux)) deltaEigenvalues(i,:) = 0.0_wp
        end if
     end do

     ! Matrix whose columns are the right eigenvectors:

     rightEigenvectors = 0.0_wp

     rightEigenvectors(1,1) = 1.0_wp
     rightEigenvectors(2,1) = velocity_(1)
     rightEigenvectors(3,1) = phiSquared / (ratioOfSpecificHeats - 1.0_wp)
     do k = 1, nSpecies
        rightEigenvectors(3+k,1) = massFraction_(k)
     end do

     rightEigenvectors(1,2) = 1.0_wp
     rightEigenvectors(2,2) = velocity_(1) + normalizedMetrics(1) * speedOfSound
     rightEigenvectors(3,2) = temperature_ + phiSquared / (ratioOfSpecificHeats - 1.0_wp) +  &
          speedOfSound * contravariantVelocity
     do k = 1, nSpecies
        rightEigenvectors(3+k,2) = massFraction_(k)
     end do

     rightEigenvectors(1,3) = 1.0_wp
     rightEigenvectors(2,3) = velocity_(1) - normalizedMetrics(1) * speedOfSound
     rightEigenvectors(3,3) = temperature_ + phiSquared / (ratioOfSpecificHeats - 1.0_wp) -  &
          speedOfSound * contravariantVelocity
     do k = 1, nSpecies
        rightEigenvectors(3+k,3) = massFraction_(k)
     end do

     do k = 1, nSpecies
        rightEigenvectors(3+k,3+k) = 1.0_wp
     end do

     ! Matrix whose rows are the left eigenvectors:

     leftEigenvectors = 0.0_wp

     leftEigenvectors(1,1) = 1.0_wp - phiSquared / speedOfSound ** 2
     leftEigenvectors(2,1) = 0.5_wp * (phiSquared / speedOfSound ** 2 -                      &
          contravariantVelocity / speedOfSound)
     leftEigenvectors(3,1) = 0.5_wp * (phiSquared / speedOfSound ** 2 +                      &
          contravariantVelocity / speedOfSound)
     do k = 1, nSpecies
        leftEigenvectors(3+k,1) = -massFraction(k)
     end do

     leftEigenvectors(1,2) = velocity_(1) / temperature_
     leftEigenvectors(2,2) = - 0.5_wp * (velocity_(1) / temperature_ -                       &
          normalizedMetrics(1) / speedOfSound)
     leftEigenvectors(3,2) = - 0.5_wp * (velocity_(1) / temperature_ +                       &
          normalizedMetrics(1) / speedOfSound)
     do k = 1, nSpecies
        leftEigenvectors(3+k,2) = 0.0_wp
     end do

     leftEigenvectors(1,3) = - 1.0_wp / temperature_
     leftEigenvectors(2,3) = 0.5_wp / temperature_
     leftEigenvectors(3,3) = 0.5_wp / temperature_
     do k = 1, nSpecies
        leftEigenvectors(3+k,3) = 0.0_wp
     end do

     do k = 1, nSpecies
        leftEigenvectors(3+k,3+k) = 1.0_wp
     end do

     ! ``Incoming'' part.
     do j = 1, 3 + nSpecies
        do i = 1, 3 + nSpecies
           incomingJacobianOfInviscidFlux(i,j) =                                             &
                rightEigenvectors(i,1) * eigenvalues(1) * leftEigenvectors(1,j) +            &
                rightEigenvectors(i,2) * eigenvalues(2) * leftEigenvectors(2,j) +            &
                rightEigenvectors(i,3) * eigenvalues(3) * leftEigenvectors(3,j)
           do k = 1, nSpecies
              incomingJacobianOfInviscidFlux(i,j) = incomingJacobianOfInviscidFlux(i,j) +    &
                   rightEigenvectors(i,3+k) * eigenvalues(3+k) * leftEigenvectors(3+k,j)
           end do
        end do
     end do

     if (present(deltaIncomingJacobianOfInviscidFlux)) then

        ! Variation of the matrix whose columns are the right eigenvectors:

        deltaRightEigenvectors = 0.0_wp

        deltaRightEigenvectors(1,1,:) = 0.0_wp
        deltaRightEigenvectors(2,1,:) = deltaVelocity(1,:)
        deltaRightEigenvectors(3,1,:) = deltaPhiSquared /                                    &
             (ratioOfSpecificHeats - 1.0_wp)
        do k = 1, nSpecies
           deltaRightEigenvectors(3+k,1,:) = deltaMassFraction(k,:)
        end do

        deltaRightEigenvectors(1,2,:) = 0.0_wp
        deltaRightEigenvectors(2,2,:) = deltaVelocity(1,:) +                                 &
             normalizedMetrics(1) * deltaSpeedOfSound
        deltaRightEigenvectors(3,2,:) = deltaTemperature +                                   &
             deltaPhiSquared / (ratioOfSpecificHeats - 1.0_wp) +                             &
             deltaSpeedOfSound * contravariantVelocity +                                     &
             speedOfSound * deltaContravariantVelocity
        do k = 1, nSpecies
           deltaRightEigenvectors(3+k,2,:) = deltaMassFraction(k,:)
        end do

        deltaRightEigenvectors(1,3,:) = 0.0_wp
        deltaRightEigenvectors(2,3,:) = deltaVelocity(1,:) -                                 &
             normalizedMetrics(1) * deltaSpeedOfSound
        deltaRightEigenvectors(3,3,:) = deltaTemperature +                                   &
             deltaPhiSquared / (ratioOfSpecificHeats - 1.0_wp) -                             &
             deltaSpeedOfSound * contravariantVelocity -                                     &
             speedOfSound * deltaContravariantVelocity
        do k = 1, nSpecies
           deltaRightEigenvectors(3+k,3,:) = deltaMassFraction(k,:)
        end do

        do k = 1, nSpecies
           deltaRightEigenvectors(3+k,3+k,:) = 0.0_wp
        end do

        ! Variation of the matrix whose rows are the left eigenvectors:

        deltaLeftEigenvectors = 0.0_wp

        temp = deltaPhiSquared / speedOfSound ** 2 -                                         &
             2.0_wp * phiSquared / speedOfSound ** 3 * deltaSpeedOfSound
        deltaLeftEigenvectors(1,1,:) = -temp
        deltaLeftEigenvectors(2,1,:) = 0.5_wp * (temp -                                      &
             deltaContravariantVelocity / speedOfSound +                                     &
             contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)
        deltaLeftEigenvectors(3,1,:) = 0.5_wp * (temp +                                      &
             deltaContravariantVelocity / speedOfSound -                                     &
             contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)
        do k = 1, nSpecies
           deltaLeftEigenvectors(3+k,1,:) = -deltaMassFraction(k,:)
        end do

        temp = deltaVelocity(1,:) / temperature_ -                                           &
             velocity_(1) / temperature_ ** 2 * deltaTemperature
        deltaLeftEigenvectors(1,2,:) = temp
        deltaLeftEigenvectors(2,2,:) = -0.5_wp * (temp +                                     &
             normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)
        deltaLeftEigenvectors(3,2,:) = -0.5_wp * (temp -                                     &
             normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)
        do k = 1, nSpecies
           deltaLeftEigenvectors(3+k,2,:) = 0.0_wp
        end do

        temp =  -1.0_wp / temperature_ ** 2 * deltaTemperature
        deltaLeftEigenvectors(1,3,:) = -temp
        deltaLeftEigenvectors(2,3,:) = 0.5_wp * temp
        deltaLeftEigenvectors(3,3,:) = 0.5_wp * temp
        do k = 1, nSpecies
           deltaLeftEigenvectors(3+k,3,:) = 0.0_wp
        end do

        do k = 1, nSpecies
           deltaLeftEigenvectors(3+k,3+k,:) = 0.0_wp
        end do

        ! Variation of the ``incoming'' part.
        do j = 1, 3 + nSpecies
           do i = 1, 3 + nSpecies
              deltaIncomingJacobianOfInviscidFlux(i,j,:) =                                   &
                   deltaRightEigenvectors(i,1,:) * eigenvalues(1) * leftEigenvectors(1,j) +  &
                   deltaRightEigenvectors(i,2,:) * eigenvalues(2) * leftEigenvectors(2,j) +  &
                   deltaRightEigenvectors(i,3,:) * eigenvalues(3) * leftEigenvectors(3,j) +  &
                   rightEigenvectors(i,1) * deltaEigenvalues(1,:) * leftEigenvectors(1,j) +  &
                   rightEigenvectors(i,2) * deltaEigenvalues(2,:) * leftEigenvectors(2,j) +  &
                   rightEigenvectors(i,3) * deltaEigenvalues(3,:) * leftEigenvectors(3,j) +  &
                   rightEigenvectors(i,1) * eigenvalues(1) * deltaLeftEigenvectors(1,j,:) +  &
                   rightEigenvectors(i,2) * eigenvalues(2) * deltaLeftEigenvectors(2,j,:) +  &
                   rightEigenvectors(i,3) * eigenvalues(3) * deltaLeftEigenvectors(3,j,:)
              do k = 1, nSpecies
                 deltaIncomingJacobianOfInviscidFlux(i,j,:) =                                &
                      deltaIncomingJacobianOfInviscidFlux(i,j,:) +                           &
                      deltaRightEigenvectors(i,3+k,:) * eigenvalues(3+k) *                   &
                      leftEigenvectors(3+k,j) + rightEigenvectors(i,3+k) *                   &
                      deltaEigenvalues(3+k,:) * leftEigenvectors(3+k,j) +                    &
                      rightEigenvectors(i,3+k) * eigenvalues(3+k) *                          &
                      deltaLeftEigenvectors(3+k,j,:)
              end do
           end do
        end do

     end if

  case (2)

     ! Eigenvalues.
     eigenvalues(1) = contravariantVelocity
     eigenvalues(2) = contravariantVelocity
     eigenvalues(3) = contravariantVelocity + speedOfSound
     eigenvalues(4) = contravariantVelocity - speedOfSound
     do k = 1, nSpecies
        eigenvalues(4+k) = contravariantVelocity
     end do
     eigenvalues = arcLength * eigenvalues

     if (present(deltaIncomingJacobianOfInviscidFlux)) then

        ! If not specified, use identity matrix for the variation of conservedVariables.
        if (present(deltaConservedVariables)) then
           deltaConservedVariables_ = deltaConservedVariables
        else
           deltaConservedVariables_ = 0.0_wp
           deltaConservedVariables_(1,1) = 1.0_wp
           deltaConservedVariables_(2,2) = 1.0_wp
           deltaConservedVariables_(3,3) = 1.0_wp
           deltaConservedVariables_(4,4) = 1.0_wp
           do k = 1, nSpecies
              deltaConservedVariables_(4+k,4+k) = 1.0_wp
           end do
        end if

        ! Compute variations of specific volume, velocity, temperature and mass fraction.
        deltaSpecificVolume = - specificVolume_ ** 2 * deltaConservedVariables_(1,:)
        deltaVelocity(1,:) = deltaSpecificVolume * conservedVariables(2) +                   &
             specificVolume_ * deltaConservedVariables_(2,:)
        deltaVelocity(2,:) = deltaSpecificVolume * conservedVariables(3) +                   &
             specificVolume_ * deltaConservedVariables_(3,:)
        deltaTemperature = ratioOfSpecificHeats *                                            &
             (deltaSpecificVolume * conservedVariables(4) +                                  &
             specificVolume_ * deltaConservedVariables_(4,:) -                               &
             (velocity_(1) * deltaVelocity(1,:) + velocity_(2) * deltaVelocity(2,:)))
        do k = 1, nSpecies
           deltaMassFraction(k,:) = deltaSpecificVolume * conservedVariables(4+k) +          &
                specificVolume_ * deltaConservedVariables_(4+k,:)
        end do

        ! Compute variations of other dependent variables.
        deltaContravariantVelocity = normalizedMetrics(1) * deltaVelocity(1,:) +             &
             normalizedMetrics(2) * deltaVelocity(2,:)
        deltaSpeedOfSound = 0.5_wp / speedOfSound *                                          &
             (ratioOfSpecificHeats - 1.0_wp) * deltaTemperature
        deltaPhiSquared = (ratioOfSpecificHeats - 1.0_wp) *                                  &
             (velocity_(1) * deltaVelocity(1,:) + velocity_(2) * deltaVelocity(2,:))

        ! Variation of matrix containing eigenvalues.
        deltaEigenvalues(1,:) = deltaContravariantVelocity
        deltaEigenvalues(2,:) = deltaContravariantVelocity
        deltaEigenvalues(3,:) = deltaContravariantVelocity + deltaSpeedOfSound
        deltaEigenvalues(4,:) = deltaContravariantVelocity - deltaSpeedOfSound
        do k = 1, nSpecies
           deltaEigenvalues(4+k,:) = deltaContravariantVelocity
        end do
        deltaEigenvalues = arcLength * deltaEigenvalues

     end if

     ! Zero-out the eigenvalues corresponding to outgoing characteristics and corresponding
     ! variations.
     do i = 1, 4 + nSpecies
        if (incomingDirection * real(eigenvalues(i), wp) < 0.0_wp) then
           eigenvalues(i) = 0.0_wp
           if (present(deltaIncomingJacobianOfInviscidFlux)) deltaEigenvalues(i,:) = 0.0_wp
        end if
     end do

     ! Matrix whose columns are the right eigenvectors:

     rightEigenvectors = 0.0_wp

     rightEigenvectors(1,1) = 1.0_wp
     rightEigenvectors(2,1) = velocity_(1)
     rightEigenvectors(3,1) = velocity_(2)
     rightEigenvectors(4,1) = phiSquared / (ratioOfSpecificHeats - 1.0_wp)
     do k = 1, nSpecies
        rightEigenvectors(4+k,1) = massFraction_(k)
     end do

     rightEigenvectors(1,2) = 0.0_wp
     rightEigenvectors(2,2) = normalizedMetrics(2) * conservedVariables(1)
     rightEigenvectors(3,2) = - normalizedMetrics(1) * conservedVariables(1)
     rightEigenvectors(4,2) = conservedVariables(1) * (normalizedMetrics(2) * velocity_(1) - &
          normalizedMetrics(1) * velocity_(2))
     do k = 1, nSpecies
        rightEigenvectors(4+k,2) = 0.0_wp
     end do

     rightEigenvectors(1,3) = 1.0_wp
     rightEigenvectors(2,3) = velocity_(1) + normalizedMetrics(1) * speedOfSound
     rightEigenvectors(3,3) = velocity_(2) + normalizedMetrics(2) * speedOfSound
     rightEigenvectors(4,3) = temperature_ + phiSquared / (ratioOfSpecificHeats - 1.0_wp) +  &
          speedOfSound * contravariantVelocity
     do k = 1, nSpecies
        rightEigenvectors(4+k,3) = massFraction_(k)
     end do

     rightEigenvectors(1,4) = 1.0_wp
     rightEigenvectors(2,4) = velocity_(1) - normalizedMetrics(1) * speedOfSound
     rightEigenvectors(3,4) = velocity_(2) - normalizedMetrics(2) * speedOfSound
     rightEigenvectors(4,4) = temperature_ + phiSquared / (ratioOfSpecificHeats - 1.0_wp) -  &
          speedOfSound * contravariantVelocity
     do k = 1, nSpecies
        rightEigenvectors(4+k,4) = massFraction_(k)
     end do

     do k = 1, nSpecies
        rightEigenvectors(4+k,4+k) = 1.0_wp
     end do

     ! Matrix whose rows are the left eigenvectors:

     leftEigenvectors = 0.0_wp

     leftEigenvectors(1,1) = 1.0_wp - phiSquared / speedOfSound ** 2
     leftEigenvectors(2,1) = - specificVolume_ * (normalizedMetrics(2) * velocity_(1) -      &
          normalizedMetrics(1) * velocity_(2))
     leftEigenvectors(3,1) = 0.5_wp * (phiSquared / speedOfSound ** 2 -                      &
          contravariantVelocity / speedOfSound)
     leftEigenvectors(4,1) = 0.5_wp * (phiSquared / speedOfSound ** 2 +                      &
          contravariantVelocity / speedOfSound)
     do k = 1, nSpecies
        leftEigenvectors(4+k,1) = -massFraction_(k)
     end do

     leftEigenvectors(1,2) = velocity_(1) / temperature_
     leftEigenvectors(2,2) = specificVolume_ * normalizedMetrics(2)
     leftEigenvectors(3,2) = - 0.5_wp * (velocity_(1) / temperature_ -                       &
          normalizedMetrics(1) / speedOfSound)
     leftEigenvectors(4,2) = - 0.5_wp * (velocity_(1) / temperature_ +                       &
          normalizedMetrics(1) / speedOfSound)
     do k = 1, nSpecies
        leftEigenvectors(4+k,2) = 0.0_wp
     end do

     leftEigenvectors(1,3) = velocity_(2) / temperature_
     leftEigenvectors(2,3) = - specificVolume_ * normalizedMetrics(1)
     leftEigenvectors(3,3) = - 0.5_wp * (velocity_(2) / temperature_ -                       &
          normalizedMetrics(2) / speedOfSound)
     leftEigenvectors(4,3) = - 0.5_wp * (velocity_(2) / temperature_ +                       &
          normalizedMetrics(2) / speedOfSound)
     do k = 1, nSpecies
        leftEigenvectors(4+k,3) = 0.0_wp
     end do

     leftEigenvectors(1,4) = - 1.0_wp / temperature_
     leftEigenvectors(2,4) = 0.0_wp
     leftEigenvectors(3,4) = 0.5_wp / temperature_
     leftEigenvectors(4,4) = 0.5_wp / temperature_
     do k = 1, nSpecies
        leftEigenvectors(4+k,4) = 0.0_wp
     end do

     do k = 1, nSpecies
        leftEigenvectors(4+k,4+k) = 1.0_wp
     end do

     ! ``Incoming'' part.
     do j = 1, 4 + nSpecies
        do i = 1, 4 + nSpecies
           incomingJacobianOfInviscidFlux(i,j) =                                             &
                rightEigenvectors(i,1) * eigenvalues(1) * leftEigenvectors(1,j) +            &
                rightEigenvectors(i,2) * eigenvalues(2) * leftEigenvectors(2,j) +            &
                rightEigenvectors(i,3) * eigenvalues(3) * leftEigenvectors(3,j) +            &
                rightEigenvectors(i,4) * eigenvalues(4) * leftEigenvectors(4,j)
           do k = 1, nSpecies
              incomingJacobianOfInviscidFlux(i,j) = incomingJacobianOfInviscidFlux(i,j) +    &
                   rightEigenvectors(i,4+k) * eigenvalues(4+k) * leftEigenvectors(4+k,j)
           end do
        end do
     end do

     if (present(deltaIncomingJacobianOfInviscidFlux)) then

        ! Variation of the matrix whose columns are the right eigenvectors:

        deltaRightEigenvectors = 0.0_wp

        deltaRightEigenvectors(1,1,:) = 0.0_wp
        deltaRightEigenvectors(2,1,:) = deltaVelocity(1,:)
        deltaRightEigenvectors(3,1,:) = deltaVelocity(2,:)
        deltaRightEigenvectors(4,1,:) = deltaPhiSquared /                                    &
             (ratioOfSpecificHeats - 1.0_wp)
        do k = 1, nSpecies
           deltaRightEigenvectors(4+k,1,:) = deltaMassFraction(k,:)
        end do

        deltaRightEigenvectors(1,2,:) = 0.0_wp
        deltaRightEigenvectors(2,2,:) = normalizedMetrics(2) * deltaConservedVariables_(1,:)
        deltaRightEigenvectors(3,2,:) = -normalizedMetrics(1) * deltaConservedVariables_(1,:)
        deltaRightEigenvectors(4,2,:) = deltaConservedVariables_(1,:) *                      &
             (normalizedMetrics(2) * velocity_(1) - normalizedMetrics(1) * velocity_(2)) +   &
             conservedVariables(1) * (normalizedMetrics(2) * deltaVelocity(1,:) -            &
             normalizedMetrics(1) * deltaVelocity(2,:))
        do k = 1, nSpecies
           deltaRightEigenvectors(4+k,2,:) = 0.0_wp
        end do

        deltaRightEigenvectors(1,3,:) = 0.0_wp
        deltaRightEigenvectors(2,3,:) = deltaVelocity(1,:) +                                 &
             normalizedMetrics(1) * deltaSpeedOfSound
        deltaRightEigenvectors(3,3,:) = deltaVelocity(2,:) +                                 &
             normalizedMetrics(2) * deltaSpeedOfSound
        deltaRightEigenvectors(4,3,:) = deltaTemperature +                                   &
             deltaPhiSquared / (ratioOfSpecificHeats - 1.0_wp) +                             &
             deltaSpeedOfSound * contravariantVelocity +                                     &
             speedOfSound * deltaContravariantVelocity
        do k = 1, nSpecies
           deltaRightEigenvectors(4+k,3,:) = deltaMassFraction(k,:)
        end do

        deltaRightEigenvectors(1,4,:) = 0.0_wp
        deltaRightEigenvectors(2,4,:) = deltaVelocity(1,:) -                                 &
             normalizedMetrics(1) * deltaSpeedOfSound
        deltaRightEigenvectors(3,4,:) = deltaVelocity(2,:) -                                 &
             normalizedMetrics(2) * deltaSpeedOfSound
        deltaRightEigenvectors(4,4,:) = deltaTemperature +                                   &
             deltaPhiSquared / (ratioOfSpecificHeats - 1.0_wp) -                             &
             deltaSpeedOfSound * contravariantVelocity -                                     &
             speedOfSound * deltaContravariantVelocity
        do k = 1, nSpecies
           deltaRightEigenvectors(4+k,4,:) = deltaMassFraction(k,:)
        end do

        do k = 1, nSpecies
           deltaRightEigenvectors(4+k,4+k,:) = 0.0_wp
        end do

        ! Variation of the matrix whose rows are the left eigenvectors:

        deltaLeftEigenvectors = 0.0_wp

        temp = deltaPhiSquared / speedOfSound ** 2 -                                         &
             2.0_wp * phiSquared / speedOfSound ** 3 * deltaSpeedOfSound
        deltaLeftEigenvectors(1,1,:) = -temp
        deltaLeftEigenvectors(2,1,:) = -deltaSpecificVolume *                                &
             (normalizedMetrics(2) * velocity_(1) - normalizedMetrics(1) * velocity_(2)) -   &
             specificVolume_ * (normalizedMetrics(2) * deltaVelocity(1,:) -                  &
             normalizedMetrics(1) * deltaVelocity(2,:))
        deltaLeftEigenvectors(3,1,:) = 0.5_wp * (temp -                                      &
             deltaContravariantVelocity / speedOfSound +                                     &
             contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)
        deltaLeftEigenvectors(4,1,:) = 0.5_wp * (temp +                                      &
             deltaContravariantVelocity / speedOfSound -                                     &
             contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)
        do k = 1, nSpecies
           deltaLeftEigenvectors(4+k,1,:) = -deltaMassFraction(k,:)
        end do
        
        temp = deltaVelocity(1,:) / temperature_ -                                           &
             velocity_(1) / temperature_ ** 2 * deltaTemperature
        deltaLeftEigenvectors(1,2,:) = temp
        deltaLeftEigenvectors(2,2,:) = deltaSpecificVolume * normalizedMetrics(2)
        deltaLeftEigenvectors(3,2,:) = -0.5_wp * (temp +                                     &
             normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)
        deltaLeftEigenvectors(4,2,:) = -0.5_wp * (temp -                                     &
             normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)

        temp = deltaVelocity(2,:) / temperature_ -                                           &
             velocity_(2) / temperature_ ** 2 * deltaTemperature
        do k = 1, nSpecies
           deltaLeftEigenvectors(4+k,2,:) = 0.0_wp
        end do

        deltaLeftEigenvectors(1,3,:) = temp
        deltaLeftEigenvectors(2,3,:) = -deltaSpecificVolume * normalizedMetrics(1)
        deltaLeftEigenvectors(3,3,:) = -0.5_wp * (temp +                                     &
             normalizedMetrics(2) / speedOfSound ** 2 * deltaSpeedOfSound)
        deltaLeftEigenvectors(4,3,:) = -0.5_wp * (temp -                                     &
             normalizedMetrics(2) / speedOfSound ** 2 * deltaSpeedOfSound)
        do k = 1, nSpecies
           deltaLeftEigenvectors(4+k,3,:) = 0.0_wp
        end do

        temp =  -1.0_wp / temperature_ ** 2 * deltaTemperature
        deltaLeftEigenvectors(1,4,:) = -temp
        deltaLeftEigenvectors(2,4,:) = 0.0_wp
        deltaLeftEigenvectors(3,4,:) = 0.5_wp * temp
        deltaLeftEigenvectors(4,4,:) = 0.5_wp * temp
        do k = 1, nSpecies
           deltaLeftEigenvectors(4+k,4,:) = 0.0_wp
        end do

        do k = 1, nSpecies
           deltaLeftEigenvectors(4+k,4+k,:) = 0.0_wp
        end do

        ! Variation of the ``incoming'' part.
        do j = 1, 4 + nSpecies
           do i = 1, 4 + nSpecies
              deltaIncomingJacobianOfInviscidFlux(i,j,:) =                                   &
                   deltaRightEigenvectors(i,1,:) * eigenvalues(1) * leftEigenvectors(1,j) +  &
                   deltaRightEigenvectors(i,2,:) * eigenvalues(2) * leftEigenvectors(2,j) +  &
                   deltaRightEigenvectors(i,3,:) * eigenvalues(3) * leftEigenvectors(3,j) +  &
                   deltaRightEigenvectors(i,4,:) * eigenvalues(4) * leftEigenvectors(4,j) +  &
                   rightEigenvectors(i,1) * deltaEigenvalues(1,:) * leftEigenvectors(1,j) +  &
                   rightEigenvectors(i,2) * deltaEigenvalues(2,:) * leftEigenvectors(2,j) +  &
                   rightEigenvectors(i,3) * deltaEigenvalues(3,:) * leftEigenvectors(3,j) +  &
                   rightEigenvectors(i,4) * deltaEigenvalues(4,:) * leftEigenvectors(4,j) +  &
                   rightEigenvectors(i,1) * eigenvalues(1) * deltaLeftEigenvectors(1,j,:) +  &
                   rightEigenvectors(i,2) * eigenvalues(2) * deltaLeftEigenvectors(2,j,:) +  &
                   rightEigenvectors(i,3) * eigenvalues(3) * deltaLeftEigenvectors(3,j,:) +  &
                   rightEigenvectors(i,4) * eigenvalues(4) * deltaLeftEigenvectors(4,j,:)
              do k = 1, nSpecies
                 deltaIncomingJacobianOfInviscidFlux(i,j,:) =                                &
                      deltaIncomingJacobianOfInviscidFlux(i,j,:) +                           &
                      deltaRightEigenvectors(i,4+k,:) * eigenvalues(4+k) *                   &
                      leftEigenvectors(4+k,j) + rightEigenvectors(i,4+k) *                   &
                      deltaEigenvalues(4+k,:) * leftEigenvectors(4+k,j) +                    &
                      rightEigenvectors(i,4+k) * eigenvalues(4+k) *                          &
                      deltaLeftEigenvectors(4+k,j,:)
              end do
           end do
        end do

     end if

  case (3)

     ! Eigenvalues.
     eigenvalues(1) = contravariantVelocity
     eigenvalues(2) = contravariantVelocity
     eigenvalues(3) = contravariantVelocity
     eigenvalues(4) = contravariantVelocity + speedOfSound
     eigenvalues(5) = contravariantVelocity - speedOfSound
     do k = 1, nSpecies
        eigenvalues(5+k) = contravariantVelocity
     end do
     eigenvalues = arcLength * eigenvalues

     if (present(deltaIncomingJacobianOfInviscidFlux)) then

        ! If not specified, use identity matrix for the variation of conservedVariables.
        if (present(deltaConservedVariables)) then
           deltaConservedVariables_ = deltaConservedVariables
        else
           deltaConservedVariables_ = 0.0_wp
           deltaConservedVariables_(1,1) = 1.0_wp
           deltaConservedVariables_(2,2) = 1.0_wp
           deltaConservedVariables_(3,3) = 1.0_wp
           deltaConservedVariables_(4,4) = 1.0_wp
           deltaConservedVariables_(5,5) = 1.0_wp
           do k = 1, nSpecies
              deltaConservedVariables_(5+k,5+k) = 1.0_wp
           end do
        end if

        ! Compute variations of specific volume, velocity, temperature and mass fraction.
        deltaSpecificVolume = - specificVolume_ ** 2 * deltaConservedVariables_(1,:)
        deltaVelocity(1,:) = deltaSpecificVolume * conservedVariables(2) +                   &
             specificVolume_ * deltaConservedVariables_(2,:)
        deltaVelocity(2,:) = deltaSpecificVolume * conservedVariables(3) +                   &
             specificVolume_ * deltaConservedVariables_(3,:)
        deltaVelocity(3,:) = deltaSpecificVolume * conservedVariables(4) +                   &
             specificVolume_ * deltaConservedVariables_(4,:)
        deltaTemperature = ratioOfSpecificHeats *                                            &
             (deltaSpecificVolume * conservedVariables(5) +                                  &
             specificVolume_ * deltaConservedVariables_(5,:) -                               &
             (velocity_(1) * deltaVelocity(1,:) + velocity_(2) * deltaVelocity(2,:) +        &
             velocity_(3) * deltaVelocity(3,:)))
        do k = 1, nSpecies
           deltaMassFraction(k,:) = deltaSpecificVolume * conservedVariables(5+k) +          &
                specificVolume_ * deltaConservedVariables_(5+k,:)
        end do

        ! Compute variations of other dependent variables.
        deltaContravariantVelocity = normalizedMetrics(1) * deltaVelocity(1,:) +             &
             normalizedMetrics(2) * deltaVelocity(2,:) +                                     &
             normalizedMetrics(3) * deltaVelocity(3,:)
        deltaSpeedOfSound = 0.5_wp / speedOfSound *                                          &
             (ratioOfSpecificHeats - 1.0_wp) * deltaTemperature
        deltaPhiSquared = (ratioOfSpecificHeats - 1.0_wp) *                                  &
             (velocity_(1) * deltaVelocity(1,:) + velocity_(2) * deltaVelocity(2,:) +        &
             velocity_(3) * deltaVelocity(3,:))

        ! Variation of matrix containing eigenvalues.
        deltaEigenvalues(1,:) = deltaContravariantVelocity
        deltaEigenvalues(2,:) = deltaContravariantVelocity
        deltaEigenvalues(3,:) = deltaContravariantVelocity
        deltaEigenvalues(4,:) = deltaContravariantVelocity + deltaSpeedOfSound
        deltaEigenvalues(5,:) = deltaContravariantVelocity - deltaSpeedOfSound
        do k = 1, nSpecies
           deltaEigenvalues(5+k,:) = deltaContravariantVelocity
        end do
        deltaEigenvalues = arcLength * deltaEigenvalues

     end if

     ! Zero-out the eigenvalues corresponding to outgoing characteristics and corresponding
     ! variations.
     do i = 1, 5 + nSpecies
        if (incomingDirection * real(eigenvalues(i), wp) < 0.0_wp) then
           eigenvalues(i) = 0.0_wp
           if (present(deltaIncomingJacobianOfInviscidFlux)) deltaEigenvalues(i,:) = 0.0_wp
        end if
     end do

     ! Matrix whose columns are the right eigenvectors:

     rightEigenvectors = 0.0_wp

     rightEigenvectors(1,1) = normalizedMetrics(1)
     rightEigenvectors(2,1) = normalizedMetrics(1) * velocity_(1)
     rightEigenvectors(3,1) = normalizedMetrics(1) * velocity_(2) +                          &
          conservedVariables(1) * normalizedMetrics(3)
     rightEigenvectors(4,1) = normalizedMetrics(1) * velocity_(3) -                          &
          conservedVariables(1) * normalizedMetrics(2)
     rightEigenvectors(5,1) = conservedVariables(1) * (normalizedMetrics(3) * velocity_(2) - &
          normalizedMetrics(2) * velocity_(3)) + phiSquared /                                &
          (ratioOfSpecificHeats - 1.0_wp) * normalizedMetrics(1)
     do k = 1, nSpecies
        rightEigenvectors(5+k,1) = massFraction_(k) * normalizedMetrics(1)
     end do

     rightEigenvectors(1,2) = normalizedMetrics(2)
     rightEigenvectors(2,2) = normalizedMetrics(2) * velocity_(1) -                          &
          conservedVariables(1) * normalizedMetrics(3)
     rightEigenvectors(3,2) = normalizedMetrics(2) * velocity_(2)
     rightEigenvectors(4,2) = normalizedMetrics(2) * velocity_(3) +                          &
          conservedVariables(1) * normalizedMetrics(1)
     rightEigenvectors(5,2) = conservedVariables(1) * (normalizedMetrics(1) * velocity_(3) - &
          normalizedMetrics(3) * velocity_(1))                                               &
          + phiSquared / (ratioOfSpecificHeats - 1.0_wp) * normalizedMetrics(2)
     do k = 1, nSpecies
        rightEigenvectors(5+k,2) = massFraction_(k) * normalizedMetrics(2)
     end do

     rightEigenvectors(1,3) = normalizedMetrics(3)
     rightEigenvectors(2,3) = normalizedMetrics(3) * velocity_(1) +                          &
          conservedVariables(1) * normalizedMetrics(2)
     rightEigenvectors(3,3) = normalizedMetrics(3) * velocity_(2) -                          &
          conservedVariables(1) * normalizedMetrics(1)
     rightEigenvectors(4,3) = normalizedMetrics(3) * velocity_(3)
     rightEigenvectors(5,3) = conservedVariables(1) * (normalizedMetrics(2) * velocity_(1) - &
          normalizedMetrics(1) * velocity_(2))                                               &
          + phiSquared / (ratioOfSpecificHeats - 1.0_wp) * normalizedMetrics(3)
     do k = 1, nSpecies
        rightEigenvectors(5+k,3) = massFraction_(k) * normalizedMetrics(3)
     end do

     rightEigenvectors(1,4) = 1.0_wp
     rightEigenvectors(2,4) = velocity_(1) + normalizedMetrics(1) * speedOfSound
     rightEigenvectors(3,4) = velocity_(2) + normalizedMetrics(2) * speedOfSound
     rightEigenvectors(4,4) = velocity_(3) + normalizedMetrics(3) * speedOfSound
     rightEigenvectors(5,4) = temperature_ + phiSquared / (ratioOfSpecificHeats - 1.0_wp) +  &
          speedOfSound * contravariantVelocity
     do k = 1, nSpecies
        rightEigenvectors(5+k,4) = massFraction_(k)
     end do

     rightEigenvectors(1,5) = 1.0_wp
     rightEigenvectors(2,5) = velocity_(1) - normalizedMetrics(1) * speedOfSound
     rightEigenvectors(3,5) = velocity_(2) - normalizedMetrics(2) * speedOfSound
     rightEigenvectors(4,5) = velocity_(3) - normalizedMetrics(3) * speedOfSound
     rightEigenvectors(5,5) = temperature_ + phiSquared / (ratioOfSpecificHeats - 1.0_wp) -  &
          speedOfSound * contravariantVelocity
     do k = 1, nSpecies
        rightEigenvectors(5+k,5) = massFraction_(k)
     end do

     do k = 1, nSpecies
        rightEigenvectors(5+k,5+k) = 1.0_wp
     end do

     ! Matrix whose rows are the left eigenvectors:

     leftEigenvectors = 0.0_wp

     leftEigenvectors(1,1) = normalizedMetrics(1) *                                          &
          (1.0_wp - phiSquared / speedOfSound ** 2) - specificVolume_ *                      &
          (normalizedMetrics(3) * velocity_(2) - normalizedMetrics(2) * velocity_(3))
     leftEigenvectors(2,1) = normalizedMetrics(2) *                                          &
          (1.0_wp - phiSquared / speedOfSound ** 2) - specificVolume_ *                      &
          (normalizedMetrics(1) * velocity_(3) - normalizedMetrics(3) * velocity_(1))
     leftEigenvectors(3,1) = normalizedMetrics(3) *                                          &
          (1.0_wp - phiSquared / speedOfSound ** 2) - specificVolume_ *                      &
          (normalizedMetrics(2) * velocity_(1) - normalizedMetrics(1) * velocity_(2))
     leftEigenvectors(4,1) = 0.5_wp * (phiSquared / speedOfSound ** 2 -                      &
          contravariantVelocity / speedOfSound)
     leftEigenvectors(5,1) = 0.5_wp * (phiSquared / speedOfSound ** 2 +                      &
          contravariantVelocity / speedOfSound)
     do k = 1, nSpecies
        leftEigenvectors(5+k,1) = -massFraction_(k)
     end do

     leftEigenvectors(1,2) = normalizedMetrics(1) * velocity_(1) / temperature_
     leftEigenvectors(2,2) = normalizedMetrics(2) * velocity_(1) / temperature_ -            &
          specificVolume_ * normalizedMetrics(3)
     leftEigenvectors(3,2) = normalizedMetrics(3) * velocity_(1) / temperature_ +            &
          specificVolume_ * normalizedMetrics(2)
     leftEigenvectors(4,2) = - 0.5_wp * (velocity_(1) / temperature_ -                       &
          normalizedMetrics(1) / speedOfSound)
     leftEigenvectors(5,2) = - 0.5_wp * (velocity_(1) / temperature_ +                       &
          normalizedMetrics(1) / speedOfSound)
     do k = 1, nSpecies
        leftEigenvectors(5+k,2) = 0.0_wp
     end do

     leftEigenvectors(1,3) = normalizedMetrics(1) * velocity_(2) / temperature_ +            &
          specificVolume_ * normalizedMetrics(3)
     leftEigenvectors(2,3) = normalizedMetrics(2) * velocity_(2) / temperature_
     leftEigenvectors(3,3) = normalizedMetrics(3) * velocity_(2) / temperature_ -            &
          specificVolume_ * normalizedMetrics(1)
     leftEigenvectors(4,3) = - 0.5_wp * (velocity_(2) / temperature_ -                       &
          normalizedMetrics(2) / speedOfSound)
     leftEigenvectors(5,3) = - 0.5_wp * (velocity_(2) / temperature_ +                       &
          normalizedMetrics(2) / speedOfSound)
     do k = 1, nSpecies
        leftEigenvectors(5+k,3) = 0.0_wp
     end do

     leftEigenvectors(1,4) = normalizedMetrics(1) * velocity_(3) / temperature_ -            &
          specificVolume_ * normalizedMetrics(2)
     leftEigenvectors(2,4) = normalizedMetrics(2) * velocity_(3) / temperature_ +            &
          specificVolume_ * normalizedMetrics(1)
     leftEigenvectors(3,4) = normalizedMetrics(3) * velocity_(3) / temperature_
     leftEigenvectors(4,4) = - 0.5_wp * (velocity_(3) / temperature_ -                       &
          normalizedMetrics(3) / speedOfSound)
     leftEigenvectors(5,4) = - 0.5_wp * (velocity_(3) / temperature_ +                       &
          normalizedMetrics(3) / speedOfSound)
     do k = 1, nSpecies
        leftEigenvectors(5+k,4) = 0.0_wp
     end do

     leftEigenvectors(1,5) = - normalizedMetrics(1) / temperature_
     leftEigenvectors(2,5) = - normalizedMetrics(2) / temperature_
     leftEigenvectors(3,5) = - normalizedMetrics(3) / temperature_
     leftEigenvectors(4,5) = 0.5_wp / temperature_
     leftEigenvectors(5,5) = 0.5_wp / temperature_
     do k = 1, nSpecies
        leftEigenvectors(5+k,5) = 0.0_wp
     end do

     do k = 1, nSpecies
        leftEigenvectors(5+k,5+k) = 1.0_wp
     end do

     ! ``Incoming'' part.
     do j = 1, 5 + nSpecies
        do i = 1, 5 + nSpecies
           incomingJacobianOfInviscidFlux(i,j) =                                             &
                rightEigenvectors(i,1) * eigenvalues(1) * leftEigenvectors(1,j) +            &
                rightEigenvectors(i,2) * eigenvalues(2) * leftEigenvectors(2,j) +            &
                rightEigenvectors(i,3) * eigenvalues(3) * leftEigenvectors(3,j) +            &
                rightEigenvectors(i,4) * eigenvalues(4) * leftEigenvectors(4,j) +            &
                rightEigenvectors(i,5) * eigenvalues(5) * leftEigenvectors(5,j)
           do k = 1, nSpecies
              incomingJacobianOfInviscidFlux(i,j) = incomingJacobianOfInviscidFlux(i,j) +    &
                   rightEigenvectors(i,5+k) * eigenvalues(5+k) * leftEigenvectors(5+k,j)
           end do
        end do
     end do

     if (present(deltaIncomingJacobianOfInviscidFlux)) then

        ! Variation of the matrix whose columns are the right eigenvectors:

        deltaRightEigenvectors = 0.0_wp

        deltaRightEigenvectors(1,1,:) = 0.0_wp
        deltaRightEigenvectors(2,1,:) = normalizedMetrics(1) * deltaVelocity(1,:)
        deltaRightEigenvectors(3,1,:) = normalizedMetrics(1) * deltaVelocity(2,:) +          &
             deltaConservedVariables_(1,:) * normalizedMetrics(3)
        deltaRightEigenvectors(4,1,:) = normalizedMetrics(1) * deltaVelocity(3,:) -          &
             deltaConservedVariables_(1,:) * normalizedMetrics(2)
        deltaRightEigenvectors(5,1,:) = deltaConservedVariables_(1,:) *                      &
             (normalizedMetrics(3) * velocity_(2) - normalizedMetrics(2) * velocity_(3)) +   &
             conservedVariables(1) * (normalizedMetrics(3) * deltaVelocity(2,:) -            &
             normalizedMetrics(2) * deltaVelocity(3,:)) + deltaPhiSquared /                  &
             (ratioOfSpecificHeats - 1.0_wp) * normalizedMetrics(1)
        do k = 1, nSpecies
           deltaRightEigenvectors(5+k,1,:) = 0.0_wp
        end do

        deltaRightEigenvectors(1,2,:) = 0.0_wp
        deltaRightEigenvectors(2,2,:) = normalizedMetrics(2) * deltaVelocity(1,:) -          &
             deltaConservedVariables_(1,:) * normalizedMetrics(3)
        deltaRightEigenvectors(3,2,:) = normalizedMetrics(2) * deltaVelocity(2,:)
        deltaRightEigenvectors(4,2,:) = normalizedMetrics(2) * deltaVelocity(3,:) +          &
             deltaConservedVariables_(1,:) * normalizedMetrics(1)
        deltaRightEigenvectors(5,2,:) = deltaConservedVariables_(1,:) *                      &
             (normalizedMetrics(1) * velocity_(3) - normalizedMetrics(3) * velocity_(1)) +   &
             conservedVariables(1) * (normalizedMetrics(1) * deltaVelocity(3,:) -            &
             normalizedMetrics(3) * deltaVelocity(1,:)) + deltaPhiSquared /                  &
             (ratioOfSpecificHeats - 1.0_wp) * normalizedMetrics(2)
        do k = 1, nSpecies
           deltaRightEigenvectors(5+k,2,:) = 0.0_wp
        end do

        deltaRightEigenvectors(1,3,:) = 0.0_wp
        deltaRightEigenvectors(2,3,:) = normalizedMetrics(3) * deltaVelocity(1,:) +          &
             deltaConservedVariables_(1,:) * normalizedMetrics(2)
        deltaRightEigenvectors(3,3,:) = normalizedMetrics(3) * deltaVelocity(2,:) -          &
             deltaConservedVariables_(1,:) * normalizedMetrics(1)
        deltaRightEigenvectors(4,3,:) = normalizedMetrics(3) * deltaVelocity(3,:)
        deltaRightEigenvectors(5,3,:) = deltaConservedVariables_(1,:) *                      &
             (normalizedMetrics(2) * velocity_(1) - normalizedMetrics(1) * velocity_(2)) +   &
             conservedVariables(1) * (normalizedMetrics(2) * deltaVelocity(1,:) -            &
             normalizedMetrics(1) * deltaVelocity(2,:)) + deltaPhiSquared /                  &
             (ratioOfSpecificHeats - 1.0_wp) * normalizedMetrics(3)
        do k = 1, nSpecies
           deltaRightEigenvectors(5+k,3,:) = 0.0_wp
        end do

        deltaRightEigenvectors(1,4,:) = 0.0_wp
        deltaRightEigenvectors(2,4,:) = deltaVelocity(1,:) +                                 &
             normalizedMetrics(1) * deltaSpeedOfSound
        deltaRightEigenvectors(3,4,:) = deltaVelocity(2,:) +                                 &
             normalizedMetrics(2) * deltaSpeedOfSound
        deltaRightEigenvectors(4,4,:) = deltaVelocity(3,:) +                                 &
             normalizedMetrics(3) * deltaSpeedOfSound
        deltaRightEigenvectors(5,4,:) = deltaTemperature +                                   &
             deltaPhiSquared / (ratioOfSpecificHeats - 1.0_wp) +                             &
             deltaSpeedOfSound * contravariantVelocity +                                     &
             speedOfSound * deltaContravariantVelocity
        do k = 1, nSpecies
           deltaRightEigenvectors(5+k,4,:) = deltaMassFraction(k,:)
        end do

        deltaRightEigenvectors(1,5,:) = 0.0_wp
        deltaRightEigenvectors(2,5,:) = deltaVelocity(1,:) -                                 &
             normalizedMetrics(1) * deltaSpeedOfSound
        deltaRightEigenvectors(3,5,:) = deltaVelocity(2,:) -                                 &
             normalizedMetrics(2) * deltaSpeedOfSound
        deltaRightEigenvectors(4,5,:) = deltaVelocity(3,:) -                                 &
             normalizedMetrics(3) * deltaSpeedOfSound
        deltaRightEigenvectors(5,5,:) = deltaTemperature +                                   &
             deltaPhiSquared / (ratioOfSpecificHeats - 1.0_wp) -                             &
             deltaSpeedOfSound * contravariantVelocity -                                     &
             speedOfSound * deltaContravariantVelocity
        do k = 1, nSpecies
           deltaRightEigenvectors(5+k,5,:) = deltaMassFraction(k,:)
        end do

        do k = 1, nSpecies
           deltaRightEigenvectors(5+k,5+k,:) = 0.0_wp
        end do

        ! Variation of the matrix whose rows are the left eigenvectors:

        deltaLeftEigenvectors = 0.0_wp

        temp = deltaPhiSquared / speedOfSound ** 2 -                                         &
             2.0_wp * phiSquared / speedOfSound ** 3 * deltaSpeedOfSound
        deltaLeftEigenvectors(1,1,:) = -normalizedMetrics(1) * temp - deltaSpecificVolume *  &
             (normalizedMetrics(3) * velocity_(2) - normalizedMetrics(2) * velocity_(3)) -   &
             specificVolume_ * (normalizedMetrics(3) * deltaVelocity(2,:) -                  &
             normalizedMetrics(2) * deltaVelocity(3,:))
        deltaLeftEigenvectors(2,1,:) = -normalizedMetrics(2) * temp - deltaSpecificVolume *  &
             (normalizedMetrics(1) * velocity_(3) - normalizedMetrics(3) * velocity_(1)) -   &
             specificVolume_ * (normalizedMetrics(1) * deltaVelocity(3,:) -                  &
             normalizedMetrics(3) * deltaVelocity(1,:))
        deltaLeftEigenvectors(3,1,:) = -normalizedMetrics(3) * temp - deltaSpecificVolume *  &
             (normalizedMetrics(2) * velocity_(1) - normalizedMetrics(1) * velocity_(2)) -   &
             specificVolume_ * (normalizedMetrics(2) * deltaVelocity(1,:) -                  &
             normalizedMetrics(1) * deltaVelocity(2,:))
        deltaLeftEigenvectors(4,1,:) = 0.5_wp * (temp - deltaContravariantVelocity /         &
             speedOfSound + contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)
        deltaLeftEigenvectors(5,1,:) = 0.5_wp * (temp + deltaContravariantVelocity /         &
             speedOfSound - contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)
        do k = 1, nSpecies
           deltaLeftEigenvectors(5+k,1,:) = -deltaMassFraction(k,:)
        end do

        temp = deltaVelocity(1,:) / temperature_ -                                           &
             velocity_(1) / temperature_ ** 2 * deltaTemperature
        deltaLeftEigenvectors(1,2,:) = normalizedMetrics(1) * temp
        deltaLeftEigenvectors(2,2,:) = normalizedMetrics(2) * temp -                         &
             deltaSpecificVolume * normalizedMetrics(3)
        deltaLeftEigenvectors(3,2,:) = normalizedMetrics(3) * temp +                         &
             deltaSpecificVolume * normalizedMetrics(2)
        deltaLeftEigenvectors(4,2,:) = -0.5_wp * (temp +                                     &
             normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)
        deltaLeftEigenvectors(5,2,:) = -0.5_wp * (temp -                                     &
             normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)
        do k = 1, nSpecies
           deltaLeftEigenvectors(5+k,2,:) = 0.0_wp
        end do

        temp = deltaVelocity(2,:) / temperature_ -                                           &
             velocity_(2) / temperature_ ** 2 * deltaTemperature
        deltaLeftEigenvectors(1,3,:) = normalizedMetrics(1) * temp +                         &
             deltaSpecificVolume * normalizedMetrics(3)
        deltaLeftEigenvectors(2,3,:) = normalizedMetrics(2) * temp
        deltaLeftEigenvectors(3,3,:) = normalizedMetrics(3) * temp -                         &
             deltaSpecificVolume * normalizedMetrics(1)
        deltaLeftEigenvectors(4,3,:) = -0.5_wp * (temp +                                     &
             normalizedMetrics(2) / speedOfSound ** 2 * deltaSpeedOfSound)
        deltaLeftEigenvectors(5,3,:) = -0.5_wp * (temp -                                     &
             normalizedMetrics(2) / speedOfSound ** 2 * deltaSpeedOfSound)
        do k = 1, nSpecies
           deltaLeftEigenvectors(5+k,3,:) = 0.0_wp
        end do

        temp = deltaVelocity(3,:) / temperature_ -                                           &
             velocity_(3) / temperature_ ** 2 * deltaTemperature
        deltaLeftEigenvectors(1,4,:) = normalizedMetrics(1) * temp -                         &
             deltaSpecificVolume * normalizedMetrics(2)
        deltaLeftEigenvectors(2,4,:) = normalizedMetrics(2) * temp +                         &
             deltaSpecificVolume * normalizedMetrics(1)
        deltaLeftEigenvectors(3,4,:) = normalizedMetrics(3) * temp
        deltaLeftEigenvectors(4,4,:) = -0.5_wp * (temp +                                     &
             normalizedMetrics(3) / speedOfSound ** 2 * deltaSpeedOfSound)
        deltaLeftEigenvectors(5,4,:) = -0.5_wp * (temp -                                     &
             normalizedMetrics(3) / speedOfSound ** 2 * deltaSpeedOfSound)
        do k = 1, nSpecies
           deltaLeftEigenvectors(5+k,4,:) = 0.0_wp
        end do

        temp =  -1.0_wp / temperature_ ** 2 * deltaTemperature
        deltaLeftEigenvectors(1,5,:) = -normalizedMetrics(1) * temp
        deltaLeftEigenvectors(2,5,:) = -normalizedMetrics(2) * temp
        deltaLeftEigenvectors(3,5,:) = -normalizedMetrics(3) * temp
        deltaLeftEigenvectors(4,5,:) = 0.5_wp * temp
        deltaLeftEigenvectors(5,5,:) = 0.5_wp * temp
        do k = 1, nSpecies
           deltaLeftEigenvectors(5+k,5,:) = 0.0_wp
        end do

        do k = 1, nSpecies
           deltaLeftEigenvectors(5+k,5+k,:) = 0.0_wp
        end do

        ! Variation of the ``incoming'' part.
        do j = 1, 5 + nSpecies
           do i = 1, 5 + nSpecies
              deltaIncomingJacobianOfInviscidFlux(i,j,:) =                                   &
                   deltaRightEigenvectors(i,1,:) * eigenvalues(1) * leftEigenvectors(1,j) +  &
                   deltaRightEigenvectors(i,2,:) * eigenvalues(2) * leftEigenvectors(2,j) +  &
                   deltaRightEigenvectors(i,3,:) * eigenvalues(3) * leftEigenvectors(3,j) +  &
                   deltaRightEigenvectors(i,4,:) * eigenvalues(4) * leftEigenvectors(4,j) +  &
                   deltaRightEigenvectors(i,5,:) * eigenvalues(5) * leftEigenvectors(5,j) +  &
                   rightEigenvectors(i,1) * deltaEigenvalues(1,:) * leftEigenvectors(1,j) +  &
                   rightEigenvectors(i,2) * deltaEigenvalues(2,:) * leftEigenvectors(2,j) +  &
                   rightEigenvectors(i,3) * deltaEigenvalues(3,:) * leftEigenvectors(3,j) +  &
                   rightEigenvectors(i,4) * deltaEigenvalues(4,:) * leftEigenvectors(4,j) +  &
                   rightEigenvectors(i,5) * deltaEigenvalues(5,:) * leftEigenvectors(5,j) +  &
                   rightEigenvectors(i,1) * eigenvalues(1) * deltaLeftEigenvectors(1,j,:) +  &
                   rightEigenvectors(i,2) * eigenvalues(2) * deltaLeftEigenvectors(2,j,:) +  &
                   rightEigenvectors(i,3) * eigenvalues(3) * deltaLeftEigenvectors(3,j,:) +  &
                   rightEigenvectors(i,4) * eigenvalues(4) * deltaLeftEigenvectors(4,j,:) +  &
                   rightEigenvectors(i,5) * eigenvalues(5) * deltaLeftEigenvectors(5,j,:)
              do k = 1, nSpecies
                 deltaIncomingJacobianOfInviscidFlux(i,j,:) =                                &
                      deltaIncomingJacobianOfInviscidFlux(i,j,:) +                           &
                      deltaRightEigenvectors(i,5+k,:) * eigenvalues(5+k) *                   &
                      leftEigenvectors(5+k,j) + rightEigenvectors(i,5+k) *                   &
                      deltaEigenvalues(5+k,:) * leftEigenvectors(5+k,j) +                    &
                      rightEigenvectors(i,5+k) * eigenvalues(5+k) *                          &
                      deltaLeftEigenvectors(5+k,j,:)
              end do
           end do
        end do

     end if

  end select !... nDimensions


end subroutine computeIncomingJacobianOfInviscidFlux

PURE_SUBROUTINE computeFirstPartialViscousJacobian(nDimensions, nSpecies,                    &
     conservedVariables, metrics, stressTensor, heatFlux, speciesFlux,                       &
     powerLawExponent, ratioOfSpecificHeats, firstPartialViscousJacobian,                    &
     specificVolume, velocity, temperature, massFraction)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: nDimensions, nSpecies
  SCALAR_TYPE, intent(in) :: conservedVariables(:), metrics(:),                              &
       stressTensor(:), heatFlux(:)
  real(SCALAR_KIND), intent(in) :: powerLawExponent, ratioOfSpecificHeats
  SCALAR_TYPE, intent(out) :: firstPartialViscousJacobian(:,:)
  SCALAR_TYPE, intent(in), optional :: specificVolume, velocity(:), temperature,             &
       massFraction(:), speciesFlux(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i,  k
  SCALAR_TYPE :: specificVolume_, velocity_(nDimensions), temperature_, phiSquared,          &
       contravariantStressTensor(nDimensions), contravariantHeatFlux, temp1, temp2,          &
       massFraction_(nSpecies), contravariantSpeciesFlux(nSpecies)

  ! Compute specific volume if it was not specified.
  if (present(specificVolume)) then
     specificVolume_ = specificVolume
  else
     specificVolume_ = 1.0_wp / conservedVariables(1)
  end if

  ! Compute velocity if it was not specified.
  if (present(velocity)) then
     velocity_ = velocity
  else
     do i = 1, nDimensions
        velocity_(i) = specificVolume_ * conservedVariables(i + 1)
     end do
  end if

  ! Compute temperature if it was not specified.
  if (present(temperature)) then
     temperature_ = temperature
  else
     temperature_ = ratioOfSpecificHeats * (specificVolume_ *                                &
          conservedVariables(nDimensions + 2) - 0.5_wp * sum(velocity_ ** 2))
  end if

  ! Compute mass fraction if it was not specified.
  if (present(massFraction) .and. nSpecies > 0) then
     massFraction_ = massFraction
  else
     do k = 1, nSpecies
        massFraction_(k) = specificVolume_ * conservedVariables(nDimensions + 2 + k)
     end do
  end if

  ! Zero-out first partial viscous Jacobian
  firstPartialViscousJacobian = 0.0_wp

  select case (nDimensions)
     
  case (1)

     ! Other dependent variables.
     phiSquared = 0.5_wp * (ratioOfSpecificHeats - 1.0_wp) * velocity_(1) ** 2
     contravariantStressTensor(1) = metrics(1) * stressTensor(1) !... not normalized.
     contravariantHeatFlux = metrics(1) * heatFlux(1) !... not normalized.
     do k = 1, nSpecies
        contravariantSpeciesFlux(k) = metrics(1) * speciesFlux(k,1) !... not normalized
     end do
     temp1 = velocity(1) * contravariantStressTensor(1) - contravariantHeatFlux

     temp2 = powerLawExponent * ratioOfSpecificHeats * specificVolume_ / temperature_ *      &
          (phiSquared / (ratioOfSpecificHeats - 1.0_wp) -                                    &
          temperature_ / ratioOfSpecificHeats)
     firstPartialViscousJacobian(1,1) = 0.0_wp
     firstPartialViscousJacobian(2,1) = temp2 * contravariantStressTensor(1)
     firstPartialViscousJacobian(3,1) = temp2 * temp1 - specificVolume_ *                    &
          (velocity(1) * contravariantStressTensor(1))
     do k = 1, nSpecies
        firstPartialViscousJacobian(3+k,1) = - temp2 * contravariantSpeciesFlux(k)
     end do
     

     temp2 = - powerLawExponent * ratioOfSpecificHeats *                                     &
          specificVolume_ / temperature_ * velocity(1)
     firstPartialViscousJacobian(1,2) = 0.0_wp
     firstPartialViscousJacobian(2,2) = temp2 * contravariantStressTensor(1)
     firstPartialViscousJacobian(3,2) = temp2 * temp1 +                                      &
          specificVolume_ * contravariantStressTensor(1)
     do k = 1, nSpecies
        firstPartialViscousJacobian(3+k,2) = - temp2 * contravariantSpeciesFlux(k)
     end do

     temp2 = powerLawExponent * ratioOfSpecificHeats * specificVolume_ / temperature_
     firstPartialViscousJacobian(1,3) = 0.0_wp
     firstPartialViscousJacobian(2,3) = temp2 * contravariantStressTensor(1)
     firstPartialViscousJacobian(3,3) = temp2 * temp1
     do k = 1, nSpecies
        firstPartialViscousJacobian(3+k,3) = - temp2 * contravariantSpeciesFlux(k)
     end do

  case (2)

     ! Other dependent variables.
     phiSquared = 0.5_wp * (ratioOfSpecificHeats - 1.0_wp) *                                 &
          (velocity_(1) ** 2 + velocity_(2) ** 2)
     contravariantStressTensor(1) = metrics(1) * stressTensor(1) +                           &
          metrics(2) * stressTensor(2) !... not normalized.
     contravariantStressTensor(2) = metrics(1) * stressTensor(3) +                           &
          metrics(2) * stressTensor(4) !... not normalized.
     contravariantHeatFlux = metrics(1) * heatFlux(1) +                                      &
          metrics(2) * heatFlux(2) !... not normalized.
     do k = 1, nSpecies
        contravariantSpeciesFlux(k) = metrics(1) * speciesFlux(k,1) +                        &
          metrics(2) * speciesFlux(k,2) !... not normalized
     end do
     temp1 = velocity(1) * contravariantStressTensor(1) +                                    &
          velocity(2) * contravariantStressTensor(2) - contravariantHeatFlux

     temp2 = powerLawExponent * ratioOfSpecificHeats * specificVolume_ / temperature_ *      &
          (phiSquared / (ratioOfSpecificHeats - 1.0_wp) -                                    &
          temperature_ / ratioOfSpecificHeats)
     firstPartialViscousJacobian(1,1) = 0.0_wp
     firstPartialViscousJacobian(2,1) = temp2 * contravariantStressTensor(1)
     firstPartialViscousJacobian(3,1) = temp2 * contravariantStressTensor(2)
     firstPartialViscousJacobian(4,1) = temp2 * temp1 - specificVolume_ *                    &
          (velocity(1) * contravariantStressTensor(1) +                                      &
          velocity(2) * contravariantStressTensor(2))
     do k = 1, nSpecies
        firstPartialViscousJacobian(4+k,1) = - temp2 * contravariantSpeciesFlux(k)
     end do

     temp2 = - powerLawExponent * ratioOfSpecificHeats *                                     &
          specificVolume_ / temperature_ * velocity(1)
     firstPartialViscousJacobian(1,2) = 0.0_wp
     firstPartialViscousJacobian(2,2) = temp2 * contravariantStressTensor(1)
     firstPartialViscousJacobian(3,2) = temp2 * contravariantStressTensor(2)
     firstPartialViscousJacobian(4,2) = temp2 * temp1 +                                      &
          specificVolume_ * contravariantStressTensor(1)
     do k = 1, nSpecies
        firstPartialViscousJacobian(4+k,2) = - temp2 * contravariantSpeciesFlux(k)
     end do

     temp2 = - powerLawExponent * ratioOfSpecificHeats *                                     &
          specificVolume_ / temperature_ * velocity(2)
     firstPartialViscousJacobian(1,3) = 0.0_wp
     firstPartialViscousJacobian(2,3) = temp2 * contravariantStressTensor(1)
     firstPartialViscousJacobian(3,3) = temp2 * contravariantStressTensor(2)
     firstPartialViscousJacobian(4,3) = temp2 * temp1 +                                      &
          specificVolume_ * contravariantStressTensor(2)
     do k = 1, nSpecies
        firstPartialViscousJacobian(4+k,3) = - temp2 * contravariantSpeciesFlux(k)
     end do

     temp2 = powerLawExponent * ratioOfSpecificHeats * specificVolume_ / temperature_
     firstPartialViscousJacobian(1,4) = 0.0_wp
     firstPartialViscousJacobian(2,4) = temp2 * contravariantStressTensor(1)
     firstPartialViscousJacobian(3,4) = temp2 * contravariantStressTensor(2)
     firstPartialViscousJacobian(4,4) = temp2 * temp1
     do k = 1, nSpecies
        firstPartialViscousJacobian(4+k,4) = - temp2 * contravariantSpeciesFlux(k)
     end do

  case (3)

     ! Other dependent variables.
     phiSquared = 0.5_wp * (ratioOfSpecificHeats - 1.0_wp) *                                 &
          (velocity_(1) ** 2 + velocity_(2) ** 2 + velocity_(3) ** 2)
     contravariantStressTensor(1) = metrics(1) * stressTensor(1) +                           &
          metrics(2) * stressTensor(2) + metrics(3) * stressTensor(3) !... not normalized.
     contravariantStressTensor(2) = metrics(1) * stressTensor(4) +                           &
          metrics(2) * stressTensor(5) + metrics(3) * stressTensor(6) !... not normalized.
     contravariantStressTensor(3) = metrics(1) * stressTensor(7) +                           &
          metrics(2) * stressTensor(8) + metrics(3) * stressTensor(9) !... not normalized.
     contravariantHeatFlux = metrics(1) * heatFlux(1) + metrics(2) * heatFlux(2) +           &
          metrics(3) * heatFlux(3) !... not normalized.
     do k = 1, nSpecies
        contravariantSpeciesFlux(k) = metrics(1) * speciesFlux(k,1) +                        &
             metrics(2) * speciesFlux(k,2) +                                                  &
             metrics(3) * speciesFlux(k,3) !... not normalized
     end do
     temp1 = velocity(1) * contravariantStressTensor(1) +                                    &
          velocity(2) * contravariantStressTensor(2) +                                       &
          velocity(3) * contravariantStressTensor(3) - contravariantHeatFlux

     temp2 = powerLawExponent * ratioOfSpecificHeats * specificVolume_ / temperature_ *      &
          (phiSquared / (ratioOfSpecificHeats - 1.0_wp) -                                    &
          temperature_ / ratioOfSpecificHeats)
     firstPartialViscousJacobian(1,1) = 0.0_wp
     firstPartialViscousJacobian(2,1) = temp2 * contravariantStressTensor(1)
     firstPartialViscousJacobian(3,1) = temp2 * contravariantStressTensor(2)
     firstPartialViscousJacobian(4,1) = temp2 * contravariantStressTensor(3)
     firstPartialViscousJacobian(5,1) = temp2 * temp1 - specificVolume_ *                    &
          (velocity(1) * contravariantStressTensor(1) +                                      &
          velocity(2) * contravariantStressTensor(2) +                                       &
          velocity(3) * contravariantStressTensor(3))
     do k = 1, nSpecies
        firstPartialViscousJacobian(5+k,1) = - temp2 * contravariantSpeciesFlux(k)
     end do

     temp2 = - powerLawExponent * ratioOfSpecificHeats *                                     &
          specificVolume_ / temperature_ * velocity(1)
     firstPartialViscousJacobian(1,2) = 0.0_wp
     firstPartialViscousJacobian(2,2) = temp2 * contravariantStressTensor(1)
     firstPartialViscousJacobian(3,2) = temp2 * contravariantStressTensor(2)
     firstPartialViscousJacobian(4,2) = temp2 * contravariantStressTensor(3)
     firstPartialViscousJacobian(5,2) = temp2 * temp1 +                                      &
          specificVolume_ * contravariantStressTensor(1)
     do k = 1, nSpecies
        firstPartialViscousJacobian(5+k,2) = - temp2 * contravariantSpeciesFlux(k)
     end do

     temp2 = - powerLawExponent * ratioOfSpecificHeats *                                     &
          specificVolume_ / temperature_ * velocity(2)
     firstPartialViscousJacobian(1,3) = 0.0_wp
     firstPartialViscousJacobian(2,3) = temp2 * contravariantStressTensor(1)
     firstPartialViscousJacobian(3,3) = temp2 * contravariantStressTensor(2)
     firstPartialViscousJacobian(4,3) = temp2 * contravariantStressTensor(3)
     firstPartialViscousJacobian(5,3) = temp2 * temp1 +                                      &
          specificVolume_ * contravariantStressTensor(2)
     do k = 1, nSpecies
        firstPartialViscousJacobian(5+k,3) = - temp2 * contravariantSpeciesFlux(k)
     end do

     temp2 = - powerLawExponent * ratioOfSpecificHeats *                                     &
          specificVolume_ / temperature_ * velocity(3)
     firstPartialViscousJacobian(1,4) = 0.0_wp
     firstPartialViscousJacobian(2,4) = temp2 * contravariantStressTensor(1)
     firstPartialViscousJacobian(3,4) = temp2 * contravariantStressTensor(2)
     firstPartialViscousJacobian(4,4) = temp2 * contravariantStressTensor(3)
     firstPartialViscousJacobian(5,4) = temp2 * temp1 +                                      &
          specificVolume_ * contravariantStressTensor(3)
     do k = 1, nSpecies
        firstPartialViscousJacobian(5+k,4) = - temp2 * contravariantSpeciesFlux(k)
     end do

     temp2 = powerLawExponent * ratioOfSpecificHeats * specificVolume_ / temperature_
     firstPartialViscousJacobian(1,5) = 0.0_wp
     firstPartialViscousJacobian(2,5) = temp2 * contravariantStressTensor(1)
     firstPartialViscousJacobian(3,5) = temp2 * contravariantStressTensor(2)
     firstPartialViscousJacobian(4,5) = temp2 * contravariantStressTensor(3)
     firstPartialViscousJacobian(5,5) = temp2 * temp1
     do k = 1, nSpecies
        firstPartialViscousJacobian(5+k,5) = - temp2 * contravariantSpeciesFlux(k)
     end do

  end select !... nDimensions

end subroutine computeFirstPartialViscousJacobian

PURE_SUBROUTINE computeSecondPartialViscousJacobian(nDimensions, nSpecies,                   &
     velocity, dynamicViscosity, secondCoefficientOfViscosity, thermalDiffusivity,           &
     massDiffusivity, jacobian, metricsAlongFirstDir, metricsAlongSecondDir,                 &
     secondPartialViscousJacobian)

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: nDimensions, nSpecies
  SCALAR_TYPE, intent(in) :: velocity(:), dynamicViscosity,                                  &
       secondCoefficientOfViscosity, thermalDiffusivity, jacobian,                           &
       metricsAlongFirstDir(:)
  SCALAR_TYPE, intent(in), optional :: massDiffusivity(:), metricsAlongSecondDir(:)
  SCALAR_TYPE, intent(out) :: secondPartialViscousJacobian(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: k
  SCALAR_TYPE :: temp1, temp2, temp3

  ! Zero-out second partial viscous Jacobian
  secondPartialViscousJacobian = 0.0_wp

  select case (nDimensions)

  case (1)

     ! Temporary variables.
     temp1 = metricsAlongFirstDir(1) * metricsAlongFirstDir(1)
     temp2 = dynamicViscosity * metricsAlongFirstDir(1) * velocity(1)
     temp3 = secondCoefficientOfViscosity * metricsAlongFirstDir(1) * velocity(1)

     secondPartialViscousJacobian(1,1) = dynamicViscosity * temp1 +                          &
          (dynamicViscosity + secondCoefficientOfViscosity) *                                &
          metricsAlongFirstDir(1) * metricsAlongFirstDir(1)
     secondPartialViscousJacobian(2,1) = dynamicViscosity * temp1 * velocity(1) +            &
          metricsAlongFirstDir(1) * temp2 + metricsAlongFirstDir(1) * temp3

     secondPartialViscousJacobian(1,2) = 0.0_wp
     secondPartialViscousJacobian(2,2) = thermalDiffusivity * temp1

     do k = 1, nSpecies
        secondPartialViscousJacobian(2+k,2+k) = massDiffusivity(k) * temp1
     end do

  case (2)
     
     ! Temporary variables.
     temp1 = metricsAlongFirstDir(1) * metricsAlongSecondDir(1) +                            &
          metricsAlongFirstDir(2) * metricsAlongSecondDir(2)
     temp2 = dynamicViscosity * (metricsAlongSecondDir(1) * velocity(1) +                    &
          metricsAlongSecondDir(2) * velocity(2))
     temp3 = secondCoefficientOfViscosity * (metricsAlongFirstDir(1) * velocity(1) +         &
          metricsAlongFirstDir(2) * velocity(2))

     secondPartialViscousJacobian(1,1) = dynamicViscosity * temp1 +                          &
          (dynamicViscosity + secondCoefficientOfViscosity) *                                &
          metricsAlongFirstDir(1) * metricsAlongSecondDir(1)
     secondPartialViscousJacobian(2,1) =                                                     &
          dynamicViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(2) +            &
          secondCoefficientOfViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(1)
     secondPartialViscousJacobian(3,1) = dynamicViscosity * temp1 * velocity(1) +            &
          metricsAlongFirstDir(1) * temp2 + metricsAlongSecondDir(1) * temp3

     secondPartialViscousJacobian(1,2) =                                                     &
          dynamicViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(1) +            &
          secondCoefficientOfViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(2)
     secondPartialViscousJacobian(2,2) = dynamicViscosity * temp1 +                          &
          (dynamicViscosity + secondCoefficientOfViscosity) *                                &
          metricsAlongFirstDir(2) * metricsAlongSecondDir(2)
     secondPartialViscousJacobian(3,2) = dynamicViscosity * temp1 * velocity(2) +            &
          metricsAlongFirstDir(2) * temp2 + metricsAlongSecondDir(2) * temp3

     secondPartialViscousJacobian(1,3) = 0.0_wp
     secondPartialViscousJacobian(2,3) = 0.0_wp
     secondPartialViscousJacobian(3,3) = thermalDiffusivity * temp1

     do k = 1, nSpecies
        secondPartialViscousJacobian(3+k,3+k) = massDiffusivity(k) * temp1
     end do

  case (3)

     ! Temporary variables.
     temp1 = metricsAlongFirstDir(1) * metricsAlongSecondDir(1) +                            &
          metricsAlongFirstDir(2) * metricsAlongSecondDir(2) +                               &
          metricsAlongFirstDir(3) * metricsAlongSecondDir(3)
     temp2 = dynamicViscosity * (metricsAlongSecondDir(1) * velocity(1) +                    &
          metricsAlongSecondDir(2) * velocity(2) + metricsAlongSecondDir(3) * velocity(3))
     temp3 = secondCoefficientOfViscosity * (metricsAlongFirstDir(1) * velocity(1) +         &
          metricsAlongFirstDir(2) * velocity(2) + metricsAlongFirstDir(3) * velocity(3))

     secondPartialViscousJacobian(1,1) = dynamicViscosity * temp1 +                          &
          (dynamicViscosity + secondCoefficientOfViscosity) *                                &
          metricsAlongFirstDir(1) * metricsAlongSecondDir(1)
     secondPartialViscousJacobian(2,1) =                                                     &
          dynamicViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(2) +            &
          secondCoefficientOfViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(1)
     secondPartialViscousJacobian(3,1) =                                                     &
          dynamicViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(3) +            &
          secondCoefficientOfViscosity * metricsAlongFirstDir(3) * metricsAlongSecondDir(1)
     secondPartialViscousJacobian(4,1) = dynamicViscosity * temp1 * velocity(1) +            &
          metricsAlongFirstDir(1) * temp2 + metricsAlongSecondDir(1) * temp3

     secondPartialViscousJacobian(1,2) =                                                     &
          dynamicViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(1) +            &
          secondCoefficientOfViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(2)
     secondPartialViscousJacobian(2,2) = dynamicViscosity * temp1 +                          &
          (dynamicViscosity + secondCoefficientOfViscosity) *                                &
          metricsAlongFirstDir(2) * metricsAlongSecondDir(2)
     secondPartialViscousJacobian(3,2) =                                                     &
          dynamicViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(3) +            &
          secondCoefficientOfViscosity * metricsAlongFirstDir(3) * metricsAlongSecondDir(2)
     secondPartialViscousJacobian(4,2) = dynamicViscosity * temp1 * velocity(2) +            &
          metricsAlongFirstDir(2) * temp2 + metricsAlongSecondDir(2) * temp3

     secondPartialViscousJacobian(1,3) =                                                     &
          dynamicViscosity * metricsAlongFirstDir(3) * metricsAlongSecondDir(1) +            &
          secondCoefficientOfViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(3)
     secondPartialViscousJacobian(2,3) =                                                     &
          dynamicViscosity * metricsAlongFirstDir(3) * metricsAlongSecondDir(2) +            &
          secondCoefficientOfViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(3)
     secondPartialViscousJacobian(3,3) = dynamicViscosity * temp1 +                          &
          (dynamicViscosity + secondCoefficientOfViscosity) *                                &
          metricsAlongFirstDir(3) * metricsAlongSecondDir(3)
     secondPartialViscousJacobian(4,3) = dynamicViscosity * temp1 * velocity(3) +            &
          metricsAlongFirstDir(3) * temp2 + metricsAlongSecondDir(3) * temp3

     secondPartialViscousJacobian(1,4) = 0.0_wp
     secondPartialViscousJacobian(2,4) = 0.0_wp
     secondPartialViscousJacobian(3,4) = 0.0_wp
     secondPartialViscousJacobian(4,4) = thermalDiffusivity * temp1

     do k = 1, nSpecies
        secondPartialViscousJacobian(4+k,4+k) = massDiffusivity(k) * temp1
     end do

  end select !... nDimensions

  ! Multiply by the Jacobian.
  secondPartialViscousJacobian = jacobian * secondPartialViscousJacobian

end subroutine computeSecondPartialViscousJacobian

PURE_SUBROUTINE computeJacobianOfSource(nDimensions, nSpecies,                               &
     conservedVariables, ratioOfSpecificHeats, combustion, jacobianOfSource,                 &
     specificVolume, velocity, temperature, massFraction)

  ! <<< Derived types >>>
  use Combustion_mod, only : t_Combustion

  ! <<< Arguments >>>
  integer, intent(in) :: nDimensions, nSpecies
  SCALAR_TYPE, intent(in) :: conservedVariables(:)
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
  SCALAR_TYPE, intent(out) :: jacobianOfSource(:,:)
  SCALAR_TYPE, intent(in), optional :: specificVolume, velocity(:), temperature,             &
       massFraction(:)
  type(t_Combustion), intent(in) :: combustion

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: k, l
  real(SCALAR_KIND) :: referenceTemperature, flameTemperature, activationTemperature,        &
       chemicalSource, H
  SCALAR_TYPE :: specificVolume_, velocity_(nDimensions), temperature_,                      &
       massFraction_(nSpecies), temp

  assert_key(nDimensions, (1, 2, 3))
  assert(nSpecies >= 0)
  assert(combustion%nReactions > 0)
  assert(size(conservedVariables) == nDimensions + 2 + nSpecies)

  ! Compute specific volume if it was not specified.
  if (present(specificVolume)) then
     specificVolume_ = specificVolume
  else
     specificVolume_ = 1.0_wp / conservedVariables(1)
  end if

  ! Compute velocity if it was not specified.
  if (present(velocity)) then
     assert(size(velocity) == nDimensions)
     velocity_ = velocity
  else
     do k = 1, nDimensions
        velocity_(k) = specificVolume_ * conservedVariables(k+1)
     end do
  end if

  ! Compute temperature if it was not specified.
  if (present(temperature)) then
     temperature_ = temperature
  else
     temperature_ = ratioOfSpecificHeats * (specificVolume_ * conservedVariables(nDimensions + 2) &
          - 0.5_wp * (sum(velocity_ ** 2)))
  end if

  ! Compute mass fraction if it was not specified.
  if (present(massFraction) .and. nSpecies > 0) then
     assert(size(massFraction) == nSpecies)
     massFraction_ = massFraction
  else
     do k = 1, nSpecies
        massFraction_(k) = specificVolume_ * conservedVariables(nDimensions+2+k)
     end do
  end if

  ! Bound the mass fractions.
  do k = 1, nSpecies
     !massFraction_(k) = max(massFraction_(k), 0.0_wp)
     !massFraction_(k) = min(massFraction_(k), 1.0_wp)
  end do

  ! Other dependent variables.
  referenceTemperature = 1.0_wp / (ratioOfSpecificHeats - 1.0_wp)
  flameTemperature = referenceTemperature / (1.0_wp - combustion%heatRelease)
  activationTemperature = combustion%zelDovich / combustion%heatRelease * flameTemperature
  chemicalSource = combustion%Damkohler * conservedVariables(1) *                            &
       massFraction_(combustion%H2) * massFraction_(combustion%O2) *                         &
       exp(- activationTemperature / temperature_)
  H = combustion%heatRelease * flameTemperature / combustion%Yfs

  ! Zero-out source Jacobian.
  jacobianOfSource = 0.0_wp

  jacobianOfSource(nDimensions+2,1) = H * chemicalSource * specificVolume_
  do k = 1, nSpecies
     jacobianOfSource(nDimensions+2+k,1) = - combustion%stoichiometricCoefficient(k) *       &
          chemicalSource * specificVolume_
  end do

  temp = activationTemperature / temperature_**2
  jacobianOfSource(nDimensions+2,nDimensions+2) = H * chemicalSource * temp
  do k = 1, nSpecies
     jacobianOfSource(nDimensions+2+k,nDimensions+2) =                                       &
          - combustion%stoichiometricCoefficient(k) * chemicalSource * temp
  end do

  do k = 1, nSpecies
     !if (massFraction_(k) <= 0.0_wp .or. massFraction_(k) > 1.0_wp) cycle
     !if (massFraction_(k) <= 0.0_wp) cycle
     if (k == combustion%H2) then
        temp = combustion%Damkohler * conservedVariables(1) * massFraction_(combustion%O2) * &
             exp(- activationTemperature / temperature_)
     else if (k == combustion%O2) then
        temp = combustion%Damkohler * conservedVariables(1) * massFraction_(combustion%H2) * &
             exp(- activationTemperature / temperature_)
     end if
     jacobianOfSource(nDimensions+2,nDimensions+2+k) = H * temp
     do l = 1, nSpecies
        jacobianOfSource(nDimensions+2+l,nDimensions+2+k) =                                  &
             - combustion%stoichiometricCoefficient(l) * temp
     end do
  end do

end subroutine computeJacobianOfSource
