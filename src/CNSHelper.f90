#include "config.h"

module CNSHelper

  implicit none
  private

  interface computeIncomingInviscidJacobian
     module procedure :: computeIncomingInviscidJacobianEqual,                               &
          computeIncomingInviscidJacobianRoe
  end interface computeIncomingInviscidJacobian

  public :: computeDependentVariables, computeTransportVariables, computeStressTensor,       &
       computeVorticityMagnitudeAndDilatation, computeCartesianInvsicidFluxes,               &
       computeCartesianViscousFluxes, computeSpectralRadius, transformFluxes,                &
       computeCfl, computeTimeStepSize, computeInviscidJacobian,                             &
       computeIncomingInviscidJacobian, computeFirstPartialViscousJacobian,                  &
       computeSecondPartialViscousJacobian

contains

  pure subroutine computeDependentVariables(nDimensions, conservedVariables,                 &
       ratioOfSpecificHeats, specificVolume,                                                 &
       velocity, pressure, temperature)

    !> Computes the requested dependent variable(s) including specific volume, velocity,
    !> pressure and temperature from the conserved state variables.

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: nDimensions
    real(SCALAR_KIND), intent(in) :: conservedVariables(:,:)
    real(SCALAR_KIND), intent(in), optional :: ratioOfSpecificHeats
    real(SCALAR_KIND), intent(out), optional :: specificVolume(:),                           &
         velocity(:,:), pressure(:), temperature(:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: ratioOfSpecificHeats_
    integer :: i

    ratioOfSpecificHeats_ = 1.4_wp
    if (present(ratioOfSpecificHeats)) ratioOfSpecificHeats_ = ratioOfSpecificHeats

    ! Specific volume.
    if (present(specificVolume)) specificVolume = 1.0_wp / conservedVariables(:,1)

    ! Velocity.
    if (present(velocity)) then
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
       if (present(velocity)) then
          pressure = (ratioOfSpecificHeats_ - 1.0_wp) *                                      &
               (conservedVariables(:,nDimensions+2) - 0.5_wp * conservedVariables(:,1) *     &
               sum(velocity ** 2, dim = 2))
       else if (present(specificVolume)) then
          pressure = (ratioOfSpecificHeats_ - 1.0_wp) *                                      &
               (conservedVariables(:,nDimensions+2) -                                        &
               0.5_wp * sum(conservedVariables(:,2:nDimensions+1) ** 2, dim = 2) *           &
               specificVolume)
       else
          pressure = (ratioOfSpecificHeats_ - 1.0_wp) *                                      &
               (conservedVariables(:,nDimensions+2) -                                        &
               0.5_wp * sum(conservedVariables(:,2:nDimensions+1) ** 2, dim = 2) /           &
               conservedVariables(:,1))
       end if
    end if

    ! Temperature.
    if (present(temperature)) then
       if (present(pressure)) then
          if (present(specificVolume)) then
             temperature = ratioOfSpecificHeats_ *                                           &
                  pressure / (ratioOfSpecificHeats_ - 1.0_wp) * specificVolume
          else
             temperature = ratioOfSpecificHeats_ *                                           &
                  pressure / (ratioOfSpecificHeats_ - 1.0_wp) / conservedVariables(:,1)
          end if
       else
          temperature = ratioOfSpecificHeats_ * (conservedVariables(:,nDimensions+2) -       &
               0.5_wp * sum(conservedVariables(:,2:nDimensions+1) ** 2, dim = 2) /           &
               conservedVariables(:,1)) / conservedVariables(:,1)
       end if
    end if

  end subroutine computeDependentVariables

  pure subroutine computeTransportVariables(temperature, powerLawExponent,                   &
       bulkViscosityRatio, ratioOfSpecificHeats, reynoldsNumberInverse,                      &
       prandtlNumberInverse, dynamicViscosity, secondCoefficientOfViscosity,                 &
       thermalDiffusivity)

    !> Computes the requested transport coefficient(s) including the dynamic viscosity,
    !> second coefficient of viscosity and the thermal conductivity assuming a power law
    !> dependence on temperature with exponent `powerLawExponent`.

    implicit none

    ! <<< Arguments >>>
    real(SCALAR_KIND), intent(in) :: temperature(:)
    real(SCALAR_KIND), intent(in) :: powerLawExponent, ratioOfSpecificHeats,                 &
         reynoldsNumberInverse
    real(SCALAR_KIND), intent(in), optional :: bulkViscosityRatio, prandtlNumberInverse
    real(SCALAR_KIND), intent(out), optional :: dynamicViscosity(:),                         &
         secondCoefficientOfViscosity(:), thermalDiffusivity(:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND

    if (powerLawExponent <= 0.0_wp) then !... handle powerLawExponent = 0 separately.

       ! Dynamic viscosity.
       if (present(dynamicViscosity)) dynamicViscosity = reynoldsNumberInverse

       ! Second coefficient of viscosity.
       if (present(secondCoefficientOfViscosity)) secondCoefficientOfViscosity =             &
            (bulkViscosityRatio - 2.0_wp / 3.0_wp) * reynoldsNumberInverse

       ! Thermal diffusivity.
       if (present(thermalDiffusivity))                                                      &
            thermalDiffusivity = reynoldsNumberInverse * prandtlNumberInverse

    else

       ! Dynamic viscosity.
       if (present(dynamicViscosity)) dynamicViscosity = ((ratioOfSpecificHeats - 1.0_wp) *  &
            temperature) ** powerLawExponent * reynoldsNumberInverse

       ! Second coefficient of viscosity.
       if (present(secondCoefficientOfViscosity)) then
          if (present(dynamicViscosity)) then
             secondCoefficientOfViscosity =                                                  &
                  (bulkViscosityRatio - 2.0_wp / 3.0_wp) * dynamicViscosity
          else
             secondCoefficientOfViscosity = (bulkViscosityRatio - 2.0_wp / 3.0_wp) *         &
                  ((ratioOfSpecificHeats - 1.0_wp) * temperature) ** powerLawExponent *      &
                  reynoldsNumberInverse
          end if
       end if

       ! Thermal diffusivity.
       if (present(thermalDiffusivity)) then
          if (present(dynamicViscosity)) then
             thermalDiffusivity = dynamicViscosity * prandtlNumberInverse
          else
             thermalDiffusivity =                                                            &
                  ((ratioOfSpecificHeats - 1.0_wp) * temperature) ** powerLawExponent *      &
                  reynoldsNumberInverse * prandtlNumberInverse
          end if
       end if

    end if

  end subroutine computeTransportVariables

  pure subroutine computeLocalRoeAverage(nDimensions,                                        &
       conservedVariablesL, conservedVariablesR, ratioOfSpecificHeats, roeAverage,           &
       deltaRoeAverage, deltaConservedVariablesL)

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: nDimensions
    real(SCALAR_KIND), intent(in) :: conservedVariablesL(:), conservedVariablesR(:)
    real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
    real(SCALAR_KIND), intent(out) :: roeAverage(:)
    real(SCALAR_KIND), intent(out), optional :: deltaRoeAverage(:,:)
    real(SCALAR_KIND), intent(in), optional :: deltaConservedVariablesL(:,:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, nUnknowns
    real(wp) :: sqrtDensityL, sqrtDensityR, specificVolumeL, specificVolumeR,                &
         enthalpyL, enthalpyR
    real(wp), allocatable :: deltaConservedVariablesL_(:,:), deltaSqrtDensity(:),            &
         deltaSpecificVolume(:), deltaEnthalpy(:)

    nUnknowns = size(conservedVariablesL)

    if (present(deltaRoeAverage)) then

       allocate(deltaConservedVariablesL_(nUnknowns, nUnknowns))
       allocate(deltaSqrtDensity(nUnknowns))
       allocate(deltaSpecificVolume(nUnknowns))
       allocate(deltaEnthalpy(nUnknowns))

       if (present(deltaConservedVariablesL)) then
          deltaConservedVariablesL_ = deltaConservedVariablesL
       else
          deltaConservedVariablesL_ = 0.0_wp
          do i = 1, nUnknowns
             deltaConservedVariablesL_(i,i) = 1.0_wp
          end do
       end if

    end if

    sqrtDensityL = sqrt(conservedVariablesL(1))
    sqrtDensityR = sqrt(conservedVariablesR(1))
    if (present(deltaRoeAverage))                                                            &
         deltaSqrtDensity = 0.5_wp / sqrtDensityL * deltaConservedVariablesL_(1,:)

    specificVolumeL = 1.0_wp / conservedVariablesL(1)
    specificVolumeR = 1.0_wp / conservedVariablesR(1)
    if (present(deltaRoeAverage))                                                            &
         deltaSpecificVolume = -specificVolumeL ** 2 * deltaConservedVariablesL_(1,:)

    enthalpyL = ratioOfSpecificHeats * conservedVariablesL(nDimensions+2) -                  &
         0.5_wp * (ratioOfSpecificHeats - 1.0_wp) * specificVolumeL *                        &
         sum(conservedVariablesL(2:nDimensions+1) ** 2)
    enthalpyR = ratioOfSpecificHeats * conservedVariablesR(nDimensions+2) -                  &
         0.5_wp * (ratioOfSpecificHeats - 1.0_wp) * specificVolumeR *                        &
         sum(conservedVariablesR(2:nDimensions+1) ** 2)
    if (present(deltaRoeAverage)) then
       deltaEnthalpy = ratioOfSpecificHeats * deltaConservedVariablesL_(nDimensions+2,:) -   &
            0.5_wp * (ratioOfSpecificHeats - 1.0_wp) * deltaSpecificVolume *                 &
            sum(conservedVariablesL(2:nDimensions+1) ** 2)
       do i = 1, nDimensions
          deltaEnthalpy = deltaEnthalpy - (ratioOfSpecificHeats - 1.0_wp) *                  &
               specificVolumeL * conservedVariablesL(i+1) * deltaConservedVariablesL_(i+1,:)
       end do
    end if

    roeAverage(1) = sqrtDensityL * sqrtDensityR
    do i = 1, nDimensions
       roeAverage(i+1) = (sqrtDensityR * conservedVariablesL(i+1) +                          &
            sqrtDensityL * conservedVariablesR(i+1)) / (sqrtDensityL + sqrtDensityR)
    end do

    roeAverage(nDimensions+2) = (sqrtDensityR * enthalpyL + sqrtDensityL * enthalpyR) /      &
         (sqrtDensityL + sqrtDensityR)

    if (present(deltaRoeAverage)) then

       deltaRoeAverage(1,:) = deltaSqrtDensity * sqrtDensityR
       do i = 1, nDimensions
          deltaRoeAverage(i+1,:) = (sqrtDensityR * deltaConservedVariablesL_(i+1,:) +        &
               deltaSqrtDensity * (conservedVariablesR(i+1) - roeAverage(i+1))) /            &
               (sqrtDensityL + sqrtDensityR)
       end do

       deltaRoeAverage(nDimensions+2,:) = (sqrtDensityR * deltaEnthalpy +                    &
            deltaSqrtDensity * (enthalpyR - roeAverage(nDimensions+2))) /                    &
            (sqrtDensityL + sqrtDensityR)

    end if

    roeAverage(nDimensions+2) = (roeAverage(nDimensions+2) +                                 &
         0.5_wp * (ratioOfSpecificHeats - 1.0_wp) / roeAverage(1) *                          &
         sum(roeAverage(2:nDimensions+1) ** 2)) / ratioOfSpecificHeats

    if (present(deltaRoeAverage)) then

       deltaRoeAverage(nDimensions+2,:) = deltaRoeAverage(nDimensions+2,:) -                 &
            0.5_wp * (ratioOfSpecificHeats - 1.0_wp) / roeAverage(1) ** 2 *                  &
            deltaRoeAverage(1,:) * sum(roeAverage(2:nDimensions+1) ** 2)
       do i = 1, nDimensions
          deltaRoeAverage(nDimensions+2,:) = deltaRoeAverage(nDimensions+2,:) +              &
               (ratioOfSpecificHeats - 1.0_wp) / roeAverage(1) * roeAverage(i+1) *           &
               deltaRoeAverage(i+1,:)
       end do
       deltaRoeAverage(nDimensions+2,:) = deltaRoeAverage(nDimensions+2,:) /                 &
            ratioOfSpecificHeats

    end if

    SAFE_DEALLOCATE(deltaEnthalpy)
    SAFE_DEALLOCATE(deltaSpecificVolume)
    SAFE_DEALLOCATE(deltaSqrtDensity)
    SAFE_DEALLOCATE(deltaConservedVariablesL_)

  end subroutine computeLocalRoeAverage

  pure subroutine computeStressTensor(nDimensions, velocityGradient, dynamicViscosity,       &
       secondCoefficientOfViscosity, stressTensor)

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: nDimensions
    real(SCALAR_KIND), intent(inout) :: velocityGradient(:,:)
    real(SCALAR_KIND), intent(in) :: dynamicViscosity(:), secondCoefficientOfViscosity(:)
    real(SCALAR_KIND), intent(out), optional :: stressTensor(:,:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp), allocatable :: divergenceOfVelocity(:)

    if (present(stressTensor)) then !... out-of-place computation.

       select case (nDimensions)

       case (1)
          stressTensor(:,1) = (2.0_wp * dynamicViscosity +                                   &
               secondCoefficientOfViscosity) * velocityGradient(:,1)

       case (2)
          stressTensor(:,1) = 2.0_wp * dynamicViscosity * velocityGradient(:,1) +            &
               secondCoefficientOfViscosity * (velocityGradient(:,1) + velocityGradient(:,4))
          stressTensor(:,2) = dynamicViscosity * (velocityGradient(:,2) +                    &
               velocityGradient(:,3))
          stressTensor(:,3) = stressTensor(:,2)
          stressTensor(:,4) = 2.0_wp * dynamicViscosity * velocityGradient(:,4) +            &
               secondCoefficientOfViscosity * (velocityGradient(:,1) + velocityGradient(:,4))

       case (3)
          stressTensor(:,1) = 2.0_wp * dynamicViscosity * velocityGradient(:,1) +            &
               secondCoefficientOfViscosity * (velocityGradient(:,1) +                       &
               velocityGradient(:,5) + velocityGradient(:,9))
          stressTensor(:,2) = dynamicViscosity * (velocityGradient(:,2) +                    &
               velocityGradient(:,4))
          stressTensor(:,3) = dynamicViscosity * (velocityGradient(:,3) +                    &
               velocityGradient(:,7))
          stressTensor(:,4) = stressTensor(:,2)
          stressTensor(:,5) = 2.0_wp * dynamicViscosity * velocityGradient(:,5) +            &
               secondCoefficientOfViscosity * (velocityGradient(:,1) +                       &
               velocityGradient(:,5) + velocityGradient(:,9))
          stressTensor(:,6) = dynamicViscosity * (velocityGradient(:,6) +                    &
               velocityGradient(:,8))
          stressTensor(:,7) = stressTensor(:,3)
          stressTensor(:,8) = stressTensor(:,6)
          stressTensor(:,9) = 2.0_wp * dynamicViscosity * velocityGradient(:,9) +            &
               secondCoefficientOfViscosity * (velocityGradient(:,1) +                       &
               velocityGradient(:,5) + velocityGradient(:,9))

       end select

    else !... in-place computation.

       select case (nDimensions)

       case (1)
          velocityGradient(:,1) = (2.0_wp * dynamicViscosity +                               &
               secondCoefficientOfViscosity) * velocityGradient(:,1)

       case (2)
          allocate(divergenceOfVelocity(size(velocityGradient, 1)))
          divergenceOfVelocity = secondCoefficientOfViscosity *                              &
               (velocityGradient(:,1) + velocityGradient(:,4))
          velocityGradient(:,1) = 2.0_wp * dynamicViscosity * velocityGradient(:,1) +        &
               divergenceOfVelocity
          velocityGradient(:,2) = dynamicViscosity *                                         &
               (velocityGradient(:,2) + velocityGradient(:,3))
          velocityGradient(:,3) = velocityGradient(:,2)
          velocityGradient(:,4) = 2.0_wp * dynamicViscosity * velocityGradient(:,4) +        &
               divergenceOfVelocity

       case (3)
          allocate(divergenceOfVelocity(size(velocityGradient, 1)))
          divergenceOfVelocity = secondCoefficientOfViscosity *                              &
               (velocityGradient(:,1) + velocityGradient(:,5) + velocityGradient(:,9))
          velocityGradient(:,1) = 2.0_wp * dynamicViscosity * velocityGradient(:,1) +        &
               divergenceOfVelocity
          velocityGradient(:,2) = dynamicViscosity *                                         &
               (velocityGradient(:,2) + velocityGradient(:,4))
          velocityGradient(:,3) = dynamicViscosity *                                         &
               (velocityGradient(:,3) + velocityGradient(:,7))
          velocityGradient(:,4) = velocityGradient(:,2)
          velocityGradient(:,5) = 2.0_wp * dynamicViscosity * velocityGradient(:,5) +        &
               divergenceOfVelocity
          velocityGradient(:,6) = dynamicViscosity *                                         &
               (velocityGradient(:,6) + velocityGradient(:,8))
          velocityGradient(:,7) = velocityGradient(:,3)
          velocityGradient(:,8) = velocityGradient(:,6)
          velocityGradient(:,9) = 2.0_wp * dynamicViscosity * velocityGradient(:,9) +        &
               divergenceOfVelocity

       end select

    end if

    SAFE_DEALLOCATE(divergenceOfVelocity)

  end subroutine computeStressTensor

  pure subroutine computeVorticityMagnitudeAndDilatation(nDimensions,                        &
       velocityGradient, vorticityMagnitude, dilatation)

    !> Computes the vorticity magnitude and dilatation from the velocity gradient.

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: nDimensions
    real(SCALAR_KIND), intent(in) :: velocityGradient(:,:)
    real(SCALAR_KIND), intent(out), optional :: vorticityMagnitude(:), dilatation(:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND

    select case (nDimensions)

    case (1)
       if (present(dilatation)) dilatation = velocityGradient(:,1)
       if (present(vorticityMagnitude)) vorticityMagnitude = 0.0_wp

    case (2)
       if (present(dilatation)) dilatation = velocityGradient(:,1) + velocityGradient(:,4)
       if (present(vorticityMagnitude))                                                      &
            vorticityMagnitude = abs(velocityGradient(:,3) - velocityGradient(:,2))

    case (3)
       if (present(dilatation))                                                              &
            dilatation = velocityGradient(:,1) + velocityGradient(:,5) + velocityGradient(:,9)
       if (present(vorticityMagnitude))                                                      &
            vorticityMagnitude = sqrt(abs(velocityGradient(:,8) -                            &
            velocityGradient(:,6)) ** 2 + abs(velocityGradient(:,3) -                        &
            velocityGradient(:,7)) ** 2 + abs(velocityGradient(:,4) -                        &
            velocityGradient(:,2)) ** 2)

    end select

  end subroutine computeVorticityMagnitudeAndDilatation

  pure subroutine computeCartesianInvsicidFluxes(nDimensions, conservedVariables,            &
       velocity, pressure, inviscidFluxes)

    !> Computes the Cartesian form of the inviscid fluxes. `velocity` and `pressure` must
    !> be consistent with `conservedVariables`, otherwise the result is unpredictable.

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: nDimensions
    real(SCALAR_KIND), intent(in) :: conservedVariables(:,:), velocity(:,:), pressure(:)
    real(SCALAR_KIND), intent(out) :: inviscidFluxes(:,:,:)

    select case (nDimensions)

    case (1)
       inviscidFluxes(:,1,1) = conservedVariables(:,2)
       inviscidFluxes(:,2,1) = conservedVariables(:,2) * velocity(:,1) + pressure
       inviscidFluxes(:,3,1) = velocity(:,1) * (conservedVariables(:,3) + pressure)

    case (2)
       inviscidFluxes(:,1,1) = conservedVariables(:,2)
       inviscidFluxes(:,2,1) = conservedVariables(:,2) * velocity(:,1) + pressure
       inviscidFluxes(:,3,1) = conservedVariables(:,2) * velocity(:,2)
       inviscidFluxes(:,4,1) = velocity(:,1) * (conservedVariables(:,4) + pressure)
       inviscidFluxes(:,1,2) = conservedVariables(:,3)
       inviscidFluxes(:,2,2) = inviscidFluxes(:,3,1)
       inviscidFluxes(:,3,2) = conservedVariables(:,3) * velocity(:,2) + pressure
       inviscidFluxes(:,4,2) = velocity(:,2) * (conservedVariables(:,4) + pressure)

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

    end select

  end subroutine computeCartesianInvsicidFluxes

  pure subroutine computeCartesianViscousFluxes(nDimensions, velocity,                       &
       stressTensor, heatFlux, viscousFluxes)

    !> Computes the Cartesian form of the viscous fluxes.

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: nDimensions
    real(SCALAR_KIND), intent(in) :: velocity(:,:), stressTensor(:,:), heatFlux(:,:)
    real(SCALAR_KIND), intent(out) :: viscousFluxes(:,:,:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND

    select case (nDimensions)

    case (1)
       viscousFluxes(:,1,1) = 0.0_wp
       viscousFluxes(:,2,1) = stressTensor(:,1)
       viscousFluxes(:,3,1) = velocity(:,1) * stressTensor(:,1) - heatFlux(:,1)

    case (2)
       viscousFluxes(:,1,1) = 0.0_wp
       viscousFluxes(:,2,1) = stressTensor(:,1)
       viscousFluxes(:,3,1) = stressTensor(:,3)
       viscousFluxes(:,4,1) = velocity(:,1) * stressTensor(:,1) +                            &
            velocity(:,2) * stressTensor(:,3) - heatFlux(:,1)
       viscousFluxes(:,1,2) = 0.0_wp
       viscousFluxes(:,2,2) = stressTensor(:,2)
       viscousFluxes(:,3,2) = stressTensor(:,4)
       viscousFluxes(:,4,2) = velocity(:,1) * stressTensor(:,2) +                            &
            velocity(:,2) * stressTensor(:,4) - heatFlux(:,2)

    case (3)
       viscousFluxes(:,1,1) = 0.0_wp
       viscousFluxes(:,2,1) = stressTensor(:,1)
       viscousFluxes(:,3,1) = stressTensor(:,4)
       viscousFluxes(:,4,1) = stressTensor(:,7)
       viscousFluxes(:,5,1) = velocity(:,1) * stressTensor(:,1) +                            &
            velocity(:,2) * stressTensor(:,4) +                                              &
            velocity(:,3) * stressTensor(:,7) - heatFlux(:,1)
       viscousFluxes(:,1,2) = 0.0_wp
       viscousFluxes(:,2,2) = stressTensor(:,2)
       viscousFluxes(:,3,2) = stressTensor(:,5)
       viscousFluxes(:,4,2) = stressTensor(:,8)
       viscousFluxes(:,5,2) = velocity(:,1) * stressTensor(:,2) +                            &
            velocity(:,2) * stressTensor(:,5) +                                              &
            velocity(:,3) * stressTensor(:,8) - heatFlux(:,2)
       viscousFluxes(:,1,3) = 0.0_wp
       viscousFluxes(:,2,3) = stressTensor(:,3)
       viscousFluxes(:,3,3) = stressTensor(:,6)
       viscousFluxes(:,4,3) = stressTensor(:,9)
       viscousFluxes(:,5,3) = velocity(:,1) * stressTensor(:,3) +                            &
            velocity(:,2) * stressTensor(:,6) +                                              &
            velocity(:,3) * stressTensor(:,9) - heatFlux(:,3)

    end select

  end subroutine computeCartesianViscousFluxes

  pure subroutine computeSpectralRadius(nDimensions, ratioOfSpecificHeats, velocity,         &
       temperature, metrics, spectralRadius, isDomainCurvilinear)

    !> Compute the spectral radii along all the directions.

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: nDimensions
    real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
    real(SCALAR_KIND), intent(in) :: velocity(:,:), temperature(:), metrics(:,:)
    real(SCALAR_KIND), intent(out) :: spectralRadius(:,:)
    logical, intent(in), optional :: isDomainCurvilinear

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    logical :: isDomainCurvilinear_

    isDomainCurvilinear_ = .true.
    if (present(isDomainCurvilinear)) isDomainCurvilinear_ = isDomainCurvilinear

    ! Temporary storage for speed of sound.
    spectralRadius(:,nDimensions) = sqrt((ratioOfSpecificHeats - 1.0_wp) * temperature)

    select case (nDimensions)

    case (1)
       spectralRadius(:,1) = abs(metrics(:,1) * velocity(:,1)) +                             &
            spectralRadius(:,1) * abs(metrics(:,1))

    case (2)
       if (isDomainCurvilinear_) then
          spectralRadius(:,1) = abs(metrics(:,1) * velocity(:,1) +                           &
               metrics(:,2) * velocity(:,2)) +                                               &
               spectralRadius(:,2) * sqrt(metrics(:,1) ** 2 + metrics(:,2) ** 2)
          spectralRadius(:,2) = abs(metrics(:,3) * velocity(:,1) +                           &
               metrics(:,4) * velocity(:,2)) +                                               &
               spectralRadius(:,2) * sqrt(metrics(:,3) ** 2 + metrics(:,4) ** 2)
       else
          spectralRadius(:,1) = abs(metrics(:,1) * velocity(:,1)) +                          &
               spectralRadius(:,2) * abs(metrics(:,1))
          spectralRadius(:,2) = abs(metrics(:,4) * velocity(:,2)) +                          &
               spectralRadius(:,2) * abs(metrics(:,4))
       end if

    case (3)
       if (isDomainCurvilinear_) then
          spectralRadius(:,1) = abs(metrics(:,1) * velocity(:,1) +                           &
               metrics(:,2) * velocity(:,2) + metrics(:,3) * velocity(:,3)) +                &
               spectralRadius(:,3) * sqrt(metrics(:,1) ** 2 +                                &
               metrics(:,2) ** 2 + metrics(:,3) ** 2)
          spectralRadius(:,2) = abs(metrics(:,4) * velocity(:,1) +                           &
               metrics(:,5) * velocity(:,2) + metrics(:,6) * velocity(:,3)) +                &
               spectralRadius(:,3) * sqrt(metrics(:,4) ** 2 +                                &
               metrics(:,5) ** 2 + metrics(:,6) ** 2)
          spectralRadius(:,3) = abs(metrics(:,7) * velocity(:,1) +                           &
               metrics(:,8) * velocity(:,2) + metrics(:,9) * velocity(:,3)) +                &
               spectralRadius(:,3) * sqrt(metrics(:,7) ** 2 +                                &
               metrics(:,8) ** 2 + metrics(:,9) ** 2)
       else
          spectralRadius(:,1) = abs(metrics(:,1) * velocity(:,1)) +                          &
               spectralRadius(:,3) * abs(metrics(:,1))
          spectralRadius(:,2) = abs(metrics(:,5) * velocity(:,2)) +                          &
               spectralRadius(:,3) * abs(metrics(:,5))
          spectralRadius(:,3) = abs(metrics(:,9) * velocity(:,3)) +                          &
               spectralRadius(:,3) * abs(metrics(:,9))
       end if

    end select

  end subroutine computeSpectralRadius

  pure subroutine transformFluxes(nDimensions, fluxes, metrics,                              &
       transformedFluxes, isDomainCurvilinear)

    !> Transforms fluxes from Cartesian form to contravariant form. If the domain is not
    !> curvilinear as specified by `isDomainCurvilinear`, the arguments `fluxes` and
    !> `transformedFluxes` may reference the same array.

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: nDimensions
    real(SCALAR_KIND), intent(in) :: fluxes(:,:,:), metrics(:,:)
    real(SCALAR_KIND), intent(out) :: transformedFluxes(:,:,:)
    logical, intent(in), optional :: isDomainCurvilinear

    ! <<< Local variables >>>
    logical :: isDomainCurvilinear_
    integer :: i

    isDomainCurvilinear_ = .true.
    if (present(isDomainCurvilinear)) isDomainCurvilinear_ = isDomainCurvilinear

    select case (nDimensions)

    case (1)
       do i = 1, size(fluxes, 2)
          transformedFluxes(:,i,1) = metrics(:,1) * fluxes(:,i,1)
       end do

    case (2)
       if (isDomainCurvilinear_) then
          do i = 1, size(fluxes, 2)
             transformedFluxes(:,i,1) = metrics(:,1) * fluxes(:,i,1) +                       &
                  metrics(:,2) * fluxes(:,i,2)
             transformedFluxes(:,i,2) = metrics(:,3) * fluxes(:,i,1) +                       &
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
             transformedFluxes(:,i,1) = metrics(:,1) * fluxes(:,i,1) +                       &
                  metrics(:,2) * fluxes(:,i,2) + metrics(:,3) * fluxes(:,i,3)
             transformedFluxes(:,i,2) = metrics(:,4) * fluxes(:,i,1) +                       &
                  metrics(:,5) * fluxes(:,i,2) + metrics(:,6) * fluxes(:,i,3)
             transformedFluxes(:,i,3) = metrics(:,7) * fluxes(:,i,1) +                       &
                  metrics(:,8) * fluxes(:,i,2) + metrics(:,9) * fluxes(:,i,3)
          end do
       else
          do i = 1, size(fluxes, 2)
             transformedFluxes(:,i,1) = metrics(:,1) * fluxes(:,i,1)
             transformedFluxes(:,i,2) = metrics(:,5) * fluxes(:,i,2)
             transformedFluxes(:,i,3) = metrics(:,9) * fluxes(:,i,3)
          end do
       end if

    end select

  end subroutine transformFluxes

  pure function computeCfl(nDimensions, iblank, jacobian, metrics,                           &
       velocity, temperature, timeStepSize, ratioOfSpecificHeats,                            &
       dynamicViscosity, thermalDiffusivity) result(cfl)

    !> Computes the CFL number.

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: nDimensions, iblank(:)
    real(SCALAR_KIND), intent(in) :: jacobian(:), metrics(:,:), velocity(:,:), temperature(:)
    real(SCALAR_KIND), intent(in) :: timeStepSize, ratioOfSpecificHeats
    real(SCALAR_KIND), intent(in), optional :: dynamicViscosity(:), thermalDiffusivity(:)

    ! <<< Result >>>
    real(SCALAR_KIND) :: cfl

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j
    real(wp) :: localSpeedOfSound, localWaveSpeed

    cfl = 0.0_wp

    ! Advection.
    do i = 1, size(iblank)
       if (iblank(i) == 0) cycle ! ... skip hole points.
       localSpeedOfSound = sqrt((ratioOfSpecificHeats - 1.0_wp) * temperature(i))
       localWaveSpeed = 0.0_wp
       do j = 1, nDimensions
          localWaveSpeed = localWaveSpeed +                                                  &
               sum(metrics(i,1+nDimensions*(j-1):nDimensions*j) ** 2)
       end do
       localWaveSpeed = localSpeedOfSound * sqrt(localWaveSpeed)
       do j = 1, nDimensions
          localWaveSpeed = localWaveSpeed + abs(dot_product(velocity(i,:),                   &
               metrics(i,1+nDimensions*(j-1):nDimensions*j)))
       end do
       localWaveSpeed = jacobian(i) * localWaveSpeed
       cfl = max(cfl, localWaveSpeed * timeStepSize)
    end do

    ! Diffusion.
    if (present(dynamicViscosity)) then
       do i = 1, size(iblank)
          if (iblank(i) == 0) cycle !... skip hole points.
          localWaveSpeed = 0.0_wp
          do j = 1, nDimensions
             localWaveSpeed = localWaveSpeed +                                               &
                  sum(metrics(i,1+nDimensions*(j-1):nDimensions*j) ** 2)
          end do
          localWaveSpeed = jacobian(i) ** 2 * sum(metrics(i,:) ** 2) *                       &
               max(2.0_wp * dynamicViscosity(i), thermalDiffusivity(i))
          cfl = max(cfl, localWaveSpeed * timeStepSize)
       end do
    end if

  end function computeCfl

  pure function computeTimeStepSize(nDimensions, iblank, jacobian,                           &
       metrics, velocity, temperature, cfl,                                                  &
       ratioOfSpecificHeats, dynamicViscosity,                                               &
       thermalDiffusivity) result(timeStepSize)

    !> Computes the time step size that leads to a specified CFL number.

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: nDimensions, iblank(:)
    real(SCALAR_KIND), intent(in) :: jacobian(:), metrics(:,:), velocity(:,:), temperature(:)
    real(SCALAR_KIND), intent(in) :: cfl, ratioOfSpecificHeats
    real(SCALAR_KIND), intent(in), optional :: dynamicViscosity(:), thermalDiffusivity(:)

    ! <<< Result >>>
    real(SCALAR_KIND) :: timeStepSize

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j
    real(wp) :: localSpeedOfSound, localWaveSpeed

    timeStepSize = huge(0.0_wp)

    ! Advection.
    do i = 1, size(iblank)
       if (iblank(i) == 0) cycle ! ... skip hole points.
       localSpeedOfSound = sqrt((ratioOfSpecificHeats - 1.0_wp) * temperature(i))
       localWaveSpeed = 0.0_wp
       do j = 1, nDimensions
          localWaveSpeed = localWaveSpeed +                                                  &
               sum(metrics(i,1+nDimensions*(j-1):nDimensions*j) ** 2)
       end do
       localWaveSpeed = localSpeedOfSound * sqrt(localWaveSpeed)
       do j = 1, nDimensions
          localWaveSpeed = localWaveSpeed + abs(dot_product(velocity(i,:),                   &
               metrics(i,1+nDimensions*(j-1):nDimensions*j)))
       end do
       localWaveSpeed = jacobian(i) * localWaveSpeed
       timeStepSize = min(timeStepSize, cfl / localWaveSpeed)
    end do

    ! Diffusion.
    if (present(dynamicViscosity)) then
       do i = 1, size(iblank)
          if (iblank(i) == 0) cycle !... skip hole points.
          localWaveSpeed = 0.0_wp
          do j = 1, nDimensions
             localWaveSpeed = localWaveSpeed +                                               &
                  sum(metrics(i,1+nDimensions*(j-1):nDimensions*j) ** 2)
          end do
          localWaveSpeed = jacobian(i) ** 2 * sum(metrics(i,:) ** 2) *                       &
               max(2.0_wp * dynamicViscosity(i), thermalDiffusivity(i))
          timeStepSize = min(timeStepSize, cfl / localWaveSpeed)
       end do
    end if

  end function computeTimeStepSize

  pure subroutine computeInviscidJacobian(nDimensions, ijkIndex, conservedVariables,         &
       metricsAlongNormalDirection, ratioOfSpecificHeats, jacobianOfInviscidFlux,            &
       deltaJacobianOfInviscidFlux, deltaConservedVariables, specificVolume, velocity,       &
       temperature)

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: nDimensions, ijkIndex
    real(SCALAR_KIND), intent(in) :: conservedVariables(:,:), metricsAlongNormalDirection(:,:)
    real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
    real(SCALAR_KIND), intent(out) :: jacobianOfInviscidFlux(:,:)
    real(SCALAR_KIND), intent(out), optional :: deltaJacobianOfInviscidFlux(:,:,:)
    real(SCALAR_KIND), intent(in), optional :: deltaConservedVariables(:,:),                 &
         specificVolume(:), velocity(:,:), temperature(:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, nUnknowns
    real(wp), allocatable :: localConservedVariables(:), localVelocity(:),                   &
         localMetricsAlongNormalDirection(:), deltaConservedVariables_(:,:),                 &
         deltaSpecificVolume(:), deltaVelocity(:,:), deltaTemperature(:)
    real(wp) :: localSpecificVolume, localTemperature

    nUnknowns = size(conservedVariables, 2)

    allocate(localConservedVariables(nUnknowns))
    allocate(localVelocity(nDimensions))
    allocate(localMetricsAlongNormalDirection(nDimensions))

    localConservedVariables = conservedVariables(ijkIndex,:)
    localMetricsAlongNormalDirection = metricsAlongNormalDirection(ijkIndex,:)

    if (present(specificVolume)) then
       localSpecificVolume = specificVolume(ijkIndex)
    else
       localSpecificVolume = 1.0_wp / localConservedVariables(1)
    end if

    if (present(velocity)) then
       localVelocity = velocity(ijkIndex,:)
    else
       localVelocity = localSpecificVolume * localConservedVariables(2:nDimensions+1)
    end if

    if (present(temperature)) then
       localTemperature = temperature(ijkIndex)
    else
       localTemperature = ratioOfSpecificHeats * (localSpecificVolume *                      &
            localConservedVariables(nDimensions+2) - 0.5_wp * sum(localVelocity ** 2))
    end if

    if (present(deltaJacobianOfInviscidFlux)) then

       if (present(deltaConservedVariables)) then
          allocate(deltaConservedVariables_(nUnknowns, nUnknowns),                           &
               source = deltaConservedVariables)
       else
          allocate(deltaConservedVariables_(nUnknowns, nUnknowns), source = 0.0_wp)
          do i = 1, nUnknowns
             deltaConservedVariables_(i,i) = 1.0_wp
          end do
       end if

       allocate(deltaSpecificVolume(nUnknowns))
       allocate(deltaVelocity(nDimensions, nUnknowns))
       allocate(deltaTemperature(nUnknowns))

       deltaSpecificVolume = -localSpecificVolume ** 2 * deltaConservedVariables_(1,:)
       do i = 1, nDimensions
          deltaVelocity(i,:) = deltaSpecificVolume * localConservedVariables(i+1) +          &
               specificVolume * deltaConservedVariables_(i+1,:)
       end do

       deltaTemperature = ratioOfSpecificHeats * (deltaSpecificVolume *                      &
            localConservedVariables(nDimensions+2) + localSpecificVolume *                   &
            deltaConservedVariables_(nDimensions+2,:) - matmul(localVelocity, deltaVelocity))

       select case (nDimensions)
       case (1)
          call computeInviscidJacobian1D(localVelocity, localTemperature,                    &
               localMetricsAlongNormalDirection, ratioOfSpecificHeats,                       &
               jacobianOfInviscidFlux, deltaJacobianOfInviscidFlux, deltaVelocity,           &
               deltaTemperature)
       case (2)
          call computeInviscidJacobian2D(localVelocity, localTemperature,                    &
               localMetricsAlongNormalDirection, ratioOfSpecificHeats,                       &
               jacobianOfInviscidFlux, deltaJacobianOfInviscidFlux, deltaVelocity,           &
               deltaTemperature)
       case (3)
          call computeInviscidJacobian3D(localVelocity, localTemperature,                    &
               localMetricsAlongNormalDirection, ratioOfSpecificHeats,                       &
               jacobianOfInviscidFlux, deltaJacobianOfInviscidFlux, deltaVelocity,           &
               deltaTemperature)
       end select

    else

       select case (nDimensions)
       case (1)
          call computeInviscidJacobian1D(localVelocity, localTemperature,                    &
               localMetricsAlongNormalDirection, ratioOfSpecificHeats, jacobianOfInviscidFlux)
       case (2)
          call computeInviscidJacobian2D(localVelocity, localTemperature,                    &
               localMetricsAlongNormalDirection, ratioOfSpecificHeats, jacobianOfInviscidFlux)
       case (3)
          call computeInviscidJacobian3D(localVelocity, localTemperature,                    &
               localMetricsAlongNormalDirection, ratioOfSpecificHeats, jacobianOfInviscidFlux)
       end select

    end if

    SAFE_DEALLOCATE(deltaTemperature)
    SAFE_DEALLOCATE(deltaVelocity)
    SAFE_DEALLOCATE(deltaSpecificVolume)
    SAFE_DEALLOCATE(deltaConservedVariables_)
    SAFE_DEALLOCATE(localMetricsAlongNormalDirection)
    SAFE_DEALLOCATE(localVelocity)
    SAFE_DEALLOCATE(localConservedVariables)

  end subroutine computeInviscidJacobian

  pure subroutine computeInviscidJacobian1D(velocity, temperature, metrics,                  &
       ratioOfSpecificHeats, jacobianOfInviscidFlux, deltaJacobianOfInviscidFlux,            &
       deltaVelocity, deltaTemperature)

    ! <<< Arguments >>>
    real(SCALAR_KIND), intent(in) :: velocity(1), temperature, metrics(1),                   &
         ratioOfSpecificHeats
    real(SCALAR_KIND), intent(out) :: jacobianOfInviscidFlux(3,3)
    real(SCALAR_KIND), intent(out), optional :: deltaJacobianOfInviscidFlux(3,3,3)
    real(SCALAR_KIND), intent(in), optional :: deltaVelocity(1,3), deltaTemperature(3)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: contravariantVelocity, phiSquared, deltaContravariantVelocity(3),            &
         deltaPhiSquared(3)

    contravariantVelocity = metrics(1) * velocity(1) !... not normalized.

    phiSquared = 0.5_wp * (ratioOfSpecificHeats - 1.0_wp) * (velocity(1) ** 2)

    jacobianOfInviscidFlux(1,1) = 0.0_wp
    jacobianOfInviscidFlux(2,1) = phiSquared * metrics(1) -                                  &
         contravariantVelocity * velocity(1)
    jacobianOfInviscidFlux(3,1) = contravariantVelocity * ((ratioOfSpecificHeats - 2.0_wp) / &
         (ratioOfSpecificHeats - 1.0_wp) * phiSquared - temperature)

    jacobianOfInviscidFlux(1,2) = metrics(1)
    jacobianOfInviscidFlux(2,2) = contravariantVelocity -                                    &
         (ratioOfSpecificHeats - 2.0_wp) * velocity(1) * metrics(1)
    jacobianOfInviscidFlux(3,2) = (temperature + phiSquared /                                &
         (ratioOfSpecificHeats - 1.0_wp)) * metrics(1) - (ratioOfSpecificHeats - 1.0_wp) *   &
         contravariantVelocity * velocity(1)

    jacobianOfInviscidFlux(1,3) = 0.0_wp
    jacobianOfInviscidFlux(2,3) = (ratioOfSpecificHeats - 1.0_wp) * metrics(1)
    jacobianOfInviscidFlux(3,3) = ratioOfSpecificHeats * contravariantVelocity

    if (present(deltaJacobianOfInviscidFlux)) then

       deltaContravariantVelocity = metrics(1) * deltaVelocity(1,:)

       deltaPhiSquared = (ratioOfSpecificHeats - 1.0_wp) * (velocity(1) * deltaVelocity(1,:))

       deltaJacobianOfInviscidFlux(1,1,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(2,1,:) = deltaPhiSquared * metrics(1) -                   &
            deltaContravariantVelocity * velocity(1) - contravariantVelocity *               &
            deltaVelocity(1,:)
       deltaJacobianOfInviscidFlux(3,1,:) = deltaContravariantVelocity *                     &
            ((ratioOfSpecificHeats - 2.0_wp) / (ratioOfSpecificHeats - 1.0_wp) *             &
            phiSquared - temperature) + contravariantVelocity *                              &
            ((ratioOfSpecificHeats - 2.0_wp) / (ratioOfSpecificHeats - 1.0_wp) *             &
            deltaPhiSquared - deltaTemperature)

       deltaJacobianOfInviscidFlux(1,2,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(2,2,:) = deltaContravariantVelocity -                     &
            (ratioOfSpecificHeats - 2.0_wp) * deltaVelocity(1,:) * metrics(1)
       deltaJacobianOfInviscidFlux(3,2,:) = (deltaTemperature + deltaPhiSquared /            &
            (ratioOfSpecificHeats - 1.0_wp)) * metrics(1) -                                  &
            (ratioOfSpecificHeats - 1.0_wp) * (deltaContravariantVelocity * velocity(1) +    &
            contravariantVelocity * deltaVelocity(1,:))

       deltaJacobianOfInviscidFlux(1,3,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(2,3,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(3,3,:) = ratioOfSpecificHeats * deltaContravariantVelocity

    end if

  end subroutine computeInviscidJacobian1D

  pure subroutine computeInviscidJacobian2D(velocity, temperature, metrics,                  &
       ratioOfSpecificHeats, jacobianOfInviscidFlux, deltaJacobianOfInviscidFlux,            &
       deltaVelocity, deltaTemperature)

    ! <<< Arguments >>>
    real(SCALAR_KIND), intent(in) :: velocity(2), temperature, metrics(2),                   &
         ratioOfSpecificHeats
    real(SCALAR_KIND), intent(out) :: jacobianOfInviscidFlux(4,4)
    real(SCALAR_KIND), intent(out), optional :: deltaJacobianOfInviscidFlux(4,4,4)
    real(SCALAR_KIND), intent(in), optional :: deltaVelocity(2,4), deltaTemperature(4)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: contravariantVelocity, phiSquared, deltaContravariantVelocity(4),            &
         deltaPhiSquared(4)

    contravariantVelocity = metrics(1) * velocity(1) +                                       &
         metrics(2) * velocity(2) !... not normalized.

    phiSquared = 0.5_wp * (ratioOfSpecificHeats - 1.0_wp) *                                  &
         (velocity(1) ** 2 + velocity(2) ** 2)

    jacobianOfInviscidFlux(1,1) = 0.0_wp
    jacobianOfInviscidFlux(2,1) = phiSquared * metrics(1) -                                  &
         contravariantVelocity * velocity(1)
    jacobianOfInviscidFlux(3,1) = phiSquared * metrics(2) -                                  &
         contravariantVelocity * velocity(2)
    jacobianOfInviscidFlux(4,1) = contravariantVelocity * ((ratioOfSpecificHeats - 2.0_wp) / &
         (ratioOfSpecificHeats - 1.0_wp) * phiSquared - temperature)

    jacobianOfInviscidFlux(1,2) = metrics(1)
    jacobianOfInviscidFlux(2,2) = contravariantVelocity -                                    &
         (ratioOfSpecificHeats - 2.0_wp) * velocity(1) * metrics(1)
    jacobianOfInviscidFlux(3,2) = velocity(2) * metrics(1) -                                 &
         (ratioOfSpecificHeats - 1.0_wp) * velocity(1) * metrics(2)
    jacobianOfInviscidFlux(4,2) = (temperature + phiSquared /                                &
         (ratioOfSpecificHeats - 1.0_wp)) * metrics(1) - (ratioOfSpecificHeats - 1.0_wp) *   &
         contravariantVelocity * velocity(1)

    jacobianOfInviscidFlux(1,3) = metrics(2)
    jacobianOfInviscidFlux(2,3) = velocity(1) * metrics(2) -                                 &
         (ratioOfSpecificHeats - 1.0_wp) * velocity(2) * metrics(1)
    jacobianOfInviscidFlux(3,3) = contravariantVelocity - (ratioOfSpecificHeats - 2.0_wp) *  &
         velocity(2) * metrics(2)
    jacobianOfInviscidFlux(4,3) = (temperature + phiSquared /                                &
         (ratioOfSpecificHeats - 1.0_wp)) * metrics(2) - (ratioOfSpecificHeats - 1.0_wp) *   &
         contravariantVelocity * velocity(2)

    jacobianOfInviscidFlux(1,4) = 0.0_wp
    jacobianOfInviscidFlux(2,4) = (ratioOfSpecificHeats - 1.0_wp) * metrics(1)
    jacobianOfInviscidFlux(3,4) = (ratioOfSpecificHeats - 1.0_wp) * metrics(2)
    jacobianOfInviscidFlux(4,4) = ratioOfSpecificHeats * contravariantVelocity

    if (present(deltaJacobianOfInviscidFlux)) then

       deltaContravariantVelocity = metrics(1) * deltaVelocity(1,:) + metrics(2) *           &
            deltaVelocity(2,:)

       deltaPhiSquared = (ratioOfSpecificHeats - 1.0_wp) *                                   &
            (velocity(1) * deltaVelocity(1,:) + velocity(2) * deltaVelocity(2,:))

       deltaJacobianOfInviscidFlux(1,1,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(2,1,:) = deltaPhiSquared * metrics(1) -                   &
            deltaContravariantVelocity * velocity(1) -                                       &
            contravariantVelocity * deltaVelocity(1,:)
       deltaJacobianOfInviscidFlux(3,1,:) = deltaPhiSquared * metrics(2) -                   &
            deltaContravariantVelocity * velocity(2) -                                       &
            contravariantVelocity * deltaVelocity(2,:)
       deltaJacobianOfInviscidFlux(4,1,:) = deltaContravariantVelocity *                     &
            ((ratioOfSpecificHeats - 2.0_wp) / (ratioOfSpecificHeats - 1.0_wp) *             &
            phiSquared - temperature) + contravariantVelocity *                              &
            ((ratioOfSpecificHeats - 2.0_wp) / (ratioOfSpecificHeats - 1.0_wp) *             &
            deltaPhiSquared - deltaTemperature)

       deltaJacobianOfInviscidFlux(1,2,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(2,2,:) = deltaContravariantVelocity -                     &
            (ratioOfSpecificHeats - 2.0_wp) * deltaVelocity(1,:) * metrics(1)
       deltaJacobianOfInviscidFlux(3,2,:) = deltaVelocity(2,:) * metrics(1) -                &
            (ratioOfSpecificHeats - 1.0_wp) * deltaVelocity(1,:) * metrics(2)
       deltaJacobianOfInviscidFlux(4,2,:) = (deltaTemperature + deltaPhiSquared /            &
            (ratioOfSpecificHeats - 1.0_wp)) * metrics(1) -                                  &
            (ratioOfSpecificHeats - 1.0_wp) * (deltaContravariantVelocity * velocity(1) +    &
            contravariantVelocity * deltaVelocity(1,:))

       deltaJacobianOfInviscidFlux(1,3,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(2,3,:) = deltaVelocity(1,:) * metrics(2) -                &
            (ratioOfSpecificHeats - 1.0_wp) * deltaVelocity(2,:) * metrics(1)
       deltaJacobianOfInviscidFlux(3,3,:) = deltaContravariantVelocity -                     &
            (ratioOfSpecificHeats - 2.0_wp) * deltaVelocity(2,:) * metrics(2)
       deltaJacobianOfInviscidFlux(4,3,:) = (deltaTemperature + deltaPhiSquared /            &
            (ratioOfSpecificHeats - 1.0_wp)) * metrics(2) -                                  &
            (ratioOfSpecificHeats - 1.0_wp) * (deltaContravariantVelocity * velocity(2) +    &
            contravariantVelocity * deltaVelocity(2,:))

       deltaJacobianOfInviscidFlux(1,4,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(2,4,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(3,4,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(4,4,:) = ratioOfSpecificHeats * deltaContravariantVelocity

    end if

  end subroutine computeInviscidJacobian2D

  pure subroutine computeInviscidJacobian3D(velocity, temperature, metrics,                  &
       ratioOfSpecificHeats, jacobianOfInviscidFlux, deltaJacobianOfInviscidFlux,            &
       deltaVelocity, deltaTemperature)

    ! <<< Arguments >>>
    real(SCALAR_KIND), intent(in) :: velocity(3), temperature, metrics(3),                   &
         ratioOfSpecificHeats
    real(SCALAR_KIND), intent(out) :: jacobianOfInviscidFlux(5,5)
    real(SCALAR_KIND), intent(out), optional :: deltaJacobianOfInviscidFlux(5,5,5)
    real(SCALAR_KIND), intent(in), optional :: deltaVelocity(3,5), deltaTemperature(5)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: contravariantVelocity, phiSquared, deltaContravariantVelocity(5),            &
         deltaPhiSquared(5)

    contravariantVelocity = metrics(1) * velocity(1) + metrics(2) * velocity(2) +            &
         metrics(3) * velocity(3) !... not normalized.

    phiSquared = 0.5_wp * (ratioOfSpecificHeats - 1.0_wp) *                                  &
         (velocity(1) ** 2 + velocity(2) ** 2 + velocity(3) ** 2)

    jacobianOfInviscidFlux(1,1) = 0.0_wp
    jacobianOfInviscidFlux(2,1) = phiSquared * metrics(1) -                                  &
         contravariantVelocity * velocity(1)
    jacobianOfInviscidFlux(3,1) = phiSquared * metrics(2) -                                  &
         contravariantVelocity * velocity(2)
    jacobianOfInviscidFlux(4,1) = phiSquared * metrics(3) -                                  &
         contravariantVelocity * velocity(3)
    jacobianOfInviscidFlux(5,1) = contravariantVelocity * ((ratioOfSpecificHeats - 2.0_wp) / &
         (ratioOfSpecificHeats - 1.0_wp) * phiSquared - temperature)

    jacobianOfInviscidFlux(1,2) = metrics(1)
    jacobianOfInviscidFlux(2,2) = contravariantVelocity -                                    &
         (ratioOfSpecificHeats - 2.0_wp) * velocity(1) * metrics(1)
    jacobianOfInviscidFlux(3,2) = velocity(2) * metrics(1) -                                 &
         (ratioOfSpecificHeats - 1.0_wp) * velocity(1) * metrics(2)
    jacobianOfInviscidFlux(4,2) = velocity(3) * metrics(1) -                                 &
         (ratioOfSpecificHeats - 1.0_wp) * velocity(1) * metrics(3)
    jacobianOfInviscidFlux(5,2) = (temperature + phiSquared /                                &
         (ratioOfSpecificHeats - 1.0_wp)) * metrics(1) - (ratioOfSpecificHeats - 1.0_wp) *   &
         contravariantVelocity * velocity(1)

    jacobianOfInviscidFlux(1,3) = metrics(2)
    jacobianOfInviscidFlux(2,3) = velocity(1) * metrics(2) -                                 &
         (ratioOfSpecificHeats - 1.0_wp) * velocity(2) * metrics(1)
    jacobianOfInviscidFlux(3,3) = contravariantVelocity - (ratioOfSpecificHeats - 2.0_wp) *  &
         velocity(2) * metrics(2)
    jacobianOfInviscidFlux(4,3) = velocity(3) * metrics(2) -                                 &
         (ratioOfSpecificHeats - 1.0_wp) * velocity(2) * metrics(3)
    jacobianOfInviscidFlux(5,3) = (temperature + phiSquared /                                &
         (ratioOfSpecificHeats - 1.0_wp)) * metrics(2) - (ratioOfSpecificHeats - 1.0_wp) *   &
         contravariantVelocity * velocity(2)

    jacobianOfInviscidFlux(1,4) = metrics(3)
    jacobianOfInviscidFlux(2,4) = velocity(1) * metrics(3) -                                 &
         (ratioOfSpecificHeats - 1.0_wp) * velocity(3) * metrics(1)
    jacobianOfInviscidFlux(3,4) = velocity(2) * metrics(3) -                                 &
         (ratioOfSpecificHeats - 1.0_wp) * velocity(3) * metrics(2)
    jacobianOfInviscidFlux(4,4) = contravariantVelocity - (ratioOfSpecificHeats - 2.0_wp) *  &
         velocity(3) * metrics(3)
    jacobianOfInviscidFlux(5,4) = (temperature + phiSquared /                                &
         (ratioOfSpecificHeats - 1.0_wp)) * metrics(3) - (ratioOfSpecificHeats - 1.0_wp) *   &
         contravariantVelocity * velocity(3)

    jacobianOfInviscidFlux(1,5) = 0.0_wp
    jacobianOfInviscidFlux(2,5) = (ratioOfSpecificHeats - 1.0_wp) * metrics(1)
    jacobianOfInviscidFlux(3,5) = (ratioOfSpecificHeats - 1.0_wp) * metrics(2)
    jacobianOfInviscidFlux(4,5) = (ratioOfSpecificHeats - 1.0_wp) * metrics(3)
    jacobianOfInviscidFlux(5,5) = ratioOfSpecificHeats * contravariantVelocity

    if (present(deltaJacobianOfInviscidFlux)) then

       deltaContravariantVelocity =                                                          &
            metrics(1) * deltaVelocity(1,:) +                                                &
            metrics(2) * deltaVelocity(2,:) +                                                &
            metrics(3) * deltaVelocity(3,:)

       deltaPhiSquared = (ratioOfSpecificHeats - 1.0_wp) *                                   &
            (velocity(1) * deltaVelocity(1,:) + velocity(2) * deltaVelocity(2,:) +           &
            velocity(3) * deltaVelocity(3,:))

       deltaJacobianOfInviscidFlux(1,1,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(2,1,:) = deltaPhiSquared * metrics(1) -                   &
            deltaContravariantVelocity * velocity(1) - contravariantVelocity *               &
            deltaVelocity(1,:)
       deltaJacobianOfInviscidFlux(3,1,:) = deltaPhiSquared * metrics(2) -                   &
            deltaContravariantVelocity * velocity(2) - contravariantVelocity *               &
            deltaVelocity(2,:)
       deltaJacobianOfInviscidFlux(4,1,:) = deltaPhiSquared * metrics(3) -                   &
            deltaContravariantVelocity * velocity(3) - contravariantVelocity *               &
            deltaVelocity(3,:)
       deltaJacobianOfInviscidFlux(5,1,:) = deltaContravariantVelocity *                     &
            ((ratioOfSpecificHeats - 2.0_wp) / (ratioOfSpecificHeats - 1.0_wp) *             &
            phiSquared - temperature) + contravariantVelocity *                              &
            ((ratioOfSpecificHeats - 2.0_wp) / (ratioOfSpecificHeats - 1.0_wp) *             &
            deltaPhiSquared - deltaTemperature)

       deltaJacobianOfInviscidFlux(1,2,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(2,2,:) = deltaContravariantVelocity -                     &
            (ratioOfSpecificHeats - 2.0_wp) * deltaVelocity(1,:) * metrics(1)
       deltaJacobianOfInviscidFlux(3,2,:) = deltaVelocity(2,:) * metrics(1) -                &
            (ratioOfSpecificHeats - 1.0_wp) * deltaVelocity(1,:) * metrics(2)
       deltaJacobianOfInviscidFlux(4,2,:) = deltaVelocity(3,:) * metrics(1) -                &
            (ratioOfSpecificHeats - 1.0_wp) * deltaVelocity(1,:) * metrics(3)
       deltaJacobianOfInviscidFlux(5,2,:) = (deltaTemperature + deltaPhiSquared /            &
            (ratioOfSpecificHeats - 1.0_wp)) * metrics(1) -                                  &
            (ratioOfSpecificHeats - 1.0_wp) * (deltaContravariantVelocity * velocity(1) +    &
            contravariantVelocity * deltaVelocity(1,:))

       deltaJacobianOfInviscidFlux(1,3,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(2,3,:) = deltaVelocity(1,:) * metrics(2) -                &
            (ratioOfSpecificHeats - 1.0_wp) * deltaVelocity(2,:) * metrics(1)
       deltaJacobianOfInviscidFlux(3,3,:) = deltaContravariantVelocity -                     &
            (ratioOfSpecificHeats - 2.0_wp) * deltaVelocity(2,:) * metrics(2)
       deltaJacobianOfInviscidFlux(4,3,:) = deltaVelocity(3,:) * metrics(2) -                &
            (ratioOfSpecificHeats - 1.0_wp) * deltaVelocity(2,:) * metrics(3)
       deltaJacobianOfInviscidFlux(5,3,:) = (deltaTemperature + deltaPhiSquared /            &
            (ratioOfSpecificHeats - 1.0_wp)) * metrics(2) -                                  &
            (ratioOfSpecificHeats - 1.0_wp) * (deltaContravariantVelocity * velocity(2) +    &
            contravariantVelocity * deltaVelocity(2,:))

       deltaJacobianOfInviscidFlux(1,4,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(2,4,:) = deltaVelocity(1,:) * metrics(3) -                &
            (ratioOfSpecificHeats - 1.0_wp) * deltaVelocity(3,:) * metrics(1)
       deltaJacobianOfInviscidFlux(3,4,:) = deltaVelocity(2,:) * metrics(3) -                &
            (ratioOfSpecificHeats - 1.0_wp) * deltaVelocity(3,:) * metrics(2)
       deltaJacobianOfInviscidFlux(4,4,:) = deltaContravariantVelocity -                     &
            (ratioOfSpecificHeats - 2.0_wp) * deltaVelocity(3,:) * metrics(3)
       deltaJacobianOfInviscidFlux(5,4,:) = (deltaTemperature + deltaPhiSquared /            &
            (ratioOfSpecificHeats - 1.0_wp)) * metrics(3) -                                  &
            (ratioOfSpecificHeats - 1.0_wp) * (deltaContravariantVelocity * velocity(3) +    &
            contravariantVelocity * deltaVelocity(3,:))

       deltaJacobianOfInviscidFlux(1,5,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(2,5,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(3,5,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(4,5,:) = 0.0_wp
       deltaJacobianOfInviscidFlux(5,5,:) = ratioOfSpecificHeats * deltaContravariantVelocity

    end if

  end subroutine computeInviscidJacobian3D

  pure subroutine computeIncomingInviscidJacobianRoe(nDimensions, ijkIndexL, ijkIndexR,      &
       conservedVariablesL, conservedVariablesR, metricsAlongNormalDirection,                &
       ratioOfSpecificHeats, incomingDirection, incomingJacobianOfInviscidFlux,              &
       deltaIncomingJacobianOfInviscidFlux, deltaConservedVariablesL)

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: nDimensions, ijkIndexL, ijkIndexR
    real(SCALAR_KIND), intent(in) :: conservedVariablesL(:,:), conservedVariablesR(:,:),     &
         metricsAlongNormalDirection(:,:)
    real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
    integer, intent(in) :: incomingDirection
    real(SCALAR_KIND), intent(out) :: incomingJacobianOfInviscidFlux(:,:)
    real(SCALAR_KIND), intent(out), optional :: deltaIncomingJacobianOfInviscidFlux(:,:,:)
    real(SCALAR_KIND), intent(in), optional :: deltaConservedVariablesL(:,:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, nUnknowns
    real(wp), allocatable :: localConservedVariablesL(:), localConservedVariablesR(:),       &
         localRoeAverage(:), localVelocity(:), localMetricsAlongNormalDirection(:),          &
         deltaConservedVariablesL_(:,:), deltaRoeAverage(:,:), deltaSpecificVolume(:),       &
         deltaVelocity(:,:), deltaTemperature(:)
    real(wp) :: localSpecificVolume, localTemperature

    nUnknowns = size(conservedVariablesL, 2)

    allocate(localConservedVariablesL(nUnknowns))
    allocate(localConservedVariablesR(nUnknowns))
    allocate(localRoeAverage(nUnknowns))
    allocate(localVelocity(nDimensions))
    allocate(localMetricsAlongNormalDirection(nDimensions))

    localConservedVariablesL = conservedVariablesL(ijkIndexL,:)
    localConservedVariablesR = conservedVariablesR(ijkIndexR,:)
    localMetricsAlongNormalDirection = metricsAlongNormalDirection(ijkIndexL,:)

    if (present(deltaIncomingJacobianOfInviscidFlux)) then

       allocate(deltaRoeAverage(nUnknowns, nUnknowns))

       if (present(deltaConservedVariablesL)) then
          allocate(deltaConservedVariablesL_(nUnknowns, nUnknowns),                          &
               source = deltaConservedVariablesL)
       else
          allocate(deltaConservedVariablesL_(nUnknowns, nUnknowns), source = 0.0_wp)
          do i = 1, nUnknowns
             deltaConservedVariablesL_(i,i) = 1.0_wp
          end do
       end if

       call computeLocalRoeAverage(nDimensions, localConservedVariablesL,                    &
            localConservedVariablesR, ratioOfSpecificHeats, localRoeAverage,                 &
            deltaRoeAverage, deltaConservedVariablesL_)

    else

       call computeLocalRoeAverage(nDimensions, localConservedVariablesL,                    &
            localConservedVariablesR, ratioOfSpecificHeats, localRoeAverage)

    end if

    localSpecificVolume = 1.0_wp / localRoeAverage(1)
    localVelocity = localSpecificVolume * localRoeAverage(2:nDimensions+1)
    localTemperature = ratioOfSpecificHeats * (localSpecificVolume *                         &
            localRoeAverage(nDimensions+2) - 0.5_wp * sum(localVelocity ** 2))

    if (present(deltaIncomingJacobianOfInviscidFlux)) then

       allocate(deltaSpecificVolume(nUnknowns))
       allocate(deltaVelocity(nDimensions, nUnknowns))
       allocate(deltaTemperature(nUnknowns))

       deltaSpecificVolume = -localSpecificVolume ** 2 * deltaRoeAverage(1,:)
       do i = 1, nDimensions
          deltaVelocity(i,:) = deltaSpecificVolume * localRoeAverage(i+1) +                  &
               localSpecificVolume * deltaRoeAverage(i+1,:)
       end do

       deltaTemperature = ratioOfSpecificHeats * (deltaSpecificVolume *                      &
            localRoeAverage(nDimensions+2) + localSpecificVolume *                           &
            deltaRoeAverage(nDimensions+2,:) - matmul(localVelocity, deltaVelocity))

       select case (nDimensions)
       case (1)
          call computeIncomingInviscidJacobian1D(localVelocity, localTemperature,            &
               localMetricsAlongNormalDirection, ratioOfSpecificHeats, incomingDirection,    &
               incomingJacobianOfInviscidFlux, deltaIncomingJacobianOfInviscidFlux,          &
               deltaVelocity, deltaTemperature)
       case (2)
          call computeIncomingInviscidJacobian2D(localRoeAverage, localSpecificVolume,       &
               localVelocity, localTemperature, localMetricsAlongNormalDirection,            &
               ratioOfSpecificHeats, incomingDirection, incomingJacobianOfInviscidFlux,      &
               deltaIncomingJacobianOfInviscidFlux, deltaRoeAverage, deltaSpecificVolume,    &
               deltaVelocity, deltaTemperature)
       case (3)
          call computeIncomingInviscidJacobian3D(localRoeAverage, localSpecificVolume,       &
               localVelocity, localTemperature, localMetricsAlongNormalDirection,            &
               ratioOfSpecificHeats, incomingDirection, incomingJacobianOfInviscidFlux,      &
               deltaIncomingJacobianOfInviscidFlux, deltaRoeAverage, deltaSpecificVolume,    &
               deltaVelocity, deltaTemperature)
       end select

    else

       select case (nDimensions)
       case (1)
          call computeIncomingInviscidJacobian1D(localVelocity, localTemperature,            &
               localMetricsAlongNormalDirection, ratioOfSpecificHeats, incomingDirection,    &
               incomingJacobianOfInviscidFlux)
       case (2)
          call computeIncomingInviscidJacobian2D(localRoeAverage, localSpecificVolume,       &
               localVelocity, localTemperature, localMetricsAlongNormalDirection,            &
               ratioOfSpecificHeats, incomingDirection, incomingJacobianOfInviscidFlux)
       case (3)
          call computeIncomingInviscidJacobian3D(localRoeAverage, localSpecificVolume,       &
               localVelocity, localTemperature, localMetricsAlongNormalDirection,            &
               ratioOfSpecificHeats, incomingDirection, incomingJacobianOfInviscidFlux)
       end select

    end if

    SAFE_DEALLOCATE(deltaTemperature)
    SAFE_DEALLOCATE(deltaVelocity)
    SAFE_DEALLOCATE(deltaSpecificVolume)
    SAFE_DEALLOCATE(deltaConservedVariablesL_)
    SAFE_DEALLOCATE(deltaRoeAverage)
    SAFE_DEALLOCATE(localMetricsAlongNormalDirection)
    SAFE_DEALLOCATE(localVelocity)
    SAFE_DEALLOCATE(localRoeAverage)
    SAFE_DEALLOCATE(localConservedVariablesR)
    SAFE_DEALLOCATE(localConservedVariablesL)

  end subroutine computeIncomingInviscidJacobianRoe

  pure subroutine computeIncomingInviscidJacobianEqual(nDimensions, ijkIndex,                &
       conservedVariables, metricsAlongNormalDirection, ratioOfSpecificHeats,                &
       incomingDirection, incomingJacobianOfInviscidFlux,                                    &
       deltaIncomingJacobianOfInviscidFlux, deltaConservedVariables, specificVolume,         &
       velocity, temperature)

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: nDimensions, ijkIndex
    real(SCALAR_KIND), intent(in) :: conservedVariables(:,:), metricsAlongNormalDirection(:,:)
    real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
    integer, intent(in) :: incomingDirection
    real(SCALAR_KIND), intent(out) :: incomingJacobianOfInviscidFlux(:,:)
    real(SCALAR_KIND), intent(out), optional :: deltaIncomingJacobianOfInviscidFlux(:,:,:)
    real(SCALAR_KIND), intent(in), optional :: deltaConservedVariables(:,:),                 &
         specificVolume(:), velocity(:,:), temperature(:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, nUnknowns
    real(wp), allocatable :: localConservedVariables(:), localVelocity(:),                   &
         localMetricsAlongNormalDirection(:), deltaConservedVariables_(:,:),                 &
         deltaSpecificVolume(:), deltaVelocity(:,:), deltaTemperature(:)
    real(wp) :: localSpecificVolume, localTemperature

    nUnknowns = size(conservedVariables, 2)

    allocate(localConservedVariables(nUnknowns))
    allocate(localVelocity(nDimensions))
    allocate(localMetricsAlongNormalDirection(nDimensions))

    localConservedVariables = conservedVariables(ijkIndex,:)
    localMetricsAlongNormalDirection = metricsAlongNormalDirection(ijkIndex,:)

    if (present(specificVolume)) then
       localSpecificVolume = specificVolume(ijkIndex)
    else
       localSpecificVolume = 1.0_wp / localConservedVariables(1)
    end if

    if (present(velocity)) then
       localVelocity = velocity(ijkIndex,:)
    else
       localVelocity = localSpecificVolume * localConservedVariables(2:nDimensions+1)
    end if

    if (present(temperature)) then
       localTemperature = temperature(ijkIndex)
    else
       localTemperature = ratioOfSpecificHeats * (localSpecificVolume *                      &
            localConservedVariables(nDimensions+2) - 0.5_wp * sum(localVelocity ** 2))
    end if

    if (present(deltaIncomingJacobianOfInviscidFlux)) then

       if (present(deltaConservedVariables)) then
          allocate(deltaConservedVariables_(nUnknowns, nUnknowns),                           &
               source = deltaConservedVariables)
       else
          allocate(deltaConservedVariables_(nUnknowns, nUnknowns), source = 0.0_wp)
          do i = 1, nUnknowns
             deltaConservedVariables_(i,i) = 1.0_wp
          end do
       end if

       allocate(deltaSpecificVolume(nUnknowns))
       allocate(deltaVelocity(nDimensions, nUnknowns))
       allocate(deltaTemperature(nUnknowns))

       deltaSpecificVolume = -localSpecificVolume ** 2 * deltaConservedVariables_(1,:)
       do i = 1, nDimensions
          deltaVelocity(i,:) = deltaSpecificVolume * localConservedVariables(i+1) +          &
               localSpecificVolume * deltaConservedVariables_(i+1,:)
       end do

       deltaTemperature = ratioOfSpecificHeats * (deltaSpecificVolume *                      &
            localConservedVariables(nDimensions+2) + localSpecificVolume *                   &
            deltaConservedVariables_(nDimensions+2,:) - matmul(localVelocity, deltaVelocity))

       select case (nDimensions)
       case (1)
          call computeIncomingInviscidJacobian1D(localVelocity, localTemperature,            &
               localMetricsAlongNormalDirection, ratioOfSpecificHeats, incomingDirection,    &
               incomingJacobianOfInviscidFlux, deltaIncomingJacobianOfInviscidFlux,          &
               deltaVelocity, deltaTemperature)
       case (2)
          call computeIncomingInviscidJacobian2D(localConservedVariables,                    &
               localSpecificVolume, localVelocity, localTemperature,                         &
               localMetricsAlongNormalDirection, ratioOfSpecificHeats, incomingDirection,    &
               incomingJacobianOfInviscidFlux, deltaIncomingJacobianOfInviscidFlux,          &
               deltaConservedVariables, deltaSpecificVolume, deltaVelocity, deltaTemperature)
       case (3)
          call computeIncomingInviscidJacobian3D(localConservedVariables,                    &
               localSpecificVolume, localVelocity, localTemperature,                         &
               localMetricsAlongNormalDirection, ratioOfSpecificHeats, incomingDirection,    &
               incomingJacobianOfInviscidFlux, deltaIncomingJacobianOfInviscidFlux,          &
               deltaConservedVariables, deltaSpecificVolume, deltaVelocity, deltaTemperature)
       end select

    else

       select case (nDimensions)
       case (1)
          call computeIncomingInviscidJacobian1D(localVelocity, localTemperature,            &
               localMetricsAlongNormalDirection, ratioOfSpecificHeats, incomingDirection,    &
               incomingJacobianOfInviscidFlux)
       case (2)
          call computeIncomingInviscidJacobian2D(localConservedVariables,                    &
               localSpecificVolume, localVelocity, localTemperature,                         &
               localMetricsAlongNormalDirection, ratioOfSpecificHeats, incomingDirection,    &
               incomingJacobianOfInviscidFlux)
       case (3)
          call computeIncomingInviscidJacobian3D(localConservedVariables,                    &
               localSpecificVolume, localVelocity, localTemperature,                         &
               localMetricsAlongNormalDirection, ratioOfSpecificHeats, incomingDirection,    &
               incomingJacobianOfInviscidFlux)
       end select

    end if

    SAFE_DEALLOCATE(deltaTemperature)
    SAFE_DEALLOCATE(deltaVelocity)
    SAFE_DEALLOCATE(deltaSpecificVolume)
    SAFE_DEALLOCATE(deltaConservedVariables_)
    SAFE_DEALLOCATE(localMetricsAlongNormalDirection)
    SAFE_DEALLOCATE(localVelocity)
    SAFE_DEALLOCATE(localConservedVariables)

  end subroutine computeIncomingInviscidJacobianEqual

  pure subroutine computeIncomingInviscidJacobian1D(velocity, temperature, metrics,          &
       ratioOfSpecificHeats, incomingDirection, incomingJacobianOfInviscidFlux,              &
       deltaIncomingJacobianOfInviscidFlux, deltaVelocity, deltaTemperature)

    implicit none

    ! <<< Arguments >>>
    real(SCALAR_KIND), intent(in) :: velocity(1), temperature, metrics(1)
    real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
    integer, intent(in) :: incomingDirection
    real(SCALAR_KIND), intent(out) :: incomingJacobianOfInviscidFlux(3,3)
    real(SCALAR_KIND), intent(out), optional :: deltaIncomingJacobianOfInviscidFlux(3,3,3)
    real(SCALAR_KIND), intent(in), optional :: deltaVelocity(1,3), deltaTemperature(3)

    ! <<< Local variables >>>
    integer :: i, j
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: arcLength, normalizedMetrics(1), contravariantVelocity, speedOfSound,        &
         phiSquared, rightEigenvectors(3,3), eigenvalues(3), leftEigenvectors(3,3),          &
         deltaContravariantVelocity(3), deltaSpeedOfSound(3), deltaPhiSquared(3),            &
         deltaRightEigenvectors(3,3,3), deltaEigenvalues(3,3), deltaLeftEigenvectors(3,3,3), &
         temp(3)

    ! Normalize the metrics.
    arcLength = abs(metrics(1))
    normalizedMetrics = metrics / arcLength

    ! Other dependent variables.
    contravariantVelocity = normalizedMetrics(1) * velocity(1)
    speedOfSound = sqrt((ratioOfSpecificHeats - 1.0_wp) * temperature)
    phiSquared = 0.5_wp * (ratioOfSpecificHeats - 1.0_wp) * velocity(1) ** 2

    ! Eigenvalues.
    eigenvalues(1) = contravariantVelocity
    eigenvalues(2) = contravariantVelocity + speedOfSound
    eigenvalues(3) = contravariantVelocity - speedOfSound
    eigenvalues = arcLength * eigenvalues

    if (present(deltaIncomingJacobianOfInviscidFlux)) then

       ! Compute variations of other dependent variables.
       deltaContravariantVelocity = normalizedMetrics(1) * deltaVelocity(1,:)
       deltaSpeedOfSound = 0.5_wp / speedOfSound *                                           &
            (ratioOfSpecificHeats - 1.0_wp) * deltaTemperature
       deltaPhiSquared = (ratioOfSpecificHeats - 1.0_wp) *                                   &
            (velocity(1) * deltaVelocity(1,:))

       ! Variation of matrix containing eigenvalues.
       deltaEigenvalues(1,:) = deltaContravariantVelocity
       deltaEigenvalues(2,:) = deltaContravariantVelocity + deltaSpeedOfSound
       deltaEigenvalues(3,:) = deltaContravariantVelocity - deltaSpeedOfSound
       deltaEigenvalues = arcLength * deltaEigenvalues

    end if

    ! Zero-out the eigenvalues corresponding to outgoing characteristics and corresponding
    ! variations.
    do i = 1, 3
       if (incomingDirection * eigenvalues(i) < 0.0_wp) then
          eigenvalues(i) = 0.0_wp
          if (present(deltaIncomingJacobianOfInviscidFlux)) deltaEigenvalues(i,:) = 0.0_wp
       end if
    end do

    ! Matrix whose columns are the right eigenvectors:

    rightEigenvectors(1,1) = 1.0_wp
    rightEigenvectors(2,1) = velocity(1)
    rightEigenvectors(3,1) = phiSquared / (ratioOfSpecificHeats - 1.0_wp)

    rightEigenvectors(1,2) = 1.0_wp
    rightEigenvectors(2,2) = velocity(1) + normalizedMetrics(1) * speedOfSound
    rightEigenvectors(3,2) = temperature + phiSquared / (ratioOfSpecificHeats - 1.0_wp) +    &
         speedOfSound * contravariantVelocity

    rightEigenvectors(1,3) = 1.0_wp
    rightEigenvectors(2,3) = velocity(1) - normalizedMetrics(1) * speedOfSound
    rightEigenvectors(3,3) = temperature + phiSquared / (ratioOfSpecificHeats - 1.0_wp) -    &
         speedOfSound * contravariantVelocity

    ! Matrix whose rows are the left eigenvectors:

    leftEigenvectors(1,1) = 1.0_wp - phiSquared / speedOfSound ** 2
    leftEigenvectors(2,1) = 0.5_wp * (phiSquared / speedOfSound ** 2 -                       &
         contravariantVelocity / speedOfSound)
    leftEigenvectors(3,1) = 0.5_wp * (phiSquared / speedOfSound ** 2 +                       &
         contravariantVelocity / speedOfSound)

    leftEigenvectors(1,2) = velocity(1) / temperature
    leftEigenvectors(2,2) = - 0.5_wp * (velocity(1) / temperature -                          &
         normalizedMetrics(1) / speedOfSound)
    leftEigenvectors(3,2) = - 0.5_wp * (velocity(1) / temperature +                          &
         normalizedMetrics(1) / speedOfSound)

    leftEigenvectors(1,3) = - 1.0_wp / temperature
    leftEigenvectors(2,3) = 0.5_wp / temperature
    leftEigenvectors(3,3) = 0.5_wp / temperature

    ! ``Incoming'' part.
    do j = 1, 3
       do i = 1, 3
          incomingJacobianOfInviscidFlux(i,j) =                                              &
               rightEigenvectors(i,1) * eigenvalues(1) * leftEigenvectors(1,j) +             &
               rightEigenvectors(i,2) * eigenvalues(2) * leftEigenvectors(2,j) +             &
               rightEigenvectors(i,3) * eigenvalues(3) * leftEigenvectors(3,j)
       end do
    end do

    if (present(deltaIncomingJacobianOfInviscidFlux)) then

       ! Variation of the matrix whose columns are the right eigenvectors:

       deltaRightEigenvectors(1,1,:) = 0.0_wp
       deltaRightEigenvectors(2,1,:) = deltaVelocity(1,:)
       deltaRightEigenvectors(3,1,:) = deltaPhiSquared /                                     &
            (ratioOfSpecificHeats - 1.0_wp)

       deltaRightEigenvectors(1,2,:) = 0.0_wp
       deltaRightEigenvectors(2,2,:) = deltaVelocity(1,:) +                                  &
            normalizedMetrics(1) * deltaSpeedOfSound
       deltaRightEigenvectors(3,2,:) = deltaTemperature +                                    &
            deltaPhiSquared / (ratioOfSpecificHeats - 1.0_wp) +                              &
            deltaSpeedOfSound * contravariantVelocity +                                      &
            speedOfSound * deltaContravariantVelocity

       deltaRightEigenvectors(1,3,:) = 0.0_wp
       deltaRightEigenvectors(2,3,:) = deltaVelocity(1,:) -                                  &
            normalizedMetrics(1) * deltaSpeedOfSound
       deltaRightEigenvectors(3,3,:) = deltaTemperature +                                    &
            deltaPhiSquared / (ratioOfSpecificHeats - 1.0_wp) -                              &
            deltaSpeedOfSound * contravariantVelocity -                                      &
            speedOfSound * deltaContravariantVelocity

       ! Variation of the matrix whose rows are the left eigenvectors:

       temp = deltaPhiSquared / speedOfSound ** 2 -                                          &
            2.0_wp * phiSquared / speedOfSound ** 3 * deltaSpeedOfSound
       deltaLeftEigenvectors(1,1,:) = -temp
       deltaLeftEigenvectors(2,1,:) = 0.5_wp * (temp -                                       &
            deltaContravariantVelocity / speedOfSound +                                      &
            contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)
       deltaLeftEigenvectors(3,1,:) = 0.5_wp * (temp +                                       &
            deltaContravariantVelocity / speedOfSound -                                      &
            contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)

       temp = deltaVelocity(1,:) / temperature -                                             &
            velocity(1) / temperature ** 2 * deltaTemperature
       deltaLeftEigenvectors(1,2,:) = temp
       deltaLeftEigenvectors(2,2,:) = -0.5_wp * (temp +                                      &
            normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)
       deltaLeftEigenvectors(3,2,:) = -0.5_wp * (temp -                                      &
            normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)

       temp =  -1.0_wp / temperature ** 2 * deltaTemperature
       deltaLeftEigenvectors(1,3,:) = -temp
       deltaLeftEigenvectors(2,3,:) = 0.5_wp * temp
       deltaLeftEigenvectors(3,3,:) = 0.5_wp * temp

       ! Variation of the ``incoming'' part.
       do j = 1, 3
          do i = 1, 3
             deltaIncomingJacobianOfInviscidFlux(i,j,:) =                                    &
                  deltaRightEigenvectors(i,1,:) * eigenvalues(1) * leftEigenvectors(1,j) +   &
                  deltaRightEigenvectors(i,2,:) * eigenvalues(2) * leftEigenvectors(2,j) +   &
                  deltaRightEigenvectors(i,3,:) * eigenvalues(3) * leftEigenvectors(3,j) +   &
                  rightEigenvectors(i,1) * deltaEigenvalues(1,:) * leftEigenvectors(1,j) +   &
                  rightEigenvectors(i,2) * deltaEigenvalues(2,:) * leftEigenvectors(2,j) +   &
                  rightEigenvectors(i,3) * deltaEigenvalues(3,:) * leftEigenvectors(3,j) +   &
                  rightEigenvectors(i,1) * eigenvalues(1) * deltaLeftEigenvectors(1,j,:) +   &
                  rightEigenvectors(i,2) * eigenvalues(2) * deltaLeftEigenvectors(2,j,:) +   &
                  rightEigenvectors(i,3) * eigenvalues(3) * deltaLeftEigenvectors(3,j,:)
          end do
       end do

    end if

  end subroutine computeIncomingInviscidJacobian1D

  pure subroutine computeIncomingInviscidJacobian2D(conservedVariables, specificVolume,      &
       velocity, temperature, metrics, ratioOfSpecificHeats, incomingDirection,              &
       incomingJacobianOfInviscidFlux, deltaIncomingJacobianOfInviscidFlux,                  &
       deltaConservedVariables, deltaSpecificVolume, deltaVelocity, deltaTemperature)

    implicit none

    ! <<< Arguments >>>
    real(SCALAR_KIND), intent(in) :: conservedVariables(4), specificVolume, velocity(2),     &
         temperature, metrics(2)
    real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
    integer, intent(in) :: incomingDirection
    real(SCALAR_KIND), intent(out) :: incomingJacobianOfInviscidFlux(4,4)
    real(SCALAR_KIND), intent(out), optional :: deltaIncomingJacobianOfInviscidFlux(4,4,4)
    real(SCALAR_KIND), intent(in), optional :: deltaConservedVariables(4,4),                 &
         deltaSpecificVolume(4), deltaVelocity(2,4), deltaTemperature(4)

    ! <<< Local variables >>>
    integer :: i, j
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: arcLength, normalizedMetrics(2), contravariantVelocity, speedOfSound,        &
         phiSquared, rightEigenvectors(4,4), eigenvalues(4), leftEigenvectors(4,4),          &
         deltaContravariantVelocity(4), deltaSpeedOfSound(4), deltaPhiSquared(4),            &
         deltaRightEigenvectors(4,4,4), deltaEigenvalues(4,4), deltaLeftEigenvectors(4,4,4), &
         temp(4)

    ! Normalize the metrics.
    arcLength = sqrt(metrics(1) ** 2 + metrics(2) ** 2)
    normalizedMetrics = metrics / arcLength

    ! Other dependent variables.
    contravariantVelocity = normalizedMetrics(1) * velocity(1) +                             &
         normalizedMetrics(2) * velocity(2)
    speedOfSound = sqrt((ratioOfSpecificHeats - 1.0_wp) * temperature)
    phiSquared = 0.5_wp * (ratioOfSpecificHeats - 1.0_wp) *                                  &
         (velocity(1) ** 2 + velocity(2) ** 2)

    ! Eigenvalues.
    eigenvalues(1) = contravariantVelocity
    eigenvalues(2) = contravariantVelocity
    eigenvalues(3) = contravariantVelocity + speedOfSound
    eigenvalues(4) = contravariantVelocity - speedOfSound
    eigenvalues = arcLength * eigenvalues

    if (present(deltaIncomingJacobianOfInviscidFlux)) then

       ! Compute variations of other dependent variables.
       deltaContravariantVelocity = normalizedMetrics(1) * deltaVelocity(1,:) +              &
            normalizedMetrics(2) * deltaVelocity(2,:)
       deltaSpeedOfSound = 0.5_wp / speedOfSound *                                           &
            (ratioOfSpecificHeats - 1.0_wp) * deltaTemperature
       deltaPhiSquared = (ratioOfSpecificHeats - 1.0_wp) *                                   &
            (velocity(1) * deltaVelocity(1,:) + velocity(2) * deltaVelocity(2,:))

       ! Variation of matrix containing eigenvalues.
       deltaEigenvalues(1,:) = deltaContravariantVelocity
       deltaEigenvalues(2,:) = deltaContravariantVelocity
       deltaEigenvalues(3,:) = deltaContravariantVelocity + deltaSpeedOfSound
       deltaEigenvalues(4,:) = deltaContravariantVelocity - deltaSpeedOfSound
       deltaEigenvalues = arcLength * deltaEigenvalues

    end if

    ! Zero-out the eigenvalues corresponding to outgoing characteristics and corresponding
    ! variations.
    do i = 1, 4
       if (incomingDirection * eigenvalues(i) < 0.0_wp) then
          eigenvalues(i) = 0.0_wp
          if (present(deltaIncomingJacobianOfInviscidFlux)) deltaEigenvalues(i,:) = 0.0_wp
       end if
    end do

    ! Matrix whose columns are the right eigenvectors:

    rightEigenvectors(1,1) = 1.0_wp
    rightEigenvectors(2,1) = velocity(1)
    rightEigenvectors(3,1) = velocity(2)
    rightEigenvectors(4,1) = phiSquared / (ratioOfSpecificHeats - 1.0_wp)

    rightEigenvectors(1,2) = 0.0_wp
    rightEigenvectors(2,2) = normalizedMetrics(2) * conservedVariables(1)
    rightEigenvectors(3,2) = - normalizedMetrics(1) * conservedVariables(1)
    rightEigenvectors(4,2) = conservedVariables(1) * (normalizedMetrics(2) * velocity(1) -   &
         normalizedMetrics(1) * velocity(2))

    rightEigenvectors(1,3) = 1.0_wp
    rightEigenvectors(2,3) = velocity(1) + normalizedMetrics(1) * speedOfSound
    rightEigenvectors(3,3) = velocity(2) + normalizedMetrics(2) * speedOfSound
    rightEigenvectors(4,3) = temperature + phiSquared / (ratioOfSpecificHeats - 1.0_wp) +    &
         speedOfSound * contravariantVelocity

    rightEigenvectors(1,4) = 1.0_wp
    rightEigenvectors(2,4) = velocity(1) - normalizedMetrics(1) * speedOfSound
    rightEigenvectors(3,4) = velocity(2) - normalizedMetrics(2) * speedOfSound
    rightEigenvectors(4,4) = temperature + phiSquared / (ratioOfSpecificHeats - 1.0_wp) -    &
         speedOfSound * contravariantVelocity

    ! Matrix whose rows are the left eigenvectors:

    leftEigenvectors(1,1) = 1.0_wp - phiSquared / speedOfSound ** 2
    leftEigenvectors(2,1) = - specificVolume * (normalizedMetrics(2) * velocity(1) -         &
         normalizedMetrics(1) * velocity(2))
    leftEigenvectors(3,1) = 0.5_wp * (phiSquared / speedOfSound ** 2 -                       &
         contravariantVelocity / speedOfSound)
    leftEigenvectors(4,1) = 0.5_wp * (phiSquared / speedOfSound ** 2 +                       &
         contravariantVelocity / speedOfSound)

    leftEigenvectors(1,2) = velocity(1) / temperature
    leftEigenvectors(2,2) = specificVolume * normalizedMetrics(2)
    leftEigenvectors(3,2) = - 0.5_wp * (velocity(1) / temperature -                          &
         normalizedMetrics(1) / speedOfSound)
    leftEigenvectors(4,2) = - 0.5_wp * (velocity(1) / temperature +                          &
         normalizedMetrics(1) / speedOfSound)

    leftEigenvectors(1,3) = velocity(2) / temperature
    leftEigenvectors(2,3) = - specificVolume * normalizedMetrics(1)
    leftEigenvectors(3,3) = - 0.5_wp * (velocity(2) / temperature -                          &
         normalizedMetrics(2) / speedOfSound)
    leftEigenvectors(4,3) = - 0.5_wp * (velocity(2) / temperature +                          &
         normalizedMetrics(2) / speedOfSound)

    leftEigenvectors(1,4) = - 1.0_wp / temperature
    leftEigenvectors(2,4) = 0.0_wp
    leftEigenvectors(3,4) = 0.5_wp / temperature
    leftEigenvectors(4,4) = 0.5_wp / temperature

    ! ``Incoming'' part.
    do j = 1, 4
       do i = 1, 4
          incomingJacobianOfInviscidFlux(i,j) =                                              &
               rightEigenvectors(i,1) * eigenvalues(1) * leftEigenvectors(1,j) +             &
               rightEigenvectors(i,2) * eigenvalues(2) * leftEigenvectors(2,j) +             &
               rightEigenvectors(i,3) * eigenvalues(3) * leftEigenvectors(3,j) +             &
               rightEigenvectors(i,4) * eigenvalues(4) * leftEigenvectors(4,j)
       end do
    end do

    if (present(deltaIncomingJacobianOfInviscidFlux)) then

       ! Variation of the matrix whose columns are the right eigenvectors:

       deltaRightEigenvectors(1,1,:) = 0.0_wp
       deltaRightEigenvectors(2,1,:) = deltaVelocity(1,:)
       deltaRightEigenvectors(3,1,:) = deltaVelocity(2,:)
       deltaRightEigenvectors(4,1,:) = deltaPhiSquared /                                     &
            (ratioOfSpecificHeats - 1.0_wp)

       deltaRightEigenvectors(1,2,:) = 0.0_wp
       deltaRightEigenvectors(2,2,:) = normalizedMetrics(2) * deltaConservedVariables(1,:)
       deltaRightEigenvectors(3,2,:) = -normalizedMetrics(1) * deltaConservedVariables(1,:)
       deltaRightEigenvectors(4,2,:) = deltaConservedVariables(1,:) *                        &
            (normalizedMetrics(2) * velocity(1) - normalizedMetrics(1) * velocity(2)) +      &
            conservedVariables(1) * (normalizedMetrics(2) * deltaVelocity(1,:) -             &
            normalizedMetrics(1) * deltaVelocity(2,:))

       deltaRightEigenvectors(1,3,:) = 0.0_wp
       deltaRightEigenvectors(2,3,:) = deltaVelocity(1,:) +                                  &
            normalizedMetrics(1) * deltaSpeedOfSound
       deltaRightEigenvectors(3,3,:) = deltaVelocity(2,:) +                                  &
            normalizedMetrics(2) * deltaSpeedOfSound
       deltaRightEigenvectors(4,3,:) = deltaTemperature +                                    &
            deltaPhiSquared / (ratioOfSpecificHeats - 1.0_wp) +                              &
            deltaSpeedOfSound * contravariantVelocity +                                      &
            speedOfSound * deltaContravariantVelocity

       deltaRightEigenvectors(1,4,:) = 0.0_wp
       deltaRightEigenvectors(2,4,:) = deltaVelocity(1,:) -                                  &
            normalizedMetrics(1) * deltaSpeedOfSound
       deltaRightEigenvectors(3,4,:) = deltaVelocity(2,:) -                                  &
            normalizedMetrics(2) * deltaSpeedOfSound
       deltaRightEigenvectors(4,4,:) = deltaTemperature +                                    &
            deltaPhiSquared / (ratioOfSpecificHeats - 1.0_wp) -                              &
            deltaSpeedOfSound * contravariantVelocity -                                      &
            speedOfSound * deltaContravariantVelocity

       ! Variation of the matrix whose rows are the left eigenvectors:

       temp = deltaPhiSquared / speedOfSound ** 2 -                                          &
            2.0_wp * phiSquared / speedOfSound ** 3 * deltaSpeedOfSound
       deltaLeftEigenvectors(1,1,:) = -temp
       deltaLeftEigenvectors(2,1,:) = -deltaSpecificVolume *                                 &
            (normalizedMetrics(2) * velocity(1) - normalizedMetrics(1) * velocity(2)) -      &
            specificVolume * (normalizedMetrics(2) * deltaVelocity(1,:) -                    &
            normalizedMetrics(1) * deltaVelocity(2,:))
       deltaLeftEigenvectors(3,1,:) = 0.5_wp * (temp -                                       &
            deltaContravariantVelocity / speedOfSound +                                      &
            contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)
       deltaLeftEigenvectors(4,1,:) = 0.5_wp * (temp +                                       &
            deltaContravariantVelocity / speedOfSound -                                      &
            contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)

       temp = deltaVelocity(1,:) / temperature -                                             &
            velocity(1) / temperature ** 2 * deltaTemperature
       deltaLeftEigenvectors(1,2,:) = temp
       deltaLeftEigenvectors(2,2,:) = deltaSpecificVolume * normalizedMetrics(2)
       deltaLeftEigenvectors(3,2,:) = -0.5_wp * (temp +                                      &
            normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)
       deltaLeftEigenvectors(4,2,:) = -0.5_wp * (temp -                                      &
            normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)

       temp = deltaVelocity(2,:) / temperature -                                             &
            velocity(2) / temperature ** 2 * deltaTemperature
       deltaLeftEigenvectors(1,3,:) = temp
       deltaLeftEigenvectors(2,3,:) = -deltaSpecificVolume * normalizedMetrics(1)
       deltaLeftEigenvectors(3,3,:) = -0.5_wp * (temp +                                      &
            normalizedMetrics(2) / speedOfSound ** 2 * deltaSpeedOfSound)
       deltaLeftEigenvectors(4,3,:) = -0.5_wp * (temp -                                      &
            normalizedMetrics(2) / speedOfSound ** 2 * deltaSpeedOfSound)

       temp =  -1.0_wp / temperature ** 2 * deltaTemperature
       deltaLeftEigenvectors(1,4,:) = -temp
       deltaLeftEigenvectors(2,4,:) = 0.0_wp
       deltaLeftEigenvectors(3,4,:) = 0.5_wp * temp
       deltaLeftEigenvectors(4,4,:) = 0.5_wp * temp

       ! Variation of the ``incoming'' part.
       do j = 1, 4
          do i = 1, 4
             deltaIncomingJacobianOfInviscidFlux(i,j,:) =                                    &
                  deltaRightEigenvectors(i,1,:) * eigenvalues(1) * leftEigenvectors(1,j) +   &
                  deltaRightEigenvectors(i,2,:) * eigenvalues(2) * leftEigenvectors(2,j) +   &
                  deltaRightEigenvectors(i,3,:) * eigenvalues(3) * leftEigenvectors(3,j) +   &
                  deltaRightEigenvectors(i,4,:) * eigenvalues(4) * leftEigenvectors(4,j) +   &
                  rightEigenvectors(i,1) * deltaEigenvalues(1,:) * leftEigenvectors(1,j) +   &
                  rightEigenvectors(i,2) * deltaEigenvalues(2,:) * leftEigenvectors(2,j) +   &
                  rightEigenvectors(i,3) * deltaEigenvalues(3,:) * leftEigenvectors(3,j) +   &
                  rightEigenvectors(i,4) * deltaEigenvalues(4,:) * leftEigenvectors(4,j) +   &
                  rightEigenvectors(i,1) * eigenvalues(1) * deltaLeftEigenvectors(1,j,:) +   &
                  rightEigenvectors(i,2) * eigenvalues(2) * deltaLeftEigenvectors(2,j,:) +   &
                  rightEigenvectors(i,3) * eigenvalues(3) * deltaLeftEigenvectors(3,j,:) +   &
                  rightEigenvectors(i,4) * eigenvalues(4) * deltaLeftEigenvectors(4,j,:)
          end do
       end do

    end if

  end subroutine computeIncomingInviscidJacobian2D

  pure subroutine computeIncomingInviscidJacobian3D(conservedVariables, specificVolume,      &
       velocity, temperature, metrics, ratioOfSpecificHeats, incomingDirection,              &
       incomingJacobianOfInviscidFlux, deltaIncomingJacobianOfInviscidFlux,                  &
       deltaConservedVariables, deltaSpecificVolume, deltaVelocity, deltaTemperature)

    implicit none

    ! <<< Arguments >>>
    real(SCALAR_KIND), intent(in) :: conservedVariables(5), specificVolume, velocity(3),     &
         temperature, metrics(3)
    real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
    integer, intent(in) :: incomingDirection
    real(SCALAR_KIND), intent(out) :: incomingJacobianOfInviscidFlux(5,5)
    real(SCALAR_KIND), intent(out), optional :: deltaIncomingJacobianOfInviscidFlux(5,5,5)
    real(SCALAR_KIND), intent(in), optional :: deltaConservedVariables(5,5),                 &
         deltaSpecificVolume(5), deltaVelocity(3,5), deltaTemperature(5)

    ! <<< Local variables >>>
    integer :: i, j
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: arcLength, normalizedMetrics(3), contravariantVelocity, speedOfSound,        &
         phiSquared, rightEigenvectors(5,5), eigenvalues(5), leftEigenvectors(5,5),          &
         deltaContravariantVelocity(5), deltaSpeedOfSound(5), deltaPhiSquared(5),            &
         deltaRightEigenvectors(5,5,5), deltaEigenvalues(5,5), deltaLeftEigenvectors(5,5,5), &
         temp(5)

    ! Normalize the metrics.
    arcLength = sqrt(metrics(1) ** 2 + metrics(2) ** 2 + metrics(3) ** 2)
    normalizedMetrics = metrics / arcLength

    ! Other dependent variables.
    contravariantVelocity = normalizedMetrics(1) * velocity(1) +                             &
         normalizedMetrics(2) * velocity(2) + normalizedMetrics(3) * velocity(3)
    speedOfSound = sqrt((ratioOfSpecificHeats - 1.0_wp) * temperature)
    phiSquared = 0.5_wp * (ratioOfSpecificHeats - 1.0_wp) *                                  &
         (velocity(1) ** 2 + velocity(2) ** 2 + velocity(3) ** 2)

    ! Eigenvalues.
    eigenvalues(1) = contravariantVelocity
    eigenvalues(2) = contravariantVelocity
    eigenvalues(3) = contravariantVelocity
    eigenvalues(4) = contravariantVelocity + speedOfSound
    eigenvalues(5) = contravariantVelocity - speedOfSound
    eigenvalues = arcLength * eigenvalues

    if (present(deltaIncomingJacobianOfInviscidFlux)) then

       ! Compute variations of other dependent variables.
       deltaContravariantVelocity = normalizedMetrics(1) * deltaVelocity(1,:) +              &
            normalizedMetrics(2) * deltaVelocity(2,:) +                                      &
            normalizedMetrics(3) * deltaVelocity(3,:)
       deltaSpeedOfSound = 0.5_wp / speedOfSound *                                           &
            (ratioOfSpecificHeats - 1.0_wp) * deltaTemperature
       deltaPhiSquared = (ratioOfSpecificHeats - 1.0_wp) *                                   &
            (velocity(1) * deltaVelocity(1,:) + velocity(2) * deltaVelocity(2,:) +           &
            velocity(3) * deltaVelocity(3,:))

       ! Variation of matrix containing eigenvalues.
       deltaEigenvalues(1,:) = deltaContravariantVelocity
       deltaEigenvalues(2,:) = deltaContravariantVelocity
       deltaEigenvalues(3,:) = deltaContravariantVelocity
       deltaEigenvalues(4,:) = deltaContravariantVelocity + deltaSpeedOfSound
       deltaEigenvalues(5,:) = deltaContravariantVelocity - deltaSpeedOfSound
       deltaEigenvalues = arcLength * deltaEigenvalues

    end if

    ! Zero-out the eigenvalues corresponding to outgoing characteristics and corresponding
    ! variations.
    do i = 1, 5
       if (incomingDirection * eigenvalues(i) < 0.0_wp) then
          eigenvalues(i) = 0.0_wp
          if (present(deltaIncomingJacobianOfInviscidFlux)) deltaEigenvalues(i,:) = 0.0_wp
       end if
    end do

    ! Matrix whose columns are the right eigenvectors:

    rightEigenvectors(1,1) = normalizedMetrics(1)
    rightEigenvectors(2,1) = normalizedMetrics(1) * velocity(1)
    rightEigenvectors(3,1) = normalizedMetrics(1) * velocity(2) +                            &
         conservedVariables(1) * normalizedMetrics(3)
    rightEigenvectors(4,1) = normalizedMetrics(1) * velocity(3) -                            &
         conservedVariables(1) * normalizedMetrics(2)
    rightEigenvectors(5,1) = conservedVariables(1) * (normalizedMetrics(3) * velocity(2) -   &
         normalizedMetrics(2) * velocity(3)) + phiSquared /                                  &
         (ratioOfSpecificHeats - 1.0_wp) * normalizedMetrics(1)

    rightEigenvectors(1,2) = normalizedMetrics(2)
    rightEigenvectors(2,2) = normalizedMetrics(2) * velocity(1) -                            &
         conservedVariables(1) * normalizedMetrics(3)
    rightEigenvectors(3,2) = normalizedMetrics(2) * velocity(2)
    rightEigenvectors(4,2) = normalizedMetrics(2) * velocity(3) +                            &
         conservedVariables(1) * normalizedMetrics(1)
    rightEigenvectors(5,2) = conservedVariables(1) * (normalizedMetrics(1) * velocity(3) -   &
         normalizedMetrics(3) * velocity(1))                                                 &
         + phiSquared / (ratioOfSpecificHeats - 1.0_wp) * normalizedMetrics(2)

    rightEigenvectors(1,3) = normalizedMetrics(3)
    rightEigenvectors(2,3) = normalizedMetrics(3) * velocity(1) +                            &
         conservedVariables(1) * normalizedMetrics(2)
    rightEigenvectors(3,3) = normalizedMetrics(3) * velocity(2) -                            &
         conservedVariables(1) * normalizedMetrics(1)
    rightEigenvectors(4,3) = normalizedMetrics(3) * velocity(3)
    rightEigenvectors(5,3) = conservedVariables(1) * (normalizedMetrics(2) * velocity(1) -   &
         normalizedMetrics(1) * velocity(2))                                                 &
         + phiSquared / (ratioOfSpecificHeats - 1.0_wp) * normalizedMetrics(3)

    rightEigenvectors(1,4) = 1.0_wp
    rightEigenvectors(2,4) = velocity(1) + normalizedMetrics(1) * speedOfSound
    rightEigenvectors(3,4) = velocity(2) + normalizedMetrics(2) * speedOfSound
    rightEigenvectors(4,4) = velocity(3) + normalizedMetrics(3) * speedOfSound
    rightEigenvectors(5,4) = temperature + phiSquared / (ratioOfSpecificHeats - 1.0_wp) +    &
         speedOfSound * contravariantVelocity

    rightEigenvectors(1,5) = 1.0_wp
    rightEigenvectors(2,5) = velocity(1) - normalizedMetrics(1) * speedOfSound
    rightEigenvectors(3,5) = velocity(2) - normalizedMetrics(2) * speedOfSound
    rightEigenvectors(4,5) = velocity(3) - normalizedMetrics(3) * speedOfSound
    rightEigenvectors(5,5) = temperature + phiSquared / (ratioOfSpecificHeats - 1.0_wp) -    &
         speedOfSound * contravariantVelocity

    ! Matrix whose rows are the left eigenvectors:

    leftEigenvectors(1,1) = normalizedMetrics(1) * (1.0_wp - phiSquared / speedOfSound ** 2) &
         - specificVolume * (normalizedMetrics(3) * velocity(2) -                            &
         normalizedMetrics(2) * velocity(3))
    leftEigenvectors(2,1) = normalizedMetrics(2) * (1.0_wp - phiSquared / speedOfSound ** 2) &
         - specificVolume * (normalizedMetrics(1) * velocity(3) -                            &
         normalizedMetrics(3) * velocity(1))
    leftEigenvectors(3,1) = normalizedMetrics(3) * (1.0_wp - phiSquared / speedOfSound ** 2) &
         - specificVolume * (normalizedMetrics(2) * velocity(1) -                            &
         normalizedMetrics(1) * velocity(2))
    leftEigenvectors(4,1) = 0.5_wp * (phiSquared / speedOfSound ** 2 -                       &
         contravariantVelocity / speedOfSound)
    leftEigenvectors(5,1) = 0.5_wp * (phiSquared / speedOfSound ** 2 +                       &
         contravariantVelocity / speedOfSound)

    leftEigenvectors(1,2) = normalizedMetrics(1) * velocity(1) / temperature
    leftEigenvectors(2,2) = normalizedMetrics(2) * velocity(1) / temperature -               &
         specificVolume * normalizedMetrics(3)
    leftEigenvectors(3,2) = normalizedMetrics(3) * velocity(1) / temperature +               &
         specificVolume * normalizedMetrics(2)
    leftEigenvectors(4,2) = - 0.5_wp * (velocity(1) / temperature -                          &
         normalizedMetrics(1) / speedOfSound)
    leftEigenvectors(5,2) = - 0.5_wp * (velocity(1) / temperature +                          &
         normalizedMetrics(1) / speedOfSound)

    leftEigenvectors(1,3) = normalizedMetrics(1) * velocity(2) / temperature +               &
         specificVolume * normalizedMetrics(3)
    leftEigenvectors(2,3) = normalizedMetrics(2) * velocity(2) / temperature
    leftEigenvectors(3,3) = normalizedMetrics(3) * velocity(2) / temperature -               &
         specificVolume * normalizedMetrics(1)
    leftEigenvectors(4,3) = - 0.5_wp * (velocity(2) / temperature -                          &
         normalizedMetrics(2) / speedOfSound)
    leftEigenvectors(5,3) = - 0.5_wp * (velocity(2) / temperature +                          &
         normalizedMetrics(2) / speedOfSound)

    leftEigenvectors(1,4) = normalizedMetrics(1) * velocity(3) / temperature -               &
         specificVolume * normalizedMetrics(2)
    leftEigenvectors(2,4) = normalizedMetrics(2) * velocity(3) / temperature +               &
         specificVolume * normalizedMetrics(1)
    leftEigenvectors(3,4) = normalizedMetrics(3) * velocity(3) / temperature
    leftEigenvectors(4,4) = - 0.5_wp * (velocity(3) / temperature -                          &
         normalizedMetrics(3) / speedOfSound)
    leftEigenvectors(5,4) = - 0.5_wp * (velocity(3) / temperature +                          &
         normalizedMetrics(3) / speedOfSound)

    leftEigenvectors(1,5) = - normalizedMetrics(1) / temperature
    leftEigenvectors(2,5) = - normalizedMetrics(2) / temperature
    leftEigenvectors(3,5) = - normalizedMetrics(3) / temperature
    leftEigenvectors(4,5) = 0.5_wp / temperature
    leftEigenvectors(5,5) = 0.5_wp / temperature

    ! ``Incoming'' part.
    do j = 1, 5
       do i = 1, 5
          incomingJacobianOfInviscidFlux(i,j) =                                              &
               rightEigenvectors(i,1) * eigenvalues(1) * leftEigenvectors(1,j) +             &
               rightEigenvectors(i,2) * eigenvalues(2) * leftEigenvectors(2,j) +             &
               rightEigenvectors(i,3) * eigenvalues(3) * leftEigenvectors(3,j) +             &
               rightEigenvectors(i,4) * eigenvalues(4) * leftEigenvectors(4,j) +             &
               rightEigenvectors(i,5) * eigenvalues(5) * leftEigenvectors(5,j)
       end do
    end do

    if (present(deltaIncomingJacobianOfInviscidFlux)) then

       ! Variation of the matrix whose columns are the right eigenvectors:

       deltaRightEigenvectors(1,1,:) = 0.0_wp
       deltaRightEigenvectors(2,1,:) = normalizedMetrics(1) * deltaVelocity(1,:)
       deltaRightEigenvectors(3,1,:) = normalizedMetrics(1) * deltaVelocity(2,:) +           &
            deltaConservedVariables(1,:) * normalizedMetrics(3)
       deltaRightEigenvectors(4,1,:) = normalizedMetrics(1) * deltaVelocity(3,:) -           &
            deltaConservedVariables(1,:) * normalizedMetrics(2)
       deltaRightEigenvectors(5,1,:) = deltaConservedVariables(1,:) *                        &
            (normalizedMetrics(3) * velocity(2) - normalizedMetrics(2) * velocity(3)) +      &
            conservedVariables(1) * (normalizedMetrics(3) * deltaVelocity(2,:) -             &
            normalizedMetrics(2) * deltaVelocity(3,:)) + deltaPhiSquared /                   &
            (ratioOfSpecificHeats - 1.0_wp) * normalizedMetrics(1)

       deltaRightEigenvectors(1,2,:) = 0.0_wp
       deltaRightEigenvectors(2,2,:) = normalizedMetrics(2) * deltaVelocity(1,:) -           &
            deltaConservedVariables(1,:) * normalizedMetrics(3)
       deltaRightEigenvectors(3,2,:) = normalizedMetrics(2) * deltaVelocity(2,:)
       deltaRightEigenvectors(4,2,:) = normalizedMetrics(2) * deltaVelocity(3,:) +           &
            deltaConservedVariables(1,:) * normalizedMetrics(1)
       deltaRightEigenvectors(5,2,:) = deltaConservedVariables(1,:) *                        &
            (normalizedMetrics(1) * velocity(3) - normalizedMetrics(3) * velocity(1)) +      &
            conservedVariables(1) * (normalizedMetrics(1) * deltaVelocity(3,:) -             &
            normalizedMetrics(3) * deltaVelocity(1,:)) + deltaPhiSquared /                   &
            (ratioOfSpecificHeats - 1.0_wp) * normalizedMetrics(2)

       deltaRightEigenvectors(1,3,:) = 0.0_wp
       deltaRightEigenvectors(2,3,:) = normalizedMetrics(3) * deltaVelocity(1,:) +           &
            deltaConservedVariables(1,:) * normalizedMetrics(2)
       deltaRightEigenvectors(3,3,:) = normalizedMetrics(3) * deltaVelocity(2,:) -           &
            deltaConservedVariables(1,:) * normalizedMetrics(1)
       deltaRightEigenvectors(4,3,:) = normalizedMetrics(3) * deltaVelocity(3,:)
       deltaRightEigenvectors(5,3,:) = deltaConservedVariables(1,:) *                        &
            (normalizedMetrics(2) * velocity(1) - normalizedMetrics(1) * velocity(2)) +      &
            conservedVariables(1) * (normalizedMetrics(2) * deltaVelocity(1,:) -             &
            normalizedMetrics(1) * deltaVelocity(2,:)) + deltaPhiSquared /                   &
            (ratioOfSpecificHeats - 1.0_wp) * normalizedMetrics(3)

       deltaRightEigenvectors(1,4,:) = 0.0_wp
       deltaRightEigenvectors(2,4,:) = deltaVelocity(1,:) +                                  &
            normalizedMetrics(1) * deltaSpeedOfSound
       deltaRightEigenvectors(3,4,:) = deltaVelocity(2,:) +                                  &
            normalizedMetrics(2) * deltaSpeedOfSound
       deltaRightEigenvectors(4,4,:) = deltaVelocity(3,:) +                                  &
            normalizedMetrics(3) * deltaSpeedOfSound
       deltaRightEigenvectors(5,4,:) = deltaTemperature +                                    &
            deltaPhiSquared / (ratioOfSpecificHeats - 1.0_wp) +                              &
            deltaSpeedOfSound * contravariantVelocity +                                      &
            speedOfSound * deltaContravariantVelocity

       deltaRightEigenvectors(1,5,:) = 0.0_wp
       deltaRightEigenvectors(2,5,:) = deltaVelocity(1,:) -                                  &
            normalizedMetrics(1) * deltaSpeedOfSound
       deltaRightEigenvectors(3,5,:) = deltaVelocity(2,:) -                                  &
            normalizedMetrics(2) * deltaSpeedOfSound
       deltaRightEigenvectors(4,5,:) = deltaVelocity(3,:) -                                  &
            normalizedMetrics(3) * deltaSpeedOfSound
       deltaRightEigenvectors(5,5,:) = deltaTemperature +                                    &
            deltaPhiSquared / (ratioOfSpecificHeats - 1.0_wp) -                              &
            deltaSpeedOfSound * contravariantVelocity -                                      &
            speedOfSound * deltaContravariantVelocity

       ! Variation of the matrix whose rows are the left eigenvectors:

       temp = deltaPhiSquared / speedOfSound ** 2 -                                          &
            2.0_wp * phiSquared / speedOfSound ** 3 * deltaSpeedOfSound
       deltaLeftEigenvectors(1,1,:) = -normalizedMetrics(1) * temp - deltaSpecificVolume *   &
            (normalizedMetrics(3) * velocity(2) - normalizedMetrics(2) * velocity(3)) -      &
            specificVolume * (normalizedMetrics(3) * deltaVelocity(2,:) -                    &
            normalizedMetrics(2) * deltaVelocity(3,:))
       deltaLeftEigenvectors(2,1,:) = -normalizedMetrics(2) * temp - deltaSpecificVolume *   &
            (normalizedMetrics(1) * velocity(3) - normalizedMetrics(3) * velocity(1)) -      &
            specificVolume * (normalizedMetrics(1) * deltaVelocity(3,:) -                    &
            normalizedMetrics(3) * deltaVelocity(1,:))
       deltaLeftEigenvectors(3,1,:) = -normalizedMetrics(3) * temp - deltaSpecificVolume *   &
            (normalizedMetrics(2) * velocity(1) - normalizedMetrics(1) * velocity(2)) -      &
            specificVolume * (normalizedMetrics(2) * deltaVelocity(1,:) -                    &
            normalizedMetrics(1) * deltaVelocity(2,:))
       deltaLeftEigenvectors(4,1,:) = 0.5_wp * (temp - deltaContravariantVelocity /          &
            speedOfSound + contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)
       deltaLeftEigenvectors(5,1,:) = 0.5_wp * (temp + deltaContravariantVelocity /          &
            speedOfSound - contravariantVelocity / speedOfSound ** 2 * deltaSpeedOfSound)

       temp = deltaVelocity(1,:) / temperature -                                             &
            velocity(1) / temperature ** 2 * deltaTemperature
       deltaLeftEigenvectors(1,2,:) = normalizedMetrics(1) * temp
       deltaLeftEigenvectors(2,2,:) = normalizedMetrics(2) * temp -                          &
            deltaSpecificVolume * normalizedMetrics(3)
       deltaLeftEigenvectors(3,2,:) = normalizedMetrics(3) * temp +                          &
            deltaSpecificVolume * normalizedMetrics(2)
       deltaLeftEigenvectors(4,2,:) = -0.5_wp * (temp +                                      &
            normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)
       deltaLeftEigenvectors(5,2,:) = -0.5_wp * (temp -                                      &
            normalizedMetrics(1) / speedOfSound ** 2 * deltaSpeedOfSound)

       temp = deltaVelocity(2,:) / temperature -                                             &
            velocity(2) / temperature ** 2 * deltaTemperature
       deltaLeftEigenvectors(1,3,:) = normalizedMetrics(1) * temp +                          &
            deltaSpecificVolume * normalizedMetrics(3)
       deltaLeftEigenvectors(2,3,:) = normalizedMetrics(2) * temp
       deltaLeftEigenvectors(3,3,:) = normalizedMetrics(3) * temp -                          &
            deltaSpecificVolume * normalizedMetrics(1)
       deltaLeftEigenvectors(4,3,:) = -0.5_wp * (temp +                                      &
            normalizedMetrics(2) / speedOfSound ** 2 * deltaSpeedOfSound)
       deltaLeftEigenvectors(5,3,:) = -0.5_wp * (temp -                                      &
            normalizedMetrics(2) / speedOfSound ** 2 * deltaSpeedOfSound)

       temp = deltaVelocity(3,:) / temperature -                                             &
            velocity(3) / temperature ** 2 * deltaTemperature
       deltaLeftEigenvectors(1,4,:) = normalizedMetrics(1) * temp -                          &
            deltaSpecificVolume * normalizedMetrics(2)
       deltaLeftEigenvectors(2,4,:) = normalizedMetrics(2) * temp +                          &
            deltaSpecificVolume * normalizedMetrics(1)
       deltaLeftEigenvectors(3,4,:) = normalizedMetrics(3) * temp
       deltaLeftEigenvectors(4,4,:) = -0.5_wp * (temp +                                      &
            normalizedMetrics(3) / speedOfSound ** 2 * deltaSpeedOfSound)
       deltaLeftEigenvectors(5,4,:) = -0.5_wp * (temp -                                      &
            normalizedMetrics(3) / speedOfSound ** 2 * deltaSpeedOfSound)

       temp =  -1.0_wp / temperature ** 2 * deltaTemperature
       deltaLeftEigenvectors(1,5,:) = -normalizedMetrics(1) * temp
       deltaLeftEigenvectors(2,5,:) = -normalizedMetrics(2) * temp
       deltaLeftEigenvectors(3,5,:) = -normalizedMetrics(3) * temp
       deltaLeftEigenvectors(4,5,:) = 0.5_wp * temp
       deltaLeftEigenvectors(5,5,:) = 0.5_wp * temp

       ! Variation of the ``incoming'' part.
       do j = 1, 5
          do i = 1, 5
             deltaIncomingJacobianOfInviscidFlux(i,j,:) =                                    &
                  deltaRightEigenvectors(i,1,:) * eigenvalues(1) * leftEigenvectors(1,j) +   &
                  deltaRightEigenvectors(i,2,:) * eigenvalues(2) * leftEigenvectors(2,j) +   &
                  deltaRightEigenvectors(i,3,:) * eigenvalues(3) * leftEigenvectors(3,j) +   &
                  deltaRightEigenvectors(i,4,:) * eigenvalues(4) * leftEigenvectors(4,j) +   &
                  deltaRightEigenvectors(i,5,:) * eigenvalues(5) * leftEigenvectors(5,j) +   &
                  rightEigenvectors(i,1) * deltaEigenvalues(1,:) * leftEigenvectors(1,j) +   &
                  rightEigenvectors(i,2) * deltaEigenvalues(2,:) * leftEigenvectors(2,j) +   &
                  rightEigenvectors(i,3) * deltaEigenvalues(3,:) * leftEigenvectors(3,j) +   &
                  rightEigenvectors(i,4) * deltaEigenvalues(4,:) * leftEigenvectors(4,j) +   &
                  rightEigenvectors(i,5) * deltaEigenvalues(5,:) * leftEigenvectors(5,j) +   &
                  rightEigenvectors(i,1) * eigenvalues(1) * deltaLeftEigenvectors(1,j,:) +   &
                  rightEigenvectors(i,2) * eigenvalues(2) * deltaLeftEigenvectors(2,j,:) +   &
                  rightEigenvectors(i,3) * eigenvalues(3) * deltaLeftEigenvectors(3,j,:) +   &
                  rightEigenvectors(i,4) * eigenvalues(4) * deltaLeftEigenvectors(4,j,:) +   &
                  rightEigenvectors(i,5) * eigenvalues(5) * deltaLeftEigenvectors(5,j,:)
          end do
       end do

    end if

  end subroutine computeIncomingInviscidJacobian3D

  pure subroutine computeFirstPartialViscousJacobian(nDimensions, ijkIndex,                  &
       conservedVariables, metricsAlongNormalDirection, stressTensor, heatFlux,              &
       powerLawExponent, ratioOfSpecificHeats, firstPartialViscousJacobian, specificVolume,  &
       velocity, temperature)

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: nDimensions, ijkIndex
    real(SCALAR_KIND), intent(in) :: conservedVariables(:,:),                                &
         metricsAlongNormalDirection(:,:), stressTensor(:,:), heatFlux(:,:)
    real(SCALAR_KIND), intent(in) :: powerLawExponent, ratioOfSpecificHeats
    real(SCALAR_KIND), intent(out) :: firstPartialViscousJacobian(:,:)
    real(SCALAR_KIND), intent(in), optional :: specificVolume(:), velocity(:,:),             &
         temperature(:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: nUnknowns
    real(wp), allocatable :: localVelocity(:), localMetricsAlongNormalDirection(:),          &
         localStressTensor(:), localHeatFlux(:)
    real(wp) :: localSpecificVolume, localTemperature

    nUnknowns = size(conservedVariables, 2)

    allocate(localVelocity(nDimensions))
    allocate(localMetricsAlongNormalDirection(nDimensions))
    allocate(localStressTensor(nDimensions ** 2))
    allocate(localHeatFlux(nDimensions))

    localMetricsAlongNormalDirection = metricsAlongNormalDirection(ijkIndex,:)
    localStressTensor = stressTensor(ijkIndex,:)
    localHeatFlux = heatFlux(ijkIndex,:)

    if (present(specificVolume)) then
       localSpecificVolume = specificVolume(ijkIndex)
    else
       localSpecificVolume = 1.0_wp / conservedVariables(ijkIndex, 1)
    end if

    if (present(velocity)) then
       localVelocity = velocity(ijkIndex,:)
    else
       localVelocity = localSpecificVolume * conservedVariables(ijkIndex, 2:nDimensions+1)
    end if

    if (present(temperature)) then
       localTemperature = temperature(ijkIndex)
    else
       localTemperature = ratioOfSpecificHeats * (localSpecificVolume *                      &
            conservedVariables(ijkIndex, nDimensions+2) - 0.5_wp * sum(localVelocity ** 2))
    end if

    select case (nDimensions)
    case (1)
       call computeFirstPartialViscousJacobian1D(localSpecificVolume, localVelocity,         &
            localTemperature, localMetricsAlongNormalDirection, localStressTensor,           &
            localHeatFlux, powerLawExponent, ratioOfSpecificHeats,                           &
            firstPartialViscousJacobian)
    case (2)
       call computeFirstPartialViscousJacobian2D(localSpecificVolume, localVelocity,         &
            localTemperature, localMetricsAlongNormalDirection, localStressTensor,           &
            localHeatFlux, powerLawExponent, ratioOfSpecificHeats,                           &
            firstPartialViscousJacobian)
    case (3)
       call computeFirstPartialViscousJacobian3D(localSpecificVolume, localVelocity,         &
            localTemperature, localMetricsAlongNormalDirection, localStressTensor,           &
            localHeatFlux, powerLawExponent, ratioOfSpecificHeats,                           &
            firstPartialViscousJacobian)
    end select

    SAFE_DEALLOCATE(localHeatFlux)
    SAFE_DEALLOCATE(localStressTensor)
    SAFE_DEALLOCATE(localMetricsAlongNormalDirection)
    SAFE_DEALLOCATE(localVelocity)

  end subroutine computeFirstPartialViscousJacobian

  pure subroutine computeFirstPartialViscousJacobian1D(specificVolume, velocity,             &
       temperature, metrics, stressTensor, heatFlux, powerLawExponent, ratioOfSpecificHeats, &
       firstPartialViscousJacobian)

    implicit none

    ! <<< Arguments >>>
    real(SCALAR_KIND), intent(in) :: specificVolume, velocity(1), temperature, metrics(1),   &
         stressTensor(1), heatFlux(1), powerLawExponent, ratioOfSpecificHeats
    real(SCALAR_KIND), intent(out) :: firstPartialViscousJacobian(3,3)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: phiSquared, contravariantStressTensor(1), contravariantHeatFlux, temp1, temp2

    ! Other dependent variables.
    phiSquared = 0.5_wp * (ratioOfSpecificHeats - 1.0_wp) * velocity(1) ** 2
    contravariantStressTensor(1) = metrics(1) * stressTensor(1) !... not normalized.
    contravariantHeatFlux = metrics(1) * heatFlux(1) !... not normalized.
    temp1 = velocity(1) * contravariantStressTensor(1) - contravariantHeatFlux

    temp2 = powerLawExponent * ratioOfSpecificHeats * specificVolume / temperature *         &
         (phiSquared / (ratioOfSpecificHeats - 1.0_wp) - temperature / ratioOfSpecificHeats)
    firstPartialViscousJacobian(1,1) = 0.0_wp
    firstPartialViscousJacobian(2,1) = temp2 * contravariantStressTensor(1)
    firstPartialViscousJacobian(3,1) = temp2 * temp1 - specificVolume *                      &
         (velocity(1) * contravariantStressTensor(1))

    temp2 = - powerLawExponent * ratioOfSpecificHeats *                                      &
         specificVolume / temperature * velocity(1)
    firstPartialViscousJacobian(1,2) = 0.0_wp
    firstPartialViscousJacobian(2,2) = temp2 * contravariantStressTensor(1)
    firstPartialViscousJacobian(3,2) = temp2 * temp1 +                                       &
         specificVolume * contravariantStressTensor(1)

    temp2 = powerLawExponent * ratioOfSpecificHeats * specificVolume / temperature
    firstPartialViscousJacobian(1,3) = 0.0_wp
    firstPartialViscousJacobian(2,3) = temp2 * contravariantStressTensor(1)
    firstPartialViscousJacobian(3,3) = temp2 * temp1

  end subroutine computeFirstPartialViscousJacobian1D

  pure subroutine computeFirstPartialViscousJacobian2D(specificVolume, velocity,             &
       temperature, metrics, stressTensor, heatFlux, powerLawExponent, ratioOfSpecificHeats, &
       firstPartialViscousJacobian)

    implicit none

    ! <<< Arguments >>>
    real(SCALAR_KIND), intent(in) :: specificVolume, velocity(2), temperature, metrics(2),   &
         stressTensor(4), heatFlux(2), powerLawExponent, ratioOfSpecificHeats
    real(SCALAR_KIND), intent(out) :: firstPartialViscousJacobian(4,4)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: phiSquared, contravariantStressTensor(2), contravariantHeatFlux, temp1, temp2

    ! Other dependent variables.
    phiSquared = 0.5_wp * (ratioOfSpecificHeats - 1.0_wp) *                                  &
         (velocity(1) ** 2 + velocity(2) ** 2)
    contravariantStressTensor(1) = metrics(1) * stressTensor(1) +                            &
         metrics(2) * stressTensor(2) !... not normalized.
    contravariantStressTensor(2) = metrics(1) * stressTensor(3) +                            &
         metrics(2) * stressTensor(4) !... not normalized.
    contravariantHeatFlux = metrics(1) * heatFlux(1) +                                       &
         metrics(2) * heatFlux(2) !... not normalized.
    temp1 = velocity(1) * contravariantStressTensor(1) +                                     &
         velocity(2) * contravariantStressTensor(2) - contravariantHeatFlux

    temp2 = powerLawExponent * ratioOfSpecificHeats * specificVolume / temperature *         &
         (phiSquared / (ratioOfSpecificHeats - 1.0_wp) - temperature / ratioOfSpecificHeats)
    firstPartialViscousJacobian(1,1) = 0.0_wp
    firstPartialViscousJacobian(2,1) = temp2 * contravariantStressTensor(1)
    firstPartialViscousJacobian(3,1) = temp2 * contravariantStressTensor(2)
    firstPartialViscousJacobian(4,1) = temp2 * temp1 - specificVolume *                      &
         (velocity(1) * contravariantStressTensor(1) +                                       &
         velocity(2) * contravariantStressTensor(2))

    temp2 = - powerLawExponent * ratioOfSpecificHeats *                                      &
         specificVolume / temperature * velocity(1)
    firstPartialViscousJacobian(1,2) = 0.0_wp
    firstPartialViscousJacobian(2,2) = temp2 * contravariantStressTensor(1)
    firstPartialViscousJacobian(3,2) = temp2 * contravariantStressTensor(2)
    firstPartialViscousJacobian(4,2) = temp2 * temp1 +                                       &
         specificVolume * contravariantStressTensor(1)

    temp2 = - powerLawExponent * ratioOfSpecificHeats *                                      &
         specificVolume / temperature * velocity(2)
    firstPartialViscousJacobian(1,3) = 0.0_wp
    firstPartialViscousJacobian(2,3) = temp2 * contravariantStressTensor(1)
    firstPartialViscousJacobian(3,3) = temp2 * contravariantStressTensor(2)
    firstPartialViscousJacobian(4,3) = temp2 * temp1 +                                       &
         specificVolume * contravariantStressTensor(2)

    temp2 = powerLawExponent * ratioOfSpecificHeats * specificVolume / temperature
    firstPartialViscousJacobian(1,4) = 0.0_wp
    firstPartialViscousJacobian(2,4) = temp2 * contravariantStressTensor(1)
    firstPartialViscousJacobian(3,4) = temp2 * contravariantStressTensor(2)
    firstPartialViscousJacobian(4,4) = temp2 * temp1

  end subroutine computeFirstPartialViscousJacobian2D

  pure subroutine computeFirstPartialViscousJacobian3D(specificVolume, velocity,             &
       temperature, metrics, stressTensor, heatFlux, powerLawExponent, ratioOfSpecificHeats, &
       firstPartialViscousJacobian)

    implicit none

    ! <<< Arguments >>>
    real(SCALAR_KIND), intent(in) :: specificVolume, velocity(3), temperature, metrics(3),   &
         stressTensor(9), heatFlux(3), powerLawExponent, ratioOfSpecificHeats
    real(SCALAR_KIND), intent(out) :: firstPartialViscousJacobian(5,5)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: phiSquared, contravariantStressTensor(3), contravariantHeatFlux, temp1, temp2

    ! Other dependent variables.
    phiSquared = 0.5_wp * (ratioOfSpecificHeats - 1.0_wp) *                                  &
         (velocity(1) ** 2 + velocity(2) ** 2 + velocity(3) ** 2)
    contravariantStressTensor(1) = metrics(1) * stressTensor(1) +                            &
         metrics(2) * stressTensor(2) + metrics(3) * stressTensor(3) !... not normalized.
    contravariantStressTensor(2) = metrics(1) * stressTensor(4) +                            &
         metrics(2) * stressTensor(5) + metrics(3) * stressTensor(6) !... not normalized.
    contravariantStressTensor(3) = metrics(1) * stressTensor(7) +                            &
         metrics(2) * stressTensor(8) + metrics(3) * stressTensor(9) !... not normalized.
    contravariantHeatFlux = metrics(1) * heatFlux(1) + metrics(2) * heatFlux(2) +            &
         metrics(3) * heatFlux(3) !... not normalized.
    temp1 = velocity(1) * contravariantStressTensor(1) +                                     &
         velocity(2) * contravariantStressTensor(2) +                                        &
         velocity(3) * contravariantStressTensor(3) - contravariantHeatFlux

    temp2 = powerLawExponent * ratioOfSpecificHeats * specificVolume / temperature *         &
         (phiSquared / (ratioOfSpecificHeats - 1.0_wp) - temperature / ratioOfSpecificHeats)
    firstPartialViscousJacobian(1,1) = 0.0_wp
    firstPartialViscousJacobian(2,1) = temp2 * contravariantStressTensor(1)
    firstPartialViscousJacobian(3,1) = temp2 * contravariantStressTensor(2)
    firstPartialViscousJacobian(4,1) = temp2 * contravariantStressTensor(3)
    firstPartialViscousJacobian(5,1) = temp2 * temp1 - specificVolume *                      &
         (velocity(1) * contravariantStressTensor(1) +                                       &
         velocity(2) * contravariantStressTensor(2) +                                        &
         velocity(3) * contravariantStressTensor(3))

    temp2 = - powerLawExponent * ratioOfSpecificHeats *                                      &
         specificVolume / temperature * velocity(1)
    firstPartialViscousJacobian(1,2) = 0.0_wp
    firstPartialViscousJacobian(2,2) = temp2 * contravariantStressTensor(1)
    firstPartialViscousJacobian(3,2) = temp2 * contravariantStressTensor(2)
    firstPartialViscousJacobian(4,2) = temp2 * contravariantStressTensor(3)
    firstPartialViscousJacobian(5,2) = temp2 * temp1 +                                       &
         specificVolume * contravariantStressTensor(1)

    temp2 = - powerLawExponent * ratioOfSpecificHeats *                                      &
         specificVolume / temperature * velocity(2)
    firstPartialViscousJacobian(1,3) = 0.0_wp
    firstPartialViscousJacobian(2,3) = temp2 * contravariantStressTensor(1)
    firstPartialViscousJacobian(3,3) = temp2 * contravariantStressTensor(2)
    firstPartialViscousJacobian(4,3) = temp2 * contravariantStressTensor(3)
    firstPartialViscousJacobian(5,3) = temp2 * temp1 +                                       &
         specificVolume * contravariantStressTensor(2)

    temp2 = - powerLawExponent * ratioOfSpecificHeats *                                      &
         specificVolume / temperature * velocity(3)
    firstPartialViscousJacobian(1,4) = 0.0_wp
    firstPartialViscousJacobian(2,4) = temp2 * contravariantStressTensor(1)
    firstPartialViscousJacobian(3,4) = temp2 * contravariantStressTensor(2)
    firstPartialViscousJacobian(4,4) = temp2 * contravariantStressTensor(3)
    firstPartialViscousJacobian(5,4) = temp2 * temp1 +                                       &
         specificVolume * contravariantStressTensor(3)

    temp2 = powerLawExponent * ratioOfSpecificHeats * specificVolume / temperature
    firstPartialViscousJacobian(1,5) = 0.0_wp
    firstPartialViscousJacobian(2,5) = temp2 * contravariantStressTensor(1)
    firstPartialViscousJacobian(3,5) = temp2 * contravariantStressTensor(2)
    firstPartialViscousJacobian(4,5) = temp2 * contravariantStressTensor(3)
    firstPartialViscousJacobian(5,5) = temp2 * temp1

  end subroutine computeFirstPartialViscousJacobian3D

  pure subroutine computeSecondPartialViscousJacobian(nDimensions, ijkIndex, velocity,       &
       dynamicViscosity, secondCoefficientOfViscosity, thermalDiffusivity, jacobian,         &
       metricsAlongNormalDirection, metricsAlongOtherDirection, secondPartialViscousJacobian)

    implicit none

    ! <<< Arguments >>>
    integer, intent(in) :: nDimensions, ijkIndex
    real(SCALAR_KIND), intent(in) :: velocity(:,:), dynamicViscosity(:),                     &
         secondCoefficientOfViscosity(:), thermalDiffusivity(:), jacobian(:),                &
         metricsAlongNormalDirection(:,:), metricsAlongOtherDirection(:,:)
    real(SCALAR_KIND), intent(out) :: secondPartialViscousJacobian(:,:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp), allocatable :: localVelocity(:), localMetricsAlongNormalDirection(:),          &
         localMetricsAlongOtherDirection(:)

    allocate(localVelocity(nDimensions))
    allocate(localMetricsAlongNormalDirection(nDimensions))
    allocate(localMetricsAlongOtherDirection(nDimensions))

    localVelocity = velocity(ijkIndex,:)
    localMetricsAlongNormalDirection = metricsAlongNormalDirection(ijkIndex,:)
    localMetricsAlongOtherDirection = metricsAlongOtherDirection(ijkIndex,:)

    select case (nDimensions)
    case (1)
       call computeSecondPartialViscousJacobian1D(localVelocity, dynamicViscosity(ijkIndex), &
            secondCoefficientOfViscosity(ijkIndex), thermalDiffusivity(ijkIndex),            &
            jacobian(ijkIndex), localMetricsAlongNormalDirection,                            &
            secondPartialViscousJacobian)
    case (2)
       call computeSecondPartialViscousJacobian2D(localVelocity, dynamicViscosity(ijkIndex), &
            secondCoefficientOfViscosity(ijkIndex), thermalDiffusivity(ijkIndex),            &
            jacobian(ijkIndex), localMetricsAlongNormalDirection,                            &
            localMetricsAlongOtherDirection, secondPartialViscousJacobian)
    case (3)
       call computeSecondPartialViscousJacobian3D(localVelocity, dynamicViscosity(ijkIndex), &
            secondCoefficientOfViscosity(ijkIndex), thermalDiffusivity(ijkIndex),            &
            jacobian(ijkIndex), localMetricsAlongNormalDirection,                            &
            localMetricsAlongOtherDirection, secondPartialViscousJacobian)
    end select

    SAFE_DEALLOCATE(localMetricsAlongOtherDirection)
    SAFE_DEALLOCATE(localMetricsAlongNormalDirection)
    SAFE_DEALLOCATE(localVelocity)

  end subroutine computeSecondPartialViscousJacobian

  pure subroutine computeSecondPartialViscousJacobian1D(velocity, dynamicViscosity,          &
       secondCoefficientOfViscosity, thermalDiffusivity, jacobian, metrics,                  &
       secondPartialViscousJacobian)

    implicit none

    ! <<< Arguments >>>
    real(SCALAR_KIND), intent(in) :: velocity(1), dynamicViscosity,                          &
         secondCoefficientOfViscosity, thermalDiffusivity, jacobian, metrics(1)
    real(SCALAR_KIND), intent(out) :: secondPartialViscousJacobian(2,2)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: temp1, temp2, temp3

    ! Temporary variables.
    temp1 = metrics(1) * metrics(1)
    temp2 = dynamicViscosity * metrics(1) * velocity(1)
    temp3 = secondCoefficientOfViscosity * metrics(1) * velocity(1)

    secondPartialViscousJacobian(1,1) = dynamicViscosity * temp1 +                           &
         (dynamicViscosity + secondCoefficientOfViscosity) *                                 &
         metrics(1) * metrics(1)
    secondPartialViscousJacobian(2,1) = dynamicViscosity * temp1 * velocity(1) +             &
         metrics(1) * temp2 + metrics(1) * temp3

    secondPartialViscousJacobian(1,2) = 0.0_wp
    secondPartialViscousJacobian(2,2) = thermalDiffusivity * temp1

    ! Multiply by the Jacobian.
    secondPartialViscousJacobian = jacobian * secondPartialViscousJacobian

  end subroutine computeSecondPartialViscousJacobian1D

  pure subroutine computeSecondPartialViscousJacobian2D(velocity, dynamicViscosity,          &
       secondCoefficientOfViscosity, thermalDiffusivity, jacobian, metricsAlongFirstDir,     &
       metricsAlongSecondDir, secondPartialViscousJacobian)

    implicit none

    ! <<< Arguments >>>
    real(SCALAR_KIND), intent(in) :: velocity(2), dynamicViscosity,                          &
         secondCoefficientOfViscosity, thermalDiffusivity, jacobian,                         &
         metricsAlongFirstDir(2), metricsAlongSecondDir(2)
    real(SCALAR_KIND), intent(out) :: secondPartialViscousJacobian(3,3)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: temp1, temp2, temp3

    ! Temporary variables.
    temp1 = metricsAlongFirstDir(1) * metricsAlongSecondDir(1) +                             &
         metricsAlongFirstDir(2) * metricsAlongSecondDir(2)
    temp2 = dynamicViscosity * (metricsAlongSecondDir(1) * velocity(1) +                     &
         metricsAlongSecondDir(2) * velocity(2))
    temp3 = secondCoefficientOfViscosity * (metricsAlongFirstDir(1) * velocity(1) +          &
         metricsAlongFirstDir(2) * velocity(2))

    secondPartialViscousJacobian(1,1) = dynamicViscosity * temp1 +                           &
         (dynamicViscosity + secondCoefficientOfViscosity) *                                 &
         metricsAlongFirstDir(1) * metricsAlongSecondDir(1)
    secondPartialViscousJacobian(2,1) =                                                      &
         dynamicViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(2) +             &
         secondCoefficientOfViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(1)
    secondPartialViscousJacobian(3,1) = dynamicViscosity * temp1 * velocity(1) +             &
         metricsAlongFirstDir(1) * temp2 + metricsAlongSecondDir(1) * temp3

    secondPartialViscousJacobian(1,2) =                                                      &
         dynamicViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(1) +             &
         secondCoefficientOfViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(2)
    secondPartialViscousJacobian(2,2) = dynamicViscosity * temp1 +                           &
         (dynamicViscosity + secondCoefficientOfViscosity) *                                 &
         metricsAlongFirstDir(2) * metricsAlongSecondDir(2)
    secondPartialViscousJacobian(3,2) = dynamicViscosity * temp1 * velocity(2) +             &
         metricsAlongFirstDir(2) * temp2 + metricsAlongSecondDir(2) * temp3

    secondPartialViscousJacobian(1,3) = 0.0_wp
    secondPartialViscousJacobian(2,3) = 0.0_wp
    secondPartialViscousJacobian(3,3) = thermalDiffusivity * temp1

    ! Multiply by the Jacobian.
    secondPartialViscousJacobian = jacobian * secondPartialViscousJacobian

  end subroutine computeSecondPartialViscousJacobian2D

  pure subroutine computeSecondPartialViscousJacobian3D(velocity, dynamicViscosity,          &
       secondCoefficientOfViscosity, thermalDiffusivity, jacobian, metricsAlongFirstDir,     &
       metricsAlongSecondDir, secondPartialViscousJacobian)

    implicit none

    ! <<< Arguments >>>
    real(SCALAR_KIND), intent(in) :: velocity(3), dynamicViscosity,                          &
         secondCoefficientOfViscosity, thermalDiffusivity, jacobian,                         &
         metricsAlongFirstDir(3), metricsAlongSecondDir(3)
    real(SCALAR_KIND), intent(out) :: secondPartialViscousJacobian(4,4)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    real(wp) :: temp1, temp2, temp3

    ! Temporary variables.
    temp1 = metricsAlongFirstDir(1) * metricsAlongSecondDir(1) +                             &
         metricsAlongFirstDir(2) * metricsAlongSecondDir(2) +                                &
         metricsAlongFirstDir(3) * metricsAlongSecondDir(3)
    temp2 = dynamicViscosity * (metricsAlongSecondDir(1) * velocity(1) +                     &
         metricsAlongSecondDir(2) * velocity(2) + metricsAlongSecondDir(3) * velocity(3))
    temp3 = secondCoefficientOfViscosity * (metricsAlongFirstDir(1) * velocity(1) +          &
         metricsAlongFirstDir(2) * velocity(2) + metricsAlongFirstDir(3) * velocity(3))

    secondPartialViscousJacobian(1,1) = dynamicViscosity * temp1 +                           &
         (dynamicViscosity + secondCoefficientOfViscosity) *                                 &
         metricsAlongFirstDir(1) * metricsAlongSecondDir(1)
    secondPartialViscousJacobian(2,1) =                                                      &
         dynamicViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(2) +             &
         secondCoefficientOfViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(1)
    secondPartialViscousJacobian(3,1) =                                                      &
         dynamicViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(3) +             &
         secondCoefficientOfViscosity * metricsAlongFirstDir(3) * metricsAlongSecondDir(1)
    secondPartialViscousJacobian(4,1) = dynamicViscosity * temp1 * velocity(1) +             &
         metricsAlongFirstDir(1) * temp2 + metricsAlongSecondDir(1) * temp3

    secondPartialViscousJacobian(1,2) =                                                      &
         dynamicViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(1) +             &
         secondCoefficientOfViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(2)
    secondPartialViscousJacobian(2,2) = dynamicViscosity * temp1 +                           &
         (dynamicViscosity + secondCoefficientOfViscosity) *                                 &
         metricsAlongFirstDir(2) * metricsAlongSecondDir(2)
    secondPartialViscousJacobian(3,2) =                                                      &
         dynamicViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(3) +             &
         secondCoefficientOfViscosity * metricsAlongFirstDir(3) * metricsAlongSecondDir(2)
    secondPartialViscousJacobian(4,2) = dynamicViscosity * temp1 * velocity(2) +             &
         metricsAlongFirstDir(2) * temp2 + metricsAlongSecondDir(2) * temp3

    secondPartialViscousJacobian(1,3) =                                                      &
         dynamicViscosity * metricsAlongFirstDir(3) * metricsAlongSecondDir(1) +             &
         secondCoefficientOfViscosity * metricsAlongFirstDir(1) * metricsAlongSecondDir(3)
    secondPartialViscousJacobian(2,3) =                                                      &
         dynamicViscosity * metricsAlongFirstDir(3) * metricsAlongSecondDir(2) +             &
         secondCoefficientOfViscosity * metricsAlongFirstDir(2) * metricsAlongSecondDir(3)
    secondPartialViscousJacobian(3,3) = dynamicViscosity * temp1 +                           &
         (dynamicViscosity + secondCoefficientOfViscosity) *                                 &
         metricsAlongFirstDir(3) * metricsAlongSecondDir(3)
    secondPartialViscousJacobian(4,3) = dynamicViscosity * temp1 * velocity(3) +             &
         metricsAlongFirstDir(3) * temp2 + metricsAlongSecondDir(3) * temp3

    secondPartialViscousJacobian(1,4) = 0.0_wp
    secondPartialViscousJacobian(2,4) = 0.0_wp
    secondPartialViscousJacobian(3,4) = 0.0_wp
    secondPartialViscousJacobian(4,4) = thermalDiffusivity * temp1

    ! Multiply by the Jacobian.
    secondPartialViscousJacobian = jacobian * secondPartialViscousJacobian

  end subroutine computeSecondPartialViscousJacobian3D

end module CNSHelper
