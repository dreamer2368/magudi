#include "config.h"

module CNSHelper

  implicit none
  public

  interface

     pure subroutine computeDependentVariables(nDimensions, nSpecies, conservedVariables,    &
          ratioOfSpecificHeats, specificVolume,                                              &
          velocity, pressure, temperature, massFraction)

       !> Computes the requested dependent variable(s) including specific volume, velocity,
       !> pressure, temperature, and mass fraction from the conserved state variables.

       integer, intent(in) :: nDimensions, nSpecies
       SCALAR_TYPE, intent(in) :: conservedVariables(:,:)

       real(SCALAR_KIND), intent(in), optional :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(out), optional :: specificVolume(:),                              &
            velocity(:,:), pressure(:), temperature(:), massFraction(:,:)

     end subroutine computeDependentVariables

  end interface

  interface

     pure subroutine computeTransportVariables(temperature, powerLawExponent,                &
          bulkViscosityRatio, ratioOfSpecificHeats, reynoldsNumberInverse,                   &
          prandtlNumberInverse, schmidtNumberInverse, dynamicViscosity,                      &
          secondCoefficientOfViscosity, thermalDiffusivity)

       !> Computes the requested transport coefficient(s) including the dynamic viscosity,
       !> second coefficient of viscosity and the thermal conductivity assuming a power law
       !> dependence on temperature with exponent `powerLawExponent`.

       SCALAR_TYPE, intent(in) :: temperature(:)
       real(SCALAR_KIND), intent(in) :: powerLawExponent,                                    &
            ratioOfSpecificHeats, reynoldsNumberInverse

       real(SCALAR_KIND), intent(in), optional :: bulkViscosityRatio, prandtlNumberInverse,  &
            schmidtNumberInverse
       SCALAR_TYPE, intent(out), optional :: dynamicViscosity(:),                            &
            secondCoefficientOfViscosity(:),                                                 &
            thermalDiffusivity(:)

     end subroutine computeTransportVariables

  end interface

  interface

     pure subroutine computeRoeAverage(nDimensions, conservedVariablesL,                     &
          conservedVariablesR, ratioOfSpecificHeats, roeAverage)

       integer, intent(in) :: nDimensions
       SCALAR_TYPE, intent(in) :: conservedVariablesL(:), conservedVariablesR(:)
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(out) :: roeAverage(:)

     end subroutine computeRoeAverage

  end interface

  interface

     pure subroutine computeStressTensor(nDimensions, nSpecies, velocityGradient,            &
       dynamicViscosity, secondCoefficientOfViscosity, stressTensor, massFraction)

       integer, intent(in) :: nDimensions, nSpecies
       SCALAR_TYPE, intent(inout) :: velocityGradient(:,:)
       SCALAR_TYPE, intent(in) :: dynamicViscosity(:), secondCoefficientOfViscosity(:)

       SCALAR_TYPE, intent(out), optional :: stressTensor(:,:)

       SCALAR_TYPE, intent(inout), optional :: massFraction(:,:)

     end subroutine computeStressTensor

  end interface

  interface

     pure subroutine computeVorticityMagnitudeAndDilatation(nDimensions,                     &
          velocityGradient, vorticityMagnitude, dilatation)

       !> Computes the vorticity magnitude and dilatation from the velocity gradient.

       integer, intent(in) :: nDimensions
       SCALAR_TYPE, intent(in) :: velocityGradient(:,:)
       SCALAR_TYPE, intent(out), optional :: vorticityMagnitude(:), dilatation(:)

     end subroutine computeVorticityMagnitudeAndDilatation

  end interface

  interface

     pure subroutine computeCartesianInvsicidFluxes(nDimensions, nSpecies,                   &
       conservedVariables, velocity, pressure, massFraction, inviscidFluxes)

       !> Computes the Cartesian form of the inviscid fluxes. `velocity`, `pressure`, and
       !> `mass fraction` must be consistent with `conservedVariables`, otherwise the result
       !> is unpredictable.

       integer, intent(in) :: nDimensions, nSpecies
       SCALAR_TYPE, intent(in) :: conservedVariables(:,:), velocity(:,:), pressure(:)
       SCALAR_TYPE, intent(out) :: inviscidFluxes(:,:,:)
       SCALAR_TYPE, intent(in), optional :: massFraction(:,:)

     end subroutine computeCartesianInvsicidFluxes

  end interface

  interface

     pure subroutine computeCartesianViscousFluxes(nDimensions, nSpecies, velocity,          &
       massFraction,stressTensor, heatFlux, speciesFlux, viscousFluxes)

       !> Computes the Cartesian form of the viscous fluxes.

       integer, intent(in) :: nDimensions, nSpecies
       SCALAR_TYPE, intent(in) :: velocity(:,:), stressTensor(:,:), heatFlux(:,:)
       SCALAR_TYPE, intent(in), optional :: massFraction(:,:), speciesFlux(:,:,:)
       SCALAR_TYPE, intent(out) :: viscousFluxes(:,:,:)

     end subroutine computeCartesianViscousFluxes

  end interface

  interface

     pure subroutine computeSpectralRadius(nDimensions, nSpecies, ratioOfSpecificHeats,      &
       velocity, temperature, massFraction, metrics, spectralRadius, isDomainCurvilinear)

       !> Compute the spectral radii along all the directions.

       integer, intent(in) :: nDimensions, nSpecies
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(in) :: velocity(:,:), temperature(:), metrics(:,:)
       SCALAR_TYPE, intent(out) :: spectralRadius(:,:)
       SCALAR_TYPE, intent(in), optional :: massFraction(:,:)
       logical, intent(in), optional :: isDomainCurvilinear

     end subroutine computeSpectralRadius

  end interface

  interface

     pure subroutine transformFluxes(nDimensions, nSpecies, fluxes, metrics,                 &
          transformedFluxes, isDomainCurvilinear)

       !> Transforms fluxes from Cartesian form to contravariant form. If the domain is not
       !> curvilinear as specified by `isDomainCurvilinear`, the arguments `fluxes` and
       !> `transformedFluxes` may reference the same array.

       integer, intent(in) :: nDimensions, nSpecies
       SCALAR_TYPE, intent(in) :: fluxes(:,:,:), metrics(:,:)
       SCALAR_TYPE, intent(out) :: transformedFluxes(:,:,:)
       logical, intent(in), optional :: isDomainCurvilinear

     end subroutine transformFluxes

  end interface

  interface

     pure function computeCfl(nDimensions, iblank, jacobian, metrics,                        &
          velocity, temperature, timeStepSize, ratioOfSpecificHeats,                         &
          dynamicViscosity, thermalDiffusivity) result(cfl)

       !> Computes the CFL number.

       integer, intent(in) :: nDimensions, iblank(:)
       SCALAR_TYPE, intent(in) :: jacobian(:), metrics(:,:), velocity(:,:), temperature(:)
       real(SCALAR_KIND), intent(in) :: timeStepSize, ratioOfSpecificHeats
       SCALAR_TYPE, intent(in), optional :: dynamicViscosity(:), thermalDiffusivity(:)

       real(SCALAR_KIND) :: cfl

     end function computeCfl

  end interface

  interface

     pure function computeTimeStepSize(nDimensions, iblank, jacobian,                        &
          metrics, velocity, temperature, cfl,                                               &
          ratioOfSpecificHeats, dynamicViscosity,                                            &
          thermalDiffusivity) result(timeStepSize)

       !> Computes the time step size that leads to a specified CFL number.

       integer, intent(in) :: nDimensions, iblank(:)
       SCALAR_TYPE, intent(in) :: jacobian(:), metrics(:,:),                                 &
            velocity(:,:), temperature(:)
       real(SCALAR_KIND), intent(in) :: cfl, ratioOfSpecificHeats
       SCALAR_TYPE, intent(in), optional :: dynamicViscosity(:), thermalDiffusivity(:)

       real(SCALAR_KIND) :: timeStepSize

     end function computeTimeStepSize

  end interface

  interface
     pure subroutine computeJacobianOfInviscidFlux(nDimensions, nSpecies,                    &
          conservedVariables, metrics, ratioOfSpecificHeats, jacobianOfInviscidFlux,         &
          deltaConservedVariables, specificVolume, velocity, temperature,                    &
          massFraction, deltaJacobianOfInviscidFlux)

       integer, intent(in) :: nDimensions, nSpecies
       SCALAR_TYPE, intent(in) :: conservedVariables(nDimensions + nSpecies + 2),            &
            metrics(nDimensions)
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(out) :: jacobianOfInviscidFlux                                    &
            (nDimensions + nSpecies + 2,nDimensions + nSpecies + 2)

       SCALAR_TYPE, intent(in), optional :: deltaConservedVariables                          &
            (nDimensions + nSpecies + 2,nDimensions + nSpecies + 2),                         &
            specificVolume, velocity(nDimensions), temperature, massFraction(nSpecies)
       SCALAR_TYPE, intent(out), optional :: deltaJacobianOfInviscidFlux                     &
            (nDimensions + nSpecies + 2,nDimensions + nSpecies + 2,nDimensions + nSpecies + 2)

     end subroutine computeJacobianOfInviscidFlux

  end interface

  interface

     pure subroutine computeIncomingJacobianOfInviscidFlux(nDimensions, nSpecies,            &
          conservedVariables, metrics, ratioOfSpecificHeats, incomingDirection,              &
          incomingJacobianOfInviscidFlux, deltaIncomingJacobianOfInviscidFlux,               &
          deltaConservedVariables, specificVolume, velocity, temperature, massFraction)

       integer, intent(in) :: nDimensions, nSpecies
       SCALAR_TYPE, intent(in) :: conservedVariables(nDimensions + nSpecies + 2),            &
            metrics(nDimensions)
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
       integer, intent(in) :: incomingDirection
       SCALAR_TYPE, intent(out) :: incomingJacobianOfInviscidFlux                            &
            (nDimensions + nSpecies + 2,nDimensions + nSpecies + 2)

       SCALAR_TYPE, intent(out), optional :: deltaIncomingJacobianOfInviscidFlux             &
            (nDimensions + nSpecies + 2,nDimensions + nSpecies + 2,nDimensions + nSpecies + 2)
       SCALAR_TYPE, intent(in), optional :: deltaConservedVariables                          &
            (nDimensions + nSpecies + 2,nDimensions + nSpecies + 2), specificVolume,         &
            velocity(nDimensions), temperature, massFraction(nSpecies)

     end subroutine computeIncomingJacobianOfInviscidFlux

  end interface

  interface

     pure subroutine computeFirstPartialViscousJacobian(nDimensions, nSpecies,               &
          conservedVariables, metrics, stressTensor, heatFlux, speciesFlux,                  &
          powerLawExponent, ratioOfSpecificHeats, firstPartialViscousJacobian,               &
          specificVolume, velocity, temperature, massFraction)

       integer, intent(in) :: nDimensions, nSpecies
       SCALAR_TYPE, intent(in) :: conservedVariables(nDimensions + nSpecies + 2),            &
            metrics(nDimensions), stressTensor(nDimensions), heatFlux(nDimensions)
       real(SCALAR_KIND), intent(in) :: powerLawExponent, ratioOfSpecificHeats
       SCALAR_TYPE, intent(out) :: firstPartialViscousJacobian                               &
            (nDimensions + nSpecies + 2,nDimensions + nSpecies + 2)

       SCALAR_TYPE, intent(in), optional :: specificVolume, velocity(nDimensions),           &
            temperature, massFraction(nSpecies), speciesFlux(nDimensions,nSpecies)

     end subroutine computeFirstPartialViscousJacobian
     
  end interface

  interface

     pure subroutine computeSecondPartialViscousJacobian(nDimensions, nSpecies,              &
          velocity, dynamicViscosity, secondCoefficientOfViscosity,                          &
          thermalDiffusivity, jacobian, metricsAlongFirstDir,                                &
          metricsAlongSecondDir, secondPartialViscousJacobian)

       integer, intent(in) :: nDimensions, nSpecies
       SCALAR_TYPE, intent(in) :: velocity(nDimensions), dynamicViscosity,                   &
            secondCoefficientOfViscosity, thermalDiffusivity, jacobian,                      &
            metricsAlongFirstDir(nDimensions), metricsAlongSecondDir(nDimensions)
       SCALAR_TYPE, intent(out) :: secondPartialViscousJacobian(nDimensions+1,nDimensions+1)

     end subroutine computeSecondPartialViscousJacobian

  end interface

end module CNSHelper
