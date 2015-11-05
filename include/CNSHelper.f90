#include "config.h"

module CNSHelper

  implicit none
  public

  interface

     pure subroutine computeDependentVariables(nDimensions, nSpecies, conservedVariables,    &
          equationOfState, ratioOfSpecificHeats, molecularWeightInverse, specificVolume,     &
          velocity, pressure, temperature, massFraction)

       !> Computes the requested dependent variable(s) including specific volume, velocity,
       !> pressure, temperature, and mass fraction from the conserved state variables.

       integer, intent(in) :: nDimensions, nSpecies
       integer, intent(in), optional :: equationOfState
       SCALAR_TYPE, intent(in) :: conservedVariables(:,:)

       real(SCALAR_KIND), intent(in), optional :: ratioOfSpecificHeats,                      &
            molecularWeightInverse(:)
       SCALAR_TYPE, intent(out), optional :: specificVolume(:),                              &
            velocity(:,:), pressure(:), temperature(:), massFraction(:,:)

     end subroutine computeDependentVariables

  end interface

  interface

     pure subroutine computeTransportVariables(nSpecies, temperature, powerLawExponent,      &
          bulkViscosityRatio, ratioOfSpecificHeats, reynoldsNumberInverse,                   &
          prandtlNumberInverse, schmidtNumberInverse, dynamicViscosity,                      &
          secondCoefficientOfViscosity, thermalDiffusivity, massDiffusivity)

       !> Computes the requested transport coefficient(s) including the dynamic viscosity,
       !> second coefficient of viscosity and the thermal conductivity assuming a power law
       !> dependence on temperature with exponent `powerLawExponent`.

       integer, intent(in) :: nSpecies
       SCALAR_TYPE, intent(in) :: temperature(:)
       real(SCALAR_KIND), intent(in) :: powerLawExponent,                                    &
            ratioOfSpecificHeats, reynoldsNumberInverse

       real(SCALAR_KIND), intent(in), optional :: bulkViscosityRatio, prandtlNumberInverse,  &
            schmidtNumberInverse(:)
       SCALAR_TYPE, intent(out), optional :: dynamicViscosity(:),                            &
            secondCoefficientOfViscosity(:), thermalDiffusivity(:), massDiffusivity(:,:)

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

     pure subroutine computeStressTensor(nDimensions, velocityGradient, dynamicViscosity,   &
          secondCoefficientOfViscosity, stressTensor)

       integer, intent(in) :: nDimensions
       SCALAR_TYPE, intent(inout) :: velocityGradient(:,:)
       SCALAR_TYPE, intent(in) :: dynamicViscosity(:), secondCoefficientOfViscosity(:)

       SCALAR_TYPE, intent(out), optional :: stressTensor(:,:)

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
       conservedVariables, velocity, pressure, inviscidFluxes)

       !> Computes the Cartesian form of the inviscid fluxes. `velocity`, `pressure`, and
       !> `mass fraction` must be consistent with `conservedVariables`, otherwise the result
       !> is unpredictable.

       integer, intent(in) :: nDimensions, nSpecies
       SCALAR_TYPE, intent(in) :: conservedVariables(:,:), velocity(:,:), pressure(:)
       SCALAR_TYPE, intent(out) :: inviscidFluxes(:,:,:)

     end subroutine computeCartesianInvsicidFluxes

  end interface

  interface

     pure subroutine computeCartesianViscousFluxes(nDimensions, nSpecies, velocity,          &
          stressTensor, heatFlux, viscousFluxes, massFraction, speciesFlux, enthalpyFlux)

       !> Computes the Cartesian form of the viscous fluxes.

       integer, intent(in) :: nDimensions, nSpecies
       SCALAR_TYPE, intent(in) :: velocity(:,:), stressTensor(:,:), heatFlux(:,:)
       SCALAR_TYPE, intent(in), optional :: massFraction(:,:), speciesFlux(:,:,:),           &
            enthalpyFlux(:,:)
       SCALAR_TYPE, intent(out) :: viscousFluxes(:,:,:)

     end subroutine computeCartesianViscousFluxes

  end interface

  interface

     pure subroutine computeSpectralRadius(nDimensions, ratioOfSpecificHeats,                &
          specificVolume, velocity, pressure, metrics, spectralRadius, isDomainCurvilinear)

       !> Compute the spectral radii along all the directions.

       integer, intent(in) :: nDimensions
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(in) :: specificVolume(:), velocity(:,:), pressure(:), metrics(:,:)
       SCALAR_TYPE, intent(out) :: spectralRadius(:,:)
       logical, intent(in), optional :: isDomainCurvilinear

     end subroutine computeSpectralRadius

  end interface

  interface

     pure subroutine transformFluxes(nDimensions, fluxes, metrics,                           &
          transformedFluxes, isDomainCurvilinear)

       !> Transforms fluxes from Cartesian form to contravariant form. If the domain is not
       !> curvilinear as specified by `isDomainCurvilinear`, the arguments `fluxes` and
       !> `transformedFluxes` may reference the same array.

       integer, intent(in) :: nDimensions
       SCALAR_TYPE, intent(in) :: fluxes(:,:,:), metrics(:,:)
       SCALAR_TYPE, intent(out) :: transformedFluxes(:,:,:)
       logical, intent(in), optional :: isDomainCurvilinear

     end subroutine transformFluxes

  end interface

  interface

     pure function computeCfl(nDimensions, iblank, jacobian, metrics,                        &
          velocity, pressure, specificVolume, timeStepSize, ratioOfSpecificHeats,            &
          dynamicViscosity, thermalDiffusivity) result(cfl)

       !> Computes the CFL number.

       integer, intent(in) :: nDimensions, iblank(:)
       SCALAR_TYPE, intent(in) :: jacobian(:), metrics(:,:), velocity(:,:), pressure(:),     &
            specificVolume(:)
       real(SCALAR_KIND), intent(in) :: timeStepSize, ratioOfSpecificHeats
       SCALAR_TYPE, intent(in), optional :: dynamicViscosity(:), thermalDiffusivity(:)

       real(SCALAR_KIND) :: cfl

     end function computeCfl

  end interface

  interface

     pure function computeTimeStepSize(nDimensions, iblank, jacobian,                        &
          metrics, velocity, pressure, specificVolume, cfl,                                  &
          ratioOfSpecificHeats, dynamicViscosity,                                            &
          thermalDiffusivity) result(timeStepSize)

       !> Computes the time step size that leads to a specified CFL number.

       integer, intent(in) :: nDimensions, iblank(:)
       SCALAR_TYPE, intent(in) :: jacobian(:), metrics(:,:),                                 &
            velocity(:,:), pressure(:), specificVolume(:)
       real(SCALAR_KIND), intent(in) :: cfl, ratioOfSpecificHeats
       SCALAR_TYPE, intent(in), optional :: dynamicViscosity(:), thermalDiffusivity(:)

       real(SCALAR_KIND) :: timeStepSize

     end function computeTimeStepSize

  end interface

  interface

     pure subroutine computeJacobianOfInviscidFlux(nDimensions, nSpecies, conservedVariables,&
          metrics, ratioOfSpecificHeats, jacobianOfInviscidFlux, deltaConservedVariables,    &
          specificVolume, velocity, pressure, massFraction, deltaJacobianOfInviscidFlux)

       integer, intent(in) :: nDimensions, nSpecies
       SCALAR_TYPE, intent(in) :: conservedVariables(:), metrics(:)
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(out) :: jacobianOfInviscidFlux(:,:)

       SCALAR_TYPE, intent(in), optional :: deltaConservedVariables, specificVolume,         &
            velocity(:), pressure, massFraction(:)
       SCALAR_TYPE, intent(out), optional :: deltaJacobianOfInviscidFlux(:,:,:)

     end subroutine computeJacobianOfInviscidFlux

  end interface

  interface

     pure subroutine computeIncomingJacobianOfInviscidFlux(nDimensions, nSpecies,            &
          conservedVariables, metrics, ratioOfSpecificHeats, incomingDirection,              &
          incomingJacobianOfInviscidFlux, deltaIncomingJacobianOfInviscidFlux,               &
          deltaConservedVariables, specificVolume, velocity, pressure, massFraction)

       integer, intent(in) :: nDimensions, nSpecies, incomingDirection
       SCALAR_TYPE, intent(in) :: conservedVariables(:), metrics(:)
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(out) :: incomingJacobianOfInviscidFlux(:,:)

       SCALAR_TYPE, intent(out), optional :: deltaIncomingJacobianOfInviscidFlux(:,:,:)
       SCALAR_TYPE, intent(in), optional :: deltaConservedVariables(:,:), specificVolume,    &
            velocity(:), pressure, massFraction(:)

     end subroutine computeIncomingJacobianOfInviscidFlux

  end interface

  interface

     pure subroutine computeFirstPartialViscousJacobian(nDimensions, nSpecies,               &
          equationOfState, conservedVariables, metrics, stressTensor, heatFlux,              &
          enthalpyFlux, speciesFlux, powerLawExponent, ratioOfSpecificHeats,                 &
          firstPartialViscousJacobian, specificVolume, velocity, temperature, massFraction,  &
          molecularWeightInverse)

       integer, intent(in) :: nDimensions, nSpecies, equationOfState
       SCALAR_TYPE, intent(in) :: conservedVariables(:), metrics(:),                         &
            stressTensor(:), heatFlux(:)
       real(SCALAR_KIND), intent(in) :: powerLawExponent, ratioOfSpecificHeats
       SCALAR_TYPE, intent(out) :: firstPartialViscousJacobian(:,:)

       SCALAR_TYPE, intent(in), optional :: specificVolume, velocity(:), temperature,        &
            massFraction(:), enthalpyFlux(:), speciesFlux(:,:), molecularWeightInverse(:)

     end subroutine computeFirstPartialViscousJacobian
     
  end interface

  interface

     pure subroutine computeSecondPartialViscousJacobian(nDimensions, nSpecies,              &
          equationOfState, velocity, dynamicViscosity, secondCoefficientOfViscosity,         &
          thermalDiffusivity, massDiffusivity, jacobian, metricsAlongFirstDir,               &
          metricsAlongSecondDir, secondPartialViscousJacobian)

       integer, intent(in) :: nDimensions, nSpecies, equationOfState
       SCALAR_TYPE, intent(in) :: velocity(:), dynamicViscosity,                             &
            secondCoefficientOfViscosity, thermalDiffusivity, massDiffusivity(:),            &
            jacobian, metricsAlongFirstDir(:)
       SCALAR_TYPE, intent(in), optional :: metricsAlongSecondDir(:)
       SCALAR_TYPE, intent(out) :: secondPartialViscousJacobian(:,:)

     end subroutine computeSecondPartialViscousJacobian

  end interface

  interface

     pure subroutine computeJacobianOfSource(nDimensions, nSpecies, equationOfState,         &
          conservedVariables, ratioOfSpecificHeats, combustion, jacobianOfSource,            &
          specificVolume, velocity, temperature, massFraction)

       use Combustion_mod, only : t_Combustion

       integer, intent(in) :: nDimensions, nSpecies, equationOfState
       SCALAR_TYPE, intent(in) :: conservedVariables(:)
       real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(out) :: jacobianOfSource(:,:)
       type(t_Combustion), intent(in) :: combustion

       SCALAR_TYPE, intent(in), optional :: specificVolume, velocity(:), temperature,        &
            massFraction(:)

     end subroutine computeJacobianOfSource

  end interface

end module CNSHelper
