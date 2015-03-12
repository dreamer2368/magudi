#include "config.h"

module CNSHelper

  implicit none

  interface

     subroutine computeDependentVariables(nDimensions, conservedVariables,                   &
          ratioOfSpecificHeats, specificVolume,                                              &
          velocity, pressure, temperature)

       !> Computes the requested dependent variable(s) including specific volume, velocity,
       !> pressure and temperature from the conserved state variables.

       integer, intent(in) :: nDimensions
       SCALAR_TYPE, intent(in) :: conservedVariables(:,:)

       real(SCALAR_KIND), intent(in), optional :: ratioOfSpecificHeats
       SCALAR_TYPE, intent(out), optional :: specificVolume(:),                              &
            velocity(:,:), pressure(:), temperature(:)

     end subroutine computeDependentVariables

  end interface

  interface

     subroutine computeTransportVariables(temperature, powerLawExponent,                     &
          ratioOfSpecificHeats, reynoldsNumber, prandtlNumber,                               &
          dynamicViscosity, secondCoefficientOfViscosity,                                    &
          thermalDiffusivity)

       !> Computes the requested transport coefficient(s) including the dynamic viscosity,
       !> second coefficient of viscosity and the thermal conductivity assuming a power law
       !> dependence on temperature with exponent `powerLawExponent`.

       SCALAR_TYPE, intent(in) :: temperature(:)
       real(SCALAR_KIND), intent(in) :: powerLawExponent, ratioOfSpecificHeats,              &
            reynoldsNumber, prandtlNumber

       SCALAR_TYPE, intent(out), optional :: dynamicViscosity(:),                            &
            secondCoefficientOfViscosity(:),                                                 &
            thermalDiffusivity(:)

     end subroutine computeTransportVariables

  end interface

  interface

     subroutine computeStressTensor(nDimensions, velocityGradient, dynamicViscosity,         &
          secondCoefficientOfViscosity, stressTensor)

       integer, intent(in) :: nDimensions
       SCALAR_TYPE, intent(inout) :: velocityGradient(:,:)
       SCALAR_TYPE, intent(in) :: dynamicViscosity(:), secondCoefficientOfViscosity(:)

       SCALAR_TYPE, intent(out), optional :: stressTensor(:,:)

     end subroutine computeStressTensor

  end interface

  interface

     subroutine computeVorticityMagnitudeAndDilatation(nDimensions,                          &
          velocityGradient, vorticityMagnitude, dilatation)

       !> Computes the vorticity magnitude and dilatation from the velocity gradient.

       integer, intent(in) :: nDimensions
       SCALAR_TYPE, intent(in) :: velocityGradient(:,:)
       SCALAR_TYPE, intent(out), optional :: vorticityMagnitude(:), dilatation(:)

     end subroutine computeVorticityMagnitudeAndDilatation

  end interface

  interface

     subroutine computeCartesianInvsicidFluxes(nDimensions, conservedVariables,              &
          velocity, pressure, inviscidFluxes)

       !> Computes the Cartesian form of the inviscid fluxes. `velocity` and `pressure` must
       !> be consistent with `conservedVariables`, otherwise the result is unpredictable.

       integer, intent(in) :: nDimensions
       SCALAR_TYPE, intent(in) :: conservedVariables(:,:), velocity(:,:), pressure(:)
       SCALAR_TYPE, intent(out) :: inviscidFluxes(:,:,:)

     end subroutine computeCartesianInvsicidFluxes

  end interface

  interface

     subroutine computeCartesianViscousFluxes(nDimensions, velocity,                         &
          stressTensor, heatFlux, viscousFluxes)

       !> Computes the Cartesian form of the viscous fluxes.

       integer, intent(in) :: nDimensions
       SCALAR_TYPE, intent(in) :: velocity(:,:), stressTensor(:,:), heatFlux(:,:)
       SCALAR_TYPE, intent(out) :: viscousFluxes(:,:,:)

     end subroutine computeCartesianViscousFluxes

  end interface

  interface

     subroutine transformFluxes(nDimensions, fluxes, metrics,                                &
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

     function computeCfl(nDimensions, iblank, jacobian, metrics,                             &
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

     function computeTimeStepSize(nDimensions, iblank, jacobian,                             &
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

     subroutine computeJacobianOfInviscidFlux1D(conservedVariables, metrics,                 &
          ratioOfSpecificHeats, jacobianOfInviscidFlux, deltaConservedVariables,             &
          specificVolume, velocity, temperature, deltaJacobianOfInviscidFlux)

       SCALAR_TYPE, intent(in) :: conservedVariables(3), metrics(1), ratioOfSpecificHeats
       SCALAR_TYPE, intent(out) :: jacobianOfInviscidFlux(3,3)

       SCALAR_TYPE, intent(in), optional :: deltaConservedVariables(3,3), specificVolume,    &
            velocity(1), temperature
       SCALAR_TYPE, intent(out), optional :: deltaJacobianOfInviscidFlux(3,3,3)

     end subroutine computeJacobianOfInviscidFlux1D

  end interface

  interface

     subroutine computeJacobianOfInviscidFlux2D(conservedVariables, metrics,                 &
          ratioOfSpecificHeats, jacobianOfInviscidFlux, deltaConservedVariables,             &
          specificVolume, velocity, temperature, deltaJacobianOfInviscidFlux)

       SCALAR_TYPE, intent(in) :: conservedVariables(4), metrics(2), ratioOfSpecificHeats
       SCALAR_TYPE, intent(out) :: jacobianOfInviscidFlux(4,4)

       SCALAR_TYPE, intent(in), optional :: deltaConservedVariables(4,4), specificVolume,    &
            velocity(2), temperature
       SCALAR_TYPE, intent(out), optional :: deltaJacobianOfInviscidFlux(4,4,4)

     end subroutine computeJacobianOfInviscidFlux2D

  end interface

  interface

     subroutine computeJacobianOfInviscidFlux3D(conservedVariables, metrics,                 &
          ratioOfSpecificHeats, jacobianOfInviscidFlux, deltaConservedVariables,             &
          specificVolume, velocity, temperature, deltaJacobianOfInviscidFlux)

       SCALAR_TYPE, intent(in) :: conservedVariables(5), metrics(3), ratioOfSpecificHeats
       SCALAR_TYPE, intent(out) :: jacobianOfInviscidFlux(5,5)

       SCALAR_TYPE, intent(in), optional :: deltaConservedVariables(5,5), specificVolume,    &
            velocity(3), temperature
       SCALAR_TYPE, intent(out), optional :: deltaJacobianOfInviscidFlux(5,5,5)

     end subroutine computeJacobianOfInviscidFlux3D

  end interface

  interface

     subroutine computeIncomingJacobianOfInviscidFlux1D(conservedVariables, metrics,         &
          ratioOfSpecificHeats, incomingDirection, incomingJacobianOfInviscidFlux,           &
          deltaConservedVariables, specificVolume, velocity, temperature,                    &
          deltaIncomingJacobianOfInviscidFlux)

       SCALAR_TYPE, intent(in) :: conservedVariables(3), metrics(1), ratioOfSpecificHeats
       integer, intent(in) :: incomingDirection
       SCALAR_TYPE, intent(out) :: incomingJacobianOfInviscidFlux(3,3)
       SCALAR_TYPE, intent(in), optional :: deltaConservedVariables(3,3),                    &
            specificVolume, velocity(1), temperature
       SCALAR_TYPE, intent(out), optional :: deltaIncomingJacobianOfInviscidFlux(3,3,3)

     end subroutine computeIncomingJacobianOfInviscidFlux1D

  end interface

  interface

     subroutine computeIncomingJacobianOfInviscidFlux2D(conservedVariables, metrics,         &
          ratioOfSpecificHeats, incomingDirection, incomingJacobianOfInviscidFlux,           &
          deltaConservedVariables, specificVolume, velocity, temperature,                    &
          deltaIncomingJacobianOfInviscidFlux)

       SCALAR_TYPE, intent(in) :: conservedVariables(4), metrics(2), ratioOfSpecificHeats
       integer, intent(in) :: incomingDirection
       SCALAR_TYPE, intent(out) :: incomingJacobianOfInviscidFlux(4,4)
       SCALAR_TYPE, intent(in), optional :: deltaConservedVariables(4,4),                    &
            specificVolume, velocity(2), temperature
       SCALAR_TYPE, intent(out), optional :: deltaIncomingJacobianOfInviscidFlux(4,4,4)

     end subroutine computeIncomingJacobianOfInviscidFlux2D

  end interface

  interface

     subroutine computeIncomingJacobianOfInviscidFlux3D(conservedVariables, metrics,         &
          ratioOfSpecificHeats, incomingDirection, incomingJacobianOfInviscidFlux,           &
          deltaConservedVariables, specificVolume, velocity, temperature,                    &
          deltaIncomingJacobianOfInviscidFlux)

       SCALAR_TYPE, intent(in) :: conservedVariables(5), metrics(3), ratioOfSpecificHeats
       integer, intent(in) :: incomingDirection
       SCALAR_TYPE, intent(out) :: incomingJacobianOfInviscidFlux(5,5)
       SCALAR_TYPE, intent(in), optional :: deltaConservedVariables(5,5),                    &
            specificVolume, velocity(3), temperature
       SCALAR_TYPE, intent(out), optional :: deltaIncomingJacobianOfInviscidFlux(5,5,5)

     end subroutine computeIncomingJacobianOfInviscidFlux3D

  end interface

  interface

     subroutine computeFirstPartialViscousJacobian1D(conservedVariables, metrics,            &
          stressTensor, heatFlux, powerLawExponent, ratioOfSpecificHeats,                    &
          firstPartialViscousJacobian, specificVolume, velocity, temperature)

       SCALAR_TYPE, intent(in) :: conservedVariables(3), metrics(1), stressTensor(1),        &
            heatFlux(1), powerLawExponent, ratioOfSpecificHeats
       SCALAR_TYPE, intent(out) :: firstPartialViscousJacobian(3,3)
       SCALAR_TYPE, intent(in), optional :: specificVolume, velocity(1), temperature

     end subroutine computeFirstPartialViscousJacobian1D

  end interface

  interface

     subroutine computeFirstPartialViscousJacobian2D(conservedVariables, metrics,            &
          stressTensor, heatFlux, powerLawExponent, ratioOfSpecificHeats,                    &
          firstPartialViscousJacobian, specificVolume, velocity, temperature)

       SCALAR_TYPE, intent(in) :: conservedVariables(4), metrics(2), stressTensor(4),        &
            heatFlux(2), powerLawExponent, ratioOfSpecificHeats
       SCALAR_TYPE, intent(out) :: firstPartialViscousJacobian(4,4)
       SCALAR_TYPE, intent(in), optional :: specificVolume, velocity(2), temperature

     end subroutine computeFirstPartialViscousJacobian2D

  end interface

  interface

     subroutine computeFirstPartialViscousJacobian3D(conservedVariables, metrics,            &
          stressTensor, heatFlux, powerLawExponent, ratioOfSpecificHeats,                    &
          firstPartialViscousJacobian, specificVolume, velocity, temperature)

       SCALAR_TYPE, intent(in) :: conservedVariables(5), metrics(3), stressTensor(9),        &
            heatFlux(3), powerLawExponent, ratioOfSpecificHeats
       SCALAR_TYPE, intent(out) :: firstPartialViscousJacobian(5,5)
       SCALAR_TYPE, intent(in), optional :: specificVolume, velocity(3), temperature

     end subroutine computeFirstPartialViscousJacobian3D

  end interface

  interface

     subroutine computeSecondPartialViscousJacobian1D(velocity, dynamicViscosity,            &
          secondCoefficientOfViscosity, thermalDiffusivity, jacobian, metrics,               &
          secondPartialViscousJacobian)

       SCALAR_TYPE, intent(in) :: velocity(1), dynamicViscosity,                             &
            secondCoefficientOfViscosity, thermalDiffusivity, jacobian, metrics
       SCALAR_TYPE, intent(out) :: secondPartialViscousJacobian(2,2)

     end subroutine computeSecondPartialViscousJacobian1D

  end interface

  interface

     subroutine computeSecondPartialViscousJacobian2D(velocity, dynamicViscosity,            &
          secondCoefficientOfViscosity, thermalDiffusivity, jacobian, metricsAlongFirstDir,  &
          metricsAlongSecondDir, secondPartialViscousJacobian)

       SCALAR_TYPE, intent(in) :: velocity(2), dynamicViscosity,                             &
            secondCoefficientOfViscosity, thermalDiffusivity, jacobian,                      &
            metricsAlongFirstDir(2), metricsAlongSecondDir(2)
       SCALAR_TYPE, intent(out) :: secondPartialViscousJacobian(3,3)

     end subroutine computeSecondPartialViscousJacobian2D

  end interface

  interface

     subroutine computeSecondPartialViscousJacobian3D(velocity, dynamicViscosity,            &
          secondCoefficientOfViscosity, thermalDiffusivity, jacobian, metricsAlongFirstDir,  &
          metricsAlongSecondDir, secondPartialViscousJacobian)

       SCALAR_TYPE, intent(in) :: velocity(3), dynamicViscosity,                             &
            secondCoefficientOfViscosity, thermalDiffusivity, jacobian,                      &
            metricsAlongFirstDir(3), metricsAlongSecondDir(3)
       SCALAR_TYPE, intent(out) :: secondPartialViscousJacobian(4,4)

     end subroutine computeSecondPartialViscousJacobian3D

  end interface

end module CNSHelper
