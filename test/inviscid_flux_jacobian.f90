#include "config.h"

program inviscid_flux_jacobian

  use CNSHelper
  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random

  implicit none

  integer, parameter :: wp = SCALAR_KIND, n = 1024, maxSpecies = 6
  real(wp), parameter :: ratioOfSpecificHeats = 1.4_wp
  logical :: success
  integer :: i, j, k, nDimensions, nSpecies
  SCALAR_TYPE, allocatable :: conservedVariables1(:,:), deltaPrimitiveVariables(:,:),        &
       deltaConservedVariables(:,:),  conservedVariables2(:,:), metrics(:,:),                &
       specificVolume(:), velocity(:,:), pressure(:), temperature(:), massFraction(:,:),     &
       fluxes1(:,:,:), fluxes2(:,:,:), fluxes3(:,:,:), localConservedVariables(:),           &
       localVelocity(:), localMassFraction(:), localMetricsAlongDirection(:),                &
       jacobianOfInviscidFlux(:,:)
  logical :: isDomainCurvilinear
  real(wp), allocatable :: stepSizes(:), errorHistory(:), convergenceHistory(:)

  interface

     real(SCALAR_KIND) function meanTrimmed(a)
       real(SCALAR_KIND), intent(inout) :: a(:)
     end function meanTrimmed

     subroutine sort(a)
       real(SCALAR_KIND), intent(inout) :: a(:)
     end subroutine sort

  end interface

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  nSpecies = random(0, maxSpecies)

  print *
  print *, 'Number of species:',nSpecies
  print *

  do nDimensions = 1, 3

     allocate(conservedVariables1(n, nDimensions + 2 + nSpecies))
     allocate(deltaPrimitiveVariables(n, nDimensions + 2 + nSpecies))
     allocate(deltaConservedVariables(n, nDimensions + 2 + nSpecies))
     allocate(conservedVariables2(n, nDimensions + 2 + nSpecies))
     allocate(metrics(n, nDimensions ** 2))
     allocate(specificVolume(n))
     allocate(velocity(n, nDimensions))
     allocate(pressure(n))
     allocate(temperature(n))
     allocate(massFraction(n, nSpecies))
     allocate(fluxes1(n, nDimensions + 2 + nSpecies, nDimensions))
     allocate(fluxes2(n, nDimensions + 2 + nSpecies, nDimensions))
     allocate(fluxes3(n, nDimensions + 2 + nSpecies, nDimensions))
     allocate(localConservedVariables(nDimensions + 2 + nSpecies))
     allocate(localVelocity(nDimensions))
     allocate(localMassFraction(nSpecies))
     allocate(localMetricsAlongDirection(nDimensions))
     allocate(jacobianOfInviscidFlux(nDimensions + 2 + nSpecies, nDimensions + 2 + nSpecies))

     isDomainCurvilinear = (random(0, 1) == 0)

     do i = 1, n

        conservedVariables1(i,1) = random(0.01_wp, 10.0_wp)
        do j = 1, nDimensions
           conservedVariables1(i,j+1) = conservedVariables1(i,1) * random(-10.0_wp, 10.0_wp)
        end do
        conservedVariables1(i, nDimensions + 2) = conservedVariables1(i,1) *                 &
             random(0.01_wp, 10.0_wp) / ratioOfSpecificHeats +                               &
             0.5_wp / conservedVariables1(i,1) *                                             &
             sum(conservedVariables1(i,2:nDimensions+1) ** 2)
        do k = 1, nSpecies
           conservedVariables1(i, nDimensions + 2 + k) = conservedVariables1(i,1) *          &
                random(0.0_wp, 1.0_wp)
        end do

        assert(conservedVariables1(i,1) > 0.0_wp)

        if (isDomainCurvilinear) then
           do j = 1, nDimensions ** 2
              metrics(i,j) = random(-1.0_wp, 1.0_wp)
           end do
        else
           metrics(i,:) = 0.0_wp
           do j = 1, nDimensions
              metrics(i,j+nDimensions*(j-1)) = random(-1.0_wp, 1.0_wp)
           end do
        end if

        do j = 1, nDimensions + 2 + nSpecies
           deltaPrimitiveVariables(i,j) = random(-1.0_wp, 1.0_wp)
        end do

     end do

     call computeDependentVariables(nDimensions, nSpecies, conservedVariables1,              &
          ratioOfSpecificHeats = ratioOfSpecificHeats, specificVolume = specificVolume,      &
          velocity = velocity, pressure = pressure, temperature = temperature,               &
          massFraction = massFraction)
     assert(all(specificVolume > 0.0_wp))
     assert(all(temperature > 0.0_wp))
     assert(all(massFraction > 0.0_wp))
     call computeCartesianInvsicidFluxes(nDimensions, nSpecies, conservedVariables1,         &
          velocity, pressure, fluxes1)
     call transformFluxes(nDimensions, fluxes1, metrics, fluxes2, isDomainCurvilinear)

     allocate(stepSizes(100))
     stepSizes(1) = 0.01_wp
     do i = 2, size(stepSizes)
        stepSizes(i) = stepSizes(i-1) / 1.2_wp
     end do

     allocate(errorHistory(size(stepSizes)))
     allocate(convergenceHistory(size(stepSizes) - 1))

     do i = 1, size(stepSizes)

        deltaConservedVariables(:,1) = deltaPrimitiveVariables(:,1)
        do j = 1, nDimensions
           deltaConservedVariables(:,j+1) = conservedVariables1(:,j+1) /                     &
                conservedVariables1(:,1) * deltaPrimitiveVariables(:,1) +                    &
                conservedVariables1(:,1) * deltaPrimitiveVariables(:,j+1)
        end do
        deltaConservedVariables(:,nDimensions+2) = conservedVariables1(:,nDimensions+2) /    &
             conservedVariables1(:,1) * deltaPrimitiveVariables(:,1) +                       &
             sum(conservedVariables1(:,2:nDimensions+1) *                                    &
             deltaPrimitiveVariables(:,2:nDimensions+1), dim = 2) +                          &
             conservedVariables1(:,1) / ratioOfSpecificHeats *                               &
             deltaPrimitiveVariables(:,nDimensions+2)
        do k = 1, nSpecies
           deltaConservedVariables(:,nDimensions+2+k) =                                      &
                conservedVariables1(:,nDimensions+2+k) / conservedVariables1(:,1) *          &
                deltaPrimitiveVariables(:,1) + conservedVariables1(:,1) *                    &
                deltaPrimitiveVariables(:,nDimensions+2+k)
        end do

        conservedVariables2 = conservedVariables1 + stepSizes(i) * deltaConservedVariables
        assert(all(conservedVariables2(:,1) > 0.0_wp))

        call computeDependentVariables(nDimensions, nSpecies, conservedVariables2,           &
             ratioOfSpecificHeats =  ratioOfSpecificHeats, specificVolume = specificVolume,  &
             velocity = velocity, pressure = pressure, temperature = temperature,            &
             massFraction = massFraction)

        assert(all(specificVolume > 0.0_wp))
        assert(all(temperature > 0.0_wp))
        call computeCartesianInvsicidFluxes(nDimensions, nSpecies, conservedVariables2,      &
             velocity, pressure, fluxes1)
        call transformFluxes(nDimensions, fluxes1, metrics, fluxes3, isDomainCurvilinear)

        errorHistory(i) = 0.0_wp

        do j = 1, n

           localConservedVariables = conservedVariables1(j,:)
           localVelocity = velocity(j,:)
           localMassFraction = massFraction(j,:)

           do k = 1, nDimensions

              localMetricsAlongDirection = metrics(j,1+nDimensions*(k-1):nDimensions*k)

              call computeJacobianOfInviscidFlux(nDimensions, nSpecies,                      &
                   localConservedVariables, localMetricsAlongDirection,                      &
                   ratioOfSpecificHeats, jacobianOfInviscidFlux,                             &
                   specificVolume = specificVolume(j), velocity = localVelocity,             &
                   temperature = temperature(j), massFraction = localMassFraction)

              errorHistory(i) = max(errorHistory(i),                                         &
                   maxval(abs(fluxes3(j,:,k) - fluxes2(j,:,k) - stepSizes(i) *               &
                   matmul(jacobianOfInviscidFlux, deltaConservedVariables(j,:)))))

           end do

        end do

        if (i > 1) then
           convergenceHistory(i-1) = log(errorHistory(i) / errorHistory(i-1)) /              &
                log(stepSizes(i) / stepSizes(i-1))
           if (convergenceHistory(i-1) < 0.0_wp) exit
        end if

     end do

     if (i > 2) then
        call sort(convergenceHistory(:i-2))
        success = success .and. (nint(meanTrimmed(convergenceHistory(:i-2))) == 2)
     else
        success = .false.
     end if

     SAFE_DEALLOCATE(convergenceHistory)
     SAFE_DEALLOCATE(errorHistory)
     SAFE_DEALLOCATE(stepSizes)

     SAFE_DEALLOCATE(jacobianOfInviscidFlux)
     SAFE_DEALLOCATE(localMetricsAlongDirection)
     SAFE_DEALLOCATE(localMassFraction)
     SAFE_DEALLOCATE(localVelocity)
     SAFE_DEALLOCATE(localConservedVariables)
     SAFE_DEALLOCATE(fluxes3)
     SAFE_DEALLOCATE(fluxes2)
     SAFE_DEALLOCATE(fluxes1)
     SAFE_DEALLOCATE(massFraction)
     SAFE_DEALLOCATE(temperature)
     SAFE_DEALLOCATE(pressure)
     SAFE_DEALLOCATE(velocity)
     SAFE_DEALLOCATE(specificVolume)
     SAFE_DEALLOCATE(metrics)
     SAFE_DEALLOCATE(conservedVariables2)
     SAFE_DEALLOCATE(deltaConservedVariables)
     SAFE_DEALLOCATE(deltaPrimitiveVariables)
     SAFE_DEALLOCATE(conservedVariables1)

     if (.not. success) exit

  end do

  call cleanupErrorHandler()

  if (.not. success) stop -1
  stop 0

end program inviscid_flux_jacobian

subroutine sort(a)

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(inout) :: a(:)

  ! <<< Local variables >>>
  integer :: i, j
  real(SCALAR_KIND) :: temp

  do i = 2, size(a)
     j = i - 1
     temp = a(i)
     do while (a(j) > temp)
        a(j+1) = a(j)
        j = j - 1
        if (j < 1) exit
     end do
     a(j+1) = temp
  end do

end subroutine sort

real(SCALAR_KIND) function median(a)

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(in) :: a(:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: n

  n = size(a)
  if (mod(n, 2) == 1) then
     median = a((n + 1) / 2)
  else
     median = 0.5_wp * (a(n / 2) + a(n / 2 + 1))
  end if

end function median

real(SCALAR_KIND) function meanTrimmed(a)

  ! <<< Arguments >>>
  real(SCALAR_KIND), intent(inout) :: a(:)

  ! <<< Scalar variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, n
  real(wp) :: firstQuartile, thirdQuartile

  interface
     real(SCALAR_KIND) function median(a)
       real(SCALAR_KIND), intent(in) :: a(:)
     end function median
  end interface

  n = size(a)

  if (mod(n, 2) == 0) then
     firstQuartile = median(a(1:n/2))
     thirdQuartile = median(a(n/2+1:n))
  else
     firstQuartile = median(a(1:(n-1)/2))
     thirdQuartile = median(a((n+1)/2+1:n))
  end if

  meanTrimmed = 0.0_wp
  n = 0

  do i = 1, size(a)
     if (a(i) >= firstQuartile .and. a(i) <= thirdQuartile) then
        meanTrimmed = meanTrimmed + a(i)
        n = n + 1
     end if
  end do

  if (n == 0) return
  meanTrimmed = meanTrimmed / real(n, wp)

end function meanTrimmed
