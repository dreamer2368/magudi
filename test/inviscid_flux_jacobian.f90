#include "config.h"

program inviscid_flux_jacobian

  use CNSHelper
  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random

  implicit none

  integer, parameter :: wp = SCALAR_KIND, n = 1024
  real(wp), parameter :: ratioOfSpecificHeats = 1.4_wp
  logical :: success
  integer :: i, j, k, nDimensions
  SCALAR_TYPE, allocatable :: conservedVariables1(:,:), deltaConservedVariables(:,:),        &
       conservedVariables2(:,:), metrics(:,:), specificVolume(:), velocity(:,:),             &
       pressure(:), temperature(:), fluxes1(:,:,:), fluxes2(:,:,:), fluxes3(:,:,:),          &
       localConservedVariables(:), localVelocity(:), localMetricsAlongDirection(:),          &
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

  do nDimensions = 1, 3

     allocate(conservedVariables1(n, nDimensions + 2))
     allocate(deltaConservedVariables(n, nDimensions + 2))
     allocate(conservedVariables2(n, nDimensions + 2))
     allocate(metrics(n, nDimensions ** 2))
     allocate(specificVolume(n))
     allocate(velocity(n, nDimensions))
     allocate(pressure(n))
     allocate(temperature(n))
     allocate(fluxes1(n, nDimensions + 2, nDimensions))
     allocate(fluxes2(n, nDimensions + 2, nDimensions))
     allocate(fluxes3(n, nDimensions + 2, nDimensions))
     allocate(localConservedVariables(nDimensions + 2))
     allocate(localVelocity(nDimensions))
     allocate(localMetricsAlongDirection(nDimensions))
     allocate(jacobianOfInviscidFlux(nDimensions + 2, nDimensions + 2))

     isDomainCurvilinear = (random(0, 2) == 0)

     do i = 1, n

        conservedVariables1(i,1) = random(0.01_wp, 10.0_wp)
        do j = 1, nDimensions
           conservedVariables1(i,j+1) = conservedVariables1(i,1) * random(-10.0_wp, 10.0_wp)
        end do
        conservedVariables1(i, nDimensions + 2) = conservedVariables1(i,1) *                 &
             (random(0.01_wp, 10.0_wp) / ratioOfSpecificHeats +                              &
             0.5_wp * sum(conservedVariables1(i,2:nDimensions+1) ** 2))

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

        deltaConservedVariables(i,1) = random(-1.0_wp, 1.0_wp)
        do j = 1, nDimensions
           deltaConservedVariables(i,j+1) =                                                  &
                deltaConservedVariables(i,1) * random(-1.0_wp, 1.0_wp)
        end do
        deltaConservedVariables(i, nDimensions + 2) = deltaConservedVariables(i,1) *         &
             (random(-1.0_wp, 1.0_wp) / ratioOfSpecificHeats +                               &
             0.5_wp * sum(deltaConservedVariables(i,2:nDimensions+1) ** 2))

     end do

     call computeDependentVariables(nDimensions, conservedVariables1,                        &
          ratioOfSpecificHeats, specificVolume = specificVolume, velocity = velocity,        &
          pressure = pressure, temperature = temperature)
     call computeCartesianInvsicidFluxes(nDimensions, conservedVariables1,                   &
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

        conservedVariables2 = conservedVariables1 + stepSizes(i) * deltaConservedVariables

        call computeDependentVariables(nDimensions, conservedVariables2,                     &
             ratioOfSpecificHeats, specificVolume = specificVolume, velocity = velocity,     &
             pressure = pressure, temperature = temperature)
        call computeCartesianInvsicidFluxes(nDimensions, conservedVariables2,                &
             velocity, pressure, fluxes1)
        call transformFluxes(nDimensions, fluxes1, metrics, fluxes3, isDomainCurvilinear)

        errorHistory(i) = 0.0_wp

        do j = 1, n

           localConservedVariables = conservedVariables1(j,:)
           localVelocity = velocity(j,:)

           do k = 1, nDimensions

              localMetricsAlongDirection = metrics(j,1+nDimensions*(k-1):nDimensions*k)

              select case (nDimensions)
              case (1)
                 call computeJacobianOfInviscidFlux1D(localConservedVariables,               &
                      localMetricsAlongDirection, ratioOfSpecificHeats,                      &
                      jacobianOfInviscidFlux, specificVolume = specificVolume(j),            &
                      velocity = localVelocity, temperature = temperature(j))
              case (2)
                 call computeJacobianOfInviscidFlux2D(localConservedVariables,               &
                      localMetricsAlongDirection, ratioOfSpecificHeats,                      &
                      jacobianOfInviscidFlux, specificVolume = specificVolume(j),            &
                      velocity = localVelocity, temperature = temperature(j))
              case (3)
                 call computeJacobianOfInviscidFlux3D(localConservedVariables,               &
                      localMetricsAlongDirection, ratioOfSpecificHeats,                      &
                      jacobianOfInviscidFlux, specificVolume = specificVolume(j),            &
                      velocity = localVelocity, temperature = temperature(j))
              end select

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
     SAFE_DEALLOCATE(localVelocity)
     SAFE_DEALLOCATE(localConservedVariables)
     SAFE_DEALLOCATE(fluxes3)
     SAFE_DEALLOCATE(fluxes2)
     SAFE_DEALLOCATE(fluxes1)
     SAFE_DEALLOCATE(temperature)
     SAFE_DEALLOCATE(pressure)
     SAFE_DEALLOCATE(velocity)
     SAFE_DEALLOCATE(specificVolume)
     SAFE_DEALLOCATE(metrics)
     SAFE_DEALLOCATE(conservedVariables2)
     SAFE_DEALLOCATE(deltaConservedVariables)
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
