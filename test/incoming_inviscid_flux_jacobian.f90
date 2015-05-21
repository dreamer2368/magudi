#include "config.h"

program incoming_inviscid_flux_jacobian

  use CNSHelper
  use ErrorHandler, only : initializeErrorHandler, cleanupErrorHandler
  use RandomNumber, only : initializeRandomNumberGenerator, random

  implicit none

  integer, parameter :: wp = SCALAR_KIND, n = 2 ** 10, maxSpecies = 6
  real(wp), parameter :: ratioOfSpecificHeats = 1.4_wp
  logical :: success
  integer :: i, j, k, nDimensions, nSpecies
  SCALAR_TYPE, allocatable :: conservedVariables(:,:), metrics(:,:), specificVolume(:),      &
       velocity(:,:), temperature(:), massFraction(:,:), localConservedVariables(:),         &
       localVelocity(:), localMassFraction(:), localMetricsAlongDirection(:),                &
       jacobianOfInviscidFlux(:,:), incomingJacobianOfInviscidFlux1(:,:),                    &
       incomingJacobianOfInviscidFlux2(:,:)
  logical :: isDomainCurvilinear

  interface

     subroutine lapackComputeIncomingJacobianOfInviscidFlux(jacobian,                        &
          incomingDirection, incomingJacobian)

       SCALAR_TYPE, intent(in) :: jacobian(:,:)
       integer, intent(in) :: incomingDirection
       SCALAR_TYPE, intent(out) :: incomingJacobian(:,:)

     end subroutine lapackComputeIncomingJacobianOfInviscidFlux

  end interface

  success = .true.

  call initializeErrorHandler()
  call initializeRandomNumberGenerator()

  nSpecies = random(0, maxSpecies)

  print *
  print *, 'Number of species:',nSpecies
  print *

  do nDimensions = 3, 3

     allocate(conservedVariables(n, nDimensions + 2 + nSpecies))
     allocate(metrics(n, nDimensions ** 2))
     allocate(specificVolume(n))
     allocate(velocity(n, nDimensions))
     allocate(temperature(n))
     allocate(massFraction(n, nSpecies))
     allocate(localConservedVariables(nDimensions + 2 + nSpecies))
     allocate(localVelocity(nDimensions))
     allocate(localMassFraction(nSpecies))
     allocate(localMetricsAlongDirection(nDimensions))
     allocate(jacobianOfInviscidFlux                                                         &
          (nDimensions + 2 + nSpecies, nDimensions + 2 + nSpecies))
     allocate(incomingJacobianOfInviscidFlux1                                                &
          (nDimensions + 2 + nSpecies, nDimensions + 2 + nSpecies))
     allocate(incomingJacobianOfInviscidFlux2                                                &
          (nDimensions + 2 + nSpecies, nDimensions + 2 + nSpecies))

     isDomainCurvilinear = (random(0, 1) == 0)

     do i = 1, n

        conservedVariables(i,1) = random(0.1_wp, 4.0_wp)
        do j = 1, nDimensions
           conservedVariables(i,j+1) = conservedVariables(i,1) * random(-4.0_wp, 4.0_wp)
        end do
        conservedVariables(i, nDimensions + 2) = conservedVariables(i,1) *                   &
             random(0.1_wp, 4.0_wp) / ratioOfSpecificHeats +                                 &
             0.5_wp / conservedVariables(i,1) *                                              &
             sum(conservedVariables(i,2:nDimensions+1) ** 2)
        do k = 1, nSpecies
           conservedVariables(i, nDimensions + 2 + k) = conservedVariables(i,1) *            &
                random(0.0_wp, 1.0_wp)
        end do

        assert(conservedVariables(i,1) > 0.0_wp)

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

     end do

     call computeDependentVariables(nDimensions, nSpecies, conservedVariables,               &
          ratioOfSpecificHeats, specificVolume = specificVolume,                             &
          velocity = velocity, temperature = temperature, massFraction = massFraction)
     assert(all(specificVolume > 0.0_wp))
     assert(all(temperature > 0.0_wp))
     assert(all(massFraction >= 0.0_wp))
     assert(all(massFraction <= 1.0_wp))

     do i = 1, n

        localConservedVariables = conservedVariables(i,:)
        localVelocity = velocity(i,:)
        localMassFraction = massFraction(i,:)

        do j = 1, nDimensions

           localMetricsAlongDirection = metrics(i,1+nDimensions*(j-1):nDimensions*j)

           call computeJacobianOfInviscidFlux(nDimensions, nSpecies,                         &
                localConservedVariables, localMetricsAlongDirection,                         &
                ratioOfSpecificHeats, jacobianOfInviscidFlux,                                &
                specificVolume = specificVolume(i), velocity = localVelocity,                &
                temperature = temperature(i), massFraction = localMassFraction)
           call computeIncomingJacobianOfInviscidFlux(nDimensions, nSpecies,              &
                localConservedVariables, localMetricsAlongDirection,                      &
                ratioOfSpecificHeats, 1, incomingJacobianOfInviscidFlux1,                 &
                specificVolume = specificVolume(i), velocity = localVelocity,             &
                temperature = temperature(i), massFraction = localMassFraction)

           call lapackComputeIncomingJacobianOfInviscidFlux(jacobianOfInviscidFlux, 1,       &
                incomingJacobianOfInviscidFlux2)

           success = success .and.                                                           &
                (maxval(abs(incomingJacobianOfInviscidFlux1 -                                &
                incomingJacobianOfInviscidFlux2)) <= sqrt(epsilon(0.0_wp)))
           if (.not. success) exit

        end do

        if (.not. success) exit

     end do

     SAFE_DEALLOCATE(incomingJacobianOfInviscidFlux2)
     SAFE_DEALLOCATE(incomingJacobianOfInviscidFlux1)
     SAFE_DEALLOCATE(jacobianOfInviscidFlux)
     SAFE_DEALLOCATE(localMetricsAlongDirection)
     SAFE_DEALLOCATE(localMassFraction)
     SAFE_DEALLOCATE(localVelocity)
     SAFE_DEALLOCATE(localConservedVariables)
     SAFE_DEALLOCATE(massFraction)
     SAFE_DEALLOCATE(temperature)
     SAFE_DEALLOCATE(velocity)
     SAFE_DEALLOCATE(specificVolume)
     SAFE_DEALLOCATE(metrics)
     SAFE_DEALLOCATE(conservedVariables)

     if (.not. success) exit

  end do

  call cleanupErrorHandler()

  if (.not. success) stop -1
  stop 0

end program incoming_inviscid_flux_jacobian

subroutine lapackComputeIncomingJacobianOfInviscidFlux(jacobianMatrix, incomingDirection,    &
     incomingJacobianMatrix)

  ! <<< Arguments >>>
  SCALAR_TYPE, intent(in) :: jacobianMatrix(:,:)
  integer, intent(in) :: incomingDirection
  SCALAR_TYPE, intent(out) :: incomingJacobianMatrix(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, LWORK, INFO
  integer, parameter :: LWMAX = 1024
  SCALAR_TYPE, allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
  SCALAR_TYPE :: WORK(LWMAX)
  integer, allocatable :: IPIV(:)

  allocate(WR(size(jacobianMatrix, 1)), WI(size(jacobianMatrix, 1)))
  allocate(VL(size(jacobianMatrix, 1), size(jacobianMatrix, 2)))
  allocate(VR(size(jacobianMatrix, 1), size(jacobianMatrix, 2)))
  allocate(IPIV(size(jacobianMatrix, 1)))

#ifdef SCALAR_TYPE_IS_real64_iso
  LWORK = -1
  call DGEEV('N', 'V', size(jacobianMatrix, 1), jacobianMatrix, size(jacobianMatrix, 2), WR, &
       WI, VL, size(jacobianMatrix, 1), VR, size(jacobianMatrix, 1), WORK, LWORK, INFO)
  LWORK = min(LWMAX, int(WORK(1)))
  call DGEEV('N', 'V', size(jacobianMatrix, 1), jacobianMatrix, size(jacobianMatrix, 2), WR, &
       WI, VL, size(jacobianMatrix, 1), VR, size(jacobianMatrix, 1), WORK, LWORK, INFO)
  VL = VR
  call DGETRF(size(jacobianMatrix, 1), size(jacobianMatrix, 2),                              &
       VL, size(jacobianMatrix, 1), IPIV, INFO)
  LWORK = -1
  call DGETRI(size(jacobianMatrix, 2), VL, size(jacobianMatrix, 1), IPIV, WORK, LWORK, INFO)
  LWORK = min(LWMAX, int(WORK(1)))
  call DGETRI(size(jacobianMatrix, 2), VL, size(jacobianMatrix, 1), IPIV, WORK, LWORK, INFO)
#endif

  do i = 1, size(jacobianMatrix, 1)
     if (incomingDirection * real(WR(i), wp) < 0.0_wp) WR(i) = 0.0_wp
  end do

  do j = 1, size(incomingJacobianMatrix, 2)
     do i = 1, size(incomingJacobianMatrix, 1)
        incomingJacobianMatrix(i,j) = sum(VR(i,:) * WR * VL(:,j))
     end do
  end do

  SAFE_DEALLOCATE(IPIV)
  SAFE_DEALLOCATE(VR)
  SAFE_DEALLOCATE(VL)
  SAFE_DEALLOCATE(WI)
  SAFE_DEALLOCATE(WR)

end subroutine lapackComputeIncomingJacobianOfInviscidFlux
