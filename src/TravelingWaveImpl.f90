#include "config.h"

module TravelingWaveImpl

  implicit none
  public

contains

  subroutine setupTravelingWave(region)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Enumerations >>>
    use State_enum, only : QOI_DUMMY_FUNCTION

    ! <<< Internal modules >>>
    use CNSHelper, only : computeDependentVariables
    use InputHelper, only : getOption, getRequiredOption


    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: region

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND

    SAFE_DEALLOCATE(region%params%buffer)
    allocate(region%params%buffer(1,1))
    region%params%buffer = getOption("traveling_wave/initial_speed",0.0_wp)

  end subroutine setupTravelingWave

  subroutine computeTravelingWaveResidual(region, residual)

    ! <<< Derived types >>>
    use Region_mod, only : t_Region
    use RegionVector_mod, only : t_RegionVector

    ! <<< Internal modules >>>
    use CNSHelper, only: transformFluxes

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: region
    class(t_RegionVector) :: residual

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, nDimensions, nUnknowns
    SCALAR_TYPE, allocatable :: F(:,:), fluxes1(:,:,:), fluxes2(:,:,:)

    ! call residual%set(region) !... seems redundant

    residual%params = 0.0_wp

    nUnknowns = region%solverOptions%nUnknowns
    nDimensions = nUnknowns - 2
    assert_key(nDimensions, (1, 2, 3))

    do i = 1, size(region%states)
      assert(size(region%states(i)%conservedVariables,1)==region%grids(i)%nGridPoints)
      assert(size(region%states(i)%conservedVariables,2)==nUnknowns)
      allocate(fluxes1(region%grids(i)%nGridPoints, nUnknowns, nDimensions))
      allocate(fluxes2(region%grids(i)%nGridPoints, nUnknowns, nDimensions))

      fluxes1 = 0.0_wp
      fluxes1(:,:,1) = region%states(i)%conservedVariables

      call transformFluxes(nDimensions, fluxes1, region%grids(i)%metrics,       &
                           fluxes2, region%grids(i)%isCurvilinear)

      SAFE_DEALLOCATE(fluxes1)

      do j = 1, nDimensions
         call region%grids(i)%firstDerivative(j)%apply(fluxes2(:,:,j),          &
                                                       region%grids(i)%localSize)
      end do

      residual%states(i)%conservedVariables = sum(fluxes2, dim = 3)
      do j = 1, nUnknowns
        residual%states(i)%conservedVariables(:,j) =                            &
        residual%states(i)%conservedVariables(:,j) * region%grids(i)%jacobian(:,1)
      end do
      residual%states(i)%conservedVariables =                                   &
            - region%params%buffer(1,1) * residual%states(i)%conservedVariables &
            - region%states(i)%rightHandSide

      SAFE_DEALLOCATE(fluxes2)
    end do

  end subroutine computeTravelingWaveResidual

end module TravelingWaveImpl
