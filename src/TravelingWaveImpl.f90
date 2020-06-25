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
    region%params%buffer = getOption("traveling_wave/initial_speed",1.0_wp)

    SAFE_DEALLOCATE(region%paramWeights%buffer)
    allocate(region%paramWeights%buffer(1,1))
    region%paramWeights%buffer = getOption("traveling_wave/parameter_weight",1.0_wp)

    region%params%buffer = region%params%buffer / region%paramWeights%buffer

  end subroutine setupTravelingWave

  subroutine computeTravelingWaveResidual(region, residual)

    ! <<< Derived types >>>
    use Region_mod, only : t_Region
    use RegionVector_mod, only : t_RegionVector

    ! <<< Enumerations >>>
    use Region_enum, only : FORWARD

    ! <<< Internal modules >>>
    use CNSHelper, only: transformFluxes

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: region
    class(t_RegionVector) :: residual

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, nDimensions, nUnknowns
    SCALAR_TYPE :: convectionSpeed
    SCALAR_TYPE, allocatable :: F(:,:), fluxes1(:,:,:), fluxes2(:,:,:)

    convectionSpeed = region%params%buffer(1,1) * region%paramWeights%buffer(1,1)

    call region%computeRhs(FORWARD,1,1)

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
            - convectionSpeed * residual%states(i)%conservedVariables           &
            - region%states(i)%rightHandSide

      SAFE_DEALLOCATE(fluxes2)
    end do

  end subroutine computeTravelingWaveResidual

  function computeTravelingWaveJacobian(region, input, mode) result(output)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_mod, only : t_Region
    use RegionVector_mod

    ! <<< Enumerations >>>
    use State_enum, only : QOI_ADJOINT_STATE
    use Region_enum, only : ADJOINT, LINEARIZED

    ! <<< Internal modules >>>
    use CNSHelper, only: transformFluxes
    use RegionVectorImpl, only: loadRegionVector

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: region
    class(t_RegionVector), intent(in) :: input
    integer, intent(in) :: mode
    type(t_RegionVector) :: output

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, direction, index, nDimensions, nUnknowns, ierror
    SCALAR_TYPE, allocatable :: fluxes1(:,:,:), fluxes2(:,:,:)
    SCALAR_TYPE :: convectionSpeed, paramWeight
    type(t_RegionVector) :: temp

    nUnknowns = region%solverOptions%nUnknowns
    nDimensions = nUnknowns - 2
    assert_key(nDimensions, (1, 2, 3))

    direction = 1
    convectionSpeed = region%params%buffer(1,1) * region%paramWeights%buffer(1,1)
    paramWeight = region%paramWeights%buffer(1,1)

    call output%set(input)

    call loadRegionVector(region,input,QOI_ADJOINT_STATE)
    call region%computeRhs(mode,1,1)
    !NOTE: current adjoint rhs includes additional negative sign.
    if (mode==ADJOINT) then
      do i = 1, size(region%states)
        region%states(i)%rightHandSide = - region%states(i)%rightHandSide
      end do
    end if

    do i = 1, size(region%states)
      assert(size(region%states(i)%rightHandSide,1)==region%grids(i)%nGridPoints)
      assert(size(region%states(i)%rightHandSide,2)==nUnknowns)
      allocate(fluxes1(region%grids(i)%nGridPoints, nUnknowns, nDimensions))

      select case(mode)
      case(ADJOINT)
        do j = 1, nDimensions
          fluxes1(:,:,j) = region%states(i)%adjointVariables
          call region%grids(i)%adjointFirstDerivative(j)%apply(fluxes1(:,:,j),    &
                                                        region%grids(i)%localSize)
        end do

        output%states(i)%conservedVariables = 0.0_wp
        do k = 1, nUnknowns
          do j = 1, nDimensions
            index = direction + (j-1) * nDimensions
            output%states(i)%conservedVariables(:,k) =                            &
                                          output%states(i)%conservedVariables(:,k)&
                               + region%grids(i)%metrics(:,index) * fluxes1(:,k,j)
          end do
          output%states(i)%conservedVariables(:,k) =                              &
          output%states(i)%conservedVariables(:,k) * region%grids(i)%jacobian(:,1)
        end do
      case(LINEARIZED)
        fluxes1 = 0.0_wp
        fluxes1(:,:,1) = region%states(i)%adjointVariables

        allocate(fluxes2(region%grids(i)%nGridPoints, nUnknowns, nDimensions))
        call transformFluxes(nDimensions, fluxes1, region%grids(i)%metrics,       &
                             fluxes2, region%grids(i)%isCurvilinear)
        do j = 1, nDimensions
          call region%grids(i)%firstDerivative(j)%apply(fluxes2(:,:,j),          &
                                                        region%grids(i)%localSize)
        end do
        output%states(i)%conservedVariables = sum(fluxes2, dim = 3)
        SAFE_DEALLOCATE(fluxes2)

        do j = 1, nUnknowns
          output%states(i)%conservedVariables(:,j) =                                &
            output%states(i)%conservedVariables(:,j) * region%grids(i)%jacobian(:,1)
        end do
      end select
      output%states(i)%conservedVariables =                                     &
              - convectionSpeed * output%states(i)%conservedVariables           &
              - region%states(i)%rightHandSide
      SAFE_DEALLOCATE(fluxes1)
    end do

    call temp%set(region)
    temp%params = 0.0_wp

    output%params = 0.0_wp
    do i = 1, size(region%states)
      allocate(fluxes1(region%grids(i)%nGridPoints, nUnknowns, nDimensions))
      allocate(fluxes2(region%grids(i)%nGridPoints, nUnknowns, nDimensions))

      fluxes1 = 0.0_wp
      fluxes1(:,:,1) = region%states(i)%conservedVariables

      call transformFluxes(nDimensions, fluxes1, region%grids(i)%metrics,       &
                           fluxes2, region%grids(i)%isCurvilinear)

      do j = 1, nDimensions
        call region%grids(i)%firstDerivative(j)%apply(fluxes2(:,:,j),          &
                                                      region%grids(i)%localSize)
      end do

      temp%states(i)%conservedVariables = sum(fluxes2, dim = 3)
      do j = 1, nUnknowns
        temp%states(i)%conservedVariables(:,j) =                                &
          temp%states(i)%conservedVariables(:,j) * region%grids(i)%jacobian(:,1)
      end do

      output%states(i)%conservedVariables = output%states(i)%conservedVariables &
            - input%params(1) * paramWeight * temp%states(i)%conservedVariables

      output%params = output%params - paramWeight *                             &
          region%grids(i)%computeInnerProduct(region%states(i)%adjointVariables,&
                                              temp%states(i)%conservedVariables)
    end do

    if (region%commGridMasters /= MPI_COMM_NULL)                                  &
         call MPI_Allreduce(MPI_IN_PLACE, output%params, 1,                       &
         SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

    do i = 1, size(region%grids)
       call MPI_Bcast(output%params, 1, SCALAR_TYPE_MPI, 0, region%grids(i)%comm, ierror)
    end do

    call temp%cleanup()

  end function computeTravelingWaveJacobian

  subroutine saveTravelingWave(region,filename)

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE

    ! <<< Internal modules >>>
    use PLOT3DHelper

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: region
    character(len = *), intent(in) :: filename

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i

    do i = 1, size(region%states)
      region%states(i)%plot3dAuxiliaryData(2) = region%params%buffer(1,1)
    end do

    call region%saveData(QOI_FORWARD_STATE,filename)

  end subroutine saveTravelingWave

  subroutine loadTravelingWave(region,filename)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Region_mod, only : t_Region

    ! <<< Enumerations >>>
    use State_enum, only : QOI_FORWARD_STATE

    ! <<< Internal modules >>>
    use PLOT3DHelper

    implicit none

    ! <<< Arguments >>>
    class(t_Region) :: region
    character(len = *), intent(in) :: filename

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i

    call region%loadData(QOI_FORWARD_STATE,filename)

    region%params%buffer(1,1) = region%states(1)%plot3dAuxiliaryData(2)

  end subroutine loadTravelingWave

end module TravelingWaveImpl
