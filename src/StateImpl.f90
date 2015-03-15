#include "config.h"

module StateImpl

  implicit none
  public

contains

  subroutine allocateData(this, simulationFlags, nGridPoints, nDimensions)

    ! <<< Derived types >>>
    use State_type
    use SimulationFlags_type

    ! <<< Arguments >>>
    type(t_State) :: this
    type(t_SimulationFlags), intent(in) :: simulationFlags
    integer, intent(in) :: nGridPoints, nDimensions

    allocate(this%conservedVariables(nGridPoints, this%nUnknowns))
    allocate(this%rightHandSide(nGridPoints, this%nUnknowns))
    allocate(this%specificVolume(nGridPoints, 1))
    allocate(this%velocity(nGridPoints, nDimensions))
    allocate(this%pressure(nGridPoints, 1))
    allocate(this%temperature(nGridPoints, 1))

    if (simulationFlags%useTargetState) then
       allocate(this%targetState(nGridPoints, this%nUnknowns))
    end if

    if (simulationFlags%viscosityOn) then

       allocate(this%dynamicViscosity(nGridPoints, 1))
       allocate(this%secondCoefficientOfViscosity(nGridPoints, 1))
       allocate(this%thermalDiffusivity(nGridPoints, 1))

       if (simulationFlags%repeatFirstDerivative) then
          allocate(this%stressTensor(nGridPoints, nDimensions ** 2))
          allocate(this%heatFlux(nGridPoints, nDimensions))
       else
          allocate(this%velocityGradient(nGridPoints, nDimensions ** 2))
       end if

    end if

    if (.not. simulationFlags%predictionOnly) then
       allocate(this%adjointVariables(nGridPoints, this%nUnknowns))
    end if

  end subroutine allocateData

end module StateImpl

subroutine setupState(this, grid, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_type
  use State_type
  use SolverOptions_type
  use SimulationFlags_type

  ! <<< Private members >>>
  use StateImpl, only : allocateData

  ! <<< Public members >>>
  use State_mod, only : cleanupState, makeQuiescent

  ! <<< Internal modules >>>
  use InputHelper, only : getOption
  use AcousticSource_mod, only : setupAcousticSource

  implicit none

  ! <<< Arguments >>>
  type(t_State) :: this
  type(t_Grid) :: grid
  type(t_SimulationFlags), intent(in), optional :: simulationFlags
  type(t_SolverOptions), intent(in), optional :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, n, nDimensions, ierror
  type(t_SimulationFlags) :: simulationFlags_
  real(SCALAR_KIND) :: ratioOfSpecificHeats, temp(3)
  character(len = STRING_LENGTH) :: key

  call cleanupState(this)

  call MPI_Cartdim_get(grid%comm, nDimensions, ierror)

  call initializeSimulationFlags(simulationFlags_)
  if (present(simulationFlags)) simulationFlags_ = simulationFlags

  this%nUnknowns = nDimensions + 2
  call allocateData(this, simulationFlags_, grid%nGridPoints, nDimensions)

  ratioOfSpecificHeats = 1.4_wp
  if (present(solverOptions)) ratioOfSpecificHeats = solverOptions%ratioOfSpecificHeats
  call makeQuiescent(this, nDimensions, ratioOfSpecificHeats)

  n = min(getOption("number_of_acoustic_sources", 0), 99)
  if (n > 0) then
     allocate(this%acousticSources(n))
     do i = 1, n
        write(key, '(A,I2.2,A)') "acoustic_source", i, "/"
        temp(1) = getOption(trim(key) // "x", 0.0_wp)
        temp(2) = getOption(trim(key) // "y", 0.0_wp)
        temp(3) = getOption(trim(key) // "z", 0.0_wp)
        call setupAcousticSource(this%acousticSources(i), temp,                              &
             getOption(trim(key) // "amplitude", 1.0_wp),                                    &
             getOption(trim(key) // "frequency", 1.0_wp),                                    &
             getOption(trim(key) // "radius", 1.0_wp),                                       &
             getOption(trim(key) // "phase", 0.0_wp))
     end do
  end if

end subroutine setupState

subroutine cleanupState(this)

  ! <<< Derived types >>>
  use State_type

  implicit none

  ! <<< Arguments >>>
  type(t_State) :: this

  SAFE_DEALLOCATE(this%acousticSources)
  SAFE_DEALLOCATE(this%conservedVariables)
  SAFE_DEALLOCATE(this%targetState)
  SAFE_DEALLOCATE(this%adjointVariables)
  SAFE_DEALLOCATE(this%rightHandSide)
  SAFE_DEALLOCATE(this%specificVolume)
  SAFE_DEALLOCATE(this%velocity)
  SAFE_DEALLOCATE(this%pressure)
  SAFE_DEALLOCATE(this%temperature)
  SAFE_DEALLOCATE(this%dynamicViscosity)
  SAFE_DEALLOCATE(this%secondCoefficientOfViscosity)
  SAFE_DEALLOCATE(this%thermalDiffusivity)
  SAFE_DEALLOCATE(this%velocityGradient)
  SAFE_DEALLOCATE(this%stressTensor)
  SAFE_DEALLOCATE(this%heatFlux)
  SAFE_DEALLOCATE(this%meanPressure)

  this%nUnknowns = 0

end subroutine cleanupState

subroutine loadStateData(this, grid, quantityOfInterest, filename, offset, success)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_type
  use State_type

  ! <<< Internal modules >>>
  use PLOT3DHelper

  implicit none

  ! <<< Arguments >>>
  type(t_State) :: this
  type(t_Grid) :: grid
  integer, intent(in) :: quantityOfInterest
  character(len = *), intent(in) :: filename
  integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
  logical, intent(out) :: success

  select case (quantityOfInterest)
  case (QOI_FORWARD_STATE, QOI_TARGET_STATE, QOI_ADJOINT_STATE)
     call plot3dReadSingleAuxiliarySolutionData(grid%comm, trim(filename),                   &
          offset, this%plot3dAuxiliaryData, success)
  end select

  select case (quantityOfInterest)
  case (QOI_FORWARD_STATE)
     call plot3dReadSingleSolution(grid%comm, trim(filename), offset,                        &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%conservedVariables, success)
  case (QOI_TARGET_STATE)
     call plot3dReadSingleSolution(grid%comm, trim(filename), offset,                        &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%targetState, success)
  case (QOI_ADJOINT_STATE)
     call plot3dReadSingleSolution(grid%comm, trim(filename), offset,                        &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%adjointVariables, success)
  case (QOI_MEAN_PRESSURE)
     call plot3dReadSingleFunction(grid%comm, trim(filename), offset,                        &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%meanPressure, success)
  end select

end subroutine loadStateData

subroutine saveStateData(this, grid, quantityOfInterest, filename, offset, success)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_type
  use State_type

  ! <<< Internal modules >>>
  use Grid_mod, only : computeGradient
  use CNSHelper, only : computeVorticityMagnitudeAndDilatation
  use PLOT3DHelper

  implicit none

  ! <<< Arguments >>>
  type(t_State) :: this
  type(t_Grid) :: grid
  integer, intent(in) :: quantityOfInterest
  character(len = *), intent(in) :: filename
  integer(kind = MPI_OFFSET_KIND), intent(inout) :: offset
  logical, intent(out) :: success

  ! <<< Local variables >>>
  integer :: nDimensions, ierror
  SCALAR_TYPE, allocatable :: f(:,:)

  if (quantityOfInterest == QOI_VORTICITY_DILATATION) then

     call MPI_Cartdim_get(grid%comm, nDimensions, ierror)

     if (.not. allocated(this%velocityGradient) .and.                                        &
          .not. allocated(this%stressTensor)) then
        allocate(f(grid%nGridPoints, max(2, nDimensions ** 2)))
     else
        allocate(f(grid%nGridPoints, 2))
     end if

     if (allocated(this%velocityGradient)) then
        call computeVorticityMagnitudeAndDilatation(nDimensions, this%velocityGradient,      &
             dilatation = f(:,1), vorticityMagnitude = f(:,2))
     else if (allocated(this%stressTensor)) then
        call computeGradient(grid, this%velocity, this%stressTensor)
        call computeVorticityMagnitudeAndDilatation(nDimensions, this%stressTensor,          &
             dilatation = f(:,1), vorticityMagnitude = f(:,2))
     else
        call computeGradient(grid, this%velocity, f)
        call computeVorticityMagnitudeAndDilatation(nDimensions, f,                          &
             dilatation = f(:,1), vorticityMagnitude = f(:,2))
     end if

  end if

  select case (quantityOfInterest)
  case (QOI_FORWARD_STATE, QOI_TARGET_STATE, QOI_ADJOINT_STATE)
     call plot3dWriteSingleAuxiliarySolutionData(grid%comm, trim(filename),                  &
          offset, this%plot3dAuxiliaryData, success)
  end select

  select case (quantityOfInterest)
  case (QOI_FORWARD_STATE)
     call plot3dWriteSingleSolution(grid%comm, trim(filename), offset,                       &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%conservedVariables, success)
  case (QOI_TARGET_STATE)
     call plot3dWriteSingleSolution(grid%comm, trim(filename), offset,                       &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%targetState, success)
  case (QOI_ADJOINT_STATE)
     call plot3dWriteSingleSolution(grid%comm, trim(filename), offset,                       &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%adjointVariables, success)
  case (QOI_MEAN_PRESSURE)
     call plot3dWriteSingleFunction(grid%comm, trim(filename), offset,                       &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%meanPressure, success)
  case (QOI_VORTICITY_DILATATION)
     call plot3dWriteSingleFunction(grid%comm, trim(filename), offset,                       &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize, f(:,1:2), success)
  end select

  SAFE_DEALLOCATE(f)

end subroutine saveStateData

function getFileType(quantityOfInterest) result(fileType)

  ! <<< Derived types >>>
  use State_type
  use PLOT3DDescriptor_type

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: quantityOfInterest

  ! <<< Result >>>
  integer :: fileType

  fileType = -1

  select case (quantityOfInterest)
  case (QOI_FORWARD_STATE, QOI_TARGET_STATE, QOI_ADJOINT_STATE)
     fileType = PLOT3D_SOLUTION_FILE
  case (QOI_MEAN_PRESSURE)
     fileType = PLOT3D_FUNCTION_FILE
  case (QOI_VORTICITY_DILATATION)
     fileType = PLOT3D_FUNCTION_FILE
  end select

end function getFileType

function getNumberOfScalars(quantityOfInterest, nDimensions) result(nScalars)

  use State_type

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: quantityOfInterest, nDimensions

  ! <<< Result >>>
  integer :: nScalars

  select case (quantityOfInterest)
  case (QOI_MEAN_PRESSURE)
     nScalars = 1
  case (QOI_VORTICITY_DILATATION)
     nScalars = 2
  end select

end function getNumberOfScalars

subroutine makeQuiescent(this, nDimensions, ratioOfSpecificHeats, conservedVariables)

  ! <<< Derived types >>>
  use State_type

  implicit none

  ! <<< Arguments >>>
  type(t_State) :: this
  integer, intent(in) :: nDimensions
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
  SCALAR_TYPE, intent(out), optional :: conservedVariables(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  if (present(conservedVariables)) then
     conservedVariables(:,1) = 1.0_wp
     conservedVariables(:,2:nDimensions+1) = 0.0_wp
     conservedVariables(:,nDimensions+2) = 1.0_wp / ratioOfSpecificHeats /                   &
          (ratioOfSpecificHeats - 1.0_wp)
  else
     this%conservedVariables(:,1) = 1.0_wp
     this%conservedVariables(:,2:nDimensions+1) = 0.0_wp
     this%conservedVariables(:,nDimensions+2) = 1.0_wp / ratioOfSpecificHeats /              &
          (ratioOfSpecificHeats - 1.0_wp)
  end if

end subroutine makeQuiescent

subroutine updateState(this, grid, time, simulationFlags, solverOptions, conservedVariables)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_type
  use State_type
  use SolverOptions_type
  use SimulationFlags_type

  ! <<< Internal modules >>>
  use Grid_mod
  use CNSHelper
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  type(t_State) :: this
  type(t_Grid) :: grid
  real(SCALAR_KIND), intent(in) :: time
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  SCALAR_TYPE, intent(in), optional :: conservedVariables(:,:)

  ! <<< Local variables >>>
  integer :: i, nDimensions, ierror

  call startTiming("updateState")

  call MPI_Cartdim_get(grid%comm, nDimensions, ierror)

  if (present(conservedVariables)) then
     call computeDependentVariables(nDimensions, conservedVariables,                         &
          solverOptions%ratioOfSpecificHeats, this%specificVolume(:,1), this%velocity,       &
          this%pressure(:,1), this%temperature(:,1))
  else
     call computeDependentVariables(nDimensions, this%conservedVariables,                    &
          solverOptions%ratioOfSpecificHeats, this%specificVolume(:,1), this%velocity,       &
          this%pressure(:,1), this%temperature(:,1))
  end if

  if (simulationFlags%viscosityOn) then

     call computeTransportVariables(this%temperature(:,1), solverOptions%powerLawExponent,   &
          solverOptions%bulkViscosityRatio, solverOptions%ratioOfSpecificHeats,              &
          solverOptions%reynoldsNumberInverse, solverOptions%prandtlNumberInverse,           &
          this%dynamicViscosity(:,1), this%secondCoefficientOfViscosity(:,1),                &
          this%thermalDiffusivity(:,1))

     if (simulationFlags%repeatFirstDerivative) then

        call computeGradient(grid, this%velocity, this%stressTensor)
        call computeStressTensor(nDimensions, this%stressTensor, this%dynamicViscosity(:,1), &
             this%secondCoefficientOfViscosity(:,1))

        call computeGradient(grid, this%temperature(:,1), this%heatFlux)
        do i = 1, nDimensions
           this%heatFlux(:,i) = - this%thermalDiffusivity(:,1) * this%heatFlux(:,i)
        end do

     else

        call computeGradient(grid, this%velocity, this%velocityGradient)

     end if

  end if

  if (simulationFlags%useConstantCfl) then

     this%cfl = solverOptions%cfl

     if (simulationFlags%viscosityOn) then
        this%timeStepSize = computeTimeStepSize(nDimensions, grid%iblank,                    &
             grid%jacobian(:,1), grid%metrics, this%velocity, this%temperature(:,1),         &
             this%cfl, solverOptions%ratioOfSpecificHeats, this%dynamicViscosity(:,1),       &
             this%thermalDiffusivity(:,1))
     else
        this%timeStepSize = computeTimeStepSize(nDimensions, grid%iblank,                    &
             grid%jacobian(:,1), grid%metrics, this%velocity, this%temperature(:,1),         &
             this%cfl, solverOptions%ratioOfSpecificHeats)
     end if

  else

     this%timeStepSize = solverOptions%timeStepSize

     if (simulationFlags%viscosityOn) then
        this%cfl = computeCfl(nDimensions, grid%iblank, grid%jacobian(:,1), grid%metrics,    &
             this%velocity, this%temperature(:,1), this%timeStepSize,                        &
             solverOptions%ratioOfSpecificHeats, this%dynamicViscosity(:,1),                 &
             this%thermalDiffusivity(:,1))
     else
        this%cfl = computeCfl(nDimensions, grid%iblank, grid%jacobian(:,1), grid%metrics,    &
             this%velocity, this%temperature(:,1), this%timeStepSize,                        &
             solverOptions%ratioOfSpecificHeats)
     end if
     this%cfl = real(0.0, SCALAR_KIND)

  end if

  call endTiming("updateState")

end subroutine updateState

subroutine computeRhsForward(this, grid, patches, time, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_type
  use Patch_type
  use State_type
  use Region_type, only : FORWARD
  use SolverOptions_type
  use PatchDescriptor_type
  use SimulationFlags_type

  ! <<< Internal modules >>>
  use CNSHelper
  use Patch_mod, only : collectAtPatch, addFarFieldPenalty, addWallPenalty
  use MPITimingsHelper, only : startTiming, endTiming
  use StencilOperator_mod, only : applyOperator

  implicit none

  ! <<< Arguments >>>
  type(t_State) :: this
  type(t_Grid) :: grid
  type(t_Patch), allocatable :: patches(:)
  real(SCALAR_KIND), intent(in) :: time
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer :: i, nDimensions, ierror
  SCALAR_TYPE, allocatable :: fluxes1(:,:,:), fluxes2(:,:,:)

  call startTiming("computeRhsForward")

  call MPI_Cartdim_get(grid%comm, nDimensions, ierror)

  allocate(fluxes1(grid%nGridPoints, this%nUnknowns, nDimensions))
  allocate(fluxes2(grid%nGridPoints, this%nUnknowns, nDimensions))

  ! Compute Cartesian form of inviscid fluxes.
  call computeCartesianInvsicidFluxes(nDimensions, this%conservedVariables,                  &
       this%velocity, this%pressure(:,1), fluxes1)

  ! Compute Cartesian form of viscous fluxes if viscous terms are included and computed using
  ! repeated first derivatives.
  if (simulationFlags%viscosityOn .and. simulationFlags%repeatFirstDerivative) then
     call computeCartesianViscousFluxes(nDimensions, this%velocity,                          &
          this%stressTensor, this%heatFlux, fluxes2)
     if (allocated(patches)) then
        do i = 1, size(patches)
           if (patches(i)%gridIndex /= grid%index .or.                                       &
                (patches(i)%patchType /= SAT_FAR_FIELD .and.                                 &
                patches(i)%patchType /= SAT_BLOCK_INTERFACE)) cycle
           call collectAtPatch(patches(i), fluxes2,                                          &
                patches(i)%viscousFluxes) !... save viscous fluxes on patches for later use.
        end do
     end if
     fluxes1 = fluxes1 - fluxes2 !... Cartesian form of total fluxes.
  end if

  ! Transform fluxes from Cartesian to contravariant form:
  ! `fluxes1` has the Cartesian form of total fluxes... upon return, `fluxes2` has the
  ! contravariant form.
  call transformFluxes(nDimensions, fluxes1, grid%metrics, fluxes2, grid%isCurvilinear)

  SAFE_DEALLOCATE(fluxes1) !... no longer needed.

  ! Take derivatives of fluxes.
  do i = 1, nDimensions
     call applyOperator(grid%firstDerivative(i), fluxes2(:,:,i), grid%localSize)
  end do
  this%rightHandSide = this%rightHandSide - sum(fluxes2, dim = 3) !... update right-hand side.

  ! Add dissipation if required.
  if (simulationFlags%dissipationOn) then
     do i = 1, nDimensions
        fluxes2(:,:,i) = this%conservedVariables !... dissipation based on high-order even
                                                 !... derivatives of conserved variables.
        call applyOperator(grid%dissipation(i), fluxes2(:,:,i), grid%localSize)
     end do
     this%rightHandSide = this%rightHandSide +                                               &
          solverOptions%dissipationAmount * sum(fluxes2, dim = 3) !... update right-hand side.
  end if

  SAFE_DEALLOCATE(fluxes2) !... no longer needed

  call endTiming("computeRhsForward")

end subroutine computeRhsForward

subroutine computeRhsAdjoint(this, grid, patches, time, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Grid_type
  use Patch_type
  use State_type
  use Region_type, only : ADJOINT
  use SolverOptions_type
  use PatchDescriptor_type
  use SimulationFlags_type

  ! <<< Internal modules >>>
  use CNSHelper
  use Patch_mod, only : addFarFieldPenalty, addWallPenalty
  use MPITimingsHelper, only : startTiming, endTiming
  use StencilOperator_mod, only : applyOperator

  implicit none

  ! <<< Arguments >>>
  type(t_State) :: this
  type(t_Grid) :: grid
  type(t_Patch), allocatable, intent(in) :: patches(:)
  real(SCALAR_KIND), intent(in) :: time
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer :: i, j, nDimensions, ierror
  SCALAR_TYPE, allocatable :: gradientOfAdjointVariables(:,:,:), localConservedVariables(:), &
       localVelocity(:), localMetricsAlongDirection(:), localStressTensor(:),                &
       localHeatFlux(:), localJacobian1(:,:), localJacobian2(:,:)

  call startTiming("computeRhsAdjoint")

  call MPI_Cartdim_get(grid%comm, nDimensions, ierror)

  ! Partial derivatives of adjoint variables w.r.t. *computational* coordinates.
  allocate(gradientOfAdjointVariables(grid%nGridPoints, this%nUnknowns, nDimensions))
  do i = 1, nDimensions
     gradientOfAdjointVariables(:,:,i) = this%adjointVariables
     call applyOperator(grid%adjointFirstDerivative(i), gradientOfAdjointVariables(:,:,i),   &
          grid%localSize)
  end do

  allocate(localJacobian1(nDimensions + 2, nDimensions + 2))
  allocate(localConservedVariables(nDimensions + 2))
  allocate(localVelocity(nDimensions))
  allocate(localMetricsAlongDirection(nDimensions))

  if (simulationFlags%viscosityOn) then
     allocate(localJacobian2(nDimensions + 2, nDimensions + 2))
     allocate(localStressTensor(nDimensions ** 2))
     allocate(localHeatFlux(nDimensions))
  end if

  do i = 1, grid%nGridPoints

     localConservedVariables = this%conservedVariables(i,:)
     localVelocity = this%velocity(i,:)
     if (simulationFlags%viscosityOn) then
        localStressTensor = this%stressTensor(i,:)
        localHeatFlux = this%heatFlux(i,:)
     end if

     do j = 1, nDimensions

        localMetricsAlongDirection = grid%metrics(i,1+nDimensions*(j-1):nDimensions*j)

        select case (nDimensions)
        case (1)
           call computeJacobianOfInviscidFlux1D(localConservedVariables,                     &
                localMetricsAlongDirection, solverOptions%ratioOfSpecificHeats,              &
                localJacobian1, specificVolume = this%specificVolume(i,1),                   &
                velocity = localVelocity, temperature = this%temperature(i,1))
        case (2)
           call computeJacobianOfInviscidFlux2D(localConservedVariables,                     &
                localMetricsAlongDirection, solverOptions%ratioOfSpecificHeats,              &
                localJacobian1, specificVolume = this%specificVolume(i,1),                   &
                velocity = localVelocity, temperature = this%temperature(i,1))
        case (3)
           call computeJacobianOfInviscidFlux3D(localConservedVariables,                     &
                localMetricsAlongDirection, solverOptions%ratioOfSpecificHeats,              &
                localJacobian1, specificVolume = this%specificVolume(i,1),                   &
                velocity = localVelocity, temperature = this%temperature(i,1))
        end select

        if (simulationFlags%viscosityOn) then
           select case (nDimensions)
           case (1)
              call computeFirstPartialViscousJacobian1D(localConservedVariables,             &
                   localMetricsAlongDirection, localStressTensor, localHeatFlux,             &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localJacobian2, specificVolume = this%specificVolume(i,1),                &
                   velocity = localVelocity, temperature = this%temperature(i,1))
           case (2)
              call computeFirstPartialViscousJacobian2D(localConservedVariables,             &
                   localMetricsAlongDirection, localStressTensor, localHeatFlux,             &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localJacobian2, specificVolume = this%specificVolume(i,1),                &
                   velocity = localVelocity, temperature = this%temperature(i,1))
           case (3)
              call computeFirstPartialViscousJacobian3D(localConservedVariables,             &
                   localMetricsAlongDirection, localStressTensor, localHeatFlux,             &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localJacobian2, specificVolume = this%specificVolume(i,1),                &
                   velocity = localVelocity, temperature = this%temperature(i,1))
           end select
           localJacobian1 = localJacobian1 - localJacobian2
        end if

        this%rightHandSide(i,:) = this%rightHandSide(i,:) -                                  &
             matmul(transpose(localJacobian1), gradientOfAdjointVariables(i,:,j))

     end do

  end do

  if (simulationFlags%viscosityOn) then
     SAFE_DEALLOCATE(localHeatFlux)
     SAFE_DEALLOCATE(localStressTensor)
     SAFE_DEALLOCATE(localJacobian2)
  end if

  SAFE_DEALLOCATE(localMetricsAlongDirection)
  SAFE_DEALLOCATE(localVelocity)
  SAFE_DEALLOCATE(localConservedVariables)
  SAFE_DEALLOCATE(localJacobian1)

  SAFE_DEALLOCATE(gradientOfAdjointVariables)

  call endTiming("computeRhsAdjoint")

end subroutine computeRhsAdjoint

subroutine addPenaltiesForward(this, grid, patches, time, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_type
  use Patch_type
  use State_type
  use Region_type, only : FORWARD
  use SolverOptions_type
  use PatchDescriptor_type
  use SimulationFlags_type

  ! <<< Internal modules >>>
  use CNSHelper
  use Patch_mod, only : collectAtPatch, addFarFieldPenalty, addWallPenalty
  use MPITimingsHelper, only : startTiming, endTiming
  use StencilOperator_mod, only : applyOperator

  implicit none

  ! <<< Arguments >>>
  type(t_State) :: this
  type(t_Grid) :: grid
  type(t_Patch), allocatable :: patches(:)
  real(SCALAR_KIND), intent(in) :: time
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer :: i, nDimensions, ierror

  call startTiming("addPenaltiesForward")

  call MPI_Cartdim_get(grid%comm, nDimensions, ierror)

  if (allocated(patches)) then
     do i = 1, size(patches)
        if (patches(i)%gridIndex /= grid%index) cycle
        select case (patches(i)%patchType)

        case (SAT_FAR_FIELD) !... non-reflecting far-field.
           call addFarFieldPenalty(patches(i), FORWARD, this%rightHandSide, grid%iblank,     &
                nDimensions, solverOptions%ratioOfSpecificHeats,                             &
                this%conservedVariables, this%targetState)

        case (SAT_SLIP_WALL) !... impenetrable wall (for inviscid flow).
           call addWallPenalty(patches(i), FORWARD, this%rightHandSide, grid%iblank,         &
                nDimensions, solverOptions%ratioOfSpecificHeats, this%conservedVariables)

        end select
     end do
  end if

  call endTiming("addPenaltiesForward")

end subroutine addPenaltiesForward

subroutine addPenaltiesAdjoint(this, grid, patches, time, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_type
  use Patch_type
  use State_type
  use Region_type, only : ADJOINT
  use SolverOptions_type
  use PatchDescriptor_type
  use SimulationFlags_type

  ! <<< Internal modules >>>
  use CNSHelper
  use Patch_mod, only : collectAtPatch, addFarFieldPenalty, addWallPenalty
  use MPITimingsHelper, only : startTiming, endTiming
  use StencilOperator_mod, only : applyOperator

  implicit none

  ! <<< Arguments >>>
  type(t_State) :: this
  type(t_Grid) :: grid
  type(t_Patch), allocatable :: patches(:)
  real(SCALAR_KIND), intent(in) :: time
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer :: i, nDimensions, ierror

  call startTiming("addPenaltiesAdjoint")

  call MPI_Cartdim_get(grid%comm, nDimensions, ierror)

  if (allocated(patches)) then
     do i = 1, size(patches)
        if (patches(i)%gridIndex /= grid%index) cycle
        select case (patches(i)%patchType)

        case (SAT_FAR_FIELD) !... non-reflecting far-field.
           call addFarFieldPenalty(patches(i), ADJOINT, this%rightHandSide, grid%iblank,     &
                nDimensions, solverOptions%ratioOfSpecificHeats,                             &
                this%conservedVariables, this%targetState, this%adjointVariables)

        case (SAT_SLIP_WALL) !... impenetrable wall (for inviscid flow).
           call addWallPenalty(patches(i), ADJOINT, this%rightHandSide, grid%iblank,         &
                nDimensions, solverOptions%ratioOfSpecificHeats, this%conservedVariables,    &
                this%adjointVariables)

        end select
     end do
  end if

  call endTiming("addPenaltiesAdjoint")

end subroutine addPenaltiesAdjoint

subroutine addSourcesForward(this, grid, patches, time, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Grid_type
  use Patch_type
  use State_type
  use Region_type, only : FORWARD
  use SolverOptions_type
  use PatchDescriptor_type
  use SimulationFlags_type

  ! <<< Internal modules >>>
  use Patch_mod, only : addDamping, addSolenoidalExcitation
  use MPITimingsHelper, only : startTiming, endTiming
  use AcousticSource_mod, only : addAcousticSource

  implicit none

  ! <<< Arguments >>>
  type(t_State) :: this
  type(t_Grid) :: grid
  type(t_Patch), allocatable, intent(in) :: patches(:)
  real(SCALAR_KIND), intent(in) :: time
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer :: i

  call startTiming("addSourcesForward")

  if (allocated(this%acousticSources)) then
     do i = 1, size(this%acousticSources)
        call addAcousticSource(this%acousticSources(i), time, grid%coordinates,              &
             grid%iblank, this%rightHandSide)
     end do
  end if

  if (allocated(patches)) then
     do i = 1, size(patches)
        if (patches(i)%gridIndex /= grid%index) cycle
        select case (patches(i)%patchType)

        case (SPONGE) !... damp towards target state in sponge zones.
           call addDamping(patches(i), FORWARD, this%rightHandSide, grid%iblank,             &
                this%conservedVariables, this%targetState)

        case (SOLENOIDAL_EXCITATION) !... momentum forcing to excite a 2D shear layer.
           call addSolenoidalExcitation(patches(i), grid%coordinates,                        &
                grid%iblank, time, this%rightHandSide)

        end select
     end do
  end if

  call endTiming("addSourcesForward")

end subroutine addSourcesForward

subroutine addSourcesAdjoint(this, grid, patches, time, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Grid_type
  use Patch_type
  use State_type
  use Region_type, only : ADJOINT
  use SolverOptions_type
  use PatchDescriptor_type
  use SimulationFlags_type

  ! <<< Internal modules >>>
  use Patch_mod, only : addDamping
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  type(t_State) :: this
  type(t_Grid) :: grid
  type(t_Patch), allocatable, intent(in) :: patches(:)
  real(SCALAR_KIND), intent(in) :: time
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer :: i

  call startTiming("addSourcesAdjoint")

  if (allocated(patches)) then
     do i = 1, size(patches)
        if (patches(i)%gridIndex /= grid%index) cycle
        select case (patches(i)%patchType)

        case (SPONGE) !... adjoint of damping source term added to forward RHS.
           call addDamping(patches(i), ADJOINT, this%rightHandSide, grid%iblank,             &
                this%adjointVariables)

        end select
     end do
  end if

  call endTiming("addSourcesAdjoint")

end subroutine addSourcesAdjoint

subroutine updatePatches(this, grid, patches, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_type
  use Patch_type
  use State_type
  use SolverOptions_type
  use PatchDescriptor_type
  use SimulationFlags_type

  ! <<< Public members >>>
  use State_mod, only : updateState

  ! <<< Internal modules >>>
  use Grid_mod, only : computeNormalizedCurveLengths
  use CNSHelper, only : computeCartesianViscousFluxes
  use State_mod, only : updateState
  use Patch_mod, only : collectAtPatch, updateSolenoidalExcitationStrength
  use MPITimingsHelper, only : startTiming, endTiming
  use AcousticSource_mod, only : addAcousticSource

  implicit none

  ! <<< Arguments >>>
  type(t_State) :: this
  type(t_Grid) :: grid
  type(t_Patch), allocatable :: patches(:)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, nDimensions, ierror
  logical :: flag
  real(wp), allocatable :: targetViscousFluxes(:,:,:)

  call startTiming("updatePatches")

  call MPI_Cartdim_get(grid%comm, nDimensions, ierror)

  if (simulationFlags%viscosityOn) then

     flag = .false.
     if (allocated(patches)) flag = any(patches(:)%gridIndex == grid%index .and.             &
          patches(:)%patchType == SAT_FAR_FIELD)
     call MPI_Allreduce(MPI_IN_PLACE, flag, 1, MPI_LOGICAL, MPI_LOR, grid%comm, ierror)

     if (flag) then
        allocate(targetViscousFluxes(grid%nGridPoints, nDimensions + 2, nDimensions))
        call updateState(this, grid, 0.0_wp, simulationFlags,                                &
             solverOptions, this%targetState)
        call computeCartesianViscousFluxes(nDimensions, this%velocity,                       &
             this%stressTensor, this%heatFlux, targetViscousFluxes)
     end if

  end if

  if (allocated(patches)) then
     do i = 1, size(patches)
        if (patches(i)%gridIndex /= grid%index) cycle

        if (simulationFlags%viscosityOn) then

           select case (patches(i)%patchType)

           case (SAT_FAR_FIELD)
              call collectAtPatch(patches(i), targetViscousFluxes,                           &
                   patches(i)%targetViscousFluxes)

           end select

        end if

        select case (patches(i)%patchType)

        case (SAT_FAR_FIELD, SAT_SLIP_WALL, SAT_ISOTHERMAL_WALL,                             &
             SAT_ADIABATIC_WALL, SAT_BLOCK_INTERFACE)
           call collectAtPatch(patches(i), grid%metrics, patches(i)%metrics)
           patches(i)%inviscidPenaltyAmount = patches(i)%inviscidPenaltyAmount /             &
                grid%firstDerivative(abs(patches(i)%normalDirection))%normBoundary(1)
           patches(i)%viscousPenaltyAmount = patches(i)%viscousPenaltyAmount /               &
                grid%firstDerivative(abs(patches(i)%normalDirection))%normBoundary(1)

        case (SOLENOIDAL_EXCITATION)
           call updateSolenoidalExcitationStrength(patches(i), grid%coordinates, grid%iblank)

        end select

     end do
  end if

  SAFE_DEALLOCATE(targetViscousFluxes)

  call endTiming("updatePatches")

end subroutine updatePatches
