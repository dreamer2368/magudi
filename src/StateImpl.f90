#include "config.h"

module StateImpl

  implicit none
  public

contains

  subroutine allocateData(this, simulationFlags, solverOptions, nGridPoints, nDimensions)

    ! <<< Derived types >>>
    use State_mod, only : t_State
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Enumerations >>>
    use SolverOptions_enum, only : SOUND

    ! <<< Arguments >>>
    class(t_State) :: this
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
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
       select case (solverOptions%costFunctionalType)
       case (SOUND)
          allocate(this%meanPressure(nGridPoints, 1))
       end select
    end if

  end subroutine allocateData

end module StateImpl

subroutine setupState(this, grid, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Private members >>>
  use StateImpl, only : allocateData

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  class(t_Grid) :: grid
  type(t_SimulationFlags), pointer, intent(in), optional :: simulationFlags
  type(t_SolverOptions), pointer, intent(in), optional :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, n, nDimensions, ierror
  type(t_SimulationFlags) :: simulationFlags_
  type(t_SolverOptions) :: solverOptions_
  real(SCALAR_KIND) :: ratioOfSpecificHeats, temp(3)
  character(len = STRING_LENGTH) :: key

  call this%cleanup()

  call MPI_Cartdim_get(grid%comm, nDimensions, ierror)
  assert_key(nDimensions, (1, 2, 3))

  if (present(simulationFlags)) then
     simulationFlags_ = simulationFlags
  else
     call simulationFlags_%initialize()
  end if

  if (present(solverOptions)) then
     solverOptions_ = solverOptions
  else
     call solverOptions_%initialize(simulationFlags_, grid%comm)
  end if

  assert(grid%nGridPoints > 0)
  this%nUnknowns = nDimensions + 2
  call allocateData(this, simulationFlags_, solverOptions_, grid%nGridPoints, nDimensions)

  ratioOfSpecificHeats = 1.4_wp
  if (present(solverOptions)) then
     assert(solverOptions%ratioOfSpecificHeats > 1.0_wp)
     ratioOfSpecificHeats = solverOptions%ratioOfSpecificHeats
  end if
  call this%makeQuiescent(nDimensions, ratioOfSpecificHeats)

  n = min(getOption("number_of_acoustic_sources", 0), 99)
  if (n > 0) then
     allocate(this%acousticSources(n))
     do i = 1, n
        write(key, '(A,I2.2,A)') "acoustic_source", i, "/"
        temp(1) = getOption(trim(key) // "x", 0.0_wp)
        temp(2) = getOption(trim(key) // "y", 0.0_wp)
        temp(3) = getOption(trim(key) // "z", 0.0_wp)
        call this%acousticSources(i)%setup(temp,                                             &     
             getOption(trim(key) // "amplitude", 1.0_wp),                                    &
             getOption(trim(key) // "frequency", 1.0_wp),                                    &
             getOption(trim(key) // "radius", 1.0_wp),                                       &
             getOption(trim(key) // "phase", 0.0_wp))
     end do
  end if

end subroutine setupState

subroutine cleanupState(this)

  ! <<< Derived types >>>
  use State_mod, only : t_State

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

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
  this%adjointForcingFactor = 1.0_wp

end subroutine cleanupState

subroutine loadStateData(this, grid, quantityOfInterest, filename, offset, success)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State

  ! <<< Enumerations >>>
  use State_enum

  ! <<< Internal modules >>>
  use PLOT3DHelper

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  class(t_Grid) :: grid
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
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State

  ! <<< Enumerations >>>
  use State_enum

  ! <<< Internal modules >>>
  use CNSHelper, only : computeVorticityMagnitudeAndDilatation
  use PLOT3DHelper

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  class(t_Grid) :: grid
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
        call grid%computeGradient(this%velocity, this%stressTensor)
        call computeVorticityMagnitudeAndDilatation(nDimensions, this%stressTensor,          &
             dilatation = f(:,1), vorticityMagnitude = f(:,2))
     else
        call grid%computeGradient(this%velocity, f)
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
  use PLOT3DDescriptor_type

  ! <<< Enumerations >>>
  use State_enum

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

  ! <<< Enumerations >>>
  use State_enum

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: quantityOfInterest, nDimensions

  ! <<< Result >>>
  integer :: nScalars

  assert_key(nDimensions, (1, 2, 3))

  select case (quantityOfInterest)
  case (QOI_MEAN_PRESSURE)
     nScalars = 1
  case (QOI_VORTICITY_DILATATION)
     nScalars = 2
  end select

end function getNumberOfScalars

subroutine makeQuiescent(this, nDimensions, ratioOfSpecificHeats, conservedVariables)

  ! <<< Derived types >>>
  use State_mod, only : t_State

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  integer, intent(in) :: nDimensions
  real(SCALAR_KIND), intent(in) :: ratioOfSpecificHeats
  SCALAR_TYPE, intent(out), optional :: conservedVariables(:,:)

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND

  assert(ratioOfSpecificHeats > 1.0_wp)

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

subroutine updateState(this, grid, simulationFlags, solverOptions, conservedVariables)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use CNSHelper
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  class(t_Grid) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  SCALAR_TYPE, intent(in), optional :: conservedVariables(:,:)

  ! <<< Local variables >>>
  integer :: i, nDimensions, ierror

  call startTiming("updateState")

  call MPI_Cartdim_get(grid%comm, nDimensions, ierror)
  assert_key(nDimensions, (1, 2, 3))

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

        call grid%computeGradient(this%velocity, this%stressTensor)
        call computeStressTensor(nDimensions, this%stressTensor, this%dynamicViscosity(:,1), &
             this%secondCoefficientOfViscosity(:,1))

        call grid%computeGradient(this%temperature(:,1), this%heatFlux)
        do i = 1, nDimensions
           this%heatFlux(:,i) = - this%thermalDiffusivity(:,1) * this%heatFlux(:,i)
        end do

     else

        call grid%computeGradient(this%velocity, this%velocityGradient)

     end if

  end if

  call endTiming("updateState")

end subroutine updateState

function computeStateCfl(this, grid, simulationFlags, solverOptions) result(cfl)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use CNSHelper, only : computeCfl

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  class(t_Grid) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Result >>>
  real(SCALAR_KIND) :: cfl

  ! <<< Local variables >>>
  integer :: nDimensions, ierror

  call MPI_Cartdim_get(grid%comm, nDimensions, ierror)
  assert_key(nDimensions, (1, 2, 3))

  assert(allocated(grid%iblank))
  assert(size(grid%iblank) > 0)
  assert(allocated(grid%metrics))
  assert(all(shape(grid%metrics) > 0))
  assert(allocated(grid%jacobian))
  assert(size(grid%jacobian, 1) > 0)
  assert(size(grid%jacobian, 2) == 1)
  assert(allocated(this%velocity))
  assert(all(shape(this%velocity) > 0))
  assert(allocated(this%temperature))
  assert(size(this%temperature, 1) > 0)
  assert(size(this%temperature, 2) == 1)

  if (simulationFlags%useConstantCfl) then
     cfl = solverOptions%cfl
  else
     if (simulationFlags%viscosityOn) then

        assert(allocated(this%dynamicViscosity))
        assert(size(this%dynamicViscosity, 1) > 0)
        assert(size(this%dynamicViscosity, 2) == 1)
        assert(allocated(this%thermalDiffusivity))
        assert(size(this%thermalDiffusivity, 1) > 0)
        assert(size(this%thermalDiffusivity, 2) == 1)

        cfl = computeCfl(nDimensions, grid%iblank, grid%jacobian(:,1), grid%metrics,         &
             this%velocity, this%temperature(:,1), solverOptions%timeStepSize,               &
             solverOptions%ratioOfSpecificHeats, this%dynamicViscosity(:,1),                 &
             this%thermalDiffusivity(:,1))
     else
        cfl = computeCfl(nDimensions, grid%iblank, grid%jacobian(:,1), grid%metrics,         &
             this%velocity, this%temperature(:,1), solverOptions%timeStepSize,               &
             solverOptions%ratioOfSpecificHeats)
     end if
  end if

end function computeStateCfl

function computeStateTimeStepSize(this, grid, simulationFlags,                               &
     solverOptions) result(timeStepSize)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use CNSHelper, only : computeTimeStepSize

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  class(t_Grid) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Result >>>
  real(SCALAR_KIND) :: timeStepSize

  ! <<< Local variables >>>
  integer :: nDimensions, ierror

  call MPI_Cartdim_get(grid%comm, nDimensions, ierror)
  assert_key(nDimensions, (1, 2, 3))

  assert(allocated(grid%iblank))
  assert(size(grid%iblank) > 0)
  assert(allocated(grid%metrics))
  assert(all(shape(grid%metrics) > 0))
  assert(allocated(grid%jacobian))
  assert(size(grid%jacobian, 1) > 0)
  assert(size(grid%jacobian, 2) == 1)
  assert(allocated(this%velocity))
  assert(all(shape(this%velocity) > 0))
  assert(allocated(this%temperature))
  assert(size(this%temperature, 1) > 0)
  assert(size(this%temperature, 2) == 1)

  if (simulationFlags%useConstantCfl) then
     if (simulationFlags%viscosityOn) then

        assert(allocated(this%dynamicViscosity))
        assert(size(this%dynamicViscosity, 1) > 0)
        assert(size(this%dynamicViscosity, 2) == 1)
        assert(allocated(this%thermalDiffusivity))
        assert(size(this%thermalDiffusivity, 1) > 0)
        assert(size(this%thermalDiffusivity, 2) == 1)

        timeStepSize = computeTimeStepSize(nDimensions, grid%iblank, grid%jacobian(:,1),     &
             grid%metrics, this%velocity, this%temperature(:,1), solverOptions%cfl,          &
             solverOptions%ratioOfSpecificHeats, this%dynamicViscosity(:,1),                 &
             this%thermalDiffusivity(:,1))
     else
        timeStepSize = computeTimeStepSize(nDimensions, grid%iblank, grid%jacobian(:,1),     &
             grid%metrics, this%velocity, this%temperature(:,1), solverOptions%cfl,          &
             solverOptions%ratioOfSpecificHeats)
     end if
  else
     timeStepSize = solverOptions%timeStepSize
  end if

end function computeStateTimeStepSize

subroutine computeRhsForward(this, grid, patches, time, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use Patch_type, only : t_Patch
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_type, only : SAT_FAR_FIELD, SAT_BLOCK_INTERFACE
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use CNSHelper
  use Patch_mod, only : collectAtPatch
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  class(t_Grid) :: grid
  type(t_Patch), allocatable :: patches(:)
  real(SCALAR_KIND), intent(in) :: time
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer :: i, nDimensions, ierror
  SCALAR_TYPE, allocatable :: fluxes1(:,:,:), fluxes2(:,:,:)

  call startTiming("computeRhsForward")

  call MPI_Cartdim_get(grid%comm, nDimensions, ierror)
  assert_key(nDimensions, (1, 2, 3))

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
           if (patches(i)%gridIndex /= grid%index) cycle
           select case (patches(i)%patchType)
           case (SAT_FAR_FIELD, SAT_BLOCK_INTERFACE)
              call collectAtPatch(patches(i), fluxes2,                                       &
                   patches(i)%viscousFluxes) !... save viscous fluxes on patch for later use.
           end select
        end do
     end if
     fluxes1 = fluxes1 - fluxes2 !... Cartesian form of total fluxes.
  end if

  ! Transform fluxes from Cartesian to contravariant form: `fluxes1` has the Cartesian form of
  ! total fluxes... upon return, `fluxes2` has the contravariant form.
  call transformFluxes(nDimensions, fluxes1, grid%metrics, fluxes2, grid%isCurvilinear)

  SAFE_DEALLOCATE(fluxes1) !... no longer needed.

  ! Take derivatives of fluxes.
  do i = 1, nDimensions
     call grid%firstDerivative(i)%apply(fluxes2(:,:,i), grid%localSize)
  end do
  this%rightHandSide = this%rightHandSide - sum(fluxes2, dim = 3) !... update right-hand side.

  ! Add dissipation if required.
  if (simulationFlags%dissipationOn) then
     do i = 1, nDimensions
        fluxes2(:,:,i) = this%conservedVariables !... dissipation based on high-order even
                                                 !... derivatives of conserved variables.
        call grid%dissipation(i)%apply(fluxes2(:,:,i), grid%localSize)
     end do
     this%rightHandSide = this%rightHandSide +                                               &
          solverOptions%dissipationAmount * sum(fluxes2, dim = 3) !... update right-hand side.
  end if

  SAFE_DEALLOCATE(fluxes2) !... no longer needed

  ! Collect conserved variables on SAT block interfaces.
  if (allocated(patches)) then
     do i = 1, size(patches)
        if (patches(i)%gridIndex /= grid%index) cycle
        select case (patches(i)%patchType)
        case (SAT_BLOCK_INTERFACE)
           call collectAtPatch(patches(i), this%conservedVariables,                          &
                patches(i)%conservedVariables)
        end select
     end do
  end if

  call endTiming("computeRhsForward")

end subroutine computeRhsForward

subroutine computeRhsAdjoint(this, grid, patches, time, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use Patch_type, only : t_Patch
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use CNSHelper
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  class(t_Grid) :: grid
  type(t_Patch), allocatable, intent(in) :: patches(:)
  real(SCALAR_KIND), intent(in) :: time
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, ierror
  SCALAR_TYPE, allocatable :: temp1(:,:,:), temp2(:,:),                                      &
       localFluxJacobian1(:,:), localFluxJacobian2(:,:), localConservedVariables(:),         &
       localVelocity(:), localMetricsAlongDirection1(:), localMetricsAlongDirection2(:),     &
       localStressTensor(:), localHeatFlux(:), localAdjointDiffusion(:,:)

  call startTiming("computeRhsAdjoint")

  call MPI_Cartdim_get(grid%comm, nDimensions, ierror)
  assert_key(nDimensions, (1, 2, 3))

  allocate(temp1(grid%nGridPoints, this%nUnknowns, nDimensions))

  ! Add dissipation if required.
  if (simulationFlags%dissipationOn) then
     do i = 1, nDimensions
        temp1(:,:,i) = this%adjointVariables
        call grid%dissipation(i)%apply(temp1(:,:,i), grid%localSize)
     end do
     this%rightHandSide = this%rightHandSide -                                               &
          solverOptions%dissipationAmount * sum(temp1, dim = 3) !... update right-hand side.
  end if

  ! Partial derivatives of adjoint variables w.r.t. *computational* coordinates.
  do i = 1, nDimensions
     temp1(:,:,i) = this%adjointVariables
     call grid%adjointFirstDerivative(i)%apply(temp1(:,:,i), grid%localSize)
  end do

  allocate(localFluxJacobian1(this%nUnknowns, this%nUnknowns))
  allocate(localConservedVariables(this%nUnknowns))
  allocate(localVelocity(nDimensions))
  allocate(localMetricsAlongDirection1(nDimensions))

  if (simulationFlags%viscosityOn) then
     allocate(localFluxJacobian2(this%nUnknowns, this%nUnknowns))
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

        localMetricsAlongDirection1 = grid%metrics(i,1+nDimensions*(j-1):nDimensions*j)

        select case (nDimensions)
        case (1)
           call computeJacobianOfInviscidFlux1D(localConservedVariables,                     &
                localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                localFluxJacobian1, specificVolume = this%specificVolume(i,1),               &
                velocity = localVelocity, temperature = this%temperature(i,1))
        case (2)
           call computeJacobianOfInviscidFlux2D(localConservedVariables,                     &
                localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                localFluxJacobian1, specificVolume = this%specificVolume(i,1),               &
                velocity = localVelocity, temperature = this%temperature(i,1))
        case (3)
           call computeJacobianOfInviscidFlux3D(localConservedVariables,                     &
                localMetricsAlongDirection1, solverOptions%ratioOfSpecificHeats,             &
                localFluxJacobian1, specificVolume = this%specificVolume(i,1),               &
                velocity = localVelocity, temperature = this%temperature(i,1))
        end select !... nDimensions

        if (simulationFlags%viscosityOn) then
           select case (nDimensions)
           case (1)
              call computeFirstPartialViscousJacobian1D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian2, specificVolume = this%specificVolume(i,1),            &
                   velocity = localVelocity, temperature = this%temperature(i,1))
           case (2)
              call computeFirstPartialViscousJacobian2D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian2, specificVolume = this%specificVolume(i,1),            &
                   velocity = localVelocity, temperature = this%temperature(i,1))
           case (3)
              call computeFirstPartialViscousJacobian3D(localConservedVariables,             &
                   localMetricsAlongDirection1, localStressTensor, localHeatFlux,            &
                   solverOptions%powerLawExponent, solverOptions%ratioOfSpecificHeats,       &
                   localFluxJacobian2, specificVolume = this%specificVolume(i,1),            &
                   velocity = localVelocity, temperature = this%temperature(i,1))
           end select
           localFluxJacobian1 = localFluxJacobian1 - localFluxJacobian2
        end if

        this%rightHandSide(i,:) = this%rightHandSide(i,:) +                                  &
             matmul(transpose(localFluxJacobian1), temp1(i,:,j))

     end do !... j = 1, nDimensions

  end do !... i = 1, grid%nGridPoints

  SAFE_DEALLOCATE(localConservedVariables)
  SAFE_DEALLOCATE(localFluxJacobian1)

  if (simulationFlags%viscosityOn) then

     SAFE_DEALLOCATE(localHeatFlux)
     SAFE_DEALLOCATE(localStressTensor)
     SAFE_DEALLOCATE(localFluxJacobian2)

     allocate(temp2(grid%nGridPoints, this%nUnknowns - 1))

     allocate(localMetricsAlongDirection2(nDimensions))
     allocate(localFluxJacobian2(this%nUnknowns - 1, this%nUnknowns - 1))
     allocate(localAdjointDiffusion(this%nUnknowns - 1, nDimensions))

     do i = 1, grid%nGridPoints

        localVelocity = this%velocity(i,:)
        localAdjointDiffusion = 0.0_wp

        do j = 1, nDimensions

           localMetricsAlongDirection1 = grid%metrics(i,1+nDimensions*(j-1):nDimensions*j)

           do k = 1, nDimensions

              localMetricsAlongDirection2 = grid%metrics(i,1+nDimensions*(k-1):nDimensions*k)

              select case (nDimensions)
              case (1)
                 call computeSecondPartialViscousJacobian1D(localVelocity,                   &
                      this%dynamicViscosity(i,1), this%secondCoefficientOfViscosity(i,1),    &
                      this%thermalDiffusivity(i,1), grid%jacobian(i,1),                      &
                      localMetricsAlongDirection1(1), localFluxJacobian2)
              case (2)
                 call computeSecondPartialViscousJacobian2D(localVelocity,                   &
                      this%dynamicViscosity(i,1), this%secondCoefficientOfViscosity(i,1),    &
                      this%thermalDiffusivity(i,1), grid%jacobian(i,1),                      &
                      localMetricsAlongDirection2, localMetricsAlongDirection1,              &
                      localFluxJacobian2)
              case (3)
                 call computeSecondPartialViscousJacobian3D(localVelocity,                   &
                      this%dynamicViscosity(i,1), this%secondCoefficientOfViscosity(i,1),    &
                      this%thermalDiffusivity(i,1), grid%jacobian(i,1),                      &
                      localMetricsAlongDirection2, localMetricsAlongDirection1,              &
                      localFluxJacobian2)
              end select !... nDimensions

              localAdjointDiffusion(:,j) = localAdjointDiffusion(:,j) + &
                   matmul(transpose(localFluxJacobian2), temp1(i,2:this%nUnknowns,k))

           end do !... k = 1, nDimensions

        end do !... j = 1, nDimensions

        do j = 1, nDimensions
           temp1(i,2:this%nUnknowns,j) = localAdjointDiffusion(:,j)
        end do

     end do !... i = 1, grid%nGridPoints

     do i = 1, nDimensions
        call grid%adjointFirstDerivative(i)%apply(temp1(:,2:this%nUnknowns,i), grid%localSize)
     end do
     temp2 = sum(temp1(:,2:this%nUnknowns,:), dim = 3)

     temp2(:,nDimensions+1) = solverOptions%ratioOfSpecificHeats *                           &
          this%specificVolume(:,1) * temp2(:,nDimensions+1)
     do i = 1, nDimensions
        temp2(:,i) = this%specificVolume(:,1) * temp2(:,i) -                                 &
             this%velocity(:,i) * temp2(:,nDimensions+1)
     end do

     this%rightHandSide(:,2:this%nUnknowns) = this%rightHandSide(:,2:this%nUnknowns) - temp2
     this%rightHandSide(:,1) = this%rightHandSide(:,1) +                                     &
          this%specificVolume(:,1) * this%conservedVariables(:,nDimensions+2) *              &
          temp2(:,nDimensions+1) + sum(this%velocity * temp2(:,1:nDimensions), dim = 2)

     SAFE_DEALLOCATE(temp2)

     SAFE_DEALLOCATE(localAdjointDiffusion)
     SAFE_DEALLOCATE(localFluxJacobian2)
     SAFE_DEALLOCATE(localMetricsAlongDirection2)

  end if !... simulationFlags%viscosityOn

  SAFE_DEALLOCATE(localMetricsAlongDirection1)
  SAFE_DEALLOCATE(localVelocity)

  SAFE_DEALLOCATE(temp1)

  call endTiming("computeRhsAdjoint")

end subroutine computeRhsAdjoint

subroutine addPenaltiesForward(this, grid, patches, time, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use Patch_type, only : t_Patch
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_type
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD

  ! <<< Internal modules >>>
  use CNSHelper
  use Patch_mod, only : collectAtPatch, addFarFieldPenalty, addWallPenalty
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  class(t_Grid) :: grid
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

        case (SAT_SLIP_WALL, SAT_ISOTHERMAL_WALL, SAT_ADIABATIC_WALL)
           call addWallPenalty(patches(i), FORWARD, this%rightHandSide, grid%iblank,         &
                nDimensions, solverOptions%ratioOfSpecificHeats, this%conservedVariables)

        end select
     end do !... i = 1, size(patches)
  end if

  call endTiming("addPenaltiesForward")

end subroutine addPenaltiesForward

subroutine addPenaltiesAdjoint(this, grid, patches, time, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use Patch_type, only : t_Patch
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_type
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : ADJOINT

  ! <<< Internal modules >>>
  use CNSHelper
  use Patch_mod, only : collectAtPatch, addFarFieldPenalty, addWallPenalty
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  class(t_Grid) :: grid
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
     end do !... i = 1, size(patches)
  end if

  call endTiming("addPenaltiesAdjoint")

end subroutine addPenaltiesAdjoint

subroutine addSourcesForward(this, grid, patches, time)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use Patch_type, only : t_Patch
  use State_mod, only : t_State
  use PatchDescriptor_type, only : SPONGE, SOLENOIDAL_EXCITATION

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD

  ! <<< Internal modules >>>
  use Patch_mod, only : addDamping, addSolenoidalExcitation
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  class(t_Grid) :: grid
  type(t_Patch), allocatable, intent(in) :: patches(:)
  real(SCALAR_KIND), intent(in) :: time

  ! <<< Local variables >>>
  integer :: i

  call startTiming("addSourcesForward")

  if (allocated(this%acousticSources)) then
     do i = 1, size(this%acousticSources)
        call this%acousticSources(i)%add(time, grid%coordinates,                             &
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
     end do !... i = 1, size(patches)
  end if

  call endTiming("addSourcesForward")

end subroutine addSourcesForward

subroutine addSourcesAdjoint(this, grid, patches, time)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use Patch_type, only : t_Patch
  use State_mod, only : t_State
  use PatchDescriptor_type, only : SPONGE

  ! <<< Enumerations >>>
  use Region_enum, only : ADJOINT

  ! <<< Internal modules >>>
  use Patch_mod, only : addDamping
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  class(t_Grid) :: grid
  type(t_Patch), allocatable, intent(in) :: patches(:)
  real(SCALAR_KIND), intent(in) :: time

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
     end do !... i = 1, size(patches)
  end if

  call endTiming("addSourcesAdjoint")

end subroutine addSourcesAdjoint

subroutine updatePatches(this, grid, patches, simulationFlags, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use Patch_type, only : t_Patch
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_type
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Internal modules >>>
  use CNSHelper, only : computeCartesianViscousFluxes, computeDependentVariables
  use Patch_mod, only : collectAtPatch, updateSolenoidalExcitationStrength
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  class(t_Grid) :: grid
  type(t_Patch), allocatable :: patches(:)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer :: i, nDimensions, ierror
  logical :: flag
  SCALAR_TYPE, allocatable :: targetViscousFluxes(:,:,:), targetTemperature(:)

  call startTiming("updatePatches")

  call MPI_Cartdim_get(grid%comm, nDimensions, ierror)

  if (simulationFlags%viscosityOn) then

     flag = .false.
     if (allocated(patches)) flag = any(patches(:)%gridIndex == grid%index .and.             &
          patches(:)%patchType == SAT_ISOTHERMAL_WALL)
     call MPI_Allreduce(MPI_IN_PLACE, flag, 1, MPI_LOGICAL, MPI_LOR, grid%comm, ierror)

     if (flag) then
        allocate(targetTemperature(grid%nGridPoints))
        call computeDependentVariables(nDimensions, this%targetState,                        &
             solverOptions%ratioOfSpecificHeats, temperature = targetTemperature)
     end if

     if (allocated(patches)) then
        do i = 1, size(patches)
           if (patches(i)%gridIndex /= grid%index .or.                                       &
                patches(i)%patchType /= SAT_ISOTHERMAL_WALL) cycle
           call collectAtPatch(patches(i), targetTemperature, patches(i)%wallTemperature)
        end do
     end if

     SAFE_DEALLOCATE(targetTemperature)

     flag = .false.
     if (allocated(patches)) flag = any(patches(:)%gridIndex == grid%index .and.             &
          patches(:)%patchType == SAT_FAR_FIELD)
     call MPI_Allreduce(MPI_IN_PLACE, flag, 1, MPI_LOGICAL, MPI_LOR, grid%comm, ierror)

     if (flag) then
        allocate(targetViscousFluxes(grid%nGridPoints, nDimensions + 2, nDimensions))
        call this%update(grid, simulationFlags, solverOptions, this%targetState)
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

           if (simulationFlags%viscosityOn) then
              patches(i)%viscousPenaltyAmount = patches(i)%viscousPenaltyAmount /            &
                   grid%firstDerivative(abs(patches(i)%normalDirection))%normBoundary(1) *   &
                   solverOptions%reynoldsNumberInverse
           end if

        case (SOLENOIDAL_EXCITATION)
           call updateSolenoidalExcitationStrength(patches(i), grid%coordinates, grid%iblank)

        end select

     end do !... i = 1, size(patches)
  end if

  SAFE_DEALLOCATE(targetViscousFluxes)
  SAFE_DEALLOCATE(targetTemperature)

  call endTiming("updatePatches")

end subroutine updatePatches

subroutine addAdjointForcing(this, grid, patch, solverOptions)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use Patch_type, only : t_Patch
  use State_mod, only : t_State
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_type

  ! <<< Enumerations >>>
  use SolverOptions_enum

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  class(t_Grid) :: grid
  class(t_Patch) :: patch
  type(t_SolverOptions), intent(in) :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, gridIndex, patchIndex, ierror
  SCALAR_TYPE :: temp

  assert(patch%patchType == CONTROL_TARGET)
  assert(patch%gridIndex == grid%index)

  assert_key(solverOptions%costFunctionalType, ( \
  SOUND, \
  LIFT,  \
  DRAG))

  call MPI_Cartdim_get(grid%comm, nDimensions, ierror)
  assert_key(nDimensions, (1, 2, 3))

  do k = patch%offset(3) + 1, patch%offset(3) + patch%patchSize(3)
     do j = patch%offset(2) + 1, patch%offset(2) + patch%patchSize(2)
        do i = patch%offset(1) + 1, patch%offset(1) + patch%patchSize(1)
           gridIndex = i - patch%gridOffset(1) + patch%gridLocalSize(1) *                    &
                (j - 1 - patch%gridOffset(2) + patch%gridLocalSize(2) *                      &
                (k - 1 - patch%gridOffset(3)))
           if (grid%iblank(gridIndex) == 0) cycle
           patchIndex = i - patch%offset(1) + patch%patchSize(1) *                           &
                (j - 1 - patch%offset(2) + patch%patchSize(2) *                              &
                (k - 1 - patch%offset(3)))

           select case (solverOptions%costFunctionalType)

           case (SOUND)
              temp = - 2.0_wp * (solverOptions%ratioOfSpecificHeats - 1.0_wp) *              &
                   this%adjointForcingFactor * grid%targetMollifier(gridIndex, 1) *          &
                   (this%pressure(gridIndex, 1) - this%meanPressure(gridIndex, 1))
              this%rightHandSide(gridIndex,1) = this%rightHandSide(gridIndex,1) +            &
                   0.5_wp * sum(this%velocity(gridIndex,:) ** 2) * temp
              this%rightHandSide(gridIndex,2:nDimensions+1) =                                &
                   this%rightHandSide(gridIndex,2:nDimensions+1) -                           &
                   this%velocity(gridIndex,:) * temp
              this%rightHandSide(gridIndex,nDimensions+2) =                                  &
                   this%rightHandSide(gridIndex,nDimensions+2) + temp

           end select !... solverOptions%costFunctionalType

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%patchSize(1)
     end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%patchSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%patchSize(3)

end subroutine addAdjointForcing
