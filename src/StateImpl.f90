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

    ! <<< Arguments >>>
    class(t_State) :: this
    type(t_SimulationFlags), intent(in) :: simulationFlags
    type(t_SolverOptions), intent(in) :: solverOptions
    integer, intent(in) :: nGridPoints, nDimensions

    allocate(this%conservedVariables(nGridPoints, solverOptions%nUnknowns))
    allocate(this%rightHandSide(nGridPoints, solverOptions%nUnknowns))
    allocate(this%specificVolume(nGridPoints, 1))
    allocate(this%velocity(nGridPoints, nDimensions))
    allocate(this%pressure(nGridPoints, 1))
    allocate(this%temperature(nGridPoints, 1))

    if (simulationFlags%useTargetState) then
       allocate(this%targetState(nGridPoints, solverOptions%nUnknowns))
    end if

    if (simulationFlags%computeTimeAverage) then
       allocate(this%timeAverage(nGridPoints, solverOptions%nUnknowns))
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

    if (simulationFlags%enableAdjoint) then
       allocate(this%adjointVariables(nGridPoints, solverOptions%nUnknowns))
    end if

    if (simulationFlags%enableIBM) then
      ! Allocate levelset arrays
      allocate(this%levelset(nGridPoints, 1))
      allocate(this%levelsetNormal(nGridPoints, nDimensions))
      ! allocate(this%levelsetCurvature(grid%nGridPoints, 1))
      ! allocate(this%indicatorFunction(grid%nGridPoints, 1))
      allocate(this%ibmDissipation(nGridPoints, solverOptions%nUnknowns))
      allocate(this%nDotGradRho(nGridPoints, 1))
      allocate(this%uDotGradRho(nGridPoints, 1))
      allocate(this%objectVelocity(nGridPoints, nDimensions))
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
  type(t_SimulationFlags), intent(in), optional :: simulationFlags
  type(t_SolverOptions), intent(in), optional :: solverOptions

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, n
  type(t_SimulationFlags) :: simulationFlags_
  type(t_SolverOptions) :: solverOptions_
  real(SCALAR_KIND) :: ratioOfSpecificHeats, temp(3)
  character(len = STRING_LENGTH) :: key

  call this%cleanup()

  assert_key(grid%nDimensions, (1, 2, 3))

  if (present(simulationFlags)) then
     simulationFlags_ = simulationFlags
  else
     call simulationFlags_%initialize()
  end if

  if (present(solverOptions)) then
     solverOptions_ = solverOptions
  else
     call solverOptions_%initialize(grid%nDimensions, simulationFlags_, grid%comm)
  end if

  this%nSpecies = solverOptions_%nSpecies

  assert(grid%nGridPoints > 0)
  call allocateData(this, simulationFlags_, solverOptions_,                                  &
       grid%nGridPoints, grid%nDimensions)

  ratioOfSpecificHeats = 1.4_wp
  if (present(solverOptions)) then
     assert(solverOptions%ratioOfSpecificHeats > 1.0_wp)
     ratioOfSpecificHeats = solverOptions%ratioOfSpecificHeats
  end if
  call this%makeQuiescent(grid%nDimensions, ratioOfSpecificHeats)

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
  SAFE_DEALLOCATE(this%timeAverage)

  ! deallocate levelset arrays
  SAFE_DEALLOCATE(this%levelset)
  SAFE_DEALLOCATE(this%levelsetNormal)
  ! SAFE_DEALLOCATE(this%levelsetCurvature)
  ! SAFE_DEALLOCATE(this%indicatorFunction)
  ! SAFE_DEALLOCATE(this%primitiveGridNorm)
  SAFE_DEALLOCATE(this%ibmDissipation)
  SAFE_DEALLOCATE(this%nDotGradRho)
  SAFE_DEALLOCATE(this%uDotGradRho)
  SAFE_DEALLOCATE(this%objectVelocity)

  this%adjointForcingFactor = 1.0_wp

end subroutine cleanupState

subroutine loadStateData(this, grid, quantityOfInterest, filename,                           &
     offset, success, speciesFilename, speciesFileOffset)

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
  character(len = *), intent(in), optional :: speciesFilename
  integer(kind = MPI_OFFSET_KIND), intent(inout), optional :: speciesFileOffset

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: nDimensions, nUnknowns

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nUnknowns = size(this%conservedVariables, 2)
  assert(nUnknowns == nDimensions + 2 + this%nSpecies)

#ifdef DEBUG
  if (quantityOfInterest == QOI_DUMMY_FUNCTION) then
     assert(associated(this%dummyFunction))
     assert(size(this%dummyFunction, 1) == grid%nGridPoints)
     assert(size(this%dummyFunction, 2) > 0)
  end if
#endif

#ifdef DEBUG
  if (this%nSpecies > 0) then
     assert(present(speciesFilename))
     assert(present(speciesFileOffset))
  end if
#endif

  select case (quantityOfInterest)
  case (QOI_FORWARD_STATE, QOI_TARGET_STATE, QOI_ADJOINT_STATE,                              &
       QOI_RIGHT_HAND_SIDE, QOI_TIME_AVERAGED_STATE)
     call plot3dReadSingleAuxiliarySolutionData(grid%comm, trim(filename),                   &
          offset, this%plot3dAuxiliaryData, success)
  end select

  select case (quantityOfInterest)

  case (QOI_FORWARD_STATE)
     call plot3dReadSingleSolution(grid%comm, trim(filename), offset,                        &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%conservedVariables(:,1:nDimensions+2), success)
     if (this%nSpecies > 0)                                                                  &
          call plot3dReadSingleFunction(grid%comm, trim(speciesFilename),                    &
          speciesFileOffset, grid%mpiDerivedTypeScalarSubarray, grid%globalSize,             &
          this%conservedVariables(:,nDimensions+3:), success)

  case (QOI_TARGET_STATE)
     call plot3dReadSingleSolution(grid%comm, trim(filename), offset,                        &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%targetState(:,1:nDimensions+2), success)
     if (this%nSpecies > 0)                                                                  &
          call plot3dReadSingleFunction(grid%comm, trim(speciesFilename),                    &
          speciesFileOffset, grid%mpiDerivedTypeScalarSubarray, grid%globalSize,             &
          this%targetState(:,nDimensions+3:), success)

  case (QOI_ADJOINT_STATE)
     call plot3dReadSingleSolution(grid%comm, trim(filename), offset,                        &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%adjointVariables(:,1:nDimensions+2), success)
     if (this%nSpecies > 0)                                                                  &
          call plot3dReadSingleFunction(grid%comm, trim(speciesFilename),                    &
          speciesFileOffset, grid%mpiDerivedTypeScalarSubarray, grid%globalSize,             &
          this%adjointVariables(:,nDimensions+3:), success)

  case (QOI_RIGHT_HAND_SIDE)
     call plot3dReadSingleSolution(grid%comm, trim(filename), offset,                        &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%rightHandSide(:,1:nDimensions+2), success)
     if (this%nSpecies > 0)                                                                  &
          call plot3dReadSingleFunction(grid%comm, trim(speciesFilename),                    &
          speciesFileOffset, grid%mpiDerivedTypeScalarSubarray, grid%globalSize,             &
          this%rightHandSide(:,nDimensions+3:), success)

  case (QOI_TIME_AVERAGED_STATE)
     call plot3dReadSingleSolution(grid%comm, trim(filename), offset,                        &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%timeAverage(:,1:nDimensions+2), success)
     if (this%nSpecies > 0)                                                                  &
          call plot3dReadSingleFunction(grid%comm, trim(speciesFilename),                    &
          speciesFileOffset, grid%mpiDerivedTypeScalarSubarray, grid%globalSize,             &
          this%timeAverage(:,nDimensions+3:), success)

  case (QOI_DUMMY_FUNCTION)
     assert(associated(this%dummyFunction))
     assert(size(this%dummyFunction, 1) == grid%nGridPoints)
     assert(size(this%dummyFunction, 2) > 0)
     call plot3dReadSingleFunction(grid%comm, trim(filename), offset,                        &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%dummyFunction, success)
     nullify(this%dummyFunction)
  end select

end subroutine loadStateData

subroutine saveStateData(this, grid, quantityOfInterest, filename,                           &
     offset, success, speciesFilename, speciesFileOffset)

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
  character(len = *), intent(in), optional :: speciesFilename
  integer(kind = MPI_OFFSET_KIND), intent(inout), optional :: speciesFileOffset

  ! <<< Local variables >>>
  integer :: nDimensions, nUnknowns

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  nUnknowns = size(this%conservedVariables, 2)
  assert(nUnknowns == nDimensions + 2 + this%nSpecies)

#ifdef DEBUG
  if (quantityOfInterest == QOI_DUMMY_FUNCTION) then
     assert(associated(this%dummyFunction))
     assert(size(this%dummyFunction, 1) == grid%nGridPoints)
     assert(size(this%dummyFunction, 2) > 0)
  end if
#endif

#ifdef DEBUG
  if (this%nSpecies > 0) then
     assert(present(speciesFilename))
     assert(present(speciesFileOffset))
  end if
#endif

  select case (quantityOfInterest)
  case (QOI_FORWARD_STATE, QOI_TARGET_STATE, QOI_RIGHT_HAND_SIDE,                            &
       QOI_TIME_AVERAGED_STATE, QOI_ADJOINT_STATE)
     this%plot3dAuxiliaryData(4) = this%time
     call plot3dWriteSingleAuxiliarySolutionData(grid%comm, trim(filename),                  &
          offset, this%plot3dAuxiliaryData, success)
  end select

  select case (quantityOfInterest)

  case (QOI_FORWARD_STATE)
     call plot3dWriteSingleSolution(grid%comm, trim(filename), offset,                       &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%conservedVariables(:,1:nDimensions+2), success)
     if (this%nSpecies > 0)                                                                  &
          call plot3dWriteSingleFunction(grid%comm, trim(speciesFilename),                   &
          speciesFileOffset, grid%mpiDerivedTypeScalarSubarray, grid%globalSize,             &
          this%conservedVariables(:,nDimensions+3:), success)

  case (QOI_TARGET_STATE)
     call plot3dWriteSingleSolution(grid%comm, trim(filename), offset,                       &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%targetState(:,1:nDimensions+2), success)
     if (this%nSpecies > 0)                                                                  &
          call plot3dWriteSingleFunction(grid%comm, trim(speciesFilename),                   &
          speciesFileOffset, grid%mpiDerivedTypeScalarSubarray, grid%globalSize,             &
          this%targetState(:,nDimensions+3:), success)

  case (QOI_ADJOINT_STATE)
     call plot3dWriteSingleSolution(grid%comm, trim(filename), offset,                       &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%adjointVariables(:,1:nDimensions+2), success)
     if (this%nSpecies > 0)                                                                  &
          call plot3dWriteSingleFunction(grid%comm, trim(speciesFilename),                   &
          speciesFileOffset, grid%mpiDerivedTypeScalarSubarray, grid%globalSize,             &
          this%adjointVariables(:,nDimensions+3:), success)

  case (QOI_RIGHT_HAND_SIDE)
     call plot3dWriteSingleSolution(grid%comm, trim(filename), offset,                       &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%rightHandSide(:,1:nDimensions+2), success)
     if (this%nSpecies > 0)                                                                  &
          call plot3dWriteSingleFunction(grid%comm, trim(speciesFilename),                   &
          speciesFileOffset, grid%mpiDerivedTypeScalarSubarray, grid%globalSize,             &
          this%rightHandSide(:,nDimensions+3:), success)

  case (QOI_TIME_AVERAGED_STATE)
     call plot3dWriteSingleSolution(grid%comm, trim(filename), offset,                       &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%timeAverage(:,1:nDimensions+2), success)
     if (this%nSpecies > 0)                                                                  &
          call plot3dWriteSingleFunction(grid%comm, trim(speciesFilename),                   &
          speciesFileOffset, grid%mpiDerivedTypeScalarSubarray, grid%globalSize,             &
          this%timeAverage(:,nDimensions+3:), success)

  case (QOI_DUMMY_FUNCTION)
     call plot3dWriteSingleFunction(grid%comm, trim(filename), offset,                       &
          grid%mpiDerivedTypeScalarSubarray, grid%globalSize,                                &
          this%dummyFunction, success)
     nullify(this%dummyFunction)

  end select

end subroutine saveStateData

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
  integer :: i, nDimensions

  call startTiming("updateState")

  nDimensions = grid%nDimensions
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
  integer :: nDimensions

  nDimensions = grid%nDimensions
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
  integer :: nDimensions

  nDimensions = grid%nDimensions
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

subroutine addSources(this, mode, grid)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  implicit none

  ! <<< Arguments >>>
  class(t_State) :: this
  integer, intent(in) :: mode
  class(t_Grid) :: grid

  ! <<< Local variables >>>
  integer :: i

  call startTiming("addSources")

  if (mode == FORWARD .and. allocated(this%acousticSources)) then
     do i = 1, size(this%acousticSources)
        call this%acousticSources(i)%add(this%time, grid%coordinates,                        &
             grid%iblank, this%rightHandSide)
     end do
  end if

  call endTiming("addSources")

end subroutine addSources
