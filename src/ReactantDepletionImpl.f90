#include "config.h"

subroutine setupReactantDepletion(this, region)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use ReactantDepletion_mod, only : t_ReactantDepletion

  ! <<< Enumerations >>>
  use State_enum, only : QOI_DUMMY_FUNCTION

  ! <<< Internal modules >>>
  use CNSHelper, only : computeDependentVariables
  use InputHelper, only : getOption, getRequiredOption
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_ReactantDepletion) :: this
  class(t_Region) :: region

  ! <<< Local variables >>>
  character(len = STRING_LENGTH) :: message, reactant

  assert(allocated(region%states))
  assert(size(region%states) > 0)
  assert(region%solverOptions%nSpecies > 0)

  call getRequiredOption("deplete_reactant", reactant, region%comm)
  select case (trim(reactant))
  case ("H2")
     this%reactant = region%states(1)%combustion%H2
  case ("O2")
     this%reactant = region%states(1)%combustion%O2
  case default
     write(message, '(A)') "WARNING, unknown reactant!"
     call gracefulExit(region%comm, message)
  end select

  assert(this%reactant > 0)
  assert(this%reactant <= region%solverOptions%nSpecies)

  call this%cleanup()

  call this%setupBase(region%simulationFlags, region%solverOptions)

end subroutine setupReactantDepletion

subroutine cleanupReactantDepletion(this)

  ! <<< Derived types >>>
  use ReactantDepletion_mod, only : t_ReactantDepletion

  implicit none

  ! <<< Arguments >>>
  class(t_ReactantDepletion) :: this

  call this%cleanupBase()

end subroutine cleanupReactantDepletion

function computeReactantDepletion(this, region) result(instantaneousFunctional)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Region_mod, only : t_Region
  use ReactantDepletion_mod, only : t_ReactantDepletion

  ! <<< Arguments >>>
  class(t_ReactantDepletion) :: this
  class(t_Region), intent(in) :: region

  ! <<< Result >>>
  SCALAR_TYPE :: instantaneousFunctional

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, ierror
  SCALAR_TYPE, allocatable :: F(:,:)

  assert(allocated(region%grids))
  assert(allocated(region%states))
  assert(size(region%grids) == size(region%states))

  assert(this%reactant > 0)
  assert(this%reactant <= region%solverOptions%nSpecies)

  instantaneousFunctional = 0.0_wp

  do i = 1, size(region%grids)

     assert(region%grids(i)%nGridPoints > 0)
     assert(allocated(region%grids(i)%targetMollifier))
     assert(size(region%grids(i)%targetMollifier, 1) == region%grids(i)%nGridPoints)
     assert(size(region%grids(i)%targetMollifier, 2) == 1)
     assert(allocated(region%states(i)%massFraction))
     assert(size(region%states(i)%massFraction, 1) == region%grids(i)%nGridPoints)
     assert(size(region%states(i)%massFraction, 2) == region%solverOptions%nSpecies)

     allocate(F(region%grids(i)%nGridPoints, 1))
     F(:,1) = region%states(i)%massFraction(:, this%reactant)
     instantaneousFunctional = instantaneousFunctional +                                     &
          region%grids(i)%computeInnerProduct(F, F, region%grids(i)%targetMollifier(:,1))
     SAFE_DEALLOCATE(F)

  end do

  if (region%commGridMasters /= MPI_COMM_NULL)                                               &
       call MPI_Allreduce(MPI_IN_PLACE, instantaneousFunctional, 1,                          &
       SCALAR_TYPE_MPI, MPI_SUM, region%commGridMasters, ierror)

  do i = 1, size(region%grids)
     call MPI_Bcast(instantaneousFunctional, 1, SCALAR_TYPE_MPI,                             &
          0, region%grids(i)%comm, ierror)
  end do

  this%cachedValue = instantaneousFunctional

end function computeReactantDepletion

subroutine computeReactantDepletionAdjointForcing(this, simulationFlags, solverOptions,      &
     grid, state, patch)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use ReactantDepletion_mod, only : t_ReactantDepletion
  use SolverOptions_mod, only : t_SolverOptions
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_ReactantDepletion) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state
  class(t_CostTargetPatch) :: patch

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, gridIndex, patchIndex
  SCALAR_TYPE :: F

  nDimensions = grid%nDimensions
  assert_key(nDimensions, (1, 2, 3))
  assert(solverOptions%nSpecies > 0)

  assert(this%reactant > 0)
  assert(this%reactant <= solverOptions%nSpecies)

  do k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)
     do j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
        do i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
           gridIndex = i - patch%gridOffset(1) + patch%gridLocalSize(1) *                    &
                (j - 1 - patch%gridOffset(2) + patch%gridLocalSize(2) *                      &
                (k - 1 - patch%gridOffset(3)))
           if (grid%iblank(gridIndex) == 0) cycle
           patchIndex = i - patch%offset(1) + patch%localSize(1) *                           &
                (j - 1 - patch%offset(2) + patch%localSize(2) *                              &
                (k - 1 - patch%offset(3)))

           F = - 2.0_wp * grid%targetMollifier(gridIndex, 1) *                               &
                state%massFraction(gridIndex, this%reactant) *                               &
                state%specificVolume(gridIndex, 1)

           patch%adjointForcing(patchIndex,:) = 0.0_wp
           patch%adjointForcing(patchIndex,1) = - F *                                        &
                state%massFraction(gridIndex, this%reactant)
           patch%adjointForcing(patchIndex,nDimensions+2+this%reactant) = F

        end do !... i = patch%offset(1) + 1, patch%offset(1) + patch%localSize(1)
     end do !... j = patch%offset(2) + 1, patch%offset(2) + patch%localSize(2)
  end do !... k = patch%offset(3) + 1, patch%offset(3) + patch%localSize(3)

end subroutine computeReactantDepletionAdjointForcing

function isReactantDepletionPatchValid(this, patchDescriptor, gridSize, normalDirection,     &
     extent, simulationFlags, message) result(isPatchValid)

  ! <<< Derived types >>>
  use ReactantDepletion_mod, only : t_ReactantDepletion
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_ReactantDepletion) :: this
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  integer, intent(in) :: gridSize(:), normalDirection, extent(6)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Result >>>
  logical :: isPatchValid

  ! <<< Local variables >>>
  integer :: i, n

  isPatchValid = .false.

  n = size(gridSize)

  do i = 1, size(gridSize)
     if (extent((i-1)*2+1) == extent((i-1)*2+2)) n = n - 1
  end do

  if (n /= size(gridSize)) then
     write(message, '(2(A,I0.0),A)') "Expected a ", size(gridSize),                          &
          "D patch, but extent represents a ", n, "D patch!"
     return
  end if

  isPatchValid = .true.

end function isReactantDepletionPatchValid
