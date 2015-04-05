#include "config.h"

subroutine setupCostTargetPatch(this, index, comm, patchDescriptor,                          &
     grid, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SolverOptions_mod, only : t_SolverOptions
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_CostTargetPatch) :: this
  integer, intent(in) :: index, comm
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  class(t_Grid), intent(in) :: grid
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions

  call this%setupBase(index, comm, patchDescriptor, grid, simulationFlags, solverOptions)

  this%penaltyInPhysicalCoordinates = .true.

  if (this%nPatchPoints > 0) then
     allocate(this%norm(this%nPatchPoints, 1))
     allocate(this%adjointForcing(this%nPatchPoints, solverOptions%nUnknowns))
  end if

end subroutine setupCostTargetPatch

subroutine cleanupCostTargetPatch(this)

  ! <<< Derived types >>>
  use CostTargetPatch_mod, only : t_CostTargetPatch

  implicit none

  ! <<< Arguments >>>
  class(t_CostTargetPatch) :: this

  call this%cleanupBase()

  SAFE_DEALLOCATE(this%norm)
  SAFE_DEALLOCATE(this%adjointForcing)

end subroutine cleanupCostTargetPatch

subroutine addAdjointForcing(this, mode, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Enumerations >>>
  use Region_enum, only : FORWARD, ADJOINT

  implicit none

  ! <<< Arguments >>>
  class(t_CostTargetPatch) :: this
  integer, intent(in) :: mode
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State) :: state

  ! <<< Local variables >>>
  integer :: i, j, k, l, gridIndex, patchIndex

  assert_key(mode, (FORWARD, ADJOINT))
  assert(this%gridIndex == grid%index)
  assert(all(grid%offset == this%gridOffset))
  assert(all(grid%localSize == this%gridLocalSize))

  if (mode == FORWARD) return

  call startTiming("addCostTargetPenalty")

  do l = 1, solverOptions%nUnknowns
     do k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
        do j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
           do i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
              gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                   &
                   (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                     &
                   (k - 1 - this%gridOffset(3)))
              if (grid%iblank(gridIndex) == 0) cycle
              patchIndex = i - this%offset(1) + this%patchSize(1) *                          &
                   (j - 1 - this%offset(2) + this%patchSize(2) *                             &
                   (k - 1 - this%offset(3)))

              state%rightHandSide(gridIndex,l) = state%rightHandSide(gridIndex,l) +          &
                   grid%targetMollifier(gridIndex,1) *          &
                   this%adjointForcing(patchIndex,l)

           end do !... i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
        end do !... j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
     end do !... k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
  end do !... l = 1, solverOptions%nUnknowns

  call endTiming("addCostTargetPenalty")

end subroutine addAdjointForcing

function verifyCostTargetPatchUsage(this, patchDescriptor, gridSize, normalDirection,        &
     extent, simulationFlags, success, message) result(isPatchUsed)

  ! <<< Derived types >>>
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use PatchDescriptor_mod, only : t_PatchDescriptor
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_CostTargetPatch) :: this
  type(t_PatchDescriptor), intent(in) :: patchDescriptor
  integer, intent(in) :: gridSize(:), normalDirection, extent(6)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  logical, intent(out) :: success
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Result >>>
  logical :: isPatchUsed

  ! <<< Local variables >>>
  integer :: i

  isPatchUsed = .false.

  success = .false.

  do i = 1, size(gridSize)
     if (extent((i-1)*2+1) < 0 .or. extent((i-1)*2+2) > gridSize(i) .or.                     &
          extent((i-1)*2+1) > extent((i-1)*2+2)) then
        write(message, '(A)') "Invalid extent!"
        return
     end if
  end do

  success = .true.

  isPatchUsed = .true.
  if (simulationFlags%predictionOnly) isPatchUsed = .false.

end function verifyCostTargetPatchUsage

subroutine updateCostTargetPatch(this, simulationFlags, solverOptions, grid, state)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use State_mod, only : t_State
  use CostTargetPatch_mod, only : t_CostTargetPatch
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  implicit none

  ! <<< Arguments >>>
  class(t_CostTargetPatch) :: this
  type(t_SimulationFlags), intent(in) :: simulationFlags
  type(t_SolverOptions), intent(in) :: solverOptions
  class(t_Grid), intent(in) :: grid
  class(t_State), intent(in) :: state

  ! Nothing to be done for this patch type...

end subroutine updateCostTargetPatch

function computeScalarInnerProductOnPatch(this, grid, f, g, weight) result(innerProduct)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use CostTargetPatch_mod, only : t_CostTargetPatch

  implicit none

  ! <<< Arguments >>>
  class(t_CostTargetPatch) :: this
  class(t_Grid) :: grid
  SCALAR_TYPE, intent(in) :: f(:), g(:)
  SCALAR_TYPE, intent(in), optional :: weight(:)

  ! <<< Result >>>
  SCALAR_TYPE :: innerProduct

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, gridIndex, patchIndex, ierror

  assert(grid%nGridPoints > 0)
  assert(size(grid%iblank) == grid%nGridPoints)
  assert(size(f) == grid%nGridPoints)
  assert(size(g) == grid%nGridPoints)
#ifdef DEBUG
  if (present(weight)) then
     assert(size(weight) == grid%nGridPoints)
  end if
  if (this%nPatchPoints > 0) then
     assert(allocated(this%norm))
     assert(size(this%norm) == this%nPatchPoints)
  end if
#endif
  assert(this%comm /= MPI_COMM_NULL)

  innerProduct = 0.0_wp

  do k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
     do j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
        do i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
           gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                      &
                (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                        &
                (k - 1 - this%gridOffset(3)))
           if (grid%iblank(gridIndex) == 0) cycle
           patchIndex = i - this%offset(1) + this%patchSize(1) *                             &
                (j - 1 - this%offset(2) + this%patchSize(2) *                                &
                (k - 1 - this%offset(3)))

           if (present(weight)) then
              innerProduct = innerProduct + f(gridIndex) *                                   &
                   this%norm(patchIndex, 1) * g(gridIndex) * weight(gridIndex)
           else
              innerProduct = innerProduct + f(gridIndex) *                                   &
                   this%norm(patchIndex, 1) * g(gridIndex)
           end if

        end do !... i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
  end do !... k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)

#ifdef SCALAR_TYPE_IS_binary128_IEEE754
  assert(allocated(this%mpiReduceBuffer))
  call MPI_Allgather(innerProduct, 1, SCALAR_TYPE_MPI, this%mpiReduceBuffer,                 &
       1, SCALAR_TYPE_MPI, this%comm, ierror)
  innerProduct = sum(this%mpiReduceBuffer)
#else
  call MPI_Allreduce(MPI_IN_PLACE, innerProduct, 1, SCALAR_TYPE_MPI,                         &
       MPI_SUM, this%comm, ierror)
#endif

end function computeScalarInnerProductOnPatch

function computeVectorInnerProductOnPatch(this, grid, f, g, weight) result(innerProduct)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use CostTargetPatch_mod, only : t_CostTargetPatch

  implicit none

  ! <<< Arguments >>>
  class(t_CostTargetPatch) :: this
  class(t_Grid) :: grid
  SCALAR_TYPE, intent(in) :: f(:,:), g(:,:)
  SCALAR_TYPE, intent(in), optional :: weight(:)

  ! <<< Result >>>
  SCALAR_TYPE :: innerProduct

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, l, gridIndex, patchIndex, ierror

  assert(grid%nGridPoints > 0)
  assert(size(grid%iblank) == grid%nGridPoints)
  assert(size(f, 1) == grid%nGridPoints)
  assert(size(g, 1) == grid%nGridPoints)
  assert(size(f, 2) == size(g, 2))
#ifdef DEBUG
  if (present(weight)) then
     assert(size(weight) == grid%nGridPoints)
  end if
  if (this%nPatchPoints > 0) then
     assert(allocated(this%norm))
     assert(size(this%norm) == this%nPatchPoints)
  end if
#endif
  assert(this%comm /= MPI_COMM_NULL)

  innerProduct = 0.0_wp

  do l = 1, size(f, 2)
     do k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
        do j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
           do i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
              gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                   &
                   (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                     &
                   (k - 1 - this%gridOffset(3)))
              if (grid%iblank(gridIndex) == 0) cycle
              patchIndex = i - this%offset(1) + this%patchSize(1) *                          &
                   (j - 1 - this%offset(2) + this%patchSize(2) *                             &
                   (k - 1 - this%offset(3)))

              if (present(weight)) then
                 innerProduct = innerProduct + f(gridIndex, l) *                             &
                      this%norm(patchIndex, 1) * g(gridIndex, l) *                           &
                      weight(gridIndex)
              else
                 innerProduct = innerProduct + f(gridIndex, l) *                             &
                      this%norm(patchIndex, 1) * g(gridIndex, l)
              end if

           end do !... i = this%offset(1) + 1, this%offset(1) + this%patchSize(1)
        end do !... j = this%offset(2) + 1, this%offset(2) + this%patchSize(2)
     end do !... k = this%offset(3) + 1, this%offset(3) + this%patchSize(3)
  end do !... l = 1, size(f, 2)

#ifdef SCALAR_TYPE_IS_binary128_IEEE754
  assert(allocated(this%mpiReduceBuffer))
  call MPI_Allgather(innerProduct, 1, SCALAR_TYPE_MPI, this%mpiReduceBuffer,                 &
       1, SCALAR_TYPE_MPI, this%comm, ierror)
  innerProduct = sum(this%mpiReduceBuffer)
#else
  call MPI_Allreduce(MPI_IN_PLACE, innerProduct, 1, SCALAR_TYPE_MPI,                         &
       MPI_SUM, this%comm, ierror)
#endif

end function computeVectorInnerProductOnPatch
