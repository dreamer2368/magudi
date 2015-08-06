#include "config.h"

module Patch_mod

#ifndef NDEBUG
  use ErrorHandler, only : assertImpl
#endif

  use MPI, only : MPI_COMM_NULL, MPI_DATATYPE_NULL

  implicit none
  private

  type, abstract, public :: t_Patch

     character(len = 20) :: name
     character(len = STRING_LENGTH) :: patchType
     integer :: index, gridIndex, normalDirection, iMin, iMax, jMin, jMax, kMin, kMax,       &
          globalSize(3), localSize(3), offset(3), nPatchPoints = 0,                          &
          gridLocalSize(3), gridOffset(3), comm = MPI_COMM_NULL,                             &
          mpiScalarSubarrayType = MPI_DATATYPE_NULL
     integer, allocatable :: mpiAllScalarSubarrayTypes(:)
#ifdef SCALAR_TYPE_IS_binary128_IEEE754
     real(SCALAR_KIND), allocatable :: mpiReduceBuffer(:)
#endif
     real(SCALAR_KIND), allocatable :: norm(:,:)

   contains

     procedure, non_overridable, pass :: setupBase
     procedure, non_overridable, pass :: cleanupBase
     generic :: collect => collectScalar, collectVector, collectTensor
     generic :: disperse => disperseScalar, disperseVector, disperseTensor
     generic :: collectMultiply => collectMultiplyScalar, collectMultiplyVector,             &
          collectMultiplyTensor
     generic :: disperseAdd => disperseAddScalar, disperseAddVector, disperseAddTensor
     generic :: gather => gatherScalar, gatherVector, gatherTensor
     generic :: scatter => scatterScalar, scatterVector, scatterTensor
     generic :: computeInnerProduct => computeScalarInnerProduct, computeVectorInnerProduct

     procedure(setup), pass, deferred :: setup
     procedure(cleanup), pass, deferred :: cleanup
     procedure(updateRhs), pass, deferred :: updateRhs

     procedure, private, pass :: collectScalar
     procedure, private, pass :: collectVector
     procedure, private, pass :: collectTensor

     procedure, private, pass :: disperseScalar
     procedure, private, pass :: disperseVector
     procedure, private, pass :: disperseTensor

     procedure, private, pass :: collectMultiplyScalar
     procedure, private, pass :: collectMultiplyVector
     procedure, private, pass :: collectMultiplyTensor

     procedure, private, pass :: disperseAddScalar
     procedure, private, pass :: disperseAddVector
     procedure, private, pass :: disperseAddTensor

     procedure, private, pass :: gatherScalar
     procedure, private, pass :: gatherVector
     procedure, private, pass :: gatherTensor

     procedure, private, pass :: scatterScalar
     procedure, private, pass :: scatterVector
     procedure, private, pass :: scatterTensor

     procedure, private, pass :: computeScalarInnerProduct
     procedure, private, pass :: computeVectorInnerProduct

  end type t_Patch

  abstract interface

     subroutine setup(this, name, comm, grid, state, extent,                                 &
          normalDirection, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Patch

       class(t_Patch) :: this
       character(len = *), intent(in) :: name
       integer, intent(in) :: comm
       class(t_Grid), intent(in) :: grid
       class(t_State), intent(in) :: state
       integer, intent(in) :: extent(6), normalDirection
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setup

  end interface

  abstract interface

     subroutine cleanup(this)

       import :: t_Patch

       class(t_Patch) :: this

     end subroutine cleanup

  end interface

  abstract interface

     subroutine updateRhs(this, mode, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_Patch

       class(t_Patch) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine updateRhs

  end interface

contains

  subroutine setupBase(this, name, comm, grid, extent, normalDirection)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use SolverOptions_mod, only : t_SolverOptions
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Internal modules >>>
    use InputHelper, only : getOption
    use ErrorHandler, only : gracefulExit

    implicit none

    ! <<< Arguments >>>
    class(t_Patch) :: this
    character(len = *), intent(in) :: name
    integer, intent(in) :: comm
    class(t_Grid), intent(in) :: grid
    integer, intent(in) :: extent(6), normalDirection

    ! <<< Local variables >>>
    logical :: isExtentValid
    character(len = STRING_LENGTH) :: message
    integer :: i, extent_(6), offset(3), procRank, numProcs, ierror
    integer, allocatable :: allLocalSizes(:,:), allOffsets(:,:)

    assert_key(grid%nDimensions, (1, 2, 3))

    assert(all(grid%offset >= 0))
    assert(all(grid%localSize > 0))

    assert(abs(normalDirection) <= grid%nDimensions)
    assert(len_trim(name) > 0)
    assert(all(extent /= 0))

    call this%cleanupBase()

    this%name = trim(name)
    if (comm /= MPI_COMM_NULL) call MPI_Comm_dup(comm, this%comm, ierror)
    this%gridIndex = grid%index

    extent_ = extent

    ! Negative values indicate counting backwards from the end.
    do i = 1, grid%nDimensions
       if (extent_(1+2*(i-1)) < 0)                                                           &
            extent_(1+2*(i-1)) = extent_(1+2*(i-1)) + grid%globalSize(i) + 1
       if (extent_(2+2*(i-1)) < 0)                                                           &
            extent_(2+2*(i-1)) = extent_(2+2*(i-1)) + grid%globalSize(i) + 1
    end do

    isExtentValid = .true.

    ! Check that extent_ describes a part of the grid.
    do i = 1, grid%nDimensions
       isExtentValid = isExtentValid .and. (extent_(2+2*(i-1)) >= extent_(1+2*(i-1)))
       isExtentValid = isExtentValid .and.                                                   &
            (extent_(2+2*(i-1)) - extent_(1+2*(i-1)) + 1 <= grid%globalSize(i))
    end do
    do while (i <= 3) !... reset for direction > number of dimensions.
       extent_(1+2*(i-1)) = 1
       extent_(2+2*(i-1)) = 1
       i = i + 1
    end do

    ! Fail if the extent is invalid.
    if (.not. isExtentValid) then
       write(message, '(3A,I0.0,A,6(1X,I0.0),A)') "Patch '", trim(this%name), "' on grid ",  &
            this%gridIndex, " has an invalid extent:", extent, "!"
       call gracefulExit(grid%comm, message)
    end if

    this%iMin = extent_(1)
    this%iMax = extent_(2)
    this%jMin = extent_(3)
    this%jMax = extent_(4)
    this%kMin = extent_(5)
    this%kMax = extent_(6)

    this%globalSize(1) = this%iMax - this%iMin + 1
    this%globalSize(2) = this%jMax - this%jMin + 1
    this%globalSize(3) = this%kMax - this%kMin + 1

    if (this%comm /= MPI_COMM_NULL) then

       ! Zero-based index of first point on the patch belonging to the ``current'' process
       ! (this value has no meaning if the patch lies outside the part of the grid belonging
       ! to the ``current'' process).
       this%offset(1) = max(this%iMin, grid%offset(1) + 1)
       this%offset(2) = max(this%jMin, grid%offset(2) + 1)
       this%offset(3) = max(this%kMin, grid%offset(3) + 1)
       this%offset = min(this%offset, grid%offset + grid%localSize) - 1

       ! Extent of the patch belonging to the ``current'' process (this value has no meaning
       ! if the patch lies outside the part of the grid belonging to the ``current'' process).
       this%localSize(1) = max(this%iMax, grid%offset(1) + 1)
       this%localSize(2) = max(this%jMax, grid%offset(2) + 1)
       this%localSize(3) = max(this%kMax, grid%offset(3) + 1)
       this%localSize = min(this%localSize, grid%offset + grid%localSize)
       this%localSize = this%localSize - this%offset

    else

       ! Reset size and offset if the patch lies outside the part of the grid belonging to the
       ! ``current'' process.
       this%offset = 0
       this%localSize = 0

    end if

    this%nPatchPoints = product(this%localSize)
    this%gridIndex = grid%index
    this%normalDirection = normalDirection

    this%gridLocalSize = grid%localSize
    this%gridOffset = grid%offset

#ifdef SCALAR_TYPE_IS_binary128_IEEE754
    if (this%comm /= MPI_COMM_NULL) then
       call MPI_Comm_size(this%comm, numProcs, ierror)
       allocate(this%mpiReduceBuffer(numProcs))
    end if
#endif

    if (this%comm /= MPI_COMM_NULL) then

       call MPI_Comm_rank(this%comm, procRank, ierror)
       call MPI_Comm_size(this%comm, numProcs, ierror)

       ! Patch derived subarray type.
       offset = this%offset - extent_(1::2) + 1
       call MPI_Type_create_subarray(3, this%globalSize, this%localSize, offset,             &
            MPI_ORDER_FORTRAN, SCALAR_TYPE_MPI, this%mpiScalarSubarrayType, ierror)
       call MPI_Type_commit(this%mpiScalarSubarrayType, ierror)

       ! Gather local size and offset from all patch processes.
       allocate(allLocalSizes(3, numProcs), source = 0)
       allocate(allOffsets(3, numProcs), source = 0)
       call MPI_Gather(this%localSize, 3, MPI_INTEGER, allLocalSizes,                        &
            3, MPI_INTEGER, 0, this%comm, ierror)
       call MPI_Gather(this%offset, 3, MPI_INTEGER, allOffsets,                              &
            3, MPI_INTEGER, 0, this%comm, ierror)

       if (procRank == 0) then

          ! Allocate patch subarray derived types for all processes on patch master process.
          allocate(this%mpiAllScalarSubarrayTypes(numProcs))

          ! Derived types describing subarrays on patches.
          do i = 1, numProcs
             allOffsets(:,i) = allOffsets(:,i) - extent_(1::2) + 1
             call MPI_Type_create_subarray(3, this%globalSize, allLocalSizes(:,i),           &
                  allOffsets(:,i), MPI_ORDER_FORTRAN, SCALAR_TYPE_MPI,                       &
                  this%mpiAllScalarSubarrayTypes(i), ierror)
             call MPI_Type_commit(this%mpiAllScalarSubarrayTypes(i), ierror)
          end do

       end if

       SAFE_DEALLOCATE(allOffsets)
       SAFE_DEALLOCATE(allLocalSizes)

    end if

    if (this%nPatchPoints > 0 .and. allocated(grid%norm)) then
       allocate(this%norm(this%nPatchPoints, 1))
       call this%collect(grid%norm, this%norm)
    end if

  end subroutine setupBase

  subroutine cleanupBase(this)

    ! <<< External modules >>>
    use MPI

    implicit none

    ! <<< Arguments >>>
    class(t_Patch) :: this

    ! <<< Local variables >>>
    integer :: i, ierror

    if (allocated(this%mpiAllScalarSubarrayTypes)) then
       do i = 1, size(this%mpiAllScalarSubarrayTypes)
          if (this%mpiAllScalarSubarrayTypes(i) /= MPI_DATATYPE_NULL)                        &
               call MPI_Type_free(this%mpiAllScalarSubarrayTypes(i), ierror)
       end do
    end if
    SAFE_DEALLOCATE(this%mpiAllScalarSubarrayTypes)

    SAFE_DEALLOCATE(this%norm)

    if (this%mpiScalarSubarrayType /= MPI_DATATYPE_NULL)                                     &
         call MPI_Type_free(this%mpiScalarSubarrayType, ierror)
    this%mpiScalarSubarrayType = MPI_DATATYPE_NULL

    if (this%comm /= MPI_COMM_NULL .and. this%comm /= MPI_COMM_WORLD)                        &
         call MPI_Comm_free(this%comm, ierror)
    this%comm = MPI_COMM_NULL

  end subroutine cleanupBase

  subroutine collectScalar(this, gridArray, patchArray)

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in) :: gridArray(:)
    real(SCALAR_KIND), intent(out) :: patchArray(:)

    ! <<< Local variables >>>
    integer :: i, j, k, patchIndex, localIndex

    call startTiming("collect")

    do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
       do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
          do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
             patchIndex = i - this%offset(1) +                                               &
                  this%localSize(1) * (j - 1 - this%offset(2) +                              &
                  this%localSize(2) * (k - 1 - this%offset(3)))
             localIndex = i - this%gridOffset(1) +                                           &
                  this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                      &
                  this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
             patchArray(patchIndex) = gridArray(localIndex)
          end do
       end do
    end do

    call endTiming("collect")

  end subroutine collectScalar

  subroutine collectVector(this, gridArray, patchArray)

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in) :: gridArray(:,:)
    real(SCALAR_KIND), intent(out) :: patchArray(:,:)

    ! <<< Local variables >>>
    integer :: i, j, k, l, patchIndex, localIndex

    call startTiming("collect")

    do l = 1, size(patchArray, 2)
       do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
          do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
             do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                patchIndex = i - this%offset(1) +                                            &
                     this%localSize(1) * (j - 1 - this%offset(2) +                           &
                     this%localSize(2) * (k - 1 - this%offset(3)))
                localIndex = i - this%gridOffset(1) +                                        &
                     this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                   &
                     this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
                if (localIndex > size(gridArray, 1)) then
                   print *, "Indices: ", i, j, k
                   print *, "Offset: ", this%offset
                   print *, "Local size: ", this%localSize
                   print *, "Grid offset: ", this%gridOffset
                   print *, "Grid local size: ", this%gridLocalSize
                end if
                patchArray(patchIndex,l) = gridArray(localIndex,l)
             end do
          end do
       end do
    end do

    call endTiming("collect")

  end subroutine collectVector

  subroutine collectTensor(this, gridArray, patchArray)

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in) :: gridArray(:,:,:)
    real(SCALAR_KIND), intent(out) :: patchArray(:,:,:)

    ! <<< Local variables >>>
    integer :: i, j, k, l, m, patchIndex, localIndex

    call startTiming("collect")

    do m = 1, size(patchArray, 3)
       do l = 1, size(patchArray, 2)
          do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
             do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
                do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                   patchIndex = i - this%offset(1) +                                         &
                        this%localSize(1) * (j - 1 - this%offset(2) +                        &
                        this%localSize(2) * (k - 1 - this%offset(3)))
                   localIndex = i - this%gridOffset(1) +                                     &
                        this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                &
                        this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
                   patchArray(patchIndex,l,m) = gridArray(localIndex,l,m)
                end do
             end do
          end do
       end do
    end do

    call endTiming("collect")

  end subroutine collectTensor

  subroutine collectMultiplyScalar(this, gridArray, patchArray)

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in) :: gridArray(:)
    real(SCALAR_KIND), intent(out) :: patchArray(:)

    ! <<< Local variables >>>
    integer :: i, j, k, patchIndex, localIndex

    call startTiming("collectMultiply")

    do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
       do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
          do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
             patchIndex = i - this%offset(1) +                                               &
                  this%localSize(1) * (j - 1 - this%offset(2) +                              &
                  this%localSize(2) * (k - 1 - this%offset(3)))
             localIndex = i - this%gridOffset(1) +                                           &
                  this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                      &
                  this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
             patchArray(patchIndex) = patchArray(patchIndex) * gridArray(localIndex)
          end do
       end do
    end do

    call endTiming("collectMultiply")

  end subroutine collectMultiplyScalar

  subroutine collectMultiplyVector(this, gridArray, patchArray)

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in) :: gridArray(:,:)
    real(SCALAR_KIND), intent(out) :: patchArray(:,:)

    ! <<< Local variables >>>
    integer :: i, j, k, l, patchIndex, localIndex

    call startTiming("collect")

    do l = 1, size(patchArray, 2)
       do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
          do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
             do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                patchIndex = i - this%offset(1) +                                            &
                     this%localSize(1) * (j - 1 - this%offset(2) +                           &
                     this%localSize(2) * (k - 1 - this%offset(3)))
                localIndex = i - this%gridOffset(1) +                                        &
                     this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                   &
                     this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
                patchArray(patchIndex,l) = patchArray(patchIndex,l) * gridArray(localIndex,l)
             end do
          end do
       end do
    end do

    call endTiming("collectMultiply")

  end subroutine collectMultiplyVector

  subroutine collectMultiplyTensor(this, gridArray, patchArray)

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in) :: gridArray(:,:,:)
    real(SCALAR_KIND), intent(out) :: patchArray(:,:,:)

    ! <<< Local variables >>>
    integer :: i, j, k, l, m, patchIndex, localIndex

    call startTiming("collectMultiply")

    do m = 1, size(patchArray, 3)
       do l = 1, size(patchArray, 2)
          do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
             do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
                do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                   patchIndex = i - this%offset(1) +                                         &
                        this%localSize(1) * (j - 1 - this%offset(2) +                        &
                        this%localSize(2) * (k - 1 - this%offset(3)))
                   localIndex = i - this%gridOffset(1) +                                     &
                        this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                &
                        this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
                   patchArray(patchIndex,l,m) = patchArray(patchIndex,l,m) *                 &
                        gridArray(localIndex,l,m)
                end do
             end do
          end do
       end do
    end do

    call endTiming("collectMultiply")

  end subroutine collectMultiplyTensor

  subroutine disperseScalar(this, patchArray, gridArray)

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in) :: patchArray(:)
    real(SCALAR_KIND), intent(out) :: gridArray(:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, patchIndex, localIndex

    call startTiming("disperse")

    gridArray = 0.0_wp

    do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
       do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
          do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
             patchIndex = i - this%offset(1) +                                               &
                  this%localSize(1) * (j - 1 - this%offset(2) +                              &
                  this%localSize(2) * (k - 1 - this%offset(3)))
             localIndex = i - this%gridOffset(1) +                                           &
                  this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                      &
                  this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
             gridArray(localIndex) = patchArray(patchIndex)
          end do
       end do
    end do

    call endTiming("disperse")

  end subroutine disperseScalar

  subroutine disperseVector(this, patchArray, gridArray)

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in) :: patchArray(:,:)
    real(SCALAR_KIND), intent(out) :: gridArray(:,:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, l, patchIndex, localIndex

    call startTiming("disperse")

    gridArray = 0.0_wp

    do l = 1, size(patchArray, 2)
       do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
          do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
             do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                patchIndex = i - this%offset(1) +                                            &
                     this%localSize(1) * (j - 1 - this%offset(2) +                           &
                     this%localSize(2) * (k - 1 - this%offset(3)))
                localIndex = i - this%gridOffset(1) +                                        &
                     this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                   &
                     this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
                gridArray(localIndex,l) = patchArray(patchIndex,l)
             end do
          end do
       end do
    end do

    call endTiming("disperse")

  end subroutine disperseVector

  subroutine disperseTensor(this, patchArray, gridArray)

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in) :: patchArray(:,:,:)
    real(SCALAR_KIND), intent(out) :: gridArray(:,:,:)

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, l, m, patchIndex, localIndex

    call startTiming("disperse")

    gridArray = 0.0_wp

    do m = 1, size(patchArray, 3)
       do l = 1, size(patchArray, 2)
          do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
             do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
                do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                   patchIndex = i - this%offset(1) +                                         &
                        this%localSize(1) * (j - 1 - this%offset(2) +                        &
                        this%localSize(2) * (k - 1 - this%offset(3)))
                   localIndex = i - this%gridOffset(1) +                                     &
                        this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                &
                        this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
                   gridArray(localIndex,l,m) = patchArray(patchIndex,l,m)
                end do
             end do
          end do
       end do
    end do

    call endTiming("disperse")

  end subroutine disperseTensor

  subroutine disperseAddScalar(this, patchArray, gridArray)

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in) :: patchArray(:)
    real(SCALAR_KIND), intent(inout) :: gridArray(:)

    ! <<< Local variables >>>
    integer :: i, j, k, patchIndex, localIndex

    call startTiming("disperseAdd")

    do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
       do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
          do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
             patchIndex = i - this%offset(1) +                                               &
                  this%localSize(1) * (j - 1 - this%offset(2) +                              &
                  this%localSize(2) * (k - 1 - this%offset(3)))
             localIndex = i - this%gridOffset(1) +                                           &
                  this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                      &
                  this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
             gridArray(localIndex) = gridArray(localIndex) + patchArray(patchIndex)
          end do
       end do
    end do

    call endTiming("disperseAdd")

  end subroutine disperseAddScalar

  subroutine disperseAddVector(this, patchArray, gridArray)

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in) :: patchArray(:,:)
    real(SCALAR_KIND), intent(inout) :: gridArray(:,:)

    ! <<< Local variables >>>
    integer :: i, j, k, l, patchIndex, localIndex

    call startTiming("disperseAdd")

    do l = 1, size(patchArray, 2)
       do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
          do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
             do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                patchIndex = i - this%offset(1) +                                            &
                     this%localSize(1) * (j - 1 - this%offset(2) +                           &
                     this%localSize(2) * (k - 1 - this%offset(3)))
                localIndex = i - this%gridOffset(1) +                                        &
                     this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                   &
                     this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
                gridArray(localIndex,l) = gridArray(localIndex,l) + patchArray(patchIndex,l)
             end do
          end do
       end do
    end do

    call endTiming("disperseAdd")

  end subroutine disperseAddVector

  subroutine disperseAddTensor(this, patchArray, gridArray)

    ! <<< Internal modules >>>
    use MPITimingsHelper, only : startTiming, endTiming

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in) :: patchArray(:,:,:)
    real(SCALAR_KIND), intent(inout) :: gridArray(:,:,:)

    ! <<< Local variables >>>
    integer :: i, j, k, l, m, patchIndex, localIndex

    call startTiming("disperseAdd")

    do m = 1, size(patchArray, 3)
       do l = 1, size(patchArray, 2)
          do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
             do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
                do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                   patchIndex = i - this%offset(1) +                                         &
                        this%localSize(1) * (j - 1 - this%offset(2) +                        &
                        this%localSize(2) * (k - 1 - this%offset(3)))
                   localIndex = i - this%gridOffset(1) +                                     &
                        this%gridLocalSize(1) * (j - 1 - this%gridOffset(2) +                &
                        this%gridLocalSize(2) * (k - 1 - this%gridOffset(3)))
                   gridArray(localIndex,l,m) = gridArray(localIndex,l,m) +                   &
                        patchArray(patchIndex,l,m)
                end do
             end do
          end do
       end do
    end do

    call endTiming("disperseAdd")

  end subroutine disperseAddTensor

  subroutine gatherScalar(this, localArray, globalArray)

    ! <<< External modules >>>
    use MPI

    implicit none

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in) :: localArray(:)
    real(SCALAR_KIND), allocatable :: globalArray(:)

    ! <<< Local variables >>>
    integer :: i, mpiRequest, procRank, numProcs, ierror

    if (this%comm == MPI_COMM_NULL) return

    assert(this%nPatchPoints > 0)
    assert(size(localArray) == this%nPatchPoints)

    call MPI_Comm_rank(this%comm, procRank, ierror)
    call MPI_Comm_size(this%comm, numProcs, ierror)

#ifndef NDEBUG
    if (procRank == 0) then
       assert(allocated(globalArray))
       assert(size(globalArray) == product(this%globalSize))
    end if
#endif

    call MPI_Isend(localArray, size(localArray), SCALAR_TYPE_MPI,                            &
         0, procRank, this%comm, mpiRequest, ierror)

    if (procRank == 0) then
       do i = 0, numProcs - 1
          call MPI_Recv(globalArray, 1, this%mpiAllScalarSubarrayTypes(i+1),                 &
               i, i, this%comm, MPI_STATUS_IGNORE, ierror)
       end do
    end if

    call MPI_Wait(mpiRequest, MPI_STATUS_IGNORE, ierror)

  end subroutine gatherScalar

  subroutine gatherVector(this, localArray, globalArray)

    ! <<< External modules >>>
    use MPI

    implicit none

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in) :: localArray(:,:)
    real(SCALAR_KIND), allocatable :: globalArray(:,:)

    ! <<< Local variables >>>
    integer :: i, mpiRequest, procRank, numProcs, ierror

    if (this%comm == MPI_COMM_NULL) return

    assert(this%nPatchPoints > 0)
    assert(size(localArray, 1) == this%nPatchPoints)
    assert(size(localArray, 2) > 0)

    call MPI_Comm_rank(this%comm, procRank, ierror)
    call MPI_Comm_size(this%comm, numProcs, ierror)

#ifndef NDEBUG
    if (procRank == 0) then
       assert(allocated(globalArray))
       assert(size(globalArray, 1) == product(this%globalSize))
       assert(size(globalArray, 2) == size(localArray, 2))
    end if
#endif

    call MPI_Isend(localArray, size(localArray), SCALAR_TYPE_MPI, 0, procRank,               &
         this%comm, mpiRequest, ierror)

    if (procRank == 0) then
       do i = 0, numProcs - 1
          call MPI_Recv(globalArray, size(globalArray, 2),                                   &
               this%mpiAllScalarSubarrayTypes(i+1), i, i,                                    &
               this%comm, MPI_STATUS_IGNORE, ierror)
       end do
    end if

    call MPI_Wait(mpiRequest, MPI_STATUS_IGNORE, ierror)

  end subroutine gatherVector

  subroutine gatherTensor(this, localArray, globalArray)

    ! <<< External modules >>>
    use MPI

    implicit none

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in) :: localArray(:,:,:)
    real(SCALAR_KIND), allocatable :: globalArray(:,:,:)

    ! <<< Local variables >>>
    integer :: i, nComponents, mpiRequest, procRank, numProcs, ierror

    if (this%comm == MPI_COMM_NULL) return

    assert(this%nPatchPoints > 0)
    assert(size(localArray, 1) == this%nPatchPoints)
    assert(size(localArray, 2) > 0)
    assert(size(localArray, 3) > 0)

    call MPI_Comm_rank(this%comm, procRank, ierror)
    call MPI_Comm_size(this%comm, numProcs, ierror)

#ifndef NDEBUG
    if (procRank == 0) then
       assert(allocated(globalArray))
       assert(size(globalArray, 1) == product(this%globalSize))
       assert(size(globalArray, 2) == size(localArray, 2))
       assert(size(globalArray, 3) == size(localArray, 3))
    end if
#endif

    call MPI_Isend(localArray, size(localArray), SCALAR_TYPE_MPI,                            &
         0, procRank, this%comm, mpiRequest, ierror)

    if (procRank == 0) then
       nComponents = size(globalArray, 2) * size(globalArray, 3)
       do i = 0, numProcs - 1
          call MPI_Recv(globalArray, nComponents, this%mpiAllScalarSubarrayTypes(i+1),       &
               i, i, this%comm, MPI_STATUS_IGNORE, ierror)
       end do
    end if

    call MPI_Wait(mpiRequest, MPI_STATUS_IGNORE, ierror)

  end subroutine gatherTensor

  subroutine scatterScalar(this, globalArray, localArray)

    ! <<< External modules >>>
    use MPI

    implicit none

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in), allocatable :: globalArray(:)
    real(SCALAR_KIND), intent(out) :: localArray(:)

    ! <<< Local variables >>>
    integer :: i, mpiRequest, procRank, numProcs, ierror

    if (this%comm == MPI_COMM_NULL) return

    assert(this%nPatchPoints > 0)
    assert(size(localArray) == this%nPatchPoints)

    call MPI_Comm_rank(this%comm, procRank, ierror)
    call MPI_Comm_size(this%comm, numProcs, ierror)

#ifndef NDEBUG
    if (procRank == 0) then
       assert(allocated(globalArray))
       assert(size(globalArray) == product(this%globalSize))
    end if
#endif

    call MPI_Irecv(localArray, size(localArray), SCALAR_TYPE_MPI,                            &
         0, procRank, this%comm, mpiRequest, ierror)

    if (procRank == 0) then
       do i = 0, numProcs - 1
          call MPI_Send(globalArray, 1, this%mpiAllScalarSubarrayTypes(i+1),                 &
               i, i, this%comm, ierror)
       end do
    end if

    call MPI_Wait(mpiRequest, MPI_STATUS_IGNORE, ierror)

  end subroutine scatterScalar

  subroutine scatterVector(this, globalArray, localArray)

    ! <<< External modules >>>
    use MPI

    implicit none

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in), allocatable :: globalArray(:,:)
    real(SCALAR_KIND), intent(out) :: localArray(:,:)

    ! <<< Local variables >>>
    integer :: i, mpiRequest, procRank, numProcs, ierror

    if (this%comm == MPI_COMM_NULL) return

    assert(this%nPatchPoints > 0)
    assert(size(localArray, 1) == this%nPatchPoints)
    assert(size(localArray, 2) > 0)

    call MPI_Comm_rank(this%comm, procRank, ierror)
    call MPI_Comm_size(this%comm, numProcs, ierror)

#ifndef NDEBUG
    if (procRank == 0) then
       assert(allocated(globalArray))
       assert(size(globalArray, 1) == product(this%globalSize))
       assert(size(globalArray, 2) == size(localArray, 2))
    end if
#endif

    call MPI_Irecv(localArray, size(localArray), SCALAR_TYPE_MPI,                            &
         0, procRank, this%comm, mpiRequest, ierror)

    if (procRank == 0) then
       do i = 0, numProcs - 1
          call MPI_Send(globalArray, size(globalArray, 2),                                   &
               this%mpiAllScalarSubarrayTypes(i+1), i, i, this%comm, ierror)
       end do
    end if

    call MPI_Wait(mpiRequest, MPI_STATUS_IGNORE, ierror)

  end subroutine scatterVector

  subroutine scatterTensor(this, globalArray, localArray)

    ! <<< External modules >>>
    use MPI

    implicit none

    ! <<< Arguments >>>
    class(t_Patch) :: this
    real(SCALAR_KIND), intent(in), allocatable :: globalArray(:,:,:)
    real(SCALAR_KIND), intent(out) :: localArray(:,:,:)

    ! <<< Local variables >>>
    integer :: i, nComponents, mpiRequest, procRank, numProcs, ierror

    if (this%comm == MPI_COMM_NULL) return

    assert(this%nPatchPoints > 0)
    assert(size(localArray, 1) == this%nPatchPoints)
    assert(size(localArray, 2) > 0)

    call MPI_Comm_rank(this%comm, procRank, ierror)
    call MPI_Comm_size(this%comm, numProcs, ierror)

#ifndef NDEBUG
    if (procRank == 0) then
       assert(allocated(globalArray))
       assert(size(globalArray, 1) == product(this%globalSize))
       assert(size(globalArray, 2) == size(localArray, 2))
       assert(size(globalArray, 3) == size(localArray, 3))
    end if
#endif

    call MPI_Irecv(localArray, size(localArray), SCALAR_TYPE_MPI,                            &
         0, procRank, this%comm, mpiRequest, ierror)

    if (procRank == 0) then
       nComponents = size(globalArray, 2) * size(globalArray, 3)
       do i = 0, numProcs - 1
          call MPI_Send(globalArray, nComponents, this%mpiAllScalarSubarrayTypes(i+1),       &
               i, i, this%comm, ierror)
       end do
    end if

    call MPI_Wait(mpiRequest, MPI_STATUS_IGNORE, ierror)

  end subroutine scatterTensor

  function computeScalarInnerProduct(this, grid, f, g, weight) result(innerProduct)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid

    implicit none

    ! <<< Arguments >>>
    class(t_Patch) :: this
    class(t_Grid) :: grid
    real(SCALAR_KIND), intent(in) :: f(:), g(:)
    real(SCALAR_KIND), intent(in), optional :: weight(:)

    ! <<< Result >>>
    real(SCALAR_KIND) :: innerProduct

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, gridIndex, patchIndex, ierror

    assert(size(f) == grid%nGridPoints)
    assert(size(g) == grid%nGridPoints)
#ifndef NDEBUG
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

    do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
       do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
          do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
             gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                    &
                  (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                      &
                  (k - 1 - this%gridOffset(3)))
             if (grid%iblank(gridIndex) == 0) cycle
             patchIndex = i - this%offset(1) + this%localSize(1) *                           &
                  (j - 1 - this%offset(2) + this%localSize(2) *                              &
                  (k - 1 - this%offset(3)))

             if (present(weight)) then
                innerProduct = innerProduct + f(gridIndex) *                                 &
                     this%norm(patchIndex, 1) * g(gridIndex) * weight(gridIndex)
             else
                innerProduct = innerProduct + f(gridIndex) *                                 &
                     this%norm(patchIndex, 1) * g(gridIndex)
             end if

          end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
       end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
    end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)

#ifdef SCALAR_TYPE_IS_binary128_IEEE754
    assert(allocated(grid%mpiReduceBuffer))
    call MPI_Allgather(innerProduct, 1, SCALAR_TYPE_MPI, grid%mpiReduceBuffer,               &
         1, SCALAR_TYPE_MPI, grid%comm, ierror)
    innerProduct = sum(grid%mpiReduceBuffer)
#else
    call MPI_Allreduce(MPI_IN_PLACE, innerProduct, 1, SCALAR_TYPE_MPI,                       &
         MPI_SUM, grid%comm, ierror)
#endif

  end function computeScalarInnerProduct

  function computeVectorInnerProduct(this, grid, f, g, weight) result(innerProduct)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid

    implicit none

    ! <<< Arguments >>>
    class(t_Patch) :: this
    class(t_Grid) :: grid
    real(SCALAR_KIND), intent(in) :: f(:,:), g(:,:)
    real(SCALAR_KIND), intent(in), optional :: weight(:)

    ! <<< Result >>>
    real(SCALAR_KIND) :: innerProduct

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, l, gridIndex, patchIndex, ierror

    assert(size(f, 1) == grid%nGridPoints)
    assert(size(g, 1) == grid%nGridPoints)
    assert(size(f, 2) == size(g, 2))
#ifndef NDEBUG
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
       do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
          do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
             do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
                gridIndex = i - this%gridOffset(1) + this%gridLocalSize(1) *                 &
                     (j - 1 - this%gridOffset(2) + this%gridLocalSize(2) *                   &
                     (k - 1 - this%gridOffset(3)))
                if (grid%iblank(gridIndex) == 0) cycle
                patchIndex = i - this%offset(1) + this%localSize(1) *                        &
                     (j - 1 - this%offset(2) + this%localSize(2) *                           &
                     (k - 1 - this%offset(3)))

                if (present(weight)) then
                   innerProduct = innerProduct + f(gridIndex, l) *                           &
                        this%norm(patchIndex, 1) * g(gridIndex, l) *                         &
                        weight(gridIndex)
                else
                   innerProduct = innerProduct + f(gridIndex, l) *                           &
                        this%norm(patchIndex, 1) * g(gridIndex, l)
                end if

             end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
          end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
       end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
    end do !... l = 1, size(f, 2)

#ifdef SCALAR_TYPE_IS_binary128_IEEE754
    assert(allocated(grid%mpiReduceBuffer))
    call MPI_Allgather(innerProduct, 1, SCALAR_TYPE_MPI, grid%mpiReduceBuffer,               &
         1, SCALAR_TYPE_MPI, grid%comm, ierror)
    innerProduct = sum(grid%mpiReduceBuffer)
#else
    call MPI_Allreduce(MPI_IN_PLACE, innerProduct, 1, SCALAR_TYPE_MPI,                       &
         MPI_SUM, grid%comm, ierror)
#endif

  end function computeVectorInnerProduct

end module Patch_mod
