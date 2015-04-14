#include "config.h"

module GridImpl

  implicit none
  public

  integer, parameter, public ::                                                              &
       NONE    = 0,                                                                          &
       PLANE   = 1,                                                                          &
       OVERLAP = 2

contains

  subroutine allocateData(this, simulationFlags)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid
    use SimulationFlags_mod, only : t_SimulationFlags

    ! <<< Arguments >>>
    class(t_Grid) :: this
    type(t_SimulationFlags), intent(in) :: simulationFlags

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
#ifdef SCALAR_TYPE_IS_binary128_IEEE754
    integer :: nProcs
#endif

    allocate(this%firstDerivative(this%nDimensions))
    if (.not. simulationFlags%repeatFirstDerivative)                                         &
         allocate(this%secondDerivative(this%nDimensions))
    if (simulationFlags%dissipationOn) allocate(this%dissipation(this%nDimensions))
    if (.not. simulationFlags%predictionOnly)                                                &
         allocate(this%adjointFirstDerivative(this%nDimensions))

    allocate(this%iblank(this%nGridPoints), source = 1)
    allocate(this%coordinates(this%nGridPoints, this%nDimensions))
    allocate(this%jacobian(this%nGridPoints, 1))
    allocate(this%metrics(this%nGridPoints, this%nDimensions ** 2))
    allocate(this%norm(this%nGridPoints, 1))
    this%norm = 1.0_wp

    if (.not. simulationFlags%predictionOnly) then
       allocate(this%targetMollifier(this%nGridPoints, 1))
       allocate(this%controlMollifier(this%nGridPoints, 1))
    end if

#ifdef SCALAR_TYPE_IS_binary128_IEEE754
    call MPI_Comm_size(this%comm, nProcs, ierror)
    allocate(this%mpiReduceBuffer(nProcs))
#endif

  end subroutine allocateData

  subroutine makeUnitCube(this)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use Grid_mod, only : t_Grid

    ! <<< Arguments >>>
    class(t_Grid) :: this

    ! <<< Local variables >>>
    integer, parameter :: wp = SCALAR_KIND
    integer :: i, j, k, l
    real(wp) :: h
    real(wp), allocatable :: unitInterval(:)

    assert_key(this%nDimensions, (1, 2, 3))

    do l = 1, this%nDimensions

       SAFE_DEALLOCATE(unitInterval)
       allocate(unitInterval(this%globalSize(l)))

       if (this%periodicityType(l) == PLANE) then
          h = 1.0_wp / real(this%globalSize(l), wp)
       else
          h = 1.0_wp / real(this%globalSize(l) - 1, wp)
       end if

       do i = 1, this%globalSize(l)
          unitInterval(i) = real(i - 1, wp) * h
       end do

       select case (l)

       case (1)
          do k = 1, this%localSize(3)
             do j = 1, this%localSize(2)
                do i = 1, this%localSize(1)
                   this%coordinates(i + this%localSize(1) * (j - 1 +                         &
                        this%localSize(2) * (k - 1)), l) = unitInterval(i + this%offset(l))
                end do
             end do
          end do

       case (2)
          do k = 1, this%localSize(3)
             do j = 1, this%localSize(2)
                do i = 1, this%localSize(1)
                   this%coordinates(i + this%localSize(1) * (j - 1 +                         &
                        this%localSize(2) * (k - 1)), l) = unitInterval(j + this%offset(l))
                end do
             end do
          end do

       case (3)
          do k = 1, this%localSize(3)
             do j = 1, this%localSize(2)
                do i = 1, this%localSize(1)
                   this%coordinates(i + this%localSize(1) * (j - 1 +                         &
                        this%localSize(2) * (k - 1)), l) = unitInterval(k + this%offset(l))
                end do
             end do
          end do

       end select

       SAFE_DEALLOCATE(unitInterval)

    end do !... l = 1, this%nDimensions

  end subroutine makeUnitCube

end module GridImpl

subroutine setupGrid(this, index, globalSize, comm, processDistribution,                     &
     periodicityType, periodicLength, simulationFlags)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Private members >>>
  use GridImpl

  ! <<< Internal modules >>>
  use MPIHelper, only : pigeonhole
  use InputHelper, only : getOption, getRequiredOption

  implicit none

  ! <<< Arguments >>>
  class(t_Grid) :: this
  integer, intent(in) :: index, globalSize(:)
  integer, intent(in), optional :: comm, processDistribution(:), periodicityType(:)
  real(SCALAR_KIND), intent(in), optional :: periodicLength(:)
  type(t_SimulationFlags), optional :: simulationFlags

  ! <<< Local variables >>>
  type(t_SimulationFlags) :: simulationFlags_
  integer :: i, comm_, procRank, nProcs, ierror
  character(len = STRING_LENGTH) :: key, val
  integer, allocatable :: processDistribution_(:), processCoordinates(:)
  logical :: isPeriodic(3)

  assert(size(globalSize) > 0 .and. size(globalSize) <= 3)
  assert(index > 0)
  assert(all(globalSize > 0))

  if (present(simulationFlags)) then
     simulationFlags_ = simulationFlags
  else
     call simulationFlags_%initialize()
  end if

  ! Clean slate.
  call this%cleanup()
  this%index = index
  this%globalSize = 1
  this%globalSize(1:size(globalSize)) = globalSize

  ! Get periodicity information from input, if not specified:

  this%periodicityType = NONE

  if (.not. present(periodicityType)) then

     do i = 1, size(globalSize)
        write(key, '(A,I3.3,A,I1.1,A)') "grid", this%index, "/dir", i, "/"
        val = getOption(trim(key) // "periodicity_type", "")
        this%periodicityType(i) = NONE
        if (trim(val) == "PLANE") then
           this%periodicityType(i) = PLANE
        else if (trim(val) == "OVERLAP") then
           this%periodicityType(i) = OVERLAP
        end if
        if (this%periodicityType(i) == PLANE)                                                &
             call getRequiredOption(trim(key) // "periodic_length",                          &
             this%periodicLength(i), this%comm)
     end do !... i =  1, size(globalSize)

  else

     assert(size(periodicityType) == size(globalSize))
     this%periodicityType(1:size(periodicityType)) = periodicityType
     assert_key(this%periodicityType(1), (NONE, PLANE, OVERLAP))
     assert_key(this%periodicityType(2), (NONE, PLANE, OVERLAP))
     assert_key(this%periodicityType(3), (NONE, PLANE, OVERLAP))

     if (any(periodicityType == PLANE)) then

        assert(present(periodicLength))
        assert(size(periodicLength) == size(globalSize))

        where (this%periodicityType(1:size(periodicityType)) == PLANE)
           this%periodicLength(1:size(periodicityType)) = periodicLength
        end where

     end if

  end if

  ! Find the number of processes in the grid communicator if specified. Otherwise,
  ! use `MPI_COMM_WORLD`.
  comm_ = MPI_COMM_WORLD
  if (present(comm)) comm_ = comm
  call MPI_Comm_size(comm_, nProcs, ierror)

  ! Generate a default process distribution. If one was specified, override the default.
  if (present(processDistribution)) then
     assert(size(processDistribution) == size(globalSize))
     assert(all(processDistribution >= 0))
     allocate(processDistribution_(size(globalSize)), source = processDistribution)
  else
     allocate(processDistribution_(size(globalSize)), source = 0)
  end if
  call MPI_Dims_create(nProcs, size(globalSize), processDistribution_, ierror)

  ! Create a Cartesian communicator.
  isPeriodic = (this%periodicityType /= NONE)
  call MPI_Cart_create(comm_, size(globalSize), processDistribution_,                        &
       isPeriodic(1:size(globalSize)), .true., this%comm, ierror)
  call MPI_Cartdim_get(this%comm, this%nDimensions, ierror)

  ! Find process coordinates in Cartesian topology.
  allocate(processCoordinates(size(globalSize)), source = 0)
  call MPI_Comm_rank(this%comm, procRank, ierror)
  call MPI_Cart_coords(this%comm, procRank, size(globalSize), processCoordinates, ierror)

  ! Find the size of the part of the grid distributed to the current process.
  this%offset = 0
  this%localSize = 1
  do i = 1, size(globalSize)
     call pigeonhole(this%globalSize(i), processDistribution_(i), processCoordinates(i),     &
          this%offset(i), this%localSize(i))
  end do
  this%nGridPoints = product(this%localSize)

  ! Derived types for MPI I/O.
  call MPI_Type_create_subarray(3, this%globalSize, this%localSize, this%offset,             &
       MPI_ORDER_FORTRAN, SCALAR_TYPE_MPI, this%mpiDerivedTypeScalarSubarray, ierror)
  call MPI_Type_commit(this%mpiDerivedTypeScalarSubarray, ierror)
  call MPI_Type_create_subarray(3, this%globalSize, this%localSize, this%offset,             &
       MPI_ORDER_FORTRAN, MPI_INTEGER, this%mpiDerivedTypeIntegerSubarray, ierror)
  call MPI_Type_commit(this%mpiDerivedTypeIntegerSubarray, ierror)

  ! Is the grid curvilinear/rectangular?
  this%isCurvilinear = simulationFlags_%isDomainCurvilinear
  write(key, '(A,I3.3,A)') "grid", this%index, "/curvilinear"
  this%isCurvilinear = getOption(key, this%isCurvilinear)

  ! Allocate grid data.
  call allocateData(this, simulationFlags_)
  call makeUnitCube(this)

  SAFE_DEALLOCATE(processCoordinates)
  SAFE_DEALLOCATE(processDistribution_)

end subroutine setupGrid

subroutine cleanupGrid(this)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid

  implicit none

  ! <<< Arguments >>>
  class(t_Grid) :: this

  ! <<< Local variables >>>
  integer :: i, ierror

  if (allocated(this%firstDerivative)) then
     do i = 1, size(this%firstDerivative)
        call this%firstDerivative(i)%cleanup()
     end do
  end if
  SAFE_DEALLOCATE(this%firstDerivative)

  if (allocated(this%secondDerivative)) then
     do i = 1, size(this%secondDerivative)
        call this%secondDerivative(i)%cleanup()
     end do
  end if
  SAFE_DEALLOCATE(this%secondDerivative)

  if (allocated(this%dissipation)) then
     do i = 1, size(this%dissipation)
        call this%dissipation(i)%cleanup()
     end do
  end if
  SAFE_DEALLOCATE(this%dissipation)

  if (allocated(this%adjointFirstDerivative)) then
     do i = 1, size(this%adjointFirstDerivative)
        call this%adjointFirstDerivative(i)%cleanup()
     end do
  end if
  SAFE_DEALLOCATE(this%adjointFirstDerivative)

  SAFE_DEALLOCATE(this%iblank)
  SAFE_DEALLOCATE(this%coordinates)
  SAFE_DEALLOCATE(this%jacobian)
  SAFE_DEALLOCATE(this%metrics)
  SAFE_DEALLOCATE(this%norm)

  SAFE_DEALLOCATE(this%targetMollifier)
  SAFE_DEALLOCATE(this%controlMollifier)

#ifdef SCALAR_TYPE_IS_binary128_IEEE754
  SAFE_DEALLOCATE(this%mpiReduceBuffer)
#endif

  if (this%mpiDerivedTypeIntegerSubarray /= MPI_DATATYPE_NULL)                               &
       call MPI_Type_free(this%mpiDerivedTypeIntegerSubarray, ierror)
  if (this%mpiDerivedTypeScalarSubarray /= MPI_DATATYPE_NULL)                                &
       call MPI_Type_free(this%mpiDerivedTypeScalarSubarray, ierror)
  if (this%comm /= MPI_COMM_NULL .and. this%comm /= MPI_COMM_WORLD)                          &
       call MPI_Comm_free(this%comm, ierror)

  this%comm = MPI_COMM_NULL
  this%mpiDerivedTypeScalarSubarray = MPI_DATATYPE_NULL
  this%mpiDerivedTypeIntegerSubarray = MPI_DATATYPE_NULL

  this%nGridPoints = 0

end subroutine cleanupGrid

subroutine loadGridData(this, quantityOfInterest, filename, offsetInBytes, success)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid

  ! <<< Enumerations >>>
  use Grid_enum

  ! <<< Internal modules >>>
  use PLOT3DHelper, only : plot3dReadSingleGrid, plot3dReadSingleFunction

  ! <<< Arguments >>>
  class(t_Grid) :: this
  integer, intent(in) :: quantityOfInterest
  character(len = *), intent(in) :: filename
  integer(kind = MPI_OFFSET_KIND), intent(inout) :: offsetInBytes
  logical, intent(out) :: success

  assert(len_trim(filename) > 0)
  assert(offsetInBytes >= 0)

  select case (quantityOfInterest)
  case (QOI_GRID)
     call plot3dReadSingleGrid(this%comm, trim(filename), offsetInBytes,                     &
          this%mpiDerivedTypeScalarSubarray, this%mpiDerivedTypeIntegerSubarray,             &
          this%globalSize, this%coordinates, this%iblank, success)
  case (QOI_JACOBIAN)
     call plot3dReadSingleFunction(this%comm, trim(filename), offsetInBytes,                 &
          this%mpiDerivedTypeScalarSubarray, this%globalSize, this%jacobian, success)
  case (QOI_METRICS)
     call plot3dReadSingleFunction(this%comm, trim(filename), offsetInBytes,                 &
          this%mpiDerivedTypeScalarSubarray, this%globalSize, this%metrics, success)
  case (QOI_TARGET_MOLLIFIER)
     call plot3dReadSingleFunction(this%comm, trim(filename), offsetInBytes,                 &
          this%mpiDerivedTypeScalarSubarray, this%globalSize,                                &
          this%targetMollifier, success)
  case (QOI_CONTROL_MOLLIFIER)
     call plot3dReadSingleFunction(this%comm, trim(filename), offsetInBytes,                 &
          this%mpiDerivedTypeScalarSubarray, this%globalSize,                                &
          this%controlMollifier, success)
  end select

end subroutine loadGridData

subroutine saveGridData(this, quantityOfInterest, filename, offsetInBytes, success)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid

  ! <<< Enumerations >>>
  use Grid_enum

  ! <<< Internal modules >>>
  use PLOT3DHelper, only : plot3dWriteSingleGrid, plot3dWriteSingleFunction

  ! <<< Arguments >>>
  class(t_Grid) :: this
  integer, intent(in) :: quantityOfInterest
  character(len = *), intent(in) :: filename
  integer(kind = MPI_OFFSET_KIND), intent(inout) :: offsetInBytes
  logical, intent(out) :: success

  assert(len_trim(filename) > 0)
  assert(offsetInBytes >= 0)

  select case (quantityOfInterest)
  case (QOI_GRID)
     call plot3dWriteSingleGrid(this%comm, trim(filename), offsetInBytes,                    &
          this%mpiDerivedTypeScalarSubarray, this%mpiDerivedTypeIntegerSubarray,             &
          this%globalSize, this%coordinates, this%iblank, success)
  case (QOI_JACOBIAN)
     call plot3dWriteSingleFunction(this%comm, trim(filename), offsetInBytes,                &
          this%mpiDerivedTypeScalarSubarray, this%globalSize, this%jacobian, success)
  case (QOI_METRICS)
     call plot3dWriteSingleFunction(this%comm, trim(filename), offsetInBytes,                &
          this%mpiDerivedTypeScalarSubarray, this%globalSize, this%metrics, success)
  case (QOI_TARGET_MOLLIFIER)
     call plot3dWriteSingleFunction(this%comm, trim(filename), offsetInBytes,                &
          this%mpiDerivedTypeScalarSubarray, this%globalSize,                                &
          this%targetMollifier, success)
  case (QOI_CONTROL_MOLLIFIER)
     call plot3dWriteSingleFunction(this%comm, trim(filename), offsetInBytes,                &
          this%mpiDerivedTypeScalarSubarray, this%globalSize,                                &
          this%controlMollifier, success)
  end select

end subroutine saveGridData

subroutine setupSpatialDiscretization(this, simulationFlags, solverOptions)

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid
  use SolverOptions_mod, only : t_SolverOptions
  use SimulationFlags_mod, only : t_SimulationFlags

  ! <<< Private members >>>
  use GridImpl, only : OVERLAP

  ! <<< Internal modules >>>
  use InputHelper, only : getOption
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_Grid) :: this
  type(t_SimulationFlags), intent(in), optional :: simulationFlags
  type(t_SolverOptions), intent(in), optional :: solverOptions

  ! <<< Local variables >>>
  type(t_SimulationFlags) :: simulationFlags_
  type(t_SolverOptions) :: solverOptions_
  integer :: i
  character(len = STRING_LENGTH) :: key, val

  if (present(simulationFlags)) then
     simulationFlags_ = simulationFlags
  else
     call simulationFlags_%initialize()
  end if

  if (present(solverOptions)) then
     solverOptions_ = solverOptions
  else
     call solverOptions_%initialize(this%nDimensions, simulationFlags_, this%comm)
  end if

  do i = 1, this%nDimensions

     write(key, '(A,I3.3,A,I1.1,A)') "grid", this%index, "/dir", i, "/"

     ! First derivative operators.
     if (this%globalSize(i) > 1) then
        val = getOption("defaults/first_derivative_scheme",                                  &
             trim(solverOptions_%discretizationType))
        val = getOption(trim(key) // "first_derivative_scheme", trim(val))
        val = trim(val) // " first derivative"
     else
        val = "null matrix"
     end if
     call this%firstDerivative(i)%setup(trim(val))
     call this%firstDerivative(i)%update(this%comm, i, this%periodicityType(i) == OVERLAP)

     ! Adjoint first derivative operator.
     if (allocated(this%adjointFirstDerivative)) then
        if (simulationFlags_%useContinuousAdjoint .or. trim(val) == "null matrix") then
           call this%adjointFirstDerivative(i)%setup(trim(val))
           this%adjointFirstDerivative(i)%rhsInterior = -this%firstDerivative(i)%rhsInterior
           this%adjointFirstDerivative(i)%rhsBoundary1 = -this%firstDerivative(i)%rhsBoundary1
           this%adjointFirstDerivative(i)%rhsBoundary2 = -this%firstDerivative(i)%rhsBoundary2
        else
           call this%firstDerivative(i)%getAdjoint(this%adjointFirstDerivative(i))
        end if
        call this%adjointFirstDerivative(i)%update(this%comm,                                &
             i, this%periodicityType(i) == OVERLAP)
     end if

     ! Second derivative operators.
     if (allocated(this%secondDerivative)) then
        if (this%globalSize(i) > 1) then
           val = getOption("defaults/second_derivative_scheme",                              &
                trim(solverOptions_%discretizationType))
           val = getOption(trim(key) // "second_derivative_scheme", trim(val))
           val = trim(val) // " second derivative"
        else
           val = "null matrix"
        end if
        call this%secondDerivative(i)%setup(val)
        call this%secondDerivative(i)%update(this%comm, i, this%periodicityType(i) == OVERLAP)
     end if

     ! Artificial dissipation operators.
     if (allocated(this%dissipation)) then
        if (this%globalSize(i) > 1) then
           val = getOption("defaults/artificial_dissipation_scheme",                         &
                trim(solverOptions_%discretizationType))
           val = getOption(trim(key) // "artificial_dissipation_scheme", trim(val))
           val = trim(val) // " dissipation"
        else
           val = "null matrix"
        end if
        call this%dissipation(i)%setup(val)
        call this%dissipation(i)%update(this%comm, i, this%periodicityType(i) == OVERLAP)
     end if

  end do !... i = 1, this%nDimensions

end subroutine setupSpatialDiscretization

subroutine computeCoordinateDerivatives(this, direction, coordinateDerivatives)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid

  ! <<< Private members >>>
  use GridImpl, only : PLANE

  ! <<< Internal modules >>>
  use MPIHelper, only : fillGhostPoints

  implicit none

  ! <<< Arguments >>>
  class(t_Grid) :: this
  integer, intent(in) :: direction
  SCALAR_TYPE, intent(out) :: coordinateDerivatives(:,:)

  ! <<< Local variables >>>
  integer :: i, j, k, l, nDimensions, gridSizeWithGhostPoints(3), numGhostPointsBegin(3),    &
       processDistribution(3), processCoordinates(3), ierror
  logical :: isPeriodic(3)
  SCALAR_TYPE, allocatable :: xWithGhostPoints(:,:,:,:)

  ! Get information from the grid communicator.
  nDimensions = this%nDimensions
  call MPI_Cart_get(this%comm, nDimensions, processDistribution,                             &
       isPeriodic, processCoordinates, ierror)

  ! Straightforward if periodicity type is not `PLANE`.
  if (this%periodicityType(direction) /= PLANE) then
     coordinateDerivatives = this%coordinates
     call this%firstDerivative(direction)%apply(coordinateDerivatives, this%localSize)
     return
  end if

  ! Special hack require for periodicity type `PLANE`:

  ! Find the grid size including ghost points, and the number of ghost points at the
  ! beginning (which is also the offset of the physical grid points).
  gridSizeWithGhostPoints = this%localSize
  gridSizeWithGhostPoints(direction) = gridSizeWithGhostPoints(direction) +                  &
       sum(this%firstDerivative(direction)%nGhost)
  numGhostPointsBegin = 0
  numGhostPointsBegin(direction) = this%firstDerivative(direction)%nGhost(1)

  ! Allocate an array to hold both physical and ghost point coordinates.
  SAFE_DEALLOCATE(xWithGhostPoints)
  allocate(xWithGhostPoints(gridSizeWithGhostPoints(1), gridSizeWithGhostPoints(2),          &
       gridSizeWithGhostPoints(3), nDimensions))

  ! Copy coordinates at physical grid points to the ghost array.
  do l = 1, nDimensions
     do k = 1, this%localSize(3)
        do j = 1, this%localSize(2)
           do i = 1, this%localSize(1)
              xWithGhostPoints(i + numGhostPointsBegin(1),                                   &
                   j + numGhostPointsBegin(2), k + numGhostPointsBegin(3), l) =              &
                   this%coordinates(i + this%localSize(1) * (j - 1 +                         &
                   this%localSize(2) * (k - 1)), l)
           end do
        end do
     end do
  end do

  ! Exchange ghost points between MPI processes.
  call fillGhostPoints(this%comm, xWithGhostPoints, direction,                               &
       this%firstDerivative(direction)%nGhost,                                               &
       this%firstDerivative(direction)%periodicOffset)

  ! At the first process, subtract the periodic length from coordinates received from
  ! the last process
  if (processCoordinates(direction) == 0) then
     select case (direction)
     case (1)
        xWithGhostPoints(1:numGhostPointsBegin(1),:,:,1) =                                   &
             xWithGhostPoints(1:numGhostPointsBegin(1),:,:,1) -                              &
             this%periodicLength(1)
     case (2)
        xWithGhostPoints(:,1:numGhostPointsBegin(2),:,2) =                                   &
             xWithGhostPoints(:,1:numGhostPointsBegin(2),:,2) -                              &
             this%periodicLength(2)
     case (3)
        xWithGhostPoints(:,:,1:numGhostPointsBegin(3),3) =                                   &
             xWithGhostPoints(:,:,1:numGhostPointsBegin(3),3) -                              &
             this%periodicLength(3)
     end select
  end if

  ! At the last process, add the periodic length from coordinates received from the
  ! first process
  if (processCoordinates(direction) == processDistribution(direction) - 1) then
     select case (direction)
     case (1)
        xWithGhostPoints(gridSizeWithGhostPoints(1) + 1 -                                    &
             numGhostPointsBegin(1) : gridSizeWithGhostPoints(1),:,:,1)                      &
             = xWithGhostPoints(gridSizeWithGhostPoints(1) + 1 -                             &
             numGhostPointsBegin(1) : gridSizeWithGhostPoints(1),:,:,1) +                    &
             this%periodicLength(1)
     case (2)
        xWithGhostPoints(:,gridSizeWithGhostPoints(2) + 1 -                                  &
             numGhostPointsBegin(2) : gridSizeWithGhostPoints(2),:,2)                        &
             = xWithGhostPoints(:,gridSizeWithGhostPoints(2) + 1 -                           &
             numGhostPointsBegin(2) : gridSizeWithGhostPoints(2),:,2) +                      &
             this%periodicLength(2)
     case (3)
        xWithGhostPoints(:,:,gridSizeWithGhostPoints(3) + 1 -                                &
             numGhostPointsBegin(3) : gridSizeWithGhostPoints(3),3)                          &
             = xWithGhostPoints(:,:,gridSizeWithGhostPoints(3) + 1 -                         &
             numGhostPointsBegin(3) : gridSizeWithGhostPoints(3),3) +                        &
             this%periodicLength(3)
     end select
  end if

  ! Stencil needs to be applied only at interior points.
  call this%firstDerivative(direction)%applyAtInteriorPoints(xWithGhostPoints,               &
       coordinateDerivatives, this%localSize)

  SAFE_DEALLOCATE(xWithGhostPoints)

end subroutine computeCoordinateDerivatives

subroutine updateGrid(this, hasNegativeJacobian, errorMessage)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid

  ! <<< Private members >>>
  use GridImpl, only : PLANE

  ! <<< Internal modules >>>
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  class(t_Grid) :: this
  logical, intent(out), optional :: hasNegativeJacobian
  character(len = STRING_LENGTH), intent(out), optional :: errorMessage

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, nDimensions, ierror
  SCALAR_TYPE, allocatable :: jacobianMatrixInverse(:,:), F(:,:)
  logical :: hasNegativeJacobian_
  character(len = STRING_LENGTH) :: message
  SCALAR_TYPE :: jacobianOutsideRange

  assert_key(this%nDimensions, (1, 2, 3))
  assert(this%nGridPoints > 0)
  assert(all(this%localSize > 0) .and. product(this%localSize) == this%nGridPoints)

  nDimensions = this%nDimensions

  allocate(jacobianMatrixInverse(this%nGridPoints, nDimensions ** 2))
  do i = 1, nDimensions
     call this%computeCoordinateDerivatives(i,                                               &
          jacobianMatrixInverse(:,(i-1)*nDimensions+1:i*nDimensions))
  end do

  ! Zero out `jacobianMatrixInverse` at hole points.
  assert(allocated(this%iblank))
  assert(size(this%iblank) == this%nGridPoints)
  do i = 1, nDimensions ** 2
     where (this%iblank == 0) jacobianMatrixInverse(:,i) = 0.0_wp
  end do

  assert(allocated(this%metrics))
  assert(size(this%metrics, 1) == this%nGridPoints)
  assert(size(this%metrics, 2) == nDimensions ** 2)

  assert(allocated(this%jacobian))
  assert(size(this%jacobian) == this%nGridPoints)

  assert(allocated(this%firstDerivative))
  assert(size(this%firstDerivative) == nDimensions)

  select case (nDimensions)

  case (1)
     this%jacobian(:,1) = jacobianMatrixInverse(:,1)
     this%metrics(:,1) = 1.0_wp

  case (2)

     if (this%isCurvilinear) then
        this%jacobian(:,1) = (jacobianMatrixInverse(:,1) * jacobianMatrixInverse(:,4) -      &
             jacobianMatrixInverse(:,2) * jacobianMatrixInverse(:,3))
        this%metrics(:,1) = jacobianMatrixInverse(:,4)
        this%metrics(:,2) = - jacobianMatrixInverse(:,3)
        this%metrics(:,3) = - jacobianMatrixInverse(:,2)
        this%metrics(:,4) = jacobianMatrixInverse(:,1)
     else
        this%jacobian(:,1) = jacobianMatrixInverse(:,1) * jacobianMatrixInverse(:,4)
        this%metrics(:,1) = jacobianMatrixInverse(:,4)
        this%metrics(:,2) = 0.0_wp
        this%metrics(:,3) = 0.0_wp
        this%metrics(:,4) = jacobianMatrixInverse(:,1)
     end if

  case (3)

     if (this%isCurvilinear) then
        this%jacobian(:,1) =                                                                 &
             jacobianMatrixInverse(:,1) *                                                    &
             (jacobianMatrixInverse(:,5) * jacobianMatrixInverse(:,9) -                      &
             jacobianMatrixInverse(:,8) * jacobianMatrixInverse(:,6)) +                      &
             jacobianMatrixInverse(:,4) *                                                    &
             (jacobianMatrixInverse(:,8) * jacobianMatrixInverse(:,3) -                      &
             jacobianMatrixInverse(:,2) * jacobianMatrixInverse(:,9)) +                      &
             jacobianMatrixInverse(:,7) *                                                    &
             (jacobianMatrixInverse(:,2) * jacobianMatrixInverse(:,6) -                      &
             jacobianMatrixInverse(:,5) * jacobianMatrixInverse(:,3))
     else
        this%jacobian(:,1) = jacobianMatrixInverse(:,1) *                                    &
             jacobianMatrixInverse(:,5) * jacobianMatrixInverse(:,9)
     end if

     if (any(this%periodicityType == PLANE)) then

        if (this%isCurvilinear) then
           this%metrics(:,1) = jacobianMatrixInverse(:,5) * jacobianMatrixInverse(:,9) -     &
                jacobianMatrixInverse(:,8) * jacobianMatrixInverse(:,6)
           this%metrics(:,4) = jacobianMatrixInverse(:,7) * jacobianMatrixInverse(:,6) -     &
                jacobianMatrixInverse(:,4) * jacobianMatrixInverse(:,9)
           this%metrics(:,7) = jacobianMatrixInverse(:,4) * jacobianMatrixInverse(:,8) -     &
                jacobianMatrixInverse(:,7) * jacobianMatrixInverse(:,5)
           this%metrics(:,2) = jacobianMatrixInverse(:,8) * jacobianMatrixInverse(:,3) -     &
                jacobianMatrixInverse(:,2) * jacobianMatrixInverse(:,9)
           this%metrics(:,5) = jacobianMatrixInverse(:,1) * jacobianMatrixInverse(:,9) -     &
                jacobianMatrixInverse(:,7) * jacobianMatrixInverse(:,3)
           this%metrics(:,8) = jacobianMatrixInverse(:,7) * jacobianMatrixInverse(:,2) -     &
                jacobianMatrixInverse(:,1) * jacobianMatrixInverse(:,8)
           this%metrics(:,3) = jacobianMatrixInverse(:,2) * jacobianMatrixInverse(:,6) -     &
                jacobianMatrixInverse(:,5) * jacobianMatrixInverse(:,3)
           this%metrics(:,6) = jacobianMatrixInverse(:,4) * jacobianMatrixInverse(:,3) -     &
                jacobianMatrixInverse(:,1) * jacobianMatrixInverse(:,6)
           this%metrics(:,9) = jacobianMatrixInverse(:,1) * jacobianMatrixInverse(:,5) -     &
                jacobianMatrixInverse(:,4) * jacobianMatrixInverse(:,2)
        else
           this%metrics(:,1) = jacobianMatrixInverse(:,5) * jacobianMatrixInverse(:,9)
           this%metrics(:,2) = 0.0_wp
           this%metrics(:,3) = 0.0_wp
           this%metrics(:,4) = 0.0_wp
           this%metrics(:,5) = jacobianMatrixInverse(:,1) * jacobianMatrixInverse(:,9)
           this%metrics(:,6) = 0.0_wp
           this%metrics(:,7) = 0.0_wp
           this%metrics(:,8) = 0.0_wp
           this%metrics(:,9) = jacobianMatrixInverse(:,1) * jacobianMatrixInverse(:,5)
        end if

     else

        assert(allocated(this%coordinates))
        assert(size(this%coordinates, 1) == this%nGridPoints)
        assert(size(this%coordinates, 2) == nDimensions)

        allocate(F(this%nGridPoints, 1))

        F(:,1) = jacobianMatrixInverse(:,5) * this%coordinates(:,3)
        call this%firstDerivative(3)%apply(F, this%localSize)
        this%metrics(:,1) = F(:,1)
        if (this%isCurvilinear) then
           F(:,1) = jacobianMatrixInverse(:,8) * this%coordinates(:,3)
           call this%firstDerivative(2)%apply(F, this%localSize)
           this%metrics(:,1) = this%metrics(:,1) - F(:,1)
        end if

        if (.not. this%isCurvilinear) then
           this%metrics(:,2) = 0.0_wp
        else
           F(:,1) = jacobianMatrixInverse(:,8) * this%coordinates(:,3)
           call this%firstDerivative(1)%apply(F, this%localSize)
           this%metrics(:,2) = F(:,1)
           F(:,1) = jacobianMatrixInverse(:,2) * this%coordinates(:,3)
           call this%firstDerivative(3)%apply(F, this%localSize)
           this%metrics(:,2) = this%metrics(:,2) - F(:,1)
        end if

        if (.not. this%isCurvilinear) then
           this%metrics(:,3) = 0.0_wp
        else
           F(:,1) = jacobianMatrixInverse(:,2) * this%coordinates(:,3)
           call this%firstDerivative(2)%apply(F, this%localSize)
           this%metrics(:,3) = F(:,1)
           F(:,1) = jacobianMatrixInverse(:,5) * this%coordinates(:,3)
           call this%firstDerivative(1)%apply(F, this%localSize)
           this%metrics(:,3) = this%metrics(:,3) - F(:,1)
        end if

        if (.not. this%isCurvilinear) then
           this%metrics(:,4) = 0.0_wp
        else
           F(:,1) = jacobianMatrixInverse(:,6) * this%coordinates(:,1)
           call this%firstDerivative(3)%apply(F, this%localSize)
           this%metrics(:,4) = F(:,1)
           F(:,1) = jacobianMatrixInverse(:,9) * this%coordinates(:,1)
           call this%firstDerivative(2)%apply(F, this%localSize)
           this%metrics(:,4) = this%metrics(:,4) - F(:,1)
        end if

        F(:,1) = jacobianMatrixInverse(:,9) * this%coordinates(:,1)
        call this%firstDerivative(1)%apply(F, this%localSize)
        this%metrics(:,5) = F(:,1)
        if (this%isCurvilinear) then
           F(:,1) = jacobianMatrixInverse(:,3) * this%coordinates(:,1)
           call this%firstDerivative(3)%apply(F, this%localSize)
           this%metrics(:,5) = this%metrics(:,5) - F(:,1)
        end if

        if (.not. this%isCurvilinear) then
           this%metrics(:,6) = 0.0_wp
        else
           F(:,1) = jacobianMatrixInverse(:,3) * this%coordinates(:,1)
           call this%firstDerivative(2)%apply(F, this%localSize)
           this%metrics(:,6) = F(:,1)
           F(:,1) = jacobianMatrixInverse(:,6) * this%coordinates(:,1)
           call this%firstDerivative(1)%apply(F, this%localSize)
           this%metrics(:,6) = this%metrics(:,6) - F(:,1)
        end if

        if (.not. this%isCurvilinear) then
           this%metrics(:,7) = 0.0_wp
        else
           F(:,1) = jacobianMatrixInverse(:,4) * this%coordinates(:,2)
           call this%firstDerivative(3)%apply(F, this%localSize)
           this%metrics(:,7) = F(:,1)
           F(:,1) = jacobianMatrixInverse(:,7) * this%coordinates(:,2)
           call this%firstDerivative(2)%apply(F, this%localSize)
           this%metrics(:,7) = this%metrics(:,7) - F(:,1)
        end if

        if (.not. this%isCurvilinear) then
           this%metrics(:,8) = 0.0_wp
        else
           F(:,1) = jacobianMatrixInverse(:,7) * this%coordinates(:,2)
           call this%firstDerivative(1)%apply(F, this%localSize)
           this%metrics(:,8) = F(:,1)
           F(:,1) = jacobianMatrixInverse(:,1) * this%coordinates(:,2)
           call this%firstDerivative(3)%apply(F, this%localSize)
           this%metrics(:,8) = this%metrics(:,8) - F(:,1)
        end if

        F(:,1) = jacobianMatrixInverse(:,1) * this%coordinates(:,2)
        call this%firstDerivative(2)%apply(F, this%localSize)
        this%metrics(:,9) = F(:,1)
        if (this%isCurvilinear) then
           F(:,1) = jacobianMatrixInverse(:,4) * this%coordinates(:,2)
           call this%firstDerivative(1)%apply(F, this%localSize)
           this%metrics(:,9) = this%metrics(:,9) - F(:,1)
        end if

        SAFE_DEALLOCATE(F)

        do i = 1, nDimensions ** 2
           where (this%iblank == 0) this%metrics(:,i) = 0.0_wp
        end do

     end if

  end select

  where (this%iblank == 0) this%jacobian(:,1) = 1.0_wp

  SAFE_DEALLOCATE(jacobianMatrixInverse)

  hasNegativeJacobian_ = .not. this%isVariableWithinRange(this%jacobian(:,1),                &
       jacobianOutsideRange, i, j, k, minValue = 0.0_wp)
  call MPI_Allreduce(MPI_IN_PLACE, hasNegativeJacobian_, 1,                                  &
       MPI_LOGICAL, MPI_LOR, this%comm, ierror)
  if (hasNegativeJacobian_) then
     write(message, '(A,4(I0.0,A),3(ES11.4,A))') "Jacobian on grid ", this%index, " at (",   &
          i, ", ", j, ", ", k, "): ", jacobianOutsideRange, " is non-positive!"
  end if

  if (present(hasNegativeJacobian)) then
     hasNegativeJacobian = hasNegativeJacobian_
     if (present(errorMessage)) then
        errorMessage = message
        return
     else if (hasNegativeJacobian) then
        call gracefulExit(this%comm, message)
     end if
  else if (hasNegativeJacobian_) then
     call gracefulExit(this%comm, message)
  end if

  ! Update the norm matrix, using the norm from the first derivative operator.
  assert(allocated(this%norm))
  assert(size(this%norm) == this%nGridPoints)
  this%norm = 1.0_wp
  do i = 1, nDimensions
     call this%firstDerivative(i)%applyNorm(this%norm, this%localSize)
  end do
  this%norm = this%norm * this%jacobian

  this%jacobian = 1.0_wp / this%jacobian

end subroutine updateGrid

function computeScalarInnerProduct(this, f, g, weight) result(innerProduct)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid

  implicit none

  ! <<< Arguments >>>
  class(t_Grid) :: this
  SCALAR_TYPE, intent(in) :: f(:), g(:)
  SCALAR_TYPE, intent(in), optional :: weight(:)

  ! <<< Result >>>
  SCALAR_TYPE :: innerProduct

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: ierror

  assert(this%nGridPoints > 0)
  assert(size(f) == this%nGridPoints)
  assert(size(g) == this%nGridPoints)
#ifdef DEBUG
  if (present(weight)) assert(size(weight) == this%nGridPoints)
#endif
  assert(allocated(this%norm))
  assert(size(this%norm) == this%nGridPoints)

  innerProduct = 0.0_wp

  if (present(weight)) then
     innerProduct = innerProduct + sum(f * this%norm(:,1) * g * weight)
  else
     innerProduct = innerProduct + sum(f * this%norm(:,1) * g)
  end if

#ifdef SCALAR_TYPE_IS_binary128_IEEE754
  assert(allocated(this%mpiReduceBuffer))
  call MPI_Allgather(innerProduct, 1, SCALAR_TYPE_MPI, this%mpiReduceBuffer,                 &
       1, SCALAR_TYPE_MPI, this%comm, ierror)
  innerProduct = sum(this%mpiReduceBuffer)
#else
  call MPI_Allreduce(MPI_IN_PLACE, innerProduct, 1, SCALAR_TYPE_MPI,                         &
       MPI_SUM, this%comm, ierror)
#endif

end function computeScalarInnerProduct

function computeVectorInnerProduct(this, f, g, weight) result(innerProduct)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid

  implicit none

  ! <<< Arguments >>>
  class(t_Grid) :: this
  SCALAR_TYPE, intent(in) :: f(:,:), g(:,:)
  SCALAR_TYPE, intent(in), optional :: weight(:)

  ! <<< Result >>>
  SCALAR_TYPE :: innerProduct

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, ierror

  assert(this%nGridPoints > 0)
  assert(size(f, 1) == this%nGridPoints)
  assert(size(g, 1) == this%nGridPoints)
  assert(size(g, 2) == size(f, 2))
#ifdef DEBUG
  if (present(weight)) assert(size(weight) == this%nGridPoints)
#endif
  assert(allocated(this%norm))
  assert(size(this%norm) == this%nGridPoints)

  innerProduct = 0.0_wp

  do i = 1, size(f, 2)
     if (present(weight)) then
        innerProduct = innerProduct + sum(f(:,i) * this%norm(:,1) * g(:,i) * weight)
     else
        innerProduct = innerProduct + sum(f(:,i) * this%norm(:,1) * g(:,i))
     end if
  end do

#ifdef SCALAR_TYPE_IS_binary128_IEEE754
  assert(allocated(this%mpiReduceBuffer))
  call MPI_Allgather(innerProduct, 1, SCALAR_TYPE_MPI, this%mpiReduceBuffer,                 &
       1, SCALAR_TYPE_MPI, this%comm, ierror)
  innerProduct = sum(this%mpiReduceBuffer)
#else
  call MPI_Allreduce(MPI_IN_PLACE, innerProduct, 1, SCALAR_TYPE_MPI,                         &
       MPI_SUM, this%comm, ierror)
#endif

end function computeVectorInnerProduct

subroutine computeGradientOfScalar(this, f, gradF)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  ! <<< Arguments >>>
  class(t_Grid) :: this
  SCALAR_TYPE, intent(in) :: f(:)
  SCALAR_TYPE, intent(out) :: gradF(:,:)

  ! <<< Local variables >>>
  integer :: nDimensions
  SCALAR_TYPE, allocatable :: temp(:,:)

  call startTiming("computeGradient")

  nDimensions = this%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  assert(this%nGridPoints > 0)
  assert(all(this%localSize > 0) .and. product(this%localSize) == this%nGridPoints)

  assert(size(f) == this%nGridPoints)
  assert(size(gradF, 1) == size(f))
  assert(size(gradF, 2) == nDimensions)

  assert(allocated(this%firstDerivative))
  assert(size(this%firstDerivative) == nDimensions)
  assert(allocated(this%metrics))
  assert(size(this%metrics, 1) == this%nGridPoints)
  assert(size(this%metrics, 2) == nDimensions ** 2)
  assert(allocated(this%jacobian))
  assert(size(this%jacobian) == this%nGridPoints)

  if (this%isCurvilinear .and. nDimensions > 1)                                              &
       allocate(temp(this%nGridPoints, nDimensions - 1))

  select case (nDimensions)

  case (1)
     gradF(:,1) = f
     call this%firstDerivative(1)%apply(gradF(:,1:1), this%localSize)
     gradF(:,1) = this%jacobian(:,1) * this%metrics(:,1) * gradF(:,1)

  case (2)

     gradF(:,1) = f
     call this%firstDerivative(1)%apply(gradF(:,1:1), this%localSize)

     gradF(:,2) = f
     call this%firstDerivative(2)%apply(gradF(:,2:2), this%localSize)

     if (this%isCurvilinear) then
        temp(:,1) = gradF(:,1)
        gradF(:,1) = this%jacobian(:,1) * (this%metrics(:,1) *                               &
             gradF(:,1) + this%metrics(:,3) * gradF(:,2))
        gradF(:,2) = this%jacobian(:,1) * (this%metrics(:,2) *                               &
             temp(:,1) + this%metrics(:,4) * gradF(:,2))
     else
        gradF(:,1) = this%jacobian(:,1) * this%metrics(:,1) * gradF(:,1)
        gradF(:,2) = this%jacobian(:,1) * this%metrics(:,4) * gradF(:,2)
     end if

  case (3)

     gradF(:,1) = f
     call this%firstDerivative(1)%apply(gradF(:,1:1), this%localSize)

     gradF(:,2) = f
     call this%firstDerivative(2)%apply(gradF(:,2:2), this%localSize)

     gradF(:,3) = f
     call this%firstDerivative(3)%apply(gradF(:,3:3), this%localSize)

     if (this%isCurvilinear) then
        temp(:,1:2) = gradF(:,1:2)
        gradF(:,1) = this%jacobian(:,1) * (this%metrics(:,1) * gradF(:,1) +                  &
             this%metrics(:,4) * gradF(:,2) + this%metrics(:,7) * gradF(:,3))
        gradF(:,2) = this%jacobian(:,1) * (this%metrics(:,2) * temp(:,1) +                   &
             this%metrics(:,5) * gradF(:,2) + this%metrics(:,8) * gradF(:,3))
        gradF(:,3) = this%jacobian(:,1) * (this%metrics(:,3) * temp(:,1) +                   &
             this%metrics(:,6) * temp(:,2) + this%metrics(:,9) * gradF(:,3))
     else
        gradF(:,1) = this%jacobian(:,1) * this%metrics(:,1) * gradF(:,1)
        gradF(:,2) = this%jacobian(:,1) * this%metrics(:,5) * gradF(:,2)
        gradF(:,3) = this%jacobian(:,1) * this%metrics(:,9) * gradF(:,3)
     end if

  end select

  SAFE_DEALLOCATE(temp)

  call endTiming("computeGradient")

end subroutine computeGradientOfScalar

subroutine computeGradientOfVector(this, f, gradF)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid

  ! <<< Internal modules >>>
  use MPITimingsHelper, only : startTiming, endTiming

  ! <<< Arguments >>>
  class(t_Grid) :: this
  SCALAR_TYPE, intent(in) :: f(:,:)
  SCALAR_TYPE, intent(out) :: gradF(:,:)

  ! <<< Local variables >>>
  integer :: nDimensions
  SCALAR_TYPE, allocatable :: temp(:,:)

  call startTiming("computeGradient")

  nDimensions = this%nDimensions
  assert_key(nDimensions, (1, 2, 3))

  assert(this%nGridPoints > 0)
  assert(all(this%localSize > 0) .and. product(this%localSize) == this%nGridPoints)

  assert(size(f, 1) == this%nGridPoints)
  assert(size(f, 2) == nDimensions)
  assert(size(gradF, 1) == size(f, 1))
  assert(size(gradF, 2) == nDimensions ** 2)

  assert(allocated(this%firstDerivative))
  assert(size(this%firstDerivative) == nDimensions)
  assert(allocated(this%metrics))
  assert(size(this%metrics, 1) == this%nGridPoints)
  assert(size(this%metrics, 2) == nDimensions ** 2)
  assert(allocated(this%jacobian))
  assert(size(this%jacobian) == this%nGridPoints)

  select case (nDimensions)

  case (1)

     gradF = f
     call this%firstDerivative(1)%apply(gradF(:,1:1), this%localSize)
     gradF(:,1) = this%jacobian(:,1) * this%metrics(:,1) * gradF(:,1)

  case (2)

     gradF(:,1:2) = f(:,1:2)
     call this%firstDerivative(1)%apply(gradF(:,1:2), this%localSize)

     gradF(:,3:4) = f(:,1:2)
     call this%firstDerivative(2)%apply(gradF(:,3:4), this%localSize)

     if (this%isCurvilinear) then

        allocate(temp(this%nGridPoints, 2), source = gradF(:,1:2))

        gradF(:,1) = this%jacobian(:,1) * (this%metrics(:,1) * gradF(:,1) +                  &
             this%metrics(:,3) * gradF(:,3))
        gradF(:,2) = this%jacobian(:,1) * (this%metrics(:,2) * temp(:,1) +                   &
             this%metrics(:,4) * gradF(:,3))
        gradF(:,3) = this%jacobian(:,1) * (this%metrics(:,1) * temp(:,2) +                   &
             this%metrics(:,3) * gradF(:,4))
        gradF(:,4) = this%jacobian(:,1) * (this%metrics(:,2) * temp(:,2) +                   &
             this%metrics(:,4) * gradF(:,4))

     else

        allocate(temp(this%nGridPoints, 1), source = gradF(:,2:2))

        gradF(:,1) = this%jacobian(:,1) * this%metrics(:,1) * gradF(:,1)
        gradF(:,2) = this%jacobian(:,1) * this%metrics(:,4) * gradF(:,3)
        gradF(:,3) = this%jacobian(:,1) * this%metrics(:,1) * temp(:,1)
        gradF(:,4) = this%jacobian(:,1) * this%metrics(:,4) * gradF(:,4)

     end if

  case (3)

     gradF(:,1:3) = f(:,1:3)
     call this%firstDerivative(1)%apply(gradF(:,1:3), this%localSize)

     gradF(:,4:6) = f(:,1:3)
     call this%firstDerivative(2)%apply(gradF(:,4:6), this%localSize)

     gradF(:,7:9) = f(:,1:3)
     call this%firstDerivative(3)%apply(gradF(:,7:9), this%localSize)

     if (this%isCurvilinear) then

        allocate(temp(this%nGridPoints, 3), source = gradF(:,1:3))

        gradF(:,1) = this%jacobian(:,1) * (this%metrics(:,1) * gradF(:,1) +                  &
             this%metrics(:,4) * gradF(:,4) + this%metrics(:,7) * gradF(:,7))
        gradF(:,2) = this%jacobian(:,1) * (this%metrics(:,2) * temp(:,1) +                   &
             this%metrics(:,5) * gradF(:,4) + this%metrics(:,8) * gradF(:,7))
        gradF(:,3) = this%jacobian(:,1) * (this%metrics(:,3) * temp(:,1) +                   &
             this%metrics(:,6) * gradF(:,4) + this%metrics(:,9) * gradF(:,7))
        gradF(:,4) = this%jacobian(:,1) * (this%metrics(:,1) * temp(:,2) +                   &
             this%metrics(:,4) * gradF(:,5) + this%metrics(:,7) * gradF(:,8))
        gradF(:,7) = this%jacobian(:,1) * (this%metrics(:,1) * temp(:,3) +                   &
             this%metrics(:,4) * gradF(:,6) + this%metrics(:,7) * gradF(:,9))

        temp(:,1) = gradF(:,8)
        gradF(:,8) = this%jacobian(:,1) * (this%metrics(:,2) * temp(:,3) +                   &
             this%metrics(:,5) * gradF(:,6) + this%metrics(:,8) * gradF(:,9))
        gradF(:,9) = this%jacobian(:,1) * (this%metrics(:,3) * temp(:,3) +                   &
             this%metrics(:,6) * gradF(:,6) + this%metrics(:,9) * gradF(:,9))

        temp(:,3) = gradF(:,5)
        gradF(:,5) = this%jacobian(:,1) * (this%metrics(:,2) * temp(:,2) +                   &
             this%metrics(:,5) * gradF(:,5) + this%metrics(:,8) * temp(:,1))
        gradF(:,6) = this%jacobian(:,1) * (this%metrics(:,3) * temp(:,2) +                   &
             this%metrics(:,6) * temp(:,3) + this%metrics(:,9) * temp(:,1))

     else

        allocate(temp(this%nGridPoints, 1))

        gradF(:,1) = this%jacobian(:,1) * this%metrics(:,1) * gradF(:,1)
        gradF(:,5) = this%jacobian(:,1) * this%metrics(:,5) * gradF(:,5)
        gradF(:,9) = this%jacobian(:,1) * this%metrics(:,9) * gradF(:,9)

        temp(:,1) = gradF(:,2)
        gradF(:,2) = this%jacobian(:,1) * this%metrics(:,5) * gradF(:,4)
        gradF(:,4) = this%jacobian(:,1) * this%metrics(:,1) * temp(:,1)

        temp(:,1) = gradF(:,3)
        gradF(:,3) = this%jacobian(:,1) * this%metrics(:,9) * gradF(:,7)
        gradF(:,7) = this%jacobian(:,1) * this%metrics(:,1) * temp(:,1)

        temp(:,1) = gradF(:,6)
        gradF(:,6) = this%jacobian(:,1) * this%metrics(:,9) * gradF(:,8)
        gradF(:,8) = this%jacobian(:,1) * this%metrics(:,5) * temp(:,1)

     end if

  end select

  SAFE_DEALLOCATE(temp)

  call endTiming("computeGradient")

end subroutine computeGradientOfVector

subroutine findMinimum(this, f, fMin, iMin, jMin, kMin)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid

  implicit none

  ! <<< Arguments >>>
  class(t_Grid) :: this
  SCALAR_TYPE, intent(in) :: f(:)
  SCALAR_TYPE, intent(out) :: fMin
  integer, intent(out), optional :: iMin, jMin, kMin

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, minIndex(3), nProcs, ierror
  real(wp), allocatable :: minValues(:)
  integer, allocatable :: minIndices(:,:)
  SCALAR_TYPE :: a
  real(wp) :: minValue

  assert(this%nGridPoints > 0)
  assert(all(this%localSize > 0) .and. product(this%localSize) == this%nGridPoints)
  assert(all(this%offset >= 0))

  assert(size(f) == this%nGridPoints)

  call MPI_Comm_size(this%comm, nProcs, ierror)
  assert(nProcs > 0)

  allocate(minValues(nProcs), minIndices(3, nProcs))

  minValue = huge(0.0_wp)
  do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
     do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
        do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
           a = f(i - this%offset(1) + this%localSize(1) * (j - 1 - this%offset(2) +          &
                this%localSize(2) * (k - 1 - this%offset(3))))
#ifdef SCALAR_IS_COMPLEX
           if (real(a, wp) < minValue) then
              minIndex = (/ i, j, k /)
              minValue = real(a, wp)
           end if
#else
           if (a < minValue) then
              minIndex = (/ i, j, k /)
              minValue = a
           end if
#endif
        end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)

  if (present(iMin) .or. present(jMin) .or. present(kMin)) then
     call MPI_Allgather(minIndex, 3, MPI_INTEGER, minIndices, 3, MPI_INTEGER,                &
          this%comm, ierror)
  end if
  call MPI_Allgather(minValue, 1, REAL_TYPE_MPI, minValues,                                  &
       1, REAL_TYPE_MPI, this%comm, ierror)

#ifdef SCALAR_IS_COMPLEX
  i = minloc(real(minValues, wp), 1)
#else
  i = minloc(minValues, 1)
#endif
  fMin = minValues(i)
  if (present(iMin)) iMin = minIndices(1,i)
  if (present(jMin)) jMin = minIndices(2,i)
  if (present(kMin)) kMin = minIndices(3,i)

  SAFE_DEALLOCATE(minValues)
  SAFE_DEALLOCATE(minIndices)

end subroutine findMinimum

subroutine findMaximum(this, f, fMax, iMax, jMax, kMax)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid

  implicit none

  ! <<< Arguments >>>
  class(t_Grid) :: this
  SCALAR_TYPE, intent(in) :: f(:)
  SCALAR_TYPE, intent(out) :: fMax
  integer, intent(out), optional :: iMax, jMax, kMax

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, j, k, maxIndex(3), nProcs, ierror
  real(wp), allocatable :: maxValues(:)
  integer, allocatable :: maxIndices(:,:)
  SCALAR_TYPE :: a
  real(wp) :: maxValue

  assert(this%nGridPoints > 0)
  assert(all(this%localSize > 0) .and. product(this%localSize) == this%nGridPoints)
  assert(all(this%offset >= 0))

  assert(size(f) == this%nGridPoints)

  call MPI_Comm_size(this%comm, nProcs, ierror)
  assert(nProcs > 0)

  allocate(maxValues(nProcs), maxIndices(3, nProcs))

  maxValue = - huge(0.0_wp)
  do k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
     do j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
        do i = this%offset(1) + 1, this%offset(1) + this%localSize(1)
           a = f(i - this%offset(1) + this%localSize(1) * (j - 1 - this%offset(2) +          &
                this%localSize(2) * (k - 1 - this%offset(3))))
#ifdef SCALAR_IS_COMPLEX
           if (real(a, wp) > maxValue) then
              maxIndex = (/ i, j, k /)
              maxValue = real(a, wp)
           end if
#else
           if (a > maxValue) then
              maxIndex = (/ i, j, k /)
              maxValue = a
           end if
#endif
        end do !... k = this%offset(3) + 1, this%offset(3) + this%localSize(3)
     end do !... j = this%offset(2) + 1, this%offset(2) + this%localSize(2)
  end do !... i = this%offset(1) + 1, this%offset(1) + this%localSize(1)

  if (present(iMax) .or. present(jMax) .or. present(kMax)) then
     call MPI_Allgather(maxIndex, 3, MPI_INTEGER, maxIndices, 3, MPI_INTEGER,                &
          this%comm, ierror)
  end if
  call MPI_Allgather(maxValue, 1, REAL_TYPE_MPI, maxValues,                                  &
       1, REAL_TYPE_MPI, this%comm, ierror)

#ifdef SCALAR_IS_COMPLEX
  i = maxloc(real(maxValues, wp), 1)
#else
  i = maxloc(maxValues, 1)
#endif
  fMax = maxValues(i)
  if (present(iMax)) iMax = maxIndices(1,i)
  if (present(jMax)) jMax = maxIndices(2,i)
  if (present(kMax)) kMax = maxIndices(3,i)

  SAFE_DEALLOCATE(maxValues)
  SAFE_DEALLOCATE(maxIndices)

end subroutine findMaximum

function isVariableWithinRange(this, f, fOutsideRange,                                       &
     iOutsideRange, jOutsideRange, kOutsideRange, minValue, maxValue)

  ! <<< External modules >>>
  use MPI

  ! <<< Derived types >>>
  use Grid_mod, only : t_Grid

  implicit none

  ! <<< Arguments >>>
  class(t_Grid) :: this
  SCALAR_TYPE, intent(in) :: f(:)
  SCALAR_TYPE, intent(out) :: fOutsideRange
  integer, intent(out) :: iOutsideRange, jOutsideRange, kOutsideRange
  real(SCALAR_KIND), intent(in), optional :: minValue, maxValue
  logical :: isVariableWithinRange

  ! <<< Local variables >>>
  SCALAR_TYPE :: fMin, fMax

  assert(this%nGridPoints > 0)
  assert(size(f) == this%nGridPoints)

  isVariableWithinRange = .true.
  if (.not. present(minValue) .and. .not. present(maxValue)) return

  if (present(minValue)) then
     call this%findMinimum(f, fMin, iOutsideRange, jOutsideRange, kOutsideRange)
     if (real(fMin, SCALAR_KIND) <= minValue) then
        fOutsideRange = fMin
        isVariableWithinRange = .false.
     end if
  end if

  if (present(maxValue)) then
     call this%findMaximum(f, fMax, iOutsideRange, jOutsideRange, kOutsideRange)
     if (real(fMax, SCALAR_KIND) >= maxValue) then
        fOutsideRange = fMax
        isVariableWithinRange = .false.
     end if
  end if

end function isVariableWithinRange
