#include "config.h"

program append_gradient

  use MPI
  use, intrinsic :: iso_fortran_env, only : output_unit

  use Region_mod, only : t_Region
  use Solver_mod, only : t_Solver
  use Patch_mod, only : t_Patch

  use Controller_mod, only : t_Controller
  use Functional_mod, only : t_Functional
  use ActuatorPatch_mod, only : t_ActuatorPatch
  use TimeIntegrator_mod, only : t_TimeIntegrator

  ! <<< Enumerations >>>
  use Grid_enum
  use State_enum
  use Region_enum, only : FORWARD

  use InputHelper, only : parseInputFile, getOption, getRequiredOption
  use ErrorHandler
  use PLOT3DHelper, only : plot3dDetectFormat, plot3dErrorMessage
  use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

  implicit none

  class(t_TimeIntegrator), pointer :: timeIntegrator => null()
  class(t_Controller), pointer :: controller => null()
  class(t_Functional), pointer :: functional => null()
  class(t_Patch), pointer :: patch => null()
  integer, parameter :: wp = SCALAR_KIND
  integer :: i,j,k,numIter, procRank, numProcs, ierror
  character(len = STRING_LENGTH) :: filename, outputPrefix, message, outputFilename
  integer(kind = MPI_OFFSET_KIND) :: outputGradientFileOffset
  integer :: firstTimesteps, secondTimesteps
  logical :: success
  integer, dimension(:,:), allocatable :: globalGridSizes
  type(t_Region) :: region
  type(t_Solver) :: solver

  ! Initialize MPI.
  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, procRank, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierror)

  call initializeErrorHandler()

  if (command_argument_count() > 1) then
     write(message, '(A)') "Usage: magudi [FILE]"
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)') "High-performance Fortran-based adjoint optimization tool."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     write(message, '(A)')                                                                   &
          "FILE is an optional restart file used if running in prediction-only mode."
     call writeAndFlush(MPI_COMM_WORLD, output_unit, message)
     call cleanupErrorHandler()
     call MPI_Finalize(ierror)
     stop -1
  end if

  call startTiming("total")

  ! Parse options from the input file.
  filename = PROJECT_NAME // ".inp"
  call parseInputFile(filename)

  outputPrefix = getOption("output_prefix", PROJECT_NAME)

  ! Verify that the grid file is in valid PLOT3D format and fetch the grid dimensions:
  ! `globalGridSizes(i,j)` is the number of grid points on grid `j` along dimension `i`.
  call getRequiredOption("grid_file", filename)
  call plot3dDetectFormat(MPI_COMM_WORLD, filename, success,                                 &
       globalGridSizes = globalGridSizes)
  if (.not. success) call gracefulExit(MPI_COMM_WORLD, plot3dErrorMessage)

  ! Setup the region.
  call region%setup(MPI_COMM_WORLD, globalGridSizes)

  ! Read the grid file.
  call getRequiredOption("grid_file", filename)
  call region%loadData(QOI_GRID, filename)

  ! Update the grids by computing the Jacobian, metrics, and norm.
  do i = 1, size(region%grids)
     call region%grids(i)%update()
  end do
  call MPI_Barrier(region%comm, ierror)

  ! Write out some useful information.
  call region%reportGridDiagnostics()

  ! Save the Jacobian and normalized metrics.
  write(filename, '(2A)') trim(outputPrefix), ".Jacobian.f"
  call region%saveData(QOI_JACOBIAN, filename)
  write(filename, '(2A)') trim(outputPrefix), ".metrics.f"
  call region%saveData(QOI_METRICS, filename)

  ! Initialize the solver.
  call solver%setup(region, outputPrefix = outputPrefix)

  ! Save the control and target mollifier if using code-generated values.
  if (region%simulationFlags%enableController) then
     filename = getOption("control_mollifier_file", "")
     if (len_trim(filename) == 0) call region%saveData(QOI_CONTROL_MOLLIFIER,                &
          trim(outputPrefix) // ".control_mollifier.f")
  end if
  if (region%simulationFlags%enableFunctional) then
     filename = getOption("target_mollifier_file", "")
     if (len_trim(filename) == 0) call region%saveData(QOI_TARGET_MOLLIFIER,                 &
          trim(outputPrefix) // ".target_mollifier.f")
  end if

  ! Main code logic.

  ! Connect to the previously allocated time integrator.
  call solver%timeIntegratorFactory%connect(timeIntegrator)
  assert(associated(timeIntegrator))

  ! Connect to the previously allocated controller.
  if (region%simulationFlags%enableController) then
     call solver%controllerFactory%connect(controller)
     assert(associated(controller))
  end if

  ! Connect to the previously allocated functional.
  if (region%simulationFlags%enableFunctional) then
     call solver%functionalFactory%connect(functional)
     assert(associated(functional))
     functional%runningTimeQuadrature = 0.0_wp
  end if

  ! Setup residual manager if this is a steady-state simulation.
  if (region%simulationFlags%steadyStateSimulation)                                          &
       call solver%residualManager%setup("", region)

  ! Call controller hooks before time marching starts.
  if (region%simulationFlags%enableController .and.                                          &
       controller%controllerSwitch) then
!     controller%onsetTime = startTime
     controller%duration = solver%nTimesteps * region%solverOptions%timeStepSize
     call controller%hookBeforeTimemarch(region, FORWARD)
  end if

  ! Reset probes.
  if (solver%probeInterval > 0) call region%resetProbes()

  !Determine size of substeps, number of iteration for loading/saving
  call getRequiredOption("merge_gradient/first_timesteps", firstTimesteps)
  call getRequiredOption("merge_gradient/second_timesteps", secondTimesteps)
  write(message,'(A,I8.8)') 'First iteration timesteps: ',firstTimesteps
  call writeAndFlush(region%comm, output_unit, message)
  write(message,'(A,I8.8)') 'Second iteration timesteps: ',secondTimesteps
  call writeAndFlush(region%comm, output_unit, message)

  do i = 1, size(region%patchFactories)
     call region%patchFactories(i)%connect(patch)
     if (.not. associated(patch)) cycle
     do j = 1, size(region%states)
        if (patch%gridIndex /= region%grids(j)%index) cycle
        select type (patch)
         class is (t_ActuatorPatch)
           numIter = firstTimesteps*timeIntegrator%nStages/size(patch%gradientBuffer,3)
           outputFilename = 'merged_'//trim(patch%gradientFilename)
           patch%iGradientBuffer = size(patch%gradientBuffer,3)
           patch%gradientFileOffset = 0
        end select
     end do
  end do

  write(message,'(A,I8.8)') 'First iteration: ',numIter
  call writeAndFlush(region%comm, output_unit, message)

  !Loading/saving iteration
  outputGradientFileOffset = int(0, MPI_OFFSET_KIND)
  do k = 1, numIter
    do i = 1, size(region%patchFactories)
       call region%patchFactories(i)%connect(patch)
       if (.not. associated(patch)) cycle
       do j = 1, size(region%states)
          if (patch%gridIndex /= region%grids(j)%index) cycle
          select type (patch)
           class is (t_ActuatorPatch)
             call loadInputGradient(patch)
             call saveOutputGradient(patch,outputFilename,                                   &
                                     outputGradientFileOffset)
          end select
       end do
    end do
  end do

  write(message,'(A)') 'First reading/saving is done!'
  call writeAndFlush(region%comm, output_unit, message)

  !Additional loading/saving iteration setup
  do i = 1, size(region%patchFactories)
     call region%patchFactories(i)%connect(patch)
     if (.not. associated(patch)) cycle
     do j = 1, size(region%states)
        if (patch%gridIndex /= region%grids(j)%index) cycle
        select type (patch)
         class is (t_ActuatorPatch)
           numIter = secondTimesteps*timeIntegrator%nStages/size(patch%gradientBuffer,3)
           patch%gradientFilename = 'after_'//trim(patch%gradientFilename)
           patch%iGradientBuffer = size(patch%gradientBuffer,3)
           patch%gradientFileOffset = 0
        end select
     end do
  end do

  write(message,'(A,I8.8)') 'Second iteration: ',numIter
  call writeAndFlush(region%comm, output_unit, message)

  !Additional Loading/saving iteration
  do k = 1, numIter
    do i = 1, size(region%patchFactories)
       call region%patchFactories(i)%connect(patch)
       if (.not. associated(patch)) cycle
       do j = 1, size(region%states)
          if (patch%gridIndex /= region%grids(j)%index) cycle
          select type (patch)
           class is (t_ActuatorPatch)
             call loadInputGradient(patch)
             call saveOutputGradient(patch,outputFilename,                                   &
                                     outputGradientFileOffset)
          end select
       end do
    end do
  end do

  ! Call controller hooks after time marching ends.
  if (region%simulationFlags%enableController .and.                                          &
       controller%controllerSwitch)                                                          &
       call controller%hookAfterTimemarch(region, FORWARD)

  call solver%residualManager%cleanup()

  call solver%cleanup()
  call region%cleanup()

  call endTiming("total")
  call reportTimings()
  call cleanupTimers()

  ! Finalize MPI.
  call cleanupErrorHandler()
  call MPI_Finalize(ierror)

contains

  subroutine loadInputGradient(this)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use ActuatorPatch_mod, only : t_ActuatorPatch

    ! <<< Arguments >>>
    class(t_ActuatorPatch) :: this

    ! <<< Local variables >>>
    integer :: arrayOfSizes(5), arrayOfSubsizes(5), arrayOfStarts(5),                          &
         mpiScalarSubarrayType, mpiFileHandle, dataSize, ierror
    integer(kind = MPI_OFFSET_KIND) :: nBytesToRead
    character(len = STRING_LENGTH) :: message

    if (this%comm == MPI_COMM_NULL) return

    arrayOfSizes(1:3) = this%globalSize
    arrayOfSizes(4) = size(this%gradientBuffer, 2)
    arrayOfSizes(5) = this%iGradientBuffer
    arrayOfSubsizes(1:3) = this%localSize
    arrayOfSubsizes(4) = size(this%gradientBuffer, 2)
    arrayOfSubsizes(5) = this%iGradientBuffer
    arrayOfStarts(1:3) = this%offset - this%extent(1::2) + 1
    arrayOfStarts(4:5) = 0
    call MPI_Type_create_subarray(5, arrayOfSizes, arrayOfSubsizes, arrayOfStarts,             &
         MPI_ORDER_FORTRAN, SCALAR_TYPE_MPI, mpiScalarSubarrayType, ierror)
    call MPI_Type_commit(mpiScalarSubarrayType, ierror)

    call MPI_File_open(this%comm, trim(this%gradientFilename) // char(0), MPI_MODE_RDONLY,     &
         MPI_INFO_NULL, mpiFileHandle, ierror)

    assert(this%gradientFileOffset >= 0)

    call MPI_File_set_view(mpiFileHandle, this%gradientFileOffset, SCALAR_TYPE_MPI,            &
         mpiScalarSubarrayType, "native", MPI_INFO_NULL, ierror)

    dataSize = this%nPatchPoints * size(this%gradientBuffer, 2) * this%iGradientBuffer
    call MPI_File_read_all(mpiFileHandle, this%gradientBuffer, dataSize,                       &
         SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)
    call MPI_File_close(mpiFileHandle, ierror)

    nBytesToRead = SIZEOF_SCALAR * product(int(this%globalSize, MPI_OFFSET_KIND)) *            &
         size(this%gradientBuffer, 2) * this%iGradientBuffer
    if (this%gradientFileOffset + nBytesToRead <= int(0, MPI_OFFSET_KIND)) then
 !      this%iGradientBuffer = int(this%gradientFileOffset / (SIZEOF_SCALAR * &
 !           product(int(this%globalSize, MPI_OFFSET_KIND)) * size(this%gradientBuffer, 2)))
 !      this%gradientFileOffset = 0
    else
       this%gradientFileOffset = this%gradientFileOffset + nBytesToRead
    end if
!   print *, 'reading is being done..',nBytesToRead,this%gradientFileOffset
    write(message,'(A,I16.8,I16.8)') 'reading is being done..',nBytesToRead,this%gradientFileOffset
    call writeAndFlush(this%comm, output_unit, message)

    call MPI_Type_free(mpiScalarSubarrayType, ierror)

  end subroutine loadInputGradient

  subroutine saveOutputGradient(this,outputFilename,outputGradientFileOffset)

    ! <<< External modules >>>
    use MPI

    ! <<< Derived types >>>
    use ActuatorPatch_mod, only : t_ActuatorPatch

    ! <<< Arguments >>>
    class(t_ActuatorPatch) :: this
    character(len = STRING_LENGTH), intent(in) :: outputFilename
    integer(kind = MPI_OFFSET_KIND), intent(inout) :: outputGradientFileOffset

    ! <<< Local variables >>>
    integer :: arrayOfSizes(5), arrayOfSubsizes(5), arrayOfStarts(5),                          &
         mpiScalarSubarrayType, mpiFileHandle, dataSize, ierror
    integer(kind = MPI_OFFSET_KIND) :: nBytesToWrite
    character(len = STRING_LENGTH) :: message

    if (this%comm == MPI_COMM_NULL .or. this%iGradientBuffer == 0) return

    arrayOfSizes(1:3) = this%globalSize
    arrayOfSizes(4) = size(this%gradientBuffer, 2)
    arrayOfSizes(5) = this%iGradientBuffer
    arrayOfSubsizes(1:3) = this%localSize
    arrayOfSubsizes(4) = size(this%gradientBuffer, 2)
    arrayOfSubsizes(5) = this%iGradientBuffer
    arrayOfStarts(1:3) = this%offset - this%extent(1::2) + 1
    arrayOfStarts(4:5) = 0
    call MPI_Type_create_subarray(5, arrayOfSizes, arrayOfSubsizes, arrayOfStarts,             &
         MPI_ORDER_FORTRAN, SCALAR_TYPE_MPI, mpiScalarSubarrayType, ierror)
    call MPI_Type_commit(mpiScalarSubarrayType, ierror)

    call MPI_File_open(this%comm, trim(outputFilename) // char(0), MPI_MODE_WRONLY,            &
         MPI_INFO_NULL, mpiFileHandle, ierror)

    assert(outputGradientFileOffset >= 0)

    call MPI_File_set_view(mpiFileHandle, outputGradientFileOffset, SCALAR_TYPE_MPI,           &
         mpiScalarSubarrayType, "native", MPI_INFO_NULL, ierror)

    dataSize = this%nPatchPoints * size(this%gradientBuffer, 2) * this%iGradientBuffer
    call MPI_File_write_all(mpiFileHandle, this%gradientBuffer, dataSize,                      &
         SCALAR_TYPE_MPI, MPI_STATUS_IGNORE, ierror)

    nBytesToWrite = SIZEOF_SCALAR * product(int(this%globalSize, MPI_OFFSET_KIND)) *            &
         size(this%gradientBuffer, 2) * this%iGradientBuffer
    outputGradientFileOffset = outputGradientFileOffset + nBytesToWrite
! print *, 'writing is being done..',nBytesToWrite,outputGradientFileOffset
    write(message,'(A,I16.8,I16.8)') 'writing is being done..',nBytesToWrite,outputGradientFileOffset
    call writeAndFlush(this%comm, output_unit, message)

    call MPI_File_close(mpiFileHandle, ierror)

    call MPI_Type_free(mpiScalarSubarrayType, ierror)

  end subroutine saveOutputGradient

end program append_gradient
