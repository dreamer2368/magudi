#include "config.h"

module ActuatorPatch_mod

  use MPI, only : MPI_OFFSET_KIND

  use Patch_mod, only : t_Patch

  implicit none

  integer, parameter, private :: wp = SCALAR_KIND

  type, extends(t_Patch), public :: t_ActuatorPatch

    !bufferOffsetIndex starts from 0 at the beginning of the file. Used only for checkpointing.
    !referenceTimestep is the number of timesteps between intermediate start timestep and initial condition.
     integer :: iGradientBuffer = 0, iControlForcingBuffer = 0,                               &
                bufferOffsetIndex = -1, forwardReferenceTimestep = -1, adjointReferenceTimestep = -1
     integer(kind = MPI_OFFSET_KIND) :: gradientFileOffset = int(0, MPI_OFFSET_KIND),         &
                                        controlForcingFileOffset = int(0, MPI_OFFSET_KIND),   &
                                        controlForcingFileSize = int(0, MPI_OFFSET_KIND)
     character(len = STRING_LENGTH) :: gradientFilename, controlForcingFilename
     SCALAR_TYPE, allocatable :: controlForcing(:,:),                                         &
                                 gradientBuffer(:,:,:), controlForcingBuffer(:,:,:)

   contains

     procedure, pass :: setup => setupActuatorPatch
     procedure, pass :: cleanup => cleanupActuatorPatch
     procedure, pass :: verifyUsage => verifyActuatorPatchUsage
     procedure, pass :: updateRhs => updateActuatorPatch
     procedure, pass :: loadForcing => loadActuatorForcing
     procedure, pass :: pinpointForcing => pinpointActuatorForcing
     procedure, pass :: saveGradient => saveActuatorGradient

  end type t_ActuatorPatch

  interface

     subroutine setupActuatorPatch(this, index, comm, patchDescriptor,                       &
          grid, simulationFlags, solverOptions)

       use Grid_mod, only : t_Grid
       use SolverOptions_mod, only : t_SolverOptions
       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ActuatorPatch

       class(t_ActuatorPatch) :: this
       integer, intent(in) :: index, comm
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       class(t_Grid), intent(in) :: grid
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions

     end subroutine setupActuatorPatch

  end interface

  interface

     subroutine cleanupActuatorPatch(this)

       import :: t_ActuatorPatch

       class(t_ActuatorPatch) :: this

     end subroutine cleanupActuatorPatch

  end interface

  interface

     function verifyActuatorPatchUsage(this, patchDescriptor, gridSize, normalDirection,     &
          extent, simulationFlags, success, message) result(isPatchUsed)

       use PatchDescriptor_mod, only : t_PatchDescriptor
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ActuatorPatch

       class(t_ActuatorPatch) :: this
       type(t_PatchDescriptor), intent(in) :: patchDescriptor
       integer, intent(in) :: gridSize(:), normalDirection, extent(6)
       type(t_SimulationFlags), intent(in) :: simulationFlags
       logical, intent(out) :: success
       character(len = STRING_LENGTH), intent(out) :: message

       logical :: isPatchUsed

     end function verifyActuatorPatchUsage

  end interface

  interface

     subroutine updateActuatorPatch(this, mode, simulationFlags, solverOptions, grid, state)

       use Grid_mod, only : t_Grid
       use State_mod, only : t_State
       use SolverOptions_mod, only : t_SolverOptions
       use SimulationFlags_mod, only : t_SimulationFlags

       import :: t_ActuatorPatch

       class(t_ActuatorPatch) :: this
       integer, intent(in) :: mode
       type(t_SimulationFlags), intent(in) :: simulationFlags
       type(t_SolverOptions), intent(in) :: solverOptions
       class(t_Grid), intent(in) :: grid
       class(t_State) :: state

     end subroutine updateActuatorPatch

  end interface

  interface

     subroutine loadActuatorForcing(this)

       import :: t_ActuatorPatch

       class(t_ActuatorPatch) :: this

     end subroutine loadActuatorForcing

  end interface

  interface

     subroutine pinpointActuatorForcing(this,bufferOffsetIndex)

       import :: t_ActuatorPatch

       class(t_ActuatorPatch) :: this
       integer, intent(in) :: bufferOffsetIndex

     end subroutine pinpointActuatorForcing

  end interface

  interface

     subroutine saveActuatorGradient(this)

       import :: t_ActuatorPatch

       class(t_ActuatorPatch) :: this

     end subroutine saveActuatorGradient

  end interface

end module ActuatorPatch_mod
