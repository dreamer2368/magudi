#include "config.h"

module ResidualManager_mod

  integer, parameter :: wp = SCALAR_KIND

  type, public :: t_ResidualManager

     real(SCALAR_KIND), allocatable :: residuals(:), tolerances(:)
     integer :: reportInterval = 0
     logical :: hasSimulationConverged = .false.

   contains

     procedure, pass :: setup => setupResidualManager
     procedure, pass :: cleanup => cleanupResidualManager
     procedure, pass :: compute => computeResiduals
     procedure, pass :: writeToFile => writeResidualsToFile

  end type t_ResidualManager

  interface

     subroutine setupResidualManager(this, prefix, region)

       use Region_mod, only : t_Region

       import :: t_ResidualManager

       class(t_ResidualManager) :: this
       character(len = *), intent(in) :: prefix
       class(t_Region) :: region

     end subroutine setupResidualManager

  end interface

  interface

     subroutine cleanupResidualManager(this)

       import :: t_ResidualManager

       class(t_ResidualManager) :: this

     end subroutine cleanupResidualManager

  end interface

  interface

     subroutine computeResiduals(this, region)

       use Region_mod, only : t_Region

       import :: t_ResidualManager

       class(t_ResidualManager) :: this
       class(t_Region) :: region

     end subroutine computeResiduals

  end interface

  interface

     subroutine writeResidualsToFile(this, comm, filename, timestep, time, append)

       import :: t_ResidualManager

       class(t_ResidualManager) :: this
       integer, intent(in) :: comm
       character(len = *), intent(in) :: filename
       integer, intent(in) :: timestep
       real(SCALAR_KIND), intent(in) :: time
       logical, intent(in), optional :: append

     end subroutine writeResidualsToFile

  end interface

end module ResidualManager_mod
