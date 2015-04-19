#include "config.h"

module InterfaceHelper

  implicit none
  public

  interface

     subroutine readPatchInterfaceInformation(region)

       use Region_mod, only : t_Region

       class(t_Region) :: region

     end subroutine readPatchInterfaceInformation

  end interface

  interface

     subroutine exchangeInterfaceData(region, mode)

       use Region_mod, only : t_Region

       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine exchangeInterfaceData

  end interface

  interface

     subroutine checkInterfaceContinuity(region, tolerance, success)

       use Region_mod, only : t_Region

       class(t_Region) :: region
       real(SCALAR_KIND), intent(in) :: tolerance
       logical, intent(out) :: success

     end subroutine checkInterfaceContinuity

  end interface

end module InterfaceHelper
