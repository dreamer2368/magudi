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

     subroutine exchangeInterfaceData(region)

       use Region_mod, only : t_Region

       class(t_Region) :: region

     end subroutine exchangeInterfaceData

  end interface

end module InterfaceHelper
