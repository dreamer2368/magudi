#include "config.h"

module RegionVector_mod

  implicit none

  type, private :: t_InternalState
    SCALAR_TYPE, allocatable :: conservedVariables(:,:)
  end type

  type, public :: t_RegionVector

    type(t_InternalState), allocatable :: states(:)
    SCALAR_TYPE, allocatable :: params(:)

  contains

    procedure, pass :: cleanup => cleanupRegionVector

  end type t_RegionVector

  interface

     subroutine cleanupRegionVector(this)

       import :: t_RegionVector

       class(t_RegionVector) :: this

     end subroutine cleanupRegionVector

  end interface

  interface

     function addRegionVector(a,b) result(r)

       import :: t_RegionVector

       class(t_RegionVector), intent(in) :: a, b
       type(t_RegionVector) :: r

     end function addRegionVector

  end interface

  interface

     function subtractRegionVector(a,b) result(r)

       import :: t_RegionVector

       class(t_RegionVector), intent(in) :: a, b
       type(t_RegionVector) :: r

     end function subtractRegionVector

  end interface

  interface

    function multiplyRegionVector(a,b) result(r)

      import :: t_RegionVector

      SCALAR_TYPE, intent(in) :: a
      class(t_RegionVector), intent(in) :: b
      type(t_RegionVector) :: r

    end function multiplyRegionVector

  end interface

  interface

    function innerProduct(a,b,region) result(r)

      use Region_mod, only : t_Region
      import :: t_RegionVector

      class(t_RegionVector), intent(in) :: a, b
      class(t_Region), intent(in) :: region
      SCALAR_TYPE :: r

    end function innerProduct

  end interface

end module RegionVector_mod
