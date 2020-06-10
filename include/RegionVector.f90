#include "config.h"

module RegionVector_mod

  implicit none

  public :: operator(*), operator(+), operator(-)

  type, private :: t_InternalState
    SCALAR_TYPE, allocatable :: conservedVariables(:,:)
  end type

  type, public :: t_RegionVector

    type(t_InternalState), allocatable :: states(:)
    SCALAR_TYPE, allocatable :: params(:)

  contains

    procedure, pass :: cleanup => cleanupRegionVector
    procedure, pass, private :: setFromRegion, setFromRegionVector
    generic, public :: set => setFromRegion, setFromRegionVector

  end type t_RegionVector

  interface

     subroutine cleanupRegionVector(this)

       import :: t_RegionVector

       class(t_RegionVector) :: this

     end subroutine cleanupRegionVector

  end interface

  interface

    subroutine setFromRegionVector(this,x)

      import :: t_RegionVector

      class(t_RegionVector), intent(inout) :: this
      class(t_RegionVector), intent(in) :: x

    end subroutine setFromRegionVector

  end interface

  interface

    subroutine setFromRegion(this,x)

      use Region_mod, only : t_Region
      import :: t_RegionVector

      class(t_RegionVector), intent(inout) :: this
      class(t_Region), intent(in) :: x

    end subroutine setFromRegion

  end interface

  interface assignment (=)

     subroutine copyRegionVector(b, a)

       import :: t_RegionVector

       class(t_RegionVector), intent(inout) :: b
       class(t_RegionVector), intent(in) :: a

     end subroutine copyRegionVector

  end interface

  interface operator (+)

     function addRegionVector(a,b) result(r)

       import :: t_RegionVector

       class(t_RegionVector), intent(in) :: a, b
       type(t_RegionVector) :: r

     end function addRegionVector

  end interface

  interface operator (-)

     function subtractRegionVector(a,b) result(r)

       import :: t_RegionVector

       class(t_RegionVector), intent(in) :: a, b
       type(t_RegionVector) :: r

     end function subtractRegionVector

  end interface

  interface operator(*)

    function multiplyRegionVector(a,b) result(r)

      import :: t_RegionVector

      class(t_RegionVector), intent(in) :: a
      SCALAR_TYPE, intent(in) :: b
      type(t_RegionVector) :: r

    end function multiplyRegionVector

  end interface

end module RegionVector_mod
