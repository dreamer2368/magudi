#include "config.h"

module MPITimingsHelper

  implicit none
  public

  interface

     subroutine startTiming(name)

       character(len = *), intent(in) :: name

     end subroutine startTiming

  end interface

  interface

     subroutine endTiming(name)

       character(len = *), intent(in) :: name

     end subroutine endTiming

  end interface

  interface

     subroutine reportTimings(comm, outputUnit)
       
       integer, intent(in), optional :: comm, outputUnit
       
     end subroutine reportTimings
     
  end interface

  interface

     subroutine cleanupTimers()

     end subroutine cleanupTimers

  end interface

end module MPITimingsHelper
