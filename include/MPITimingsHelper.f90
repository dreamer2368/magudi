#include "config.h"

module MPITimingsHelper

  implicit none
  public

  interface

     subroutine startTiming(name)

       !> Starts the timer labeled `name`. Creates a new timer if no timer with label `name`
       !> exists. The timers are internally stored as a sorted binary tree to minimize the
       !> overhead of starting and stopping them.

       character(len = *), intent(in) :: name

     end subroutine startTiming

  end interface

  interface

     subroutine endTiming(name)

       !> Stops the timer labeled `name`. Does nothing if no timer lebeled `name` exists. The
       !> number of calls is incremented by 1, and the accumulated time is incremented by the
       !> time elapsed since the most recent call to `startTiming` with `name` as the
       !> argument.

       character(len = *), intent(in) :: name

     end subroutine endTiming

  end interface

  interface

     subroutine reportTimings(comm, outputUnit)

       !> Reports timings. This subroutine must be called collectively by all processes in the
       !> MPI communicator `comm`. The master process in `comm` gathers the timers from all
       !> other processes. The timers are then sorted by their labels. Timers from different
       !> processes with the same label are clubbed together by adding their number of calls
       !> and accumulated time. The total program time is counted as the time elapsed since
       !> the first call to `startTiming`.

       integer, intent(in), optional :: comm, outputUnit

     end subroutine reportTimings

  end interface

  interface

     subroutine cleanupTimers()

       !> Frees memory that was used to internally store the timers.

     end subroutine cleanupTimers

  end interface

end module MPITimingsHelper
