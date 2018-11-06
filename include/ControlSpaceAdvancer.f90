#include "config.h"

module ControlSpaceAdvancer

  implicit none
  public

  interface
    subroutine ZAXPY(comm,ZFilename,A,XFilename,YFilename)
      use MPI
      use, intrinsic :: iso_fortran_env, only : output_unit

      use InputHelper, only : parseInputFile, getOption, getRequiredOption
      use ErrorHandler
      use MPITimingsHelper, only : startTiming, endTiming, reportTimings, cleanupTimers

      ! <<< Arguments >>>
      integer, intent(in) :: comm
      SCALAR_TYPE, intent(in) :: A
      character(len = *), intent(in) :: ZFilename, XFilename
      character(len = *), intent(in), optional :: YFilename

    end subroutine ZAXPY
  end interface

end module ControlSpaceAdvancer
