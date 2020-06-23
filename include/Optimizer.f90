#include "config.h"

module Optimizer_enum

  implicit none
  public

  integer, parameter ::                                                                      &
       VERIFY = 301,                                                                         &
       NLCG   = 302,                                                                         &
       NEWTON = 303,                                                                         &
       GMRES  = 304,                                                                         &
       BICGSTAB = 305

end module Optimizer_enum

module Optimizer_mod

  use RegionVector_mod, only : t_RegionVector

  implicit none

  type, public :: t_Optimizer

     integer :: numGrids, numParams, saveInterval, reportInterval, nTimesteps
     type(t_RegionVector) :: residual, base, grad, conjGrad, prevGrad, prevCG,  &
                             x, r, rcg, Ap, s, As
     SCALAR_TYPE :: initialStep, goldenRatio, linminTol, cgTol
     SCALAR_TYPE :: bracket(3,2), gg1

     integer :: maxGMRES, maxRestart
     SCALAR_TYPE, allocatable :: H(:,:), sn(:), cs(:), beta(:)
     type(t_RegionVector), allocatable :: Q(:)
     character(len=STRING_LENGTH) :: outputPrefix
     logical :: verbose

   contains

     procedure, pass :: setup => setupOptimizer
     procedure, pass :: cleanup => cleanupOptimizer
     procedure, pass :: verifyAdjoint
     procedure, pass :: runNLCG
     ! procedure, pass :: runCGS
     procedure, pass :: runGMRES
     procedure, pass :: runBICGSTAB
     procedure, pass :: printBracket
     procedure, pass :: showProgress => showProgressOptimizer

  end type t_Optimizer

  interface

     subroutine setupOptimizer(this, region, mode)

       use Region_mod, only : t_Region
       import :: t_Optimizer

       implicit none

       class(t_Optimizer) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode

     end subroutine setupOptimizer

  end interface

  interface

     subroutine cleanupOptimizer(this)

       import :: t_Optimizer

       class(t_Optimizer) :: this

     end subroutine cleanupOptimizer

  end interface

  interface
     subroutine verifyAdjoint(this, region)

       use Region_mod, only : t_Region

       import :: t_Optimizer

       class(t_Optimizer) :: this
       class(t_Region) :: region

     end subroutine verifyAdjoint

  end interface

  interface

     subroutine runNLCG(this, region, restartFilename)

       use Region_mod, only : t_Region
       import :: t_Optimizer

       class(t_Optimizer) :: this
       class(t_Region) :: region
       character(len = *), intent(in), optional :: restartFilename

     end subroutine runNLCG

  end interface

  interface

     subroutine runGMRES(this, region, restartFilename)

       use Region_mod, only : t_Region
       import :: t_Optimizer

       class(t_Optimizer) :: this
       class(t_Region) :: region
       character(len = *), intent(in) :: restartFilename

     end subroutine runGMRES

  end interface

  interface

     subroutine runBICGSTAB(this, region, restartFilename)

       use Region_mod, only : t_Region
       import :: t_Optimizer

       class(t_Optimizer) :: this
       class(t_Region) :: region
       character(len = *), intent(in) :: restartFilename

     end subroutine runBICGSTAB

  end interface

  interface

     subroutine printBracket(this)

       import :: t_Optimizer

       implicit none

       class(t_Optimizer) :: this

     end subroutine printBracket

  end interface

  interface

     subroutine showProgressOptimizer(this, region, mode, step,                 &
                                      scalars, outputFilename, append)

       use Region_mod, only : t_Region
       import :: t_Optimizer

       implicit none

       class(t_Optimizer) :: this
       class(t_Region) :: region
       integer, intent(in) :: mode, step
       real(SCALAR_KIND), intent(in) :: scalars(:)
       character(len=*), intent(in) :: outputFilename
       logical, intent(in), optional :: append

     end subroutine showProgressOptimizer

  end interface

end module Optimizer_mod
