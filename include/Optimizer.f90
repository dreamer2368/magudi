#include "config.h"

module Optimier_enum

  implicit none
  public

  integer, parameter ::                                                                      &
       VERIFY = 301,                                                                         &
       NLCG   = 302,                                                                         &
       NEWTON = 303

end module Optimier_enum

module Optimizer_mod

  use RegionVector_mod, only : t_RegionVector

  implicit none

  type, public :: t_Optimizer

     integer :: numGrids, numParams, saveInterval, reportInterval, nTimesteps
     type(t_RegionVector) :: residual, base, grad, conjGrad, prevGrad, prevCG
     SCALAR_TYPE :: initialStep, goldenRatio, linminTol, cgTol
     SCALAR_TYPE :: bracket(3,2), gg1
     character(len=STRING_LENGTH) :: outputPrefix
     logical :: verbose

   contains

     procedure, pass :: setup => setupOptimizer
     procedure, pass :: cleanup => cleanupOptimizer
     procedure, pass :: verifyAdjoint
     procedure, pass :: runNLCG
     procedure, pass :: printBracket
     procedure, pass :: showProgress => showProgressOptimizer
     ! procedure, pass :: linmin => linminConjugateGradient
     ! procedure, pass :: frprmn => frprmnConjugateGradient

  end type t_Optimizer

  interface

     subroutine setupOptimizer(this, region)

       use Region_mod, only : t_Region
       import :: t_Optimizer

       implicit none

       class(t_Optimizer) :: this
       class(t_Region) :: region

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

     subroutine printBracket(this)

       import :: t_Optimizer

       implicit none

       class(t_Optimizer) :: this

     end subroutine printBracket

  end interface

  interface

     subroutine showProgressOptimizer(this, region, costFunctional,             &
                                     costSensitivity, outputFilename, append)

       use Region_mod, only : t_Region
       import :: t_Optimizer

       implicit none

       class(t_Optimizer) :: this
       class(t_Region) :: region
       real(SCALAR_KIND), intent(in) :: costFunctional, costSensitivity
       character(len=*), intent(in) :: outputFilename
       logical, intent(in) :: append

     end subroutine showProgressOptimizer

  end interface

end module Optimizer_mod
