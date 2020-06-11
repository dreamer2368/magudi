#include "config.h"

module Optimizer_mod

  use RegionVector_mod, only : t_RegionVector

  implicit none

  type, public :: t_Optimizer

     integer :: numGrids, numParams
     type(t_RegionVector) :: residual, base, grad, conjGrad, prevGrad, prevCG
     SCALAR_TYPE :: initialStep, goldenRatio, linminTol, cgTol
     SCALAR_TYPE :: bracket(3), step

   contains

     procedure, pass :: setup => setupOptimizer
     procedure, pass :: cleanup => cleanupOptimizer
     procedure, pass :: verifyNLCG
     procedure, pass :: mnbrak => mnbrakConjugateGradient
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
     subroutine verifyNLCG(this, region)

       use Region_mod, only : t_Region

       import :: t_Optimizer

       class(t_Optimizer) :: this
       class(t_Region) :: region

     end subroutine verifyNLCG

  end interface

  interface

     subroutine mnbrakConjugateGradient(this, solver, region)

       use Solver_mod, only : t_Solver
       use Region_mod, only : t_Region
       import :: t_Optimizer

       class(t_Optimizer) :: this
       class(t_Solver) :: solver
       class(t_Region) :: region

     end subroutine mnbrakConjugateGradient

  end interface

  ! interface
  !
  !    subroutine linminConjugateGradient(this, solver, region)
  !
  !      import :: t_Optimizer
  !      use Solver_mod, only : t_Solver
  !      use Region_mod, only : t_Region
  !
  !      class(t_Optimizer) :: this
  !      class(t_Solver) :: solver
  !      class(t_Region) :: region
  !
  !    end subroutine linminConjugateGradient
  !
  ! end interface
  !
  ! interface
  !
  !    subroutine frprmnConjugateGradient(this, solver, region)
  !
  !      import :: t_Optimizer
  !      use Solver_mod, only : t_Solver
  !      use Region_mod, only : t_Region
  !
  !      class(t_Optimizer) :: this
  !      class(t_Solver) :: solver
  !      class(t_Region) :: region
  !
  !    end subroutine frprmnConjugateGradient
  !
  ! end interface

end module Optimizer_mod
