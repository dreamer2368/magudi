#include "config.h"

subroutine setupOptimizer(this, solver, region)

  ! <<< Derived types >>>
  use Optimizer_mod, only : t_Optimizer
  use Solver_mod, only : t_Solver
  use Region_mod, only : t_Region
  use Functional_mod, only : t_Functional

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this
  class(t_Solver) :: solver
  class(t_Region) :: region

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  integer :: i, bufferSize
  class(t_Functional), pointer :: functional => null()

  call this%cleanup()

  this%goldenRatio = 0.5_wp * ( 1.0_wp + sqrt(5.0_wp) )
  this%initialStep = getOption("optimization/initial_",0.1_wp)
  this%linminTol = getOption("optimization/linmin_tolerance",0.1_wp)
  this%cgTol = getOption("optimization/cg_tolerance",(1.0_wp)**(-5))

  this%numGrids = size(region%grids)
  assert(region%simulationFlags%enableFunctional)
  call solver%functionalFactory%connect(functional)
  assert(associated(functional))
  this%numParams = 0
  if (allocated(functional%params)) this%numParams = size(functional%params)

  allocate(this%base%states(this%numGrids))
  allocate(this%grad%states(this%numGrids))
  allocate(this%conjGrad%states(this%numGrids))
  allocate(this%prevGrad%states(this%numGrids))
  allocate(this%prevCG%states(this%numGrids))
  do i = 1, this%numGrids
    allocate(this%base%states(i)%conservedVariables,MOLD=region%states(i)%conservedVariables)
    allocate(this%grad%states(i)%conservedVariables,MOLD=region%states(i)%conservedVariables)
    allocate(this%conjGrad%states(i)%conservedVariables,MOLD=region%states(i)%conservedVariables)
    allocate(this%prevGrad%states(i)%conservedVariables,MOLD=region%states(i)%conservedVariables)
    allocate(this%prevCG%states(i)%conservedVariables,MOLD=region%states(i)%conservedVariables)
  end do
  allocate(this%base%params(this%numParams))
  allocate(this%grad%params(this%numParams))
  allocate(this%conjGrad%params(this%numParams))
  allocate(this%prevGrad%params(this%numParams))
  allocate(this%prevCG%params(this%numParams))

end subroutine setupOptimizer

subroutine cleanupOptimizer(this)

  ! <<< Derived types >>>
  use Optimizer_mod, only : t_Optimizer
  use RegionVector_mod, only : t_RegionVector

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this

  call this%base%cleanup()
  call this%grad%cleanup()
  call this%conjGrad%cleanup()
  call this%prevGrad%cleanup()
  call this%prevCG%cleanup()

end subroutine cleanupOptimizer

subroutine mnbrakConjugateGradient(this, solver, region)

  ! <<< Derived types >>>
  use Optimizer_mod, only : t_Optimizer
  use Solver_mod, only : t_Solver
  use Region_mod, only : t_Region

  implicit none

  ! <<< Arguments >>>
  class(t_Optimizer) :: this
  class(t_Solver) :: solver
  class(t_Region) :: region

  ! <<< Local variables >>>

end subroutine mnbrakConjugateGradient
