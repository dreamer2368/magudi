#include "config.h"

module EnsightHelper

  implicit none
  public

  type, public :: t_Ensight

     integer :: nOutputTimes, fileview, dataSize, gDataSize
     real(SCALAR_KIND), allocatable :: outputTimes(:)
     real(KIND=4), allocatable :: buffer1_sp(:), buffer3_sp(:,:)
     character(len = STRING_LENGTH) :: directory, filename

   contains

     procedure, public, pass :: setup => setupEnsight
     procedure, public, pass :: output => outputEnsight

  end type t_Ensight

  interface

     subroutine setupEnsight(this, grid, time)

       !> Sets up the EnSight files, including timing and file types.

       use Grid_mod, only : t_Grid

       import :: t_Ensight

       class(t_Ensight) :: this
       class(t_Grid) :: grid
       real(SCALAR_KIND), intent(in) :: time

     end subroutine setupEnsight

  end interface

  interface

     subroutine outputEnsight(this, state, grid, mode, time, nSpecies)

       !> Writes the primitive variables (density, velocity, temperature, mass fraction)
       !> or adjoint variables in EnSight Gold format for visualization purposes.

       use State_mod, only : t_State
       use Grid_mod, only : t_Grid

       import :: t_Ensight

       class(t_Ensight) :: this
       class(t_State) :: state
       class(t_Grid) :: grid
       integer, intent(in) :: mode
       integer, intent(in), optional :: nSpecies
       real(SCALAR_KIND), intent(in) :: time

     end subroutine outputEnsight

  end interface

end module EnsightHelper
