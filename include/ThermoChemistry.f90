#include "config.h"

module ThermoChemistry

  implicit none
  public

  interface

     subroutine getDensity(comm, species, density)

       integer, intent(in) :: comm
       character(len = *), intent(in) :: species

       real(SCALAR_KIND), intent(out) :: density

     end subroutine getDensity

  end interface

  interface

     subroutine getMolecularWeight(comm, species, molecularWeight)

       integer, intent(in) :: comm
       character(len = *), intent(in) :: species

       real(SCALAR_KIND), intent(out) :: molecularWeight

     end subroutine getMolecularWeight

  end interface

end module ThermoChemistry
