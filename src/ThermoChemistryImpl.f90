#include "config.h"

module ThermoChemistryImpl

  implicit none
  public

end module ThermoChemistryImpl

subroutine getDensity(comm, species, density)

  ! <<< Internal modules >>>
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: species
  real(SCALAR_KIND), intent(out) :: density

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: message

  ! Get density at STP [kg/m^3]
  select case (trim(species))
  case('AIR')
     density = 1.205_wp
  case ('Ar', 'ARGON')
     density = 1.661_wp
  case ('CH4', 'METHANE')
     density = 0.668_wp
  case ('CO', 'CARBON_MONOXIDE')
     density = 1.165_wp
  case ('CO2', 'CARBON_DIOXIDE')
     density = 1.842_wp
  case ('H')
     density = 1.00794_wp
  case ('H2', 'HYDROGEN')
     density = 0.0899_wp
  case('H20', 'WATER')
     density = 18.01528_wp
  case ('N')
     density = 14.0067_wp
  case ('N2', 'NITROGEN')
     density = 28.0134_wp
  case ('O')
     density = 15.9994_wp
  case ('O2', 'OXYGEN')
     density = 31.9988_wp
  case default
     write(message, '(3A)') "Unknown species: ", trim(species), "!"
     call gracefulExit(comm, message)
  end select

end subroutine getDensity


subroutine getMolecularWeight(comm, species, molecularWeight)

  ! <<< Internal modules >>>
  use ErrorHandler, only : gracefulExit

  implicit none

  ! <<< Arguments >>>
  integer, intent(in) :: comm
  character(len = *), intent(in) :: species
  real(SCALAR_KIND), intent(out) :: molecularWeight

  ! <<< Local variables >>>
  integer, parameter :: wp = SCALAR_KIND
  character(len = STRING_LENGTH) :: message

  ! Get molar mass at STP [g/mol]
  select case (trim(species))
  case('AIR')
     molecularWeight = 28.851_wp
  case ('Ar', 'ARGON')
     molecularWeight = 39.948_wp
  case ('CH4', 'METHANE')
     molecularWeight = 16.04_wp
  case ('CO', 'CARBON_MONOXIDE')
     molecularWeight = 28.01_wp
  case ('CO2', 'CARBON_DIOXIDE')
     molecularWeight = 44.0095_wp
  case ('H')
     molecularWeight = 1.00794_wp
  case ('H2', 'HYDROGEN')
     molecularWeight = 2.01588_wp
  case('H20', 'WATER')
     molecularWeight = 18.01528_wp
  case ('N')
     molecularWeight = 14.0067_wp
  case ('N2', 'NITROGEN')
     molecularWeight = 28.0134_wp
  case ('O')
     molecularWeight = 15.9994_wp
  case ('O2' ,'OXYGEN')
     molecularWeight = 31.9988_wp
  case default
     write(message, '(3A)') "Unknown species: ", trim(species), "!"
     call gracefulExit(comm, message)
  end select

end subroutine getMolecularWeight
