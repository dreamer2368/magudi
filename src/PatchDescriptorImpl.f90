#include "config.h"

subroutine parsePatchType(patchTypeString, patchType)

  ! <<< Derived types >>>
  use PatchDescriptor_type

  implicit none

  ! <<< Arguments >>>
  character(len = *), intent(in) :: patchTypeString
  integer, intent(out) :: patchType

  patchType = -1

  if (trim(patchTypeString) == "SPONGE") then
     patchType = SPONGE
  else if (trim(patchTypeString) == "ACTUATOR") then
     patchType = ACTUATOR
  else if (trim(patchTypeString) == "SOLENOIDAL_EXCITATION") then
     patchType = SOLENOIDAL_EXCITATION
  else if (trim(patchTypeString) == "SAT_FAR_FIELD") then
     patchType = SAT_FAR_FIELD
  else if (trim(patchTypeString) == "SAT_SLIP_WALL") then
     patchType = SAT_SLIP_WALL
  else if (trim(patchTypeString) == "SAT_ISOTHERMAL_WALL") then
     patchType = SAT_ISOTHERMAL_WALL
  else if (trim(patchTypeString) == "SAT_ADIABATIC_WALL") then
     patchType = SAT_ADIABATIC_WALL
  else if (trim(patchTypeString) == "SAT_BLOCK_INTERFACE") then
     patchType = SAT_BLOCK_INTERFACE
  end if

end subroutine parsePatchType

subroutine validatePatchDescriptor(this, globalGridSizes,                                    &
     simulationFlags, errorCode, message)

  ! <<< Derived types >>>
  use PatchDescriptor_type
  use SimulationFlags_type

  ! <<< Arguments >>>
  type(t_PatchDescriptor) :: this
  integer, intent(in) :: globalGridSizes(:,:)
  type(t_SimulationFlags), intent(in) :: simulationFlags
  integer, intent(out) :: errorCode
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Local variables >>>
  integer :: i, extent(6)
  logical :: isExtentValid

  extent = (/ this%iMin, this%iMax,                                                          &
       this%jMin, this%jMax,                                                                 &
       this%kMin, this%kMax /)

  ! Check if the grid index is valid.
  if (this%gridIndex <= 0 .or. this%gridIndex > size(globalGridSizes, 2)) then
     write(message, '(3A,I0.0,A)') "Patch '", trim(this%name),                               &
          "' has an invalid grid index: ", this%gridIndex, "!"
     errorCode = 2
     return
  end if

  ! Check if the normal direction is valid.
  select case (this%patchType)
  case (ACTUATOR, SOLENOIDAL_EXCITATION)
  case default
     if (abs(this%normalDirection) > size(globalGridSizes, 1) .or.                           &
          this%normalDirection == 0) then
        write(message, '(3A,2(I0.0,A))') "Patch '", trim(this%name), "' on grid ",           &
             this%gridIndex, " has an invalid normal direction: ",                           &
             this%normalDirection, "!"
        errorCode = 2
        return
     end if
  end select

  ! Extents must be non-zero.
  isExtentValid = all(extent /= 0)

  ! Negative values indicate counting backwards from the end.
  do i = 1, size(globalGridSizes, 1)
     if (extent(1+2*(i-1)) < 0) extent(1+2*(i-1)) = extent(1+2*(i-1))                        &
          + globalGridSizes(i, this%gridIndex) + 1
     if (extent(2+2*(i-1)) < 0) extent(2+2*(i-1)) = extent(2+2*(i-1))                        &
          + globalGridSizes(i, this%gridIndex) + 1
  end do

  ! Check that extent describes a part of the grid.
  do i = 1, size(globalGridSizes, 1)
     isExtentValid = isExtentValid .and. (extent(2+2*(i-1)) >= extent(1+2*(i-1)))
     isExtentValid = isExtentValid .and. (extent(2+2*(i-1))                                  &
          - extent(1+2*(i-1)) + 1 <= globalGridSizes(i, this%gridIndex))
  end do
  do while (i <= 3) !... reset for direction > number of dimensions.
     extent(1+2*(i-1)) = 1
     extent(2+2*(i-1)) = 1
     i = i + 1
  end do

  ! Fail if the extent is invalid.
  if (.not. isExtentValid) then
     write(message, '(3A,I0.0,A,6(I0.0,1X))') "Patch '", trim(this%name), "' on grid ",      &
          this%gridIndex, " has an invalid extent: ",                                        &
          this%iMin, this%iMax, this%jMin, this%jMax, this%kMin, this%kMax, "!"
     errorCode = 2
     return
  end if

  ! Check if patch spans more than one grid point along normal direction.
  i = abs(this%normalDirection)
  select case (this%patchType)
  case (SPONGE, ACTUATOR, SOLENOIDAL_EXCITATION)
  case default
     if (extent(1+2*(i-1)) /= extent(2+2*(i-1))) then
        write(message, '(3A,I0.0,A)') "Patch '", trim(this%name), "' on grid ",              &
             this%gridIndex, " is not allowed to extend across more &
             &than 1 grid point along the normal direction!"
        errorCode = 2
        return
     end if
  end select

  ! Check if sponge zones are aligned with computational boundaries.
  if (this%patchType == SPONGE) then
     i = abs(this%normalDirection)
     if ((this%normalDirection > 0 .and. .not. extent(1+2*(i-1)) == 1) .or.                  &
          (this%normalDirection < 0 .and.                                                    &
          .not. extent(2+2*(i-1)) == globalGridSizes(i, this%gridIndex))) then
        write(message, '(2(A,I0.0),A)') "Sponge patch '", trim(this%name), "' on grid ",     &
             this%gridIndex, " is not aligned with a computational boundary!"
        errorCode = 2
        return
     end if
  end if

  ! Check if target state is being used if the patch requires it.
  select case (this%patchType)
  case (SPONGE, SAT_FAR_FIELD)
     if (.not. simulationFlags%useTargetState) then
        write(message, '(3A,I0.0,A)') "Not using patch '", trim(this%name), "' on grid ",    &
             this%gridIndex, " because there is no target state!"
        errorCode = 1
        return
     end if
  end select

  ! Check if domain is two-dimensional for using patch `SOLENOIDAL_EXCITATION`.
  select case (this%patchType)
  case (SOLENOIDAL_EXCITATION)
     if (size(globalGridSizes, 1) /= 2) then
        write(message, '(3A,I0.0,A)') "Not using patch '", trim(this%name), "' on grid ",    &
             this%gridIndex, " because the domain is not two-dimensional!"
        errorCode = 1
        return
     end if
  end select

  ! Check if an `ACTUATOR` patch will be used.
  select case (this%patchType)
  case (ACTUATOR)
     if (simulationFlags%predictionOnly) then
        write(message, '(3A,I0.0,A)') "Not using patch '", trim(this%name), "' on grid ",    &
             this%gridIndex, " because control actuation has been disabled!"
        errorCode = 1
        return
     end if
  end select

  ! Copy over extent to iMin, iMax, etc.
  this%iMin = extent(1)
  this%iMax = extent(2)
  this%jMin = extent(3)
  this%jMax = extent(4)
  this%kMin = extent(5)
  this%kMax = extent(6)

  errorCode = 0

end subroutine validatePatchDescriptor

subroutine validatePatchesConnectivity(patchDescriptors, errorCode, message)

  ! <<< Derived types >>>
  use PatchDescriptor_type

  ! <<< Internal modules >>>
  use InputHelper, only : getOption

  implicit none

  ! <<< Arguments >>>
  type(t_PatchDescriptor), intent(in) :: patchDescriptors(:)
  integer, intent(out) :: errorCode
  character(len = STRING_LENGTH), intent(out) :: message

  ! <<< Local variables >>>
  integer :: i, j
  character(len = STRING_LENGTH) :: key, val
  logical, allocatable :: foundPair(:)

  errorCode = 0

  allocate(foundPair(size(patchDescriptors)), source = .false.)

  do i = 1, size(patchDescriptors)
     select case (patchDescriptors(i)%patchType)
     case (SAT_BLOCK_INTERFACE)

        write(key, '(3A)') "patches/", trim(patchDescriptors(i)%name), "/conforms_with"
        val = getOption(key, "")

        if (len_trim(val) > 0) then

           do j = 1, size(patchDescriptors)
              if (i == j) cycle
              if (trim(patchDescriptors(j)%name) == trim(val)) then

                 if (patchDescriptors(j)%patchType /= SAT_BLOCK_INTERFACE) then
                    write(message, '(2(3A,I0.0),A)') "Patch '",                              &
                         trim(patchDescriptors(i)%name), "' on grid ",                       &
                         patchDescriptors(i)%gridIndex, " conforms with patch '",            &
                         trim(patchDescriptors(j)%name), "' on grid ",                       &
                         patchDescriptors(j)%gridIndex, ", which is not an interface patch!"
                    errorCode = 2
                    return
                 end if

                 if (patchDescriptors(j)%normalDirection /=                                  &
                      -patchDescriptors(i)%normalDirection) then
                    write(message, '(2(3A,I0.0),A)') "Patch '",                              &
                         trim(patchDescriptors(i)%name), "' on grid ",                       &
                         patchDescriptors(i)%gridIndex, " conforms with patch '",            &
                         trim(patchDescriptors(j)%name), "' on grid ",                       &
                         patchDescriptors(j)%gridIndex, &
                         ", but the dot product of their unit normals is not -1!"
                    errorCode = 2
                    return
                 end if

                 if (patchDescriptors(j)%gridIndex == patchDescriptors(i)%gridIndex) then
                    write(message, '(2(3A,I0.0),A)') "Patch '",                              &
                         trim(patchDescriptors(i)%name), "' on grid ",                       &
                         patchDescriptors(i)%gridIndex, " conforms with patch '",            &
                         trim(patchDescriptors(j)%name), "' also on grid ",                  &
                         patchDescriptors(j)%gridIndex, "!"
                    errorCode = 2
                    return
                 end if

                 if ((patchDescriptors(j)%iMax - patchDescriptors(j)%iMin /=                 &
                      patchDescriptors(i)%iMax - patchDescriptors(i)%iMin) .or.              &
                      (patchDescriptors(j)%jMax - patchDescriptors(j)%jMin /=                &
                      patchDescriptors(i)%jMax - patchDescriptors(i)%jMin) .or.              &
                      (patchDescriptors(j)%kMax - patchDescriptors(j)%kMin /=                &
                      patchDescriptors(i)%kMax - patchDescriptors(i)%kMin)) then
                    write(message, '(2(3A,I0.0),A)') "Patch '",                              &
                         trim(patchDescriptors(i)%name), "' on grid ",                       &
                         patchDescriptors(i)%gridIndex, " does not conform with patch '",    &
                         trim(patchDescriptors(j)%name), "' on grid ",                       &
                         patchDescriptors(j)%gridIndex, "!"
                    errorCode = 2
                    return
                 end if

                 foundPair(i) = .true.
                 foundPair(j) = .true.

              end if
           end do

        end if

     end select
  end do

  do i = 1, size(patchDescriptors)
     select case (patchDescriptors(i)%patchType)
     case (SAT_BLOCK_INTERFACE)
        if (.not. foundPair(i)) then
           write(message, '(3A,I0.0,A)')                                                     &
                "Unable to find a patch that conforms with interface patch '",               &
                trim(patchDescriptors(i)%name), "' on grid ",                                &
                patchDescriptors(i)%gridIndex, "!"
           errorCode = 2
           return
        end if
     end select
  end do

end subroutine validatePatchesConnectivity
