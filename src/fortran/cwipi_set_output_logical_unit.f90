!
!********************************************************************************
!
! Set output logical unit
!
! parameters
!    output_listing      <-- Output listing file
!
!********************************************************************************
!
subroutine cwipi_set_output_logical_unit (outputUnit)

  use cwipi_printfort

  implicit none

  integer :: outputUnit
  
  ifile = outputUnit

  call cwipi_set_output_listing_f(outputUnit)

end subroutine cwipi_set_output_logical_unit
