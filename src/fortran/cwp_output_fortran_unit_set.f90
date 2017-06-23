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
subroutine cwp_output_fortran_unit_set (outputUnit)

  use cwipi_printfort

  implicit none

  integer :: outputUnit
  
  ifile = outputUnit

  call cwipi_set_output_listing_f(outputUnit)

end subroutine cwp_output_fortran_unit_set
