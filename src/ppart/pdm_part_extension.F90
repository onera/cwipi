!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2020  ONERA
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 3 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------

#include "pdm_configf.h"

module pdm_part_extension

  use pdm

  implicit none

  integer, parameter :: PDM_EXTEND_FROM_FACE = 0
  integer, parameter :: PDM_EXTEND_FROM_EDGE = 1
  integer, parameter :: PDM_EXTEND_FROM_VTX  = 2

interface

!>
!!
!! \brief Compute a part extension structure
!!
!! \param [in]   part_ext          PDM_part_extension_t
!!
!!

subroutine PDM_part_extension_compute (part_ext) &
bind (c, name='PDM_part_extension_compute')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: part_ext

end subroutine PDM_part_extension_compute



!>
!!
!! \brief Free a part extension structure
!!
!! \param [in]   part_ext          PDM_part_extension_t
!!
!!

subroutine PDM_part_extension_free (part_ext) &
bind (c, name='PDM_part_extension_free')
  use iso_c_binding
  implicit none

  type(c_ptr), value :: part_ext

end subroutine PDM_part_extension_free

end interface


contains


!>
!!
!! \brief Initialize a part extension structure
!!
!! \param [out]   part_ext          Pointer to a new \ref PDM_part_extension_t object
!! \param [in]    n_part            Number of partitions
!! \param [in]    extend_type       Type of extension
!! \param [in]    depth             Depth of extension
!! \param [in]    comm              MPI communicator
!! \param [in]    owner             Ownership
!!
!!

subroutine PDM_part_extension_create (part_ext,    &
                                      n_part,      &
                                      extend_type, &
                                      depth,       &
                                      comm,        &
                                      owner)
  use iso_c_binding
  implicit none

  type(c_ptr)         :: part_ext
  integer, intent(in) :: n_part
  integer, intent(in) :: extend_type
  integer, intent(in) :: depth
  integer, intent(in) :: comm
  integer, intent(in) :: owner

  integer(c_int)      :: c_comm

  interface
    function PDM_part_extension_create_c (n_part,      &
                                          extend_type, &
                                          depth,       &
                                          comm,        &
                                          owner)       &
    result (part_ext)                                  &
    bind (c, name='PDM_part_extension_create')
      use iso_c_binding
      implicit none

      integer(c_int), value :: n_part
      integer(c_int), value :: extend_type
      integer(c_int), value :: depth
      integer(c_int), value :: comm
      integer(c_int), value :: owner
      type(c_ptr)           :: part_ext

    end function PDM_part_extension_create_c
  end interface

  c_comm = PDM_MPI_Comm_f2c(comm)

  part_ext = PDM_part_extension_create_c (n_part,      &
                                          extend_type, &
                                          depth,       &
                                          comm,        &
                                          owner)

end subroutine PDM_part_extension_create



end module pdm_part_extension
