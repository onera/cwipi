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

module pdm_gnum

  use pdm

  implicit none

  interface

    function PDM_gnum_create (dim,         &
                              n_part,      &
                              merge,       &
                              tolerance,   &
                              fcomm,       &
                              owner)       &
                              result (gen_gnum) &
      bind (c, name = 'PDM_gnum_create_cf')

      use iso_c_binding
      implicit none

      integer(c_int), value :: dim
      integer(c_int), value :: n_part
      integer(c_int), value :: merge
      real(c_double), value :: tolerance
      integer(c_int), value :: fComm
      integer(c_int), value :: owner

      type (c_ptr) :: gen_gnum

    end function PDM_gnum_create

    subroutine PDM_gnum_set_from_coords (gen_gnum,   &
                                         i_part,     &
                                         n_elts,     &
                                         coords,     &
                                         char_length)&
      bind (c, name = 'PDM_gnum_set_from_coords')

      use iso_c_binding
      implicit none

      type(c_ptr), value :: gen_gnum

      integer(c_int), value :: i_part
      integer(c_int), value :: n_elts
      type(c_ptr), value :: coords
      type(c_ptr), value :: char_length

    end subroutine PDM_gnum_set_from_coords

    subroutine PDM_gnum_set_from_parents (gen_gnum,   &
                                         i_part,     &
                                         n_elts,     &
                                         parent_gnum)&
      bind (c, name = 'PDM_gnum_set_from_parents')

      use iso_c_binding
      implicit none

      type(c_ptr), value :: gen_gnum

      integer(c_int), value :: i_part
      integer(c_int), value :: n_elts
      type(c_ptr), value :: parent_gnum

    end subroutine PDM_gnum_set_from_parents

    subroutine PDM_gnum_compute (gen_gnum) &
      bind (c, name = 'PDM_gnum_compute')

      use iso_c_binding
      implicit none

      type(c_ptr), value :: gen_gnum

    end subroutine PDM_gnum_compute

    function PDM_gnum_get (gen_gnum, i_part) &
      result (g_nums) &
      bind (c, name = 'PDM_gnum_get')

      use iso_c_binding
      implicit none

      type(c_ptr) :: g_nums
      type(c_ptr), value :: gen_gnum
      integer(c_int), value :: i_part

    end function PDM_gnum_get

    subroutine PDM_gnum_free (gen_gnum)     &
      bind (c, name = 'PDM_gnum_free')
      use iso_c_binding
      implicit none

      type (c_ptr)  , value :: gen_gnum

    end subroutine PDM_gnum_free

  end interface

end module pdm_gnum
