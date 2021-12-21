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

module pdm_dcube_gen

  use pdm

  implicit none

  interface

    function PDM_dcube_gen_init (fcomm,          &
                                  n_vtx_seg,     &
                                  length,        &
                                  zero_x,        &
                                  zero_y,        &
                                  zero_z,        &
                                  owner)         &
                                  result (dcube) &
      bind (c, name = 'PDM_dcube_gen_init_cf')

      use iso_c_binding
      implicit none

      integer(c_int), value :: fComm
      integer(c_int), value :: owner

#ifdef PDM_LONG_G_NUM
      integer (c_long),  value :: n_vtx_seg
#else
      integer (c_int),  value :: n_vtx_seg
#endif
      real(c_double), value :: length
      real(c_double), value :: zero_x, zero_y, zero_z

      type (c_ptr) :: dcube

    end function PDM_dcube_gen_init

    subroutine PDM_dcube_gen_dim_get (dcube,           &
                                      n_face_group,    &
                                      dn_cell,         &
                                      dn_face,         &
                                      dn_vtx,          &
                                      sface_vtx,       &
                                      sface_group)     &
      bind (c, name = 'PDM_dcube_gen_dim_get')

      use iso_c_binding
      implicit none

      type(c_ptr), value :: dcube

      integer(c_int) :: n_face_group
      integer(c_int) :: dn_cell
      integer(c_int) :: dn_face
      integer(c_int) :: dn_vtx
      integer(c_int) :: sface_vtx
      integer(c_int) :: sface_group

    end subroutine PDM_dcube_gen_dim_get

    subroutine PDM_dcube_gen_data_get (dcube,          &
                                       dface_cell,      &
                                       dface_vtx_idx,   &
                                       dface_vtx,       &
                                       dvtx_coord,      &
                                       dface_group_idx, &
                                       dface_group)     &
      bind (c, name = 'PDM_dcube_gen_data_get')

      use iso_c_binding
      implicit none

      type   (c_ptr), value :: dcube
      type   (c_ptr) :: dface_cell
      type   (c_ptr) :: dface_vtx_idx
      type   (c_ptr) :: dface_vtx
      type   (c_ptr) :: dvtx_coord
      type   (c_ptr) :: dface_group_idx
      type   (c_ptr) :: dface_group

    end subroutine PDM_dcube_gen_data_get

    subroutine PDM_dcube_gen_free (dcube)     &
      bind (c, name = 'PDM_dcube_gen_free')
      use iso_c_binding
      implicit none

      type (c_ptr)  , value :: dcube

    end subroutine PDM_dcube_gen_free

  end interface

end module pdm_dcube_gen
