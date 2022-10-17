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


module pdm_sphere_surf_gen

  use pdm
  use pdm_pointer_array

  implicit none

  contains


  subroutine PDM_sphere_surf_icosphere_gen_part(comm,           &
                                                n,              &
                                                x_center,       &
                                                y_center,       &
                                                z_center,       &
                                                radius,         &
                                                n_part,         &
                                                part_method,    &
                                                pn_vtx,         &
                                                pvtx_coord,     &
                                                pvtx_ln_to_gn,  &
                                                pn_face,        &
                                                pface_vtx_idx,  &
                                                pface_vtx,      &
                                                pface_ln_to_gn)
    use iso_c_binding
    implicit none

    integer,              intent(in)  :: comm
    integer(pdm_g_num_s), intent(in)  :: n
    double precision,     intent(in)  :: x_center
    double precision,     intent(in)  :: y_center
    double precision,     intent(in)  :: z_center
    double precision,     intent(in)  :: radius
    integer,              intent(in)  :: n_part
    integer,              intent(in)  :: part_method
    integer(pdm_l_num_s), pointer     :: pn_vtx(:)
    type(PDM_pointer_array_t)         :: pvtx_coord
    type(PDM_pointer_array_t)         :: pvtx_ln_to_gn
    integer(pdm_l_num_s), pointer     :: pn_face(:)
    type(PDM_pointer_array_t)         :: pface_vtx_idx
    type(PDM_pointer_array_t)         :: pface_vtx
    type(PDM_pointer_array_t)         :: pface_ln_to_gn

    integer(c_int)                    :: c_comm
    type(c_ptr)                       :: c_pn_vtx            = C_NULL_PTR
    type(c_ptr)                       :: c_pn_face           = C_NULL_PTR
    type(c_ptr)                       :: cc_pvtx_coord       = C_NULL_PTR
    type(c_ptr)                       :: cc_pvtx_ln_to_gn    = C_NULL_PTR
    type(c_ptr)                       :: cc_pface_vtx_idx    = C_NULL_PTR
    type(c_ptr)                       :: cc_pface_vtx        = C_NULL_PTR
    type(c_ptr)                       :: cc_pface_ln_to_gn   = C_NULL_PTR
    type(c_ptr), pointer              :: c_pvtx_coord(:)     => null()
    type(c_ptr), pointer              :: c_pvtx_ln_to_gn(:)  => null()
    type(c_ptr), pointer              :: c_pface_vtx_idx(:)  => null()
    type(c_ptr), pointer              :: c_pface_vtx(:)      => null()
    type(c_ptr), pointer              :: c_pface_ln_to_gn(:) => null()
    integer                           :: ipart

    interface

      subroutine PDM_sphere_surf_icosphere_gen_part_c(comm,           &
                                                      n,              &
                                                      x_center,       &
                                                      y_center,       &
                                                      z_center,       &
                                                      radius,         &
                                                      n_part,         &
                                                      part_method,    &
                                                      pn_vtx,         &
                                                      pvtx_coord,     &
                                                      pvtx_ln_to_gn,  &
                                                      pn_face,        &
                                                      pface_vtx_idx,  &
                                                      pface_vtx,      &
                                                      pface_ln_to_gn) &
      bind (c, name='PDM_sphere_surf_icosphere_gen_part')

        use iso_c_binding
        implicit none

        integer(c_int),  value :: comm
#ifdef PDM_LONG_G_NUM
        integer(c_long), value :: n
#else
        integer(c_int),  value :: n
#endif
        real(c_double),  value :: x_center
        real(c_double),  value :: y_center
        real(c_double),  value :: z_center
        real(c_double),  value :: radius
        integer(c_int),  value :: n_part
        integer(c_int),  value :: part_method
        type(c_ptr)            :: pn_vtx
        type(c_ptr)            :: pvtx_coord
        type(c_ptr)            :: pvtx_ln_to_gn
        type(c_ptr)            :: pn_face
        type(c_ptr)            :: pface_vtx_idx
        type(c_ptr)            :: pface_vtx
        type(c_ptr)            :: pface_ln_to_gn

      end subroutine PDM_sphere_surf_icosphere_gen_part_c
    end interface

    c_comm = PDM_MPI_Comm_f2c(comm)

    call PDM_sphere_surf_icosphere_gen_part_c(c_comm,            &
                                              n,                 &
                                              x_center,          &
                                              y_center,          &
                                              z_center,          &
                                              radius,            &
                                              n_part,            &
                                              part_method,       &
                                              c_pn_vtx,          &
                                              cc_pvtx_coord,     &
                                              cc_pvtx_ln_to_gn,  &
                                              c_pn_face,         &
                                              cc_pface_vtx_idx,  &
                                              cc_pface_vtx,      &
                                              cc_pface_ln_to_gn)

    call c_f_pointer(c_pn_vtx, &
                     pn_vtx,   &
                     [n_part])

    call c_f_pointer(c_pn_face, &
                     pn_face,   &
                     [n_part])


    call PDM_pointer_array_create(pvtx_coord,     &
                                  n_part,         &
                                  PDM_TYPE_DOUBLE)

    call PDM_pointer_array_create(pvtx_ln_to_gn,  &
                                  n_part,         &
                                  PDM_TYPE_G_NUM)

    call PDM_pointer_array_create(pface_vtx_idx,  &
                                  n_part,         &
                                  PDM_TYPE_INT)

    call PDM_pointer_array_create(pface_vtx,      &
                                  n_part,         &
                                  PDM_TYPE_INT)

    call PDM_pointer_array_create(pface_ln_to_gn, &
                                  n_part,         &
                                  PDM_TYPE_G_NUM)

    call c_f_pointer(cc_pvtx_coord, &
                     c_pvtx_coord,  &
                     [n_part])

    call c_f_pointer(cc_pvtx_ln_to_gn, &
                     c_pvtx_ln_to_gn,  &
                     [n_part])

    call c_f_pointer(cc_pface_vtx_idx, &
                     c_pface_vtx_idx,  &
                     [n_part])

    call c_f_pointer(cc_pface_vtx, &
                     c_pface_vtx,  &
                     [n_part])

    call c_f_pointer(cc_pface_ln_to_gn, &
                     c_pface_ln_to_gn,  &
                     [n_part])

    do ipart = 1, n_part
      call PDM_pointer_array_part_set_from_cptr(pvtx_coord,          &
                                                ipart-1,             &
                                                c_pvtx_coord(ipart), &
                                                3*pn_vtx(ipart))

      call PDM_pointer_array_part_set_from_cptr(pvtx_ln_to_gn,          &
                                                ipart-1,                &
                                                c_pvtx_ln_to_gn(ipart), &
                                                pn_vtx(ipart))

      call PDM_pointer_array_part_set_from_cptr(pface_vtx_idx,          &
                                                ipart-1,                &
                                                c_pface_vtx_idx(ipart), &
                                                pn_face(ipart))

      call PDM_pointer_array_part_set_from_cptr(pface_vtx,          &
                                                ipart-1,            &
                                                c_pface_vtx(ipart), &
                                                pn_face(ipart))

      call PDM_pointer_array_part_set_from_cptr(pface_ln_to_gn,          &
                                                ipart-1,                 &
                                                c_pface_ln_to_gn(ipart), &
                                                pn_face(ipart))
    end do

  end subroutine PDM_sphere_surf_icosphere_gen_part

end module pdm_sphere_surf_gen
