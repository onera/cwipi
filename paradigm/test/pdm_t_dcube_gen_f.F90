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

program testf

  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE  
  use mpi
#endif  
  use pdm_dcube_gen
  use iso_c_binding
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE  
  include "mpif.h"
#endif  

  integer :: code
  integer :: i_rank
  integer :: n_rank

  ! integer (kind = pdm_g_num_s), pointer :: block_distrib_index(:)
  ! type(c_ptr)                           :: cptr_block_distrib_index
  ! integer (kind = pdm_g_num_s), pointer :: gnum_elt(:)
  ! type(c_ptr), pointer                  :: cptr_gnum_elt(:)
  ! type(c_ptr)                           :: cptr_cptr_gnum_elt

  integer(c_int), parameter :: fComm = MPI_COMM_WORLD

  type(c_ptr)              :: dcube

  integer(pdm_g_num_s), parameter :: n_vtx_seg = 10
  double precision, parameter     :: length = 5.
  double precision, parameter     :: zero_x = 1.
  double precision, parameter     :: zero_y = 1.
  double precision, parameter     :: zero_z = 1.

  integer (c_int) :: n_face_group
  integer (c_int) :: dn_cell
  integer (c_int) :: dn_face
  integer (c_int) :: dn_vtx
  integer (c_int) :: sface_vtx
  integer (c_int) :: sface_group

  type(c_ptr) :: cptr_dface_cell
  type(c_ptr) :: cptr_dface_vtx_idx
  type(c_ptr) :: cptr_dface_vtx
  type(c_ptr) :: cptr_dvtx_coord
  type(c_ptr) :: cptr_dface_group_idx
  type(c_ptr) :: cptr_dface_group

  integer (kind = pdm_g_num_s), pointer :: dface_cell(:)

  integer (kind = pdm_l_num_s), pointer :: dface_vtx_idx(:)
  integer (kind = pdm_g_num_s), pointer :: dface_vtx(:)

  double precision, pointer :: dvtx_coord(:)

  integer (kind = pdm_l_num_s), pointer :: dface_group_idx(:)
  integer (kind = pdm_g_num_s), pointer :: dface_group(:)

  n_face_group = -1
  dn_cell = -1
  dn_face = -1
  dn_vtx = -1
  sface_vtx = -1
  sface_group = -1

  call mpi_init(code)
  call mpi_comm_rank(mpi_comm_world, i_rank, code)
  call mpi_comm_size(mpi_comm_world, n_rank, code)

  dcube = pdm_dcube_gen_init(fComm, n_vtx_seg, length, zero_x, zero_y, zero_z, PDM_OWNERSHIP_KEEP)

  call PDM_dcube_gen_dim_get (dcube,           &
                              n_face_group,    &
                              dn_cell,         &
                              dn_face,         &
                              dn_vtx,          &
                              sface_vtx,       &
                              sface_group)

  write(*,*) "n_face_group = ", n_face_group
  write(*,*) "dn_cell      = ", dn_cell
  write(*,*) "dn_face      = ", dn_face
  write(*,*) "dn_vtx       = ", dn_vtx
  write(*,*) "sface_vtx    = ", sface_vtx
  write(*,*) "sface_group  = ", sface_group

  call PDM_dcube_gen_data_get (dcube,                &
                               cptr_dface_cell,      &
                               cptr_dface_vtx_idx,   &
                               cptr_dface_vtx,       &
                               cptr_dvtx_coord,      &
                               cptr_dface_group_idx, &
                               cptr_dface_group)

  call c_f_pointer(cptr_dface_cell     , dface_cell     , [2*dn_face])

  write(*,*) dface_cell(1), dface_cell(2)

  call c_f_pointer(cptr_dface_vtx_idx  , dface_vtx_idx  , [dn_face+1])
  write(*,*) dface_vtx_idx(1), dface_vtx_idx(2)
  call c_f_pointer(cptr_dface_vtx      , dface_vtx      , [dface_vtx_idx(dn_face+1)])
  write(*,*) dface_vtx(1), dface_vtx(2)
  call c_f_pointer(cptr_dvtx_coord     , dvtx_coord     , [3*dn_vtx])
  write(*,*) dvtx_coord(1), dvtx_coord(2)

  call c_f_pointer(cptr_dface_group_idx, dface_group_idx, [n_face_group+1])
  write(*,*) "dface_group_idx(sface_group+1)  = ", dface_group_idx(n_face_group+1)
  call c_f_pointer(cptr_dface_group    , dface_group    , [dface_group_idx(n_face_group+1)])


  call pdm_dcube_gen_free(dcube)

  call mpi_finalize(code)

end program testf
