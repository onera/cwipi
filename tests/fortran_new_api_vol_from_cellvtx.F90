!-----------------------------------------------------------------------------
! This file is part of the CWIPI library.
!
! Copyright (C) 2024  ONERA
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
#include "cwipi_configf.h"

program fortran_new_api_vol_from_cellvtx

#ifdef CWP_HAVE_FORTRAN_MPI_MODULE
    use mpi
#endif
    use cwp

    implicit none

#ifndef CWP_HAVE_FORTRAN_MPI_MODULE
    include "mpif.h"
#endif

  !--------------------------------------------------------------------
  integer, parameter          :: comm = MPI_COMM_WORLD
  integer                     :: i_rank, err
  logical                     :: is_code1
  character(len = 5), pointer :: code_names(:)         => null()
  character(len = 5), pointer :: coupled_code_names(:) => null()
  integer,            pointer :: intra_comms(:)        => null()
  character(len = 99)         :: coupling_name

  integer                     :: n_vtx
  double precision,   pointer :: vtx_coord(:,:) => null()
  integer(c_long),    pointer :: vtx_gnum(:)    => null()
  integer                     :: n_cell
  integer(c_int),     pointer :: cell_vtx_idx(:) => null()
  integer(c_int),     pointer :: cell_vtx(:)     => null()
  integer(c_long),    pointer :: cell_gnum(:)    => null()

  character(len = 99)         :: field_name
  double precision,   pointer :: send_val(:) => null()
  double precision,   pointer :: recv_val(:) => null()

  integer                     :: i, j, k, i_vtx
  double precision            :: x, max_err
  !--------------------------------------------------------------------


  ! Initialize MPI
  call mpi_init(err)
  call mpi_comm_rank(comm, i_rank, err)

  is_code1 = (mod(i_rank, 2) == 0)

  allocate(code_names(1),         &
           coupled_code_names(1), &
           intra_comms(1))

  if (is_code1) then
    code_names        (1) = "code1"
    coupled_code_names(1) = "code2"
  else
    code_names        (1) = "code2"
    coupled_code_names(1) = "code1"
  endif

  ! Initialize CWIPI
  call cwp_init(comm,          &
                1,             &
                code_names,    &
                CWP_STATUS_ON, &
                intra_comms)

  ! Create coupling environment
  coupling_name = "fortran_new_api_vol_from_cellvtx"

  call cwp_cpl_create(code_names(1),                                         &
                      coupling_name,                                         &
                      coupled_code_names(1),                                 &
                      CWP_INTERFACE_VOLUME,                                  &
                      CWP_COMM_PAR_WITH_PART,                                &
                      CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, &
                      1,                                                     &
                      CWP_DYNAMIC_MESH_STATIC,                               &
                      CWP_TIME_EXCH_USER_CONTROLLED)

  call cwp_visu_set(code_names(1),           &
                    coupling_name,           &
                    1,                       &
                    CWP_VISU_FORMAT_ENSIGHT, &
                    "text")

  ! Define interface mesh
  n_vtx = 27
  allocate(vtx_coord(3, n_vtx))
  i_vtx = 1
  do k = 0,2
    do j = 0,2
      do i = 0,2
        vtx_coord(1,i_vtx) = i - 1
        vtx_coord(2,i_vtx) = j - 1
        vtx_coord(3,i_vtx) = k - 1
        i_vtx = i_vtx + 1
      enddo
    enddo
  enddo

  if (is_code1) then
    ! Rotate 90° around z-axis
    do i_vtx = 1,n_vtx
      x = vtx_coord(1,i_vtx)
      vtx_coord(1,i_vtx) = vtx_coord(2,i_vtx)
      vtx_coord(2,i_vtx) = x
    enddo
  endif

  n_cell = 22
  allocate(cell_vtx_idx(n_cell+1), cell_vtx(110))
  cell_vtx_idx = [0, 8, 16, 22, 28, 34, 40, 45, 50, 55, 60, 65, 70, 74, 78, 82, 86, 90, 94, 98, 102, 106, 110]
  cell_vtx     = [ &
    1, 2, 5, 4, 10, 11, 14, 13, &
    4, 5, 8, 7, 13, 14, 17, 16, &
    2, 3, 5, 11, 12, 14,        &
    3, 6, 5, 12, 15, 14,        &
    5, 6, 9, 14, 15, 18,        &
    5, 9, 8, 14, 18, 17,        &
    10, 11, 14, 13, 20,         &
    13, 14, 23, 22, 20,         &
    10, 13, 22, 19, 20,         &
    13, 14, 17, 16, 26,         &
    13, 22, 23, 14, 26,         &
    13, 16, 25, 22, 26,         &
    11, 12, 14, 20,             &
    14, 23, 20, 24,             &
    14, 15, 24, 12,             &
    14, 12, 24, 20,             &
    12, 24, 20, 21,             &
    14, 15, 18, 24,             &
    14, 23, 24, 26,             &
    14, 18, 17, 26,             &
    14, 18, 26, 24,             &
    18, 26, 24, 27              &
  ]

  call cwp_mesh_interf_vtx_set(code_names(1), &
                               coupling_name, &
                               0,             &
                               n_vtx,         &
                               vtx_coord,     &
                               vtx_gnum)

  call cwp_mesh_interf_from_cellvtx_set(code_names(1), &
                                        coupling_name, &
                                        0,             &
                                        n_cell,        &
                                        cell_vtx_idx,  &
                                        cell_vtx,      &
                                        cell_gnum)

  call cwp_mesh_interf_finalize(code_names(1), &
                                coupling_name)

  ! Define fields
  field_name = "field"

  allocate(send_val(n_vtx), &
           recv_val(n_vtx))

  send_val = vtx_coord(1,:) + vtx_coord(2,:) + vtx_coord(3,:)

  call cwp_field_create(code_names(1),                &
                        coupling_name,                &
                        field_name,                   &
                        CWP_DOUBLE,                   &
                        CWP_FIELD_STORAGE_INTERLACED, &
                        1,                            &
                        CWP_DOF_LOCATION_NODE,        &
                        CWP_FIELD_EXCH_SENDRECV,      &
                        CWP_STATUS_ON)

  call cwp_field_data_set(code_names(1),        &
                          coupling_name,        &
                          field_name,           &
                          0,                    &
                          CWP_FIELD_MAP_SOURCE, &
                          send_val)

  call cwp_field_data_set(code_names(1),        &
                          coupling_name,        &
                          field_name,           &
                          0,                    &
                          CWP_FIELD_MAP_TARGET, &
                          recv_val)

  call cwp_time_step_beg(code_names(1), &
                         0.d0)

  ! Compute spatial interpolation weights
  call cwp_spatial_interp_weights_compute(code_names(1), &
                                          coupling_name)

  ! Exchange interpolated fields
  call cwp_field_issend(code_names(1), &
                        coupling_name, &
                        field_name)

  call cwp_field_irecv(code_names(1), &
                       coupling_name, &
                       field_name)

  call cwp_field_wait_issend(code_names(1), &
                             coupling_name, &
                             field_name)

  call cwp_field_wait_irecv(code_names(1), &
                             coupling_name, &
                             field_name)

  call cwp_time_step_end(code_names(1))

  ! Check interpolated field
  max_err = maxval(abs(send_val - recv_val))

  if (max_err > 1.d-14) then
    print *, "Max. interpolation error is too high:", max_err
    stop
  endif

  deallocate(send_val, &
             recv_val)

  deallocate(vtx_coord,    &
             cell_vtx_idx, &
             cell_vtx)

  deallocate(code_names,         &
             coupled_code_names, &
             intra_comms)


  ! Finalize CWIPI
  call cwp_finalize()

  if (i_rank == 0) then
    print *, "The End"
  endif

  ! Finalize MPI
  call mpi_finalize(err)

  end program fortran_new_api_vol_from_cellvtx
