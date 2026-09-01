!-----------------------------------------------------------------------------
! This file is part of the CWIPI library.
!
! Copyright (C) 2021-2023  ONERA
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

program testf
#ifdef CWP_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use cwp
  use, intrinsic :: iso_c_binding

  implicit none

#ifndef CWP_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  type my_mesh
    integer(c_int)            :: pn_vtx
    double precision, pointer :: pvtx_coord(:,:) => null()
    double precision, pointer :: field(:)        => null()
    integer(c_int),   pointer :: pface_vtx(:)    => null()
  end type my_mesh

  ! --------------------------------------------------------------------
  integer, parameter                 :: comm = MPI_COMM_WORLD
  integer                            :: ierr
  integer                            :: i_rank
  integer                            :: n_rank

  ! CMD args
  logical                            :: verbose        = .false.
  integer(c_int)                     :: spatial_interp_algo = CWP_SPATIAL_INTERP_FROM_IDENTITY

  type(my_mesh), allocatable         :: mesh(:)
  integer(c_int)                     :: id_block

  integer(c_long),           pointer :: vtx_ln_to_gn(:)  => null()

  integer(c_int)                     :: n_code
  integer,                   pointer :: code_id(:)           => null()
  character(len=5),          pointer :: code_name(:)         => null()
  character(len=5),          pointer :: coupled_code_name(:) => null()
  integer(c_int)                     :: is_active_rank = CWP_STATUS_ON
  integer(c_int),            pointer :: intra_comms(:)       => null()
  character(len=99)                  :: coupling_name


  integer(c_int)                     :: visu_status, map_type, exch_type
  integer(c_int)                     :: stride
  character(len=99)                  :: field_name

  character                          :: strnum
  integer                            :: i, k
  integer                            :: iiunit = 13
  ! --------------------------------------------------------------------

  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(comm, i_rank, ierr)
  call MPI_Comm_size(comm, n_rank, ierr)

  if (verbose) then
    write (strnum, '(i1)') i_rank
    open(unit=iiunit, file='fortran_new_api_point_'//strnum//'.log')
  endif

  n_code = 2
  allocate(code_id(n_code),           &
           code_name(n_code),         &
           coupled_code_name(n_code), &
           intra_comms(n_code),       &
           mesh(n_code))

  do i = 1,n_code
    code_id(i) = i
    write (strnum, '(i1)') i
    code_name        (i) = "code" // strnum
    write (strnum, '(i1)') mod(i,2)+1
    coupled_code_name(i) = "code" // strnum
    if (verbose) then
      write(iiunit, *) "running ", code_name(i), ", coupled with ", coupled_code_name(i)
    endif
  enddo

  call CWP_output_fortran_unit_set(iiunit)

  call CWP_Init(comm,           &
                n_code,         &
                code_name,      &
                is_active_rank, &
                intra_comms)

  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "CWIPI Init OK"
  endif

  coupling_name = "fortran_new_api_point"
  do i = 1, n_code
    call CWP_Cpl_create(code_name(i),                  &
                        coupling_name,                 &
                        coupled_code_name(i),          &
                        CWP_INTERFACE_POINT,           &
                        CWP_COMM_PAR_WITH_PART,        &
                        spatial_interp_algo,           &
                        1,                             &
                        CWP_DYNAMIC_MESH_STATIC,       &
                        CWP_TIME_EXCH_USER_CONTROLLED)
  enddo

  ! this must be performed in 2 separate loops if the intra comms do overlap
  do i = 1, n_code
    call CWP_Visu_set(code_name(i),            &
                      coupling_name,           &
                      1,                       &
                      CWP_VISU_FORMAT_ENSIGHT, &
                      "text")
  enddo

  do i = 1, n_code
    call CWP_Cpl_barrier(code_name(i), &
                         coupling_name)
  enddo

  if (i_rank == 0) then
    print *, "Create coupling OK"
  endif

  ! Define interface mesh
  do i = 1, n_code

    id_block = CWP_Mesh_interf_block_add(code_name(i),   &
                                         coupling_name,  &
                                         CWP_BLOCK_NODE)

    if (i_rank == i-1) then
      mesh(i)%pn_vtx = 1
      allocate(mesh(i)%pvtx_coord(3,1))
      mesh(i)%pvtx_coord = 0.0
    else
      mesh(i)%pn_vtx = 0
      allocate(mesh(i)%pvtx_coord(3,0))
    endif

    call CWP_Mesh_interf_vtx_set(code_name(i),       &
                                 coupling_name,      &
                                 0,                  &
                                 mesh(i)%pn_vtx,     &
                                 mesh(i)%pvtx_coord, &
                                 vtx_ln_to_gn)

    allocate(mesh(i)%pface_vtx(mesh(i)%pn_vtx))
    do k = 1, mesh(i)%pn_vtx
      mesh(i)%pface_vtx(k) = k
    enddo

    call CWP_Mesh_interf_block_std_set(code_name(i),      &
                                       coupling_name,     &
                                       0,                 &
                                       id_block,          &
                                       mesh(i)%pn_vtx,    &
                                       mesh(i)%pface_vtx, &
                                       vtx_ln_to_gn);

    call CWP_Mesh_interf_finalize(code_name(i),  &
                                  coupling_name)

  enddo

  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "Set mesh OK"
  endif

  ! Create fields
  visu_status = CWP_STATUS_ON
  stride = 2

  do i = 1, n_code
    allocate(mesh(i)%field(mesh(i)%pn_vtx * stride))
    if (code_id(i) == 1) then
      exch_type = CWP_FIELD_EXCH_SEND
      map_type  = CWP_FIELD_MAP_SOURCE
      mesh(i)%field(k::stride) = mesh(i)%pvtx_coord(k,:)
    else
      exch_type = CWP_FIELD_EXCH_RECV
      map_type  = CWP_FIELD_MAP_TARGET
    endif

    field_name="exchanged_field"
    call CWP_Field_create(code_name(i),                 &
                          coupling_name,                &
                          field_name,                   &
                          CWP_DOUBLE,                   &
                          CWP_FIELD_STORAGE_INTERLACED, &
                          stride,                       &
                          CWP_DOF_LOCATION_CELL_CENTER, &
                          exch_type,                    &
                          visu_status)

    call CWP_Field_data_set(code_name(i),  &
                            coupling_name, &
                            field_name,    &
                            0,             &
                            map_type,      &
                            mesh(i)%field)

  enddo

  do i = 1, n_code
    call CWP_Time_step_beg(code_name(i), &
                           0.d0)
  enddo

  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "Create fields OK"
  endif

  do i = 1, n_code
    ! Compute spatial interpolation weights
    call CWP_Spatial_interp_weights_compute(code_name(i),  &
                                            coupling_name)
  enddo

  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "Interpolation weights computation OK"
  endif

  ! Exchange interpolated field
  do i = 1, n_code
    if (code_id(i) == 1) then
      call CWP_Field_issend(code_name(i),  &
                            coupling_name, &
                            field_name)
    else
      call CWP_Field_irecv(code_name(i),  &
                           coupling_name, &
                           field_name)
    endif
  enddo

  do i = 1, n_code
    if (code_id(i) == 1) then
      call CWP_Field_wait_issend(code_name(i),  &
                                 coupling_name, &
                                 field_name)
    else
      call CWP_Field_wait_irecv(code_name(i),  &
                                coupling_name, &
                                field_name)
    endif
  enddo

  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "Exchange interpolated field OK"
  endif

  if (verbose) then
    close(iiunit)
  endif

  do i = 1, n_code
    call CWP_Time_step_end(code_name(i))
  enddo

  do i = 1, n_code
    call CWP_Field_del(code_name(i),  &
                       coupling_name, &
                       field_name)
    deallocate(mesh(i)%field)
  enddo

  ! Delete interface mesh
  do i = 1, n_code
    call CWP_Mesh_interf_del(code_name(i),  &
                             coupling_name)
    deallocate(mesh(i)%pvtx_coord)
    deallocate(mesh(i)%pface_vtx)
  enddo

  ! Delete coupling
  do i = 1, n_code
    call CWP_Cpl_Del(code_name(i),  &
                     coupling_name)
  enddo

  ! Free memory
  deallocate(code_id,           &
             code_name,         &
             coupled_code_name, &
             intra_comms,       &
             mesh)

  call MPI_Barrier(MPI_comm_world, ierr)
  if (i_rank == 0) then
    print *, "End"
  endif

  ! Finalize
  call CWP_Finalize()
  call MPI_Finalize(ierr)

end program testf
