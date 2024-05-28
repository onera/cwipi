#include "cwipi_configf.h"

program fortran_new_api_empty_parts

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
  integer,            pointer  :: intra_comms(:)       => null()
  character(len = 99)         :: coupling_name

  integer                     :: n_vtx_seg
  integer                     :: n_vtx
  double precision,   pointer :: vtx_coord(:,:)
  integer(c_long),    pointer :: vtx_ln_to_gn(:)
  integer                     :: n_elt
  integer(c_int),     pointer :: connec(:)
  integer(c_long),    pointer :: elt_ln_to_gn(:)
  integer                     :: id_block

  character(len = 99)         :: field_name
  double precision,   pointer :: send_data(:)
  double precision,   pointer :: recv_data(:)

  integer                     :: i
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
  coupling_name = "fortran_new_api_empty_parts"

  call cwp_cpl_create(code_names(1),                                         &
                      coupling_name,                                         &
                      coupled_code_names(1),                                 &
                      CWP_INTERFACE_SURFACE,                                 &
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
  n_vtx_seg = 3
  if (is_code1) then
    n_vtx_seg = 4
  endif
  call generate_mesh(intra_comms(1), &
                     n_vtx_seg,      &
                     n_vtx,          &
                     vtx_coord,      &
                     vtx_ln_to_gn,   &
                     n_elt,          &
                     connec,         &
                     elt_ln_to_gn)

  print *, "rank", i_rank, " has", n_vtx, " vtx and", n_elt, " elt"

  call cwp_mesh_interf_vtx_set(code_names(1), &
                               coupling_name, &
                               0,             &
                               n_vtx,         &
                               vtx_coord,     &
                               vtx_ln_to_gn)

  id_block = cwp_mesh_interf_block_add(code_names(1),       &
                                       coupling_name,       &
                                       CWP_BLOCK_FACE_QUAD4)

  call cwp_mesh_interf_block_std_set(code_names(1), &
                                     coupling_name, &
                                     0,             &
                                     id_block,      &
                                     n_elt,         &
                                     connec,        &
                                     elt_ln_to_gn)

  call cwp_mesh_interf_finalize(code_names(1), &
                                coupling_name)

  ! Define fields
  field_name = "field"

  allocate(send_data(n_elt), &
           recv_data(n_vtx))

  if (is_code1) then
    do i = 1, n_elt
      send_data(i) = real(elt_ln_to_gn(i), kind=8)
    enddo

    call cwp_field_create(code_names(1),                &
                          coupling_name,                &
                          field_name,                   &
                          CWP_DOUBLE,                   &
                          CWP_FIELD_STORAGE_INTERLACED, &
                          1,                            &
                          CWP_DOF_LOCATION_CELL_CENTER, &
                          CWP_FIELD_EXCH_SEND,          &
                          CWP_STATUS_ON)

    call cwp_field_data_set(code_names(1),        &
                            coupling_name,        &
                            field_name,           &
                            0,                    &
                            CWP_FIELD_MAP_SOURCE, &
                            send_data)
  else
    call cwp_field_create(code_names(1),                &
                          coupling_name,                &
                          field_name,                   &
                          CWP_DOUBLE,                   &
                          CWP_FIELD_STORAGE_INTERLACED, &
                          1,                            &
                          CWP_DOF_LOCATION_NODE,        &
                          CWP_FIELD_EXCH_RECV,          &
                          CWP_STATUS_ON)

    call cwp_field_data_set(code_names(1),        &
                            coupling_name,        &
                            field_name,           &
                            0,                    &
                            CWP_FIELD_MAP_TARGET, &
                            recv_data)
  endif


  call cwp_time_step_beg(code_names(1), &
                         0.d0)

  ! Compute spatial interpolation weights
  call cwp_spatial_interp_weights_compute(code_names(1), &
                                          coupling_name)

  ! Exchange interpolated fields
  if (is_code1) then
    call cwp_field_issend(code_names(1), &
                          coupling_name, &
                          field_name)
  else
    call cwp_field_irecv(code_names(1), &
                         coupling_name, &
                         field_name)
  endif

  if (is_code1) then
    call cwp_field_wait_issend(code_names(1), &
                               coupling_name, &
                               field_name)
  else
    call cwp_field_wait_irecv(code_names(1), &
                              coupling_name, &
                              field_name)
  endif

  call cwp_time_step_end(code_names(1))


  if (i_rank == 0) then
    print *, "End :)"
  endif


  ! Free memory
  deallocate(code_names,         &
             coupled_code_names, &
             intra_comms,        &
             send_data,          &
             recv_data)

  if (associated(vtx_coord)) then
    deallocate(vtx_coord)
  endif
  if (associated(vtx_ln_to_gn)) then
    deallocate(vtx_ln_to_gn)
  endif
  if (associated(connec)) then
    deallocate(connec)
  endif
  if (associated(elt_ln_to_gn)) then
    deallocate(elt_ln_to_gn)
  endif


  ! Finalize CWIPI
  call cwp_finalize()

  ! Finalize MPI
  call mpi_finalize(err)

contains

subroutine generate_mesh(intra_comm,   &
                         n_vtx_seg,    &
                         n_vtx,        &
                         vtx_coord,    &
                         vtx_ln_to_gn, &
                         n_elt,        &
                         connec,       &
                         elt_ln_to_gn)
  implicit none

  integer, intent(in)       :: intra_comm
  integer, intent(in)       :: n_vtx_seg
  integer, intent(out)      :: n_vtx
  double precision, pointer :: vtx_coord(:,:)
  integer(c_long),  pointer :: vtx_ln_to_gn(:)
  integer, intent(out)      :: n_elt
  integer(c_int),   pointer :: connec(:)
  integer(c_long),  pointer :: elt_ln_to_gn(:)

  double precision          :: step
  integer                   :: i_rank_intra, err
  integer                   :: i, j, k

  call mpi_comm_rank(intra_comm, i_rank_intra, err)

  if (i_rank_intra == 0) then
    n_vtx = n_vtx_seg * n_vtx_seg
    n_elt = (n_vtx_seg - 1) * (n_vtx_seg - 1)

    allocate(vtx_coord(3, n_vtx), &
             vtx_ln_to_gn(n_vtx), &
             connec(4*n_elt),     &
             elt_ln_to_gn(n_elt))

    step = 1.d0 / real(n_vtx_seg - 1, kind=8)

    ! Vertices
    k = 1
    do j = 1, n_vtx_seg
      do i = 1, n_vtx_seg
        vtx_coord(1,k) = (i - 1) * step
        vtx_coord(2,k) = (j - 1) * step
        vtx_coord(3,k) = 0
        vtx_ln_to_gn(k) = k
        k = k + 1
      enddo
    enddo

    ! Elements
    k = 1
    do j = 1, n_vtx_seg-1
      do i = 1, n_vtx_seg-1
        connec(4*(k-1)+1) = n_vtx_seg*(j-1) + (i-1) + 1
        connec(4*(k-1)+2) = n_vtx_seg*(j-1) + i     + 1
        connec(4*(k-1)+3) = n_vtx_seg*j     + i     + 1
        connec(4*(k-1)+4) = n_vtx_seg*j     + (i-1) + 1
        elt_ln_to_gn(k) = k
        k = k + 1
      enddo
    enddo

  else
    n_vtx = 0
    n_elt = 0

    vtx_coord    => null()
    vtx_ln_to_gn => null()
    connec       => null()
    elt_ln_to_gn => null()
  endif

end subroutine generate_mesh

end program fortran_new_api_empty_parts
