#include "cwipi_configf.h"

program fortran_new_api_callback_sendrecv_sol

#ifdef CWP_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif

  use cwp
  use pdm_generate_mesh

  implicit none

#ifndef CWP_HAVE_FORTRAN_MPI_MODULE
    include "mpif.h"
#endif


  !--------------------------------------------------------------------
  integer                                     :: ierr
  integer                                     :: i_rank, n_rank
  integer(c_long), parameter                  :: n_vtx_seg = 10
  integer(c_int),  parameter                  :: n_part    = 1

  ! Coupling
  integer                                     :: n_code
  character(len = 99),                pointer :: code_names(:)         => null()
  integer,                            pointer :: is_active_rank(:)     => null()
  double precision,                   pointer :: time_init(:)          => null()
  integer,                            pointer :: intra_comms(:)        => null()
  character(len = 99),                pointer :: coupled_code_names(:) => null()
  character(len = 99)                         :: coupling_name

  ! Mesh
  integer(c_int)                              :: n_vtx = 0
  double precision,                   pointer :: vtx_coord(:)   => null()
  integer(c_long),                    pointer :: vtx_g_num(:)   => null()
  integer(c_int)                              :: n_elt = 0
  integer(c_int),                     pointer :: elt_vtx_idx(:) => null()
  integer(c_int),                     pointer :: elt_vtx(:)     => null()
  integer(c_long),                    pointer :: elt_g_num(:)   => null()
  integer(c_int)                              :: id_block

  ! Fields
  character(len = 99)                         :: field_name
  double precision,                   pointer :: send_field_data(:) => null()
  double precision,                   pointer :: recv_field_data(:) => null()

  integer                                     :: i, j
  !--------------------------------------------------------------------

  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(mpi_comm_world, i_rank, ierr)
  call MPI_Comm_size(mpi_comm_world, n_rank, ierr)

  ! Initialize CWIPI
  n_code = 1

  allocate(code_names     (n_code), &
           is_active_rank(n_code),  &
           time_init      (n_code), &
           intra_comms    (n_code))

  code_names(1)     = "codeFortran"
  is_active_rank(1) = CWP_STATUS_ON
  time_init(1)      = 0.d0

  call CWP_Init(mpi_comm_world, &
                n_code,         &
                code_names,     &
                is_active_rank, &
                time_init,      &
                intra_comms)

  ! Create the coupling
  coupling_name = "coupling_Python_Fortran"
  allocate(coupled_code_names(n_code))
  coupled_code_names(1) = "codePython"

  call CWP_Cpl_create(code_names(1),                                         &
                      coupling_name,                                         &
                      coupled_code_names(1),                                 &
                      CWP_INTERFACE_SURFACE,                                 &
                      CWP_COMM_PAR_WITH_PART,                                &
                      CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, &
                      n_part,                                                &
                      CWP_DYNAMIC_MESH_STATIC,                               &
                      CWP_TIME_EXCH_USER_CONTROLLED)

  ! Set coupling visualisation
  call CWP_Visu_set(code_names(1),           &
                    coupling_name,           &
                    1,                       &
                    CWP_VISU_FORMAT_ENSIGHT, &
                    "text")

  ! Generate mesh
  call PDM_generate_mesh_rectangle_simplified(intra_comms(1), &
                                              n_vtx_seg,      &
                                              n_vtx,          &
                                              n_elt,          &
                                              vtx_coord,      &
                                              elt_vtx_idx,    &
                                              elt_vtx)

  ! Set mesh vertices
  call CWP_Mesh_interf_vtx_set(code_names(1), &
                               coupling_name, &
                               0,             &
                               n_vtx,         &
                               vtx_coord,     &
                               vtx_g_num)

  ! Set mesh elements
  id_block = CWP_Mesh_interf_block_add(code_names(1),       &
                                       coupling_name,       &
                                       CWP_BLOCK_FACE_POLY)

  call CWP_Mesh_interf_f_poly_block_set(code_names(1), &
                                        coupling_name, &
                                        0,             &
                                        id_block,      &
                                        n_elt,         &
                                        elt_vtx_idx,   &
                                        elt_vtx,       &
                                        elt_g_num)

  ! Finalize mesh
  call CWP_Mesh_interf_finalize(code_names(1), &
                                coupling_name)


  ! Define field
  field_name = "coord_x"

  call CWP_Field_create(code_names(1),                &
                        coupling_name,                &
                        field_name,                   &
                        CWP_DOUBLE,                   &
                        CWP_FIELD_STORAGE_INTERLACED, &
                        1,                            &
                        CWP_DOF_LOCATION_CELL_CENTER, &
                        CWP_FIELD_EXCH_SENDRECV,      &
                        CWP_STATUS_ON)

  allocate(send_field_data(n_elt), &
           recv_field_data(n_elt))

  do i = 1, n_elt
    send_field_data(i) = 0.d0
    do j = elt_vtx_idx(i)+1, elt_vtx_idx(i+1)
      send_field_data(i) = send_field_data(i) + vtx_coord(3*(elt_vtx(j)-1))
    enddo
    send_field_data(i) = send_field_data(i) / (elt_vtx_idx(i+1) - elt_vtx_idx(i))
  enddo


  call CWP_Field_data_set(code_names(1),        &
                          coupling_name,        &
                          field_name,           &
                          0,                    &
                          CWP_FIELD_MAP_SOURCE, &
                          send_field_data)

  call CWP_Field_data_set(code_names(1),        &
                          coupling_name,        &
                          field_name,           &
                          0,                    &
                          CWP_FIELD_MAP_TARGET, &
                          recv_field_data)

  ! Set user-defined interpolation function
  call CWP_Interp_function_set(code_names(1), &
                               coupling_name, &
                               field_name,    &
                               my_interpolation)


  ! Spatial interpolation
  call CWP_Spatial_interp_property_set(code_names(1), &
                                       coupling_name, &
                                       "tolerance",   &
                                       "double",      &
                                       "0.1")

  call CWP_Spatial_interp_weights_compute(code_names(1), &
                                          coupling_name)


  ! Exchange interpolated fields
  call CWP_Field_irecv(code_names(1), &
                       coupling_name, &
                       field_name)

  call CWP_Field_issend(code_names(1), &
                        coupling_name, &
                        field_name)

  call CWP_Field_wait_irecv(code_names(1), &
                            coupling_name, &
                            field_name)

  call CWP_Field_wait_issend(code_names(1), &
                             coupling_name, &
                             field_name)

  ! Delete field
  call CWP_Field_Del(code_names(1), &
                     coupling_name, &
                     field_name)

  ! Delete Mesh
  call CWP_Mesh_interf_del(code_names(1), &
                           coupling_name)

  ! Delete the coupling
  call CWP_Cpl_Del(code_names(1), &
                   coupling_name)


  ! Free memory
  deallocate(code_names,         &
             is_active_rank,     &
             time_init,          &
             intra_comms,        &
             coupled_code_names, &
             send_field_data,    &
             recv_field_data)

  call pdm_fortran_free_c(c_loc(vtx_coord))
  call pdm_fortran_free_c(c_loc(elt_vtx_idx))
  call pdm_fortran_free_c(c_loc(elt_vtx))

  ! Finalize CWIPI
  call CWP_Finalize()

  ! Finalize MPI
  call MPI_Finalize(ierr)


contains

  subroutine my_interpolation(local_code_name,          &
                              cpl_id,                   &
                              field_id,                 &
                              i_part,                   &
                              spatial_interp_algorithm, &
                              storage,                  &
                              buffer_in,                &
                              buffer_out)
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = 1)   :: local_code_name, cpl_id, field_id
    integer(kind = c_int)               :: i_part
    integer(kind = c_int)               :: spatial_interp_algorithm
    integer(kind = c_int)               :: storage
    real(kind = c_double), dimension(*) :: buffer_in
    real(kind = c_double), dimension(*) :: buffer_out

    integer(c_int)                      :: n_components
    integer(c_int)                      :: n_elt_src
    integer(c_int),             pointer :: src_to_tgt_idx(:) => null()
    integer                             :: i, j, k


    n_components = CWP_Interp_field_n_components_get(local_code_name, &
                                                     cpl_id,          &
                                                     field_id)

    call CWP_Interp_src_data_get(local_code_name, &
                                 cpl_id,          &
                                 field_id,        &
                                 i_part,          &
                                 n_elt_src,       &
                                 src_to_tgt_idx)

    do i = 1, n_elt_src
      do j = src_to_tgt_idx(i)+1, src_to_tgt_idx(i+1)
        do k = 1, n_components
          buffer_out(n_components*(j-1)+k) = buffer_in(n_components*(i-1)+k)
          enddo
      enddo
    enddo

  end subroutine my_interpolation

end program fortran_new_api_callback_sendrecv_sol
