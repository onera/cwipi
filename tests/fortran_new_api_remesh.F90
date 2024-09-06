#include "cwipi_configf.h"

program fortran_new_api_remesh

#ifdef CWP_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use cwp

  implicit none

#ifndef CWP_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif


  !--------------------------------------------------------------------
  ! MPI
  integer, parameter          :: comm = MPI_COMM_WORLD
  integer                     :: i_rank                          ! Rank in COMM_WORLD
  integer                     :: err                             ! MPI error code
  integer                     :: i_rank_intra                    ! Rank in intra communicator

  ! CWIPI
  integer                     :: i_code                          ! Code identifier
  character(len = 5), pointer :: code_names(:)         => null() ! Local code name
  character(len = 5), pointer :: coupled_code_names(:) => null() ! Coupled code name
  integer,            pointer :: intra_comms(:)        => null() ! Intra communicator


  ! Coupling
  character(len = 99)         :: cpl_name                        ! Coupling name

  ! Interface mesh
  integer(c_long)             :: n_vtx_seg                       ! Number of vertices along each side of the interface (rectangle)
  integer                     :: n_vtx                           ! Local number of mesh vertices
  integer                     :: n_elt                           ! Local number of mesh elements
  double precision,   pointer :: coords(:,:)    => null()        ! Vertex coordinates
  integer(c_long),    pointer :: vtx_gnum(:)    => null()        ! Vertex global IDs
  integer(c_int),     pointer :: elt_vtx_idx(:) => null()        ! Index for Element->Vertex connectivity
  integer(c_int),     pointer :: elt_vtx(:)     => null()        ! Element->Vertex connectivity
  integer(c_long),    pointer :: elt_gnum(:)    => null()        ! Element global IDs
  integer                     :: id_block                        ! Interface mesh block identifier

  ! Fields
  character(len = 99)         :: field_name1                     ! Field1 name
  character(len = 99)         :: field_name2                     ! Field2 name
  double precision,   pointer :: send_data1(:)                   ! Source field1 pointer
  double precision,   pointer :: recv_data1(:)                   ! Target field1 pointer
  double precision,   pointer :: send_data2(:)                   ! Source field2 pointer
  double precision,   pointer :: recv_data2(:)                   ! Target field2 pointer

  ! Time loop
  integer, parameter          :: n_time_steps = 30               ! Number of time steps
  double precision            :: time                            ! Time
  integer                     :: freq_remesh                     ! Remeshing period
  integer                     :: i_step                          ! Current time step
  !--------------------------------------------------------------------


  ! Initialize MPI
  call mpi_init(err)
  call mpi_comm_rank(comm, i_rank, err)


  ! Initialize CWIPI
  allocate(code_names(1),         &
           coupled_code_names(1), &
           intra_comms(1))

  i_code = mod(i_rank, 2)

  if (i_code .eq. 0) then
    code_names        (1) = "code1"
    coupled_code_names(1) = "code2"
  else
    code_names        (1) = "code2"
    coupled_code_names(1) = "code1"
  endif

  call cwp_init(comm,          &
                1,             &
                code_names,    &
                CWP_STATUS_ON, &
                intra_comms)

  call mpi_comm_rank(intra_comms(1), i_rank_intra, err)


  ! Create coupling environment
  cpl_name = "fortran_new_api_remesh"

  call cwp_cpl_create(code_names(1),                                         &
                      cpl_name,                                              &
                      coupled_code_names(1),                                 &
                      CWP_INTERFACE_SURFACE,                                 &
                      CWP_COMM_PAR_WITH_PART,                                &
                      CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, &
                      1,                                                     &
                      CWP_DYNAMIC_MESH_VARIABLE,                             &
                      CWP_TIME_EXCH_USER_CONTROLLED)

  call cwp_visu_set(code_names(1),           &
                    cpl_name,                &
                    1,                       &
                    CWP_VISU_FORMAT_ENSIGHT, &
                    "text")


  ! Initial interface mesh
  n_vtx_seg = 10 - 3*i_code
  call cwpt_generate_mesh_rectangle_simplified(intra_comms(1), &
                                               n_vtx_seg,      &
                                               n_vtx,          &
                                               n_elt,          &
                                               coords,         &
                                               elt_vtx_idx,    &
                                               elt_vtx)

  vtx_gnum => null()
  call cwp_mesh_interf_vtx_set(code_names(1), &
                               cpl_name,      &
                               0,             &
                               n_vtx,         &
                               coords,        &
                               vtx_gnum)

  id_block = cwp_mesh_interf_block_add(code_names(1), &
                                       cpl_name,      &
                                       CWP_BLOCK_FACE_TRIA3)

  elt_gnum => null()
  call cwp_mesh_interf_block_std_set(code_names(1), &
                                     cpl_name,      &
                                     0,             &
                                     id_block,      &
                                     n_elt,         &
                                     elt_vtx,       &
                                     elt_gnum)

  call cwp_mesh_interf_finalize(code_names(1), &
                                cpl_name)


  ! Create fields
  field_name1 = "field1"
  field_name2 = "field2"

  if (i_code .eq. 0) then
    ! Code0 is node-centered
    allocate(send_data1(n_vtx), &
             recv_data2(n_vtx))

    call cwp_field_create(code_names(1),                &
                          cpl_name,                     &
                          field_name1,                  &
                          CWP_DOUBLE,                   &
                          CWP_FIELD_STORAGE_INTERLACED, &
                          1,                            &
                          CWP_DOF_LOCATION_NODE,        &
                          CWP_FIELD_EXCH_SEND,          &
                          CWP_STATUS_ON)

    call cwp_field_create(code_names(1),                &
                          cpl_name,                     &
                          field_name2,                  &
                          CWP_DOUBLE,                   &
                          CWP_FIELD_STORAGE_INTERLACED, &
                          1,                            &
                          CWP_DOF_LOCATION_NODE,        &
                          CWP_FIELD_EXCH_RECV,          &
                          CWP_STATUS_ON)
  else
    ! Code1 is cell-centered
    allocate(recv_data1(n_elt), &
             send_data2(n_elt))

    call cwp_field_create(code_names(1),                &
                          cpl_name,                     &
                          field_name1,                  &
                          CWP_DOUBLE,                   &
                          CWP_FIELD_STORAGE_INTERLACED, &
                          1,                            &
                          CWP_DOF_LOCATION_CELL_CENTER, &
                          CWP_FIELD_EXCH_RECV,          &
                          CWP_STATUS_ON)

    call cwp_field_create(code_names(1),                &
                          cpl_name,                     &
                          field_name2,                  &
                          CWP_DOUBLE,                   &
                          CWP_FIELD_STORAGE_INTERLACED, &
                          1,                            &
                          CWP_DOF_LOCATION_CELL_CENTER, &
                          CWP_FIELD_EXCH_SEND,          &
                          CWP_STATUS_ON)
  endif


  ! Time loop
  time = 0.d0
  freq_remesh = 5 + 3*i_code

  do i_step = 1, n_time_steps

    if (i_rank .eq. 0) then
      print "(a5,i0)", "Step ", i_step
    endif

    ! Begin time step
    call cwp_time_step_beg(code_names(1), &
                           time)
    time = time + 1.d0


    ! Remesh if required
    if (mod(i_step, freq_remesh) .eq. 0) then

      if (i_rank_intra .eq. 0) then
        print *, "  Remesh ", code_names(1)
      endif

      call CWP_Mesh_interf_del(code_names(1), cpl_name)

      call cwpt_fortran_free_c(c_loc(coords))
      call cwpt_fortran_free_c(c_loc(elt_vtx_idx))
      call cwpt_fortran_free_c(c_loc(elt_vtx))

      n_vtx_seg = n_vtx_seg + 2 + i_code


      call cwpt_generate_mesh_rectangle_simplified(intra_comms(1), &
                                                   n_vtx_seg,      &
                                                   n_vtx,          &
                                                   n_elt,          &
                                                   coords,         &
                                                   elt_vtx_idx,    &
                                                   elt_vtx)

      call cwp_mesh_interf_vtx_set(code_names(1), &
                                   cpl_name,      &
                                   0,             &
                                   n_vtx,         &
                                   coords,        &
                                   vtx_gnum)

      id_block = cwp_mesh_interf_block_add(code_names(1), &
                                           cpl_name,      &
                                           CWP_BLOCK_FACE_TRIA3)

      call cwp_mesh_interf_block_std_set(code_names(1), &
                                         cpl_name,      &
                                         0,             &
                                         id_block,      &
                                         n_elt,         &
                                         elt_vtx,       &
                                         elt_gnum)

      call cwp_mesh_interf_finalize(code_names(1), &
                                    cpl_name)


      ! Field data
      if (i_code .eq. 0) then
        deallocate(send_data1, recv_data2)
        allocate(send_data1(n_vtx), &
                 recv_data2(n_vtx))
      else
        deallocate(recv_data1, send_data2)
        allocate(recv_data1(n_elt), &
                 send_data2(n_elt))
      endif

    endif



    ! Set field value pointers
    if (i_code .eq. 0) then
      call eval_field1(coords, time, send_data1)
      call cwp_field_data_set(code_names(1),        &
                              cpl_name,             &
                              field_name1,          &
                              0,                    &
                              CWP_FIELD_MAP_SOURCE, &
                              send_data1)
      call cwp_field_data_set(code_names(1),        &
                              cpl_name,             &
                              field_name2,          &
                              0,                    &
                              CWP_FIELD_MAP_TARGET, &
                              recv_data2)
    else
      call eval_field2(n_elt, elt_vtx_idx, elt_vtx, coords, time, send_data2)
      call cwp_field_data_set(code_names(1),        &
                              cpl_name,             &
                              field_name1,          &
                              0,                    &
                              CWP_FIELD_MAP_TARGET, &
                              recv_data1)
      call cwp_field_data_set(code_names(1),        &
                              cpl_name,             &
                              field_name2,          &
                              0,                    &
                              CWP_FIELD_MAP_SOURCE, &
                              send_data2)
    endif


    ! Compute spatial interpolation weights
    call cwp_spatial_interp_weights_compute(code_names(1), cpl_name)


    ! Exchange fields
    if (i_code .eq. 0) then
      call cwp_field_issend(code_names(1), cpl_name, field_name1)
      call cwp_field_irecv (code_names(1), cpl_name, field_name2)
    else
      call cwp_field_irecv (code_names(1), cpl_name, field_name1)
      call cwp_field_issend(code_names(1), cpl_name, field_name2)
    endif

    ! Overlap exchanges with computations if possible...

    ! Finalize field exchange
    if (i_code .eq. 0) then
      call cwp_field_wait_issend(code_names(1), cpl_name, field_name1)
      call cwp_field_wait_irecv (code_names(1), cpl_name, field_name2)
    else
      call cwp_field_wait_irecv (code_names(1), cpl_name, field_name1)
      call cwp_field_wait_issend(code_names(1), cpl_name, field_name2)
    endif

    ! End time step
    call cwp_time_step_end(code_names(1))

  enddo


  ! Delete mesh
  call cwp_mesh_interf_del(code_names(1), cpl_name)

  ! Delete coupling
  call cwp_cpl_del(code_names(1), cpl_name)

  ! Finalize
  call cwp_finalize()

  ! Free memory
  deallocate(code_names,         &
             coupled_code_names, &
             intra_comms)

  call cwpt_fortran_free_c(c_loc(coords))
  call cwpt_fortran_free_c(c_loc(elt_vtx_idx))
  call cwpt_fortran_free_c(c_loc(elt_vtx))
  if (i_code .eq. 0) then
    deallocate(send_data1, &
               recv_data2)
  else
    deallocate(recv_data1, &
               send_data2)
  endif

  if (i_rank .eq. 0) then
    print *, "The End :)"
  endif

  call mpi_finalize(err)


contains

subroutine eval_field1(coords, &
                       time,   &
                       field)
  implicit none

  double precision, pointer :: coords(:,:)
  double precision          :: time
  double precision, pointer :: field(:)

  field(:) = cos(sqrt(coords(1,:)**2 + coords(2,:)**2) - 0.1*time)

end subroutine eval_field1


subroutine eval_field2(n_elt,       &
                       elt_vtx_idx, &
                       elt_vtx,     &
                       coords,      &
                       time,        &
                       field)
  implicit none

  integer                   :: n_elt
  integer(c_int),   pointer :: elt_vtx_idx(:)
  integer(c_int),   pointer :: elt_vtx(:)
  double precision, pointer :: coords(:,:)
  double precision          :: time
  double precision, pointer :: field(:)
  integer                   :: i_elt, i

  do i_elt = 1, n_elt
    field(i_elt) = 0.d0
    do i = elt_vtx_idx(i_elt)+1, elt_vtx_idx(i_elt+1)
      field(i_elt) = field(i_elt) + cos(coords(1, elt_vtx(i)) + 0.1*time)
    enddo
    field(i_elt) = field(i_elt) / (elt_vtx_idx(i_elt+1) - elt_vtx_idx(i_elt))
  enddo

end subroutine eval_field2

end program fortran_new_api_remesh
