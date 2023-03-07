#include "cwipi_configf.h"

program fortran_new_api_deformable_sol

#ifdef CWP_HAVE_FORTRAN_MPI_MODULE
    use mpi
#endif
    use cwp

    implicit none

#ifndef CWP_HAVE_FORTRAN_MPI_MODULE
    include "mpif.h"
#endif

  !--------------------------------------------------------------------
  integer                                     :: ierr
  integer                                     :: i_rank, n_rank

  integer                                     :: n_code
  character(len = 5),                 pointer :: code_names(:)         => null()
  integer,                            pointer :: is_coupled_rank(:)    => null()
  double precision,                   pointer :: time_init(:)          => null()
  integer,                            pointer :: intra_comms(:)        => null()

  integer                                     :: n_part
  character(len = 5),                 pointer :: coupled_code_names(:) => null()
  character(len = 99)                         :: coupling_name

  double precision, dimension(:), pointer     :: coords
  integer(c_long), pointer, dimension(:)      :: vtx_g_num => null()

  integer, pointer, dimension(:)              :: connec_idx
  integer, pointer, dimension(:)              :: connec
  integer(c_long), pointer, dimension(:)      :: elt_g_num => null()
  integer(c_int)                              :: id_block

  character(len = 99)                         :: field_name
  integer(c_int)                              :: n_components

  double precision                            :: ampl
  double precision                            :: dt
  double precision                            :: freq
  double precision                            :: omega
  double precision                            :: phi
  double precision                            :: randLevel
  double precision                            :: time
  double precision                            :: xmax, xmin, ymax, ymin
  integer                                     :: i, j, it, itdeb, itend, n2, n_partition
  integer                                     :: n_vtx_seg, n_vtx, n_elt

  double precision,                   pointer :: send_field_data(:) => null()
  double precision,                   pointer :: recv_field_data(:) => null()

  integer(c_int)                              :: n_uncomputed_tgts
  integer(c_int),                     pointer :: uncomputed_tgts(:) => null()
  !--------------------------------------------------------------------

  ! MPI Initialization :
  call MPI_Init(ierr)
  call MPI_Comm_rank(mpi_comm_world, i_rank, ierr)
  call MPI_Comm_size(mpi_comm_world, n_rank, ierr)

  ! Check running on correct number of MPI ranks :
  n_partition = 1
  do while ((2 * n_partition**2) < n_rank)
    n_partition = n_partition + 1
  enddo

  n2 = 2 * n_partition**2

  if (n2 /= n_rank) then
    if (i_rank == 0) then
      write(6,*) '      Not executed : only available if the number of processus in the form of 2 * n_partition**2'
    endif
    call exit()
  endif

  ! Initialize CWIPI :
  n_code = 1

  allocate(code_names(n_code),         &
           is_coupled_rank(n_code),    &
           time_init(n_code),          &
           intra_comms(n_code))

  code_names(1)         = "code1"
  is_coupled_rank(1)    = CWP_STATUS_ON
  time_init(1)          = 0.d0

  call CWP_Init(mpi_comm_world,  &
                n_code,          &
                code_names,      &
                is_coupled_rank, &
                time_init,       &
                intra_comms)

  ! Create the coupling :
  ! CWP_DYNAMIC_MESH_DEFORMABLE allows us to take into account the modifications
  ! to the mesh over the coupling steps.
  coupling_name     = "code1_code2";

  allocate(coupled_code_names(n_code))

  coupled_code_names(1) = "code2"

  call CWP_Cpl_create(code_names(1),                                         &
                      coupling_name,                                         &
                      coupled_code_names(1),                                 &
                      CWP_INTERFACE_SURFACE,                                 &
                      CWP_COMM_PAR_WITH_PART,                                &
                      CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, &
                      n_part,                                                &
                      CWP_DYNAMIC_MESH_DEFORMABLE,                           &
                      CWP_TIME_EXCH_USER_CONTROLLED)

  ! Set coupling visualisation:
  call CWP_Visu_set(code_names(1),           &
                    coupling_name,           &
                    1,                       &
                    CWP_VISU_FORMAT_ENSIGHT, &
                    "text")

  ! Set mesh data :
  xmin       = -10.d0
  xmax       =  10.d0
  ymin       = -10.d0
  ymax       =  10.d0
  itdeb      = 1
  itend      = 50
  freq       = 0.20d0
  ampl       = 0.012d0
  phi        = 0.1d0
  n_vtx_seg  = 10
  randLevel  = 0.1d0

  n_vtx = n_vtx_seg * n_vtx_seg
  n_elt = (n_vtx_seg - 1) * (n_vtx_seg - 1)

  allocate(coords(3 * n_vtx))
  allocate(connec_idx(n_elt + 1))
  allocate(connec(4 * n_elt))

  call grid_mesh_f(xmin, &
                   xmax, &
                   ymin, &
                   ymax, &
                   randLevel, &
                   n_vtx_seg, &
                   n_partition, &
                   coords,  &
                   connec_idx,&
                   connec,&
                   intra_comms(1))

  ! Interations :
  time  = 0.0d0
  dt    = 0.1d0
  omega = 2.0d0*acos(-1.0d0)*freq

  do it = itdeb, itend

    time = (it-itdeb)*dt

    do i = 1, n_vtx
      coords((i-1)*3+3) = (coords((i-1)*3+1)*coords((i-1)*3+1)+ &
                          coords((i-1)*3+2)*coords((i-1)*3+2))* &
                          ampl*cos(omega*time+phi)
    enddo

    ! do  i = 1, n_elt
    !   values(i) = 0.
    !   do j = connec_idx(i)+1, connec_idx(i+1)
    !     values(i) = values(i) + coords((connec(j)-1)*3+3)
    !   enddo
    !   values(i) = values(i)/ (connec_idx(i+1) - connec_idx(i))
    ! enddo

    if (it == itdeb) then

      ! Set the mesh vertices coordinates :
      call CWP_Mesh_interf_vtx_set(code_names(1), &
                                   coupling_name, &
                                   0,             &
                                   n_vtx,         &
                                   coords,        &
                                   vtx_g_num)

      ! Set the mesh polygons connectiviy :
      id_block = CWP_Mesh_interf_block_add(code_names(1),       &
                                           coupling_name,       &
                                           CWP_BLOCK_FACE_POLY)

      call CWP_Mesh_interf_f_poly_block_set(code_names(1), &
                                            coupling_name, &
                                            0,             &
                                            id_block,      &
                                            n_elt,         &
                                            connec_idx,    &
                                            connec,        &
                                            elt_g_num)

      ! Finalize mesh :
      call CWP_Mesh_interf_finalize(code_names(1), &
                                    coupling_name)

      ! Create field :
      call CWP_Field_create(code_names(1),                &
                            coupling_name,                &
                            field_name,                   &
                            CWP_DOUBLE,                   &
                            CWP_FIELD_STORAGE_INTERLACED, &
                            n_components,                 &
                            CWP_DOF_LOCATION_NODE,        &
                            CWP_FIELD_EXCH_SEND,          &
                            CWP_STATUS_ON)

      call CWP_Field_create(code_names(1),                &
                            coupling_name,                &
                            field_name,                   &
                            CWP_DOUBLE,                   &
                            CWP_FIELD_STORAGE_INTERLACED, &
                            n_components,                 &
                            CWP_DOF_LOCATION_NODE,        &
                            CWP_FIELD_EXCH_RECV,          &
                            CWP_STATUS_ON)

      ! Set interpolation property :
      call CWP_Spatial_interp_property_set(code_names(1), &
                                           coupling_name, &
                                           "tolerance",   &
                                           "double",      &
                                           "0.1")

    else

      ! Update mesh :
      ! Nothing to do since CWIPI stores the pointers of the arrays passed. Thus if the data
      ! in the pointer is changed, it is automatically in the CWIPI code.

    endif

    ! Compute interpolation weights :
    call CWP_Spatial_interp_weights_compute(code_names(1), &
                                            coupling_name)

    ! Set and exchange the field values :
    call CWP_Field_data_set(code_names(1),        &
                            coupling_name,        &
                            field_name,           &
                            0,                    &
                            CWP_FIELD_MAP_SOURCE, &
                            send_field_data)

    call CWP_Field_issend(code_names(1), &
                          coupling_name, &
                          field_name)

    call CWP_Field_data_set(code_names(1),        &
                            coupling_name,        &
                            field_name,           &
                            0,                    &
                            CWP_FIELD_MAP_TARGET, &
                            recv_field_data)

    call CWP_Field_irecv(code_names(1), &
                         coupling_name, &
                         field_name)

    call CWP_Field_wait_issend(code_names(1), &
                               coupling_name, &
                               field_name)

    call CWP_Field_wait_irecv(code_names(1), &
                              coupling_name, &
                              field_name)

    ! Check interpolation :
    n_uncomputed_tgts = CWP_N_uncomputed_tgts_get(code_names(1), &
                                                  coupling_name, &
                                                  field_name,    &
                                                  0)

    uncomputed_tgts => CWP_Uncomputed_tgts_get(code_names(1), &
                                               coupling_name, &
                                               field_name,    &
                                               0)

  enddo

  ! Delete field :
  call CWP_Field_Del(code_names(1),   &
                     coupling_name,   &
                     field_name)

  ! Delete Mesh :
  call CWP_Mesh_interf_del(code_names(1), &
                           coupling_name)

  ! Delete the coupling :
  call CWP_Cpl_Del(code_names(1), &
                   coupling_name)

  ! Finalize CWIPI :
  call CWP_Finalize()

  ! Finalize MPI :
  call MPI_Finalize(ierr)

end program fortran_new_api_deformable_sol
