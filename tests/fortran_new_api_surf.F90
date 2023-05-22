
#include "cwipi_configf.h"

program testf
#ifdef CWP_HAVE_FORTRAN_MPI_MODULE  
    use mpi
#endif
  use cwp
  use pdm_pointer_array

  implicit none

#ifndef CWP_HAVE_FORTRAN_MPI_MODULE  
    include "mpif.h"
#endif  


  type my_type
    integer(c_long), pointer :: g_num(:) => null()
    integer(c_int),  pointer :: data(:)  => null()
  end type my_type

  !--------------------------------------------------------------------
  integer                       :: ierr
  integer                       :: i_rank, n_rank

  integer,          parameter   :: n_vtx_seg  = 10
  double precision, parameter   :: rand_level = 0.25
  double precision              :: shift

  integer                       :: n_code, n_part
  character(len = 5),   pointer :: code_names(:)         => null()
  character(len = 5),   pointer :: coupled_code_names(:) => null()
  integer,              pointer :: is_coupled_rank(:)    => null()
  double precision,     pointer :: time_init(:)          => null()
  integer,              pointer :: intra_comms(:)        => null()
  character(len = 99)           :: coupling_name

  integer(c_int)                :: n_vtx, n_elt
  double precision,     pointer :: vtx_coord(:,:) => null()
  integer(c_int),       pointer :: connec_idx(:)  => null()
  integer(c_int),       pointer :: connec(:)      => null()
  integer(c_long),      pointer :: vtx_g_num(:)   => null()
  integer(c_long),      pointer :: elt_g_num(:)   => null()
  integer(c_int)                :: id_block

  character(len = 99)           :: field_name
  integer(c_int)                :: map_type, exch_type, visu_status
  double precision,     pointer :: field_data(:) => null()

  integer(c_int)                :: n_computed_tgts
  integer(c_int),       pointer :: computed_tgts(:) => null()

  integer(c_int),       target  :: toto
  integer(c_int)                :: n_param
  type(c_ptr)                   :: param_value = C_NULL_PTR
  integer(c_int),       pointer :: tata => null()

  integer(c_int)                :: n_elt2
  integer(c_int),       pointer :: connec_idx2(:) => null()
  integer(c_int),       pointer :: connec2(:)     => null()
  ! integer(c_long),      pointer :: vtx_g_num2(:)  => null()
  integer(c_long),      pointer :: elt_g_num2(:)  => null()

  integer                       :: i, ivtx, n_wrong
  double precision              :: distance

  character                     :: strnum
  logical                       :: debug = .true.

  integer(c_int)                :: n_elt_src
  integer(c_int), pointer       :: src_to_tgt_idx(:) => null()
  double precision, pointer     :: interp_weights(:) => null()
  integer(c_int)                :: n_elt_tgt
  integer(c_int)                :: n_ref_tgt
  integer(c_int), pointer       :: ref_tgt(:)        => null()
  integer(c_int), pointer       :: tgt_to_src_idx(:) => null()

  !--> list getters
  character(256), allocatable :: code_list(:)
  character(256), allocatable :: loc_code_list(:)
  character(256), allocatable :: f_param_names(:)

  !--> output file
  ! character(len = 23),   pointer :: output_file(:)         => null()
  integer :: iiunit = 9

  !--> reduce
  integer(c_int),       pointer :: res => null()
  character(len = 5),   pointer :: g_code_names(:)         => null()

  ! Global data
  character(len=99)             :: global_data_name
  integer(c_int), pointer       :: global_data(:,:) => null()

  ! Part data
  character(len=99)             :: part_data_name
  integer(c_int), pointer       :: recv_data(:) => null()
  type(PDM_pointer_array_t)     :: gnum_elt, part_data
  integer(c_int), pointer       :: n_elt_part(:) => null()
  type(my_type), allocatable    :: my_part(:)
  integer(c_int)                :: n_comp, j, request
  !--------------------------------------------------------------------


  !! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_comm_world, i_rank, ierr)
  call MPI_Comm_size(MPI_comm_world, n_rank, ierr)

  if (n_rank /= 2) then
    print *, "n_rank must be 2"
    stop
  endif

  if (debug) then
    write (strnum, '(i1)') i_rank
    open(unit=iiunit, file='fortran_new_api_surf_'//strnum//'.txt')
  endif

  !! Initialize CWIPI
  n_code = 1
  n_part = 1

  allocate(code_names(n_code),         &
           coupled_code_names(n_code), &
           is_coupled_rank(n_code),    &
           time_init(n_code),          &
           intra_comms(n_code),        &
           ! output_file(1),             &
           g_code_names(2))

  if (i_rank == 0) then
    code_names(1)         = "code1"
    coupled_code_names(1) = "code2"
    is_coupled_rank(1)    = CWP_STATUS_ON
    time_init(1)          = 0.d0
  else
    code_names(1)         = "code2"
    coupled_code_names(1) = "code1"
    is_coupled_rank(1)    = CWP_STATUS_ON
    time_init(1)          = 0.d0
  endif

  !--> output file
  ! output_file(1) = "fortran_output_file.txt"
  ! call CWP_Output_file_set(output_file(1))
  call cwp_output_fortran_unit_set(iiunit)

  call CWP_Init(MPI_comm_world,  &
                n_code,          &
                code_names,      &
                is_coupled_rank, &
                time_init,       &
                intra_comms)

  print *, "> CWP_set_toto", loc(mon_toto)
  call CWP_set_toto(mon_toto)
  print *, "> CWP_call_toto"
  call CWP_call_toto()
  print *, "OK!"

  !-->>
  if (debug) then
    write(iiunit,*) "n_code =", CWP_Codes_nb_get()
    write(iiunit,*) "n_local_code =", CWP_Loc_codes_nb_get()
    write(iiunit,*) "state =", CWP_State_get(code_names(1))
  endif
  call CWP_Properties_dump()

  !--> character array getters
  code_list     = CWP_Codes_list_get()
  loc_code_list = CWP_Loc_codes_list_get()

  !--> n_total_codes = 2
  if (debug) then
    do i=1,2
      write(iiunit,*) i_rank, " --> ", "code_list(", i, ") :", code_list(i)
    end do

    write(iiunit,*) i_rank, " --> ", "loc_code_list(", 1, ") :", loc_code_list(1)
  endif

  toto = 5
  if (code_names(1) == "code1") then
    ! call CWP_Param_lock("code1")
    call CWP_Param_add("code1", "toto", toto)
    ! call CWP_Param_unlock("code1")
  endif

  toto = 12
  if (code_names(1) == "code2") then
    ! call CWP_Param_lock("code2")
    call CWP_Param_add("code2", "toto", toto)
    ! call CWP_Param_unlock("code2")
  endif

  call MPI_Barrier(MPI_comm_world, ierr)

  n_param = CWP_Param_n_get("code1", &
                            CWP_INT)
  if (debug) then
    write(iiunit,*) "n_param =", n_param
  endif

  ! write(iiunit,*) code_names(1), ", param is? :", CWP_Param_is("code1",  &
  !                                                       "toto", &
  !                                                       CWP_INT)

  if (code_names(1) == "code1") then
    call CWP_Param_get("code1",  &
                       "toto",   &
                       CWP_INT,  &
                       param_value)
    call c_f_pointer(param_value, tata)
    if (debug) then
      write(iiunit,*) "tata =", loc(tata)
    endif
  endif
  !<<--

  ! --> get list of int parameter names
  call CWP_Param_list_get("code1",  &
                          CWP_INT,  &
                          n_param,  &
                          f_param_names)

  if (debug) then
    do i=1, n_param
     write(iiunit,*) "f_param_names : ", i, " -> ", f_param_names(i)
    end do
  endif

  !--> sum of parameters
  g_code_names(1) = "code1"
  g_code_names(2) = "code2"
  call CWP_Param_reduce(CWP_OP_SUM, &
                        "toto",     &
                        CWP_INT,    &
                        param_value,&
                        2,          &
                        g_code_names)

  call c_f_pointer(param_value, res)
  if (debug) then
    write(iiunit,*) "res of sum reduce : ", loc(res)
  endif

  !! Create a coupling
  coupling_name = "fortran_new_api_surf"
  call CWP_Cpl_create(code_names(1),                                         &
                      coupling_name,                                         &
                      coupled_code_names(1),                                 &
                      CWP_INTERFACE_SURFACE,                                 &
                      CWP_COMM_PAR_WITH_PART,                                &
                      CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, &
                      n_part,                                                &
                      CWP_DYNAMIC_MESH_STATIC,                               &
                      CWP_TIME_EXCH_USER_CONTROLLED)

  call CWP_Visu_set(code_names(1),           &
                    coupling_name,           &
                    1,                       &
                    CWP_VISU_FORMAT_ENSIGHT, &
                    "text")

  !! Define the interface mesh
  shift = 0.1*i_rank
  call gen_mesh(shift, shift+1.d0, &
                shift, shift+1.d0, &
                rand_level,        &
                n_vtx_seg,         &
                n_vtx,             &
                vtx_coord,         &
                n_elt,             &
                connec_idx,        &
                connec,            &
                vtx_g_num,         &
                elt_g_num)


  call CWP_Mesh_interf_vtx_set(code_names(1), &
                               coupling_name, &
                               0,             &
                               n_vtx,         &
                               vtx_coord,     &
                               vtx_g_num)

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

  call CWP_Mesh_interf_finalize(code_names(1), &
                                coupling_name)


  !-->>
  call CWP_Mesh_interf_f_poly_block_get(code_names(1), &
                                        coupling_name, &
                                        0,             &
                                        id_block,      &
                                        n_elt2,        &
                                        connec_idx2,   &
                                        connec2,       &
                                        elt_g_num2)
  !<<--


  !! Create fields
  field_name = "field"
  visu_status = CWP_STATUS_ON
  allocate(field_data(3*n_vtx))
  if (code_names(1) == "code1") then
    map_type  = CWP_FIELD_MAP_SOURCE
    exch_type = CWP_FIELD_EXCH_SEND
    field_data(1:3*n_vtx:3) = vtx_coord(1,1:n_vtx)
    field_data(2:3*n_vtx:3) = vtx_coord(2,1:n_vtx)
    field_data(3:3*n_vtx:3) = vtx_coord(3,1:n_vtx)
  else
    map_type  = CWP_FIELD_MAP_TARGET
    exch_type = CWP_FIELD_EXCH_RECV
  endif
  call CWP_Field_create(code_names(1),                &
                        coupling_name,                &
                        field_name,                   &
                        CWP_DOUBLE,                   &
                        CWP_FIELD_STORAGE_INTERLACED, &
                        3,                            &
                        CWP_DOF_LOCATION_NODE,        &
                        exch_type,                    &
                        visu_status)

  call CWP_Time_step_beg(code_names(1), &
                         time_init(1))

  call CWP_Field_data_set(code_names(1), &
                          coupling_name, &
                          field_name,    &
                          0,             &
                          map_type,      &
                          field_data)

  !! Set tolerance for spatial interpolation
  call CWP_Spatial_interp_property_set(code_names(1), &
                                       coupling_name, &
                                       "tolerance",   &
                                       "double",      &
                                       "1e-3")

  !! Set user-defined interpolation function
  ! if (code_names(1) == "code1") then
  !   call CWP_Interp_function_set(code_names(1), &
  !                                coupling_name, &
  !                                field_name,    &
  !                                my_interpolation)
  ! endif

  !! Compute spatial interpolation weights
  call CWP_Spatial_interp_weights_compute(code_names(1), &
                                          coupling_name)


  !! Exchange interpolated field
  if (code_names(1) == "code1") then
    call CWP_Field_issend(code_names(1), &
                          coupling_name, &
                          field_name)
  else
    call CWP_Field_irecv(code_names(1), &
                         coupling_name, &
                         field_name)
  endif



  if (code_names(1) == "code1") then
    call CWP_Field_wait_issend(code_names(1), &
                               coupling_name, &
                               field_name)


    call CWP_Interp_src_data_get(code_names(1), &
                                 coupling_name, &
                                 field_name,    &
                                 0,             &
                                 n_elt_src,     &
                                 src_to_tgt_idx)
    call CWP_Interp_location_weights_get(code_names(1), &
                                         coupling_name, &
                                         field_name,    &
                                         0,             &
                                         interp_weights)
    if (debug) then
      write(iiunit,*) "-- CWP_Interp_src_data_get & CWP_Interp_location_weights_get --"
      write(iiunit,*) "n_elt_src      : ", n_elt_src
      write(iiunit,*) "src_to_tgt_idx : ", src_to_tgt_idx
      write(iiunit,*) "interp_weights : ", interp_weights
      write(iiunit,*) ""
    endif
  else
    call CWP_Field_wait_irecv(code_names(1), &
                              coupling_name, &
                              field_name)

    call CWP_Interp_tgt_data_get(code_names(1), &
                                 coupling_name, &
                                 field_name,    &
                                 0,             &
                                 n_elt_tgt,     &
                                 n_ref_tgt,     &
                                 ref_tgt,       &
                                 tgt_to_src_idx)
    if (debug) then
      write(iiunit,*) "-- CWP_Interp_tgt_data_get --"
      write(iiunit,*) "n_elt_tgt      : ", n_elt_tgt
      write(iiunit,*) "n_ref_tgt      : ", n_ref_tgt
      write(iiunit,*) "ref_tgt        : ", ref_tgt
      write(iiunit,*) "tgt_to_src_idx : ", tgt_to_src_idx
      write(iiunit,*) ""
    endif
  endif



  !! Check
  n_wrong = 0
  if (code_names(1) == "code2") then

    n_computed_tgts = CWP_N_computed_tgts_get(code_names(1), &
                                              coupling_name, &
                                              field_name,    &
                                              0)
    if (debug) then
      write(iiunit,*) n_computed_tgts, " computed tgts /", n_vtx
    endif

    computed_tgts => CWP_computed_tgts_get(code_names(1), &
                                           coupling_name, &
                                           field_name,    &
                                           0)

    do i = 1, n_computed_tgts
      ivtx = computed_tgts(i)

      distance = sqrt(sum(vtx_coord(:,ivtx) - field_data(3*i-2:3*i)))

      if (distance > 1.d-6) then
        n_wrong = n_wrong + 1
        if (debug) then
          write(iiunit,*) "error vtx", ivtx, "distance =", distance
          write(iiunit,*) "coord =", vtx_coord(:,ivtx), " recv =", field_data(3*i-2:3*i)
        endif
      endif

    enddo

    if (debug) then
      write(iiunit,*) "n_wrong =", n_wrong
    endif
  endif



  !! Global data
  global_data_name = "chocolatine"
  allocate(global_data(3,2))
  if (code_names(1) == "code1") then
    do i = 1,2
      global_data(:,i) = [i,2*i,3*i]
    enddo
    call CWP_Global_data_issend(code_names(1),    &
                                coupling_name,    &
                                global_data_name, &
                                global_data)
  else
    call CWP_Global_data_irecv(code_names(1),    &
                               coupling_name,    &
                               global_data_name, &
                               global_data)
  endif

  ! call MPI_Barrier(MPI_comm_world, ierr)

  if (code_names(1) == "code1") then
    call CWP_Global_data_wait_issend(code_names(1), &
                                     coupling_name, &
                                     global_data_name)
  else
    call CWP_Global_data_wait_irecv(code_names(1),    &
                                    coupling_name,    &
                                    global_data_name)
  endif

  call MPI_Barrier(MPI_comm_world, ierr)

  if (debug) then
    write(iiunit,*) "-- Global data --"
    write(iiunit,*) "size = ", size(global_data,1), size(global_data,2)
    do i = 1,size(global_data,2)
      write(iiunit,*) global_data(:,i)
    enddo
  endif

  deallocate(global_data)


  !! Part data
  if (debug) then
    write(iiunit,*) "-- Part data --"
  endif
  part_data_name = "All work and no play makes Jack a dull boy"

  if (code_names(1) == "code1") then
    n_part    = 2
    exch_type = CWP_PARTDATA_SEND
  else
    n_part    = 1
    exch_type = CWP_PARTDATA_RECV
  endif


  call PDM_pointer_array_create(gnum_elt,       &
                                n_part,         &
                                PDM_TYPE_G_NUM)

  call PDM_pointer_array_create(part_data,     &
                                n_part,        &
                                PDM_TYPE_INT)

  allocate(my_part(n_part))
  if (code_names(1) == "code1") then
    allocate(n_elt_part(n_part))
    n_elt_part = [3, 2]
  else
    allocate(n_elt_part(n_part))
    n_elt_part = [5]
  endif

  n_comp = 2
  do i = 1, n_part
    allocate(my_part(i)%g_num(n_elt_part(i)), &
             my_part(i)%data (n_elt_part(i) * n_comp))
  enddo

  if (code_names(1) == "code1") then
    my_part(1)%g_num = [1, 3, 5]
    my_part(2)%g_num = [2, 4]
  else
    my_part(1)%g_num = [1, 4, 2, 5, 3]
  endif

  do i = 1, n_part
    do j = 1, n_elt_part(i)
      my_part(i)%data(2*(j-1)+1) =   my_part(i)%g_num(j)
      my_part(i)%data(2*j      ) = 2*my_part(i)%g_num(j)
    enddo

    call PDM_pointer_array_part_set(gnum_elt,         &
                                    i-1,              &
                                    my_part(i)%g_num)

    call PDM_pointer_array_part_set(part_data,        &
                                    i-1,              &
                                    my_part(i)%data)
  enddo

  call CWP_Part_data_create(code_names(1),  &
                            coupling_name,  &
                            part_data_name, &
                            exch_type,      &
                            gnum_elt,       &
                            n_elt_part,     &
                            n_part)

  if (code_names(1) == "code1") then
    call CWP_Part_data_issend(code_names(1),  &
                              coupling_name,  &
                              part_data_name, &
                              n_comp,         &
                              part_data,      &
                              request)
  else
    call CWP_Part_data_irecv(code_names(1),  &
                             coupling_name,  &
                             part_data_name, &
                             n_comp,         &
                             part_data,      &
                             request)
  endif


  if (code_names(1) == "code1") then
    call CWP_Part_data_wait_issend(code_names(1),  &
                                   coupling_name,  &
                                   part_data_name, &
                                   request)
    if (debug) then
      do i = 1, n_part
        call PDM_pointer_array_part_get(part_data, &
                                        i-1,       &
                                        recv_data)
        do j = 1, n_elt_part(i)
          write(iiunit, *) my_part(i)%g_num(j), " sends    ", recv_data(2*(j-1)+1:2*j)
        enddo
      enddo
    endif
  else
    call CWP_Part_data_wait_irecv(code_names(1),  &
                                  coupling_name,  &
                                  part_data_name, &
                                  request)

    if (debug) then
      do i = 1, n_part
        call PDM_pointer_array_part_get(part_data, &
                                        i-1,       &
                                        recv_data)
        do j = 1, n_elt_part(i)
          write(iiunit, *) my_part(i)%g_num(j), " receives ", recv_data(2*(j-1)+1:2*j)
        enddo
      enddo
    endif
  endif

  call CWP_Part_data_del(code_names(1),  &
                         coupling_name,  &
                         part_data_name, &
                         exch_type)

  do i = 1, n_part
    deallocate(my_part(i)%g_num, &
               my_part(i)%data)
  enddo
  deallocate(n_elt_part)
  deallocate(my_part)

  call PDM_pointer_array_free(part_data)
  call PDM_pointer_array_free(gnum_elt)


  if (debug) then
    close(iiunit)
  endif

  call CWP_Time_step_end(code_names(1))

  call CWP_Visu_end(code_names(1), &
                    coupling_name)

  !! Delete interface mesh
  call CWP_Mesh_interf_del(code_names(1), &
                           coupling_name)

  !! Delete coupling
  call CWP_Cpl_Del(code_names(1), &
                   coupling_name)


  !! Free memory
  deallocate(code_names,         &
             coupled_code_names, &
             is_coupled_rank,    &
             time_init,          &
             intra_comms,        &
             g_code_names)
  deallocate(vtx_coord, connec_idx, connec)
  deallocate(field_data)
  deallocate(code_list,     &
             loc_code_list, &
             f_param_names)

  !! Finalize
  call CWP_Finalize()
  call MPI_Finalize(ierr)


contains

subroutine gen_mesh(xmin, xmax, &
                    ymin, ymax, &
                    rand_level, &
                    n_vtx_seg,  &
                    n_vtx,      &
                    vtx_coord,  &
                    n_elt,      &
                    connec_idx, &
                    connec,     &
                    vtx_g_num,  &
                    elt_g_num)
  implicit none

  double precision, intent(in)  :: xmin, xmax, ymin, ymax
  double precision, intent(in)  :: rand_level
  integer,          intent(in)  :: n_vtx_seg
  integer,          intent(out) :: n_vtx
  double precision, pointer     :: vtx_coord(:,:)
  integer,          intent(out) :: n_elt
  integer(c_int),   pointer     :: connec_idx(:)
  integer(c_int),   pointer     :: connec(:)
  integer(c_long),  pointer     :: vtx_g_num(:)
  integer(c_long),  pointer     :: elt_g_num(:)

  double precision              :: step_x, step_y, r(2)
  integer                       :: i, j, k, l


  !! Vertices !!
  step_x = (xmax - xmin) / (n_vtx_seg - 1)
  step_y = (ymax - ymin) / (n_vtx_seg - 1)

  n_vtx = n_vtx_seg * n_vtx_seg
  allocate(vtx_coord(3,n_vtx))

  k = 1
  do j = 1, n_vtx_seg
    do i = 1, n_vtx_seg
      call random_number(r)
      r = 2*r - 1
      vtx_coord(1,k) = xmin + step_x*(i-1)
      if (i > 1 .and. i < n_vtx_seg) then
        vtx_coord(1,k) = vtx_coord(1,k) + step_x*rand_level*r(1)
      end if
      vtx_coord(2,k) = ymin + step_y*(j-1)
      if (j > 1 .and. j < n_vtx_seg) then
        vtx_coord(2,k) = vtx_coord(2,k) + step_x*rand_level*r(2)
      end if
      vtx_coord(3,k) = 0.d0
      k = k+1
    end do
  end do

  !! Elements !!
  n_elt = (n_vtx_seg - 1) * (n_vtx_seg - 1)
  allocate(connec_idx(n_elt+1), connec(4*n_elt))
  connec_idx(1) = 0
  k = 1
  l = 1
  do j = 1, n_vtx_seg - 1
    do i = 1, n_vtx_seg - 1
      connec(l) = i + n_vtx_seg*(j-1)
      l = l+1
      connec(l) = i + n_vtx_seg*(j-1) + 1
      l = l+1
      connec(l) = i + n_vtx_seg*j     + 1
      l = l+1
      connec(l) = i + n_vtx_seg*j
      l = l+1

      connec_idx(k+1) = connec_idx(k) + 4
      k = k+1
    end do
  end do


  !! Global numbers !!
  vtx_g_num => null()
  elt_g_num => null()
  ! allocate(vtx_g_num(n_vtx), elt_g_num(n_elt))
  ! do i = 1, n_vtx
  !   vtx_g_num(i) = i
  ! end do
  ! do i = 1, n_elt
  !   elt_g_num(i) = i
  ! end do

end subroutine gen_mesh


  subroutine my_interpolation(local_code_name,          &
                              cpl_id,                   &
                              field_id,                 &
                              i_part,                   &
                              spatial_interp_algorithm, &
                              storage,                  &
                              c_buffer_in,              &
                              c_buffer_out)             &
    bind(c)
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
    integer(kind = c_int)             :: i_part
    integer(kind = c_int)             :: spatial_interp_algorithm
    integer(kind = c_int)             :: storage
    type(c_ptr), value                :: c_buffer_in
    type(c_ptr), value                :: c_buffer_out

    integer(kind = c_int)             :: n_components
    integer(kind = c_int)             :: n_elt_src
    integer(kind = c_int), pointer    :: src_to_tgt_idx(:) => null()
    real(kind = c_double), pointer    :: buffer_in(:)      => null()
    real(kind = c_double), pointer    :: buffer_out(:)     => null()
    integer                           :: i, j, k

    print *, ">> my_interpolation"


    n_components = CWP_Interp_field_n_components_get(local_code_name, &
                                                     cpl_id,          &
                                                     field_id)

    call CWP_Interp_src_data_get(local_code_name, &
                                 cpl_id,          &
                                 field_id,        &
                                 i_part,          &
                                 n_elt_src,       &
                                 src_to_tgt_idx)

    call c_f_pointer(c_buffer_in,  buffer_in,  [n_elt_src])
    call c_f_pointer(c_buffer_out, buffer_out, [src_to_tgt_idx(n_elt_src+1)])

    do i = 1, n_elt_src
      print *, "i = ", i
      do j = src_to_tgt_idx(i)+1, src_to_tgt_idx(i+1)
        print *, "  j = ", j
        do k = 1, n_components
          buffer_out(n_components*(j-1)+k) = buffer_in(n_components*(i-1)+k)
        enddo
      enddo
    enddo

  end subroutine my_interpolation



  subroutine mon_toto(i) bind(c)
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value, intent(in) :: i
    print *, "Toto i = ", i
    ! print *, "Coucou toto!"
  end subroutine mon_toto







end program testf
