program testf
  use mpi
  use cwp

  implicit none

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
  double precision,     pointer :: vtx_coord(:)  => null()
  integer(c_int),       pointer :: connec_idx(:) => null()
  integer(c_int),       pointer :: connec(:)     => null()
  integer(c_long),      pointer :: vtx_g_num(:)  => null()
  integer(c_long),      pointer :: elt_g_num(:)  => null()
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

  integer(c_int)                :: n_vtx2, n_elt2
  integer(c_int),       pointer :: connec_idx2(:) => null()
  integer(c_int),       pointer :: connec2(:)     => null()
  integer(c_long),      pointer :: vtx_g_num2(:)  => null()
  integer(c_long),      pointer :: elt_g_num2(:)  => null()

  integer                       :: i, ivtx, n_wrong
  double precision              :: distance

  character                     :: strnum
  integer                       :: fid = 13
  logical                       :: debug = .true.

  !--> character array getters
  character(256), allocatable :: code_list(:)
  character(256), allocatable :: loc_code_list(:)
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
    open(unit=fid, file="fortran_new_api_surf_"//strnum//".log", action='write')
  endif


  !! Initialize CWIPI
  n_code = 1
  n_part = 1

  allocate(code_names(n_code),         &
           coupled_code_names(n_code), &
           is_coupled_rank(n_code),    &
           time_init(n_code),          &
           intra_comms(n_code))

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

  call CWP_Init(MPI_comm_world,  &
                n_code,          &
                code_names,      &
                is_coupled_rank, &
                time_init,       &
                intra_comms)

  !--> character array getters
  code_list     = CWP_Codes_list_get()
  loc_code_list = CWP_Loc_codes_list_get()

  !--> n_total_codes = 2
  do i=1,2
    print *, i_rank, " --> ", "code_list(", i, ") :", code_list(i)
  end do

  print *, i_rank, " --> ", "loc_code_list(", 1, ") :", loc_code_list(1)

  !-->>
  ! print *, "n_code =", CWP_Codes_nb_get()
  ! print *, "n_local_code =", CWP_Loc_codes_nb_get()
  ! print *, "state =", CWP_State_get(code_names(1))
  ! call CWP_Properties_dump()

  toto = 123456
  if (code_names(1) == "code1") then
    call CWP_Param_add("code1", "toto", toto)
  endif

  n_param = CWP_Param_n_get("code1", &
                            CWP_INT)
  print *, "n_param =", n_param

  ! print *, code_names(1), ", param is? :", CWP_Param_is("code1",  &
  !                                                       "toto", &
  !                                                       CWP_INT)

  if (code_names(1) == "code1") then
    call CWP_Param_get("code1",  &
                       "toto",   &
                       CWP_INT,  &
                       param_value)
    call c_f_pointer(param_value, tata)
    print *, "tata =", loc(tata)
  endif
  !<<--

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
    field_data(:) = vtx_coord(:)
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
  else
    call CWP_Field_wait_irecv(code_names(1), &
                              coupling_name, &
                              field_name)
  endif



  !! Check
  n_wrong = 0
  if (code_names(1) == "code2") then

    n_computed_tgts = CWP_N_computed_tgts_get(code_names(1), &
                                              coupling_name, &
                                              field_name,    &
                                              0)

    print *, n_computed_tgts, " computed tgts /", n_vtx

    computed_tgts => CWP_computed_tgts_get(code_names(1), &
                                           coupling_name, &
                                           field_name,    &
                                           0)

    do i = 1, n_computed_tgts
      ivtx = computed_tgts(i)

      distance = sqrt(sum(vtx_coord(3*ivtx-2:3*ivtx) - field_data(3*i-2:3*i)))

      if (distance > 1.d-6) then
        n_wrong = n_wrong + 1
        print *, "error vtx", ivtx, "distance =", distance
        print *, "coord =", vtx_coord(3*ivtx-2:3*ivtx), " recv =", field_data(3*i-2:3*i)
      endif

    enddo

    print *, "n_wrong =", n_wrong
  endif


  if (debug) then
    close(fid)
  endif



  !! Delete interface mesh
  call CWP_Mesh_interf_del(code_names(1), &
                           coupling_name)

  !! Delete coupling
  call CWP_Cpl_Del(code_names(1), &
                   coupling_name)

  !! Free memory
  deallocate(code_names, coupled_code_names, is_coupled_rank, time_init, intra_comms)
  deallocate(vtx_coord, connec_idx, connec)
  deallocate(field_data)

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
  double precision, pointer     :: vtx_coord(:)
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
  allocate(vtx_coord(3*n_vtx))

  k = 1
  l = 1
  do j = 1, n_vtx_seg
    do i = 1, n_vtx_seg
      call random_number(r)
      r = 2*r - 1
      vtx_coord(l) = xmin + step_x*(i-1)
      if (i > 1 .and. i < n_vtx_seg) then
        vtx_coord(l) = vtx_coord(l) + step_x*rand_level*r(1)
      end if
      l = l+1
      vtx_coord(l) = ymin + step_y*(j-1)
      if (j > 1 .and. j < n_vtx_seg) then
        vtx_coord(l) = vtx_coord(l) + step_x*rand_level*r(2)
      end if
      l = l+1
      vtx_coord(l) = 0.d0
      l = l+1
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


end program testf
