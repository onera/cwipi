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

  double precision, dimension(:), pointer     :: xyz_dest
  integer(c_long), pointer, dimension(:)      :: pts_g_num => null()

  integer, pointer, dimension(:)              :: connec_idx
  integer, pointer, dimension(:)              :: connec
  integer(c_long), pointer, dimension(:)      :: elt_g_num => null()
  integer(c_int)                              :: id_block

  character(len = 99)                         :: send_field_name
  character(len = 99)                         :: recv_field_name
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

  print *, "FORTRAN - CWP_Init : OK"

  ! Create the coupling :
  ! CWP_DYNAMIC_MESH_DEFORMABLE allows us to take into account the modifications
  ! to the mesh over the coupling steps.
  coupling_name     = "code1_code2"

  allocate(coupled_code_names(n_code))

  coupled_code_names(1) = "code2"

  n_part = 1;
  call CWP_Cpl_create(code_names(1),                                         &
                      coupling_name,                                         &
                      coupled_code_names(1),                                 &
                      CWP_INTERFACE_SURFACE,                                 &
                      CWP_COMM_PAR_WITH_PART,                                &
                      CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, &
                      n_part,                                                &
                      CWP_DYNAMIC_MESH_DEFORMABLE,                           &
                      CWP_TIME_EXCH_USER_CONTROLLED)

  print *, "FORTRAN - CWP_Cpl_create : OK"

  ! Set coupling visualisation:
  call CWP_Visu_set(code_names(1),           &
                    coupling_name,           &
                    1,                       &
                    CWP_VISU_FORMAT_ENSIGHT, &
                    "text")

  print *, "FORTRAN - CWP_Visu_set : OK"

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
  allocate(send_field_data(n_elt))
  allocate(recv_field_data(n_elt))
  allocate(xyz_dest(3 * n_elt))
  allocate(uncomputed_tgts(n_elt))

  call grid_mesh(xmin, &
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

  print *, "FORTRAN - grid_mesh : OK"

  send_field_name = "girafe"
  recv_field_name = "chinchilla"

  ! Interations :
  time  = 0.0d0
  dt    = 0.1d0
  omega = 2.0d0*acos(-1.0d0)*freq

  do it = itdeb, itend

    time = (it-itdeb)*dt

    do i = 1, n_vtx
      coords((i-1)*3+3) = (coords((i-1)*3+1)*coords((i-1)*3+1)+  &
                           coords((i-1)*3+2)*coords((i-1)*3+2))* &
                           ampl*cos(omega*time+phi)
    enddo

    do  i = 1, n_elt
      send_field_data(i) = 0.
      do j = connec_idx(i)+1, connec_idx(i+1)
        send_field_data(i) = send_field_data(i) + coords((connec(j)-1)*3+3)
      enddo
      send_field_data(i) = send_field_data(i)/ (connec_idx(i+1) - connec_idx(i))
    enddo

    if (it == itdeb) then

      ! Set the mesh vertices coordinates :
      call CWP_Mesh_interf_vtx_set(code_names(1), &
                                   coupling_name, &
                                   0,             &
                                   n_vtx,         &
                                   coords,        &
                                   vtx_g_num)

      print *, "FORTRAN - CWP_Mesh_interf_vtx_set : OK"

      ! Set the mesh polygons connectivity :
      id_block = CWP_Mesh_interf_block_add(code_names(1),       &
                                           coupling_name,       &
                                           CWP_BLOCK_FACE_POLY)

      print *, "FORTRAN - CWP_Mesh_interf_block_add : OK"

      call CWP_Mesh_interf_f_poly_block_set(code_names(1), &
                                            coupling_name, &
                                            0,             &
                                            id_block,      &
                                            n_elt,         &
                                            connec_idx,    &
                                            connec,        &
                                            elt_g_num)

      print *, "FORTRAN - CWP_Mesh_interf_f_poly_block_set : OK"

      ! Finalize mesh :
      call CWP_Mesh_interf_finalize(code_names(1), &
                                    coupling_name)

      print *, "FORTRAN - CWP_Mesh_interf_finalize : OK"

      ! Create field :
      call CWP_Field_create(code_names(1),                &
                            coupling_name,                &
                            send_field_name,              &
                            CWP_DOUBLE,                   &
                            CWP_FIELD_STORAGE_INTERLACED, &
                            n_components,                 &
                            CWP_DOF_LOCATION_USER,        &
                            CWP_FIELD_EXCH_SEND,          &
                            CWP_STATUS_ON)

      print *, "FORTRAN - CWP_Field_create send : OK"

      call CWP_Field_create(code_names(1),                &
                            coupling_name,                &
                            recv_field_name,              &
                            CWP_DOUBLE,                   &
                            CWP_FIELD_STORAGE_INTERLACED, &
                            n_components,                 &
                            CWP_DOF_LOCATION_USER,        &
                            CWP_FIELD_EXCH_RECV,          &
                            CWP_STATUS_ON)

      print *, "FORTRAN - CWP_Field_create recv : OK"

      ! Set interpolation property :
      call CWP_Spatial_interp_property_set(code_names(1), &
                                           coupling_name, &
                                           "tolerance",   &
                                           "double",      &
                                           "0.1")

      print *, "FORTRAN - CWP_Spatial_interp_property_set : OK"

    else

      ! Update user defined degrees of freedom :
      ! When CWP_DOF_LOCATION_USER, calling CWP_User_tgt_pts_set is mandatory.

      do i = 1, n_elt
         xyz_dest((i-1)*3+1) = 0.
         xyz_dest((i-1)*3+2) = 0.
         xyz_dest((i-1)*3+3) = 0.
         do j = connec_idx(i)+1, connec_idx(i+1)
             xyz_dest((i-1)*3+1) = xyz_dest((i-1)*3+1) + coords((connec(j)-1)*3+1)
             xyz_dest((i-1)*3+2) = xyz_dest((i-1)*3+2) + coords((connec(j)-1)*3+2)
             xyz_dest((i-1)*3+3) = xyz_dest((i-1)*3+3) + coords((connec(j)-1)*3+3)
         enddo
         xyz_dest((i-1)*3+1) = xyz_dest((i-1)*3+1) / (connec_idx(i+1) - connec_idx(i))
         xyz_dest((i-1)*3+2) = xyz_dest((i-1)*3+2) / (connec_idx(i+1) - connec_idx(i))
         xyz_dest((i-1)*3+3) = xyz_dest((i-1)*3+3) / (connec_idx(i+1) - connec_idx(i))
      enddo

      call CWP_User_tgt_pts_set(code_names(1), &
                                coupling_name, &
                                0,             &
                                n_elt,         &
                                xyz_dest,      &
                                pts_g_num)

      print *, "FORTRAN - CWP_User_tgt_pts_set : OK"

      ! Finalize mesh :
      call CWP_Mesh_interf_finalize(code_names(1), &
                                    coupling_name)

      print *, "FORTRAN - CWP_Mesh_interf_finalize : OK"

    endif

    ! Compute interpolation weights :
    call CWP_Spatial_interp_weights_compute(code_names(1), &
                                            coupling_name)

    print *, "FORTRAN - CWP_Spatial_interp_weights_compute : OK"

    call CWP_Field_data_set(code_names(1),        &
                            coupling_name,        &
                            recv_field_name,      &
                            0,                    &
                            CWP_FIELD_MAP_TARGET, &
                            recv_field_data)

    print *, "FORTRAN - CWP_Field_data_set recv : OK"

    call CWP_Field_irecv(code_names(1), &
                         coupling_name, &
                         recv_field_name)

    print *, "FORTRAN - CWP_Field_irecv : OK"

    ! Set and exchange the field values :
    call CWP_Field_data_set(code_names(1),        &
                            coupling_name,        &
                            send_field_name,      &
                            0,                    &
                            CWP_FIELD_MAP_SOURCE, &
                            send_field_data)

    print *, "FORTRAN - CWP_Field_data_set send : OK"

    call CWP_Field_issend(code_names(1), &
                          coupling_name, &
                          send_field_name)

    print *, "FORTRAN - CWP_Field_issend : OK"

    call CWP_Field_wait_irecv(code_names(1), &
                              coupling_name, &
                              recv_field_name)

    print *, "FORTRAN - CWP_Field_wait_irecv : OK"

    call CWP_Field_wait_issend(code_names(1), &
                               coupling_name, &
                               send_field_name)

    print *, "FORTRAN - CWP_Field_wait_issend : OK"

    ! Check interpolation :
    n_uncomputed_tgts = CWP_N_uncomputed_tgts_get(code_names(1),   &
                                                  coupling_name,   &
                                                  recv_field_name, &
                                                  0)

    print *, "FORTRAN - CWP_N_uncomputed_tgts_get : OK"

    if (n_uncomputed_tgts /= 0) then
      uncomputed_tgts => CWP_Uncomputed_tgts_get(code_names(1),   &
                                                 coupling_name,   &
                                                 recv_field_name, &
                                                 0)

      print *, "FORTRAN - CWP_Uncomputed_tgts_get : OK"
    endif

  enddo

  ! Delete field :
  call CWP_Field_Del(code_names(1),   &
                     coupling_name,   &
                     send_field_name)

  print *, "FORTRAN - CWP_Field_Del send : OK"

  call CWP_Field_Del(code_names(1),   &
                     coupling_name,   &
                     recv_field_name)

  print *, "FORTRAN - CWP_Field_Del recv : OK"

  ! Delete Mesh :
  call CWP_Mesh_interf_del(code_names(1), &
                           coupling_name)

  print *, "FORTRAN - CWP_Mesh_interf_del : OK"

  ! Delete the coupling :
  call CWP_Cpl_Del(code_names(1), &
                   coupling_name)

  print *, "FORTRAN - CWP_Cpl_Del : OK"

  ! free
  deallocate(code_names)
  deallocate(is_coupled_rank)
  deallocate(time_init)
  deallocate(intra_comms)
  deallocate(coupled_code_names)
  deallocate(coords)
  deallocate(connec_idx)
  deallocate(connec)
  deallocate(send_field_data)
  deallocate(recv_field_data)
  deallocate(xyz_dest)
  deallocate(uncomputed_tgts)

  ! Finalize CWIPI :
  call CWP_Finalize()

  print *, "FORTRAN - CWP_Finalize : OK"

  ! Finalize MPI :
  call MPI_Finalize(ierr)

contains
   ! Return a random number beteween -1 and 1
   function random01() &
     result (resultat)

     implicit none


     integer          :: signe, rsigna, rsignb, rd
     double precision :: u
     double precision :: resultat

     call random_number(u)

     rsigna = floor(32768*u) ! RAND_MAX + 1
     rsignb = floor(32768*u)
     signe = (rsigna - rsignb) / abs(rsigna - rsignb)
     rd = floor(32768*u)
     resultat = signe*(rd/32767) ! RAND_MAX

   end function random01

   ! Create mesh
   subroutine grid_mesh(xmin, &
                        xmax, &
                        ymin, &
                        ymax, &
                        randLevel, &
                        n_vtx_seg, &
                        n_partition, &
                        coords,  &
                        connec_idx,&
                        connec,&
                        comm)

     implicit none

     double precision, intent(in) :: xmax, xmin, ymax, ymin, randLevel
     integer,          intent(in) :: n_vtx_seg, n_partition, comm

     double precision, pointer, dimension(:), intent(out)  :: coords
     integer, pointer, dimension(:), intent(out)           :: connec_idx
     integer, pointer, dimension(:), intent(out)           :: connec

     integer          :: i, j, k, ii, jj, p1, p2, p3, p4
     integer          :: ierr, i_rank, n_rank
     integer          :: nBoundVerticesSeg, nBoundVertices, nVertex, nElts
     double precision :: lX, lY, rd01
     double precision :: deltaU, deltaV, u, v, randU, randV

     double precision, pointer, dimension(:) :: boundRanks, boundRank

     call MPI_Comm_rank(comm, i_rank, ierr)
     call MPI_Comm_size(comm, n_rank, ierr)

     lX = (xmax - xmin) / n_partition
     lY = (ymax - ymin) / n_partition

     ! Compute local partition bounds with random level

     nBoundVerticesSeg = n_partition + 1
     nBoundVertices = nBoundVerticesSeg * nBoundVerticesSeg

     allocate(boundRanks(3*nBoundVertices))

     if (i_rank == 0) then
       do j=1,nBoundVerticesSeg
         do i=1,nBoundVerticesSeg
           boundRanks(3 * ((j-1) * nBoundVerticesSeg + (i-1)) + 1) = xmin + (i-1) * lX
           boundRanks(3 * ((j-1) * nBoundVerticesSeg + (i-1)) + 2) = ymin + (j-1) * lY
           boundRanks(3 * ((j-1) * nBoundVerticesSeg + (i-1)) + 3) = 0.
           if (j /= 1 .and. j /= nBoundVerticesSeg) then
             rd01 =  random01()
             boundRanks(3 * ((j-1) * nBoundVerticesSeg + (i-1)) + 2) =  &
             boundRanks(3 * ((j-1) * nBoundVerticesSeg + (i-1)) + 2) &
             + rd01 * randLevel * lY
           endif
           if (i /= 1 .and. i /= nBoundVerticesSeg) then
             rd01 =  random01()
             boundRanks(3 * ((j-1) * nBoundVerticesSeg + (i-1)) + 1) = &
             boundRanks(3 * ((j-1) * nBoundVerticesSeg + (i-1)) + 1) + &
             rd01 * randLevel * lX
           endif
         enddo
       enddo
     endif

     call MPI_Bcast(boundRanks, 3 * nBoundVertices, MPI_DOUBLE_PRECISION, 0, comm)

     allocate(boundRank(3*4))

     ii = modulo(i_rank, n_partition)
     jj = i_rank / n_partition

     p1 = (nBoundVerticesSeg * jj)       + ii
     p2 = (nBoundVerticesSeg * jj)       + ii + 1
     p3 = (nBoundVerticesSeg * (jj + 1)) + ii + 1
     p4 = (nBoundVerticesSeg * (jj + 1)) + ii

     boundRank(0 * 3 + 1) = boundRanks(3 * p1 + 1)
     boundRank(0 * 3 + 2) = boundRanks(3 * p1 + 2)
     boundRank(0 * 3 + 3) = boundRanks(3 * p1 + 3)

     boundRank(1 * 3 + 1) = boundRanks(3 * p2 + 1)
     boundRank(1 * 3 + 2) = boundRanks(3 * p2 + 2)
     boundRank(1 * 3 + 3) = boundRanks(3 * p2 + 3)

     boundRank(2 * 3 + 1) = boundRanks(3 * p3 + 1)
     boundRank(2 * 3 + 2) = boundRanks(3 * p3 + 2)
     boundRank(2 * 3 + 3) = boundRanks(3 * p3 + 3)

     boundRank(3 * 3 + 1) = boundRanks(3 * p4 + 1)
     boundRank(3 * 3 + 2) = boundRanks(3 * p4 + 2)
     boundRank(3 * 3 + 3) = boundRanks(3 * p4 + 3)

     deallocate(boundRanks)

     ! Number of vertices and elements in the partition

     nVertex = n_vtx_seg * n_vtx_seg
     nElts   = (n_vtx_seg - 1) * (n_vtx_seg - 1)

     ! Define coordinates

     deltaU = 2.0/(n_vtx_seg - 1)
     deltaV = 2.0/(n_vtx_seg - 1)
     u = -1
     v = -1
     do j=1,n_vtx_seg
       do i=1,n_vtx_seg
         randU = u
         randV = v
         if ((i /= 1) .and. (j /= 1) .and. (j /= n_vtx_seg) .and. (i /= n_vtx_seg)) then
           rd01 =  random01()
           randU = randU + rd01 * randLevel * deltaU
           rd01 =  random01()
           randV = randV + rd01 * randLevel * deltaV
         endif

         coords(3 * ((j-1) * n_vtx_seg + (i-1)) + 1) = &
         0.25 * ((1 - randU - randV + randU * randV) * boundRank(0 * 3 + 1) + &
                 (1 + randU - randV - randU * randV) * boundRank(1 * 3 + 1) + &
                 (1 + randU + randV + randU * randV) * boundRank(2 * 3 + 1) + &
                 (1 - randU + randV - randU * randV) * boundRank(3 * 3 + 1) )

         coords(3 * ((j-1) * n_vtx_seg + (i-1)) + 2) = &
         0.25 * ((1 - randU - randV + randU * randV) * boundRank(0 * 3 + 2) + &
                 (1 + randU - randV - randU * randV) * boundRank(1 * 3 + 2) + &
                 (1 + randU + randV + randU * randV) * boundRank(2 * 3 + 2) + &
                 (1 - randU + randV - randU * randV) * boundRank(3 * 3 + 2) )

         coords(3 * ((j-1) * n_vtx_seg + (i-1)) + 3) = 0.

         u = u + deltaU

       enddo

       v = v + deltaV
       u = -1

     enddo

     deallocate(boundRank)

     ! Define connectivity

     connec_idx(1) = 0;
     do i=2,(nElts+1)
       connec_idx(i) = connec_idx(i-1) + 4
     enddo


     k = 0;
     do j=1,(n_vtx_seg-1)
       do i=1,(n_vtx_seg-1)
         connec(4 * k + 1) = (j - 1) * n_vtx_seg + (i-1)     + 1
         connec(4 * k + 2) = (j - 1) * n_vtx_seg + (i-1) + 1 + 1
         connec(4 * k + 3) =  j      * n_vtx_seg + (i-1) + 1 + 1
         connec(4 * k + 4) =  j      * n_vtx_seg + (i-1)     + 1
         k = k + 1
       enddo
     enddo

   end subroutine grid_mesh

end program fortran_new_api_deformable_sol
