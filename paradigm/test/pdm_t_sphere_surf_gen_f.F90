#include "pdm_configf.h"

program testf

  use iso_c_binding
  use pdm
#ifdef PDM_HAVE_FORTRAN_MPI_MODULE
  use mpi
#endif
  use pdm_io
  use pdm_writer
  use pdm_sphere_surf_gen
  use pdm_fortran

  implicit none

#ifndef PDM_HAVE_FORTRAN_MPI_MODULE
  include "mpif.h"
#endif

  !-----------------------------------------------------------
  integer, parameter                :: comm        = MPI_COMM_WORLD
  integer, parameter                :: n_part      = 1
  integer, parameter                :: part_method = 1

  integer(pdm_g_num_s)              :: n        = 5
  double precision                  :: x_center = 1.d0
  double precision                  :: y_center = 2.d0
  double precision                  :: z_center = 3.d0
  double precision                  :: radius   = 4.d0

  integer(pdm_l_num_s), pointer     :: pn_vtx(:)  => null()
  type(PDM_pointer_array_t)         :: pvtx_coord
  type(PDM_pointer_array_t)         :: pvtx_ln_to_gn
  integer(pdm_l_num_s), pointer     :: pn_face(:) => null()
  type(PDM_pointer_array_t)         :: pface_vtx_idx
  type(PDM_pointer_array_t)         :: pface_vtx
  type(PDM_pointer_array_t)         :: pface_ln_to_gn

  double precision,     pointer     :: vtx_coord(:)     => null()
  double precision,     pointer     :: vtx_coord2(:,:)  => null()
  integer(pdm_g_num_s), pointer     :: vtx_ln_to_gn(:)  => null()
  ! integer(pdm_l_num_s), pointer     :: face_vtx_idx(:)  => null()
  integer(pdm_l_num_s), pointer     :: face_vtx(:)      => null()
  integer(pdm_g_num_s), pointer     :: face_ln_to_gn(:) => null()

  ! Writer
  type(c_ptr)                       :: wrt
  integer                           :: id_geom
  integer                           :: id_block

  ! MPI
  integer                           :: code
  integer                           :: i_rank
  integer                           :: n_rank

  integer                           :: ipart
  !-----------------------------------------------------------

  call mpi_init(code)
  call mpi_comm_rank(comm, i_rank, code)
  call mpi_comm_size(comm, n_rank, code)

  if (i_rank .eq. 0) then
    write(*, *) "-- Generate mesh"
  end if
  call PDM_sphere_surf_icosphere_gen_part(comm,           &
                                          n,              &
                                          x_center,       &
                                          y_center,       &
                                          z_center,       &
                                          radius,         &
                                          n_part,         &
                                          part_method,    &
                                          pn_vtx,         &
                                          pvtx_coord,     &
                                          pvtx_ln_to_gn,  &
                                          pn_face,        &
                                          pface_vtx_idx,  &
                                          pface_vtx,      &
                                          pface_ln_to_gn)


  if (i_rank .eq. 0) then
    write(*, *) "-- Write mesh"
  end if
  call PDM_writer_create(wrt,                    &
                         "Ensight",              &
                         PDM_WRITER_FMT_ASCII,   &
                         PDM_WRITER_TOPO_CST,    &
                         PDM_WRITER_OFF,         &
                         "sphere_surf",          &
                         "sphere_surf",          &
                         comm,                   &
                         PDM_IO_KIND_MPI_SIMPLE, &
                         1.d0,                   &
                         "")

  call PDM_writer_geom_create(wrt,     &
                              id_geom, &
                              "mesh",  &
                              n_part)

  call PDM_writer_step_beg(wrt, 0.d0)

  !  Write geometry
  do ipart = 1, n_part
    call PDM_pointer_array_part_get(pvtx_ln_to_gn, &
                                    ipart-1,       &
                                    vtx_ln_to_gn)
    call PDM_pointer_array_part_get(pvtx_coord, &
                                    ipart-1,    &
                                    vtx_coord)
    allocate(vtx_coord2(3, pn_vtx(ipart)))
    vtx_coord2(1,:) = vtx_coord(1:3*pn_vtx(ipart):3)
    vtx_coord2(2,:) = vtx_coord(2:3*pn_vtx(ipart):3)
    vtx_coord2(3,:) = vtx_coord(3:3*pn_vtx(ipart):3)
    call PDM_pointer_array_part_get(pvtx_coord, &
                                    ipart-1,    &
                                    vtx_coord)

    call PDM_writer_geom_coord_set(wrt,                &
                                   id_geom,            &
                                   ipart-1,            &
                                   pn_vtx(ipart),      &
                                   vtx_coord2,         &
                                   vtx_ln_to_gn,       &
                                   PDM_OWNERSHIP_KEEP)


    call PDM_writer_geom_bloc_add(wrt,                  &
                                  id_geom,              &
                                  2,                    & ! PDM_MESH_NODAL_TRIA3
                                  PDM_OWNERSHIP_USER,   &
                                  id_block)

    call PDM_pointer_array_part_get(pface_ln_to_gn, &
                                    ipart-1,        &
                                    face_ln_to_gn)
    call PDM_pointer_array_part_get(pface_vtx, &
                                    ipart-1,   &
                                    face_vtx)

    call PDM_writer_geom_bloc_std_set(wrt,            &
                                      id_geom,        &
                                      id_block,       &
                                      ipart-1,        &
                                      pn_face(ipart), &
                                      face_vtx,       &
                                      face_ln_to_gn)
  end do

  call PDM_writer_geom_write(wrt,     &
                             id_geom)

  call PDM_writer_step_end(wrt)

  call PDM_writer_free(wrt)


  !  Free memory
  call pdm_fortran_free_c(c_loc(pn_vtx))
  call pdm_fortran_free_c(c_loc(pn_face))
  call PDM_pointer_array_free(pvtx_coord)
  call PDM_pointer_array_free(pvtx_ln_to_gn)
  call PDM_pointer_array_free(pface_vtx_idx)
  call PDM_pointer_array_free(pface_vtx)
  call PDM_pointer_array_free(pface_ln_to_gn)


  if (i_rank .eq. 0) then
    write(*, *) "-- End"
  end if

  call mpi_finalize(code)

end program testf
