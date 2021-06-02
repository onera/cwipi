!-----------------------------------------------------------------------------
! This file is part of the ParaDiGM library.
!
! Copyright (C) 2020  ONERA
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

program testf

  use pdm
  use mpi
  use pdm_mesh_location
  use iso_c_binding


  implicit none

  !
  ! About point cloud
  !

  integer, parameter :: n_point_cloud = 1
  integer, parameter :: n_part_cloud = 1
  integer, parameter :: i_point_cloud = 0
  integer, parameter :: n_points_into_cloud = 2

  double precision, pointer :: coords_cloud(:) ! pointer or allocatble, target
  integer (kind = pdm_g_num_s), pointer :: gnum_cloud(:) ! pointer or allocatble, target

  type(c_ptr)  :: cptr_coords_cloud
  type(c_ptr)  :: cptr_gnum_cloud

  !
  ! About mesh
  !

  integer, parameter :: n_part_mesh = 1

  integer, parameter :: n_cell = 2
  integer (c_int), pointer :: cell_face_idx(:)
  integer (c_int), pointer :: cell_face(:)
  integer (kind = pdm_g_num_s), pointer :: gnum_cell(:)
  type(c_ptr)  :: cptr_cell_face_idx
  type(c_ptr)  :: cptr_cell_face
  type(c_ptr)  :: cptr_gnum_cell

  integer, parameter :: n_face = 11
  integer (c_int), pointer :: face_vtx_idx(:)
  integer (c_int), pointer :: face_vtx(:)
  integer (kind = pdm_g_num_s), pointer :: gnum_face(:)
  type(c_ptr)  :: cptr_face_vtx_idx
  type(c_ptr)  :: cptr_face_vtx
  type(c_ptr)  :: cptr_gnum_face

  integer, parameter :: n_vtx = 12
  double precision, pointer :: coords_vtx(:) ! pointer or allocatble, target
  integer (kind = pdm_g_num_s), pointer :: gnum_vtx(:) ! pointer or allocatble, target
  type(c_ptr)  :: cptr_coords_vtx
  type(c_ptr)  :: cptr_gnum_vtx

  integer :: i
  integer, parameter :: i_part_cloud = 0
  integer, parameter :: i_part_mesh = 0

  type(c_ptr)              :: ml

  ! integer, parameter :: partial = 0 ! Put 1 to keep results when the subroutine closest_points_free is

  integer :: code
  integer :: i_rank
  integer :: n_rank

  !
  ! results
  !

  type(c_ptr) :: cptr_location
  type(c_ptr) :: cptr_dist2
  type(c_ptr) :: cptr_projected_coords

  integer (kind = pdm_g_num_s), pointer :: location(:)
  double precision, pointer :: dist2(:)
  double precision, pointer :: projected_coords(:)

  integer :: n_located
  integer :: n_unlocated
  integer (c_int), pointer :: located(:)
  integer (c_int), pointer :: unlocated(:)
  type(c_ptr)  :: cptr_located
  type(c_ptr)  :: cptr_unlocated

  type(c_ptr) :: cptr_elt_pts_inside_idx
  type(c_ptr) :: cptr_points_gnum
  type(c_ptr) :: cptr_points_coords
  type(c_ptr) :: cptr_points_uvw
  type(c_ptr) :: cptr_points_weights_idx
  type(c_ptr) :: cptr_points_weights
  type(c_ptr) :: cptr_points_dist2
  type(c_ptr) :: cptr_points_projected_coords

  integer (c_int), pointer :: elt_pts_inside_idx (:)
  integer (kind = pdm_g_num_s), pointer :: points_gnum (:)
  double precision, pointer :: points_coords (:)
  double precision, pointer :: points_uvw (:)
  integer (c_int), pointer :: points_weights_idx (:)
  double precision, pointer :: points_weights (:)
  double precision, pointer :: points_dist2 (:)
  double precision, pointer :: points_projected_coords (:)

  !
  ! Init
  !

  call mpi_init(code)
  call mpi_comm_rank(mpi_comm_world, i_rank, code)
  call mpi_comm_size(mpi_comm_world, n_rank, code)

  if (n_rank .ne. 1) then
    print *,'Error : 1 MPI process is mandatory'
    call mpi_finalize(code)
    stop
  end if

  !
  ! Set point cloud : 2 points
  !

  allocate(coords_cloud(3*n_points_into_cloud))
  allocate(gnum_cloud(n_points_into_cloud))

  coords_cloud(1) = 0.5
  coords_cloud(2) = 0.5
  coords_cloud(3) = 0.5

  coords_cloud(4) = 1.5
  coords_cloud(5) = 0.5
  coords_cloud(6) = 0.5

  do i = 1, n_points_into_cloud
    gnum_cloud(i) = i
  end do

  cptr_coords_cloud = c_loc (coords_cloud)
  cptr_gnum_cloud = c_loc (gnum_cloud)

  !
  ! Set mesh : 2 hexa
  !

  allocate(coords_vtx(3*n_vtx))

  coords_vtx(1) = 0.d0
  coords_vtx(2) = 0.d0
  coords_vtx(3) = 0.d0

  coords_vtx(4) = 1.d0
  coords_vtx(5) = 0.d0
  coords_vtx(6) = 0.d0

  coords_vtx(7) = 1.d0
  coords_vtx(8) = 1.d0
  coords_vtx(9) = 0.d0

  coords_vtx(10) = 0.d0
  coords_vtx(11) = 1.d0
  coords_vtx(12) = 0.d0

  coords_vtx(13) = 0.d0
  coords_vtx(14) = 0.d0
  coords_vtx(15) = 1.d0

  coords_vtx(16) = 1.d0
  coords_vtx(17) = 0.d0
  coords_vtx(18) = 1.d0

  coords_vtx(19) = 1.d0
  coords_vtx(20) = 1.d0
  coords_vtx(21) = 1.d0

  coords_vtx(22) = 0.d0
  coords_vtx(23) = 1.d0
  coords_vtx(24) = 1.d0

  coords_vtx(25) = 2.d0
  coords_vtx(26) = 0.d0
  coords_vtx(27) = 0.d0

  coords_vtx(28) = 2.d0
  coords_vtx(29) = 1.d0
  coords_vtx(30) = 0.d0

  coords_vtx(31) = 2.d0
  coords_vtx(32) = 0.d0
  coords_vtx(33) = 1.d0

  coords_vtx(34) = 2.d0
  coords_vtx(35) = 1.d0
  coords_vtx(36) = 1.d0

  allocate(gnum_vtx(3*n_vtx))

  do i = 1, n_vtx
    gnum_vtx(i) = i
  end do

  cptr_coords_vtx = c_loc (coords_vtx)
  cptr_gnum_vtx = c_loc (gnum_vtx)


  allocate(face_vtx_idx(n_face+1))
  allocate(face_vtx(4*n_face))

  face_vtx_idx(1) = 0
  do i = 1, n_face
    face_vtx_idx(i+1) = face_vtx_idx(i) + 4
  end do

  face_vtx(1) = 1
  face_vtx(2) = 4
  face_vtx(3) = 3
  face_vtx(4) = 2

  face_vtx(5) = 5
  face_vtx(6) = 6
  face_vtx(7) = 7
  face_vtx(8) = 8

  face_vtx(9) = 1
  face_vtx(10) = 5
  face_vtx(11) = 8
  face_vtx(12) = 4

  face_vtx(13) = 6
  face_vtx(14) = 2
  face_vtx(15) = 3
  face_vtx(16) = 7

  face_vtx(17) = 3
  face_vtx(18) = 4
  face_vtx(19) = 8
  face_vtx(20) = 7

  face_vtx(21) = 1
  face_vtx(22) = 2
  face_vtx(23) = 6
  face_vtx(24) = 5

  face_vtx(25) = 11
  face_vtx(26) = 12
  face_vtx(27) = 7
  face_vtx(28) = 6

  face_vtx(29) = 9
  face_vtx(30) = 2
  face_vtx(31) = 3
  face_vtx(32) = 10

  face_vtx(33) = 12
  face_vtx(34) = 10
  face_vtx(35) = 3
  face_vtx(36) = 7

  face_vtx(37) = 11
  face_vtx(38) = 6
  face_vtx(39) = 2
  face_vtx(40) = 9

  face_vtx(41) = 9
  face_vtx(42) = 10
  face_vtx(43) = 12
  face_vtx(44) = 11

  allocate(gnum_face(n_face))

  do i = 1, n_face
    gnum_face(i) = i
  end do

  cptr_face_vtx_idx = c_loc (face_vtx_idx)
  cptr_face_vtx = c_loc (face_vtx)
  cptr_gnum_face = c_loc (gnum_face)


  allocate(cell_face_idx(n_cell+1))
  allocate(cell_face(6*n_cell))

  cell_face_idx(1) = 0
  do i = 1, n_cell
    cell_face_idx(i+1) = cell_face_idx(i) + 6
  end do

  cell_face(1) = 1
  cell_face(2) = 2
  cell_face(3) = 3
  cell_face(4) = 4
  cell_face(5) = 5
  cell_face(6) = 6

  cell_face(7) = 7
  cell_face(8) = 8
  cell_face(9) = 9
  cell_face(10) = 10
  cell_face(11) = 11
  cell_face(12) = -4

  allocate(gnum_cell(n_cell))

  do i = 1, n_cell
    gnum_cell(i) = i
  end do

  cptr_cell_face_idx = c_loc (cell_face_idx)
  cptr_cell_face = c_loc (cell_face)
  cptr_gnum_cell = c_loc (gnum_cell)


  !
  ! Create a new PDM_mesh_location structure
  !   The MPI communicator and the number of point cloud are setted
  !

  ml = PDM_mesh_location_create (PDM_MESH_NATURE_MESH_SETTED, &
                                 n_point_cloud, &
                                 MPI_COMM_WORLD)

  !
  ! Set the local number partition for any point cloud
  !

  call PDM_mesh_location_n_part_cloud_set (ml, &
                                           i_point_cloud, &
                                           n_part_cloud)

  !
  ! Set point properties for any partition of point cloud
  !

  call PDM_mesh_location_cloud_set (ml, &
                                    i_point_cloud, &
                                    i_part_cloud, &
                                    n_points_into_cloud, &
                                    cptr_coords_cloud, &
                                    cptr_gnum_cloud)

  !
  ! Set mesh
  !

  call PDM_mesh_location_mesh_global_data_set (ml, &
                                               n_part_mesh)

  call PDM_mesh_location_part_set (ml, &
                                   i_part_mesh, &
                                   n_cell, &
                                   cptr_cell_face_idx, &
                                   cptr_cell_face, &
                                   cptr_gnum_cell, &
                                   n_face, &
                                   cptr_face_vtx_idx, &
                                   cptr_face_vtx, &
                                   cptr_gnum_face, &
                                   n_vtx, &
                                   cptr_coords_vtx, &
                                   cptr_gnum_vtx)


  !
  ! Set options
  !

  call PDM_mesh_location_tolerance_set (ml, 1.d-4)

  call PDM_mesh_location_method_set (ml, PDM_MESH_LOCATION_OCTREE) ! or PDM_MESH_LOCATION_DBBTREE

  !
  ! Compute
  !

  call PDM_mesh_location_compute (ml)

  !
  ! Get results
  !

  n_unlocated = PDM_mesh_location_n_unlocated_get (ml, &
                                                   i_point_cloud, &
                                                   i_part_cloud)

  n_located = PDM_mesh_location_n_located_get (ml, &
                                               i_point_cloud, &
                                               i_part_cloud)

  cptr_unlocated = PDM_mesh_location_unlocated_get (ml, &
                                                    i_point_cloud, &
                                                    i_part_cloud)

  cptr_located = PDM_mesh_location_located_get (ml, &
                                                i_point_cloud, &
                                                i_part_cloud)
  call c_f_pointer(cptr_located, located, [n_located])
  call c_f_pointer(cptr_unlocated, unlocated, [n_unlocated])

  call PDM_mesh_location_point_location_get (ml, &
                                             i_point_cloud, &
                                             i_part_cloud, &
                                             cptr_location, &
                                             cptr_dist2, &
                                             cptr_projected_coords)

  call c_f_pointer(cptr_location, location, [n_located])
  call c_f_pointer(cptr_dist2, dist2, [n_located])
  call c_f_pointer(cptr_projected_coords, projected_coords, [3*n_located])

  call PDM_mesh_location_points_in_elt_get (ml, &
                                            i_part_mesh, &
                                            i_point_cloud, &
                                            cptr_elt_pts_inside_idx, &
                                            cptr_points_gnum, &
                                            cptr_points_coords, &
                                            cptr_points_uvw, &
                                            cptr_points_weights_idx, &
                                            cptr_points_weights, &
                                            cptr_points_dist2, &
                                            cptr_points_projected_coords)

  call c_f_pointer(cptr_elt_pts_inside_idx, elt_pts_inside_idx, [n_cell + 1])
  call c_f_pointer(cptr_points_gnum, points_gnum, [elt_pts_inside_idx(n_cell+1)])
  call c_f_pointer(cptr_points_coords, points_coords, [3 * elt_pts_inside_idx(n_cell+1)])
  call c_f_pointer(cptr_points_uvw, points_uvw, [3 * elt_pts_inside_idx(n_cell+1)])
  call c_f_pointer(cptr_points_weights_idx, points_weights_idx, [elt_pts_inside_idx(n_cell+1)])
  call c_f_pointer(cptr_points_weights, points_weights, [points_weights_idx(elt_pts_inside_idx(n_cell + 1))])
  call c_f_pointer(cptr_points_dist2, points_dist2, [elt_pts_inside_idx(n_cell+1)])
  call c_f_pointer(cptr_points_projected_coords, points_projected_coords, [3 * elt_pts_inside_idx(n_cell+1)])

  call PDM_mesh_location_dump_times (ml)

  call PDM_mesh_location_free (ml, 0)

  deallocate(coords_cloud)
  deallocate(gnum_cloud)

  deallocate(cell_face_idx)
  deallocate(cell_face)
  deallocate(gnum_cell)

  deallocate(face_vtx_idx)
  deallocate(face_vtx)
  deallocate(gnum_face)

  deallocate(gnum_vtx)
  deallocate(coords_vtx)

  call mpi_finalize(code)


end program testf
