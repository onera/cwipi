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

#include "pdm_configf.h"

module pdm_mesh_location

  use pdm
  use iso_c_binding

  implicit none

  !!
  !! Enum type PDM_mesh_location_method_t
  !!

  integer(c_int), parameter :: PDM_MESH_LOCATION_OCTREE = 0
  integer(c_int), parameter :: PDM_MESH_LOCATION_DBBTREE = 1


  interface

    !>
    !!
    !! \brief Create a structure to compute the location of point clouds inta a mesh
    !!
    !! \param [in]   mesh_nature    Nature of the mesh
    !! \param [in]   n_point_cloud  Number of point cloud
    !! \param [in]   comm           MPI communicator
    !!
    !! \return     Identifier
    !!
    !!

    subroutine PDM_mesh_location_create (mesh_nature, &
                                         n_point_cloud, &
                                         fcomm, &
                                         id) &
      bind (c, name = 'PDM_mesh_location_create_cf')

      use iso_c_binding

      implicit none

      integer(c_int), value :: mesh_nature
      integer(c_int), value :: n_point_cloud
      integer(c_int), value :: fComm

      integer(c_int)        :: id

    end subroutine PDM_mesh_location_create

    !>
    !!
    !! \brief Set the number of partitions of a point cloud
    !!
    !! \param [in]   id              Identifier
    !! \param [in]   i_point_cloud   Index of point cloud
    !! \param [in]   n_part          Number of partitions
    !!
    !!

    subroutine PDM_mesh_location_n_part_cloud_set (id, &
                                                   i_point_cloud, &
                                                   n_part) &
     bind (c, name = 'PDM_mesh_location_n_part_cloud_set')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id
      integer(c_int), value :: i_point_cloud
      integer(c_int), value :: n_part

    end subroutine PDM_mesh_location_n_part_cloud_set

    !>
    !!
    !! \brief Set a point cloud
    !!
    !! \param [in]   id              Identifier
    !! \param [in]   i_point_cloud   Index of point cloud
    !! \param [in]   i_part          Index of partition
    !! \param [in]   n_points        Number of points
    !! \param [in]   coords          Point coordinates
    !! \param [in]   gnum            Point global number
    !!
    !!

    subroutine PDM_mesh_location_cloud_set (id, &
                                            i_point_cloud, &
                                            i_part, &
                                            n_points, &
                                            coords, &
                                            gnum) &
     bind (c, name = 'PDM_mesh_location_cloud_set')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id
      integer(c_int), value :: i_point_cloud
      integer(c_int), value :: i_part
      integer(c_int), value :: n_points
      type(c_ptr), value    :: coords
      type(c_ptr), value    :: gnum

    end subroutine PDM_mesh_location_cloud_set

        !>
    !!
    !! \brief Set a point cloud
    !!
    !! \param [in]   id              Identifier
    !! \param [in]   i_point_cloud   Index of point cloud
    !! \param [in]   i_part          Index of partition
    !! \param [in]   n_points        Number of points
    !! \param [in]   coords          Point coordinates
    !! \param [in]   gnum            Point global number
    !!
    !!

    subroutine PDM_mesh_location_cloud_get (id, &
                                            i_point_cloud, &
                                            i_part, &
                                            n_points, &
                                            coords, &
                                            gnum) &
     bind (c, name = 'PDM_mesh_location_cloud_get')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id
      integer(c_int), value :: i_point_cloud
      integer(c_int), value :: i_part
      integer(c_int), value :: n_points
      type(c_ptr)           :: coords
      type(c_ptr)           :: gnum

    end subroutine PDM_mesh_location_cloud_get

    !>
    !!
    !! \brief Set the mesh nodal
    !!
    !! \param [in]   id             Identifier
    !! \param [in]   mesh_nodal_id  Mesh nodal identifier
    !!
    !!

    subroutine PDM_mesh_location_shared_nodal_mesh_set (id, &
                                                        mesh_nodal_id) &
     bind (c, name = 'PDM_mesh_location_shared_nodal_mesh_set')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id
      integer(c_int), value :: mesh_nodal_id

    end subroutine PDM_mesh_location_shared_nodal_mesh_set

    !>
    !!
    !! \brief Set global data of a mesh
    !!
    !! \param [in]   id             Identifier
    !! \param [in]   n_part         Number of partition
    !!
    !!

    subroutine PDM_mesh_location_mesh_global_data_set (id, &
                                                       n_part) &
     bind (c, name = 'PDM_mesh_location_mesh_global_data_set')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id
      integer(c_int), value :: n_part

    end subroutine PDM_mesh_location_mesh_global_data_set

    !>
    !!
    !! \brief Set a part of a mesh
    !!
    !! \param [in]   id            Identifier
    !! \param [in]   i_part        Partition to define
    !! \param [in]   n_cell        Number of cells
    !! \param [in]   cell_face_idx Index in the cell -> face connectivity
    !! \param [in]   cell_face     cell -> face connectivity
    !! \param [in]   cell_ln_to_gn Local cell numbering to global cel numbering
    !! \param [in]   n_face        Number of faces
    !! \param [in]   face_vtx_idx  Index in the face -> vertex connectivity
    !! \param [in]   face_vtx      face -> vertex connectivity
    !! \param [in]   face_ln_to_gn Local face numbering to global face numbering
    !! \param [in]   n_vtx         Number of vertices
    !! \param [in]   coords        Coordinates
    !! \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
    !!
    !!

    subroutine PDM_mesh_location_part_set (id, &
                                           i_part, &
                                           n_cell, &
                                           cell_face_idx, &
                                           cell_face, &
                                           cell_ln_to_gn, &
                                           n_face, &
                                           face_vtx_idx, &
                                           face_vtx, &
                                           face_ln_to_gn, &
                                           n_vtx, &
                                           coords, &
                                           vtx_ln_to_gn) &
     bind (c, name = 'PDM_mesh_location_part_set')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id
      integer(c_int), value :: i_part
      integer(c_int), value :: n_cell
      type(c_ptr), value    :: cell_face_idx
      type(c_ptr), value    :: cell_face
      type(c_ptr), value    :: cell_ln_to_gn
      integer(c_int), value :: n_face
      type(c_ptr), value    :: face_vtx_idx
      type(c_ptr), value    :: face_vtx
      type(c_ptr), value    :: face_ln_to_gn
      integer(c_int), value :: n_vtx
      type(c_ptr), value    :: coords
      type(c_ptr), value    :: vtx_ln_to_gn

    end subroutine PDM_mesh_location_part_set

    !>
    !!
    !! \brief Set a part of a mesh (2d version)
    !!
    !! \param [in]   id            Identifier
    !! \param [in]   i_part        Partition to define
    !! \param [in]   n_cell        Number of cells
    !! \param [in]   cell_edge_idx Index in the cell -> edge connectivity
    !! \param [in]   cell_edge     cell -> edge connectivity
    !! \param [in]   cell_ln_to_gn Local cell numbering to global cel numbering
    !! \param [in]   n_edge        Number of edges
    !! \param [in]   edge_vtx_idx  Index in the edge -> vertex connectivity
    !! \param [in]   edge_vtx      edge -> vertex connectivity
    !! \param [in]   edge_ln_to_gn Local edge numbering to global edge numbering
    !! \param [in]   n_vtx         Number of vertices
    !! \param [in]   coords        Coordinates
    !! \param [in]   vtx_ln_to_gn  Local vertex numbering to global vertex numbering
    !!
    !!

    subroutine PDM_mesh_location_part_set_2d (id, &
                                              i_part, &
                                              n_cell, &
                                              cell_edge_idx, &
                                              cell_edge, &
                                              cell_ln_to_gn, &
                                              n_edge, &
                                              edge_vtx_idx, &
                                              edge_vtx, &
                                              edge_ln_to_gn, &
                                              n_vtx, &
                                              coords, &
                                              vtx_ln_to_gn) &
     bind (c, name = 'PDM_mesh_location_cloud_set_2d')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id
      integer(c_int), value :: i_part
      integer(c_int), value :: n_cell
      type(c_ptr), value    :: cell_edge_idx
      type(c_ptr), value    :: cell_edge
      type(c_ptr), value    :: cell_ln_to_gn
      integer(c_int), value :: n_edge
      type(c_ptr), value    :: edge_vtx_idx
      type(c_ptr), value    :: edge_vtx
      type(c_ptr), value    :: edge_ln_to_gn
      integer(c_int), value :: n_vtx
      type(c_ptr), value    :: coords
      type(c_ptr), value    :: vtx_ln_to_gn

    end subroutine PDM_mesh_location_part_set_2d

    !>
    !!
    !! \brief Set the tolerance for bounding boxes
    !!
    !! \param [in]   id              Identifier
    !! \param [in]   tol             Tolerance
    !!
    !!

    subroutine PDM_mesh_location_tolerance_set (id, &
                                                tol) &
     bind (c, name = 'PDM_mesh_location_tolerance_set')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id
      real(c_double), value :: tol

    end subroutine PDM_mesh_location_tolerance_set

    !>
    !!
    !! \brief Set the method for computing location
    !!
    !! \param [in]   id              Identifier
    !! \param [in]   method          Method
    !!
    !!

    subroutine PDM_mesh_location_method_set (id, &
                                             method) &
     bind (c, name = 'PDM_mesh_location_method_set')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id
      integer(c_int), value :: method

    end subroutine PDM_mesh_location_method_set

    !>
    !!
    !! \brief Compute point location
    !!
    !! \param [in]   id  Identifier
    !!
    !!

    subroutine PDM_mesh_location_compute (id) &
      bind (c, name = 'PDM_mesh_location_compute')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id

    end subroutine PDM_mesh_location_compute

    !>
    !!
    !! \brief Get mesh location
    !!
    !! \param [in]   id                    Identifier
    !! \param [in]   i_point_cloud         Current cloud
    !! \param [in]   i_part                Index of partition of the cloud
    !! \param [out]  n_points              Number of points in point cloud
    !! \param [out]  coord                 Coordinates of points in point cloud
    !! \param [out]  g_num                 Global numbers of points in point cloud
    !! \param [out]  location              The global number of the closest element if the point is located,
    !!                                     -1 otherwise
    !!
    !!

    subroutine PDM_mesh_location_get (id, &
                                      i_point_cloud, &
                                      i_part, &
                                      location, &
                                      weights_idx, &
                                      weights, &
                                      projected_coords) &
     bind (c, name = 'PDM_mesh_location_get')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id
      integer(c_int), value :: i_point_cloud
      integer(c_int), value :: i_part
      type(c_ptr)           :: location
      type(c_ptr)           :: weights_idx
      type(c_ptr)           :: weights
      type(c_ptr)           :: projected_coords


    end subroutine PDM_mesh_location_get

    !>
    !!
    !! \brief Free a locationd mesh structure
    !!
    !! \param [in]  id       Identifier
    !! \param [in]  partial  if partial is equal to 0, all data are removed.
    !!                       Otherwise, results are kept.
    !!
    !!

    subroutine PDM_mesh_location_free (id, &
                                       partial) &
     bind (c, name = 'PDM_mesh_location_free')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id
      integer(c_int), value :: partial

    end subroutine PDM_mesh_location_free


    !>
    !!
    !! \brief  Dump elapsed an CPU time
    !!
    !! \param [in]  id       Identifier
    !!
    !!

    subroutine PDM_mesh_location_dump_times (id) &
      bind (c, name = 'PDM_mesh_location_dump_times')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id

    end subroutine PDM_mesh_location_dump_times


    function PDM_mesh_location_mesh_nodal_id_get (id) &
                                                 result(mesh_nodal_id) &
      bind (c, name = 'PDM_mesh_location_mesh_nodal_id_get')

      use iso_c_binding

      implicit none

      integer(c_int), value :: id
      integer(c_int)        :: mesh_nodal_id

    end function PDM_mesh_location_mesh_nodal_id_get

  end interface

end module pdm_mesh_location
