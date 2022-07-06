!-----------------------------------------------------------------------------
! This file is part of the CWIPI library.
!
! Copyright (C) 2011  ONERA
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

module cwp
    use iso_c_binding

    ! CWP_Status_t
    enum, bind(c)
        enumerator :: &
                CWP_STATUS_OFF, &
                CWP_STATUS_ON
    end enum

    ! CWP_Type_t
    enum, bind(c)
        enumerator :: &
                CWP_DOUBLE, &
                CWP_INT, &
                CWP_CHAR
    end enum

    ! CWP_Visu_format_t
    enum, bind(c)
        enumerator :: &
                CWP_VISU_FORMAT_ENSIGHT
    end enum

    ! CWP_Comm_t
    enum, bind(c)
        enumerator :: &
                CWP_COMM_PAR_WITH_PART, &
                CWP_COMM_PAR_WITHOUT_PART, &
                CWP_COMM_SEQ, &
                CWP_COMM_INTERNAL
    end enum

    ! CWP_Time_exch_t
    enum, bind(c)
        enumerator :: &
            CWP_TIME_EXCH_EACH_TIME_STEP, &
            CWP_TIME_EXCH_N_TIME_STEP, &
            CWP_TIME_EXCH_CPL_TIME_STEP, &
            CWP_TIME_EXCH_ASYNCHRONOUS, &
            CWP_TIME_EXCH_SLAVE, &
            CWP_TIME_EXCH_MASTER
    end enum

    ! CWP_Dof_location_t
    enum, bind(c)
        enumerator :: &
                CWP_DOF_LOCATION_UNDEF, &
                CWP_DOF_LOCATION_CELL_CENTER, &
                CWP_DOF_LOCATION_NODE, &
                CWP_DOF_LOCATION_USER
    end enum

    ! CWP_Field_exch_t
    enum, bind(c)
        enumerator :: &
                CWP_FIELD_EXCH_SEND, &
                CWP_FIELD_EXCH_RECV, &
                CWP_FIELD_EXCH_SENDRECV
    end enum

    ! CWP_Field_exch_t
    enum, bind(c)
      enumerator :: &
                CWP_FIELD_MAP_SOURCE, &
                CWP_FIELD_MAP_TARGET
    end enum

    ! CWP_Field_storage_t
    enum, bind(c)
        enumerator :: &
                CWP_FIELD_STORAGE_INTERLACED, &
                CWP_FIELD_STORAGE_BLOCK
    end enum

    ! CWP_Block_t
    enum, bind(c)
        enumerator :: &
                CWP_BLOCK_NODE, &
                CWP_BLOCK_EDGE2, &
                CWP_BLOCK_FACE_TRIA3, &
                CWP_BLOCK_FACE_QUAD4, &
                CWP_BLOCK_FACE_POLY, &
                CWP_BLOCK_CELL_TETRA4, &
                CWP_BLOCK_CELL_HEXA8, &
                CWP_BLOCK_CELL_PRISM6, &
                CWP_BLOCK_CELL_PYRAM5, &
                CWP_BLOCK_CELL_POLY
    end enum

    ! CWP_Spatial_interp_t
    enum, bind(c)
        enumerator :: &
                CWP_SPATIAL_INTERP_FROM_CLOSEST_POINT_LEAST_SQUARES, &
                CWP_SPATIAL_INTERP_FROM_INTERSECTION, &
                CWP_SPATIAL_INTERP_FROM_LOCATION_DIST_CLOUD_SURF, &
                CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, &
                CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_DBBTREE
    end enum

    ! CWP_Interface_t
    enum, bind(c)
        enumerator :: &
                CWP_INTERFACE_POINT, &
                CWP_INTERFACE_LINEAR, &
                CWP_INTERFACE_SURFACE, &
                CWP_INTERFACE_VOLUME
    end enum

    ! CWP_Dynamic_mesh_t
    enum, bind(c)
        enumerator :: &
                CWP_DYNAMIC_MESH_STATIC, &
                CWP_DYNAMIC_MESH_DEFORMABLE, &
                CWP_DYNAMIC_MESH_VARIABLE
    end enum

    interface
      subroutine CWP_Init_cf(fcomm, n_code, code_names, l_code_names, is_active_rank, time_init, intra_comms) &
              bind(c, name = 'CWP_Init_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int), value :: fcomm
        integer(c_int), value :: n_code
        type(c_ptr) :: code_names
        integer(c_int), dimension(n_code) :: l_code_names
        integer(c_int), dimension(n_code) :: is_active_rank
        real(c_double), dimension(n_code) :: time_init
        type(c_ptr), value :: intra_comms
      end subroutine CWP_Init_cf

      subroutine CWP_Cpl_create_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, coupled_code_name, l_coupled_code_name, &
        entities_dim, comm_type, spatial_interp, n_part, displacement, freq) &
        bind(c, name = 'CWP_Cpl_create_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, coupled_code_name
        integer(c_int), value :: l_local_code_name, l_cpl_id, l_coupled_code_name
        integer(kind = c_int), value :: entities_dim
        integer(kind = c_int), value :: comm_type
        integer(kind = c_int), value :: spatial_interp
        integer(kind = c_int), value :: n_part
        integer(kind = c_int), value :: displacement
        integer(kind = c_int), value :: freq
      end subroutine CWP_Cpl_Create_cf

      subroutine CWP_Cpl_del_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id) &
        bind(c, name = 'CWP_Cpl_del_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Cpl_del_cf

      function CWP_N_uncomputed_tgts_get_cf (local_code_name,   &
                                             l_local_code_name, &
                                             cpl_id,            &
                                             l_cpl_id,          &
                                             field_id,          &
                                             l_field_id,        &
                                             i_part)            &
                                             result (n_uncomputed_tgts) &
        bind(c, name = 'CWP_N_uncomputed_tgts_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        character(kind = c_char, len = 1) :: field_id
        integer(kind = c_int), value :: l_field_id
        integer(c_int), value :: i_part
        integer(kind = c_int) :: n_uncomputed_tgts
      end function CWP_N_uncomputed_tgts_get_cf

      function CWP_N_computed_tgts_get_cf(local_code_name,   &
                                          l_local_code_name, &
                                          cpl_id,            &
                                          l_cpl_id,          &
                                          field_id,          &
                                          l_field_id,        &
                                          i_part)            &
                                          result (n_computed_tgts) &
            bind(c, name = 'CWP_N_computed_tgts_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        character(kind = c_char, len = 1) :: field_id
        integer(kind = c_int), value :: l_field_id
        integer(c_int), value :: i_part
        integer(kind = c_int) :: n_computed_tgts
      end function CWP_N_computed_tgts_get_cf

      function CWP_N_involved_srcs_get_cf(local_code_name,   &
              l_local_code_name, &
              cpl_id,            &
              l_cpl_id,          &
              field_id,          &
              l_field_id,        &
              i_part)            &
              result (n_involved_srcs) &
                      bind(c, name = 'CWP_N_involved_srcs_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        character(kind = c_char, len = 1) :: field_id
        integer(kind = c_int), value :: l_field_id
        integer(c_int), value :: i_part
        integer(kind = c_int) :: n_involved_srcs
      end function CWP_N_involved_srcs_get_cf

      function CWP_Computed_tgts_get_cf(local_code_name,   &
                                        l_local_code_name, &
                                        cpl_id,            &
                                        l_cpl_id,          &
                                        field_id,          &
                                        l_field_id,        &
                                        i_part)            &
                                        result (computed_tgts) &
            bind(c, name = 'CWP_Computed_tgts_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        character(kind = c_char, len = 1) :: field_id
        integer(kind = c_int), value :: l_field_id
        integer(c_int), value :: i_part
        type(c_ptr) :: computed_tgts
      end function CWP_Computed_tgts_get_cf


      function CWP_Uncomputed_tgts_get_cf (local_code_name,   &
                                          l_local_code_name, &
                                          cpl_id,            &
                                          l_cpl_id,          &
                                          field_id,          &
                                          l_field_id,        &
                                          i_part)            &
                                          result (uncomputed_tgts) &
        bind(c, name = 'CWP_Uncomputed_tgts_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        character(kind = c_char, len = 1) :: field_id
        integer(kind = c_int), value :: l_field_id
        integer(c_int), value :: i_part
        type(c_ptr) :: uncomputed_tgts
      end function CWP_Uncomputed_tgts_get_cf

      function CWP_Involved_srcs_get_cf(local_code_name,   &
              l_local_code_name, &
              cpl_id,            &
              l_cpl_id,          &
              field_id,          &
              l_field_id,        &
              i_part)            &
              result (involved_srcs) &
                      bind(c, name = 'CWP_Involved_srcs_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        character(kind = c_char, len = 1) :: field_id
        integer(kind = c_int), value :: l_field_id
        integer(c_int), value :: i_part
        type(c_ptr) :: involved_srcs
      end function CWP_Involved_srcs_get_cf

      function CWP_Computed_tgts_dist_to_spatial_interp_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id) &
            result (dists) &
        bind(c, name = 'CWP_Computed_tgts_dist_to_spatial_interp_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        type(c_ptr) :: dists
      end function CWP_Computed_tgts_dist_to_spatial_interp_get_cf

      subroutine CWP_Spatial_interp_weights_compute_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id) &
            bind(c, name = 'CWP_Spatial_interp_weights_compute_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Spatial_interp_weights_compute_cf

      subroutine CWP_Visu_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, freq, &
            format, format_option, l_format_option) &
            bind(c, name = 'CWP_Visu_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, format_option
        integer(c_int), value :: freq
        integer(c_int), value :: format
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_format_option
      end subroutine CWP_Visu_set_cf

      subroutine CWP_User_tgt_pts_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, n_pts, coord, global_num) &
              bind(c, name = 'CWP_User_tgt_pts_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        integer(c_int), value :: i_part, n_pts
        type(c_ptr), value :: coord
        type(c_ptr), value :: global_num
      end subroutine CWP_User_tgt_pts_set_cf

      subroutine CWP_Mesh_interf_finalize_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id) &
            bind(c, name = 'CWP_Mesh_interf_finalize_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_finalize_cf

      subroutine CWP_Mesh_interf_vtx_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
            i_part, n_pts, coord, global_num) &
            bind(c, name = 'CWP_Mesh_interf_vtx_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: i_part, n_pts
        type(c_ptr), value :: coord
        type(c_ptr), value :: global_num
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_vtx_set_cf

      function CWP_Mesh_interf_block_add_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, block_type) result(block_id) &
        bind(c, name = 'CWP_Mesh_interf_block_add_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: block_type
        integer(c_int) :: block_id
        integer(c_int), value :: l_local_code_name, l_cpl_id
      end function CWP_Mesh_interf_block_add_cf

      subroutine CWP_Mesh_interf_block_std_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
            i_part, block_id, n_elts, connec, global_num) &
        bind(c, name = 'CWP_Mesh_interf_block_std_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: i_part, block_id, n_elts
        type(c_ptr), value :: connec, global_num
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_block_std_set_cf

      subroutine CWP_Mesh_interf_f_poly_block_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
            i_part, block_id, n_elts, connec_idx, connec, global_num) &
        bind(c, name = 'CWP_Mesh_interf_f_poly_block_set')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: i_part, block_id, n_elts
        type(c_ptr), value :: connec_idx, connec, global_num
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_f_poly_block_set_cf

      subroutine CWP_Mesh_interf_c_poly_block_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
            i_part, block_id, n_elts, n_faces, connec_faces_idx, connec_faces, connec_cells_idx, connec_cells, global_num) &
        bind(c, name = 'CWP_Mesh_interf_c_poly_block_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: i_part, block_id, n_elts, n_faces
        type(c_ptr), value :: connec_faces_idx, connec_faces, connec_cells_idx, connec_cells, global_num
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_c_poly_block_set_cf

      subroutine CWP_Mesh_interf_del_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id) &
              bind(c, name = 'CWP_Mesh_interf_del_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_del_cf

      subroutine CWP_Mesh_interf_from_cellface_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, n_cells, &
            cell_face_idx, cell_face, n_faces, face_vtx_idx, face_vtx, parent_num) &
            bind(c, name = 'CWP_Mesh_interf_from_cellface_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: i_part, n_cells, n_faces
        type(c_ptr), value :: cell_face_idx, cell_face, face_vtx_idx, face_vtx, parent_num
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_from_cellface_set_cf

      subroutine CWP_Mesh_interf_from_faceedge_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, n_faces, &
            face_edge_idx, face_edge, n_edges, edge_vtx_idx, edge_vtx, parent_num) &
            bind(c, name = 'CWP_Mesh_interf_from_faceedge_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: i_part, n_faces, n_edges
        type(c_ptr), value :: face_edge_idx, face_edge, edge_vtx_idx, edge_vtx, parent_num
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Mesh_interf_from_faceedge_set_cf

      subroutine CWP_Field_create_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id, &
          data_type, storage, n_component, target_location, exch_type, visu_status) &
          bind(c, name = 'CWP_Field_create_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
        integer(c_int), value :: data_type, storage, n_component, target_location, exch_type, visu_status
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      end subroutine

      subroutine CWP_Field_data_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id, &
            i_part, map_type, data) &
            bind(c, name = 'CWP_Field_data_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
        integer(c_int), value :: i_part
        integer(c_int), value :: map_type
        type(c_ptr), value :: data
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      end subroutine CWP_Field_data_set_cf

      subroutine CWP_Field_issend_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, src_field_id, l_src_field_id) &
            bind(c, name = 'CWP_Field_issend_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, src_field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_src_field_id
      end subroutine CWP_Field_issend_cf

      subroutine CWP_Field_irecv_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, tgt_field_id, l_tgt_field_id) &
            bind(c, name = 'CWP_Field_irecv_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, tgt_field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_tgt_field_id
      end subroutine CWP_Field_irecv_cf

      subroutine CWP_Field_wait_issend_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, src_field_id, l_src_field_id) &
            bind(c, name = 'CWP_Field_wait_issend_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, src_field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_src_field_id
      end subroutine CWP_Field_wait_issend_cf

      subroutine CWP_Field_wait_irecv_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, tgt_field_id, l_tgt_field_id) &
            bind(c, name = 'CWP_Field_wait_irecv_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, tgt_field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_tgt_field_id
      end subroutine CWP_Field_wait_irecv_cf

      subroutine CWP_Interp_from_location_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
              src_field_id, l_src_field_id, ptInterpolationFct) &
              bind(c, name = 'CWP_Interp_from_location_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none

        abstract interface
          subroutine ptInterpolationFct ( &
                  interface_type, &
                  n_src_vtcs, &
                  n_src_std_elts, &
                  n_tgt_pts, &
                  src_vtcs_coords, &
                  src_connec_idx, &
                  src_connec, &
                  tgt_pts_coords, &
                  tgt_pts_target_location, &
                  tgt_pts_dist, &
                  tgt_pts_bary_coords_idx, &
                  tgt_pts_bary_coords, &
                  stride, &
                  src_field_location, &
                  src_field, &
                  tgt_field_location, &
                  tgt_field &
                  ) &
                  bind (c)
            use, intrinsic :: iso_c_binding
            integer(kind = c_int)               :: interface_type
            integer(kind = c_int)               :: n_src_vtcs
            integer(kind = c_int)               :: n_src_std_elts
            integer(kind = c_int)               :: n_tgt_pts
            real(kind = c_double), dimension(*) :: src_vtcs_coords
            integer(kind = c_int), dimension(*) :: src_connec_idx
            integer(kind = c_int), dimension(*) :: src_connec
            real(kind = c_double), dimension(*) :: tgt_pts_coords
            integer(kind = c_int), dimension(*) :: tgt_pts_target_location
            real(kind = c_double), dimension(*) :: tgt_pts_dist
            integer(kind = c_int), dimension(*) :: tgt_pts_bary_coords_idx
            real(kind = c_double), dimension(*) :: tgt_pts_bary_coords
            integer(kind = c_int)               :: stride
            integer(kind = c_int)               :: src_field_location
            real(kind = c_double), dimension(*) :: src_field
            integer(kind = c_int)               :: tgt_field_location
            real(kind = c_double), dimension(*) :: tgt_field
          end subroutine ptInterpolationFct
        end interface
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, src_field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_src_field_id
      end subroutine CWP_Interp_from_location_set_cf

      subroutine CWP_Finalize() &
            bind(c, name = 'CWP_Finalize')
      end subroutine CWP_Finalize
  end interface

contains


  !>
  !! \brief Initialize CWIPI.
  !!
  !! This function creates the MPI intra communicators of the codes from
  !! the \p global_comm MPI communicator that contains all code ranks. This
  !! function has to be called from all ranks contained in the \p global_comm.
  !!
  !! \param [in]  global_comm    MPI global communicator
  !! \param [in]  n_code         Number of codes on the current rank
  !! \param [in]  code_names     Names of codes on the current rank (size = \p n_code)
  !! \param [in]  is_active_rank Is current rank have to be used by CWIPI (size = \p n_code)
  !! \param [in]  time_init      Initial time (size = \p n_code)
  !! \param [out] intra_comms    MPI intra communicators of each code (size = \p n_code)
  !!
  !!

  subroutine CWP_Init (fcomm,          &
                       n_code,         &
                       code_names,     &
                       is_active_rank, &
                       time_init,      &
                       intra_comms)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int) :: fcomm
    integer(c_int) :: n_code
    character(kind = c_char, len = *), dimension(:), pointer :: code_names
    integer(c_int), dimension(n_code) :: is_active_rank
    real(c_double), dimension(n_code) :: time_init
    integer(c_int), dimension(:), pointer :: intra_comms
    integer, dimension(n_code) :: l_code_names
    integer :: i

    do i = 1, n_code
      l_code_names(i) = len(code_names(i))
    end do

    call CWP_Init_cf(fcomm, n_code,     &
                     c_loc(code_names), &
                    l_code_names,       &
                    is_active_rank,     &
                    time_init,          &
                    c_loc(intra_comms))
  end subroutine CWP_Init


  !>
  !! \brief Create a coupling object and define its properties.
  !!
  !! \param [in]  local_code_name     Local code name
  !! \param [in]  cpl_id              Coupling identifier
  !! \param [in]  coupled_code_name   Distant or local coupled code name
  !! \param [in]  comm_type           Communication type
  !! \param [in]  spatial_interp      Spatial interpolation method
  !! \param [in]  n_part              Number of interface partition
  !! \param [in]  displacement        Mesh moving status
  !! \param [in]  recv_freq_type      Type of receiving frequency
  !!
  !!

  subroutine CWP_Cpl_create (local_code_name,   &
                             cpl_id,            &
                             coupled_code_name, &
                             entities_dim,      &
                             comm_type,         &
                             spatial_interp,    &
                             n_part,            &
                             displacement,      &
                             freq)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, coupled_code_name
    integer(kind = c_int) :: entities_dim
    integer(kind = c_int) :: comm_type
    integer(kind = c_int) :: spatial_interp
    integer(kind = c_int) :: n_part
    integer(kind = c_int) :: displacement
    integer(kind = c_int) :: freq
    integer(c_int) :: l_local_code_name, l_cpl_id, l_coupled_code_name

    l_local_code_name   = len(local_code_name)
    l_cpl_id            = len(cpl_id)
    l_coupled_code_name = len(coupled_code_name)

    call CWP_Cpl_create_cf (local_code_name,    &
                            l_local_code_name,  &
                            cpl_id,             &
                            l_cpl_id,           &
                            coupled_code_name,  &
                            l_coupled_code_name,&
                            entities_dim,       &
                            comm_type,          &
                            spatial_interp,     &
                            n_part,             &
                            displacement,       &
                            freq)
  end subroutine CWP_Cpl_Create


  !>
  !!
  !! \brief Delete a coupling object.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !!
  !!

  subroutine CWP_Cpl_Del (local_code_name, &
                          cpl_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)

    call CWP_Cpl_del_cf (local_code_name, l_local_code_name, cpl_id, l_cpl_id)
  end subroutine CWP_Cpl_Del


  !>
  !!
  !! \brief Return the number of uncomputed targets.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] field_id         Field identifier
  !! \param [in] i_part           Current partition
  !!
  !! \return                Number of uncomputed targets
  !!

  function CWP_N_uncomputed_tgts_get (local_code_name, &
                                      cpl_id,          &
                                      field_id,        &
                                      i_part)          &
                                      result (n_uncomputed_tgts)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: cpl_id
    character(kind = c_char, len = *) :: field_id
    integer(c_int) :: i_part

    integer(c_int) :: n_uncomputed_tgts

    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    n_uncomputed_tgts = CWP_N_uncomputed_tgts_get_cf (local_code_name,   &
                                                      l_local_code_name, &
                                                      cpl_id,            &
                                                      l_cpl_id,          &
                                                      field_id,          &
                                                      l_field_id,        &
                                                      i_part)
  end function CWP_N_uncomputed_tgts_get


  !>
  !!
  !! \brief Return uncomputed targets.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] field_id         Field identifier
  !! \param [in] i_part           Current partition
  !!
  !! \return                Uncomputed targets
  !!

  function CWP_Uncomputed_tgts_get (local_code_name, &
                                    cpl_id,          &
                                    field_id,        &
                                    i_part)          &
                                    result (uncomputed_tgts)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: cpl_id
    character(kind = c_char, len = *) :: field_id
    integer(c_int) :: i_part

    integer(c_int), dimension(:), pointer :: uncomputed_tgts

    type(c_ptr) :: cptr_uncomputed_tgts
    integer(c_int) :: n_uncomputed_tgts
    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    n_uncomputed_tgts = CWP_N_uncomputed_tgts_get_cf (local_code_name,   &
                                                      l_local_code_name, &
                                                      cpl_id,            &
                                                      l_cpl_id,          &
                                                      field_id,          &
                                                      l_field_id,        &
                                                      i_part)

    cptr_uncomputed_tgts = CWP_Uncomputed_tgts_get_cf(local_code_name,   &
                                                      l_local_code_name, &
                                                      cpl_id,            &
                                                      l_cpl_id,          &
                                                      field_id,          &
                                                      l_field_id,        &
                                                      i_part)


    call c_f_pointer (cptr = cptr_uncomputed_tgts, &
                      fptr = uncomputed_tgts  ,    &
                      shape= [ n_uncomputed_tgts])
  end function CWP_Uncomputed_tgts_get


  !>
  !!
  !! \brief Return the number of computed targets.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] field_id         Field identifier
  !! \param [in] i_part           Current partition
  !!
  !! \return                Number of computed targets
  !!

  function CWP_N_computed_tgts_get (local_code_name, &
                                    cpl_id,          &
                                    field_id,        &
                                    i_part)          &
                                    result (n_computed_tgts)
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: cpl_id
    character(kind = c_char, len = *) :: field_id
    integer(c_int) :: i_part

    integer(c_int) :: n_computed_tgts

    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    n_computed_tgts = CWP_N_computed_tgts_get_cf (local_code_name,   &
                                                  l_local_code_name, &
                                                  cpl_id,            &
                                                  l_cpl_id,          &
                                                  field_id,          &
                                                  l_field_id,        &
                                                  i_part)
  end function CWP_N_computed_tgts_get


  !>
  !!
  !! \brief Return computed targets.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] field_id         Field identifier
  !! \param [in] i_part           Current partition
  !!
  !! \return                Computed targets
  !!

  function CWP_Computed_tgts_get (local_code_name, &
                                  cpl_id,          &
                                  field_id,        &
                                  i_part)          &
                                  result (computed_tgts)

    use, intrinsic :: iso_c_binding
    implicit none


    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: cpl_id
    character(kind = c_char, len = *) :: field_id
    integer(c_int) :: i_part

    integer(c_int), dimension(:), pointer :: computed_tgts

    type(c_ptr) :: cptr_computed_tgts
    integer(c_int) :: n_computed_tgts
    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)


    n_computed_tgts = CWP_N_computed_tgts_get_cf (local_code_name,   &
                                                  l_local_code_name, &
                                                  cpl_id,            &
                                                  l_cpl_id,          &
                                                  field_id,          &
                                                  l_field_id,        &
                                                  i_part)

    cptr_computed_tgts = CWP_Computed_tgts_get_cf(local_code_name,   &
                                                  l_local_code_name, &
                                                  cpl_id,            &
                                                  l_cpl_id,          &
                                                  field_id,          &
                                                  l_field_id,        &
                                                  i_part)


    call c_f_pointer (cptr = cptr_computed_tgts, &
                      fptr = computed_tgts  ,    &
                      shape= [n_computed_tgts])
  end function CWP_Computed_tgts_get


  !>
  !!
  !! \brief Return the number of involved sources.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] field_id         Field identifier
  !! \param [in] i_part           Current partition
  !!
  !! \return                Number of involved sources
  !!

  function CWP_N_involved_srcs_get (local_code_name, &
          cpl_id,          &
          field_id,        &
          i_part)          &
          result (n_involved_srcs)
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: cpl_id
    character(kind = c_char, len = *) :: field_id
    integer(c_int) :: i_part

    integer(c_int) :: n_involved_srcs

    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    n_involved_srcs = CWP_N_involved_srcs_get_cf (local_code_name,   &
            l_local_code_name, &
            cpl_id,            &
            l_cpl_id,          &
            field_id,          &
            l_field_id,        &
            i_part)
  end function CWP_N_involved_srcs_get


  !>
  !!
  !! \brief Return involved sources.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] field_id         Field identifier
  !! \param [in] i_part           Current partition
  !!
  !! \return                Involved sources
  !!

  function CWP_Involved_srcs_get (local_code_name, &
          cpl_id,          &
          field_id,        &
          i_part)          &
          result (involved_srcs)

    use, intrinsic :: iso_c_binding
    implicit none


    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: cpl_id
    character(kind = c_char, len = *) :: field_id
    integer(c_int) :: i_part

    integer(c_int), dimension(:), pointer :: involved_srcs

    type(c_ptr) :: cptr_involved_srcs
    integer(c_int) :: n_involved_srcs
    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)


    n_involved_srcs = CWP_N_involved_srcs_get_cf (local_code_name,   &
            l_local_code_name, &
            cpl_id,            &
            l_cpl_id,          &
            field_id,          &
            l_field_id,        &
            i_part)

    cptr_involved_srcs = CWP_Involved_srcs_get_cf(local_code_name,   &
            l_local_code_name, &
            cpl_id,            &
            l_cpl_id,          &
            field_id,          &
            l_field_id,        &
            i_part)


    call c_f_pointer (cptr = cptr_involved_srcs, &
            fptr = involved_srcs  ,    &
            shape= [n_involved_srcs])
  end function CWP_Involved_srcs_get


  !>
  !! \brief Return distance from each target to the source interface. <b>(Not implemented yet)</b>
  !!
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !!
  !! \return               Distance
  !!
  !!

  function CWP_Computed_tgts_dist_to_spatial_interp_get (local_code_name, &
                                                         cpl_id)          &
                                                         result (dists)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: l_local_code_name, l_cpl_id
    double precision, dimension(:), pointer :: dists
!!    type(c_ptr) :: cptr_dists

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    !cptr_dists = CWP_Computed_tgts_dist_to_spatial_interp_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id)

    ! TODO The types of return variables may probably not be the right ones
    print *, "CWP_Computed_tgts_dist_to_spatial_interp_get not implemented"
    allocate(dists(1))
    dists = (/0./)
  end function CWP_Computed_tgts_dist_to_spatial_interp_get


  !>
  !! \brief Compute spatial interpolation weights.
  !!
  !! \param [in]  local_code_name     Local code name
  !! \param [in]  cpl_id              Coupling identifier
  !!

  subroutine CWP_Spatial_interp_weights_compute (local_code_name, &
                                                 cpl_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Spatial_interp_weights_compute_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id)
  end subroutine CWP_Spatial_interp_weights_compute


  !>
  !! \brief Enable visualization output.
  !!
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !! \param [in]  freq             Output frequency
  !! \param [in]  format           Output format to visualize exchanged fieldsDouble
  !!                               on the coupled mesh. Choice between :
  !!                               - "EnSight Gold"
  !! \param [in]  format_option   Output options "opt1, opt2, ..."
  !!                               - text : output text files
  !!                               - binary : output binary files (default)
  !!

  subroutine CWP_Visu_set (local_code_name, &
                           cpl_id,          &
                           freq,            &
                           format,          &
                           format_option)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, format_option
    integer(c_int) :: freq
    integer(c_int) :: format
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_format_option

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_format_option = len(format_option)

    call CWP_Visu_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, freq, format, format_option, l_format_option)
  end subroutine


  !>
  !! \brief Setting user target points.
  !!
  !! This function must be called if the degrees of freedom locations are
  !! \ref CWP_DOF_LOCATION_USER
  !!
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !! \param [in]  i_part           Current partition
  !! \param [in]  n_pts            Number of points
  !! \param [in]  coord            Coordinates (size = 3 * n_pts)
  !! \param [in]  g_num            global number or NUL (size = n_pts)
  !!

  subroutine CWP_User_tgt_pts_set(local_code_name, &
                                  cpl_id, &
                                  i_part, &
                                  n_pts, &
                                  coord, &
                                  global_num)
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id
    integer(kind = c_int) :: i_part, n_pts
    double precision, dimension(:), pointer :: coord
    integer(kind = c_long), dimension(:), pointer :: global_num

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_User_tgt_pts_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, n_pts, &
            c_loc(coord), c_loc(global_num))
  end subroutine CWP_User_tgt_pts_set


  !>
  !! \brief Finalize interface mesh.
  !!
  !! This function computes the global numbers of mesh entities if they are
  !! not provided.
  !!
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !!

  subroutine CWP_Mesh_interf_finalize (local_code_name, &
                                       cpl_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_finalize_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id)
  end subroutine CWP_Mesh_interf_finalize


  !>
  !! \brief Set vertices.
  !!
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !! \param [in]  i_part           Current partition
  !! \param [in]  n_pts            Number of points
  !! \param [in]  coord            Coordinates (size = 3    !! \p n_pts)
  !! \param [in]  global_num       Pointer to parent element number (or NULL)
  !!

  subroutine CWP_Mesh_interf_vtx_set (local_code_name, &
                                      cpl_id,          &
                                      i_part,          &
                                      n_pts,           &
                                      coord,           &
                                      global_num)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: i_part, n_pts
    double precision, dimension(:), pointer :: coord
    integer(c_long), dimension(:), pointer :: global_num
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_vtx_set_cf (local_code_name,   &
                                     l_local_code_name, &
                                     cpl_id,            &
                                     l_cpl_id,          &
                                     i_part,            &
                                     n_pts,             &
                                     c_loc(coord),      &
                                     c_loc(global_num))
  end subroutine CWP_Mesh_interf_vtx_set


  !>
  !! \brief Add a connectivity block to the interface mesh.
  !!
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !! \param [in]  block_type       Block type
  !!
  !! \return block identifier
  !!

  function CWP_Mesh_interf_block_add (local_code_name, &
                                      cpl_id,          &
                                      block_type)      &
                                      result(block_id)
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: block_type
    integer(c_int) :: block_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    block_id = CWP_Mesh_interf_block_add_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, block_type)
  end function CWP_Mesh_interf_block_add


  !>
  !! \brief Set a standard block to the interface mesh.
  !!
  !! This function adds a connectivity block to the interface mesh.
  !! Definition of element connectivity is :
  !!
  !!  - edge (\ref CWP_BLOCK_EDGE2) :
  !!
  !!   \code
  !!       1 x-------x 2
  !!   \endcode
  !!
  !!  - triangle (\ref CWP_BLOCK_FACE_TRIA3):
  !!
  !!   \code
  !!       1 x-------x 3
  !!          \     /
  !!           \   /
  !!            \ /
  !!             x 2
  !!   \endcode
  !!
  !!  - quadrangle (\ref CWP_BLOCK_FACE_QUAD4) :
  !!
  !!   \code
  !!          4 x-------x 3
  !!           /       /
  !!          /       /
  !!       1 x-------x2
  !!   \endcode
  !!
  !!   - tetrahedron (\ref CWP_BLOCK_CELL_TETRA4) :
  !!
  !!   \code
  !!             x 4
  !!            /|\
  !!           / | \
  !!          /  |  \
  !!       1 x- -|- -x 3
  !!          \  |  /
  !!           \ | /
  !!            \|/
  !!             x 2
  !!   \endcode
  !!
  !!   - pyramid (\ref CWP_BLOCK_CELL_PYRAM5) :
  !!
  !!   \code
  !!              5 x
  !!               /|\
  !!              //| \
  !!             // |  \
  !!          4 x/--|---x 3
  !!           //   |  /
  !!          //    | /
  !!       1 x-------x 2
  !!   \endcode
  !!
  !!  - prism (\ref CWP_BLOCK_CELL_PRISM6) :
  !!
  !!   \code
  !!       4 x-------x 6
  !!         |\     /|
  !!         | \   / |
  !!       1 x- \-/ -x 3
  !!          \ 5x  /
  !!           \ | /
  !!            \|/
  !!             x 2
  !!   \endcode
  !!
  !!  -  hexaedron (\ref CWP_BLOCK_CELL_HEXA8) :
  !!
  !!   \code
  !!          8 x-------x 7
  !!           /|      /|
  !!          / |     / |
  !!       5 x-------x6 |
  !!         | 4x----|--x 3
  !!         | /     | /
  !!         |/      |/
  !!       1 x-------x 2
  !!   \endcode
  !!
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !! \param [in]  i_part           Partition identifier
  !! \param [in]  block_id         Block identifier
  !! \param [in]  n_elts           Number of elements
  !! \param [in]  connec           Connectivity (size = n_vertex_elt * n_elts)
  !! \param [in]  global_num       Pointer to global element number (or NULL)
  !!

  subroutine CWP_Mesh_interf_block_std_set (local_code_name, &
                                            cpl_id,          &
                                            i_part,          &
                                            block_id,        &
                                            n_elts,          &
                                            connec,          &
                                            global_num)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: i_part, block_id, n_elts
    integer(c_int), dimension(:), pointer :: connec
    integer(c_long), dimension(:), pointer :: global_num
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_block_std_set_cf (local_code_name,   &
                                           l_local_code_name, &
                                           cpl_id,            &
                                           l_cpl_id,          &
                                           i_part,            &
                                           block_id,          &
                                           n_elts,            &
                                           c_loc(connec),     &
                                           c_loc(global_num))
  end subroutine CWP_Mesh_interf_block_std_set


  !>
  !! \brief Set the connectivity of a polygon block in a interface mesh partition.
  !!
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !! \param [in]  i_part           Current partition
  !! \param [in]  block_id         Block identifier
  !! \param [in]  n_elts           Number of elements
  !! \param [in]  connec_idx       Connectivity index (\p connec_id[0] = 0 and
  !!                               size = \p n_elts + 1)
  !! \param [in]  connec           Connectivity (size = \p connec_idx[\p n_elts])
  !! \param [in]  global_num       Pointer to global element number (or NULL)
  !!

  subroutine CWP_Mesh_interf_f_poly_block_set( local_code_name, &
                                               cpl_id,          &
                                               i_part,          &
                                               block_id,        &
                                               n_elts,          &
                                               connec_idx,      &
                                               connec,          &
                                               global_num)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: i_part, block_id, n_elts
    integer(c_int), dimension(:), pointer :: connec_idx, connec
    integer(c_long), dimension(:), pointer :: global_num
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_f_poly_block_set_cf (local_code_name,   &
                                              l_local_code_name, &
                                              cpl_id,            &
                                              l_cpl_id,          &
                                              i_part,            &
                                              block_id,          &
                                              n_elts,            &
                                              c_loc(connec_idx), &
                                              c_loc(connec),     &
                                              c_loc(global_num))
  end subroutine CWP_Mesh_interf_f_poly_block_set


  !>
  !! \brief Adding a polyhedron connectivity block to the interface mesh.
  !!
  !! \param [in]  local_code_name   Local code name
  !! \param [in]  cpl_id            Coupling identifier
  !! \param [in]  i_part            Current partition
  !! \param [in]  block_id          Block identifier
  !! \param [in]  n_elts            Number of elements
  !! \param [in]  connec_cells_idx  Polyhedron to face index
  !!                                (\p src_poly_cell_face_idx[0] = 0 and
  !!                                 size = \p n_elts + 1)
  !! \param [in]  connec_cells      Polyhedron to face connectivity
  !!                                (size = \p cell_face_idx[\p n_elts])
  !! \param [in]  n_faces           Number of faces
  !! \param [in]  connec_faces_idx  Polyhedron face to vertex index
  !!                                (\p face_vertex_idx[0] = 0 and
  !!                                 size = max(\p cell_face_connec) + 1)
  !! \param [in]  connec_faces      Polyhedron face to vertex connectivity
  !!                                (size = \p face_vertex_idx[\p n_elts])
  !! \param [in]  global_num        Pointer to global element number (or NULL)
  !!

  subroutine CWP_Mesh_interf_c_poly_block_set (local_code_name, &
                                               cpl_id,          &
                                               i_part,          &
                                               block_id,        &
                                               n_elts,          &
                                               n_faces,         &
                                               connec_faces_idx,&
                                               connec_faces,    &
                                               connec_cells_idx,&
                                               connec_cells,    &
                                               global_num)
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: i_part, block_id, n_elts, n_faces
    integer(c_int), dimension(:), pointer :: connec_faces_idx, connec_faces, connec_cells_idx, connec_cells
    integer(c_long), dimension(:), pointer :: global_num
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_c_poly_block_set_cf (local_code_name,         &
                                              l_local_code_name,       &
                                              cpl_id,                  &
                                              l_cpl_id,                &
                                              i_part,                  &
                                              block_id,                &
                                              n_elts,                  &
                                              n_faces,                 &
                                              c_loc(connec_faces_idx), &
                                              c_loc(connec_faces),     &
                                              c_loc(connec_cells_idx), &
                                              c_loc(connec_cells),     &
                                              c_loc(global_num))
  end subroutine CWP_Mesh_interf_c_poly_block_set


  !>
  !! \brief Delete interface mesh.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !!
  !!

  subroutine CWP_Mesh_interf_del (local_code_name, &
                                  cpl_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)

    call CWP_Mesh_interf_del_cf (local_code_name, l_local_code_name, cpl_id, l_cpl_id)
  end subroutine CWP_Mesh_interf_del


  !>
  !! \brief Define the interface mesh from a cell to face connectivity.
  !!
  !! \param [in]  local_code_name   Local code name
  !! \param [in]  cpl_id            Coupling identifier
  !! \param [in]  i_part            Current partition
  !! \param [in]  n_cells           Number of cells
  !! \param [in]  cell_face_idx     Polyhedron to face index
  !!                                (\p src_poly_cell_face_idx[0] = 0 and
  !!                                 size = \p n_elts + 1)
  !! \param [in]  cell_face         Cell to face connectivity
  !!                                (size = \p cell_face_idx[\p n_elts])
  !! \param [in]  n_faces           Number of faces
  !! \param [in]  face_vtx_idx      Polyhedron face to vertex index
  !!                                (\p face_vtx_idx[0] = 0 and
  !!                                 size = \p n_faces + 1)
  !! \param [in]  face_vtx          Face to vertex connectivity
  !!                                (size = \p face_vtx_idx[\p n_elts])
  !! \param [in]  parent_num        Pointer to parent element number (or NULL)
  !!

  subroutine CWP_Mesh_interf_from_cellface_set (local_code_name, &
                                                cpl_id,          &
                                                i_part,          &
                                                n_cells,         &
                                                cell_face_idx,   &
                                                cell_face,       &
                                                n_faces,         &
                                                face_vtx_idx,    &
                                                face_vtx,        &
                                                parent_num)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: i_part, n_cells, n_faces
    integer(c_int), dimension(:), pointer :: cell_face_idx, cell_face, face_vtx_idx, face_vtx, parent_num
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_from_cellface_set_cf (local_code_name,      &
                                               l_local_code_name,    &
                                               cpl_id,               &
                                               l_cpl_id,             &
                                               i_part,               &
                                               n_cells,              &
                                               c_loc(cell_face_idx), &
                                               c_loc(cell_face),     &
                                               n_faces,              &
                                               c_loc(face_vtx_idx),  &
                                               c_loc(face_vtx),      &
                                               c_loc(parent_num))
  end subroutine CWP_Mesh_interf_from_cellface_set


  !>
  !! \brief Define the surface interface mesh from a face to edge connectivity.
  !!
  !! \param [in]  local_code_name   Local code name
  !! \param [in]  cpl_id            Coupling identifier
  !! \param [in]  i_part            Current partition
  !! \param [in]  n_faces           Number of cells
  !! \param [in]  face_edge_idx     Polygon to edge index
  !!                                (\p face_edge_idx[0] = 0 and
  !!                                 size =  \p n_faces + 1)
  !! \param [in]  face_edge         Face to edge connectivity
  !!                                (size = \p face_edge_idx[\p n_faces])
  !! \param [in]  n_edges           Number of faces
  !! \param [in]  edge_vtx_idx      Polyhedron face to vertex index
  !!                                (\p edge_vtx_idx[0] = 0 and
  !!                                 size = \p n_edges + 1)
  !! \param [in]  edge_vtx          Face to vertex connectivity
  !!                                (size = \p edge_vtx_idx[\p n_edges])
  !! \param [in]  parent_num        Pointer to parent element number (or NULL)
  !!

  subroutine CWP_Mesh_interf_from_faceedge_set (local_code_name, &
                                                cpl_id,          &
                                                i_part,          &
                                                n_faces,         &
                                                face_edge_idx,   &
                                                face_edge,       &
                                                n_edges,         &
                                                edge_vtx_idx,    &
                                                edge_vtx,        &
                                                parent_num)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: i_part, n_faces, n_edges
    integer(c_int), dimension(:), pointer :: face_edge_idx, face_edge, edge_vtx_idx, edge_vtx, parent_num
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_from_faceedge_set_cf (local_code_name,     &
                                               l_local_code_name,   &
                                               cpl_id,              &
                                               l_cpl_id,            &
                                               i_part,              &
                                               n_faces,             &
                                               c_loc(face_edge_idx),&
                                               c_loc(face_edge),    &
                                               n_edges,             &
                                               c_loc(edge_vtx_idx), &
                                               c_loc(edge_vtx),     &
                                               c_loc(parent_num))
  end subroutine CWP_Mesh_interf_from_faceedge_set


  !>
  !!
  !! \brief Create a new field.
  !!
  !! \param [in]  local_code_name Local code name
  !! \param [in]  cpl_id          Coupling identifier
  !! \param [in]  field_id        Field id
  !! \param [in]  data_type       Data type
  !! \param [in]  storage         Storage type
  !! \param [in]  n_component     Number of component
  !! \param [in]  target_location Target location
  !! \param [in]  exch_type       Exchange type
  !! \param [in]  visu_status     Visualization status
  !!

  subroutine CWP_Field_create (local_code_name,      &
                               cpl_id,               &
                               field_id,             &
                               data_type,            &
                               storage,              &
                               n_component,          &
                               target_location,      &
                               exch_type,            &
                               visu_status)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(c_int) :: data_type, storage, n_component, target_location, exch_type, visu_status
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_create_cf (local_code_name,    &
                              l_local_code_name,  &
                              cpl_id,             &
                              l_cpl_id,           &
                              field_id,           &
                              l_field_id,         &
                              data_type,          &
                              storage,            &
                              n_component,        &
                              target_location,    &
                              exch_type,          &
                              visu_status)
  end subroutine CWP_Field_create


  !>
  !!
  !! \brief Set field data.
  !!
  !! \param [in] local_code_name   Local code name
  !! \param [in] cpl_id            Coupling identifier
  !! \param [in] field_id          Field identifier
  !! \param [in] i_part            Current partition
  !! \param [in] data_type         Choice if data is setted for the source or the target
  !! \param [in] data              Storage array (Mapping)
  !!

  subroutine CWP_Field_data_set (local_code_name, &
                                 cpl_id,          &
                                 field_id,        &
                                 i_part,          &
                                 map_type,        &
                                 data)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(c_int) :: i_part
    integer(c_int) :: map_type
    double precision, dimension(:), pointer :: data
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_data_set_cf (local_code_name,   &
                                l_local_code_name, &
                                cpl_id,            &
                                l_cpl_id,          &
                                field_id,          &
                                l_field_id,        &
                                i_part,            &
                                map_type,          &
                                c_loc(data))
  end subroutine CWP_Field_data_set


  !>
  !! \brief Send a spatially interpolated field to the coupled code with
  !!        nonblocking communications.
  !!
  !! This function is independant of \ref CWP_Time_exch_t mode. The user has to
  !! manually check the consistency of the exchanges.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in]  cpl_id          Coupling identifier
  !! \param [in]  src_field_id    Source field id
  !!

  subroutine CWP_Field_issend (local_code_name, &
                               cpl_id,          &
                               src_field_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, src_field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_src_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_src_field_id = len(src_field_id)

    call CWP_Field_issend_cf (local_code_name,    &
                              l_local_code_name,  &
                              cpl_id,             &
                              l_cpl_id,           &
                              src_field_id,       &
                              l_src_field_id)
  end subroutine CWP_Field_issend


  !>
  !!
  !! \brief Receive a spatially interpolated field from the coupled code
  !!        with nonblocking communications.
  !!
  !! This function is independant of \ref CWP_Time_exch_t mode. The user has to
  !! manually check the consistency of the exchanges.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in]  cpl_id          Coupling identifier
  !! \param [in]  tgt_field_id    Target field id
  !!

  subroutine CWP_Field_irecv (local_code_name, &
                              cpl_id,          &
                              tgt_field_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, tgt_field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_tgt_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_tgt_field_id = len(tgt_field_id)

    call CWP_Field_irecv_cf(local_code_name,    &
                            l_local_code_name,  &
                            cpl_id,             &
                            l_cpl_id,           &
                            tgt_field_id,       &
                            l_tgt_field_id)
  end subroutine CWP_Field_irecv


  !>
  !!
  !! \brief Wait the end of an exchange related to request from \ref CWP_Field_issend.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] src_field_id     Source field id
  !!

  subroutine CWP_Field_wait_issend (local_code_name, &
                                    cpl_id,          &
                                    src_field_id)
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, src_field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_src_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_src_field_id = len(src_field_id)

    call CWP_Field_wait_issend_cf (local_code_name,    &
                                   l_local_code_name,  &
                                   cpl_id,             &
                                   l_cpl_id,           &
                                   src_field_id,       &
                                   l_src_field_id)
  end subroutine CWP_Field_wait_issend


  !>
  !!
  !! \brief Wait the end of an exchange related to request from \ref CWP_Field_irecv.
  !!
  !! This function waits the end of exchange related to request
  !! from \ref CWP_Field_irecv
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] tgt_field_id     Target field id
  !!

  subroutine CWP_Field_wait_irecv (local_code_name, &
                                   cpl_id,          &
                                   tgt_field_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, tgt_field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_tgt_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_tgt_field_id = len(tgt_field_id)

    call CWP_Field_wait_irecv_cf(local_code_name,     &
                                 l_local_code_name,   &
                                 cpl_id,              &
                                 l_cpl_id,            &
                                 tgt_field_id,        &
                                 l_tgt_field_id)
  end subroutine CWP_Field_wait_irecv


  !>
  !!
  !! \brief Setting of an user interpolation from location.
  !!
  !! This function takes into account an user interpolation function written with
  !! void (*\ref CWP_Interp_from_location_t) interface.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] src_field_id     Source field id
  !! \param [in] fct              Function
  !!

  subroutine CWP_Interp_from_location_set(local_code_name, &
                                          cpl_id, &
                                          src_field_id, &
                                          ptInterpolationFct)

    use, intrinsic :: iso_c_binding
    implicit none

    abstract interface
      subroutine ptInterpolationFct ( &
              interface_type, &
              n_src_vtcs, &
              n_src_std_elts, &
              n_tgt_pts, &
              src_vtcs_coords, &
              src_connec_idx, &
              src_connec, &
              tgt_pts_coords, &
              tgt_pts_target_location, &
              tgt_pts_dist, &
              tgt_pts_bary_coords_idx, &
              tgt_pts_bary_coords, &
              stride, &
              src_field_location, &
              src_field, &
              tgt_field_location, &
              tgt_field &
              )
        use, intrinsic :: iso_c_binding
        implicit none

        integer(kind = c_int)               :: interface_type
        integer(kind = c_int)               :: n_src_vtcs
        integer(kind = c_int)               :: n_src_std_elts
        integer(kind = c_int)               :: n_tgt_pts
        real(kind = c_double), dimension(*) :: src_vtcs_coords
        integer(kind = c_int), dimension(*) :: src_connec_idx
        integer(kind = c_int), dimension(*) :: src_connec
        real(kind = c_double), dimension(*) :: tgt_pts_coords
        integer(kind = c_int), dimension(*) :: tgt_pts_target_location
        real(kind = c_double), dimension(*) :: tgt_pts_dist
        integer(kind = c_int), dimension(*) :: tgt_pts_bary_coords_idx
        real(kind = c_double), dimension(*) :: tgt_pts_bary_coords
        integer(kind = c_int)               :: stride
        integer(kind = c_int)               :: src_field_location
        real(kind = c_double), dimension(*) :: src_field
        integer(kind = c_int)               :: tgt_field_location
        real(kind = c_double), dimension(*) :: tgt_field
      end subroutine ptInterpolationFct
    end interface

    character(kind = c_char, len = *) :: local_code_name, cpl_id, src_field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_src_field_id
    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_src_field_id = len(src_field_id)

    ! call CWP_Interp_from_location_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
    !         & src_field_id, l_src_field_id, ptInterpolationFct)
  end subroutine CWP_Interp_from_location_set
end module cwp