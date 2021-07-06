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
                CWP_TEST, &
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

        function CWP_N_uncomputed_tgts_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
                target_location, i_part) result (n_uncomputed_tgts) &
                bind(c, name = 'CWP_N_uncomputed_tgts_get_cf')
            use, intrinsic :: iso_c_binding
            implicit none

            character(kind = c_char, len = 1) :: local_code_name, cpl_id
            integer(kind = c_int), value :: l_local_code_name, l_cpl_id
            integer(c_int) :: target_location, i_part
            integer(kind = c_int) :: n_uncomputed_tgts
        end function CWP_N_uncomputed_tgts_get_cf

        function CWP_Uncomputed_tgts_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id) result (uncomputed_tgts) &
                bind(c, name = 'CWP_Uncomputed_tgts_get_cf')
            use, intrinsic :: iso_c_binding
            implicit none

            character(kind = c_char, len = 1) :: local_code_name, cpl_id
            integer(kind = c_int), value :: l_local_code_name, l_cpl_id
            type(c_ptr) :: uncomputed_tgts
        end function CWP_Uncomputed_tgts_get_cf

        function CWP_N_computed_tgts_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id) result (n_computed_tgts) &
                bind(c, name = 'CWP_N_computed_tgts_get_cf')
            use, intrinsic :: iso_c_binding
            implicit none

            character(kind = c_char, len = 1) :: local_code_name, cpl_id
            integer(kind = c_int), value :: l_local_code_name, l_cpl_id
            integer(kind = c_int) :: n_computed_tgts
        end function CWP_N_computed_tgts_get_cf

        function CWP_Computed_tgts_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id) result (computed_tgts) &
                bind(c, name = 'CWP_Computed_tgts_get_cf')
            use, intrinsic :: iso_c_binding
            implicit none

            character(kind = c_char, len = 1) :: local_code_name, cpl_id
            integer(kind = c_int), value :: l_local_code_name, l_cpl_id
            type(c_ptr) :: computed_tgts
        end function CWP_Computed_tgts_get_cf

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
            integer(c_int) :: i_part, block_id, n_elts
            type(c_ptr) :: connec, global_num
            integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        end subroutine CWP_Mesh_interf_block_std_set_cf

        subroutine CWP_Mesh_interf_f_poly_block_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
                i_part, block_id, n_elts, connec_idx, connec, global_num) &
            bind(c, name = 'CWP_Mesh_interf_f_poly_block_set')
            use, intrinsic :: iso_c_binding
            implicit none

            character(kind = c_char, len = 1) :: local_code_name, cpl_id
            integer(c_int) :: i_part, block_id, n_elts
            type(c_ptr) :: connec_idx, connec, global_num
            integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        end subroutine CWP_Mesh_interf_f_poly_block_set_cf

        subroutine CWP_Mesh_interf_c_poly_block_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
                i_part, block_id, n_elts, n_faces, connec_faces_idx, connec_faces, connec_cells_idx, connec_cells, global_num) &
            bind(c, name = 'CWP_Mesh_interf_c_poly_block_set_cf')
            use, intrinsic :: iso_c_binding
            implicit none

            character(kind = c_char, len = 1) :: local_code_name, cpl_id
            integer(c_int) :: i_part, block_id, n_elts, n_faces
            type(c_ptr) :: connec_faces_idx, connec_faces, connec_cells_idx, connec_cells, global_num
            integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        end subroutine CWP_Mesh_interf_c_poly_block_set_cf

        subroutine CWP_Mesh_interf_from_cellface_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, n_cells, &
                cell_face_idx, cell_face, n_faces, face_vtx_idx, face_vtx, parent_num) &
                bind(c, name = 'CWP_Mesh_interf_from_cellface_set_cf')
            use, intrinsic :: iso_c_binding
            implicit none

            character(kind = c_char, len = 1) :: local_code_name, cpl_id
            integer(c_int) :: i_part, n_cells, n_faces
            type(c_ptr) :: cell_face_idx, cell_face, face_vtx_idx, face_vtx, parent_num
            integer(kind = c_int) :: l_local_code_name, l_cpl_id
        end subroutine CWP_Mesh_interf_from_cellface_set_cf

        subroutine CWP_Mesh_interf_from_faceedge_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, n_faces, &
                face_edge_idx, face_edge, n_edges, edge_vtx_idx, edge_vtx, parent_num) &
                bind(c, name = 'CWP_Mesh_interf_from_faceedge_set_cf')
            use, intrinsic :: iso_c_binding
            implicit none

            character(kind = c_char, len = 1) :: local_code_name, cpl_id
            integer(c_int) :: i_part, n_faces, n_edges
            type(c_ptr) :: face_edge_idx, face_edge, edge_vtx_idx, edge_vtx, parent_num
            integer(kind = c_int) :: l_local_code_name, l_cpl_id
        end subroutine CWP_Mesh_interf_from_faceedge_set_cf

        subroutine CWP_Field_create_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id, &
                data_type, storage, n_component, target_location, exch_type, visu_status) &
                bind(c, name = 'CWP_Field_create_cf')
            use, intrinsic :: iso_c_binding
            implicit none

            character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
            integer(c_int) :: data_type, storage, n_component, target_location, exch_type, visu_status
            integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_field_id
        end subroutine

        subroutine CWP_Field_data_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id, &
                i_part, data) &
                bind(c, name = 'CWP_Field_data_set_cf')
            use, intrinsic :: iso_c_binding
            implicit none

            character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
            integer(c_int) :: i_part
            type(c_ptr) :: data
            integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_field_id
        end subroutine CWP_Field_data_set_cf

        subroutine CWP_Field_issend_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, src_field_id, l_src_field_id) &
                bind(c, name = 'CWP_Field_irecv_cf')
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
                bind(c, name = 'CWP_Field_wait_irecv_cf')
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

        subroutine CWP_Finalize() &
                bind(c, name = 'CWP_Finalize')
        end subroutine CWP_Finalize
    end interface

contains
    subroutine CWP_Init(fcomm, n_code, code_names, is_active_rank, time_init, intra_comms)
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
        call CWP_Init_cf(fcomm, n_code, c_loc(code_names), l_code_names, is_active_rank, time_init, c_loc(intra_comms))
    end subroutine CWP_Init

    subroutine CWP_Cpl_create(local_code_name, cpl_id, coupled_code_name, &
            entities_dim, comm_type, spatial_interp, n_part, displacement, freq)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id, coupled_code_name
        integer(kind = c_int), value :: entities_dim
        integer(kind = c_int), value :: comm_type
        integer(kind = c_int) :: spatial_interp
        integer(kind = c_int), value :: n_part
        integer(kind = c_int), value :: displacement
        integer(kind = c_int), value :: freq
        integer(c_int) :: l_local_code_name, l_cpl_id, l_coupled_code_name

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)
        l_coupled_code_name = len(coupled_code_name)
        call CWP_Cpl_create_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, coupled_code_name, l_coupled_code_name, &
                entities_dim, comm_type, spatial_interp, n_part, displacement, freq)
    end subroutine CWP_Cpl_Create

    subroutine CWP_Cpl_Del(local_code_name, cpl_id)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id
        integer(c_int) :: l_local_code_name, l_cpl_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)
        call CWP_Cpl_del_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id)
    end subroutine CWP_Cpl_Del

    function CWP_N_uncomputed_tgts_get(local_code_name, cpl_id, target_location, i_part) result (n_uncomputed_tgts)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id
        integer(c_int) :: target_location, i_part
        integer(c_int) :: n_uncomputed_tgts
        integer(c_int) :: l_local_code_name, l_cpl_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)
        n_uncomputed_tgts = CWP_N_uncomputed_tgts_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
                target_location, i_part)
    end function CWP_N_uncomputed_tgts_get

    function CWP_Uncomputed_tgts_get(local_code_name, cpl_id) result (uncomputed_tgts)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id
        integer(c_int) :: l_local_code_name, l_cpl_id
        integer(c_int), dimension(:), pointer :: uncomputed_tgts
        type(c_ptr) :: cptr_uncomputed_tgts

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)

        cptr_uncomputed_tgts = CWP_Uncomputed_tgts_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id)

        ! TODO The types of return variables may probably not be the right ones
        print *, "CWP_Uncomputed_tgts_get not implemented"
        allocate(uncomputed_tgts(1))
        uncomputed_tgts = (/0/)
    end function CWP_Uncomputed_tgts_get

    function CWP_N_computed_tgts_get(local_code_name, cpl_id) result (n_computed_tgts)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id
        integer(c_int) :: n_computed_tgts
        integer(c_int) :: l_local_code_name, l_cpl_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)
        n_computed_tgts = CWP_N_computed_tgts_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id)
    end function CWP_N_computed_tgts_get

    function CWP_Computed_tgts_get(local_code_name, cpl_id) result (computed_tgts)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id
        integer(c_int) :: l_local_code_name, l_cpl_id
        integer(c_int), dimension(:), pointer :: computed_tgts
        type(c_ptr) :: cptr_computed_tgts

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)

        cptr_computed_tgts = CWP_Computed_tgts_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id)

        ! TODO The types of return variables may probably not be the right ones
        print *, "CWP_Computed_tgts_get not implemented"
        allocate(computed_tgts(1))
        computed_tgts = (/0/)
    end function CWP_Computed_tgts_get

    function CWP_Computed_tgts_dist_to_spatial_interp_get(local_code_name, cpl_id) result (dists)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id
        integer(c_int) :: l_local_code_name, l_cpl_id
        double precision, dimension(:), pointer :: dists
        type(c_ptr) :: cptr_dists

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)

        cptr_dists = CWP_Computed_tgts_dist_to_spatial_interp_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id)
        ! TODO The types of return variables may probably not be the right ones
        print *, "CWP_Computed_tgts_dist_to_spatial_interp_get not implemented"
        allocate(dists(1))
        dists = (/0./)
    end function CWP_Computed_tgts_dist_to_spatial_interp_get

    subroutine CWP_Spatial_interp_weights_compute(local_code_name, cpl_id)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id
        integer(kind = c_int) :: l_local_code_name, l_cpl_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)

        call CWP_Spatial_interp_weights_compute_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id)
    end subroutine CWP_Spatial_interp_weights_compute

    subroutine CWP_Visu_set(local_code_name, cpl_id, freq, format, format_option)
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

    subroutine CWP_Mesh_interf_finalize(local_code_name, cpl_id)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id
        integer(kind = c_int) :: l_local_code_name, l_cpl_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)

        call CWP_Mesh_interf_finalize_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id)
    end subroutine CWP_Mesh_interf_finalize

    subroutine CWP_Mesh_interf_vtx_set(local_code_name, cpl_id, i_part, n_pts, coord, global_num)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id
        integer(c_int) :: i_part, n_pts
        double precision, dimension(:), pointer :: coord
        integer(c_int), dimension(:), pointer :: global_num
        integer(kind = c_int) :: l_local_code_name, l_cpl_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)

        call CWP_Mesh_interf_vtx_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, n_pts, &
                c_loc(coord), c_loc(global_num))
    end subroutine CWP_Mesh_interf_vtx_set

    function CWP_Mesh_interf_block_add(local_code_name, cpl_id, block_type) result(block_id)
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

    subroutine CWP_Mesh_interf_block_std_set(local_code_name, cpl_id, i_part, block_id, n_elts, connec, global_num)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id
        integer(c_int) :: i_part, block_id, n_elts
        integer(c_int), dimension(:), pointer :: connec, global_num
        integer(kind = c_int) :: l_local_code_name, l_cpl_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)

        call CWP_Mesh_interf_block_std_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, block_id, &
                n_elts, c_loc(connec), c_loc(global_num))
    end subroutine CWP_Mesh_interf_block_std_set

    subroutine CWP_Mesh_interf_f_poly_block_set(local_code_name, cpl_id, i_part, block_id, n_elts, connec_idx, connec, global_num)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id
        integer(c_int) :: i_part, block_id, n_elts
        integer(c_int), dimension(:), pointer :: connec_idx, connec, global_num
        integer(kind = c_int) :: l_local_code_name, l_cpl_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)

        call CWP_Mesh_interf_f_poly_block_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, block_id, &
                n_elts, c_loc(connec_idx), c_loc(connec), c_loc(global_num))
    end subroutine CWP_Mesh_interf_f_poly_block_set

    subroutine CWP_Mesh_interf_c_poly_block_set(local_code_name, cpl_id, i_part, block_id, n_elts, n_faces, &
            connec_faces_idx, connec_faces, connec_cells_idx, connec_cells, global_num)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id
        integer(c_int) :: i_part, block_id, n_elts, n_faces
        integer(c_int), dimension(:), pointer :: connec_faces_idx, connec_faces, connec_cells_idx, connec_cells, global_num
        integer(kind = c_int) :: l_local_code_name, l_cpl_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)

        call CWP_Mesh_interf_c_poly_block_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
                i_part, block_id, n_elts, n_faces, c_loc(connec_faces_idx), c_loc(connec_faces), &
                c_loc(connec_cells_idx), c_loc(connec_cells), c_loc(global_num))
    end subroutine CWP_Mesh_interf_c_poly_block_set

    subroutine CWP_Mesh_interf_from_cellface_set(local_code_name, cpl_id, i_part, n_cells, &
            cell_face_idx, cell_face, n_faces, face_vtx_idx, face_vtx, parent_num)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id
        integer(c_int) :: i_part, n_cells, n_faces
        integer(c_int), dimension(:), pointer :: cell_face_idx, cell_face, face_vtx_idx, face_vtx, parent_num
        integer(kind = c_int) :: l_local_code_name, l_cpl_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)

        call CWP_Mesh_interf_from_cellface_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, n_cells, &
                c_loc(cell_face_idx), c_loc(cell_face), n_faces, c_loc(face_vtx_idx), c_loc(face_vtx), c_loc(parent_num))
    end subroutine CWP_Mesh_interf_from_cellface_set

    subroutine CWP_Mesh_interf_from_faceedge_set(local_code_name, cpl_id, i_part, n_faces, &
            face_edge_idx, face_edge, n_edges, edge_vtx_idx, edge_vtx, parent_num)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id
        integer(c_int) :: i_part, n_faces, n_edges
        integer(c_int), dimension(:), pointer :: face_edge_idx, face_edge, edge_vtx_idx, edge_vtx, parent_num
        integer(kind = c_int) :: l_local_code_name, l_cpl_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)

        call CWP_Mesh_interf_from_faceedge_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, n_faces, &
                c_loc(face_edge_idx), c_loc(face_edge), n_edges, c_loc(edge_vtx_idx), c_loc(edge_vtx), c_loc(parent_num))
    end subroutine CWP_Mesh_interf_from_faceedge_set

    subroutine CWP_Field_create(local_code_name, cpl_id, field_id, &
            data_type, storage, n_component, target_location, exch_type, visu_status)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
        integer(c_int) :: data_type, storage, n_component, target_location, exch_type, visu_status
        integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)
        l_field_id = len(field_id)

        call CWP_Field_create_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id, &
                data_type, storage, n_component, target_location, exch_type, visu_status)
    end subroutine CWP_Field_create

    subroutine CWP_Field_data_set(local_code_name, cpl_id, field_id, i_part, data)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
        integer(c_int) :: i_part
        double precision, dimension(:), pointer :: data
        integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)
        l_field_id = len(field_id)

        call CWP_Field_data_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id, &
                i_part, c_loc(data))
    end subroutine CWP_Field_data_set

    subroutine CWP_Field_issend(local_code_name, cpl_id, src_field_id)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id, src_field_id
        integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_src_field_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)
        l_src_field_id = len(src_field_id)
        call CWP_Field_issend_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, src_field_id, l_src_field_id)
    end subroutine CWP_Field_issend

    subroutine CWP_Field_irecv(local_code_name, cpl_id, tgt_field_id)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id, tgt_field_id
        integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_tgt_field_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)
        l_tgt_field_id = len(tgt_field_id)
        call CWP_Field_irecv_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, tgt_field_id, l_tgt_field_id)
    end subroutine CWP_Field_irecv

    subroutine CWP_Field_wait_issend(local_code_name, cpl_id, src_field_id)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id, src_field_id
        integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_src_field_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)
        l_src_field_id = len(src_field_id)
        call CWP_Field_wait_issend_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, src_field_id, l_src_field_id)
    end subroutine CWP_Field_wait_issend

    subroutine CWP_Field_wait_irecv(local_code_name, cpl_id, tgt_field_id)
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = *) :: local_code_name, cpl_id, tgt_field_id
        integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_tgt_field_id

        l_local_code_name = len(local_code_name)
        l_cpl_id = len(cpl_id)
        l_tgt_field_id = len(tgt_field_id)
        call CWP_Field_wait_irecv_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, tgt_field_id, l_tgt_field_id)
    end subroutine CWP_Field_wait_irecv
end module cwp