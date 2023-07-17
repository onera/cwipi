!-----------------------------------------------------------------------------
! This file is part of the CWIPI library.
!
! Copyright (C) 2021-2023  ONERA
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
    use pdm_fortran

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
                CWP_COMM_PAR_WITHOUT_PART
                ! CWP_COMM_INTERNAL
    end enum

    ! CWP_Time_exch_t
    enum, bind(c)
        enumerator :: &
            CWP_TIME_EXCH_USER_CONTROLLED
            ! CWP_TIME_EXCH_EACH_TIME_STEP, &
            ! CWP_TIME_EXCH_N_TIME_STEP, &
            ! CWP_TIME_EXCH_CPL_TIME_STEP, &
            ! CWP_TIME_EXCH_ASYNCHRONOUS, &
            ! CWP_TIME_EXCH_SLAVE, &
            ! CWP_TIME_EXCH_MASTER
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
                CWP_FIELD_STORAGE_INTERLEAVED
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
                CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES, &
                CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES, &
                CWP_SPATIAL_INTERP_FROM_INTERSECTION, &
                CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT, &
                CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE, &
                CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE, &
                CWP_SPATIAL_INTERP_FROM_IDENTITY
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

    ! CWP_Status_t
    enum, bind(c)
        enumerator :: &
                CWP_STATUS_OFF, &
                CWP_STATUS_ON
    end enum

    ! CWP_Err_t
    enum, bind(c)
      enumerator :: &
        CWP_ERR_NO_ERROR, &
        CWP_ERR_DEFAULT
    end enum

    ! CWP_State_t
    enum, bind(c)
      enumerator :: &
        CWP_STATE_IN_PROGRESS, &
        CWP_STATE_END,         &
        CWP_STATE_OUTPUT_ERROR
    end enum

    ! CWP_Op_t
    enum, bind(c)
      enumerator :: &
        CWP_OP_MIN, &
        CWP_OP_MAX, &
        CWP_OP_SUM
    end enum

    ! CWP_PartData_exch_t
    enum, bind(c)
      enumerator :: &
        CWP_PARTDATA_SEND, &
        CWP_PARTDATA_RECV
    end enum

    interface CWP_Param_set; module procedure &
      CWP_Param_set_int_, &
      CWP_Param_set_double_, &
      CWP_Param_set_char_
    end interface CWP_Param_set

    interface CWP_Param_add; module procedure &
      CWP_Param_add_int_, &
      CWP_Param_add_double_, &
      CWP_Param_add_char_
    end interface CWP_Param_add

    interface CWP_Init ; module procedure &
        CWP_Init_
    end interface CWP_Init

    interface CWP_C_to_f_string ; module procedure &
        CWP_C_to_f_string_
    end interface CWP_C_to_f_string

    interface CWP_State_update ; module procedure &
        CWP_State_update_
    end interface CWP_State_update

    interface CWP_Time_step_beg ; module procedure &
        CWP_Time_step_beg_
    end interface CWP_Time_step_beg

    interface CWP_Time_step_end ; module procedure &
        CWP_Time_step_end_
    end interface CWP_Time_step_end

    interface CWP_User_structure_set ; module procedure &
        CWP_User_structure_set_
    end interface CWP_User_structure_set

    interface CWP_User_structure_get ; module procedure &
        CWP_User_structure_get_
    end interface CWP_User_structure_get

    interface CWP_Output_file_set ; module procedure &
        CWP_Output_file_set_
    end interface CWP_Output_file_set

    interface CWP_State_get ; module procedure &
        CWP_State_get_
    end interface CWP_State_get

    interface CWP_Cpl_create ; module procedure &
        CWP_Cpl_create_
    end interface CWP_Cpl_create

    interface CWP_Cpl_barrier ; module procedure &
        CWP_Cpl_barrier_
    end interface CWP_Cpl_barrier

    interface CWP_Cpl_Del ; module procedure &
        CWP_Cpl_Del_
    end interface CWP_Cpl_Del

    interface CWP_Computed_tgts_bcast_enable ; module procedure &
        CWP_Computed_tgts_bcast_enable_
    end interface CWP_Computed_tgts_bcast_enable

    interface CWP_Involved_srcs_bcast_enable ; module procedure &
        CWP_Involved_srcs_bcast_enable_
    end interface CWP_Involved_srcs_bcast_enable

    interface CWP_N_uncomputed_tgts_get ; module procedure &
        CWP_N_uncomputed_tgts_get_
    end interface CWP_N_uncomputed_tgts_get

    interface CWP_Uncomputed_tgts_get ; module procedure &
        CWP_Uncomputed_tgts_get_
    end interface CWP_Uncomputed_tgts_get

    interface CWP_N_computed_tgts_get ; module procedure &
        CWP_N_computed_tgts_get_
    end interface CWP_N_computed_tgts_get

    interface CWP_Computed_tgts_get ; module procedure &
        CWP_Computed_tgts_get_
    end interface CWP_Computed_tgts_get

    interface CWP_N_involved_srcs_get ; module procedure &
        CWP_N_involved_srcs_get_
    end interface CWP_N_involved_srcs_get

    interface CWP_Involved_srcs_get ; module procedure &
        CWP_Involved_srcs_get_
    end interface CWP_Involved_srcs_get

    interface CWP_Spatial_interp_weights_compute ; module procedure &
        CWP_Spatial_interp_weights_compute_
    end interface CWP_Spatial_interp_weights_compute

    interface CWP_Spatial_interp_property_set ; module procedure &
        CWP_Spatial_interp_property_set_
    end interface CWP_Spatial_interp_property_set

    interface CWP_Visu_set ; module procedure &
        CWP_Visu_set_
    end interface CWP_Visu_set

    interface CWP_User_tgt_pts_set ; module procedure &
        CWP_User_tgt_pts_set_
    end interface CWP_User_tgt_pts_set

    interface CWP_Mesh_interf_finalize ; module procedure &
        CWP_Mesh_interf_finalize_
    end interface CWP_Mesh_interf_finalize

    interface CWP_Mesh_interf_vtx_set ; module procedure &
        CWP_Mesh_interf_vtx_set_
    end interface CWP_Mesh_interf_vtx_set

    interface CWP_Mesh_interf_block_add ; module procedure &
        CWP_Mesh_interf_block_add_
    end interface CWP_Mesh_interf_block_add

    interface CWP_Mesh_interf_block_std_set ; module procedure &
        CWP_Mesh_interf_block_std_set_
    end interface CWP_Mesh_interf_block_std_set

    interface CWP_Mesh_interf_block_std_get ; module procedure &
        CWP_Mesh_interf_block_std_get_
    end interface CWP_Mesh_interf_block_std_get

    interface CWP_Mesh_interf_f_poly_block_set ; module procedure &
        CWP_Mesh_interf_f_poly_block_set_
    end interface CWP_Mesh_interf_f_poly_block_set

    interface CWP_Mesh_interf_f_poly_block_get ; module procedure &
        CWP_Mesh_interf_f_poly_block_get_
    end interface CWP_Mesh_interf_f_poly_block_get

    interface CWP_Mesh_interf_c_poly_block_set ; module procedure &
        CWP_Mesh_interf_c_poly_block_set_
    end interface CWP_Mesh_interf_c_poly_block_set

    interface CWP_Mesh_interf_c_poly_block_get ; module procedure &
        CWP_Mesh_interf_c_poly_block_get_
    end interface CWP_Mesh_interf_c_poly_block_get

    interface CWP_Mesh_interf_del ; module procedure &
        CWP_Mesh_interf_del_
    end interface CWP_Mesh_interf_del

    interface CWP_Mesh_interf_from_cellface_set ; module procedure &
        CWP_Mesh_interf_from_cellface_set_
    end interface CWP_Mesh_interf_from_cellface_set

    interface CWP_Mesh_interf_from_faceedge_set ; module procedure &
        CWP_Mesh_interf_from_faceedge_set_
    end interface CWP_Mesh_interf_from_faceedge_set

    interface CWP_Field_create ; module procedure &
        CWP_Field_create_
    end interface CWP_Field_create

    interface CWP_Field_data_set ; module procedure &
        CWP_Field_data_set_
    end interface CWP_Field_data_set

    interface CWP_Field_dof_location_get ; module procedure &
        CWP_Field_dof_location_get_
    end interface CWP_Field_dof_location_get

    interface CWP_Field_storage_get ; module procedure &
        CWP_Field_storage_get_
    end interface CWP_Field_storage_get

    interface CWP_Field_n_dof_get ; module procedure &
        CWP_Field_n_dof_get_
    end interface CWP_Field_n_dof_get

    interface CWP_Field_del ; module procedure &
        CWP_Field_del_
    end interface CWP_Field_del

    interface CWP_Field_issend ; module procedure &
        CWP_Field_issend_
    end interface CWP_Field_issend

    interface CWP_Field_irecv ; module procedure &
        CWP_Field_irecv_
    end interface CWP_Field_irecv

    interface CWP_Field_wait_issend ; module procedure &
        CWP_Field_wait_issend_
    end interface CWP_Field_wait_issend

    interface CWP_Field_wait_irecv ; module procedure &
        CWP_Field_wait_irecv_
    end interface CWP_Field_wait_irecv

    interface CWP_Field_interp_function_unset ; module procedure &
        CWP_Field_interp_function_unset_
    end interface CWP_Field_interp_function_unset

    interface CWP_Field_interp_function_set ; module procedure &
        CWP_Field_interp_function_set_
    end interface CWP_Field_interp_function_set

    interface CWP_Field_n_components_get ; module procedure &
        CWP_Field_n_components_get_
    end interface CWP_Field_n_components_get

    interface CWP_Field_src_data_properties_get ; module procedure &
        CWP_Field_src_data_properties_get_
    end interface CWP_Field_src_data_properties_get

    interface CWP_Field_tgt_data_properties_get ; module procedure &
        CWP_Field_tgt_data_properties_get_
    end interface CWP_Field_tgt_data_properties_get

    interface CWP_Field_location_weights_get ; module procedure &
        CWP_Field_location_weights_get_
    end interface CWP_Field_location_weights_get

    interface CWP_Field_location_point_data_get ; module procedure &
        CWP_Field_location_point_data_get_
    end interface CWP_Field_location_point_data_get

    interface CWP_Field_location_internal_cell_vtx_get ; module procedure &
        CWP_Field_location_internal_cell_vtx_get_
    end interface CWP_Field_location_internal_cell_vtx_get

    interface CWP_Field_intersection_volumes_get ; module procedure &
        CWP_Field_intersection_volumes_get_
    end interface CWP_Field_intersection_volumes_get

    interface CWP_Field_intersection_tgt_elt_volumes_get ; module procedure &
        CWP_Field_intersection_tgt_elt_volumes_get_
    end interface CWP_Field_intersection_tgt_elt_volumes_get

    interface CWP_Field_nearest_neighbors_distances_get ; module procedure &
        CWP_Field_nearest_neighbors_distances_get_
    end interface CWP_Field_nearest_neighbors_distances_get

    interface CWP_Field_nearest_neighbors_coord_get ; module procedure &
        CWP_Field_nearest_neighbors_coord_get_
    end interface CWP_Field_nearest_neighbors_coord_get

    interface CWP_Param_del ; module procedure &
        CWP_Param_del_
    end interface CWP_Param_del

    interface CWP_Param_n_get ; module procedure &
        CWP_Param_n_get_
    end interface CWP_Param_n_get

    interface CWP_Param_is ; module procedure &
        CWP_Param_is_
    end interface CWP_Param_is

    interface CWP_Param_get
      module procedure CWP_Param_get_int
      module procedure CWP_Param_get_double
      module procedure CWP_Param_get_char
    end interface CWP_Param_get

    interface CWP_Param_reduce
      module procedure CWP_Param_reduce_int
      module procedure CWP_Param_reduce_double
      module procedure CWP_Param_reduce_char
    end interface CWP_Param_reduce

    interface CWP_Param_lock ; module procedure &
        CWP_Param_lock_
    end interface CWP_Param_lock

    interface CWP_Param_unlock ; module procedure &
        CWP_Param_unlock_
    end interface CWP_Param_unlock

    interface CWP_Codes_list_get ; module procedure &
        CWP_Codes_list_get_
    end interface CWP_Codes_list_get

    interface CWP_Loc_codes_list_get ; module procedure &
        CWP_Loc_codes_list_get_
    end interface CWP_Loc_codes_list_get

    interface CWP_Param_list_get ; module procedure &
        CWP_Param_list_get_
    end interface CWP_Param_list_get

    interface CWP_Global_data_issend
      module procedure CWP_Global_data_issend_int
      module procedure CWP_Global_data_issend_long
      module procedure CWP_Global_data_issend_double
      module procedure CWP_Global_data_issend_complex4
      module procedure CWP_Global_data_issend_complex8
      module procedure CWP_Global_data_issend_real4
    end interface CWP_Global_data_issend

    interface CWP_Global_data_irecv
      module procedure CWP_Global_data_irecv_int
      module procedure CWP_Global_data_irecv_long
      module procedure CWP_Global_data_irecv_double
      module procedure CWP_Global_data_irecv_complex4
      module procedure CWP_Global_data_irecv_complex8
      module procedure CWP_Global_data_irecv_real4
    end interface CWP_Global_data_irecv

    interface CWP_Global_data_wait_issend
      module procedure CWP_Global_data_wait_issend_
    end interface CWP_Global_data_wait_issend

    interface CWP_Global_data_wait_irecv
      module procedure CWP_Global_data_wait_irecv_
    end interface CWP_Global_data_wait_irecv

    interface CWP_Part_data_create
      module procedure CWP_Part_data_create_
    end interface CWP_Part_data_create

    interface CWP_Part_data_del
      module procedure CWP_Part_data_del_
    end interface CWP_Part_data_del

    interface CWP_Part_data_issend
      module procedure CWP_Part_data_issend_
    end interface CWP_Part_data_issend

    interface CWP_Part_data_irecv
      module procedure CWP_Part_data_irecv_
    end interface CWP_Part_data_irecv

    interface CWP_Part_data_wait_issend
      module procedure CWP_Part_data_wait_issend_
    end interface CWP_Part_data_wait_issend

    interface CWP_Part_data_wait_irecv
      module procedure CWP_Part_data_wait_irecv_
    end interface CWP_Part_data_wait_irecv

    interface CWP_Cpl_spatial_interp_algo_get
      module procedure CWP_Cpl_spatial_interp_algo_get_
    end interface CWP_Cpl_spatial_interp_algo_get

  !
  ! Private

  private :: c_f_char_array,&
             CWP_Init_ ,&
             CWP_C_to_f_string_ ,&
             CWP_State_update_ ,&
             CWP_Time_step_beg_ ,&
             CWP_Time_step_end_ ,&
             CWP_User_structure_set_ ,&
             CWP_User_structure_get_ ,&
             CWP_Output_file_set_,&
             CWP_State_get_ ,&
             CWP_Cpl_create_ ,&
             CWP_Cpl_barrier_ ,&
             CWP_Cpl_Del_ ,&
             CWP_Computed_tgts_bcast_enable_,&
             CWP_Involved_srcs_bcast_enable_,&
             CWP_N_uncomputed_tgts_get_ ,&
             CWP_Uncomputed_tgts_get_ ,&
             CWP_N_computed_tgts_get_ ,&
             CWP_Computed_tgts_get_ ,&
             CWP_N_involved_srcs_get_ ,&
             CWP_Involved_srcs_get_ ,&
             CWP_Spatial_interp_weights_compute_ ,&
             CWP_Spatial_interp_property_set_ ,&
             CWP_Visu_set_ ,&
             CWP_User_tgt_pts_set_ ,&
             CWP_Mesh_interf_finalize_ ,&
             CWP_Mesh_interf_vtx_set_ ,&
             CWP_Mesh_interf_block_add_ ,&
             CWP_Mesh_interf_block_std_set_ ,&
             CWP_Mesh_interf_block_std_get_ ,&
             CWP_Mesh_interf_f_poly_block_set_ ,&
             CWP_Mesh_interf_f_poly_block_get_ ,&
             CWP_Mesh_interf_c_poly_block_set_ ,&
             CWP_Mesh_interf_c_poly_block_get_ ,&
             CWP_Mesh_interf_del_ ,&
             CWP_Mesh_interf_from_cellface_set_ ,&
             CWP_Mesh_interf_from_faceedge_set_ ,&
             CWP_Field_create_ ,&
             CWP_Field_data_set_ ,&
             CWP_Field_dof_location_get_ ,&
             CWP_Field_storage_get_ ,&
             CWP_Field_n_dof_get_, &
             CWP_Field_del_ ,&
             CWP_Field_issend_ ,&
             CWP_Field_irecv_ ,&
             CWP_Field_wait_issend_ ,&
             CWP_Field_wait_irecv_ ,&
             CWP_Field_interp_function_unset_ ,&
             CWP_Field_interp_function_set_ ,&
             CWP_Field_n_components_get_, &
             CWP_Field_src_data_properties_get_, &
             CWP_Field_tgt_data_properties_get_, &
             CWP_Field_location_weights_get_, &
             CWP_Field_location_point_data_get_, &
             CWP_Field_location_internal_cell_vtx_get_, &
             CWP_Field_intersection_volumes_get_, &
             CWP_Field_intersection_tgt_elt_volumes_get_, &
             CWP_Field_nearest_neighbors_distances_get_, &
             CWP_Field_nearest_neighbors_coord_get_, &
             CWP_Param_add_int_ ,&
             CWP_Param_add_double_ ,&
             CWP_Param_add_char_ ,&
             CWP_Param_set_int_ ,&
             CWP_Param_set_double_ ,&
             CWP_Param_set_char_ ,&
             CWP_Param_del_ ,&
             CWP_Param_n_get_ ,&
             CWP_Param_is_ ,&
             CWP_Param_get_int ,&
             CWP_Param_get_double ,&
             CWP_Param_get_char, &
             CWP_Param_reduce_int,&
             CWP_Param_reduce_double,&
             CWP_Param_reduce_char, &
             CWP_Param_lock_ ,&
             CWP_Param_unlock_,&
             CWP_Codes_list_get_,&
             CWP_Loc_codes_list_get_,&
             CWP_Param_list_get_, &
             CWP_Global_data_issend_int, &
             CWP_Global_data_issend_long, &
             CWP_Global_data_issend_double, &
             CWP_Global_data_issend_complex4, &
             CWP_Global_data_issend_complex8, &
             CWP_Global_data_issend_real4, &
             CWP_Global_data_irecv_int, &
             CWP_Global_data_irecv_long, &
             CWP_Global_data_irecv_double, &
             CWP_Global_data_irecv_complex4, &
             CWP_Global_data_irecv_complex8, &
             CWP_Global_data_irecv_real4, &
             CWP_Global_data_wait_issend_, &
             CWP_Global_data_wait_irecv_, &
             CWP_Part_data_create_, &
             CWP_Part_data_del_, &
             CWP_Part_data_issend_, &
             CWP_Part_data_irecv_, &
             CWP_Part_data_wait_issend_, &
             CWP_Part_data_wait_irecv_, &
             CWP_Cpl_spatial_interp_algo_get_

    interface
      subroutine CWP_Init_cf(fcomm, n_code, code_names, l_code_names, is_active_rank, intra_comms) &
              bind(c, name = 'CWP_Init_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int), value :: fcomm
        integer(c_int), value :: n_code
        type(c_ptr),    value :: code_names
        type(c_ptr),    value :: l_code_names
        integer(c_int), value :: is_active_rank
        type(c_ptr),    value :: intra_comms
      end subroutine CWP_Init_cf

      subroutine CWP_State_update_cf(local_code_name, l_local_code_name, state) &
        bind(c, name='CWP_State_update_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name
        integer(c_int), value             :: l_local_code_name, state
      end subroutine CWP_State_update_cf

      subroutine CWP_Time_step_beg_cf(local_code_name, l_local_code_name, current_time) &
        bind(c, name='CWP_Time_step_beg_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name
        integer(c_int), value             :: l_local_code_name
        real(c_double), value             :: current_time
      end subroutine CWP_Time_step_beg_cf

      subroutine CWP_Time_step_end_cf(local_code_name, l_local_code_name) &
        bind(c, name='CWP_Time_step_end_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name
        integer(c_int), value             :: l_local_code_name
      end subroutine CWP_Time_step_end_cf

      subroutine CWP_User_structure_set_cf(local_code_name, l_local_code_name, user_structure) &
        bind (c, name="CWP_User_structure_set_cf")
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name
        integer(c_int), value             :: l_local_code_name
        type(c_ptr),    value             :: user_structure
      end subroutine CWP_User_structure_set_cf

      function CWP_User_structure_get_cf(local_code_name, l_local_code_name) &
        result (user_structure)                                              &
        bind (c, name="CWP_User_structure_get_cf")
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name
        integer(c_int), value             :: l_local_code_name
        type(c_ptr)                       :: user_structure
      end function CWP_User_structure_get_cf

      subroutine CWP_Output_file_set_cf(f_output_file_name, l_output_file_name) &
        bind (c, name="CWP_Output_file_set_cf")
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: f_output_file_name
        integer(c_int), value             :: l_output_file_name
      end subroutine CWP_Output_file_set_cf

      function CWP_State_get_cf(local_code_name, l_local_code_name) &
        result (state)                                              &
        bind(c, name='CWP_State_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name
        integer(c_int), value             :: l_local_code_name
        integer(c_int)                    :: state
      end function CWP_State_get_cf

      subroutine CWP_Cpl_create_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, coupled_code_name, l_coupled_code_name, &
        entities_dim, comm_type, spatial_interp, n_part, displacement, freq) &
        bind(c, name = 'CWP_Cpl_create_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, coupled_code_name
        integer(c_int), value :: l_local_code_name, l_cpl_id, l_coupled_code_name
        integer(c_int), value :: entities_dim
        integer(c_int), value :: comm_type
        integer(c_int), value :: spatial_interp
        integer(c_int), value :: n_part
        integer(c_int), value :: displacement
        integer(c_int), value :: freq
      end subroutine CWP_Cpl_Create_cf

      subroutine CWP_Cpl_barrier_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id) &
        bind(c, name = 'CWP_Cpl_barrier_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value             :: l_local_code_name, l_cpl_id
      end subroutine CWP_Cpl_barrier_cf

      subroutine CWP_Cpl_del_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id) &
        bind(c, name = 'CWP_Cpl_del_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
      end subroutine CWP_Cpl_del_cf

      subroutine CWP_Computed_tgts_bcast_enable_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id) &
        bind(c, name = 'CWP_Computed_tgts_bcast_enable_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      end subroutine CWP_Computed_tgts_bcast_enable_cf

      subroutine CWP_Involved_srcs_bcast_enable_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id) &
        bind(c, name = 'CWP_Involved_srcs_bcast_enable_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      end subroutine CWP_Involved_srcs_bcast_enable_cf

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
        bind(c, name = 'CWP_Mesh_interf_f_poly_block_set_cf')
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
            face_edge_idx, face_edge, n_edges, edge_vtx, parent_num) &
            bind(c, name = 'CWP_Mesh_interf_from_faceedge_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: i_part, n_faces, n_edges
        type(c_ptr), value :: face_edge_idx, face_edge, edge_vtx, parent_num
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

      subroutine CWP_Field_issend_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id) &
            bind(c, name = 'CWP_Field_issend_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      end subroutine CWP_Field_issend_cf

      subroutine CWP_Field_irecv_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, tgt_field_id, l_tgt_field_id) &
            bind(c, name = 'CWP_Field_irecv_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, tgt_field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_tgt_field_id
      end subroutine CWP_Field_irecv_cf

      subroutine CWP_Field_wait_issend_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id) &
            bind(c, name = 'CWP_Field_wait_issend_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      end subroutine CWP_Field_wait_issend_cf

      subroutine CWP_Field_wait_irecv_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, tgt_field_id, l_tgt_field_id) &
            bind(c, name = 'CWP_Field_wait_irecv_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, tgt_field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_tgt_field_id
      end subroutine CWP_Field_wait_irecv_cf

      subroutine CWP_Field_interp_function_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
              field_id, l_field_id, user_interpolation_fct) &
              bind(c, name = 'CWP_Field_interp_function_set_cf')
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = 1)       :: local_code_name, cpl_id, field_id
        integer(kind = c_int),            value :: l_local_code_name, l_cpl_id, l_field_id
        type(c_funptr),                   value :: user_interpolation_fct
      end subroutine CWP_Field_interp_function_set_cf

      subroutine CWP_Field_interp_function_unset_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
              field_id, l_field_id) &
              bind(c, name = 'CWP_Field_interp_function_unset_cf')
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      end subroutine CWP_Field_interp_function_unset_cf

      function CWP_Field_n_components_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
                                             field_id, l_field_id) &
        result (n_components)                                                          &
        bind(c, name='CWP_Field_n_components_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
        integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_field_id
        integer(c_int)                    :: n_components
      end function CWP_Field_n_components_get_cf

      subroutine CWP_Field_src_data_properties_get_cf(local_code_name,   &
                                                      l_local_code_name, &
                                                      cpl_id,            &
                                                      l_cpl_id,          &
                                                      field_id,          &
                                                      l_field_id,        &
                                                      i_part,            &
                                                      n_elt_src,         &
                                                      c_src_to_tgt_idx)  &
          bind(c, name = 'CWP_Field_src_data_properties_get_cf')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
          integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_field_id
          integer(c_int), value             :: i_part
          integer(c_int)                    :: n_elt_src
          type(c_ptr)                       :: c_src_to_tgt_idx
      end subroutine CWP_Field_src_data_properties_get_cf

      subroutine CWP_Field_tgt_data_properties_get_cf(local_code_name,         &
                                                      l_local_code_name,       &
                                                      cpl_id,                  &
                                                      l_cpl_id,                &
                                                      field_id,                &
                                                      l_field_id,              &
                                                      i_part,                  &
                                                      n_elt_tgt,               &
                                                      n_referenced_tgt,        &
                                                      c_referenced_tgt,        &
                                                      c_tgt_come_from_src_idx) &
          bind(c, name = 'CWP_Field_tgt_data_properties_get_cf')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
          integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_field_id
          integer(c_int), value             :: i_part
          integer(c_int)                    :: n_elt_tgt, n_referenced_tgt
          type(c_ptr)                       :: c_referenced_tgt, c_tgt_come_from_src_idx
      end subroutine CWP_Field_tgt_data_properties_get_cf

      subroutine CWP_Field_location_weights_get_cf(local_code_name,          &
                                                   l_local_code_name,        &
                                                   cpl_id,                   &
                                                   l_cpl_id,                 &
                                                   field_id,                 &
                                                   l_field_id,               &
                                                   i_part,                   &
                                                   c_weights,                &
                                                   s_weights)                & !! TO DO: get with coupling
          bind(c, name = 'CWP_Field_location_weights_get_cf')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
          integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_field_id
          integer(c_int), value             :: i_part
          type(c_ptr)                       :: c_weights
          integer(kind = c_int)             :: s_weights
      end subroutine CWP_Field_location_weights_get_cf

      subroutine CWP_Field_location_point_data_get_cf(local_code_name,           &
                                                       l_local_code_name,        &
                                                       cpl_id,                   &
                                                       l_cpl_id,                 &
                                                       field_id,                 &
                                                       l_field_id,               &
                                                       i_part,                   &
                                                       c_points_coords,          &
                                                       c_points_uvw,             &
                                                       c_points_dist2,           &
                                                       c_points_projected_coords,&
                                                       s_size)&

          bind(c, name = 'CWP_Field_location_point_data_get_cf')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
          integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_field_id
          integer(c_int), value             :: i_part
          type(c_ptr)                       :: c_points_coords, c_points_uvw, c_points_dist2, c_points_projected_coords
          integer(kind = c_int) :: s_size
      end subroutine CWP_Field_location_point_data_get_cf

      subroutine CWP_Field_intersection_volumes_get_cf(local_code_name,      &
                                                        l_local_code_name,    &
                                                        cpl_id,               &
                                                        l_cpl_id,             &
                                                        field_id,         &
                                                        l_field_id,       &
                                                        i_part,               &
                                                        c_volumes,            &
                                                        s_volumes)            &
          bind(c, name = 'CWP_Field_intersection_volumes_get_cf')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
          integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_field_id
          integer(c_int), value             :: i_part
          type(c_ptr)                       :: c_volumes
          integer(kind = c_int)             :: s_volumes
      end subroutine CWP_Field_intersection_volumes_get_cf

      subroutine CWP_Field_intersection_tgt_elt_volumes_get_cf(local_code_name,      &
                                                                l_local_code_name,    &
                                                                cpl_id,               &
                                                                l_cpl_id,             &
                                                                field_id,         &
                                                                l_field_id,       &
                                                                i_part,               &
                                                                c_tgt_elt_volumes,    &
                                                                n_elt)                &
          bind(c, name = 'CWP_Field_intersection_tgt_elt_volumes_get_cf')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind=c_char, len = 1) :: local_code_name, cpl_id, field_id
          integer(c_int), value           :: l_local_code_name, l_cpl_id, l_field_id
          integer(c_int), value           :: i_part
          type(c_ptr)                     :: c_tgt_elt_volumes
          integer(c_int)                  :: n_elt
      end subroutine CWP_Field_intersection_tgt_elt_volumes_get_cf

      subroutine CWP_Field_nearest_neighbors_distances_get_cf(local_code_name,   &
                                                               l_local_code_name, &
                                                               cpl_id,            &
                                                               l_cpl_id,          &
                                                               field_id,          &
                                                               l_field_id,        &
                                                               i_part,            &
                                                               c_distances2,      &
                                                               s_distances2)      &
          bind(c, name = 'CWP_Field_nearest_neighbors_distances_get_cf')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
          integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_field_id
          integer(c_int), value             :: i_part
          type(c_ptr)                       :: c_distances2
          integer(kind = c_int)             :: s_distances2
      end subroutine CWP_Field_nearest_neighbors_distances_get_cf

      subroutine CWP_Field_nearest_neighbors_coord_get_cf(f_local_code_name,   &
                                                           l_local_code_name,   &
                                                           f_cpl_id,            &
                                                           l_cpl_id,            &
                                                           f_field_id,          &
                                                           l_field_id,          &
                                                           i_part,              &
                                                           c_nearest_src_coord, &
                                                           n_nearest_src_pts)   &
      bind (c, name='CWP_Field_nearest_neighbors_coord_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind=c_char, len = 1) :: f_local_code_name, f_cpl_id, f_field_id
        integer(c_int), value           :: l_local_code_name, l_cpl_id, l_field_id
        integer(c_int), value           :: i_part
        type(c_ptr)                     :: c_nearest_src_coord
        integer(c_int)                  :: n_nearest_src_pts
      end subroutine CWP_Field_nearest_neighbors_coord_get_cf

      subroutine CWP_Field_location_internal_cell_vtx_get_cf(local_code_name,   &
                                                              l_local_code_name, &
                                                              cpl_id,            &
                                                              l_cpl_id,          &
                                                              field_id,          &
                                                              l_field_id,        &
                                                              i_part,            &
                                                              cell_vtx_idx,      &
                                                              n_cell,            &
                                                              cell_vtx)          &

          bind(c, name = 'CWP_Field_location_internal_cell_vtx_get_cf')
          use, intrinsic :: iso_c_binding
          implicit none
          character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
          integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_field_id
          integer(c_int), value             :: i_part
          type(c_ptr)                       :: cell_vtx_idx, cell_vtx
          integer(kind = c_int) :: n_cell
      end subroutine CWP_Field_location_internal_cell_vtx_get_cf

      subroutine CWP_Spatial_interp_property_set_cf(local_code_name,   &
                                                    l_local_code_name, &
                                                    cpl_id,            &
                                                    l_cpl_id,          &
                                                    property_name,     &
                                                    l_property_name,   &
                                                    property_type,     &
                                                    property_value,    &
                                                    l_property_value)  &
      bind(c, name = 'CWP_Spatial_interp_property_set_cf')
        use, intrinsic :: iso_c_binding
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, property_name, property_value
        integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_property_name, property_type, l_property_value
      end subroutine CWP_Spatial_interp_property_set_cf


    !>
    !!
    !! \brief Finalize CWIPI.
    !!
    !!

    subroutine CWP_Finalize() &
          bind(c, name = 'CWP_Finalize')
    end subroutine CWP_Finalize

    !>
    !! \brief Return the number of codes known by CWIPI.
    !!
    !! \return Number of codes
    !!
    !!

    function CWP_Codes_nb_get() &
      result (n_codes)          &
      bind(c, name='CWP_Codes_nb_get')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: n_codes
    end function CWP_Codes_nb_get


    !>
    !! \brief Return the number of local codes known by CWIPI.
    !!
    !! \return Number of local codes
    !!
    !!

    function CWP_Loc_codes_nb_get() &
      result (n_local_codes)        &
      bind(c, name='CWP_Loc_codes_nb_get')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: n_local_codes
    end function CWP_Loc_codes_nb_get


    !>
    !! \brief Dump code properties.
    !!
    !!

    subroutine CWP_Properties_dump() &
      bind(c, name='CWP_Properties_dump')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine CWP_Properties_dump


    subroutine CWP_Param_add_int_cf(local_code_name,   &
                                    l_local_code_name, &
                                    param_name,        &
                                    l_param_name,      &
                                    initial_value)     &
      bind(c, name='CWP_Param_add_int_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      integer(c_int), value :: initial_value
    end subroutine CWP_Param_add_int_cf


    subroutine CWP_Param_add_double_cf(local_code_name,   &
                                       l_local_code_name, &
                                       param_name,        &
                                       l_param_name,      &
                                       initial_value)     &
      bind(c, name='CWP_Param_add_double_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      real(c_double), value :: initial_value
    end subroutine CWP_Param_add_double_cf


    subroutine CWP_Param_add_char_cf(local_code_name,   &
                                     l_local_code_name, &
                                     param_name,        &
                                     l_param_name,      &
                                     initial_value,     &
                                     l_initial_value)   &
      bind(c, name='CWP_Param_add_char_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name, l_initial_value
      character(kind = c_char, len = 1) :: initial_value
    end subroutine CWP_Param_add_char_cf

    subroutine CWP_Param_set_int_cf(local_code_name,   &
                                    l_local_code_name, &
                                    param_name,        &
                                    l_param_name,      &
                                    initial_value)     &
      bind(c, name='CWP_Param_set_int_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      integer(c_int), value :: initial_value
    end subroutine CWP_Param_set_int_cf


    subroutine CWP_Param_set_double_cf(local_code_name,   &
                                       l_local_code_name, &
                                       param_name,        &
                                       l_param_name,      &
                                       initial_value)     &
      bind(c, name='CWP_Param_set_double_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      real(c_double), value :: initial_value
    end subroutine CWP_Param_set_double_cf


    subroutine CWP_Param_set_char_cf(local_code_name,   &
                                     l_local_code_name, &
                                     param_name,        &
                                     l_param_name,      &
                                     initial_value,     &
                                     l_initial_value)   &
      bind(c, name='CWP_Param_set_char_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name, l_initial_value
      character(kind = c_char, len = 1) :: initial_value
    end subroutine CWP_Param_set_char_cf


    subroutine CWP_Param_del_cf(local_code_name,   &
                                l_local_code_name, &
                                param_name,        &
                                l_param_name,      &
                                data_type)         &
      bind(c, name='CWP_Param_del_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      integer(c_int), value :: data_type
    end subroutine CWP_Param_del_cf


    function CWP_Param_n_get_cf(local_code_name,   &
                                l_local_code_name, &
                                data_type)         &
      result (n_param)                             &
      bind(c, name='CWP_Param_n_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name
      integer(c_int), value :: l_local_code_name
      integer(c_int), value :: data_type
      integer(c_int)        :: n_param
    end function CWP_Param_n_get_cf


    function CWP_Param_is_cf(local_code_name,   &
                             l_local_code_name, &
                             param_name,        &
                             l_param_name,      &
                             data_type)         &
      result (is_param)                         &
      bind(c, name='CWP_Param_is_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      integer(c_int), value :: data_type
      integer(c_int)        :: is_param
    end function CWP_Param_is_cf


    subroutine CWP_Param_get_int_cf(local_code_name,   &
                                    l_local_code_name, &
                                    param_name,        &
                                    l_param_name,      &
                                    value)             &
      bind(c, name='CWP_Param_get_int_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      integer(c_int)        :: value
    end subroutine CWP_Param_get_int_cf

    subroutine CWP_Param_get_double_cf(local_code_name,   &
                                       l_local_code_name, &
                                       param_name,        &
                                       l_param_name,      &
                                       value)             &
      bind(c, name='CWP_Param_get_double_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      real(c_double)        :: value
    end subroutine CWP_Param_get_double_cf

    subroutine CWP_Param_get_char_cf(local_code_name,   &
                                     l_local_code_name, &
                                     param_name,        &
                                     l_param_name,      &
                                     val,               &
                                     l_value)           &
      bind(c, name='CWP_Param_get_char_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value             :: l_local_code_name, l_param_name
      type(c_ptr)                       :: val
      integer(c_int)                    :: l_value
    end subroutine CWP_Param_get_char_cf


    subroutine CWP_Param_reduce_int_cf(op,           &
                                        param_name,   &
                                        l_param_name, &
                                        res,          &
                                        n_codes,      &
                                        code_names,   &
                                        l_code_names) &

      bind(c, name='CWP_Param_reduce_int_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value             :: op, n_codes
      integer(c_int)                    :: res
      character(kind = c_char, len = 1) :: param_name
      integer(c_int), value             :: l_param_name
      type(c_ptr),    value             :: code_names
      type(c_ptr),    value             :: l_code_names
    end subroutine CWP_Param_reduce_int_cf


    subroutine CWP_Param_reduce_double_cf(op,           &
                                          param_name,   &
                                          l_param_name, &
                                          res,          &
                                          n_codes,      &
                                          code_names,   &
                                          l_code_names) &

      bind(c, name='CWP_Param_reduce_double_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value             :: op, n_codes
      real(c_double)                    :: res
      character(kind = c_char, len = 1) :: param_name
      integer(c_int), value             :: l_param_name
      type(c_ptr),    value             :: code_names
      type(c_ptr),    value             :: l_code_names
    end subroutine CWP_Param_reduce_double_cf

    subroutine CWP_Param_reduce_char_cf(op,           &
                                        param_name,   &
                                        l_param_name, &
                                        res,          &
                                        l_res,        &
                                        n_codes,      &
                                        code_names,   &
                                        l_code_names) &

      bind(c, name='CWP_Param_reduce_char_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value             :: op, n_codes
      type(c_ptr)                       :: res
      integer(c_int)                    :: l_res
      character(kind = c_char, len = 1) :: param_name
      integer(c_int), value             :: l_param_name
      type(c_ptr),    value             :: code_names
      type(c_ptr),    value             :: l_code_names
    end subroutine CWP_Param_reduce_char_cf


    subroutine CWP_Param_lock_cf(local_code_name,   &
                                 l_local_code_name) &
      bind(c, name='CWP_Param_lock_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name
      integer(c_int), value :: l_local_code_name
    end subroutine CWP_Param_lock_cf


    subroutine CWP_Param_unlock_cf(local_code_name,   &
                                   l_local_code_name) &
      bind(c, name='CWP_Param_unlock_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name
      integer(c_int), value :: l_local_code_name
    end subroutine CWP_Param_unlock_cf

    subroutine CWP_Mesh_interf_block_std_get_cf(local_code_name,   &
                                                l_local_code_name, &
                                                cpl_id,            &
                                                l_cpl_id,          &
                                                i_part,            &
                                                block_id,          &
                                                n_elts,            &
                                                c_connec,          &
                                                c_global_num,      &
                                                s_connec)          &
        bind(c, name = 'CWP_Mesh_interf_block_std_get_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name, cpl_id
        integer(c_int), value :: i_part, block_id, n_elts
        type(c_ptr)      :: c_connec, c_global_num
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id
        integer(kind = c_int) :: s_connec
    end subroutine CWP_Mesh_interf_block_std_get_cf


    subroutine CWP_Mesh_interf_f_poly_block_get_cf(local_code_name,   &
                                                   l_local_code_name, &
                                                   cpl_id,            &
                                                   l_cpl_id,          &
                                                   i_part,            &
                                                   block_id,          &
                                                   n_elts,            &
                                                   connec_idx,        &
                                                   connec,            &
                                                   global_num)        &
      bind(c, name='CWP_Mesh_interf_f_poly_block_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, cpl_id
      integer(c_int), value :: l_local_code_name, l_cpl_id
      integer(c_int), value :: i_part, block_id
      integer(c_int)        :: n_elts
      type(c_ptr)           :: connec_idx, connec, global_num
    end subroutine CWP_Mesh_interf_f_poly_block_get_cf


    subroutine CWP_Mesh_interf_c_poly_block_get_cf(local_code_name,   &
                                                   l_local_code_name, &
                                                   cpl_id,            &
                                                   l_cpl_id,          &
                                                   i_part,            &
                                                   block_id,          &
                                                   n_elts,            &
                                                   n_faces,           &
                                                   connec_faces_idx,  &
                                                   connec_faces,      &
                                                   connec_cells_idx,  &
                                                   connec_cells,      &
                                                   global_num)        &
      bind(c, name='CWP_Mesh_interf_c_poly_block_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, cpl_id
      integer(c_int), value :: l_local_code_name, l_cpl_id
      integer(c_int), value :: i_part, block_id
      integer(c_int)        :: n_elts, n_faces
      type(c_ptr)           :: connec_faces_idx, connec_faces, connec_cells_idx, connec_cells, global_num
    end subroutine CWP_Mesh_interf_c_poly_block_get_cf


    function CWP_Field_dof_location_get_cf(local_code_name,   &
                                           l_local_code_name, &
                                           cpl_id,            &
                                           l_cpl_id,          &
                                           field_id,          &
                                           l_field_id)        &
      result (dof_location)                                   &
      bind(c, name='CWP_Field_dof_location_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
      integer(c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      integer(c_int)        :: dof_location
    end function CWP_Field_dof_location_get_cf


    function CWP_Field_storage_get_cf(local_code_name,   &
                                      l_local_code_name, &
                                      cpl_id,            &
                                      l_cpl_id,          &
                                      field_id,          &
                                      l_field_id)        &
      result (storage)                                   &
      bind(c, name='CWP_Field_storage_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
      integer(c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      integer(c_int)        :: storage
    end function CWP_Field_storage_get_cf

    function CWP_Field_n_dof_get_cf(local_code_name,   &
                                    l_local_code_name, &
                                    cpl_id,            &
                                    l_cpl_id,          &
                                    field_id,          &
                                    l_field_id,        &
                                    i_part)            &
      result (n_dof)                                     &
      bind(c, name='CWP_Field_n_dof_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
      integer(c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      integer(c_int), value :: i_part
      integer(c_int)        :: n_dof
    end function CWP_Field_n_dof_get_cf

    subroutine CWP_Field_del_cf(local_code_name,   &
                                l_local_code_name, &
                                cpl_id,            &
                                l_cpl_id,          &
                                field_id,          &
                                l_field_id)        &
      bind(c, name='CWP_Field_del_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
      integer(c_int), value :: l_local_code_name, l_cpl_id, l_field_id
    end subroutine CWP_Field_del_cf

    subroutine CWP_Codes_list_get_cf(code_list, code_list_s, n_codes) &
      bind (c, name="CWP_Codes_list_get_cf")
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            :: code_list, code_list_s
      integer(c_int)         :: n_codes
    end subroutine CWP_Codes_list_get_cf

    subroutine CWP_Loc_codes_list_get_cf(loc_code_list, loc_code_list_s, n_loc_codes) &
      bind (c, name="CWP_Loc_codes_list_get_cf")
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            :: loc_code_list, loc_code_list_s
      integer(c_int)         :: n_loc_codes
    end subroutine CWP_Loc_codes_list_get_cf

    subroutine CWP_Param_list_get_cf(code_name, l_code_name, data_type, n_param, c_param_names, c_param_sizes) &
      bind (c, name="CWP_Param_list_get_cf")
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)                       :: c_param_names, c_param_sizes
      integer(c_int)                    :: n_param
      character(kind = c_char, len = 1) :: code_name
      integer(c_int), value             :: l_code_name
      integer(c_int), value             :: data_type
    end subroutine CWP_Param_list_get_cf

    subroutine CWP_Global_data_issend_cf(f_local_code_name, &
                                         l_local_code_name, &
                                         f_cpl_id,          &
                                         l_cpl_id,          &
                                         f_global_data_id,  &
                                         l_global_data_id,  &
                                         s_send_entity,     &
                                         send_stride,       &
                                         n_send_entity,     &
                                         send_data)         &
    bind (c, name='CWP_Global_data_issend_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_global_data_id
      integer(c_int),  value        :: l_local_code_name, l_cpl_id, l_global_data_id
      integer(c_long), value        :: s_send_entity
      integer(c_int),  value        :: send_stride
      integer(c_int),  value        :: n_send_entity
      type(c_ptr),     value        :: send_data
    end subroutine CWP_Global_data_issend_cf

    subroutine CWP_Global_data_irecv_cf(f_local_code_name,      &
                                        l_local_code_name,      &
                                        f_cpl_id,               &
                                        l_cpl_id,               &
                                        f_global_data_id,       &
                                        l_global_data_id,       &
                                        s_recv_entity,          &
                                        recv_stride,            &
                                        n_recv_entity,          &
                                        recv_data)              &
    bind (c, name='CWP_Global_data_irecv_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_global_data_id
      integer(c_int),  value        :: l_local_code_name, l_cpl_id, l_global_data_id
      integer(c_long), value        :: s_recv_entity
      integer(c_int),  value        :: recv_stride
      integer(c_int),  value        :: n_recv_entity
      type(c_ptr),     value        :: recv_data
    end subroutine CWP_Global_data_irecv_cf

    subroutine CWP_Global_data_wait_issend_cf(f_local_code_name, &
                                              l_local_code_name, &
                                              f_cpl_id,          &
                                              l_cpl_id,          &
                                              f_global_data_id,  &
                                              l_global_data_id)  &
    bind (c, name='CWP_Global_data_wait_issend_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_global_data_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id, l_global_data_id
    end subroutine CWP_Global_data_wait_issend_cf

    subroutine CWP_Global_data_wait_irecv_cf(f_local_code_name, &
                                             l_local_code_name, &
                                             f_cpl_id,          &
                                             l_cpl_id,          &
                                             f_global_data_id,  &
                                             l_global_data_id)  &
    bind (c, name='CWP_Global_data_wait_irecv_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_global_data_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id, l_global_data_id
    end subroutine CWP_Global_data_wait_irecv_cf

    subroutine CWP_Part_data_create_cf(f_local_code_name, &
                                       l_local_code_name, &
                                       f_cpl_id,          &
                                       l_cpl_id,          &
                                       f_part_data_id,    &
                                       l_part_data_id,    &
                                       exch_type,         &
                                       gnum_elt,          &
                                       n_elt,             &
                                       n_part)            &
    bind (c, name='CWP_Part_data_create_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_part_data_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id, l_part_data_id
      integer(c_int), value         :: exch_type
      type(c_ptr),    value         :: gnum_elt
      type(c_ptr),    value         :: n_elt
      integer(c_int), value         :: n_part
    end subroutine CWP_Part_data_create_cf


    subroutine CWP_Part_data_del_cf(f_local_code_name, &
                                    l_local_code_name, &
                                    f_cpl_id,          &
                                    l_cpl_id,          &
                                    f_part_data_id,    &
                                    l_part_data_id,    &
                                    exch_type)         &
    bind (c, name='CWP_Part_data_del_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_part_data_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id, l_part_data_id
      integer(c_int), value         :: exch_type
    end subroutine CWP_Part_data_del_cf

    subroutine CWP_Part_data_issend_cf(f_local_code_name,   &
                                       l_local_code_name,   &
                                       f_cpl_id,            &
                                       l_cpl_id,            &
                                       f_part_data_id,      &
                                       l_part_data_id,      &
                                       s_data,              &
                                       n_components,        &
                                       part1_to_part2_data) &
    bind (c, name='CWP_Part_data_issend_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_part_data_id
      integer(c_int),  value        :: l_local_code_name, l_cpl_id, l_part_data_id
      integer(c_long), value        :: s_data
      integer(c_int),  value        :: n_components
      type(c_ptr),     value        :: part1_to_part2_data
    end subroutine CWP_Part_data_issend_cf

    subroutine CWP_Part_data_irecv_cf(f_local_code_name,   &
                                      l_local_code_name,   &
                                      f_cpl_id,            &
                                      l_cpl_id,            &
                                      f_part_data_id,      &
                                      l_part_data_id,      &
                                      s_data,              &
                                      n_components,        &
                                      part2_data)          &
    bind (c, name='CWP_Part_data_irecv_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_part_data_id
      integer(c_int),  value        :: l_local_code_name, l_cpl_id, l_part_data_id
      integer(c_long), value        :: s_data
      integer(c_int),  value        :: n_components
      type(c_ptr),     value        :: part2_data
    end subroutine CWP_Part_data_irecv_cf

    subroutine CWP_Part_data_wait_issend_cf(f_local_code_name,   &
                                            l_local_code_name,   &
                                            f_cpl_id,            &
                                            l_cpl_id,            &
                                            f_part_data_id,      &
                                            l_part_data_id)      &
    bind (c, name='CWP_Part_data_wait_issend_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_part_data_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id, l_part_data_id
    end subroutine CWP_Part_data_wait_issend_cf

    subroutine CWP_Part_data_wait_irecv_cf(f_local_code_name,   &
                                           l_local_code_name,   &
                                           f_cpl_id,            &
                                           l_cpl_id,            &
                                           f_part_data_id,      &
                                           l_part_data_id)      &
    bind (c, name='CWP_Part_data_wait_irecv_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_part_data_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id, l_part_data_id
    end subroutine CWP_Part_data_wait_irecv_cf

    subroutine CWP_Part_data_n_part_get_cf(f_local_code_name, &
                                           l_local_code_name, &
                                           f_cpl_id,          &
                                           l_cpl_id,          &
                                           f_part_data_id,    &
                                           l_part_data_id,    &
                                           n_part)            &
    bind (c, name='CWP_Part_data_n_part_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_part_data_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id, l_part_data_id
      integer(c_int)                :: n_part
    end subroutine CWP_Part_data_n_part_get_cf

    subroutine CWP_Part_data_n_ref_get_cf(f_local_code_name, &
                                          l_local_code_name, &
                                          f_cpl_id,          &
                                          l_cpl_id,          &
                                          f_part_data_id,    &
                                          l_part_data_id,    &
                                          i_part,            &
                                          n_ref)             &
    bind (c, name='CWP_Part_data_n_ref_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id, f_part_data_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id, l_part_data_id
      integer(c_int), value         :: i_part
      integer(c_int)                :: n_ref
    end subroutine CWP_Part_data_n_ref_get_cf


    function CWP_Cpl_spatial_interp_algo_get_cf(f_local_code_name,      &
                                                l_local_code_name,      &
                                                f_cpl_id,               &
                                                l_cpl_id) result (algo) &
    bind (c, name='CWP_Cpl_spatial_interp_algo_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1) :: f_local_code_name, f_cpl_id
      integer(c_int), value         :: l_local_code_name, l_cpl_id
      integer(c_int)                :: algo
    end function CWP_Cpl_spatial_interp_algo_get_cf

  end interface

contains

  !> convert an array of char * in c to an array
  !  of characters of maximum size 256 in FORTRAN

  subroutine c_f_char_array(c_char_array, &
                            c_size_array, &
                            n_chars, &
                            f_char_array, &
                            free_all)

    use, intrinsic :: iso_c_binding
    implicit none

    type(c_ptr) :: c_char_array, c_size_array
    type(c_ptr), pointer :: fptr2(:) => null()
    character(c_char), pointer :: fptr(:) => null()
    integer(c_int), pointer     :: f_size_array(:)
    character(256), allocatable :: f_char_array(:)
    logical, optional           :: free_all
    integer(c_int) :: i, n_chars, strlen

    call c_f_pointer(c_char_array, fptr2, [n_chars])
    call c_f_pointer(c_size_array, f_size_array, [n_chars])

    allocate(f_char_array(n_chars))
    do i = 1, n_chars
      strlen      = f_size_array(i)
      call c_f_pointer(fptr2(i), fptr, [strlen])
      f_char_array(i) = transfer(fptr(1:strlen), f_char_array(i))
      f_char_array(i) = f_char_array(i)(1:strlen)
      if (present(free_all)) then
        if (free_all) then
          call pdm_fortran_free_c(fptr2(i))
        endif
      endif

    end do

    call pdm_fortran_free_c(c_char_array)
    call pdm_fortran_free_c(c_size_array)

  end subroutine c_f_char_array


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
  !! \param [in]  is_active_rank Current rank is available for CWIPI
  !! \param [out] intra_comms    MPI intra communicators of each code (size = \p n_code)
  !!
  !!

  subroutine CWP_Init_(fcomm,          &
                       n_code,         &
                       code_names,     &
                       is_active_rank, &
                       intra_comms)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int) :: fcomm
    integer(c_int), intent(in) :: n_code
    character(kind = c_char, len = *), dimension(n_code), target :: code_names
    integer(c_int) :: is_active_rank
    integer(c_int), dimension(:), pointer :: intra_comms
    integer, dimension(n_code), target :: l_code_names
    integer :: i

    do i = 1, n_code
      l_code_names(i) = len(code_names(i))
    end do

    call CWP_Init_cf(fcomm, n_code,       &
                     c_loc(code_names),   &
                     c_loc(l_code_names), &
                     is_active_rank,      &
                     c_loc(intra_comms))

  end subroutine CWP_Init_

  !>
  !! \brief Create a Fortran string from a C string
  !!
  !! This function creates a Fortran string from a C string. There is a string copy
  !!
  !! \param [in]  c_str       C string
  !!
  !!  \return      Fortran string
  !!

  function CWP_C_to_f_string_(c_str)  result(f_str)
    use iso_c_binding

    implicit none

    character(kind=c_char,len=1), intent(in) :: c_str(*)
    character(len=:), pointer :: f_str
    integer i, nchars


    i = 1
    do
       if (c_str(i) == c_null_char) exit
       i = i + 1
    end do
    nchars = i - 1  ! Exclude null character from Fortran string
    allocate(character(len=nchars) :: f_str)

    f_str = transfer(c_str(1:nchars), f_str)
  end function

  !>
  !! \brief Update code state.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] state            State
  !!
  !!

  subroutine CWP_State_update_(local_code_name, &
                              state)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    integer(kind = c_int), intent(in) :: state
    integer(c_int)                    :: l_local_code_name

    l_local_code_name = len(local_code_name)

    call CWP_State_update_cf(local_code_name,   &
                             l_local_code_name, &
                             state)

  end subroutine CWP_State_update_

  !>
  !! \brief Begin code time step.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] current_time     Current time
  !!
  !!

  subroutine CWP_Time_step_beg_(local_code_name, &
                             current_time)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    double precision, intent(in)      :: current_time
    integer(c_int)                    :: l_local_code_name

    l_local_code_name = len(local_code_name)

    call CWP_Time_step_beg_cf(local_code_name,   &
                            l_local_code_name, &
                            current_time)

  end subroutine CWP_Time_step_beg_

  !>
  !! \brief End code time step.
  !!
  !! \param [in] local_code_name  Local code name
  !!
  !!

  subroutine CWP_Time_step_end_(local_code_name)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    integer(c_int)                    :: l_local_code_name

    l_local_code_name = len(local_code_name)

    call CWP_Time_step_end_cf(local_code_name,   &
                              l_local_code_name)

  end subroutine CWP_Time_step_end_

  subroutine CWP_User_structure_set_(local_code_name, &
                                    user_structure)
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    type(c_ptr), value                :: user_structure
    integer(c_int)                    :: l_local_code_name

    l_local_code_name = len(local_code_name)

    call CWP_User_structure_set_cf(local_code_name,   &
                                 l_local_code_name, &
                                 user_structure)

  end subroutine CWP_User_structure_set_


  !>
  !! \brief Return the user structure associated
  !!
  !! This structure can be called into a callback
  !!
  !! \param [in] local_code_name  Local code name
  !!
  !! \return  User structure
  !!
  !!

  function CWP_User_structure_get_(local_code_name) &
    result (user_structure)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    type(c_ptr)                       :: user_structure
    integer(c_int)                    :: l_local_code_name

    l_local_code_name = len(local_code_name)

    user_structure = CWP_User_structure_get_cf(local_code_name,   &
                                               l_local_code_name)

  end function CWP_User_structure_get_

  !>
  !! \brief Writing output to Fortran file (shared by fortran and C code).
  !!
  !! This function set the file Fortran logical unit for writing output.
  !!
  !!  \param [in]  iunit        File Fortan logical unit
  !!
  !!

  subroutine cwp_output_fortran_unit_set (outputUnit)

    use, intrinsic :: iso_c_binding
    use cwp_printfort

    implicit none

    integer :: outputUnit

    ifile = outputUnit

    call cwp_set_output_listing_f(outputUnit)

  end subroutine cwp_output_fortran_unit_set

  !>
  !! \brief Define output file (in which only C code writes).
  !!
  !! \param [in] output_file_name    Output file name
  !!
  !!

  subroutine CWP_Output_file_set_ (f_output_file_name)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: f_output_file_name
    integer(c_int)                    :: l_output_file_name

    l_output_file_name = len(f_output_file_name)

    call CWP_Output_file_set_cf (f_output_file_name, &
                                 l_output_file_name)
  end subroutine CWP_Output_file_set_

  !>
  !! \brief Return code state.
  !!
  !! \param [in]  code_name    Code name
  !!
  !! \return      Code state
  !!

  function CWP_State_get_(local_code_name) &
    result (state)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    integer(c_int)                    :: state
    integer(c_int)                    :: l_local_code_name

    l_local_code_name = len(local_code_name)

    state = CWP_State_get_cf(local_code_name,   &
                             l_local_code_name)

  end function CWP_State_get_

    !>
    !! \brief Return list of codes known by CWIPI.
    !!
    !! \return List of code names
    !!

    function CWP_Codes_list_get_() &
      result (fstrings)

      use, intrinsic :: iso_c_binding
      implicit none

      type(c_ptr) :: code_list, code_list_s
      integer(c_int) :: n_codes
      character(256), allocatable :: fstrings(:)

      call CWP_Codes_list_get_cf(code_list, code_list_s, n_codes)

      call c_f_char_array(code_list, code_list_s, n_codes, fstrings)

    end function CWP_Codes_list_get_

    !>
    !! \brief Return list of local codes known by CWIPI.
    !!
    !! \return List of local code names
    !!

    function CWP_Loc_codes_list_get_() &
      result (fstrings)

      use, intrinsic :: iso_c_binding
      implicit none

      type(c_ptr) :: loc_code_list, loc_code_list_s
      integer(c_int) :: n_loc_codes
      character(256), allocatable :: fstrings(:)

      call CWP_Loc_codes_list_get_cf(loc_code_list, loc_code_list_s, n_loc_codes)

      call c_f_char_array(loc_code_list, loc_code_list_s, n_loc_codes, fstrings)

    end function CWP_Loc_codes_list_get_

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

  subroutine CWP_Cpl_create_ (local_code_name,   &
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
  end subroutine CWP_Cpl_Create_


  !>
  !! \brief MPI Barrier on the coupling communicator.
  !!
  !! \param [in]  local_code_name     Local code name
  !! \param [in]  cpl_id              Coupling identifier
  !!
  !!

  subroutine CWP_Cpl_barrier_(local_code_name, &
                              cpl_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int)                    :: l_local_code_name, l_cpl_id

    l_local_code_name   = len(local_code_name)
    l_cpl_id            = len(cpl_id)

    call CWP_Cpl_barrier_cf(local_code_name,   &
                            l_local_code_name, &
                            cpl_id,            &
                            l_cpl_id)
  end subroutine CWP_Cpl_barrier_


  !>
  !!
  !! \brief Delete a coupling object.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !!
  !!

  subroutine CWP_Cpl_Del_ (local_code_name, &
                          cpl_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)

    call CWP_Cpl_del_cf (local_code_name, l_local_code_name, cpl_id, l_cpl_id)
  end subroutine CWP_Cpl_Del_

  !>
  !!
  !! \brief Enable broadcast of the computed targets ids (in \ref CWP_COMM_PAR_WITHOUT_PART mode).
  !!
  !! This function must be called in order for the computed targets to be accessible
  !! on non-root ranks
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] i_part           Current partition
  !!

  subroutine CWP_Computed_tgts_bcast_enable_ (local_code_name, &
                                              cpl_id,          &
                                              field_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    call CWP_Computed_tgts_bcast_enable_cf (local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id)
  end subroutine CWP_Computed_tgts_bcast_enable_

  !>
  !!
  !! \brief Enable broadcast of the involved sources ids (in \ref CWP_COMM_PAR_WITHOUT_PART mode).
  !!
  !! This function must be called in order for the involved sources to be accessible
  !! on non-root ranks
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] i_part           Current partition
  !!

  subroutine CWP_Involved_srcs_bcast_enable_ (local_code_name, &
                                              cpl_id,          &
                                              field_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    call CWP_Involved_srcs_bcast_enable_cf (local_code_name, l_local_code_name, cpl_id, l_cpl_id, field_id, l_field_id)
  end subroutine CWP_Involved_srcs_bcast_enable_

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

  function CWP_N_uncomputed_tgts_get_ (local_code_name, &
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
  end function CWP_N_uncomputed_tgts_get_


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

  function CWP_Uncomputed_tgts_get_ (local_code_name, &
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
  end function CWP_Uncomputed_tgts_get_


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

  function CWP_N_computed_tgts_get_(local_code_name, &
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
  end function CWP_N_computed_tgts_get_


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

  function CWP_Computed_tgts_get_ (local_code_name, &
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
  end function CWP_Computed_tgts_get_


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

  function CWP_N_involved_srcs_get_ (local_code_name, &
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
  end function CWP_N_involved_srcs_get_


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

  function CWP_Involved_srcs_get_ (local_code_name, &
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
  end function CWP_Involved_srcs_get_


  !>
  !! \brief Return distance from each target to the source interface. <b>(Not implemented yet)</b>
  !!
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !!
  !! \return               Distance
  !!
  !!

!   function CWP_Computed_tgts_dist_to_spatial_interp_get (local_code_name, &
!                                                          cpl_id)          &
!                                                          result (dists)

!     use, intrinsic :: iso_c_binding
!     implicit none

!     character(kind = c_char, len = *) :: local_code_name, cpl_id
!     integer(c_int) :: l_local_code_name, l_cpl_id
!     double precision, dimension(:), pointer :: dists
! !!    type(c_ptr) :: cptr_dists

!     l_local_code_name = len(local_code_name)
!     l_cpl_id = len(cpl_id)

!     !cptr_dists = CWP_Computed_tgts_dist_to_spatial_interp_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id)

!     ! TODO The types of return variables may probably not be the right ones
!     print *, "CWP_Computed_tgts_dist_to_spatial_interp_get not implemented"
!     allocate(dists(1))
!     dists = (/0./)
!   end function CWP_Computed_tgts_dist_to_spatial_interp_get


  !>
  !! \brief Compute spatial interpolation weights.
  !!
  !! \param [in]  local_code_name     Local code name
  !! \param [in]  cpl_id              Coupling identifier
  !!

  subroutine CWP_Spatial_interp_weights_compute_ (local_code_name, &
                                                 cpl_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Spatial_interp_weights_compute_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id)
  end subroutine CWP_Spatial_interp_weights_compute_


  !>
  !! \brief Set a property of the spatial interpolation algorithm.
  !!
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !! \param [in]  property_name    Name of the property
  !! \param [in]  property_type    Type of the property
  !! \param [in]  property_value   Value of the property
  !!

  subroutine CWP_Spatial_interp_property_set_(local_code_name, &
                                              cpl_id,          &
                                              property_name,   &
                                              property_type,   &
                                              property_value)
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: cpl_id
    character(kind = c_char, len = *) :: property_name
    integer  (kind = c_int)           :: property_type
    character(kind = c_char, len = *) :: property_value

    call CWP_Spatial_interp_property_set_cf(local_code_name,      &
                                            len(local_code_name), &
                                            cpl_id,               &
                                            len(cpl_id),          &
                                            property_name,        &
                                            len(property_name),   &
                                            property_type,        &
                                            property_value,       &
                                            len(property_value))

  end subroutine CWP_Spatial_interp_property_set_


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

  subroutine CWP_Visu_set_ (local_code_name, &
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

  subroutine CWP_User_tgt_pts_set_(local_code_name, &
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
    double precision, dimension(:,:), pointer :: coord
    integer(kind = c_long), dimension(:), pointer :: global_num
    type(c_ptr) :: c_global_num

    if (associated(global_num)) then
      c_global_num = c_loc(global_num)
    else
      c_global_num = c_null_ptr
    endif

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_User_tgt_pts_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, i_part, n_pts, &
            c_loc(coord), c_global_num)
  end subroutine CWP_User_tgt_pts_set_


  !>
  !! \brief Finalize interface mesh.
  !!
  !! This function computes the global numbers of mesh entities if they are
  !! not provided.
  !!
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !!

  subroutine CWP_Mesh_interf_finalize_ (local_code_name, &
                                       cpl_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_finalize_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id)
  end subroutine CWP_Mesh_interf_finalize_


  !>
  !! \brief Set vertices.
  !!
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !! \param [in]  i_part           Current partition
  !! \param [in]  n_pts            Number of points
  !! \param [in]  coord            Coordinates (size = 3 * \p n_pts)
  !! \param [in]  global_num       Pointer to parent element number (or NULL)
  !!

  subroutine CWP_Mesh_interf_vtx_set_ (local_code_name, &
                                      cpl_id,          &
                                      i_part,          &
                                      n_pts,           &
                                      coord,           &
                                      global_num)

    use :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: i_part, n_pts, array_size
    double precision, dimension(:,:), pointer :: coord
    integer(c_long), dimension(:), pointer :: global_num
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    type(c_ptr) :: c_global_num

    if (associated(global_num)) then
      c_global_num = c_loc(global_num)
    else
      c_global_num = c_null_ptr
    endif

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    array_size = size(coord)
    if (modulo(array_size, 3) /= 0) then
        print *, "Error : Length of connectivity array is not a multiple of the 3"
        stop 'error'
    endif


    call CWP_Mesh_interf_vtx_set_cf (local_code_name,   &
                                     l_local_code_name, &
                                     cpl_id,            &
                                     l_cpl_id,          &
                                     i_part,            &
                                     n_pts,             &
                                     c_loc(coord),      &
                                     c_global_num)

  end subroutine CWP_Mesh_interf_vtx_set_


  !>
  !! \brief Add a connectivity block to the interface mesh.
  !!
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !! \param [in]  block_type       Block type
  !!
  !! \return block identifier
  !!

  function CWP_Mesh_interf_block_add_ (local_code_name, &
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
  end function CWP_Mesh_interf_block_add_


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

  subroutine CWP_Mesh_interf_block_std_set_ (local_code_name, &
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
    ! integer(c_int) :: array_size
    integer(c_int), dimension(:), pointer :: connec
    integer(c_long), dimension(:), pointer :: global_num
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    type(c_ptr) :: c_global_num

    if (associated(global_num)) then
      c_global_num = c_loc(global_num)
    else
      c_global_num = c_null_ptr
    endif

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    ! TODO get stride
    ! array_size = size(connec)
    ! if (modulo(array_size, 3) /= 0) then
    !     print *, "Error : Length of connectivity array is not a multiple of the element size"
    !     stop 'error'
    ! endif

    call CWP_Mesh_interf_block_std_set_cf (local_code_name,   &
                                           l_local_code_name, &
                                           cpl_id,            &
                                           l_cpl_id,          &
                                           i_part,            &
                                           block_id,          &
                                           n_elts,            &
                                           c_loc(connec),     &
                                           c_global_num)
  end subroutine CWP_Mesh_interf_block_std_set_


  !>
  !! \brief Get the properties of a standard block of the interface mesh.
  !!
  !! \param [in]   local_code_name  Local code name
  !! \param [in]   cpl_id           Coupling identifier
  !! \param [in]   i_part           Partition identifier
  !! \param [in]   block_id         Block identifier
  !! \param [out]  n_elts           Number of elements
  !! \param [out]  connec           Connectivity (size = n_vertex_elt * n_elts)
  !! \param [out]  global_num       Pointer to global element number (or NULL)
  !!

  subroutine CWP_Mesh_interf_block_std_get_(local_code_name, &
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
    integer(c_int),  dimension(:), pointer :: connec
    integer(c_long), dimension(:), pointer :: global_num
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, s_connec
    type(c_ptr)      :: c_connec, c_global_num

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)

    call CWP_Mesh_interf_block_std_get_cf(local_code_name,   &
                                          l_local_code_name, &
                                          cpl_id,            &
                                          l_cpl_id,          &
                                          i_part,            &
                                          block_id,          &
                                          n_elts,            &
                                          c_connec,          &
                                          c_global_num,      &
                                          s_connec)

    call c_f_pointer(c_global_num, global_num, [n_elts])
    call c_f_pointer(c_connec,     connec,     [s_connec])

  end subroutine CWP_Mesh_interf_block_std_get_

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

  subroutine CWP_Mesh_interf_f_poly_block_set_(local_code_name, &
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


    type(c_ptr) :: c_global_num

    if (associated(global_num)) then
      c_global_num = c_loc(global_num)
    else
      c_global_num = c_null_ptr
    endif

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
                                              c_global_num)
  end subroutine CWP_Mesh_interf_f_poly_block_set_


  !>
  !! \brief Get the properties of a polygon block of the interface mesh partition.
  !!
  !! \param [in]   local_code_name  Local code name
  !! \param [in]   cpl_id           Coupling identifier
  !! \param [in]   i_part           Current partition
  !! \param [in]   block_id         Block identifier
  !! \param [out]  n_elts           Number of elements
  !! \param [out]  connec_idx       Connectivity index (\p connec_id[0] = 0 and
  !!                               size = \p n_elts + 1)
  !! \param [out]  connec           Connectivity (size = \p connec_idx[\p n_elts])
  !! \param [out]  global_num       Pointer to global element number (or NULL)
  !!
  !!

  subroutine CWP_Mesh_interf_f_poly_block_get_(local_code_name, &
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
    type(c_ptr) :: c_connec_idx, c_connec, c_global_num

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_f_poly_block_get_cf(local_code_name,   &
                                             l_local_code_name, &
                                             cpl_id,            &
                                             l_cpl_id,          &
                                             i_part,            &
                                             block_id,          &
                                             n_elts,            &
                                             c_connec_idx,      &
                                             c_connec,          &
                                             c_global_num)

    call c_f_pointer(c_connec_idx, connec_idx, [n_elts+1])
    call c_f_pointer(c_connec,     connec,     [connec_idx(n_elts+1)])
    call c_f_pointer(c_global_num, global_num, [n_elts])

  end subroutine CWP_Mesh_interf_f_poly_block_get_


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

  subroutine CWP_Mesh_interf_c_poly_block_set_ (local_code_name, &
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

    type(c_ptr) :: c_global_num

    if (associated(global_num)) then
      c_global_num = c_loc(global_num)
    else
      c_global_num = c_null_ptr
    endif

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
                                              c_global_num)
  end subroutine CWP_Mesh_interf_c_poly_block_set_


  !>
  !! \brief Get the properties of a polyhedron block of the interface mesh partition..
  !!
  !! \param [in]  local_code_name   Local code name
  !! \param [in]  cpl_id            Coupling identifier
  !! \param [in]  i_part            Current partition
  !! \param [in]  block_id          Block identifier
  !! \param [out]  n_elts            Number of elements
  !! \param [out]  connec_cells_idx  Polyhedron to face index
  !!                                (\p src_poly_cell_face_idx[0] = 0 and
  !!                                 size = \p n_elts + 1)
  !! \param [out]  connec_cells      Polyhedron to face connectivity
  !!                                (size = \p cell_face_idx[\p n_elts])
  !! \param [out]  n_faces           Number of faces
  !! \param [out]  connec_faces_idx  Polyhedron face to vertex index
  !!                                (\p face_vertex_idx[0] = 0 and
  !!                                 size = max(\p cell_face_connec) + 1)
  !! \param [out]  connec_faces      Polyhedron face to vertex connectivity
  !!                                (size = \p face_vertex_idx[\p n_elts])
  !! \param [out]  global_num        Pointer to global element number (or NULL)
  !!
  !!

  subroutine CWP_Mesh_interf_c_poly_block_get_(local_code_name,  &
                                              cpl_id,           &
                                              i_part,           &
                                              block_id,         &
                                              n_elts,           &
                                              n_faces,          &
                                              connec_faces_idx, &
                                              connec_faces,     &
                                              connec_cells_idx, &
                                              connec_cells,     &
                                              global_num)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: i_part, block_id, n_elts, n_faces
    integer(c_int),  dimension(:), pointer :: connec_faces_idx, connec_faces, connec_cells_idx, connec_cells
    integer(c_long), dimension(:), pointer :: global_num
    integer(kind = c_int) :: l_local_code_name, l_cpl_id
    type(c_ptr) :: c_connec_faces_idx, c_connec_faces, c_connec_cells_idx, c_connec_cells, c_global_num

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)

    call CWP_Mesh_interf_c_poly_block_get_cf(local_code_name,    &
                                             l_local_code_name,  &
                                             cpl_id,             &
                                             l_cpl_id,           &
                                             i_part,             &
                                             block_id,           &
                                             n_elts,             &
                                             n_faces,            &
                                             c_connec_faces_idx, &
                                             c_connec_faces,     &
                                             c_connec_cells_idx, &
                                             c_connec_cells,     &
                                             c_global_num)

    call c_f_pointer(c_connec_faces_idx, connec_faces_idx, [n_faces+1])
    call c_f_pointer(c_connec_faces,     connec_faces,     [connec_faces_idx(n_faces+1)])
    call c_f_pointer(c_connec_cells_idx, connec_cells_idx, [n_elts+1])
    call c_f_pointer(c_connec_cells,     connec_cells,     [connec_cells_idx(n_elts+1)])
    call c_f_pointer(c_global_num, global_num, [n_elts])

  end subroutine CWP_Mesh_interf_c_poly_block_get_


  !>
  !! \brief Delete interface mesh.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !!
  !!

  subroutine CWP_Mesh_interf_del_ (local_code_name, &
                                  cpl_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: l_local_code_name, l_cpl_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)

    call CWP_Mesh_interf_del_cf (local_code_name, l_local_code_name, cpl_id, l_cpl_id)
  end subroutine CWP_Mesh_interf_del_


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
  !! \param [in]  global_num        Pointer to parent element number (or NULL)
  !!

  subroutine CWP_Mesh_interf_from_cellface_set_ (local_code_name, &
                                                cpl_id,          &
                                                i_part,          &
                                                n_cells,         &
                                                cell_face_idx,   &
                                                cell_face,       &
                                                n_faces,         &
                                                face_vtx_idx,    &
                                                face_vtx,        &
                                                global_num)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: i_part, n_cells, n_faces
    integer(c_int), dimension(:), pointer :: cell_face_idx, cell_face, face_vtx_idx, face_vtx
    integer(c_long), dimension(:), pointer :: global_num

    integer(kind = c_int) :: l_local_code_name, l_cpl_id
    type(c_ptr) :: c_global_num

    if (associated(global_num)) then
      c_global_num = c_loc(global_num)
    else
      c_global_num = c_null_ptr
    endif

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
                                               c_global_num)
  end subroutine CWP_Mesh_interf_from_cellface_set_


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
  !! \param [in]  edge_vtx          Edge to vertex connectivity
  !!                                (size = 2 * \p n_edges)
  !! \param [in]  global_num        Pointer to parent element number (or NULL)
  !!

  subroutine CWP_Mesh_interf_from_faceedge_set_ (local_code_name, &
                                                cpl_id,          &
                                                i_part,          &
                                                n_faces,         &
                                                face_edge_idx,   &
                                                face_edge,       &
                                                n_edges,         &
                                                edge_vtx,        &
                                                global_num)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: i_part, n_faces, n_edges
    integer(c_int), dimension(:), pointer :: face_edge_idx, face_edge, edge_vtx
    integer(c_long), dimension(:), pointer :: global_num
    integer(kind = c_int) :: l_local_code_name, l_cpl_id

    type(c_ptr) :: c_global_num

    if (associated(global_num)) then
      c_global_num = c_loc(global_num)
    else
      c_global_num = c_null_ptr
    endif

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
                                               c_loc(edge_vtx),     &
                                               c_global_num)

  end subroutine CWP_Mesh_interf_from_faceedge_set_


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

  subroutine CWP_Field_create_(local_code_name,      &
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
  end subroutine CWP_Field_create_


  !>
  !!
  !! \brief Set field data.
  !!
  !! \param [in] local_code_name   Local code name
  !! \param [in] cpl_id            Coupling identifier
  !! \param [in] field_id          Field identifier
  !! \param [in] i_part            Current partition
  !! \param [in] data_type         Choice if data is set for the source or the target
  !! \param [in] data              Storage array (Mapping)
  !!

  subroutine CWP_Field_data_set_(local_code_name, &
                                 cpl_id,          &
                                 field_id,        &
                                 i_part,          &
                                 map_type,        &
                                 data)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(c_int) :: i_part
    integer(c_int) :: n_dof, n_components, array_size, spatial_interp_algorithm
    integer(c_int) :: map_type
    double precision, dimension(:), pointer :: data
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    ! check data array size
    array_size = size(data)
    n_dof = CWP_Field_n_dof_get(local_code_name,   &
                                cpl_id,            &
                                field_id,          &
                                i_part)
    n_components = CWP_Field_n_components_get(local_code_name,   &
                                              cpl_id,            &
                                              field_id)
    spatial_interp_algorithm = CWP_Cpl_spatial_interp_algo_get(local_code_name,   &
                                                               cpl_id)

    if (spatial_interp_algorithm .eq. CWP_SPATIAL_INTERP_FROM_NEAREST_SOURCES_LEAST_SQUARES .or. &
        spatial_interp_algorithm .eq. CWP_SPATIAL_INTERP_FROM_NEAREST_TARGETS_LEAST_SQUARES .or. &
        spatial_interp_algorithm .eq. CWP_SPATIAL_INTERP_FROM_INTERSECTION .or.                  &
        spatial_interp_algorithm .eq. CWP_SPATIAL_INTERP_FROM_IDENTITY) then

        if (map_type .eq. CWP_FIELD_MAP_TARGET) then

            if (array_size < n_dof * n_components) then
                print *, "Error : Data array size is inconsistent (too small) with the number of components &
                          & and degree of freedom location"
                stop
            endif

        endif

    else if (spatial_interp_algorithm .eq. CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT .or. &
             spatial_interp_algorithm .eq. CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE .or.         &
             spatial_interp_algorithm .eq. CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_BOXTREE) then

        if (map_type .eq. CWP_FIELD_MAP_SOURCE) then

            if (array_size < n_dof * n_components) then
                print *, "Error : Data array size is inconsistent (too small) with the number of components &
                          & and degree of freedom location"
                stop
            endif

        endif

    endif

    call CWP_Field_data_set_cf (local_code_name,   &
                                l_local_code_name, &
                                cpl_id,            &
                                l_cpl_id,          &
                                field_id,          &
                                l_field_id,        &
                                i_part,            &
                                map_type,          &
                                c_loc(data))
  end subroutine CWP_Field_data_set_


  !>
  !!
  !! \brief Get target degrees of freedom location.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] field_id         Field identifier
  !!
  !! \return                      Location of degrees of freedom
  !!
  !!

  function CWP_Field_dof_location_get_(local_code_name, &
                                       cpl_id,          &
                                       field_id)        &
    result (dof_location)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(c_int) :: dof_location
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    dof_location = CWP_Field_dof_location_get_cf(local_code_name,   &
                                                 l_local_code_name, &
                                                 cpl_id,            &
                                                 l_cpl_id,          &
                                                 field_id,          &
                                                 l_field_id)

  end function CWP_Field_dof_location_get_


  !>
  !!
  !! \brief Get field storage type.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] field_id         Field identifier
  !!
  !! \return                      Field storage type
  !!

  function CWP_Field_storage_get_(local_code_name, &
                                  cpl_id,          &
                                  field_id)        &
    result (storage)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(c_int) :: storage
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    storage = CWP_Field_storage_get_cf(local_code_name,   &
                                       l_local_code_name, &
                                       cpl_id,            &
                                       l_cpl_id,          &
                                       field_id,          &
                                       l_field_id)

  end function CWP_Field_storage_get_

  !>
  !!
  !! \brief Get field number of degrees of freedom.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] field_id         Field identifier
  !! \param [in] n_dof            Field number of degrees of freedom
  !!
  !! \return                      Field storage type
  !!

  function CWP_Field_n_dof_get_(local_code_name, &
                                cpl_id,          &
                                field_id,        &
                                i_part)          &
    result (n_dof)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(c_int) :: i_part
    integer(c_int) :: n_dof
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    n_dof = CWP_Field_n_dof_get_cf(local_code_name,   &
                                   l_local_code_name, &
                                   cpl_id,            &
                                   l_cpl_id,          &
                                   field_id,          &
                                   l_field_id,        &
                                   i_part)

  end function CWP_Field_n_dof_get_


  !>
  !! \brief Delete a field.
  !!
  !! \param [in]  local_code_name Local code name
  !! \param [in]  cpl_id          Coupling identifier
  !! \param [in]  field_id        Field identifier
  !!
  !!

  subroutine CWP_Field_del_(local_code_name, &
                            cpl_id,          &
                            field_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    call CWP_Field_del_cf(local_code_name,   &
                          l_local_code_name, &
                          cpl_id,            &
                          l_cpl_id,          &
                          field_id,          &
                          l_field_id)

  end subroutine CWP_Field_del_

  !>
  !! \brief Send a spatially interpolated field to the coupled code with
  !!        non-blocking communications.
  !!
  !! This function is independent of \ref CWP_Time_exch_t mode. The user has to
  !! manually check the consistency of the exchanges.
  !!
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !! \param [in]  field_id         Field identifier
  !!

  subroutine CWP_Field_issend_(local_code_name, &
                               cpl_id,          &
                               field_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_issend_cf (local_code_name,    &
                              l_local_code_name,  &
                              cpl_id,             &
                              l_cpl_id,           &
                              field_id,       &
                              l_field_id)
  end subroutine CWP_Field_issend_


  !>
  !!
  !! \brief Receive a spatially interpolated field from the coupled code
  !!        with non-blocking communications.
  !!
  !! This function is independent of \ref CWP_Time_exch_t mode. The user has to
  !! manually check the consistency of the exchanges.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in]  cpl_id          Coupling identifier
  !! \param [in]  tgt_field_id    Target field id
  !!

  subroutine CWP_Field_irecv_(local_code_name, &
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
  end subroutine CWP_Field_irecv_


  !>
  !!
  !! \brief Wait the end of an exchange related to request from \ref CWP_Field_issend.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] field_id         Field identifier
  !!

  subroutine CWP_Field_wait_issend_(local_code_name, &
                                    cpl_id,          &
                                    field_id)
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_wait_issend_cf (local_code_name,    &
                                   l_local_code_name,  &
                                   cpl_id,             &
                                   l_cpl_id,           &
                                   field_id,       &
                                   l_field_id)
  end subroutine CWP_Field_wait_issend_


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

  subroutine CWP_Field_wait_irecv_ (local_code_name, &
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
  end subroutine CWP_Field_wait_irecv_

  !>
  !!
  !!  \brief Unsetting of a user interpolation.
  !!
  !!  \param [in] local_code_name  Local code name
  !!  \param [in] cpl_id           Coupling identifier
  !!  \param [in] field_id         Field identifier
  !!
  !!

  subroutine CWP_Field_interp_function_unset_ (local_code_name, &
                                         cpl_id,          &
                                         field_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_interp_function_unset_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
            & field_id, l_field_id)

  end subroutine CWP_Field_interp_function_unset_

  !>
  !!
  !! \brief Setting of a user interpolation from location.
  !!
  !! This function takes into account an user interpolation function written with
  !! void (*\ref CWP_Field_interp_function_t) interface.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] field_id         Field identifier
  !! \param [in] fct              Function
  !!

  subroutine CWP_Field_interp_function_set_(local_code_name,        &
                                      cpl_id,                 &
                                      field_id,               &
                                      user_interpolation_fct)

    use, intrinsic :: iso_c_binding
    implicit none

    interface
        subroutine user_interpolation_fct(c_local_code_name, &
                                          c_cpl_id,          &
                                          c_field_id,        &
                                          i_part,            &
                                          c_buffer_in,       &
                                          c_buffer_out)      &
          bind(c)
          use, intrinsic :: iso_c_binding
          implicit none

          character(kind=c_char,len=1) :: c_local_code_name(*)
          character(kind=c_char,len=1) :: c_cpl_id(*)
          character(kind=c_char,len=1) :: c_field_id(*)
          integer(kind=c_int), value   :: i_part
          type(c_ptr), value           :: c_buffer_in
          type(c_ptr), value           :: c_buffer_out
        end subroutine user_interpolation_fct
    end interface

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(kind = c_int)             :: l_local_code_name, l_cpl_id, l_field_id
    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    call CWP_Field_interp_function_set_cf(local_code_name,   &
                                    l_local_code_name, &
                                    cpl_id,            &
                                    l_cpl_id,          &
                                    field_id,      &
                                    l_field_id,    &
                                    c_funloc(user_interpolation_fct))
  end subroutine CWP_Field_interp_function_set_

  !>
  !!
  !! \brief Get spatial interpolation number of algorithms.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] field_id         Field identifier
  !!

  function CWP_Field_n_components_get_ (local_code_name, &
                                        cpl_id,          &
                                        field_id)        &
    result (n_components)
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    integer                :: n_components
    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    n_components = CWP_Field_n_components_get_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
                                                 & field_id, l_field_id)

  end function CWP_Field_n_components_get_

  !>
  !!
  !!  \brief Get spatial interpolation source data.
  !!
  !!  \param [in]  local_code_name  Local code name
  !!  \param [in]  cpl_id           Coupling identifier
  !!  \param [in]  field_id         Field identifier
  !!  \param [out] i_part           Partition identifier
  !!  \param [out] n_elt_src        Number of local source entities in current partition
  !!  \param [out] src_to_tgt_idx   Index for source->target mapping
  !!

  subroutine CWP_Field_src_data_properties_get_(local_code_name, &
                                                cpl_id, &
                                                field_id, &
                                                i_part, &
                                                n_elt_src, &
                                                src_to_tgt_idx)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    integer(c_int) :: i_part, n_elt_src
    integer(c_int),  dimension(:), pointer :: src_to_tgt_idx
    type(c_ptr)      :: c_src_to_tgt_idx

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_src_data_properties_get_cf(local_code_name, &
                                              l_local_code_name, &
                                              cpl_id, &
                                              l_cpl_id, &
                                              field_id, &
                                              l_field_id, &
                                              i_part, &
                                              n_elt_src, &
                                              c_src_to_tgt_idx)

    call c_f_pointer(c_src_to_tgt_idx, src_to_tgt_idx, [n_elt_src+1])

  end subroutine CWP_Field_src_data_properties_get_

  !>
  !!
  !!  \brief Get spatial interpolation target data.
  !!
  !!  \param [in]  local_code_name       Local code name
  !!  \param [in]  cpl_id                Coupling identifier
  !!  \param [in]  field_id              Field identifier
  !!  \param [out] i_part                Partition identifier
  !!  \param [out] n_elt_tgt             Number of local target entities in current partition
  !!  \param [out] n_referenced_tgt      Number of referenced target entities in current partition
  !!  \param [out] referenced_tgt        Ids of referenced target entities in current partition (1-based)
  !!  \param [out] tgt_come_from_src_idx Index for target->source mapping
  !!

  subroutine CWP_Field_tgt_data_properties_get_(local_code_name, &
                                                cpl_id, &
                                                field_id, &
                                                i_part, &
                                                n_elt_tgt, &
                                                n_referenced_tgt, &
                                                referenced_tgt, &
                                                tgt_come_from_src_idx)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    integer(c_int) :: i_part, n_elt_tgt, n_referenced_tgt
    integer(c_int),  dimension(:), pointer :: referenced_tgt, tgt_come_from_src_idx
    type(c_ptr)      :: c_referenced_tgt, c_tgt_come_from_src_idx

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_tgt_data_properties_get_cf(local_code_name, &
                                              l_local_code_name, &
                                              cpl_id, &
                                              l_cpl_id, &
                                              field_id, &
                                              l_field_id, &
                                              i_part, &
                                              n_elt_tgt, &
                                              n_referenced_tgt, &
                                              c_referenced_tgt, &
                                              c_tgt_come_from_src_idx)

    call c_f_pointer(c_referenced_tgt, referenced_tgt, [n_referenced_tgt])
    call c_f_pointer(c_tgt_come_from_src_idx, tgt_come_from_src_idx, [n_referenced_tgt+1])

  end subroutine CWP_Field_tgt_data_properties_get_

  !>
  !!
  !!  \brief Get spatial interpolation weights (location algorithm).
  !!
  !!  \param [in]  local_code_name  Local code name
  !!  \param [in]  cpl_id           Coupling identifier
  !!  \param [in]  field_id         Field identifier
  !!  \param [in]  i_part           Partition identifier
  !!  \param [out] weights          Interpolation weights
  !!

  subroutine CWP_Field_location_weights_get_(local_code_name, &
                                             cpl_id,          &
                                             field_id,        &
                                             i_part,          &
                                             weights)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    integer(c_int) :: i_part
    double precision, dimension(:), pointer :: weights
    integer(kind = c_int) :: s_weights
    type(c_ptr)           :: c_weights

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_location_weights_get_cf(local_code_name, &
                                           l_local_code_name, &
                                           cpl_id, &
                                           l_cpl_id, &
                                           field_id, &
                                           l_field_id, &
                                           i_part, &
                                           c_weights, &
                                           s_weights)

    call c_f_pointer(c_weights, weights, [s_weights])

  end subroutine CWP_Field_location_weights_get_

  !>
  !!
  !!  \brief Get spatial interpolation point data (location algorithm).
  !!
  !!  \param [in]  local_code_name          Local code name
  !!  \param [in]  cpl_id                   Coupling identifier
  !!  \param [in]  field_id         Field identifier
  !!  \param [in]  i_part                   Partition identifier
  !!  \param [out] points_coords            Cartesian coordinates of points inside local elements
  !!  \param [out] points_uvw               Parametric coordinates of points inside local elements
  !!  \param [out] points_dist2             Squared distance from points to elements
  !!  \param [out] points_projected_coords  Cartesian coordinates of projection on points on local elements
  !!

  subroutine CWP_Field_location_point_data_get_(local_code_name, &
                                                 cpl_id, &
                                                 field_id, &
                                                 i_part, &
                                                 points_coords, &
                                                 points_uvw, &
                                                 points_dist2, &
                                                 points_projected_coords)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    integer(c_int) :: i_part
    double precision, pointer :: points_coords(:,:)
    double precision, pointer :: points_uvw(:,:)
    double precision, pointer :: points_dist2(:)
    double precision, pointer :: points_projected_coords(:,:)
    type(c_ptr)      :: c_points_coords, c_points_uvw, c_points_dist2, c_points_projected_coords
    integer(kind = c_int) :: s_size

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_location_point_data_get_cf(local_code_name, &
                                               l_local_code_name, &
                                               cpl_id, &
                                               l_cpl_id, &
                                               field_id, &
                                               l_field_id, &
                                               i_part, &
                                               c_points_coords, &
                                               c_points_uvw, &
                                               c_points_dist2, &
                                               c_points_projected_coords, &
                                               s_size)

    call c_f_pointer(c_points_coords, points_coords, [3, s_size])
    call c_f_pointer(c_points_uvw, points_uvw, [3, s_size])
    call c_f_pointer(c_points_dist2, points_dist2, [s_size])
    call c_f_pointer(c_points_projected_coords, points_projected_coords, [3, s_size])

  end subroutine CWP_Field_location_point_data_get_

  !>
  !!
  !!  \brief Get spatial interpolation internal cell->vertex connectivity (location algorithm).
  !!
  !!  \param [in]  local_code_name  Local code name
  !!  \param [in]  cpl_id           Coupling identifier
  !!  \param [in]  field_id         Field identifier
  !!  \param [in]  i_part           Partition identifier
  !!  \param [out] cell_vtx_idx     Index for local cell->vertex connectivity
  !!  \param [out] cell_vtx         Local cell->vertex connectivity
  !!

  subroutine CWP_Field_location_internal_cell_vtx_get_(local_code_name, &
                                                        cpl_id, &
                                                        field_id, &
                                                        i_part, &
                                                        cell_vtx_idx, &
                                                        cell_vtx)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    integer(c_int) :: i_part
    integer(c_int),  dimension(:), pointer :: cell_vtx_idx, cell_vtx
    type(c_ptr)      :: c_cell_vtx_idx, c_cell_vtx
    integer(kind = c_int) :: n_cell

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_location_internal_cell_vtx_get_cf(local_code_name, &
                                                      l_local_code_name, &
                                                      cpl_id, &
                                                      l_cpl_id, &
                                                      field_id, &
                                                      l_field_id, &
                                                      i_part, &
                                                      c_cell_vtx_idx, &
                                                      n_cell, &
                                                      c_cell_vtx)

    call c_f_pointer(c_cell_vtx_idx, cell_vtx_idx, [n_cell+1])
    call c_f_pointer(c_cell_vtx, cell_vtx, [cell_vtx_idx(n_cell+1)])

  end subroutine CWP_Field_location_internal_cell_vtx_get_

  !>
  !!
  !!  \brief Get spatial interpolation volumes (intersection algorithm).
  !!
  !!  \param [in]  local_code_name  Local code name
  !!  \param [in]  cpl_id           Coupling identifier
  !!  \param [in]  field_id         Field identifier
  !!  \param [in]  i_part           Partition identifier
  !!  \param [out] volumes          Volumes of intersection polyhedra
  !!

  subroutine CWP_Field_intersection_volumes_get_(local_code_name, &
                                                  cpl_id, &
                                                  field_id, &
                                                  i_part, &
                                                  volumes)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    integer(c_int) :: i_part
    double precision,  dimension(:), pointer :: volumes
    integer(kind = c_int) :: s_volumes
    type(c_ptr)           :: c_volumes

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_intersection_volumes_get_cf(local_code_name, &
                                                l_local_code_name, &
                                                cpl_id, &
                                                l_cpl_id, &
                                                field_id, &
                                                l_field_id, &
                                                i_part, &
                                                c_volumes, &
                                                s_volumes)

    call c_f_pointer(c_volumes, volumes, [s_volumes])

  end subroutine CWP_Field_intersection_volumes_get_


  !>
  !!
  !!  \brief Get spatial local target elements volumes (intersection algorithm).
  !!
  !!  \param [in]  local_code_name           Local code name
  !!  \param [in]  cpl_id                    Coupling identifier
  !!  \param [in]  field_id                  Field identifier
  !!  \param [in]  i_part                    Partition identifier
  !!  \param [out] tgt_elt_volumes           Volumes of local target elements
  !!

  subroutine CWP_Field_intersection_tgt_elt_volumes_get_(local_code_name, &
                                                          cpl_id,          &
                                                          field_id,        &
                                                          i_part,          &
                                                          tgt_elt_volumes)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind=c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(c_int)                  :: l_local_code_name, l_cpl_id, l_field_id
    integer(c_int)                  :: i_part
    double precision, pointer       :: tgt_elt_volumes(:)
    type(c_ptr)                     :: c_tgt_elt_volumes = C_NULL_PTR
    integer(c_int)                  :: n_elt

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_intersection_tgt_elt_volumes_get_cf(local_code_name,   &
                                                        l_local_code_name, &
                                                        cpl_id,            &
                                                        l_cpl_id,          &
                                                        field_id,          &
                                                        l_field_id,        &
                                                        i_part,            &
                                                        c_tgt_elt_volumes, &
                                                        n_elt)

    call c_f_pointer(c_tgt_elt_volumes, tgt_elt_volumes, [n_elt])

  end subroutine CWP_Field_intersection_tgt_elt_volumes_get_

  !>
  !!
  !!  \brief Get spatial interpolation distances (nearest neighbors algorithm).
  !!
  !!  \param [in]  local_code_name  Local code name
  !!  \param [in]  cpl_id           Coupling identifier
  !!  \param [in]  field_id         Field identifier
  !!  \param [in]  i_part           Partition identifier
  !!  \param [out] distances2       Squared distances from nearest source points
  !!

  subroutine CWP_Field_nearest_neighbors_distances_get_(local_code_name, &
                                                         cpl_id,          &
                                                         field_id,        &
                                                         i_part,          &
                                                         distances2)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id
    integer(c_int) :: i_part
    double precision,  dimension(:), pointer :: distances2
    integer(kind = c_int) :: s_distances2
    type(c_ptr)           :: c_distances2

    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_field_id = len(field_id)

    call CWP_Field_nearest_neighbors_distances_get_cf(local_code_name, &
                                                       l_local_code_name, &
                                                       cpl_id, &
                                                       l_cpl_id, &
                                                       field_id, &
                                                       l_field_id, &
                                                       i_part, &
                                                       c_distances2, &
                                                       s_distances2)

    call c_f_pointer(c_distances2, distances2, [s_distances2])

  end subroutine CWP_Field_nearest_neighbors_distances_get_


  !>
  !!
  !!  \brief Get coordinates of nearest source points (nearest neighbors algorithm).
  !!
  !!  \param [in]  local_code_name    Local code name
  !!  \param [in]  cpl_id             Coupling identifier
  !!  \param [in]  field_id           Field identifier
  !!  \param [in]  i_part             Partition identifier
  !!  \param [out] nearest_src_coord  Coordinates of nearest source points
  !!

  subroutine CWP_Field_nearest_neighbors_coord_get_(local_code_name,   &
                                                     cpl_id,            &
                                                     field_id,          &
                                                     i_part,            &
                                                     nearest_src_coord)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind=c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(c_int)                  :: l_local_code_name, l_cpl_id, l_field_id
    integer(c_int), intent(in)      :: i_part
    double precision, pointer       :: nearest_src_coord(:,:)

    type(c_ptr)                     :: c_nearest_src_coord = C_NULL_PTR
    integer(c_int)                  :: n_nearest_src_pts

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    call CWP_Field_nearest_neighbors_coord_get_cf(local_code_name,     &
                                                  l_local_code_name,   &
                                                  cpl_id,              &
                                                  l_cpl_id,            &
                                                  field_id,            &
                                                  l_field_id,          &
                                                  i_part,              &
                                                  c_nearest_src_coord, &
                                                  n_nearest_src_pts)

    call c_f_pointer(c_nearest_src_coord, nearest_src_coord, [3, n_nearest_src_pts])

  end subroutine CWP_Field_nearest_neighbors_coord_get_

! /*----------------------------------------------------------------------------*
!  * Functions about all code parameters                                        *
!  *----------------------------------------------------------------------------*/


  !>
  !!
  !! \brief Add a new parameter and initialize it.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] param_name       Parameter name
  !! \param [in] data_type        Parameter type
  !! \param [in] initial_value    Initial value
  !!
  !!

  subroutine CWP_Param_add_int_(local_code_name, &
                                param_name,      &
                                initial_value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    integer(kind = c_int), intent(in) :: initial_value
    integer(kind = c_int)             :: l_local_code_name, l_param_name

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)

    call CWP_Param_add_int_cf(local_code_name,   &
                              l_local_code_name, &
                              param_name,        &
                              l_param_name,      &
                              initial_value)

  end subroutine CWP_Param_add_int_


  subroutine CWP_Param_add_double_(local_code_name, &
                                   param_name,      &
                                   initial_value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    real(kind = c_double), intent(in) :: initial_value
    integer(kind = c_int)             :: l_local_code_name, l_param_name

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)

    call CWP_Param_add_double_cf(local_code_name,   &
                                 l_local_code_name, &
                                 param_name,        &
                                 l_param_name,      &
                                 initial_value)

  end subroutine CWP_Param_add_double_

  subroutine CWP_Param_add_char_(local_code_name, &
                                 param_name,      &
                                 initial_value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    character(kind = c_char, len = *) :: initial_value
    integer(kind = c_int)             :: l_local_code_name, l_param_name, l_initial_value

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)
    l_initial_value   = len(initial_value)

    call CWP_Param_add_char_cf(local_code_name,   &
                               l_local_code_name, &
                               param_name,        &
                               l_param_name,      &
                               initial_value,     &
                               l_initial_value)

  end subroutine CWP_Param_add_char_

  !>
  !!
  !! \brief Set a parameter.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] param_name       Parameter name
  !! \param [in] data_type        Parameter type
  !! \param [in] value            Value
  !!
  !!

  subroutine CWP_Param_set_int_(local_code_name, &
                                param_name,      &
                                value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    integer(kind = c_int), intent(in) :: value
    integer(kind = c_int)             :: l_local_code_name, l_param_name

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)

    call CWP_Param_set_int_cf(local_code_name,   &
                              l_local_code_name, &
                              param_name,        &
                              l_param_name,      &
                              value)

  end subroutine CWP_Param_set_int_


  subroutine CWP_Param_set_double_(local_code_name, &
                                   param_name,      &
                                   value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    real(kind = c_double), intent(in) :: value
    integer(kind = c_int)             :: l_local_code_name, l_param_name

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)

    call CWP_Param_set_double_cf(local_code_name,   &
                                 l_local_code_name, &
                                 param_name,        &
                                 l_param_name,      &
                                 value)

  end subroutine CWP_Param_set_double_

  subroutine CWP_Param_set_char_(local_code_name, &
                                 param_name,      &
                                 value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    character(kind = c_char, len = *) :: value
    integer(kind = c_int)             :: l_local_code_name, l_param_name, l_value

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)
    l_value           = len(value)

    call CWP_Param_set_char_cf(local_code_name,   &
                               l_local_code_name, &
                               param_name,        &
                               l_param_name,      &
                               value,             &
                               l_value)

  end subroutine CWP_Param_set_char_


  !>
  !!
  !! \brief Delete a parameter.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] param_name       Parameter name
  !! \param [in] data_type        Parameter type,
  !!
  !!

  subroutine CWP_Param_del_(local_code_name, &
                            param_name,      &
                            data_type)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    integer, intent(in)               :: data_type
    integer(kind = c_int)             :: l_local_code_name, l_param_name

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)

    call CWP_Param_del_cf(local_code_name,   &
                          l_local_code_name, &
                          param_name,        &
                          l_param_name,      &
                          data_type)

  end subroutine CWP_Param_del_


! /*----------------------------------------------------------------------------*
!  * Functions about all code parameters                                        *
!  *----------------------------------------------------------------------------*/

  !>
  !!
  !! \brief Return the number of parameters for the code \p code_name.
  !!
  !! \param [in] code_name       Local or distant code name
  !! \param [in] data_type       Parameter type,
  !!
  !! return  Number of parameters
  !!
  !!

  function CWP_Param_n_get_(code_name, &
                            data_type) &
    result (n_param)
    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name
    integer, intent(in)               :: data_type
    integer(kind = c_int)             :: l_code_name
    integer                           :: n_param

    l_code_name = len(code_name)

    n_param = CWP_Param_n_get_cf(code_name,   &
                                 l_code_name, &
                                 data_type)

  end function CWP_Param_n_get_

  !>
  !!
  !! \brief Return the list of parameters for the code \p code_name.
  !!
  !! \param [in] code_name      Local or distant code name
  !! \param [in] data_type      Parameter type
  !! \param [in] n_param        Number of parameters
  !! \param [in] param_names    Parameter names
  !!

  subroutine CWP_Param_list_get_(code_name, &
                                 data_type, &
                                 n_param,   &
                                 f_param_names)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name
    integer(kind = c_int)             :: l_code_name
    integer                           :: data_type
    type(c_ptr)                       :: c_param_names, c_param_sizes
    integer(c_int)                    :: n_param
    character(256), allocatable       :: f_param_names(:)

    l_code_name  = len(code_name)

    call CWP_Param_list_get_cf(code_name, l_code_name, data_type, n_param, c_param_names, c_param_sizes)

    call c_f_char_array(c_param_names, c_param_sizes, n_param, f_param_names, .true.)

  end subroutine CWP_Param_list_get_

  !>
  !!
  !! \brief Is this \p code_name a parameter ?
  !!
  !! \param [in] code_name      Local or distant code name
  !! \param [in] param_name     Parameter name
  !! \param [in] data_type      Parameter type,
  !!
  !! return  1 : true / 0 : false
  !!
  !!

  function CWP_Param_is_(code_name,  &
                        param_name, &
                        data_type)  &
    result (is_param)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name
    character(kind = c_char, len = *) :: param_name
    integer, intent(in)               :: data_type
    integer(kind = c_int)             :: l_code_name
    integer(kind = c_int)             :: l_param_name
    integer                           :: is_param

    l_code_name  = len(code_name)
    l_param_name = len(param_name)

    is_param = CWP_Param_is_cf(code_name,    &
                               l_code_name,  &
                               param_name,   &
                               l_param_name, &
                               data_type)

  end function CWP_Param_is_


  !>
  !!
  !! \brief Return the parameter value of \p param_name on \p code_name.
  !!
  !! \param [in]  code_name  Local or distant code name
  !! \param [in]  param_name Parameter name
  !! \param [in]  data_type  Parameter type
  !! \param [out] value      Parameter value
  !!
  !!

  subroutine CWP_Param_get_int(code_name,  &
                               param_name, &
                               value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name
    character(kind = c_char, len = *) :: param_name
    integer(c_int), intent(out)       :: value
    integer(kind = c_int)             :: l_code_name
    integer(kind = c_int)             :: l_param_name

    l_code_name  = len(code_name)
    l_param_name = len(param_name)

    call CWP_Param_get_int_cf(code_name,    &
                              l_code_name,  &
                              param_name,   &
                              l_param_name, &
                              value)

  end subroutine CWP_Param_get_int


  subroutine CWP_Param_get_double(code_name,  &
                                  param_name, &
                                  value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name
    character(kind = c_char, len = *) :: param_name
    real(c_double), intent(out)       :: value
    integer(kind = c_int)             :: l_code_name
    integer(kind = c_int)             :: l_param_name

    l_code_name  = len(code_name)
    l_param_name = len(param_name)

    call CWP_Param_get_double_cf(code_name,    &
                                 l_code_name,  &
                                 param_name,   &
                                 l_param_name, &
                                 value)

  end subroutine CWP_Param_get_double

  subroutine CWP_Param_get_char(code_name,  &
                                param_name, &
                                val)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name
    character(kind = c_char, len = *) :: param_name
    integer(kind = c_int)             :: l_code_name
    integer(kind = c_int)             :: l_param_name
    character(c_char), pointer        :: val(:)
    type(c_ptr)                       :: cptr = C_NULL_PTR
    integer(c_int)                    :: l_value

    l_code_name  = len(code_name)
    l_param_name = len(param_name)

    call CWP_Param_get_char_cf(code_name,    &
                               l_code_name,  &
                               param_name,   &
                               l_param_name, &
                               cptr,         &
                               l_value)

    call c_f_pointer(cptr, val, [l_value])

  end subroutine CWP_Param_get_char


  !>
  !!
  !! \brief Return the result of a reduce operation about a parameter
  !!
  !! \param [in]  op           Operation
  !! \param [in]  param_name   Parameter name
  !! \param [in]  data_type    Parameter type,
  !! \param [out] res          Result
  !! \param [in]  n_codes      Number of codes
  !! \param [in]  code_names   Codes name
  !!

  subroutine CWP_Param_reduce_int(op,         &
                                  param_name, &
                                  res,        &
                                  n_codes,    &
                                  code_names)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                           :: op
    integer(c_int), intent(out)                                   :: res
    integer(c_int)                                                :: n_codes, i
    character(kind = c_char, len = *)                             :: param_name
    integer(kind = c_int)                                         :: l_param_name
    character(kind = c_char, len = *), dimension(n_codes), target :: code_names
    integer, dimension(n_codes), target                           :: l_code_names

    l_param_name = len(param_name)
    do i=1,n_codes
      l_code_names(i) = len(code_names(i))
    end do

    call CWP_Param_reduce_int_cf(op,                &
                                 param_name,        &
                                 l_param_name,      &
                                 res,               &
                                 n_codes,           &
                                 c_loc(code_names), &
                                 c_loc(l_code_names))

  end subroutine CWP_Param_reduce_int

  subroutine CWP_Param_reduce_double(op,         &
                                     param_name, &
                                     res,        &
                                     n_codes,    &
                                     code_names)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                           :: op
    real(c_double), intent(out)                                   :: res
    integer(c_int)                                                :: n_codes, i
    character(kind = c_char, len = *)                             :: param_name
    integer(kind = c_int)                                         :: l_param_name
    character(kind = c_char, len = *), dimension(n_codes), target :: code_names
    integer, dimension(n_codes), target                           :: l_code_names

    l_param_name = len(param_name)
    do i=1,n_codes
      l_code_names(i) = len(code_names(i))
    end do

    call CWP_Param_reduce_double_cf(op,                &
                                    param_name,        &
                                    l_param_name,      &
                                    res,               &
                                    n_codes,           &
                                    c_loc(code_names), &
                                    c_loc(l_code_names))

  end subroutine CWP_Param_reduce_double

  subroutine CWP_Param_reduce_char(op,         &
                                   param_name, &
                                   res,        &
                                   n_codes,    &
                                   code_names)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                           :: op
    character(c_char), pointer                                    :: res(:)
    integer(c_int)                                                :: n_codes, i
    character(kind = c_char, len = *)                             :: param_name
    integer(kind = c_int)                                         :: l_param_name
    character(kind = c_char, len = *), dimension(n_codes), target :: code_names
    integer, dimension(n_codes), target                           :: l_code_names
    type(c_ptr)                                                   :: cptr = C_NULL_PTR
    integer(c_int)                                                :: l_res

    l_param_name = len(param_name)
    do i=1,n_codes
      l_code_names(i) = len(code_names(i))
    end do

    call CWP_Param_reduce_char_cf(op,                &
                                  param_name,        &
                                  l_param_name,      &
                                  cptr,              &
                                  l_res,             &
                                  n_codes,           &
                                  c_loc(code_names), &
                                  c_loc(l_code_names))

    call c_f_pointer(cptr, res, [l_res])

  end subroutine CWP_Param_reduce_char



  !>
  !!
  !! \brief Lock access to local parameters from a distant code.
  !!
  !! \param [in]  code_name  Code to lock
  !!
  !!

  subroutine CWP_Param_lock_(code_name)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name
    integer(kind = c_int)             :: l_code_name

    l_code_name  = len(code_name)

    call CWP_Param_lock_cf(code_name,  &
                          l_code_name)

  end subroutine CWP_Param_lock_


  !>
  !!
  !! \brief Unlock access to local parameters from a distant code.
  !!
  !! \param [in]  code_name  Code to unlock
  !!
  !!

  subroutine CWP_Param_unlock_(code_name)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name
    integer(kind = c_int)             :: l_code_name

    l_code_name  = len(code_name)

    call CWP_Param_unlock_cf(code_name,   &
                             l_code_name)

  end subroutine CWP_Param_unlock_


  !>
  !! \brief Initiate the sending of a data array.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] global_data_id   GlobalData identifier
  !! \param [in] send_data        Pointer to data array
  !!
  !!

  subroutine CWP_Global_data_issend_int(f_local_code_name, &
                                        f_cpl_id,          &
                                        f_global_data_id,  &
                                        send_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: f_local_code_name, f_cpl_id, f_global_data_id
    integer(c_int), pointer       :: send_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: send_stride
    integer(c_int)                :: n_send_entity

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)
    l_global_data_id  = len(f_global_data_id)

    s_data        = 4
    send_stride   = size(send_data, 1)
    n_send_entity = size(send_data, 2)

    call CWP_Global_data_issend_cf(f_local_code_name, &
                                   l_local_code_name, &
                                   f_cpl_id,          &
                                   l_cpl_id,          &
                                   f_global_data_id,  &
                                   l_global_data_id,  &
                                   s_data,            &
                                   send_stride,       &
                                   n_send_entity,     &
                                   c_loc(send_data))

  end subroutine CWP_Global_data_issend_int

  subroutine CWP_Global_data_issend_long(f_local_code_name, &
                                         f_cpl_id,          &
                                         f_global_data_id,  &
                                         send_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: f_local_code_name, f_cpl_id, f_global_data_id
    integer(c_long), pointer      :: send_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: send_stride
    integer(c_int)                :: n_send_entity

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)
    l_global_data_id  = len(f_global_data_id)

    s_data        = 8
    send_stride   = size(send_data, 1)
    n_send_entity = size(send_data, 2)

    call CWP_Global_data_issend_cf(f_local_code_name, &
                                   l_local_code_name, &
                                   f_cpl_id,          &
                                   l_cpl_id,          &
                                   f_global_data_id,  &
                                   l_global_data_id,  &
                                   s_data,            &
                                   send_stride,       &
                                   n_send_entity,     &
                                   c_loc(send_data))

  end subroutine CWP_Global_data_issend_long

  subroutine CWP_Global_data_issend_double(f_local_code_name, &
                                           f_cpl_id,          &
                                           f_global_data_id,  &
                                           send_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: f_local_code_name, f_cpl_id, f_global_data_id
    double precision, pointer     :: send_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: send_stride
    integer(c_int)                :: n_send_entity

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)
    l_global_data_id  = len(f_global_data_id)

    s_data        = 8
    send_stride   = size(send_data, 1)
    n_send_entity = size(send_data, 2)

    call CWP_Global_data_issend_cf(f_local_code_name, &
                                   l_local_code_name, &
                                   f_cpl_id,          &
                                   l_cpl_id,          &
                                   f_global_data_id,  &
                                   l_global_data_id,  &
                                   s_data,            &
                                   send_stride,       &
                                   n_send_entity,     &
                                   c_loc(send_data))

  end subroutine CWP_Global_data_issend_double

  subroutine CWP_Global_data_issend_complex4(f_local_code_name, &
                                             f_cpl_id,          &
                                             f_global_data_id,  &
                                             send_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: f_local_code_name, f_cpl_id, f_global_data_id
    complex(kind = 4), pointer    :: send_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: send_stride
    integer(c_int)                :: n_send_entity

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)
    l_global_data_id  = len(f_global_data_id)

    s_data        = 8
    send_stride   = size(send_data, 1)
    n_send_entity = size(send_data, 2)

    call CWP_Global_data_issend_cf(f_local_code_name, &
                                   l_local_code_name, &
                                   f_cpl_id,          &
                                   l_cpl_id,          &
                                   f_global_data_id,  &
                                   l_global_data_id,  &
                                   s_data,            &
                                   send_stride,       &
                                   n_send_entity,     &
                                   c_loc(send_data))

  end subroutine CWP_Global_data_issend_complex4

  subroutine CWP_Global_data_issend_complex8(f_local_code_name, &
                                             f_cpl_id,          &
                                             f_global_data_id,  &
                                             send_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: f_local_code_name, f_cpl_id, f_global_data_id
    complex(kind = 8), pointer    :: send_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: send_stride
    integer(c_int)                :: n_send_entity

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)
    l_global_data_id  = len(f_global_data_id)

    s_data        = 16
    send_stride   = size(send_data, 1)
    n_send_entity = size(send_data, 2)

    call CWP_Global_data_issend_cf(f_local_code_name, &
                                   l_local_code_name, &
                                   f_cpl_id,          &
                                   l_cpl_id,          &
                                   f_global_data_id,  &
                                   l_global_data_id,  &
                                   s_data,            &
                                   send_stride,       &
                                   n_send_entity,     &
                                   c_loc(send_data))

  end subroutine CWP_Global_data_issend_complex8

  subroutine CWP_Global_data_issend_real4(f_local_code_name, &
                                          f_cpl_id,          &
                                          f_global_data_id,  &
                                          send_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: f_local_code_name, f_cpl_id, f_global_data_id
    real(kind = 4), pointer       :: send_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: send_stride
    integer(c_int)                :: n_send_entity

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)
    l_global_data_id  = len(f_global_data_id)

    s_data        = 4
    send_stride   = size(send_data, 1)
    n_send_entity = size(send_data, 2)

    call CWP_Global_data_issend_cf(f_local_code_name, &
                                   l_local_code_name, &
                                   f_cpl_id,          &
                                   l_cpl_id,          &
                                   f_global_data_id,  &
                                   l_global_data_id,  &
                                   s_data,            &
                                   send_stride,       &
                                   n_send_entity,     &
                                   c_loc(send_data))

  end subroutine CWP_Global_data_issend_real4

  !>
  !! \brief Initiate the reception of a data array.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] global_data_id   GlobalData identifier
  !! \param [in] recv_data        Pointer to data array
  !!
  !!

  subroutine CWP_Global_data_irecv_int(f_local_code_name, &
                                       f_cpl_id,          &
                                       f_global_data_id,  &
                                       recv_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: f_local_code_name, f_cpl_id, f_global_data_id
    integer(c_int), pointer       :: recv_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: recv_stride
    integer(c_int)                :: n_recv_entity

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)
    l_global_data_id  = len(f_global_data_id)

    s_data = 4
    recv_stride   = size(recv_data, 1)
    n_recv_entity = size(recv_data, 2)

    call CWP_Global_data_irecv_cf(f_local_code_name, &
                                  l_local_code_name, &
                                  f_cpl_id,          &
                                  l_cpl_id,          &
                                  f_global_data_id,  &
                                  l_global_data_id,  &
                                  s_data,            &
                                  recv_stride,       &
                                  n_recv_entity,     &
                                  c_loc(recv_data))

  end subroutine CWP_Global_data_irecv_int

  subroutine CWP_Global_data_irecv_long(f_local_code_name, &
                                        f_cpl_id,          &
                                        f_global_data_id,  &
                                        recv_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: f_local_code_name, f_cpl_id, f_global_data_id
    integer(c_long), pointer      :: recv_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: recv_stride
    integer(c_int)                :: n_recv_entity

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)
    l_global_data_id  = len(f_global_data_id)

    s_data = 8
    recv_stride   = size(recv_data, 1)
    n_recv_entity = size(recv_data, 2)

    call CWP_Global_data_irecv_cf(f_local_code_name, &
                                  l_local_code_name, &
                                  f_cpl_id,          &
                                  l_cpl_id,          &
                                  f_global_data_id,  &
                                  l_global_data_id,  &
                                  s_data,            &
                                  recv_stride,       &
                                  n_recv_entity,     &
                                  c_loc(recv_data))

  end subroutine CWP_Global_data_irecv_long

  subroutine CWP_Global_data_irecv_double(f_local_code_name, &
                                          f_cpl_id,          &
                                          f_global_data_id,  &
                                          recv_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: f_local_code_name, f_cpl_id, f_global_data_id
    double precision, pointer     :: recv_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: recv_stride
    integer(c_int)                :: n_recv_entity

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)
    l_global_data_id  = len(f_global_data_id)

    s_data = 8
    recv_stride   = size(recv_data, 1)
    n_recv_entity = size(recv_data, 2)

    call CWP_Global_data_irecv_cf(f_local_code_name, &
                                  l_local_code_name, &
                                  f_cpl_id,          &
                                  l_cpl_id,          &
                                  f_global_data_id,  &
                                  l_global_data_id,  &
                                  s_data,            &
                                  recv_stride,       &
                                  n_recv_entity,     &
                                  c_loc(recv_data))

  end subroutine CWP_Global_data_irecv_double

  subroutine CWP_Global_data_irecv_complex4(f_local_code_name, &
                                            f_cpl_id,          &
                                            f_global_data_id,  &
                                            recv_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: f_local_code_name, f_cpl_id, f_global_data_id
    complex(kind=4), pointer      :: recv_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: recv_stride
    integer(c_int)                :: n_recv_entity

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)
    l_global_data_id  = len(f_global_data_id)

    s_data = 8
    recv_stride   = size(recv_data, 1)
    n_recv_entity = size(recv_data, 2)

    call CWP_Global_data_irecv_cf(f_local_code_name, &
                                  l_local_code_name, &
                                  f_cpl_id,          &
                                  l_cpl_id,          &
                                  f_global_data_id,  &
                                  l_global_data_id,  &
                                  s_data,            &
                                  recv_stride,       &
                                  n_recv_entity,     &
                                  c_loc(recv_data))

  end subroutine CWP_Global_data_irecv_complex4

  subroutine CWP_Global_data_irecv_complex8(f_local_code_name, &
                                            f_cpl_id,          &
                                            f_global_data_id,  &
                                            recv_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: f_local_code_name, f_cpl_id, f_global_data_id
    complex(kind=8), pointer      :: recv_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: recv_stride
    integer(c_int)                :: n_recv_entity

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)
    l_global_data_id  = len(f_global_data_id)

    s_data = 16
    recv_stride   = size(recv_data, 1)
    n_recv_entity = size(recv_data, 2)

    call CWP_Global_data_irecv_cf(f_local_code_name, &
                                  l_local_code_name, &
                                  f_cpl_id,          &
                                  l_cpl_id,          &
                                  f_global_data_id,  &
                                  l_global_data_id,  &
                                  s_data,            &
                                  recv_stride,       &
                                  n_recv_entity,     &
                                  c_loc(recv_data))

  end subroutine CWP_Global_data_irecv_complex8

  subroutine CWP_Global_data_irecv_real4(f_local_code_name, &
                                         f_cpl_id,          &
                                         f_global_data_id,  &
                                         recv_data)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: f_local_code_name, f_cpl_id, f_global_data_id
    real(kind=4), pointer         :: recv_data(:,:)

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id
    integer(c_long)               :: s_data
    integer(c_int)                :: recv_stride
    integer(c_int)                :: n_recv_entity

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)
    l_global_data_id  = len(f_global_data_id)

    s_data = 8
    recv_stride   = size(recv_data, 1)
    n_recv_entity = size(recv_data, 2)

    call CWP_Global_data_irecv_cf(f_local_code_name, &
                                  l_local_code_name, &
                                  f_cpl_id,          &
                                  l_cpl_id,          &
                                  f_global_data_id,  &
                                  l_global_data_id,  &
                                  s_data,            &
                                  recv_stride,       &
                                  n_recv_entity,     &
                                  c_loc(recv_data))

  end subroutine CWP_Global_data_irecv_real4


  !>
  !! \brief Finalize the sending of a data array.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] global_data_id   GlobalData identifier
  !!
  !!

  subroutine CWP_Global_data_wait_issend_(f_local_code_name, &
                                          f_cpl_id,          &
                                          f_global_data_id)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: f_local_code_name, f_cpl_id, f_global_data_id

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)
    l_global_data_id  = len(f_global_data_id)

    call CWP_Global_data_wait_issend_cf(f_local_code_name, &
                                        l_local_code_name, &
                                        f_cpl_id,          &
                                        l_cpl_id,          &
                                        f_global_data_id,  &
                                        l_global_data_id)

  end subroutine CWP_Global_data_wait_issend_


  !>
  !! \brief Finalize the reception of a data array.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] global_data_id   GlobalData identifier
  !!
  !!

  subroutine CWP_Global_data_wait_irecv_(f_local_code_name, &
                                         f_cpl_id,          &
                                         f_global_data_id)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: f_local_code_name, f_cpl_id, f_global_data_id

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_global_data_id

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)
    l_global_data_id  = len(f_global_data_id)

    call CWP_Global_data_wait_irecv_cf(f_local_code_name, &
                                       l_local_code_name, &
                                       f_cpl_id,          &
                                       l_cpl_id,          &
                                       f_global_data_id,  &
                                       l_global_data_id)

  end subroutine CWP_Global_data_wait_irecv_


  !>
  !! \brief Create partitioned data exchange object
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] part_data_id     PartData identifier
  !! \param [in] exch_type        Exchange type
  !! \param [in] gnum_elt         Global ids
  !! \param [in] n_elt            Number of elements in partitions (size = \p n_part)
  !! \param [in] n_part           Number of partitions
  !!
  !!

  subroutine CWP_Part_data_create_(f_local_code_name, &
                                   f_cpl_id,          &
                                   f_part_data_id,    &
                                   exch_type,         &
                                   gnum_elt,          &
                                   n_elt,             &
                                   n_part)
    use, intrinsic :: iso_c_binding
    use pdm_pointer_array
    implicit none
    character(kind=c_char, len=*)     :: f_local_code_name, f_cpl_id, f_part_data_id
    integer(c_int), intent(in)        :: exch_type
    type(PDM_pointer_array_t), target :: gnum_elt
    integer(c_int), pointer           :: n_elt(:)
    integer,        intent(in)        :: n_part

    integer(c_int)                    :: l_local_code_name, l_cpl_id, l_part_data_id

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)
    l_part_data_id    = len(f_part_data_id)

    call CWP_Part_data_create_cf(f_local_code_name,    &
                                 l_local_code_name,    &
                                 f_cpl_id,             &
                                 l_cpl_id,             &
                                 f_part_data_id,       &
                                 l_part_data_id,       &
                                 exch_type,            &
                                 c_loc(gnum_elt%cptr), &
                                 c_loc(n_elt),         &
                                 n_part)

  end subroutine CWP_Part_data_create_


  !>
  !! \brief Delete partitioned data exchange object
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] part_data_id     PartData identifier
  !! \param [in] exch_type        Exchange type
  !!
  !!

  subroutine CWP_Part_data_del_(f_local_code_name, &
                                f_cpl_id,          &
                                f_part_data_id,    &
                                exch_type)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: f_local_code_name, f_cpl_id, f_part_data_id
    integer(c_int), intent(in)    :: exch_type

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_part_data_id

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)
    l_part_data_id    = len(f_part_data_id)

    call CWP_Part_data_del_cf(f_local_code_name, &
                              l_local_code_name, &
                              f_cpl_id,          &
                              l_cpl_id,          &
                              f_part_data_id,    &
                              l_part_data_id,    &
                              exch_type)

  end subroutine CWP_Part_data_del_


  !>
  !! \brief Send a data array.
  !!
  !! \param [in]  local_code_name      Local code name
  !! \param [in]  cpl_id               Coupling identifier
  !! \param [in]  part_data_id         PartData identifier
  !! \param [in]  n_components         Number of components
  !! \param [in]  part1_to_part2_data  Pointer to data array to send
  !! \param [out] request              MPI request
  !!
  !!

  subroutine CWP_Part_data_issend_(local_code_name,     &
                                   cpl_id,              &
                                   part_data_id,        &
                                   n_components,        &
                                   part1_to_part2_data)
    use, intrinsic :: iso_c_binding
    use pdm_pointer_array
    implicit none

    character(kind=c_char, len=*)     :: local_code_name, cpl_id, part_data_id
    integer(c_int), intent(in)        :: n_components
    type(PDM_pointer_array_t), target :: part1_to_part2_data

    integer(c_int)                    :: l_local_code_name, l_cpl_id, l_part_data_id
    integer(c_long)                   :: s_data

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_part_data_id    = len(part_data_id)

    s_data = part1_to_part2_data%s_data

    call CWP_Part_data_issend_cf(local_code_name,                 &
                                 l_local_code_name,               &
                                 cpl_id,                          &
                                 l_cpl_id,                        &
                                 part_data_id,                    &
                                 l_part_data_id,                  &
                                 s_data,                          &
                                 n_components,                    &
                                 c_loc(part1_to_part2_data%cptr))

  end subroutine CWP_Part_data_issend_


  !>
  !! \brief Receive a data array.
  !!
  !! \param [in]    local_code_name      Local code name
  !! \param [in]    cpl_id               Coupling identifier
  !! \param [in]    part_data_id         PartData identifier
  !! \param [in]    n_components         Number of components
  !! \param [inout] part2_data           Pointer to data array to receive
  !!
  !!

  subroutine CWP_Part_data_irecv_(local_code_name, &
                                  cpl_id,          &
                                  part_data_id,    &
                                  n_components,    &
                                  part2_data)
    use, intrinsic :: iso_c_binding
    use pdm_pointer_array
    implicit none

    character(kind=c_char, len=*)     :: local_code_name, cpl_id, part_data_id
    integer(c_int), intent(in)        :: n_components
    type(PDM_pointer_array_t), target :: part2_data

    integer(c_int)                    :: l_local_code_name, l_cpl_id, l_part_data_id
    integer(c_long)                   :: s_data
    integer(c_int)                    :: n_part, n_ref, i

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_part_data_id    = len(part_data_id)

    s_data = part2_data%s_data

    call CWP_Part_data_irecv_cf(local_code_name,        &
                                l_local_code_name,      &
                                cpl_id,                 &
                                l_cpl_id,               &
                                part_data_id,           &
                                l_part_data_id,         &
                                s_data,                 &
                                n_components,           &
                                c_loc(part2_data%cptr))

    call CWP_Part_data_n_part_get_cf(local_code_name,   &
                                     l_local_code_name, &
                                     cpl_id,            &
                                     l_cpl_id,          &
                                     part_data_id,      &
                                     l_part_data_id,    &
                                     n_part)

    ! Compute lengths
    do i = 1, n_part
      call CWP_Part_data_n_ref_get_cf(local_code_name,   &
                                      l_local_code_name, &
                                      cpl_id,            &
                                      l_cpl_id,          &
                                      part_data_id,      &
                                      l_part_data_id,    &
                                      i-1,               &
                                      n_ref)

      part2_data%length(i) = n_components * n_ref
    enddo

  end subroutine CWP_Part_data_irecv_


  !>
  !! \brief Wait of send a data array.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] part_data_id     PartData identifier
  !!
  !!

  subroutine CWP_Part_data_wait_issend_(local_code_name, &
                                        cpl_id,          &
                                        part_data_id)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name, cpl_id, part_data_id

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_part_data_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_part_data_id    = len(part_data_id)

    call CWP_Part_data_wait_issend_cf(local_code_name,   &
                                      l_local_code_name, &
                                      cpl_id,            &
                                      l_cpl_id,          &
                                      part_data_id,      &
                                      l_part_data_id)

  end subroutine CWP_Part_data_wait_issend_


  !>
  !! \brief Wait of receive a data array.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] part_data_id     PartData identifier
  !!
  !!

  subroutine CWP_Part_data_wait_irecv_(local_code_name, &
                                       cpl_id,          &
                                       part_data_id)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: local_code_name, cpl_id, part_data_id

    integer(c_int)                :: l_local_code_name, l_cpl_id, l_part_data_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_part_data_id    = len(part_data_id)

    call CWP_Part_data_wait_irecv_cf(local_code_name,   &
                                     l_local_code_name, &
                                     cpl_id,            &
                                     l_cpl_id,          &
                                     part_data_id,      &
                                     l_part_data_id)

  end subroutine CWP_Part_data_wait_irecv_




  function CWP_Cpl_spatial_interp_algo_get_(f_local_code_name, &
                                            f_cpl_id) result(algo)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char, len=*) :: f_local_code_name, f_cpl_id
    integer(c_int)                :: l_local_code_name, l_cpl_id
    integer(c_int)                :: algo

    l_local_code_name = len(f_local_code_name)
    l_cpl_id          = len(f_cpl_id)

    algo = CWP_Cpl_spatial_interp_algo_get_cf(f_local_code_name, &
                                              l_local_code_name, &
                                              f_cpl_id,          &
                                              l_cpl_id)

  end function CWP_Cpl_spatial_interp_algo_get_


end module cwp
