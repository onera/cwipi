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
                CWP_COMM_PAR_WITHOUT_PART, &
                CWP_COMM_SEQ
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

    interface CWP_State_update ; module procedure &
        CWP_State_update_
    end interface CWP_State_update

    interface CWP_Time_update ; module procedure &
        CWP_Time_update_
    end interface CWP_Time_update

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

    interface CWP_Cpl_Del ; module procedure &
        CWP_Cpl_Del_
    end interface CWP_Cpl_Del

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

    interface CWP_Field_n_component_get ; module procedure &
        CWP_Field_n_component_get_
    end interface CWP_Field_n_component_get

    interface CWP_Field_target_dof_location_get ; module procedure &
        CWP_Field_target_dof_location_get_
    end interface CWP_Field_target_dof_location_get

    interface CWP_Field_storage_get ; module procedure &
        CWP_Field_storage_get_
    end interface CWP_Field_storage_get

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

    interface CWP_Interp_from_location_unset ; module procedure &
        CWP_Interp_from_location_unset_
    end interface CWP_Interp_from_location_unset

!!$    interface CWP_Interp_from_location_set ; module procedure &
!!$        CWP_Interp_from_location_set_
!!$    end interface CWP_Interp_from_location_set

    interface CWP_Param_add_int ; module procedure &
        CWP_Param_add_int_
    end interface CWP_Param_add_int

    interface CWP_Param_add_double ; module procedure &
        CWP_Param_add_double_
    end interface CWP_Param_add_double

    interface CWP_Param_add_char ; module procedure &
        CWP_Param_add_char_
    end interface CWP_Param_add_char

    interface CWP_Param_set_int ; module procedure &
        CWP_Param_set_int_
    end interface CWP_Param_set_int

    interface CWP_Param_set_double ; module procedure &
        CWP_Param_set_double_
    end interface CWP_Param_set_double

    interface CWP_Param_set_char ; module procedure &
        CWP_Param_set_char_
    end interface CWP_Param_set_char

    interface CWP_Param_del ; module procedure &
        CWP_Param_del_
    end interface CWP_Param_del

    interface CWP_Param_n_get ; module procedure &
        CWP_Param_n_get_
    end interface CWP_Param_n_get

    interface CWP_Param_is ; module procedure &
        CWP_Param_is_
    end interface CWP_Param_is

    interface CWP_Param_get ; module procedure &
        CWP_Param_get_
    end interface CWP_Param_get

    interface CWP_Param_reduce ; module procedure &
        CWP_Param_reduce_
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

  !
  ! Private

  private :: c_f_char_array,&
             CWP_Init_ ,&
             CWP_State_update_ ,&
             CWP_Time_update_ ,&
             CWP_User_structure_set_ ,&
             CWP_User_structure_get_ ,&
             CWP_Output_file_set_,&
             CWP_State_get_ ,&
             CWP_Cpl_create_ ,&
             CWP_Cpl_Del_ ,&
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
             CWP_Field_n_component_get_ ,&
             CWP_Field_target_dof_location_get_ ,&
             CWP_Field_storage_get_ ,&
             CWP_Field_del_ ,&
             CWP_Field_issend_ ,&
             CWP_Field_irecv_ ,&
             CWP_Field_wait_issend_ ,&
             CWP_Field_wait_irecv_ ,&
             CWP_Interp_from_location_unset_ ,&
!!$             CWP_Interp_from_location_set_ ,&
             CWP_Param_add_int_ ,&
             CWP_Param_add_double_ ,&
             CWP_Param_add_char_ ,&
             CWP_Param_set_int_ ,&
             CWP_Param_set_double_ ,&
             CWP_Param_set_char_ ,&
             CWP_Param_del_ ,&
             CWP_Param_n_get_ ,&
             CWP_Param_is_ ,&
             CWP_Param_get_ ,&
             CWP_Param_reduce_,&
             CWP_Param_lock_ ,&
             CWP_Param_unlock_,&
             CWP_Codes_list_get_,&
             CWP_Loc_codes_list_get_,&
             CWP_Param_list_get_

    interface
      subroutine CWP_Init_cf(fcomm, n_code, code_names, l_code_names, is_active_rank, time_init, intra_comms) &
              bind(c, name = 'CWP_Init_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int), value             :: fcomm
        integer(c_int), value             :: n_code
        type(c_ptr),    value             :: code_names
        type(c_ptr),    value             :: l_code_names
        integer(c_int), dimension(n_code) :: is_active_rank
        real(c_double), dimension(n_code) :: time_init
        type(c_ptr),    value             :: intra_comms
      end subroutine CWP_Init_cf

      subroutine CWP_State_update_cf(local_code_name, l_local_code_name, state) &
        bind(c, name='CWP_State_update_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name
        integer(c_int), value             :: l_local_code_name, state
      end subroutine CWP_State_update_cf

      subroutine CWP_Time_update_cf(local_code_name, l_local_code_name, current_time) &
        bind(c, name='CWP_Time_update_cf')
        use, intrinsic :: iso_c_binding
        implicit none
        character(kind = c_char, len = 1) :: local_code_name
        integer(c_int), value             :: l_local_code_name
        real(c_double), value             :: current_time
      end subroutine CWP_Time_update_cf

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

      subroutine CWP_Interp_from_location_fortran_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
              src_field_id, l_src_field_id, ptInterpolationFct) &
              bind(c, name = 'CWP_Interp_from_location_fortran_set_cf')
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
      end subroutine CWP_Interp_from_location_fortran_set_cf

      subroutine CWP_Interp_from_location_fortran_unset_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
              src_field_id, l_src_field_id) &
              bind(c, name = 'CWP_Interp_from_location_fortran_unset_cf')
        use, intrinsic :: iso_c_binding
        implicit none

        character(kind = c_char, len = 1) :: local_code_name, cpl_id, src_field_id
        integer(kind = c_int), value :: l_local_code_name, l_cpl_id, l_src_field_id
      end subroutine CWP_Interp_from_location_fortran_unset_cf

      subroutine CWP_Spatial_interp_property_set_cf(local_code_name,   &
                                                    l_local_code_name, &
                                                    cpl_id,            &
                                                    l_cpl_id,          &
                                                    property_name,     &
                                                    l_property_name,   &
                                                    property_type,     &
                                                    l_property_type,   &
                                                    property_value,    &
                                                    l_property_value)  &
      bind(c, name = 'CWP_Spatial_interp_property_set_cf')
        use, intrinsic :: iso_c_binding
        character(kind = c_char, len = 1) :: local_code_name, cpl_id, property_name, property_type, property_value
        integer(kind = c_int), value      :: l_local_code_name, l_cpl_id, l_property_name, l_property_type, l_property_value
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



    subroutine CWP_Param_add_cf(local_code_name,   &
                                l_local_code_name, &
                                param_name,        &
                                l_param_name,      &
                                data_type,         &
                                initial_value)     &
      bind(c, name='CWP_Param_add_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      integer(c_int), value :: data_type
      type(c_ptr),    value :: initial_value
    end subroutine CWP_Param_add_cf


    subroutine CWP_Param_set_cf(local_code_name,   &
                                l_local_code_name, &
                                param_name,        &
                                l_param_name,      &
                                data_type,         &
                                value)             &
      bind(c, name='CWP_Param_set_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      integer(c_int), value :: data_type
      type(c_ptr),    value :: value
    end subroutine CWP_Param_set_cf


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


    subroutine CWP_Param_get_cf(local_code_name,   &
                                l_local_code_name, &
                                param_name,        &
                                l_param_name,      &
                                data_type,         &
                                value)             &
      bind(c, name='CWP_Param_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, param_name
      integer(c_int), value :: l_local_code_name, l_param_name
      integer(c_int), value :: data_type
      type(c_ptr)           :: value
    end subroutine CWP_Param_get_cf

    subroutine CWP_Param_reduce_cf(op,           &
                                   param_name,   &
                                   l_param_name, &
                                   data_type,    &
                                   res,          &
                                   n_codes,      &
                                   code_names,   &
                                   l_code_names) &

      bind(c, name='CWP_Param_reduce_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value             :: op, data_type, n_codes
      type(c_ptr)                       :: res
      character(kind = c_char, len = 1) :: param_name
      integer(c_int), value             :: l_param_name
      type(c_ptr),    value             :: code_names
      type(c_ptr),    value             :: l_code_names
    end subroutine CWP_Param_reduce_cf


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


    function CWP_Field_n_component_get_cf(local_code_name,   &
                                          l_local_code_name, &
                                          cpl_id,            &
                                          l_cpl_id,          &
                                          field_id,          &
                                          l_field_id)        &
      result (n_component)                                   &
      bind(c, name='CWP_Field_n_component_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
      integer(c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      integer(c_int)        :: n_component
    end function CWP_Field_n_component_get_cf


    function CWP_Field_target_dof_location_get_cf(local_code_name,   &
                                                  l_local_code_name, &
                                                  cpl_id,            &
                                                  l_cpl_id,          &
                                                  field_id,          &
                                                  l_field_id)        &
      result (dof_location)                                          &
      bind(c, name='CWP_Field_target_dof_location_get_cf')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind = c_char, len = 1) :: local_code_name, cpl_id, field_id
      integer(c_int), value :: l_local_code_name, l_cpl_id, l_field_id
      integer(c_int)        :: dof_location
    end function CWP_Field_target_dof_location_get_cf


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

  end interface

contains

  !> convert an array of char * in c to an array
  !  of characters of maximum size 256 in FORTRAN

  subroutine c_f_char_array(c_char_array, &
                            c_size_array, &
                            n_chars, &
                            f_char_array)

    use, intrinsic :: iso_c_binding
    implicit none

    type(c_ptr) :: c_char_array, c_size_array
    type(c_ptr), pointer :: fptr2(:) => null()
    character(c_char), pointer :: fptr(:) => null()
    integer(c_int), pointer     :: f_size_array(:)
    character(256), allocatable :: f_char_array(:)
    integer(c_int) :: i, n_chars, strlen

    call c_f_pointer(c_char_array, fptr2, [n_chars])
    call c_f_pointer(c_size_array, f_size_array, [n_chars])

    allocate(f_char_array(n_chars))
    do i = 1, n_chars
      strlen      = f_size_array(i)
      call c_f_pointer(fptr2(i), fptr, [strlen])
      f_char_array(i) = transfer(fptr(1:strlen), f_char_array(i))
      f_char_array(i) = f_char_array(i)(1:strlen)//char(0)
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
  !! \param [in]  is_active_rank Is current rank have to be used by CWIPI (size = \p n_code)
  !! \param [in]  time_init      Initial time (size = \p n_code)
  !! \param [out] intra_comms    MPI intra communicators of each code (size = \p n_code)
  !!
  !!

  subroutine CWP_Init_ (fcomm,          &
                       n_code,         &
                       code_names,     &
                       is_active_rank, &
                       time_init,      &
                       intra_comms)

    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int) :: fcomm
    integer(c_int), intent(in) :: n_code
    character(kind = c_char, len = *), dimension(n_code), target :: code_names
    integer(c_int), dimension(n_code) :: is_active_rank
    real(c_double), dimension(n_code) :: time_init
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
                     time_init,           &
                     c_loc(intra_comms))

  end subroutine CWP_Init_


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
  !! \brief Update code time.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in]  current_time Current time
  !!
  !!

  subroutine CWP_Time_update_(local_code_name, &
                             current_time)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    double precision, intent(in)      :: current_time
    integer(c_int)                    :: l_local_code_name

    l_local_code_name = len(local_code_name)

    call CWP_Time_update_cf(local_code_name,   &
                            l_local_code_name, &
                            current_time)

  end subroutine CWP_Time_update_

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
  !! \brief Writing output to fortran file (shared by fortran and C code).
  !!
  !! This function set the file fortran logical unit for writing output.
  !!
  !!  \param [in]  iunit        File fortan logical unit
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
  !! \param [in]  property_type    Type of the property ("double" or "int")
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
    character(kind = c_char, len = *) :: property_type
    character(kind = c_char, len = *) :: property_value

    call CWP_Spatial_interp_property_set_cf(local_code_name,      &
                                            len(local_code_name), &
                                            cpl_id,               &
                                            len(cpl_id),          &
                                            property_name,        &
                                            len(property_name),   &
                                            property_type,        &
                                            len(property_type),   &
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
    double precision, dimension(:), pointer :: coord
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
  !! \param [in]  coord            Coordinates (size = 3    !! \p n_pts)
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
    integer(c_int) :: i_part, n_pts
    double precision, dimension(:), pointer :: coord
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
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !! \param [in]  i_part           Partition identifier
  !! \param [in]  block_id         Block identifier
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

  subroutine CWP_Mesh_interf_f_poly_block_set_( local_code_name, &
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
  !! \param [in]  local_code_name  Local code name
  !! \param [in]  cpl_id           Coupling identifier
  !! \param [in]  i_part           Current partition
  !! \param [in]  block_id         Block identifier
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
  !! \param [in]  edge_vtx_idx      Polyhedron face to vertex index
  !!                                (\p edge_vtx_idx[0] = 0 and
  !!                                 size = \p n_edges + 1)
  !! \param [in]  edge_vtx          Face to vertex connectivity
  !!                                (size = \p edge_vtx_idx[\p n_edges])
  !! \param [in]  global_num        Pointer to parent element number (or NULL)
  !!

  subroutine CWP_Mesh_interf_from_faceedge_set_ (local_code_name, &
                                                cpl_id,          &
                                                i_part,          &
                                                n_faces,         &
                                                face_edge_idx,   &
                                                face_edge,       &
                                                n_edges,         &
                                                edge_vtx_idx,    &
                                                edge_vtx,        &
                                                global_num)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id
    integer(c_int) :: i_part, n_faces, n_edges
    integer(c_int), dimension(:), pointer :: face_edge_idx, face_edge, edge_vtx_idx, edge_vtx
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
                                               c_loc(edge_vtx_idx), &
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

  subroutine CWP_Field_create_ (local_code_name,      &
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
  !! \param [in] data_type         Choice if data is setted for the source or the target
  !! \param [in] data              Storage array (Mapping)
  !!

  subroutine CWP_Field_data_set_ (local_code_name, &
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
  end subroutine CWP_Field_data_set_


  !>
  !!
  !! \brief Get number of field components.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in] cpl_id           Coupling identifier
  !! \param [in] field_id         Field identifier
  !!
  !! \return                      number of field components
  !!
  !!

  function CWP_Field_n_component_get_(local_code_name, &
                                     cpl_id,          &
                                     field_id)        &
    result (n_component)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, field_id
    integer(c_int) :: n_component
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_field_id

    l_local_code_name = len(local_code_name)
    l_cpl_id          = len(cpl_id)
    l_field_id        = len(field_id)

    n_component = CWP_Field_n_component_get_cf(local_code_name,   &
                                               l_local_code_name, &
                                               cpl_id,            &
                                               l_cpl_id,          &
                                               field_id,          &
                                               l_field_id)

  end function CWP_Field_n_component_get_


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

  function CWP_Field_target_dof_location_get_(local_code_name, &
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

    dof_location = CWP_Field_target_dof_location_get_cf(local_code_name,   &
                                                        l_local_code_name, &
                                                        cpl_id,            &
                                                        l_cpl_id,          &
                                                        field_id,          &
                                                        l_field_id)

  end function CWP_Field_target_dof_location_get_


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
  !! \brief Delete a field.
  !!
  !! \param [in] local_code_name Local code name
  !! \param [in]  cpl_id         Coupling identifier
  !! \param [in]  field_id       Field identifier
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
  !!        nonblocking communications.
  !!
  !! This function is independant of \ref CWP_Time_exch_t mode. The user has to
  !! manually check the consistency of the exchanges.
  !!
  !! \param [in] local_code_name  Local code name
  !! \param [in]  cpl_id          Coupling identifier
  !! \param [in]  src_field_id    Source field id
  !!

  subroutine CWP_Field_issend_ (local_code_name, &
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
  end subroutine CWP_Field_issend_


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

  subroutine CWP_Field_irecv_ (local_code_name, &
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
  !! \param [in] src_field_id     Source field id
  !!

  subroutine CWP_Field_wait_issend_ (local_code_name, &
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
  !!  \brief Unsetting of an user interpolation.
  !!
  !!  \param [in] local_code_name  Local code name
  !!  \param [in] cpl_id           Coupling identifier
  !!  \param [in] src_field_id     Source field id
  !!
  !!

  subroutine CWP_Interp_from_location_unset_ (local_code_name, &
                                             cpl_id, &
                                             src_field_id)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name, cpl_id, src_field_id
    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_src_field_id
    l_local_code_name = len(local_code_name)
    l_cpl_id = len(cpl_id)
    l_src_field_id = len(src_field_id)

    call CWP_Interp_from_location_fortran_unset_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
            & src_field_id, l_src_field_id)

  end subroutine CWP_Interp_from_location_unset_

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

!!$  subroutine CWP_Interp_from_location_set_(local_code_name, &
!!$                                          cpl_id, &
!!$                                          src_field_id, &
!!$                                          ptInterpolationFct)
!!$
!!$    use, intrinsic :: iso_c_binding
!!$    implicit none
!!$
!!$    interface
!!$      subroutine ptInterpolationFct ( &
!!$              interface_type, &
!!$              n_src_vtcs, &
!!$              n_src_std_elts, &
!!$              n_tgt_pts, &
!!$              src_vtcs_coords, &
!!$              src_connec_idx, &
!!$              src_connec, &
!!$              tgt_pts_coords, &
!!$              tgt_pts_target_location, &
!!$              tgt_pts_dist, &
!!$              tgt_pts_bary_coords_idx, &
!!$              tgt_pts_bary_coords, &
!!$              stride, &
!!$              src_field_location, &
!!$              src_field, &
!!$              tgt_field_location, &
!!$              tgt_field &
!!$              )
!!$        use, intrinsic :: iso_c_binding
!!$        implicit none
!!$
!!$        integer(kind = c_int)               :: interface_type
!!$        integer(kind = c_int)               :: n_src_vtcs
!!$        integer(kind = c_int)               :: n_src_std_elts
!!$        integer(kind = c_int)               :: n_tgt_pts
!!$        real(kind = c_double), dimension(*) :: src_vtcs_coords
!!$        integer(kind = c_int), dimension(*) :: src_connec_idx
!!$        integer(kind = c_int), dimension(*) :: src_connec
!!$        real(kind = c_double), dimension(*) :: tgt_pts_coords
!!$        integer(kind = c_int), dimension(*) :: tgt_pts_target_location
!!$        real(kind = c_double), dimension(*) :: tgt_pts_dist
!!$        integer(kind = c_int), dimension(*) :: tgt_pts_bary_coords_idx
!!$        real(kind = c_double), dimension(*) :: tgt_pts_bary_coords
!!$        integer(kind = c_int)               :: stride
!!$        integer(kind = c_int)               :: src_field_location
!!$        real(kind = c_double), dimension(*) :: src_field
!!$        integer(kind = c_int)               :: tgt_field_location
!!$        real(kind = c_double), dimension(*) :: tgt_field
!!$      end subroutine ptInterpolationFct
!!$    end interface
!!$
!!$    character(kind = c_char, len = *) :: local_code_name, cpl_id, src_field_id
!!$    integer(kind = c_int) :: l_local_code_name, l_cpl_id, l_src_field_id
!!$    l_local_code_name = len(local_code_name)
!!$    l_cpl_id = len(cpl_id)
!!$    l_src_field_id = len(src_field_id)
!!$
!!$    call CWP_Interp_from_location_fortran_set_cf(local_code_name, l_local_code_name, cpl_id, l_cpl_id, &
!!$            & src_field_id, l_src_field_id, ptInterpolationFct)
!!$  end subroutine CWP_Interp_from_location_set_

! /*----------------------------------------------------------------------------*
!  * Functions about all code parameters                                        *
!  *----------------------------------------------------------------------------*/


  !>
  !!
  !! \brief Add a new parameter and intialize it.
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
    integer, intent(in)               :: initial_value
    integer(kind = c_int)             :: l_local_code_name, l_param_name

    integer, dimension(1), target :: cvalue
    integer, pointer               :: cptrvalue(:)
    integer                       :: data_type

    data_type = CWP_INT

    cvalue(1) = initial_value

    cptrvalue => cvalue

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)

    call CWP_Param_add_cf(local_code_name,   &
                          l_local_code_name, &
                          param_name,        &
                          l_param_name,      &
                          data_type,         &
                          c_loc(cptrvalue))

  end subroutine CWP_Param_add_int_


  subroutine CWP_Param_add_double_(local_code_name, &
                                  param_name,      &
                                  initial_value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    real(kind = 8), intent(in)               :: initial_value
    integer(kind = c_int)             :: l_local_code_name, l_param_name

    real(kind = 8), dimension(1), target :: cvalue
    real(kind = 8), pointer              :: cptrvalue(:)
    integer                              :: data_type

    data_type = CWP_DOUBLE

    cvalue(1) = initial_value

    cptrvalue => cvalue

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)

    call CWP_Param_add_cf(local_code_name,   &
                          l_local_code_name, &
                          param_name,        &
                          l_param_name,      &
                          data_type,         &
                          c_loc(cptrvalue))

  end subroutine CWP_Param_add_double_

  subroutine CWP_Param_add_char_(local_code_name, &
                                param_name,      &
                                initial_value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    character(len=*), intent(in)      :: initial_value
    integer(kind = c_int)             :: l_local_code_name, l_param_name

    character(len=512), dimension(1), target :: cvalue
    character(len=512), pointer            :: cptrvalue(:)
    integer                                :: data_type

    data_type = CWP_CHAR

    cvalue(1) = initial_value

    cptrvalue => cvalue

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)

    call CWP_Param_add_cf(local_code_name,   &
                          l_local_code_name, &
                          param_name,        &
                          l_param_name,      &
                          data_type,         &
                          c_loc(cptrvalue))

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
    integer, intent(in)               :: value
    integer(kind = c_int)             :: l_local_code_name, l_param_name

    integer, dimension(1), target :: cvalue
    integer, pointer               :: cptrvalue(:)

    integer                        :: data_type = CWP_INT

    cvalue(1) = value

    cptrvalue => cvalue

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)

    call CWP_Param_set_cf(local_code_name,   &
                          l_local_code_name, &
                          param_name,        &
                          l_param_name,      &
                          data_type,         &
                          c_loc(cptrvalue))

  end subroutine CWP_Param_set_int_

  subroutine CWP_Param_set_double_(local_code_name, &
                                  param_name,      &
                                  value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    real (kind = 8), intent(in)       :: value
    integer(kind = c_int)             :: l_local_code_name, l_param_name

    real(kind = 8), dimension(1), target :: cvalue
    real(kind = 8), pointer               :: cptrvalue(:)

    integer                        :: data_type = CWP_DOUBLE

    cvalue(1) = value

    cptrvalue => cvalue

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)

    call CWP_Param_set_cf(local_code_name,   &
                          l_local_code_name, &
                          param_name,        &
                          l_param_name,      &
                          data_type,         &
                          c_loc(cptrvalue))

  end subroutine CWP_Param_set_double_

  subroutine CWP_Param_set_char_(local_code_name, &
                                param_name,      &
                                value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: local_code_name
    character(kind = c_char, len = *) :: param_name
    character(len=*)                  :: value
    integer(kind = c_int)             :: l_local_code_name, l_param_name

    character(len = 512), dimension(1), target :: cvalue
    character(len = 512), pointer               :: cptrvalue(:)

    integer                        :: data_type = CWP_CHAR

    cvalue(1) = value

    cptrvalue => cvalue

    l_local_code_name = len(local_code_name)
    l_param_name      = len(param_name)

    call CWP_Param_set_cf(local_code_name,   &
                          l_local_code_name, &
                          param_name,        &
                          l_param_name,      &
                          data_type,         &
                          c_loc(cptrvalue))

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
  !! \brief eturn the list of parameters for the code \p code_name.
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

    call c_f_char_array(c_param_names, c_param_sizes, n_param, f_param_names)

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

  subroutine CWP_Param_get_(code_name,  &
                           param_name, &
                           data_type,  &
                           value)

    use, intrinsic :: iso_c_binding
    implicit none

    character(kind = c_char, len = *) :: code_name
    character(kind = c_char, len = *) :: param_name
    integer, intent(in)               :: data_type
    type(c_ptr)                       :: value
    integer(kind = c_int)             :: l_code_name
    integer(kind = c_int)             :: l_param_name

    l_code_name  = len(code_name)
    l_param_name = len(param_name)

    call CWP_Param_get_cf(code_name,    &
                          l_code_name,  &
                          param_name,   &
                          l_param_name, &
                          data_type,    &
                          value)

  end subroutine CWP_Param_get_

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

  subroutine CWP_Param_reduce_(op,         &
                               param_name, &
                               data_type,  &
                               res,        &
                               n_codes,    &
                               code_names)

    use, intrinsic :: iso_c_binding
    implicit none

    integer, intent(in)                                           :: op, data_type
    type(c_ptr)                                                   :: res
    integer(c_int)                                                :: n_codes, i
    character(kind = c_char, len = *)                             :: param_name
    integer(kind = c_int)                                         :: l_param_name
    character(kind = c_char, len = *), dimension(n_codes), target :: code_names
    integer, dimension(n_codes), target                           :: l_code_names

    l_param_name = len(param_name)
    do i=1,n_codes
      l_code_names(i) = len(code_names(i))
    end do

    call CWP_Param_reduce_cf(op,                &
                             param_name,        &
                             l_param_name,      &
                             data_type,         &
                             res,               &
                             n_codes,           &
                             c_loc(code_names), &
                             c_loc(l_code_names))

  end subroutine CWP_Param_reduce_

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


end module cwp
