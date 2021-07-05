#ifndef CWP_CF_H_
#define CWP_CF_H_

/*
  This file is part of the CWIPI library.

  Copyright (C) 2011-2017  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include "cwp.h"

void CWP_Init_cf(MPI_Fint f_global_comm, const int n_code, const char **code_names, const int *l_code_names, const int *is_active_rank, const double *time_init, MPI_Fint *intra_comms);

void CWP_Cpl_create_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id,
                       const char *coupled_code_name, int l_coupled_code_name, CWP_Interface_t entities_dim,
                       CWP_Comm_t comm_type, CWP_Spatial_interp_t spatial_interp, int n_part, CWP_Dynamic_mesh_t displacement, CWP_Time_exch_t freq);

void CWP_Cpl_del_cf(const char *f_local_code_name, int l_local_code_name, const char *f_cpl_id, int l_cpl_id);

int CWP_N_uncomputed_tgts_get_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id, CWP_Dof_location_t target_location, int i_part);

const int *CWP_Uncomputed_tgts_get_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id);

int CWP_N_computed_tgts_get_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id);

const int *CWP_Computed_tgts_get_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id);

const double *CWP_Computed_tgts_dist_to_spatial_interp_get_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id);

void CWP_Spatial_interp_weights_compute_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id);

void CWP_Visu_set_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id, int freq, CWP_Visu_format_t format, const char *format_option, int l_format_option);

void CWP_Mesh_interf_finalize_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id);

void CWP_Mesh_interf_vtx_set_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id, int i_part, int n_pts, double coord[], CWP_g_num_t global_num[]);

int CWP_Mesh_interf_block_add_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id, CWP_Block_t block_type);

void CWP_Mesh_interf_block_std_set_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id, int i_part, int block_id, int n_elts, int connec[], CWP_g_num_t global_num[]);

void CWP_Mesh_interf_h_order_block_set_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id, int i_part, int block_id, int n_elts, int order, int connec[], CWP_g_num_t global_num[]);

void CWP_Mesh_interf_f_poly_block_set_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id, int i_part, int block_id, int n_elts, int connec_idx[], int connec[], CWP_g_num_t global_num[]);

void CWP_Mesh_interf_c_poly_block_set_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id, int i_part, int block_id, int n_elts, int n_faces,
                                      int connec_faces_idx[], int connec_faces[], int connec_cells_idx[], int connec_cells[], CWP_g_num_t global_num[]);

void CWP_Mesh_interf_from_cellface_set_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id, int i_part, int n_cells,
                                          int cell_face_idx[], int cell_face[], int n_faces, int face_vtx_idx[], int face_vtx[], CWP_g_num_t parent_num[]);

void CWP_Mesh_interf_from_faceedge_set_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id, int i_part, int n_faces,
                                          int face_edge_idx[], int face_edge[], int n_edges, int edge_vtx_idx[], int edge_vtx[], CWP_g_num_t parent_num[]);

void CWP_Field_create_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id, const char *field_id, int l_field_id, CWP_Type_t data_type, CWP_Field_storage_t storage, int n_component, CWP_Dof_location_t target_location, CWP_Field_exch_t exch_type, CWP_Status_t visu_status);

void CWP_Field_data_set_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id, const char *field_id, int l_field_id, int i_part, double data[]);

void CWP_Field_exch_cf(const char *local_code_name, int l_local_code_name, const char *cpl_id, int l_cpl_id);

void CWP_Field_issend_cf(const char *f_local_code_name, int l_local_code_name, const char *f_cpl_id, int l_cpl_id, const char *f_src_field_id, int l_src_field_id);

void CWP_Field_irecv_cf(const char *f_local_code_name, int l_local_code_name, const char *f_cpl_id, int l_cpl_id, const char *f_tgt_field_id, int l_tgt_field_id);

void CWP_Field_wait_issend_cf(const char *f_local_code_name, int l_local_code_name, const char *f_cpl_id, int l_cpl_id, const char *f_src_field_id, int l_src_field_id);

void CWP_Field_wait_irecv_cf(const char *f_local_code_name, int l_local_code_name, const char *f_cpl_id, int l_cpl_id, const char *f_tgt_field_id, int l_tgt_field_id);

#endif //CWP_CF_H_
