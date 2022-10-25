#ifndef __STRUCT_H__
#define __STRUCT_H__
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

/*
  This file is inspired from OpenPALM.
  OpenPALM is a free software under the GNU Lesser General Public License.
  See: https://www.cerfacs.fr/globc/PALM_WEB/
*/

/*----------------------------------------------------------------------------
 * Standard C++ library headers
 *----------------------------------------------------------------------------*/

#include <string>
#include <map>
#include <tuple>

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "mpi.h"
#include "cwp.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Types definition
 *============================================================================*/

typedef struct t_field_settings
{
  int             i_part;
  CWP_Field_map_t map_type;
  int             n_component;
  int             n_entities;
} t_field_settings, *p_field_settings;

typedef struct t_field
{
  // Field_data_set
  double *data;

  // Field_wait_irecv
  std::map<std::tuple<std::string, std::string , std::string> , p_field_settings> field_settings;
} t_field, *p_field;

typedef struct t_coupling
{
  // Mesh_interf_vtx_set
  double *vtx_coord;
  CWP_g_num_t *vtx_global_num;

  // Mesh_interf_c_poly_block_set
  int *connec_faces_idx;
  int *connec_faces;
  int *connec_cells_idx;
  int *connec_cells;
  CWP_g_num_t *cell_global_num;

  // Mesh_interf_f_poly_block_set
  int *connec_idx;
  int *connec;
  CWP_g_num_t *elt_global_num;

  // Mesh_interf_from_faceedge_set
  int *face_edge_idx;
  int *face_edge;
  int *edge_vtx_idx;
  int *edge_vtx;
  CWP_g_num_t *face_global_num;

  // Mesh_interf_block_std_set
  int *std_connec;
  CWP_g_num_t *std_global_num;

  // User_tgt_pts_set
  double *usr_coord;
  CWP_g_num_t *usr_global_num;

  // Spatial_interp_property_set
  char *property_name;
  char *property_type;
  char *property_value;

  // field_id
  std::map<std::string, p_field> field;
} t_coupling, *p_coupling;

typedef struct t_cwp
{
  // Codes_list_get
  int          n_code_names;
  const char **code_names;

  // Loc_codes_list_get
  int          n_loc_code_names;
  const char **loc_code_names;

  // Param_list_get, Param_get
  int                             n_param_names;
  char                          **param_names;
  std::map<std::string, char *>   char_param_value;

  // Output_file_set
  FILE *output_file;

  // cpl_id
  std::map<std::string, p_coupling> coupling;
} t_cwp, *p_cwp;

typedef struct t_server_mpi
{
  MPI_Comm  global_comm;
  MPI_Comm *intra_comms;
} t_server_mpi, *p_server_mpi;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __STRUCT_H__ */
