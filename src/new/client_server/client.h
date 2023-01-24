#ifndef __CLIENT_H__
#define __CLIENT_H__
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
 *  Header for the current file
 *----------------------------------------------------------------------------*/

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

typedef struct t_client
{
  MPI_Comm comm;
  int      i_rank;
  MPI_Comm intra_comm;
  int      intra_i_rank;
  char     *code_name;
  int      server_port;
  int      flags;
  int      socket;
  int      max_msg_size;
  int      listen_socket;
  int      connected_socket;
  int      client_endianess;
  int      server_endianess;
  char     server_name[256];

}t_client,*p_client;

/*=============================================================================
 * Client CWIPI function interfaces
 *============================================================================*/

/**
 * \brief Initialize CWIPI.
 *        /!\ n_code == 1 in client-server mode
 *
 * \param [in]  n_code         Number of codes on the current rank
 * \param [in]  code_names     Names of codes on the current rank (size = \p n_code)
 * \param [in]  is_active_rank Is current rank have to be used by CWIPI (size = \p n_code)
 * \param [in]  time_init      Initial time (size = \p n_code)
 * \param [out] intra_comms    MPI intra communicators of each code (size = \p n_code)
 *
 */

void
CWP_client_Init
(
  MPI_Comm                  comm,
  char                    *config,
  const int                n_code,
  const char             **code_names,
  const CWP_Status_t      *is_active_rank,
  const double            *time_init
);

/**
 *
 * \brief Finalize CWIPI.
 *
 */

void
CWP_client_Finalize
(
 void
);

/**
 *
 * \brief Param_lock CWIPI.
 *
 * \param [in]  code_name  Code to lock
 *
 */

void
CWP_client_Param_lock
(
const char *code_name
);

/**
 *
 * \brief Param_unlock CWIPI.
 *
 * \param [in]  code_name  Code to unlock
 *
 */


void
CWP_client_Param_unlock
(
const char *code_name
);

/**
 *
 * \brief Param_add CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 * \param [in] initial_value    Initial value
 *
 */

void
CWP_client_Param_add
(
 const char        *local_code_name,
 const char        *param_name,
 const CWP_Type_t  data_type,
 void              *initial_value
);

/**
 *
 * \brief Param_get CWIPI.
 *
 * \param [in]  code_name  Local or distant code name
 * \param [in]  param_name Parameter name
 * \param [in]  data_type  Parameter type
 * \param [out] value      Parameter value
 *
 */

void
CWP_client_Param_get
(
 const char       *code_name,
 const char       *param_name,
 const CWP_Type_t  data_type,
 void             *value
);

/**
 *
 * \brief Param_set CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 * \param [in] value            Value
 *
 */

void
CWP_client_Param_set
(
 const char             *local_code_name,
 const char             *param_name,
 const CWP_Type_t        data_type,
 void                   *value
);

/**
 *
 * \brief Param_del CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 *
 */

void
CWP_client_Param_del
(
 const char       *local_code_name,
 const char       *param_name,
 const CWP_Type_t  data_type
);

/**
 *
 * \brief Param_n_get CWIPI.
 *
 * \param [in] code_name       Local or distant code name
 * \param [in] data_type       Parameter type,
 *
 * return  Number of parameters
 *
 */

int
CWP_client_Param_n_get
(
 const char             *code_name,
 const CWP_Type_t        data_type
);

/**
 *
 * \brief Param_list_get CWIPI.
 *
 * \param [in]  code_name      Local or distant code name
 * \param [in]  data_type      Parameter type,
 * \param [out] nParam         Number of parameters
 * \param [out] paramNames     Parameter names
 *
 *
 */

void
CWP_client_Param_list_get
(
 const char             *code_name,
 const CWP_Type_t        data_type,
 int                    *nParam,
 char                 ***paramNames
);

/**
 *
 * \brief Param_is CWIPI.
 *
 * \param [in] code_name      Local or distant code name
 * \param [in] param_name     Parameter name
 * \param [in] data_type      Parameter type
 *
 * return  1 : true / 0 : false
 *
 */

int
CWP_client_Param_is
(
 const char             *code_name,
 const char             *param_name,
 const CWP_Type_t        data_type
);

/**
 *
 * \brief Param_reduce CWIPI.
 *
 * The parameter name has to be the same for all codes.
 *
 * \param [in]  op           Operation
 * \param [in]  param_name   Parameter name
 * \param [in]  data_type    Parameter type,
 * \param [out] res          Result
 * \param [in]  nCode        Number of codes
 * \param       code_names   Codes name
 *
 */

void
CWP_client_Param_reduce
(
 const CWP_Op_t    op,
 const char       *param_name,
 const CWP_Type_t  data_type,
 void             *res,
 const int         nCode,
 const char      **code_names
);

/**
 * \brief Cpl_create CWIPI.
 *
 * \param [in]  local_code_name     Local code name
 * \param [in]  cpl_id              Coupling identifier
 * \param [in]  coupled_code_name   Distant or local coupled code name
 * \param [in]  comm_type           Communication type
 * \param [in]  spatial_interp      Spatial interpolation method
 * \param [in]  n_part              Number of interface partition
 * \param [in]  displacement        Mesh moving status
 * \param [in]  recv_freq_type      Type of receiving frequency
 *
 */

void
CWP_client_Cpl_create
(
 const char                *local_code_name,
 const char                *cpl_id,
 const char                *coupled_code_name,
 CWP_Interface_t            entities_dim,
 const CWP_Comm_t           comm_type,
 const CWP_Spatial_interp_t spatial_interp,
 const int                  n_part,
 const CWP_Dynamic_mesh_t   displacement,
 const CWP_Time_exch_t      recv_freq_type
);

/**
 *
 * \brief Cpl_del CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 */

void
CWP_client_Cpl_del
(
 const char *local_code_name,
 const char *cpl_id
);

/**
 * \brief State_update CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] state            State
 *
 */

void
CWP_client_State_update
(
 const char* local_code_name,
 const CWP_State_t state
);

/**
 * \brief Time_update CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in]  current_time Current time
 *
 */

void
CWP_client_Time_update
(
 const char* local_code_name,
 const double current_time
);

/**
 * \brief Output_file_set CWIPI.
 *
 * /!\ one file per code
 *
 * \param [in] output_filename    Output file directory or name
 *
 */

void
CWP_client_Output_file_set
(
  FILE *output_file
);


/**
 * \brief User_structure_set CWIPI.
 *
 * This structure can be called into a callback
 *
 * \param [in] local_code_name  Local code name
 * \param [in] user_structure   User structure
 *
 */

void
CWP_client_User_structure_set
(
 const char* local_code_name,
       void* user_structure
);


/**
 * \brief User_structure_get CWIPI.
 *
 * This structure can be called into a callback
 *
 * \param [in] local_code_name  Local code name
 *
 * \return  User structure
 *
 */

void *
CWP_client_User_structure_get
(
 const char* local_code_name
);

/**
 * \brief State_get CWIPI.
 *
 * \param [in]  code_name    Code name
 *
 * \return      Code state
 */

CWP_State_t
CWP_client_State_get
(
 const char    *code_name
);


/**
 * \brief Codes_nb_get CWIPI.
 *
 * \return Number of codes
 *
 */

int
CWP_client_Codes_nb_get
(
 void
);

/**
 * \brief Codes_list_get CWIPI.
 *
 * \return list of codes.
 */

const char **
CWP_client_Codes_list_get
(
void
);


/**
 * \brief Loc_codes_nb_get CWIPI.
 *
 * \return number of local codes.
 */

int
CWP_client_Loc_codes_nb_get
(
 void
);


/**
 * \brief Loc_codes_list_get CWIPI.
 *
 * \return list of local codes.
 */

const char **
CWP_client_Loc_codes_list_get
(
 void
);

/**
 * \brief Properties_dump CWIPI.
 *
 */

void
CWP_client_Properties_dump
(
void
);

/**
 *
 * \brief N_uncomputed_tgts_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return                Number of uncomputed targets
 */

int
CWP_client_N_uncomputed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);


/**
 *
 * \brief Uncomputed_tgts_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return                Uncomputed targets
 */

const int *
CWP_client_Uncomputed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);

/**
 *
 * \brief N_computed_tgts_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return                Number of computed targets
 */

int
CWP_client_N_computed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);

/**
 *
 * \brief Computed_tgts_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return                Computed targets
 */

const int *
CWP_client_Computed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);


/**
 *
 * \brief N_involved_srcs_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return               Number of involved sources
 */

int
CWP_client_N_involved_srcs_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);


/**
 *
 * \brief Involved_srcs_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 * \param [in] i_part           Current partition
 *
 * \return              Involved sources
 */

const int *
CWP_client_Involved_srcs_get
(
 const char *local_code_name,
 const char *cpl_id,
 const char *field_id,
 const int   i_part
);

/**
 * \brief Spatial_interp_weights_compute CWIPI.
 *
 * \param [in]  local_code_name     Local code name
 * \param [in]  cpl_id              Coupling identifier
 *
 */

void
CWP_client_Spatial_interp_weights_compute
(
 const char     *local_code_name,
 const char     *cpl_id
);


/**
 * \brief Spatial_interp_property_set CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  property_name    Name of the property
 * \param [in]  property_type    Type of the property ("double" or "int")
 * \param [in]  property_value   Value of the property
 *
 */

void
CWP_client_Spatial_interp_property_set
(
 const char     *local_code_name,
 const char     *cpl_id,
 const char     *property_name,
 const char     *property_type,
 const char     *property_value
);

/**
 * \brief Visu_set CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  freq             Output frequency
 * \param [in]  format           Output format to visualize exchanged fieldsDouble
 *                               on the coupled mesh. Choice between :
 *                               - "EnSight Gold"
 * \param [in]  format_option   Output options "opt1, opt2, ..."
 *                               - text : output text files
 *                               - binary : output binary files (default)
 */

void
CWP_client_Visu_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const int                   freq,
 const CWP_Visu_format_t     format,
 const char                 *format_option
);

/**
 * \brief SUser_tgt_pts_set CWIPI.
 *
 * This function must be called if the degrees of freedom locations are
 * \ref CWP_DOF_LOCATION_USER
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  n_pts            Number of points
 * \param [in]  coord            Coordinates (size = 3 * n_pts)
 * \param [in]  g_num            global number or NUL (size = n_pts)
 *
 */

void
CWP_client_User_tgt_pts_set
(
 const char    *local_code_name,
 const char    *cpl_id,
 const int      i_part,
 const int      n_pts,
 double         coord[],
 CWP_g_num_t    global_num[]
);

/**
 * \brief Mesh_interf_finalize CWIPI.
 *
 * This function computes the global numbers of mesh entities if they are
 * not provided.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 *
 */

void
CWP_client_Mesh_interf_finalize
(
 const char           *local_code_name,
 const char           *cpl_id
);


/**
 * \brief Mesh_interf_vtx_set CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  n_pts            Number of points
 * \param [in]  coord            Coordinates (size = 3 * \p n_pts)
 * \param [in]  global_num       Pointer to parent element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_vtx_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             n_pts,
 double                coord[],
 CWP_g_num_t           global_num[]
);


/**
 * \brief Mesh_interf_block_add CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  block_type       Block type
 *
 * \return block identifier
 */

int
CWP_client_Mesh_interf_block_add
(
 const char           *local_code_name,
 const char           *cpl_id,
 const CWP_Block_t     block_type
);


/**
 * \brief Mesh_interf_block_std_set CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  block_id         Block identifier
 * \param [in]  block_type       Block type
 * \param [in]  n_elts           Number of elements
 * \param [in]  connec           Connectivity (size = n_vertex_elt * n_elts)
 * \param [in]  global_num       Pointer to global element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_block_std_set
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const int          block_id,
 const int          n_elts,
 int                connec[],
 CWP_g_num_t        global_num[]
);


/**
 * \brief Mesh_interf_block_std_get CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  block_id         Block identifier
 * \param [out]  n_elts           Number of elements
 * \param [out]  connec           Connectivity (size = n_vertex_elt * n_elts)
 * \param [out]  global_num       Pointer to global element number (or NULL)
 */

void
CWP_client_Mesh_interf_block_std_get
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const int          block_id,
 int               *n_elts,
 int              **connec,
 CWP_g_num_t      **global_num
);

/**
  * \brief Get the standard block type
  *
  * \param [in]  local_code_name  Local code name
  * \param [in]  cpl_id           Coupling identifier
  * \param [in]  block_id    Block identifier
  *
  * \return block type
  */

CWP_Block_t
CWP_client_std_block_type_get
(
 const char             *local_code_name,
 const char             *cpl_id,
 const int               block_id
);

/**
 * \brief Mesh_interf_f_poly_block_set CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  block_id         Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  connec_idx       Connectivity index (\p connec_id[0] = 0 and
 *                               size = \p n_elts + 1)
 * \param [in]  connec           Connectivity (size = \p connec_idx[\p n_elts])
 * \param [in]  global_num       Pointer to global element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_f_poly_block_set
(
 const char             *local_code_name,
 const char             *cpl_id,
 const int               i_part,
 const int               block_id,
 const int               n_elts,
 int                     connec_idx[],
 int                     connec[],
 CWP_g_num_t             global_num[]
);


/**
 * \brief Mesh_interf_f_poly_block_get CWIPI.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  block_id         Block identifier
 * \param [out]  n_elts           Number of elements
 * \param [out]  connec_idx       Connectivity index (\p connec_id[0] = 0 and
 *                               size = \p n_elts + 1)
 * \param [out]  connec           Connectivity (size = \p connec_idx[\p n_elts])
 * \param [out]  global_num       Pointer to global element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_f_poly_block_get
(
 const char             *local_code_name,
 const char             *cpl_id,
 const int               i_part,
 const int               block_id,
 int                    *n_elts,
 int                   **connec_idx,
 int                   **connec,
 CWP_g_num_t           **global_num
);


/**
 * \brief Mesh_interf_c_poly_block_set CWIPI.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  block_id          Block identifier
 * \param [in]  n_elts            Number of elements
 * \param [in]  connec_cells_idx  Polyhedron to face index
 *                                (\p src_poly_cell_face_idx[0] = 0 and
 *                                 size = \p n_elts + 1)
 * \param [in]  connec_cells      Polyhedron to face connectivity
 *                                (size = \p cell_face_idx[\p n_elts])
 * \param [in]  n_faces           Number of faces
 * \param [in]  connec_faces_idx  Polyhedron face to vertex index
 *                                (\p face_vertex_idx[0] = 0 and
 *                                 size = max(\p cell_face_connec) + 1)
 * \param [in]  connec_faces      Polyhedron face to vertex connectivity
 *                                (size = \p face_vertex_idx[\p n_elts])
 * \param [in]  global_num        Pointer to global element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_c_poly_block_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             block_id,
 const int             n_elts,
 const int             n_faces,
 int                   connec_faces_idx[],
 int                   connec_faces[],
 int                   connec_cells_idx[],
 int                   connec_cells[],
 CWP_g_num_t           global_num[]
);


/**
 * \brief Mesh_interf_c_poly_block_get CWIPI.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  block_id          Block identifier
 * \param [out]  n_elts            Number of elements
 * \param [out]  connec_cells_idx  Polyhedron to face index
 *                                (\p src_poly_cell_face_idx[0] = 0 and
 *                                 size = \p n_elts + 1)
 * \param [out]  connec_cells      Polyhedron to face connectivity
 *                                (size = \p cell_face_idx[\p n_elts])
 * \param [out]  n_faces           Number of faces
 * \param [out]  connec_faces_idx  Polyhedron face to vertex index
 *                                (\p face_vertex_idx[0] = 0 and
 *                                 size = max(\p cell_face_connec) + 1)
 * \param [out]  connec_faces      Polyhedron face to vertex connectivity
 *                                (size = \p face_vertex_idx[\p n_elts])
 * \param [out]  global_num        Pointer to global element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_c_poly_block_get
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             block_id,
 int                  *n_elts,
 int                  *n_faces,
 int                 **connec_faces_idx,
 int                 **connec_faces,
 int                 **connec_cells_idx,
 int                 **connec_cells,
 CWP_g_num_t         **global_num
);


/**
 * \brief Mesh_interf_del CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 */

void
CWP_client_Mesh_interf_del
(
 const char *local_code_name,
 const char *cpl_id
);


/**
 * \brief Mesh_interf_from_cellface_set CWIPI.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  n_cells           Number of cells
 * \param [in]  cell_face_idx     Polyhedron to face index
 *                                (\p src_poly_cell_face_idx[0] = 0 and
 *                                 size = \p n_elts + 1)
 * \param [in]  cell_face         Cell to face connectivity
 *                                (size = \p cell_face_idx[\p n_elts])
 * \param [in]  n_faces           Number of faces
 * \param [in]  face_vtx_idx      Polyhedron face to vertex index
 *                                (\p face_vtx_idx[0] = 0 and
 *                                 size = \p n_faces + 1)
 * \param [in]  face_vtx          Face to vertex connectivity
 *                                (size = \p face_vtx_idx[\p n_elts])
 * \param [in]  global_num        Pointer to parent element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_from_cellface_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             n_cells,
 int                   cell_face_idx[],
 int                   cell_face[],
 const int             n_faces,
 int                   face_vtx_idx[],
 int                   face_vtx[],
 CWP_g_num_t           global_num[]
);


/**
 * \brief Mesh_interf_from_faceedge_set CWIPI.
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  n_faces           Number of cells
 * \param [in]  face_edge_idx     Polygon to edge index
 *                                (\p face_edge_idx[0] = 0 and
 *                                 size =  \p n_faces + 1)
 * \param [in]  face_edge         Face to edge connectivity
 *                                (size = \p face_edge_idx[\p n_faces])
 * \param [in]  n_edges           Number of faces
 * \param [in]  edge_vtx_idx      Polyhedron face to vertex index
 *                                (\p edge_vtx_idx[0] = 0 and
 *                                 size = \p n_edges + 1)
 * \param [in]  edge_vtx          Face to vertex connectivity
 *                                (size = \p edge_vtx_idx[\p n_edges])
 * \param [in]  global_num        Pointer to parent element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_from_faceedge_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             n_faces,
 int             face_edge_idx[],
 int             face_edge[],
 const int             n_edges,
 int             edge_vtx_idx[],
 int             edge_vtx[],
 CWP_g_num_t     global_num[]
);

/**
 *
 * \brief Field_create CWIPI.
 *
 * \param [in]  local_code_name Local code name
 * \param [in]  cpl_id          Coupling identifier
 * \param [in]  field_id        Field id
 * \param [in]  data_type       Data type
 * \param [in]  storage         Storage type
 * \param [in]  n_component     Number of component
 * \param [in]  target_location Target location
 * \param [in]  exch_type       Exchange type
 * \param [in]  visu_status     Visualization status
 *
 */

void
CWP_client_Field_create
(
 const char                  *local_code_name,
 const char                  *cpl_id,
 const char                  *field_id,
 const CWP_Type_t             data_type,
 const CWP_Field_storage_t    storage,
 const int                    n_component,
 const CWP_Dof_location_t     target_location,
 const CWP_Field_exch_t       exch_type,
 const CWP_Status_t           visu_status
);


/**
 *
 * \brief Field_data_set CWIPI.
 *
 * \param [in] local_code_name   Local code name
 * \param [in] cpl_id            Coupling identifier
 * \param [in] field_id          Field identifier
 * \param [in] i_part            Current partition
 * \param [in] data_type         Choice if data is setted for the source or the target
 * \param [in] data              Storage array (Mapping)
 * \param [in] n_component       Number of components of the field
 * \param [in] n_dof             Number of (cells, vtx or user)
 *
 */

void
CWP_client_Field_data_set
(
 const char              *local_code_name,
 const char              *cpl_id,
 const char              *field_id,
 const int                i_part,
 const CWP_Field_map_t    map_type,
 int                      n_entities,
 double                   data[]
);

/**
 *
 * \brief Field_n_component_get CWIPI.
 *  * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      number of field components
 *
 */

int
CWP_client_Field_n_component_get
(
 const char      *local_code_name,
 const char      *cpl_id,
 const char      *field_id
);

/**
 *
 * \brief Field_target_dof_location_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      Location of degrees of freedom
 *
 */

CWP_Dof_location_t
CWP_client_Field_target_dof_location_get
(
 const char      *local_code_name,
 const char      *cpl_id,
 const char      *field_id
);

/**
 *
 * \brief Field_storage_get CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      Field storage type
 */

CWP_Field_storage_t
CWP_client_Field_storage_get
(
 const char      *local_code_name,
 const char      *cpl_id         ,
 const char      *field_id
);

/**
 * \brief Field_del CWIPI.
 *
 * \param [in] local_code_name Local code name
 * \param [in]  cpl_id         Coupling identifier
 * \param [in]  field_id       Field identifier
 *
 */

void
CWP_client_Field_del
(
 const char      *local_code_name,
 const char      *cpl_id         ,
 const char      *field_id
);

/**
 * \brief Field_issend CWIPI.
 *
 * This function is independant of \ref CWP_Time_exch_t mode. The user has to
 * manually check the consistency of the exchanges.
 *
 * \param [in] local_code_name  Local code name
 * \param [in]  cpl_id          Coupling identifier
 * \param [in]  src_field_id    Source field id
 *
 *
 */

void
CWP_client_Field_issend
(
 const char     *local_code_name,
 const char     *cpl_id,
 const char     *src_field_id
);

/**
 *
 * \brief Field_irecv CWIPI.
 *
 * This function is independant of \ref CWP_Time_exch_t mode. The user has to
 * manually check the consistency of the exchanges.
 *
 * \param [in] local_code_name  Local code name
 * \param [in]  cpl_id          Coupling identifier
 * \param [in]  tgt_field_id    Target field id
 *
 *
 */

void
CWP_client_Field_irecv
(
 const char        *local_code_name,
 const char        *cpl_id,
 const char        *tgt_field_id
);

/**
 *
 * \brief Field_wait_issend CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] src_field_id     Source field id
 *
 */

void
CWP_client_Field_wait_issend
(
 const char  *local_code_name,
 const char  *cpl_id,
 const char  *src_field_id
);


/**
 *
 * \brief Field_wait_irecv CWIPI.
 *
 * This function waits the end of exchange related to request
 * from \ref CWP_Field_irecv
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] tgt_field_id     Target field id
 *
 */

int
CWP_client_Field_wait_irecv
(
 const char              *local_code_name,
 const char              *cpl_id,
 const char              *tgt_field_id,
 double                 **data
);

/**
 *
 * \brief Interp_from_location_unset CWIPI.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] src_field_id     Source field id
 *
 */

void
CWP_client_Interp_from_location_unset
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *src_field_id
);


/**
 *
 * \brief Interp_from_location_set CWIPI.
 *
 * This function takes into account an user interpolation function written with
 * void (*\ref CWP_Interp_from_location_t) interface.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] src_field_id     Source field id
 * \param [in] data_type        Field data type
 * \param [in] fct              Function
 *
 */

void
CWP_client_Interp_from_location_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *src_field_id,
 CWP_Interp_from_location_t  fct
);

/*=============================================================================
 * Client CWIPI function interfaces that are not implemented yet
 *============================================================================*/

/**
 * \brief Mesh_interf_h_order_block_set CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  block_id         Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  order            Element order
 * \param [in]  connec           Connectivity (size = n_vertex_elt * n_elts)
 * \param [in]  global_num       Pointer to global element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_h_order_block_set
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const int          block_id,
 const int          n_elts,
 const int          order,
 int                connec[],
 CWP_g_num_t        global_num[]
);


/**
 * \brief Mesh_interf_h_order_block_get CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  block_id         Block identifier
 * \param [out]  n_elts           Number of elements
 * \param [out]  order            Element order
 * \param [out]  connec           Connectivity (size = n_vertex_elt * n_elts)
 * \param [out]  global_num       Pointer to global element number (or NULL)
 *
 */

void
CWP_client_Mesh_interf_h_order_block_get
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const int          block_id,
 int               *n_elts,
 int               *order,
 int              **connec,
 CWP_g_num_t      **global_num
);


/**
 * \brief Cpl_trans_init CWIPI. <b>(Not implemented yet)</b>
 *
 * This function defines the translation direction. The movement is updated
 * before the end of the current time step by \ref CWP_Cpl_trans_update
 * Two coupled codes have to define the same properties. The distant code is always
 * considered as a static interface.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  vect             Direction
 *
 */

void
CWP_client_Cpl_trans_init
(
 const char      *local_code_name,
 const char      *cpl_id,
 const double     vect[3]
);


/**
 * \brief Cpl_trans_update CWIPI. <b>(Not implemented yet)</b>
 *
 * This function computes the next time step position from a relative distance. If
 * it is a known position, spatial interpolation weights are not recomputed. Otherwise,
 * spatial interpolation weights are computed from previous results.
 *
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  dist             Relative distance from previous displacement
 *
 */

void
CWP_client_Cpl_trans_update
(
 const char      *local_code_name,
 const char      *cpl_id,
 const double     dist
);


/**
 * \brief Cpl_rotation_init CWIPI. <b>(Not implemented yet)</b>
 *
 * This function defines the rotation properties. The movement is updated
 * before the end of the current time step by \ref CWP_Cpl_rotation_update.
 * Two coupled codes have to define the same properties. The distant code is always
 * considered as a static interface.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  vect             Direction
 * \param [in]  center           Center
 *
 */

void
CWP_client_Cpl_rotation_init
(
 const char      *local_code_name,
 const char      *cpl_id,
 const double     vect[3],
 const double     center[3]
);


/**
 * \brief Cpl_rotation_update CWIPI. <b>(Not implemented yet)</b>
 *
 * This function computes the next time step position from a relative angle. If
 * it is a known position, the spatial interpolation weights are not reprocessed. Otherwise,
 * the spatial interpolation weights are computed from previous results.
 *
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  local_code_name  Local code name
 * \param [in]  angle            Relative angle from previous displacement
 *
 */

void
CWP_client_Cpl_rotation_update
(
 const char      *local_code_name,
 const char      *cpl_id,
 const double     angle
);


/**
 * \brief Cpl_storage_properties_set CWIPI. <b>(Not implemented yet)</b>
 *
 * This functions activates the storage of the spatial interpolation weights in case of rotation
 * or translation of the coupling interface.
 *
 * \param [in] cpl_id              Coupling identifier
 * \param [in] local_code_name     Local code name
 * \param [in] buffer_size         Size of buffer (Mo) on each coupling
 *                                 communicator rank (same value for each)
 * \param [in] disk_storage_size   Total size of disk storage when the buffer
 *                                 is full (Mo) (Same value for each rank)
 *
 */

void
CWP_client_Cpl_storage_properties_set
(
 const char     *local_code_name,
 const char     *cpl_id,
 const int       buffer_size,
 const int       disk_storage_size
);


/**
 *
 * \brief Interp_from_intersect_set CWIPI. <b>(Not implemented yet)</b>
 *
 * This function takes into account an user interpolation function written with
 * void (*\ref CWP_Interp_from_intersect_t) interface.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] fct              Function
 *
 */

void
CWP_client_Interp_from_intersect_set
(
 const char                *local_code_name,
 const char                *cpl_id,
 CWP_Interp_from_intersect_t fct
);

/**
 *
 * \brief Interp_from_closest_pts_set CLIENT. <b>(Not implemented yet)</b>
 *
 * This function takes into account an user interpolation function written with
 *  void (*\ref CWP_Interp_from_closest_pts_t) interface.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] fct              Function
 *
 */

void
CWP_client_Interp_from_closest_pts_set
(
 const char                     *local_code_name,
 const char                     *cpl_id,
 CWP_Interp_from_closest_pts_t   fct
);



/**
 * \brief Computed_tgts_dist_to_spatial_interp_get CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 *
 * \return               Distance
 *
 */

const double *
CWP_client_Computed_tgts_dist_to_spatial_interp_get
(
 const char *local_code_name,
 const char *cpl_id
);

/**
 * \brief Recv_freq_set CWIPI. <b>(Not implemented yet)</b>
 *
 * This function set the receiving frequency. It must be used when
 * the type of receiving frequency is \ref CWP_TIME_EXCH_N_TIME_STEP
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  n_step           Frequency in steps number
 *
 */

void
CWP_client_Recv_freq_set
(
 const char      *local_code_name,
 const char      *cpl_id,
 const int        n_step
);

/**
 * \brief next_recv_time_set CWIPI. <b>(Not implemented yet)</b>
 *
 * It must be used when the type of receiving frequency is
 * \ref CWP_TIME_EXCH_ASYNCHRONOUS
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  next_time        Next receiving time
 *
 */

void
CWP_client_next_recv_time_set
(
 const char      *local_code_name,
 const char      *cpl_id,
 const double     next_time
);


/**
 * \brief Cpl_time_step_set CWIPI. <b>(Not implemented yet)</b>
 *
 * This function sets the coupling time step. It must be used when
 * the type of receiving frequency is \ref CWP_TIME_EXCH_CPL_TIME_STEP
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  next_time_step   Coupling time step
 *
 */

void
CWP_client_Cpl_time_step_set
(
 const char      *local_code_name,
 const char      *cpl_id,
 const int        next_time_step
);


/**
 * \brief Field_exch CWIPI. <b>(Not implemented yet)</b>
 *
 * This function exchanges the interpolated fields for each coupling depending
 * on mode of time exchange \ref CWP_Time_exch_t.
 *
 * \param [in] local_code_name      Local code name
 * \param [in] cpl_id               Coupling identifier
 *
 */

void
CWP_client_Field_exch
(
 const char *local_code_name,
 const char *cpl_id
);


/**
 * \brief Mesh_interf_shared_pdm_nodal CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in] local_code_name   Local code name
 * \param [in] cpl_id            Coupling identifier
 * \param [in] i_part            Current partition
 * \param [in] pdm_nodal         pdm nodal mesh
 *
 */


void
CWP_client_Mesh_interf_shared_pdm_nodal
(
 const char   *local_code_name,
 const char   *cpl_id,
 const int     i_part,
 void         *pdm_nodal
);


/**
 *
 * \brief Field_data_type_get CWIPI.  <b>(Not implemented yet)</b>
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      Field data type
 *
 */

CWP_Type_t
CWP_client_Field_data_type_get
(
 const char      *local_code_name,
 const char      *cpl_id         ,
 const char      *field_id
);


/**
 *
 * \brief Field_gradient_data_set CWIPI. <b>(Not implemented yet)</b>
 *
 * \param [in] local_code_name Local code name
 * \param [in] cpl_id          Coupling identifier
 * \param [in] field_id        Field identifier
 * \param [in] i_part          Current partition
 * \param [in] order           Order
 * \param [in] data_type       Choice if data is setted for the source or the target
 * \param [in] data            Storage array (Mapping)
 *
 */

void
CWP_client_Field_gradient_data_set
(
 const char             *local_code_name,
 const char             *cpl_id,
 const char             *field_id,
 const int               i_part,
 const int               order,
 const CWP_Field_storage_t  storage_type,
 double                  data[]
);

/*=============================================================================
 * Client function interfaces
 *============================================================================*/

/* Connect to a server */

int
CWP_client_connect
(
 MPI_Comm  comm,
 MPI_Comm  intra_comm,
 const char* server_name,
 int server_port,
 int flags
);

/* Disconnect */

int
CWP_client_disconnect
(
 void
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CLIENT_H__ */
