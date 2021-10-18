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

/*=============================================================================
 * Macro definitions
 *============================================================================*/


#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

/*============================================================================*
 *                                                                            *
 *                      Public function prototypes                            *
 *                      --------------------------                            *
 *                                                                            *
 *============================================================================*/

/*----------------------------------------------------------------------------*
 * General functions                                                          *
 *----------------------------------------------------------------------------*/


/**
 * \brief Initialize CWIPI.
 *
 * This function creates the MPI intra communicators of the codes from
 * the \p global_comm MPI communicator that contains all code ranks. This
 * function has to be called from all ranks contained in the \p global_comm.
 *
 * \param [in]  global_comm    MPI global communicator
 * \param [in]  n_code         Number of codes on the current rank
 * \param [in]  f_code_names   Names of codes on the current rank (size = \p n_code)
 * \param [in]  l_code_names   Length of code names on the current rank (size = \p n_code)
 * \param [in]  is_active_rank Is current rank have to be used by CWIPI (size = \p n_code)
 * \param [in]  time_init      Initial time (size = \p n_code)
 * \param [out] intra_comms    MPI intra communicators of each code (size = \p n_code)
 *
 */

void 
CWP_Init_cf (
  MPI_Fint f_global_comm, 
  const int n_code, 
  const char **f_code_names, 
  const int *l_code_names, 
  const int *is_active_rank, 
  const double *time_init, 
  MPI_Fint *f_intra_comms
); 

/*----------------------------------------------------------------------------*
 * Functions about current code properties                                    *
 *----------------------------------------------------------------------------*/

/**
 * \brief Update code state.
 *
 * \param [in] local_code_name   Fortran local code name
 * \param [in] l_local_code_name Length of Fortran local code name
 * \param [in] state             State
 *
 */

void
CWP_State_update_cf
(
 const char* local_code_name,
 const int l_local_code_name,
 const CWP_State_t state
);


/**
 * \brief Update code time.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] l_local_code_name Length of Fortran local code name
 * \param [in]  current_time Current time
 *
 */

void
CWP_Time_update_cf
(
 const char* local_code_name,
 const int l_local_code_name,
 const double current_time
);



/*----------------------------------------------------------------------------*
 * Functions about other code properties                               *
 *----------------------------------------------------------------------------*/


/**
 * \brief Return code state.
 *
 * \param [in]  code_name    Code name
 *
 * \return      Code state
 */

CWP_State_t
CWP_State_get_cf
(
 const char    *code_name,
 const int l_local_code_name
);


/**
 * \brief Create a coupling object and define its properties.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_coupled_code_name Distant or local coupled code name (Fortran)
 * \param [in]  l_coupled_code_name Length of Distant or local coupled code name 
 * \param [in]  comm_type           Communication type
 * \param [in]  spatial_interp      Spatial interpolation method
 * \param [in]  n_part              Number of interface partition
 * \param [in]  displacement        Mesh moving status
 * \param [in]  recv_freq_type      Type of receiving frequency
 *
 */

void 
CWP_Cpl_create_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id,
  const int l_cpl_id, 
  const char *f_coupled_code_name, 
  const int l_coupled_code_name, 
  const CWP_Interface_t entities_dim,
  const CWP_Comm_t comm_type,
  const CWP_Spatial_interp_t spatial_interp, 
  const int n_part, 
  const CWP_Dynamic_mesh_t displacement, 
  const CWP_Time_exch_t freq
); 


/**
 *
 * \brief Delete a coupling object.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 *
 */

void 
CWP_Cpl_del_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id, 
  const int l_cpl_id
);



/**
 *
 * \brief Return the number of uncomputed targets.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_field_id          Fortran Field identifier
 * \param [in]  l_field_id          Length of Fortran Field identifier
 * \param [in] i_part               Current partition
 *
 * \return                Number of uncomputed targets
 */

int 
CWP_N_uncomputed_tgts_get_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id, 
  const int l_cpl_id, 
  const char *f_field_id, 
  const int l_field_id, 
  int i_part
); 


/**
 *
 * \brief Return uncomputed targets.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_field_id          Fortran Field identifier
 * \param [in]  l_field_id          Length of Fortran Field identifier
 * \param [in]  i_part              Current partition
 *
 * \return                Uncomputed targets
 */

const int 
*CWP_Uncomputed_tgts_get_cf (
  const char *f_local_code_name,
  const int l_local_code_name,
  const char *f_cpl_id, 
  const int l_cpl_id, 
  const char *f_field_id, 
  const int l_field_id, 
  int i_part
);

/**
 *
 * \brief Return the number of computed targets. <b>(Not implemented yet)</b>
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_field_id          Fortran Field identifier
 * \param [in]  l_field_id          Length of Fortran Field identifier
 * \param [in]   i_part             Current partition
 *
 * \return                Number of computed targets
 */

int 
CWP_N_computed_tgts_get_cf (
  const char *f_local_code_name, 
  const int l_local_code_name,
  const char *f_cpl_id, 
  const int l_cpl_id, 
  const char *f_field_id, 
  const int l_field_id, 
  int i_part
);

/**
 *
 * \brief Return computed targets. <b>(Not implemented yet)</b>
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_field_id          Fortran Field identifier
 * \param [in]  l_field_id          Length of Fortran Field identifier
 * \param [in]  i_part              Current partition
 *
 * \return                Computed targets
 */

const int *
CWP_Computed_tgts_get_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id, 
  const int l_cpl_id, 
  const char *f_field_id, 
  const int l_field_id, 
  int i_part
);

/**
 * \brief Return distance from each target to the source interface. <b>(Not implemented yet)</b>
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 *
 * \return               Distance
 *
 */

const double *
CWP_Computed_tgts_dist_to_spatial_interp_get_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id,
  const int l_cpl_id
);


/**
 * \brief Compute spatial interpolation weights.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 *
 */

void 
CWP_Spatial_interp_weights_compute_cf (
  const char *f_local_code_name, 
  const int l_local_code_name,
  const char *f_cpl_id, 
  const int l_cpl_id
);

/*----------------------------------------------------------------------------*
 * Functions about visualization                                              *
 *----------------------------------------------------------------------------*/

/**
 * \brief Enable visualization output.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  freq                Output frequency
 * \param [in]  format              Output format to visualize exchanged fieldsDouble
 *                                  on the coupled mesh. Choice between :
 *                                  - "EnSight Gold"
 * \param [in]  f_format_option     Fortran Output options "opt1, opt2, ..."
 *                                  - text : output text files
 *                                  - binary : output binary files (default)
 * \param [in]  l_format_option    Length of Fortran option 
 */

void CWP_Visu_set_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id, 
  const int l_cpl_id, 
  int freq, 
  CWP_Visu_format_t format, 
  const char *f_format_option, 
  const int l_format_option
);


/**
 * \brief Finalize interface mesh.
 *
 * This function computes the global numbers of mesh entities if they are
 * not provided.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 *
 */

void 
CWP_Mesh_interf_finalize_cf
(
  const char *f_local_code_name, 
  int l_local_code_name, 
  const char *f_cpl_id, 
  int l_cpl_id
);


/**
 * \brief Set vertices.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  i_part              Current partition
 * \param [in]  n_pts               Number of points
 * \param [in]  coord               Coordinates (size = 3 * \p n_pts)
 * \param [in]  global_num          Pointer to parent element number (or NULL)
 *
 */

void 
CWP_Mesh_interf_vtx_set_cf(
  const char *f_local_code_name, 
  int l_local_code_name, 
  const char *f_cpl_id, 
  int l_cpl_id, int i_part, 
  int n_pts, double coord[], 
  CWP_g_num_t global_num[]
);



/**
 * \brief Add a connectivity block to the interface mesh.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  block_type          Block type
 *
 * \return block identifier
 */

int 
CWP_Mesh_interf_block_add_cf (
  const char *f_local_code_name, 
  int l_local_code_name,
  const char *f_cpl_id, 
  int l_cpl_id, 
  CWP_Block_t block_type
);




/**
 * \brief Set a standard block to the interface mesh.
 *
 * This function adds a connectivity block to the interface mesh.
 * Definition of element connectivity is :
 *
 *  - edge (\ref CWP_BLOCK_EDGE2) :
 *
 *   \code
 *       1 x-------x 2
 *   \endcode
 *
 *  - triangle (\ref CWP_BLOCK_FACE_TRIA3):
 *
 *   \code
 *       1 x-------x 3
 *          \     /
 *           \   /
 *            \ /
 *             x 2
 *   \endcode
 *
 *  - quadrangle (\ref CWP_BLOCK_FACE_QUAD4) :
 *
 *   \code
 *          4 x-------x 3
 *           /       /
 *          /       /
 *       1 x-------x2
 *   \endcode
 *
 *   - tetrahedron (\ref CWP_BLOCK_CELL_TETRA4) :
 *
 *   \code
 *             x 4
 *            /|\
 *           / | \
 *          /  |  \
 *       1 x- -|- -x 3
 *          \  |  /
 *           \ | /
 *            \|/
 *             x 2
 *   \endcode
 *
 *   - pyramid (\ref CWP_BLOCK_CELL_PYRAM5) :
 *
 *   \code
 *              5 x
 *               /|\
 *              //| \
 *             // |  \
 *          4 x/--|---x 3
 *           //   |  /
 *          //    | /
 *       1 x-------x 2
 *   \endcode
 *
 *  - prism (\ref CWP_BLOCK_CELL_PRISM6) :
 *
 *   \code
 *       4 x-------x 6
 *         |\     /|
 *         | \   / |
 *       1 x- \-/ -x 3
 *          \ 5x  /
 *           \ | /
 *            \|/
 *             x 2
 *   \endcode
 *
 *  -  hexaedron (\ref CWP_BLOCK_CELL_HEXA8) :
 *
 *   \code
 *          8 x-------x 7
 *           /|      /|
 *          / |     / |
 *       5 x-------x6 |
 *         | 4x----|--x 3
 *         | /     | /
 *         |/      |/
 *       1 x-------x 2
 *   \endcode
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  i_part              Partition identifier
 * \param [in]  block_id            Block identifier
 * \param [in]  n_elts              Number of elements
 * \param [in]  connec              Connectivity (size = n_vertex_elt * n_elts)
 * \param [in]  global_num          Pointer to global element number (or NULL)
 *
 */

void 
CWP_Mesh_interf_block_std_set_cf (
  const char *f_local_code_name, 
  int l_local_code_name, 
  const char *f_cpl_id, 
  int l_cpl_id, 
  int i_part, 
  int block_id,
  int n_elts, 
  int connec[], 
  CWP_g_num_t global_num[]
); 


/**
 * \brief Set the connectivity of a polygon block in a interface mesh partition.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  i_part              Current partition
 * \param [in]  block_id            Block identifier
 * \param [in]  n_elts              Number of elements
 * \param [in]  connec_idx          Connectivity index (\p connec_id[0] = 0 and
 *                                  size = \p n_elts + 1)
 * \param [in]  connec              Connectivity (size = \p connec_idx[\p n_elts])
 * \param [in]  global_num          Pointer to global element number (or NULL)
 *
 */

void 
CWP_Mesh_interf_f_poly_block_set_cf (
  const char *f_local_code_name, 
  int l_local_code_name, 
  const char *f_cpl_id, 
  int l_cpl_id, 
  int i_part, 
  int block_id, 
  int n_elts, 
  int connec_idx[], 
  int connec[], 
  CWP_g_num_t global_num[]
);


/**
 * \brief Adding a polyhedron connectivity block to the interface mesh.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  i_part              Current partition
 * \param [in]  block_id            Block identifier
 * \param [in]  n_elts              Number of elements
 * \param [in]  connec_cells_idx    Polyhedron to face index
 *                                  (\p src_poly_cell_face_idx[0] = 0 and
 *                                   size = \p n_elts + 1)
 * \param [in]  connec_cells        Polyhedron to face connectivity
 *                                  (size = \p cell_face_idx[\p n_elts])
 * \param [in]  n_faces             Number of faces
 * \param [in]  connec_faces_idx    Polyhedron face to vertex index
 *                                  (\p face_vertex_idx[0] = 0 and
 *                                   size = max(\p cell_face_connec) + 1)
 * \param [in]  connec_faces        Polyhedron face to vertex connectivity
 *                                  (size = \p face_vertex_idx[\p n_elts])
 * \param [in]  global_num          Pointer to global element number (or NULL)
 *
 */

void 
CWP_Mesh_interf_c_poly_block_set_cf (
  const char *f_local_code_name,
  int l_local_code_name,
  const char *f_cpl_id, 
  int l_cpl_id, 
  int i_part, 
  int block_id, 
  int n_elts, 
  int n_faces,
  int connec_faces_idx[], 
  int connec_faces[], 
  int connec_cells_idx[], 
  int connec_cells[], 
  CWP_g_num_t global_num[]
);



/**
 * \brief Define the interface mesh from a cell to face connectivity.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  i_part              Current partition
 * \param [in]  n_cells             Number of cells
 * \param [in]  cell_face_idx       Polyhedron to face index
 *                                  (\p src_poly_cell_face_idx[0] = 0 and
 *                                   size = \p n_elts + 1)
 * \param [in]  cell_face           Cell to face connectivity
 *                                  (size = \p cell_face_idx[\p n_elts])
 * \param [in]  n_faces             Number of faces
 * \param [in]  face_vtx_idx        Polyhedron face to vertex index
 *                                  (\p face_vtx_idx[0] = 0 and
 *                                   size = \p n_faces + 1)
 * \param [in]  face_vtx            Face to vertex connectivity
 *                                  (size = \p face_vtx_idx[\p n_elts])
 * \param [in]  parent_num          Pointer to parent element number (or NULL)
 *
 */

void 
CWP_Mesh_interf_from_cellface_set_cf (
  const char *f_local_code_name, 
  int l_local_code_name, 
  const char *f_cpl_id, 
  int l_cpl_id, 
  int i_part, 
  int n_cells,
  int cell_face_idx[], 
  int cell_face[], 
  int n_faces, 
  int face_vtx_idx[], 
  int face_vtx[], 
  CWP_g_num_t parent_num[]
);


/**
 * \brief Define the surface interface mesh from a face to edge connectivity.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  i_part              Current partition
 * \param [in]  n_faces             Number of cells
 * \param [in]  face_edge_idx       Polygon to edge index
 *                                  (\p face_edge_idx[0] = 0 and
 *                                   size =  \p n_faces + 1)
 * \param [in]  face_edge           Face to edge connectivity
 *                                  (size = \p face_edge_idx[\p n_faces])
 * \param [in]  n_edges             Number of faces
 * \param [in]  edge_vtx_idx        Polyhedron face to vertex index
 *                                  (\p edge_vtx_idx[0] = 0 and
 *                                   size = \p n_edges + 1)
 * \param [in]  edge_vtx            Face to vertex connectivity
 *                                  (size = \p edge_vtx_idx[\p n_edges])
 * \param [in]  parent_num          Pointer to parent element number (or NULL)
 *
 */

void 
CWP_Mesh_interf_from_faceedge_set_cf (
  const char *f_local_code_name, 
  int l_local_code_name, 
  const char *f_cpl_id, 
  int l_cpl_id, 
  int i_part, 
  int n_faces,
  int face_edge_idx[], 
  int face_edge[], 
  int n_edges, 
  int edge_vtx_idx[], 
  int edge_vtx[], 
  CWP_g_num_t parent_num[]
);

/*----------------------------------------------------------------------------*
 * Functions about field                                                      *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Create a new field.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_field_id          Fortran Field identifier
 * \param [in]  l_field_id          Length of Fortran Field identifier
 * \param [in]  data_type           Data type
 * \param [in]  storage             Storage type
 * \param [in]  n_component         Number of component
 * \param [in]  target_location     Target location
 * \param [in]  exch_type           Exchange type
 * \param [in]  visu_status         Visualization status
 *
 */

void 
CWP_Field_create_cf
(
  const char *f_local_code_name, 
  int l_local_code_name, 
  const char *f_cpl_id, 
  int l_cpl_id, 
  const char *f_field_id, 
  int l_field_id, 
  CWP_Type_t data_type, 
  CWP_Field_storage_t storage, 
  int n_component, 
  CWP_Dof_location_t target_location, 
  CWP_Field_exch_t exch_type, 
  CWP_Status_t visu_status
);



/**
 *
 * \brief Set field data.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  i_part              Current partition
 * \param [in]  data                Storage array (Mapping)
 *
 */

void 
CWP_Field_data_set_cf (
  const char *f_local_code_name, 
  int l_local_code_name, 
  const char *f_cpl_id, 
  int l_cpl_id, 
  const char *f_field_id, 
  int l_field_id, 
  int i_part,
  int map_type,
  double data[]
);

/**
 * \brief Send a spatially interpolated field to the coupled code with
 *        nonblocking communications.
 *
 * This function is independant of \ref CWP_Time_exch_t mode. The user has to
 * manually check the consistency of the exchanges.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_src_field_id      Fortran Source field id
 * \param [in]  l_src_field_id      Length of Source field id
 *
 *    
 */

void 
CWP_Field_issend_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id, 
  const int l_cpl_id, 
  const char *f_src_field_id, 
  const int l_src_field_id
);

/**
 *
 * \brief Receive a spatially interpolated field from the coupled code
 *        with nonblocking communications.
 *
 * This function is independant of \ref CWP_Time_exch_t mode. The user has to
 * manually check the consistency of the exchanges.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_tgt_field_id      Fortran Target field id
 * \param [in]  f_tgt_field_id      Length ofFortran Target field id
 *
 *
 */

void 
CWP_Field_irecv_cf(
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id, 
  const int l_cpl_id, 
  const char *f_tgt_field_id, 
  const int l_tgt_field_id
);

/**
 *
 * \brief Wait the end of an exchange related to request from \ref CWP_Field_issend.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_src_field_id      Fortran Source field id
 * \param [in]  l_src_field_id      Length of Source field id
 *
 */

void 
CWP_Field_wait_issend_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id,
  const int l_cpl_id, 
  const char *f_src_field_id, 
  const int l_src_field_id
);

/**
 *
 * \brief Wait the end of an exchange related to request from \ref CWP_Field_irecv.
 *
 * This function waits the end of exchange related to request
 * from \ref CWP_Field_irecv
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_tgt_field_id      Fortran Target field id
 * \param [in]  f_tgt_field_id      Length ofFortran Target field id
 *
 */

void 
CWP_Field_wait_irecv_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id, 
  const int l_cpl_id, 
  const char *f_tgt_field_id, 
  const int l_tgt_field_id
);


#endif //CWP_CF_H_
