#ifndef __CWP_H__
#define __CWP_H__
/*
  This file is part of the CWIPI library.

  Copyright (C) 2011-20  ONERA

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


/** \file cwp.h
  * \brief CWIPI new API header file
  *
  */


#include <mpi.h>
#include <stdio.h>


/*=============================================================================
 * Macro definitions
 *============================================================================*/

#if !defined (__hpux) && !defined (_AIX)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \typedef CWP_g_num_t
 * \brief Long int in cwipi
 *
 */

typedef long CWP_g_num_t;

/*============================================================================
 * Enumeration definitions
 *============================================================================*/

/**
 * \enum CWP_Type_t
 * \brief  Type of data.
 */

typedef enum {

  CWP_DOUBLE,      /*!< Double precision type */
  CWP_INT,         /*!< Integer type */
  CWP_CHAR         /*!< String type */

} CWP_Type_t;

/**
 * \enum CWP_Visu_format_t
 * \brief  List of available formats for Visualization.
 *
 */

typedef enum {

  CWP_VISU_FORMAT_ENSIGHT      /*!< Ensight visualization format */

} CWP_Visu_format_t;



/**
 * \enum CWP_Comm_t
 * \brief Communication mode
 *
 * CWP_Comm_t gives the different communication mode
 */

typedef enum {

  CWP_COMM_PAR_WITH_PART,    /*!< Parallel communication
                                  on partitioned source mesh */
  CWP_COMM_PAR_WITHOUT_PART, /*!< Parallel communication
                                  on unpartitioned source mesh defined on
                                  all processes */
  CWP_COMM_SEQ,              /*!< Parallel communication
                                  on unpartitioned source mesh defined on
                                  master processes */
  CWP_COMM_INTERNAL         /*!< Internal communication within a process */

} CWP_Comm_t;

/**
 * \enum CWP_Time_exch_t
 * \brief  Modes of time exchange used by \ref CWP_Exch function to determine when to launch interpolation/exchange.
 *
 */


typedef enum {

  CWP_TIME_EXCH_EACH_TIME_STEP,      /*!< Exchange at each time step */
  CWP_TIME_EXCH_N_TIME_STEP,         /*!< Exchange every \it n time steps  */
  CWP_TIME_EXCH_CPL_TIME_STEP,       /*!< Coupling time step        */
  CWP_TIME_EXCH_ASYNCHRONOUS,        /*!< Exchanges are asynchronous with temporal interpolation */
  CWP_TIME_EXCH_SLAVE,               /*!< Give a converged state    */
  CWP_TIME_EXCH_MASTER               /*!< Request a converged state */
} CWP_Time_exch_t;


/**
 * \enum CWP_Dof_location_t
 * \brief Modes of degrees of freedom location.
 *
 */

typedef enum {

  CWP_DOF_LOCATION_CELL_CENTER,   /*!< Field defined in cell point  */
  CWP_DOF_LOCATION_NODE,     /*!< Cell vertex field */
  CWP_DOF_LOCATION_USER           /*!< User defined field */
} CWP_Dof_location_t ;


/**
 * \enum CWP_Field_exch_t
 * \brief Modes of field exchange.
 */

typedef enum {

  CWP_FIELD_EXCH_SEND,        /*!< Send */
  CWP_FIELD_EXCH_RECV,        /*!< Receive */
  CWP_FIELD_EXCH_SENDRECV     /*!< Send and receive */

} CWP_Field_exch_t ;


/**
 * \enum CWP_Field_storage_t
 * \brief Modes of field storage.
 */

typedef enum {

  CWP_FIELD_STORAGE_INTERLACED,  /*!< Interlaced storage */
  CWP_FIELD_STORAGE_BLOCK        /*!< Block storage */

} CWP_Field_storage_t ;


/**
 * \enum CWP_Status_t
 * \brief on/off status
 *
 */

typedef enum {

  CWP_STATUS_OFF,           /*!< OFF */
  CWP_STATUS_ON             /*!< ON */

} CWP_Status_t;

/**
 * \enum CWP_Err_t
 * \brief Error codes.
 *
 */

typedef enum {

  CWP_ERR_NO_ERROR,       /*!< Output without error */
  CWP_ERR_DEFAULT         /*!< Output with default error */

} CWP_Err_t;

/**
 * \enum CWP_Block_t
 * \brief Elements taken into account
 *
 * CWP_Block_t defines elements taken into account
 */

typedef enum {

  CWP_BLOCK_NODE,          /*!< Node */
  CWP_BLOCK_EDGE2,         /*!< Edge with two nodes */
  CWP_BLOCK_FACE_TRIA3,    /*!< Triangle with three nodes */
  CWP_BLOCK_FACE_QUAD4,    /*!< Quadrangle with three nodes */
  CWP_BLOCK_FACE_POLY,     /*!< Generic polygon */
  CWP_BLOCK_CELL_TETRA4,   /*!< Tetrahedron with four nodes */
  CWP_BLOCK_CELL_HEXA8,    /*!< Hexahedron with eight nodes */
  CWP_BLOCK_CELL_PRISM6,   /*!< Prism with six nodes */
  CWP_BLOCK_CELL_PYRAM5,   /*!< Pyramid with five nodes */
  CWP_BLOCK_CELL_POLY      /*!< Generic polyhedron */

} CWP_Block_t;


/**
 * \enum CWP_Dynamic_mesh_t
 * \brief Modes of time dynamic mesh.
 *
 */

typedef enum {

  CWP_DYNAMIC_MESH_STATIC,      /*!< Static mesh*/
  CWP_DYNAMIC_MESH_DEFORMABLE,  /*!< Deformable mesh with constant topology */
  CWP_DYNAMIC_MESH_VARIABLE     /*!< Variable mesh topology */
} CWP_Dynamic_mesh_t;

/**
 * \enum CWP_Spatial_interp_t
 * \brief List of available spatial interpolation methods.
  */

typedef enum {

  CWP_SPATIAL_INTERP_FROM_CLOSEST_POINT_LEAST_SQUARES, /*!< Least squares from closest points */
  CWP_SPATIAL_INTERP_FROM_INTERSECTION,                /*!< Meshes intersection */
  CWP_SPATIAL_INTERP_FROM_LOCATION                     /*!< Location into a mesh */

} CWP_Spatial_interp_t;

/**
 * \enum CWP_Interface_t
 * \brief Coupling interfaces.
 */

typedef enum{

  CWP_INTERFACE_POINT,    /*!< Point interface */
  CWP_INTERFACE_LINEAR,   /*!< Linear interface */
  CWP_INTERFACE_SURFACE,  /*!< Surface interface */
  CWP_INTERFACE_VOLUME    /*!< Volume interface */

} CWP_Interface_t;

/**
 * \enum CWP_State_t
 * \brief Code state.
 *
 */

typedef enum {

  CWP_STATE_IN_PROGRESS,         /*!< Code running */
  CWP_STATE_END,                 /*!< Computation end */
  CWP_STATE_OUTPUT_ERROR         /*!< Output on error */

} CWP_State_t;

/**
 * \enum CWP_Op_t
 * \brief Operations on control parameters.
 *
 */

typedef enum {

  CWP_OP_MIN,            /*!< Minimum */
  CWP_OP_MAX,            /*!< Maximum */
  CWP_OP_SUM             /*!< Sum */

} CWP_Op_t;

/*============================================================================
 * User interpolation type
 *============================================================================*/

/**
 * \typedef void (*CWP_Interp_from_location_t)
 * \brief User interpolation function interface from location into a mesh.
 *
 * void (*CWP_Interp_from_location_t) defines the user interpolation
 * interface to take into account an user interpolation from location of target
 * points into the source mesh. Use \ref CWP_Interp_from_location_set to activate
 * the function.
 *
 * \param [in]  interface_type              Interface type
 * \param [in]  n_src_vtcs                  Number of source mesh vertices
 * \param [in]  n_src_std_elts              Number of source mesh standard elements
 * \param [in]  n_src_poly                  Number of source mesh polyhedra
 * \param [in]  n_tgt_pts                   Number of target points
 * \param [in]  src_vts_coords              Source Mesh vertices coordinates
 * \param [in]  src_parent_elts_num         Pointer to parent element number
 *                                          (or NULL)
 * \param [in]  src_parent_vtcs_num         Pointer to parent vertex number
 *                                          (or NULL)
 * \param [in]  src_connec_idx              Element to vertex index
 *                                          (src_connec_idx[0] = 0 and
 *                                          size = n_src_std_element + 1)
 * \param [in]  src_connec                  Element to vertex connectivity.
 *                                          (size = src_connec_idx[n_src_std_element])
 * \param [in]  src_poly_cell_face_idx      Polyhedron to face index
 *                                          (src_poly_cell_face_idx[0] = 0 and
 *                                          size = n_src_polyhedron + 1)
 * \param [in]  src_poly_cell_face_connec   Polyhedron to face connectivity
 *                                          (size = src_poly_cell_face_idx[n_src_polyhedron])
 * \param [in]  src_poly_face_vtx_idx       Polyhedron face to vertex index
 *                                          (src_poly_face_vertex_idx[0] = 0 and
 *                                          size_idx = max(src_poly_cell_face_connec) + 1)
 * \param [in]  src_poly_face_vtx_connec    Polyhedron face to vertex connectivity
 *                                          (size = src_poly_face_vertex_idx[size_iudx - 1])
 * \param [in]  tgt_pts_coords              Target points coordinates
 *                                          (size = 3 * n_tgt_pts)
 * \param [in]  tgt_pts_target_location            target points location
 *                                          (size = n_tgt_pts)
 * \param [in]  tgt_pts_dist                target points distance to location element
 *                                          (size = n_tgt_pts)
 * \param [in]  tgt_pts_bary_coords_idx     Index of Barycentric coordinates target points
 *                                          in location element
 *                                          (tgt_pts_bary_coords_idx[0] = 0 and
 *                                          size = n_tgt_pts + 1)
 * \param [in]  tgt_pts_bary_coords         Barycentric coordinates target points
 *                                          in location element
 *                                          (size = tgt_pts_bary_coords_idx[n_tgt_pts])
 * \param [in]  stride                      Number of field components
 * \param [in]  src_field_target_location          source field location
 * \param [in]  src_field                   source field
 *                                          (size depends on field type and stride)
 * \param [in]  src_field_target_location          target field location
 * \param [out] tgt_field                   target field
 *                                          (size = stride * n_tgt_pts)
 */


typedef void (*CWP_Interp_from_location_t)
(
 const int                   interface_type,
 const int                   n_src_vtcs,
 const int                   n_src_std_elts,
//const int                  n_src_poly,
 const int                   n_tgt_pts,
 const double                src_vtcs_coords[],
//const CWP_g_num_t          src_global_elts_num[],
// const CWP_g_num_t         src_global_vtcs_num[],
 const int                   src_connec_idx[],
 const int                   src_connec[],
 //const int                 src_poly_cell_face_idx[],
 //const int                 src_poly_cell_face_connec[],
 //const int                 src_poly_face_vtx_idx[],
 //onst int                  src_poly_face_vtx_connec[],
 const double                tgt_pts_coords[],
 const int                   tgt_pts_target_location[],
 const double                tgt_pts_dist[],
 const int                   tgt_pts_bary_coords_idx[],
 const double                tgt_pts_bary_coords[],
 const int                   stride,
 const CWP_Dof_location_t     src_field_location,
 const void                 *src_field,
 const CWP_Dof_location_t     tgt_field_location,
 void                       *tgt_field
);

/**
 * \typedef void (*CWP_Interp_from_intersect_t)
 * \brief User interpolation function interface from intersection between meshes. <b>(Not implemented yet)</b>
 *
 * void (*CWP_Interp_from_intersect_t) defines the user interpolation
 * interface to take into account an user interpolation from intersection
 * between source and target meshes. Use \ref CWP_Interp_from_intersect_set to activate
 * the function.
 *
 * \param [in]  interface_type              Interface type
 *
 */

typedef void (*CWP_Interp_from_intersect_t)
(
 const int interface_type
);

/**
 * \typedef void (*CWP_Interp_from_closest_pts_t)
 * \brief User interpolation function from closest points. <b>(Not implemented yet)</b>
 *
 * void (*CWP_Interp_from_closest_pts_t) defines the user interpolation
 * interface to take into account an user interpolation from <i>n</i> closest
 * points. Use \ref CWP_Interp_from_closest_pts_set to activate
 * the function.
 *
 * \param [in]  interface_type              Interface type
 *
 */

typedef void (*CWP_Interp_from_closest_pts_t)
(
 const int interface_type
);

/*=============================================================================
 * Static global variables
 *============================================================================*/

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
 * \param [in]  code_names     Names of codes on the current rank (size = \p n_code)
 * \param [in]  is_active_rank Is current rank have to be used by CWIPI (size = \p n_code)
 * \param [in]  time_init      Initial time (size = \p n_code)
 * \param [out] intra_comms    MPI intra communicators of each code (size = \p n_code)
 *
 */

void
CWP_Init
(
 const MPI_Comm           global_comm,
 const int                n_code,
 const char             **code_names,
 const CWP_Status_t      *is_active_rank,
 const double            *time_init,
 MPI_Comm                *intra_comms
);

/**
 *
 * \brief Finalize CWIPI.
 *
 */

void
CWP_Finalize
(void
);

/*----------------------------------------------------------------------------*
 * Functions about current code properties                                    *
 *----------------------------------------------------------------------------*/

/**
 * \brief Update code state.
 *
 * This function set the code state.
 *
 * \param [in] state    State
 *
 */

void
CWP_State_update
(
 const char* local_code_name,
 const CWP_State_t state
);



MPI_Comm
CWP_Connectable_comm_get
(
  char* local_code_name
);

/**
 * \brief Update code time
 *
 * This function update the code current time.
 *
 * \param [in]  current_time Current time
 *
 */

void
CWP_Time_update
(
 const char* local_code_name,
 const double current_time
);


/**
 * \brief Define output file.
 *
 * \param [in] output_file    Output file
 *
 */

void
CWP_Output_file_set
(
FILE *output_file
);


//**
// * \brief Writing output to fortran file.
// *
// * This function set the file fortran logical unit for writing output.
// *
// * \param [in]  iunit        File fortan logical unit
// *
// */

//void
//PROCF (cwp_output_fortran_unit_set, CWP_OUTPUT_FORTRAN_UNIT_SET)
//(
// int *iunit
//);

/*----------------------------------------------------------------------------*
 * Functions about other code properties                               *
 *----------------------------------------------------------------------------*/

/**
 * \brief Code state.
 *
 * This function return the state code.
 *
 * \param [in]  code_name    Code name
 *
 */

CWP_State_t
CWP_State_get
(
 const char    *code_name
);


/**
 * \brief Return the number of codes known by CWIPI.
 *
 */

int
CWP_Codes_nb_get
(
 void
);

/**
 * \brief Return the list of code names known by CWIPI.
 *
 */

const char **
CWP_Codes_list_get
(
void
);


/**
 * \brief Return the number of local codes known by CWIPI.
 *
 */

int
CWP_Loc_codes_nb_get
(
 void
);


/**
 * \brief Return the list of local code names known by CWIPI.
 *
 */

const char **
CWP_Loc_codes_list_get
(
 void
);

/*----------------------------------------------------------------------------*
 * Functions about properties                                                 *
 *----------------------------------------------------------------------------*/

/**
 * \brief Dump code properties.
 *
 * This function dump code properties.
 *
 */

void
CWP_Properties_dump
(void
);

/*----------------------------------------------------------------------------*
 * General functions about coupling                                           *
 *----------------------------------------------------------------------------*/

/**
 * \brief Create a coupling object and define its properties.
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
CWP_Cpl_create
(
 const char               *local_code_name,
 const char               *cpl_id,
 const char               *coupled_code_name,
 const CWP_Comm_t          comm_type,
 const CWP_Spatial_interp_t spatial_interp,
 const int                 n_part,
 const CWP_Dynamic_mesh_t  displacement,
 const CWP_Time_exch_t          recv_freq_type
);


/**
 * \brief Initialize a translation displacement. <b>(Not implemented yet)</b>
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
CWP_Cpl_trans_init
(
 const char      *local_code_name,
 const char      *cpl_id,
 const double     vect[3]
);


/**
 * \brief Define the next time step position in a local mesh translation. <b>(Not implemented yet)</b>
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
CWP_Cpl_trans_update
(
 const char      *local_code_name,
 const char      *cpl_id,
 const double     dist
);


/**
 * \brief Initialize a rotation displacement. <b>(Not implemented yet)</b>
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
CWP_Cpl_rotation_init
(
 const char      *local_code_name,
 const char      *cpl_id,
 const double     vect[3],
 const double     center[3]
);


/**
 * \brief Define the next time step position in a local mesh rotation. <b>(Not implemented yet)</b>
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
CWP_Cpl_rotation_update
(
 const char      *local_code_name,
 const char      *cpl_id,
 const double     angle
);


/**
 * \brief Set the storage properties of the spatial interpolation weights. <b>(Not implemented yet)</b>
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
CWP_Cpl_storage_properties_set
(
 const char     *local_code_name,
 const char     *cpl_id,
 const int       buffer_size,
 const int       disk_storage_size
);


/**
 *
 * \brief Delete a coupling object.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 */

void
CWP_Cpl_del
(
 const char *local_code_name,
 const char *cpl_id
);

/**
 *
 * \brief Return the number of uncomputed targets.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 * \return                Number of uncomputed targets
 */

int
CWP_N_uncomputed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id,
 const CWP_Dof_location_t pointsCloudLocation,
 const int  i_part
);

/**
 *
 * \brief Return uncomputed targets.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 * \return                Uncomputed targets
 */

const int *
CWP_Uncomputed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id
);

/**
 *
 * \brief Return the number of computed targets. <b>(Not implemented yet)</b>
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 * \return                Number of computed targets
 */

int
CWP_N_computed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id
);

/**
 *
 * \brief Return computed targets. <b>(Not implemented yet)</b>
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 * \return                Computed targets
 */

const int *
CWP_Computed_tgts_get
(
 const char *local_code_name,
 const char *cpl_id
);

/**
 * \brief Return distance from each target to the source interface. <b>(Not implemented yet)</b>
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 *
 * \return               Distance
 *
 */

const double *
CWP_Computed_tgts_dist_to_spatial_interp_get
(
 const char *local_code_name,
 const char *cpl_id
);

//----------------------------------------------------------------------------
// Functions about exchange frequency
//----------------------------------------------------------------------------

/**
 * \brief Setting receiving frequency.
 *
 * This function set receiving frequency. It must be used when
 * the type of receiving frequency is \ref CWP_TIME_EXCH_RELATED_N_TIME_STEP
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  n_step           Frequency in steps number
 *
 */

void
CWP_Recv_freq_set
(
 const char      *local_code_name,
 const char      *cpl_id,
 const int        n_step
);

/**
 * \brief Set the next receiving time.
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
CWP_Next_recv_time_set
(
 const char      *local_code_name,
 const char      *cpl_id,
 const double     next_time
);


/**
 * \brief Set the coupling time step. <b>(Not implemented yet)</b>
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
CWP_Cpl_time_step_set
(
 const char      *local_code_name,
 const char      *cpl_id,
 const int        next_time_step
);

/*----------------------------------------------------------------------------*
 * Functions about spatial interpolation                                      *
 *----------------------------------------------------------------------------*/

/**
 * \brief Compute spatial interpolation weights
 *
 * \param [in]  local_code_name     Local code name
 * \param [in]  cpl_id              Coupling identifier
 *
 */

void
CWP_Spatial_interp_weights_compute
(
 const char     *local_code_name,
 const char        *cpl_id
);

/**
 * \brief Set the properties of the spatial interpolation algorithm
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  fmt              Format with the syntax : "prop1, prop2, ..."
 * \param       ...              Values of each properties
 *
 */

void
CWP_spatial_interp_properties_set
(
 const char     *local_code_name,
 const char     *cpl_id,
 const char     *fmt,
 ...
);

/*----------------------------------------------------------------------------*
 * Functions about visualization                                              *
 *----------------------------------------------------------------------------*/

/**
 * \brief Enable visualization output
 *
 * This function enable visualization output.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  freq             Output frequency
 * \param [in]  format           Output format to visualize exchanged fieldsDouble
 *                               on the coupled mesh. Choice between :
 *                               - "EnSight Gold"
 *                               - "MED_ficher"
 *                               - "CGNS"
 *                               .
 * \param [in]  format_option   Output options "opt1, opt2, ..."
 *                               - text               output text files
 *                               - binary             output binary files (default)
 *                               - big_endian         force binary files
 *                                                    to big-endian
 *                               - discard_polygons   do not output polygons
 *                                                    or related values
 *                               - discard_polyhedra  do not output polyhedra
 *                                                    or related values
 *                               - divide_polygons    tesselate polygons
 *                                                    with triangles
 *                               - divide_polyhedra   tesselate polyhedra
 *                                                    with tetrahedra and pyramids
 *                                                    (adding a vertex near
 *                                                    each polyhedron's center)
 *                               .
 *
 */

void
 CWP_Visu_set
 (const char                 *local_code_name,
  const char                 *cpl_id,
  const int                   freq,
  const CWP_Visu_format_t     format,
  const char                 *format_option
 );

/*----------------------------------------------------------------------------*
 * Functions about User target points                                         *
 *----------------------------------------------------------------------------*/

/**
 * \brief Setting user target points
 *
 * This function must be called if the nature of receiving fieldsDouble
 * is \ref CWP_DOF_LOCATION_USER
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  n_pts            Number of points
 * \param [in]  coord            Coordinates (size = 3 * n_pts)
 *
 */

void
CWP_User_tgt_pts_set
(
 const char    *local_code_name,
 const char    *cpl_id,
 const int      i_part,
 const int      n_pts,
 double         coord[]
);

/*----------------------------------------------------------------------------*
 * Functions about Mesh                                                    *
 *----------------------------------------------------------------------------*/


/**
 * \brief Finalize interface mesh.
 *
 * This function computes the global numbers of mesh entities if they are
 * not provided.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 *
 */

void
CWP_Mesh_interf_finalize
(
 const char           *local_code_name,
 const char           *cpl_id
);

/**
 * \brief Set vertices.
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
CWP_Mesh_interf_vtx_set
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             n_pts,
 double                coord[],
 CWP_g_num_t           global_num[]
);

/**
 * \brief Add a connectivity block to the interface mesh.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  block_type       Block type
 *
 * \return block identifier
 */

int
CWP_Mesh_interf_block_add
(
 const char           *local_code_name,
 const char           *cpl_id,
 const CWP_Block_t     block_type
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
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  block_id         Block identifier
 * \param [in]  n_elts           Number of elements
 * \param [in]  connec           Connectivity (size = n_vertex_elt * n_elts)
 * \param [in]  global_num       Pointer to global element number (or NULL)
 *
 */

void
CWP_Mesh_interf_block_std_set
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
 * \brief Set a generic high order block to the interface mesh
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

/*
void
CWP_Mesh_interf_h_order_block_set
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
*/

/**
 * \brief Set the connectivity of a polygon block in a interface mesh partition.
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
CWP_Mesh_interf_f_poly_block_set
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
 * \brief Adding a polyhedron connectivity block to the interface mesh.
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
CWP_Mesh_interf_c_poly_block_set
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
 * \brief Delete interface mesh.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 */

void
CWP_Mesh_interf_del
(
 const char *local_code_name,
 const char *cpl_id
);


/**
 * \brief Map a PDM mesh nodal as interface mesh.
 *
 * \param [in] local_code_name   Local code name
 * \param [in] cpl_id            Coupling identifier
 * \param [in] i_part            Current partition
 * \param [in] pdm_nodal         pdm nodal mesh
 *
 */

void
CWP_Mesh_interf_shared_fvm_nodal
(
 const char   *local_code_name,
 const char   *cpl_id,
 const int     i_part,
 void         *pdm_nodal
);

/**
 * \brief Define the interface mesh from a cell to face connectivity.
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
 * \param [in]  parent_num        Pointer to parent element number (or NULL)
 *
 */

void
CWP_Mesh_interf_from_cellface_set
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
 CWP_g_num_t           parent_num[]
);


/**
 * \brief Define the surface interface mesh from a face to edge connectivity.
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
 * \param [in]  parent_num        Pointer to parent element number (or NULL)
 *
 */

void
CWP_Mesh_interf_from_faceedge_set
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
 CWP_g_num_t     parent_num[]
);

/*----------------------------------------------------------------------------*
 * Functions about field                                                      *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Create a new field.
 *
 * \param [in] local_code_name Local code name
 * \param [in]  cpl_id         Coupling identifier
 * \param [in]  field_id       Field id
 * \param [in]  data_type      Data type
 * \param [in]  storage        Storage type
 * \param [in]  n_component    Number of component
 * \param [in]  nature         Value location
 * \param [in]  exch_type      Exchange type
 * \param [in]  visu_status    Visualization status
 *
 */

void
CWP_Field_create
(
 const char                  *local_code_name,
 const char                  *cpl_id,
 const char                  *field_id,
 const CWP_Type_t             data_type,
 const CWP_Field_storage_t    storage,
 const int                    n_component,
 const CWP_Dof_location_t      value_target_location,
 const CWP_Field_exch_t       exch_type,
 const CWP_Status_t           visu_status
);


/**
 *
 * \brief Set field data.
 *
 * \param [in] local_code_name   Local code name
 * \param [in] cpl_id            Coupling identifier
 * \param [in] field_id          Field identifier
 * \param [in] i_part            Current partition
 * \param [in] data              Storage array (Mapping)
 *
 */

void
CWP_Field_data_set
(
 const char         *local_code_name,
 const char         *cpl_id,
 const char         *field_id,
 const int           i_part,
 double              data[]
);


/**
 *
 * \brief Set field gradient (optional).  <b>(Not implemented yet)</b>
 *
 * \param [in] local_code_name Local code name
 * \param [in] cpl_id          Coupling identifier
 * \param [in] field_id        Field identifier
 * \param [in] i_part          Current partition
 * \param [in] order           Order
 * \param [in] data            Storage array (Mapping)
 *
 */

void
CWP_Field_gradient_data_set
(
 const char      *local_code_name,
 const char      *cpl_id,
 const char      *field_id,
 const int        i_part,
 const int        order,
 double           data[]
);

/**
 *
 * \brief Get number of field components.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      number of field components
 *
 */

int
CWP_Field_n_component_get
(
 const char      *local_code_name,
 const char      *cpl_id,
 const char      *field_id
);

/**
 *
 * \brief Get target degrees of freedom location.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      Field nature
 *
 */

CWP_Dof_location_t
CWP_Field_target_dof_location_get
(
 const char      *local_code_name,
 const char      *cpl_id,
 const char      *field_id
);


/**
 *
 * \brief Get field data type.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      Field data type
 *
 */

CWP_Type_t
CWP_Field_data_type_get
(
 const char      *local_code_name,
 const char      *cpl_id         ,
 const char      *field_id
);

/**
 *
 * \brief Get field storage type.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      Field storage type
 */

CWP_Field_storage_t
CWP_Field_storage_get
(
 const char      *local_code_name,
 const char      *cpl_id         ,
 const char      *field_id
);

/**
 * \brief Delete a field.
 *
 * \param [in] local_code_name Local code name
 * \param [in]  cpl_id         Coupling identifier
 * \param [in]  field_id       Field identifier
 *
 */

void
CWP_Field_del
(
 const char      *local_code_name,
 const char      *cpl_id         ,
 const char      *field_id
);

/*----------------------------------------------------------------------------*
 * Functions about exchange                                                   *
 *----------------------------------------------------------------------------*/

/**
 * \brief Exchange spatially interpolated fields. <b>(Not implemented yet)</b>
 *
 * This function exchanges the interpolated fields for each coupling depending
 * on mode of time exchange \ref CWP_Time_exch_t.
 *
 * \param [in] local_code_name      Local code name
 * \param [in] cpl_id               Coupling identifier
 *
 */

void
CWP_Field_Exch
(
 const char *local_code_name,
 const char *cpl_id
);


/**
 * \brief Send a spatially interpolated field to the coupled code with
 *        nonblocking communications.
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
CWP_Field_Issend
(
 const char     *local_code_name,
 const char     *cpl_id,
 const char     *src_field_id
);

/**
 *
 * \brief Receive a spatially interpolated field from the coupled code
 *        with nonblocking communications.
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
CWP_Field_Irecv
(
 const char        *local_code_name,
 const char        *cpl_id,
 const char        *targetFieldID
);

/**
 *
 * \brief Waiting of the end of exchange related to request.
 *
 * This function waits the end of exchange related to request
 * from \ref CWP_Issend
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 */

void
CWP_Wait_issend
(
 const char  *local_code_name,
 const char  *cpl_id,
 const char  *src_field_id
);

/**
 *
 * \brief Waiting of the end of exchange related to request.
 *
 * This function waits the end of exchange related to request
 * from \ref CWP_Irecv
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 */

void
CWP_Wait_irecv
(
 const char  *local_code_name,
 const char  *cpl_id,
 const char  *distant_field_id
);


 CWP_g_num_t*
 CWP_GlobalNumGet
 (
  const char  *local_code_name,
  const char  *cpl_id,
  const int    id_block,
  const int    i_part
 );


/*----------------------------------------------------------------------------*
 * Functions about user interpolation                                         *
 *----------------------------------------------------------------------------*/

/* /TODO: associer l'interpolation utilisateur au champ plutot qu'au couplage */

/**
 *
 * \brief Setting of an user interpolation from location.
 *
 * This function takes into account an user interpolation function written with
 * void (*\ref CWP_Interp_from_location_t) interface.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] fct              Function
 *
 */

void
CWP_Interp_from_location_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 const char                 *field_id,
 CWP_Interp_from_location_t  fct
);

/**
 *
 * \brief Setting of an user interpolation from intersection. <b>(Not implemented yet)</b>
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
CWP_Interp_from_intersect_set
(
 const char                *local_code_name,
 const char                *cpl_id,
 CWP_Interp_from_intersect_t fct
);

/**
 *
 * \brief Setting of an user interpolation from closest points. <b>(Not implemented yet)</b>
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
CWP_Interp_from_closest_pts_set
(
 const char                     *local_code_name,
 const char                     *cpl_id,
 CWP_Interp_from_closest_pts_t   fct
);

/*----------------------------------------------------------------------------*
 * Functions about current code control parameters                            *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Add a control parameter
 *
 * Addition of a control parameter in the code properties.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 * \param [in] initial_value    Initial value
 *
 */

void
CWP_Param_add
(
 const char        *local_code_name,
 const char        *param_name,
 const CWP_Type_t   data_type,
 void              *initial_value
);


/**
 *
 * \brief Set a control parameter
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 * \param [in] value            Value
 *
 */

void
CWP_Param_set
(
 const char             *local_code_name,
 const char             *param_name,
 const CWP_Type_t        data_type,
 void                   *value
);


/**
 *
 * \brief Removing a local int control parameter
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type,
 *
 */

void
CWP_Param_del
(
 const char       *local_code_name,
 const char       *param_name,
 const CWP_Type_t  data_type
);


/*----------------------------------------------------------------------------*
 * Functions about all code parameters                                        *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Return the number of parameters of a code
 *
 * \param [in] code_name       Local or distant code name
 * \param [in] data_type       Parameter type,
 *
 * return  Number of parameters
 *
 */

int
CWP_Param_n_get
(
 const char             *code_name,
 const CWP_Type_t        data_type
);

/**
 *
 * \brief Return the parameter list of a code
 *
 * \param [in]  code_name      Local or distant code name
 * \param [in]  data_type      Parameter type,
 * \param [out] nParam         Number of parameters
 * \param [out] paramNames     Parameter names
 *
 *
 */

void
CWP_Param_list_get
(
 const char             *code_name,
 const CWP_Type_t        data_type,
 int                    *nParam,
 char                 ***paramNames
);

/**
 *
 * \brief Is a parameter ?
 *
 * \param [in] code_name      Local or distant code name
 * \param [in] param_name     Parameter name
 * \param [in] data_type      Parameter type,
 *
 * return  1 : true / 0 : false
 *
 */

int
CWP_Param_is
(
 const char             *code_name,
 const char             *param_name,
 const CWP_Type_t        data_type
);

/**
 *
 * \brief Return the value of a int control parameter from "code_name" code
 *
 * \param [in]  code_name  Local or distant code name
 * \param [in]  param_name Parameter name
 * \param [in]  data_type  Parameter type
 * \param [out] value      Parameter value
 *
 */

void
CWP_Param_get
(
 const char       *code_name,
 const char       *name,
 const CWP_Type_t  data_type,
 void             *value
);

/**
 *
 * \brief Return the operation result on a control parameter
 *        (same name in all codes)
 *
 * \param [in]  op        Operation
 * \param [in]  name      Parameter name
 * \param [in]  data_type Parameter type,
 * \param [out] res       Result
 * \param [in]  nCode     Number of codes
 * \param       ...       Codes name
 *
 */

void
CWP_Param_reduce
(
 const CWP_Op_t    op,
 const char       *name,
 const CWP_Type_t  data_type,
 void             *res,
 const int         nCode,
 ...
);

/**
 *
 * \brief Lock access to local parameters from a distant code
 *
 * \param [in]  code_name  Code to lock
 *
 */

void
CWP_Param_lock
(
const char *code_name
);

/**
 *
 * \brief unlock access to local parameters from a distant code
 *
 * \param [in]  code_name  Code to unlock
 *
 */

void
CWP_Param_unlock
(
const char *code_name
);


void
CWP_surf_gen_init
(char* genName,
  int nx, int ny, int nPart, MPI_Comm* comm, double prop, double width, double randomVar
);

void
CWP_surf_gen_compute
(char* genName
);



void
CWP_surf_gen_by_block_get
( char* genName, int i_part,
  int* nVtx , double** coords, CWP_g_num_t** vtxGnum, int* nElts,
  int* nTri , int** eltsConnecTri , CWP_g_num_t** eltsGnumTri,
  int* nQuad, int** eltsConnecQuad, CWP_g_num_t** eltsGnumQuad,
  int* nPoly, int** eltsConnecPolyIndex, int** eltsConnecPoly, CWP_g_num_t** eltsGnumPoly
);

void
CWP_surf_gen_one_connectivity_get
( char* genName, int i_part,
  int* nVtx , double** coords, CWP_g_num_t** vtxGnum,
  int* nElts, int** eltsConnecIndex, int** eltsConnec, CWP_g_num_t** eltsGnum
);

void
CWP_surf_face_edge_get
( char* genName, int i_part,
  int* nVtx , double** coords, CWP_g_num_t** vtxGnum,
  int* nFace, int** faceEdgeIdx, int** faceEdge,
  int* nEdge, int** edgeVtxIdx, int** edgeVtx,
  CWP_g_num_t** faceLNToGN
);


void
CWP_surf_gen_tri_field_get
( char* genName, int i_part,
  double** field
);

void
CWP_surf_gen_quad_field_get
( char* genName, int i_part,
  double** field
);

void
CWP_surf_gen_poly_field_get
( char* genName, int i_part,
  double** field
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CWP_H__ */
