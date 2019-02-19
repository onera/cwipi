#ifndef __CWP_H__
#define __CWP_H__
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


/** \file cwp.h
  * \brief CWIPI new API header file 
  * 
  */


#include <mpi.h>
#include <stdio.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/
#include <fvmc_parall.h>
#include <fvmc_nodal.h>



/*=============================================================================
 * Macro definitions
 *============================================================================*/

#if !defined (__hpux) && !defined (_AIX) 
/**
  *  \def PROCF(x, y)
  *  \brief Macro to rename function for Fortran.
  *  \param x input value.
  *  \param y input value.
  *  \returns ?? 
  */
#define PROCF(x, y) x##_
#else
/**
  *  \def PROCF(x, y)
  *  \brief Macro to rename function for Fortran.
  *  \param x input value.
  *  \param y input value.
  *  \returns ?? 
  */
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
 * \brief  Type of parameters               
 *
 * CWP_Type_t gives the available parameter types                
 */

typedef enum {

  CWP_DOUBLE,      /*!< Double precision type */
  CWP_INT,         /*!< Integer type */ 
  CWP_CHAR         /*!< String type */ 

} CWP_Type_t;


/**
 * \enum CWP_Comm_t
 * \brief Communication mode
 *
 * CWP_Comm_t gives the different communication mode 
 */

typedef enum {

  CWP_COMM_PAR_WITH_PART,    /*!< Parallel communcation 
                                  on partitioned source */
  CWP_COMM_PAR_WITHOUT_PART, /*!< Parallel communcation 
                                  on unpartitioned source defined on 
                                  all process */
  CWP_COMM_SEQ,              /*!< Parallel communcation 
                                  on unpartitioned source defined on 
                                  master process */
  CWP_COMM_INTERNAL         /*!< Internal communcation within a process */

} CWP_Comm_t;

/**
 * \enum CWP_Freq_t
 * \brief  Echange frequency 
 *
 * CWP_Freq_t describes the different ways to define coupling data 
 */

typedef enum {

  CWP_FREQ_EACH_TIME_STEP,      /*!< Exchange at each time step */
  CWP_FREQ_RELATED_N_TIME_STEP, /*!< Exchange frequency is related to 
                                     the number of time steps  */ 
  CWP_FREQ_CPL_TIME_STEP,       /*!< Coupling time step        */ 
  CWP_FREQ_ASYNCHRONOUS,        /*!< Exchanges are asynchronous */ 
  CWP_FREQ_SLAVE,               /*!< Give a converged state    */ 
  CWP_FREQ_MASTER               /*!< Request a converged state */ 

} CWP_Freq_t;


/**
 * \enum CWP_Field_value_t
 * \brief Field nature
 *
 * CWP_Field_value_t gives different nature 
 */

typedef enum {

  CWP_FIELD_VALUE_CELL_MEAN,   /*!< Cell mean value */
  CWP_FIELD_VALUE_CELL_POINT,  /*!< Field defined in cell point  */
  CWP_FIELD_VALUE_NODE,        /*!< Cell vertex field */
  CWP_FIELD_VALUE_USER         /*!< User defined field */ 

} CWP_Field_value_t ;


/**
 * \enum CWP_Field_exch_t
 * \brief Field exchange type
 *
 * CWP_Field_exch_t gives different type of field exchange
 */

typedef enum {

  CWP_FIELD_EXCH_SEND,        /*!< Send */
  CWP_FIELD_EXCH_RECV,        /*!< Receive */
  CWP_FIELD_EXCH_SENDRECV    /*!< Send and receive */ 

} CWP_Field_exch_t ;


/**
 * \enum CWP_Field_storage_t
 * \brief Field storage
 *
 * CWP_Field_storage_t gives dfieferent nature 
 */

typedef enum {

  CWP_FIELD_STORAGE_INTERLACED,  /*!< Interlaced storage */
  CWP_FIELD_STORAGE_BLOCK        /*!< Block storage */

} CWP_Field_storage_t ;


/* TODO: Voir si on met cet attribut sur FIELd */

/**
 * \enum CWP_Interpolation_t
 * \brief Interpolation type
 *
 * CWP_Interpolation_t gives the different ways to interpolate
 */

typedef enum {

  CWP_INTERPOLATION_DEFAULT,  /*!< Default interpolation */
  CWP_INTERPOLATION_USER      /*!< User interpolation */

} CWP_Interpolation_t;

/**
 * \enum CWP_Status_t
 * \brief Status on/off
 *
 * CWP_Status_t defines the different status
 */

typedef enum {

  CWP_STATUS_OFF,           /*!< output with error */
  CWP_STATUS_ON             /*!< Output without error */

} CWP_Status_t;

/**
 * \enum CWP_Err_t
 * \brief Error codes
 *
 * CWP_Err_t defines the different error codes 
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
 * \enum CWP_Displacement_t
 * \brief Active support moving
 *
 * CWP_Displacement_t actives support moving (mesh or point cloud)  
 */

typedef enum {

  CWP_DISPLACEMENT_STATIC,       /*!< Static */ 
  CWP_DISPLACEMENT_UNSPECIFIED,  /*!< Unspecified displacement */ 
  CWP_DISPLACEMENT_TRANSLATION,  /*!< Translation */ 
  CWP_DISPLACEMENT_ROTATION      /*!< Rotation */ 

} CWP_Displacement_t;

/**
 * \enum CWP_Geom_t
 * \brief Geomtric algorithms
 *
 * CWP_Geom_t gives different geometric algorithm on which interpolation 
 * method is based 
 */

typedef enum {

  CWP_GEOM_CLOSEST_POINT, /*!< Closest points */
  CWP_GEOM_INTERSECTION,  /*!< Meshes intersection */
  CWP_GEOM_LOCATION       /*!< Location into a mesh */

} CWP_Geom_t;

/**
 * \enum CWP_Interface_t
 * \brief Coupling interfaces
 *
 * CWP_Interface_t gives different coupling interfaces 
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
 * CWP_State_t gives different code states. 
 * 
 */

typedef enum {

  CWP_STATE_IN_PROGRESS,         /*!< Code running */
  CWP_STATE_END,                 /*!< Computation end */
  CWP_STATE_OUTPUT_ERROR         /*!< Output one error */

} CWP_State_t;

/**
 * \enum CWP_Op_t
 * \brief Operations on crontrol parameters
 *
 */

typedef enum {

  CWP_OP_MIN,            /*!< output with error */
  CWP_OP_MAX,            /*!< Output without error */
  CWP_OP_SUM             /*!< Output without error */

} CWP_Op_t;

/*============================================================================
 * User interpolation type
 *============================================================================*/

/**
 * \typedef void (*CWP_Interp_from_location_t)
 * \brief User interpolation function from location into a mesh.
 *
 * void (*CWP_Interp_from_location_t) defines the user interpolation 
 * interface to take into account an user interpolation from location of target 
 * points into the source mesh.
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
 * \param [in]  tgt_pts_location            target points location
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
 * \param [in]  src_field_location          source field location
 * \param [in]  src_field                   source field
 *                                          (size depends on field type and stride)
 * \param [in]  src_field_location          target field location
 * \param [out] tgt_field                   target field
 *                                          (size = stride * n_tgt_pts)
 */

typedef void (*CWP_Interp_from_location_t)
(
 const int                   interface_type,
 const int                   n_src_vtcs,
 const int                   n_src_std_elts,
 const int                   n_src_poly,
 const int                   n_tgt_pts,
 const double                src_vtcs_coords[],
 const CWP_g_num_t            src_parent_elts_num[],
 const CWP_g_num_t            src_parent_vtcs_num[],
 const int                   src_connec_idx[],
 const int                   src_connec[],
 const int                   src_poly_cell_face_idx[],
 const int                   src_poly_cell_face_connec[],
 const int                   src_poly_face_vtx_idx[],
 const int                   src_poly_face_vtx_connec[],
 const double                tgt_pts_coords[],
 const int                   tgt_pts_location[],
 const float                 tgt_pts_dist[],
 const int                   tgt_pts_bary_coords_idx[],
 const double                tgt_pts_bary_coords[],
 const int                   stride,
 const CWP_Field_value_t     src_field_location,
 const void                 *src_field,
 const CWP_Field_value_t     tgt_field_location,
 void                       *tgt_field
);

/**
 * \typedef void (*CWP_Interp_from_intersec_t)
 * \brief User interpolation function from intersection between meshes <b>(Not implemented yet)</b>  
 *
 * void (*CWP_Interp_from_intersec_t) defines the user interpolation 
 * interface to take into account an user interpolation from intersection 
 * between source and target meshes 
 *
 * \param [in]  interface_type              Interface type
 *
 */

typedef void (*CWP_Interp_from_intersec_t)
(
 const int interface_type
);

/**
 * \typedef void (*CWP_Interp_from_closest_pts_t)
 * \brief User interpolation function from closest points <b>(Not implemented yet)</b>
 *
 * void (*CWP_Interp_from_closest_pts_t) defines the user interpolation 
 * interface to take into account an user interpolation from <i>n</i> closest
 * point 
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
 * This function create the MPI intra communicator for this code from
 * the MPI inter communicator that contains all code process. It is a
 * synchronization point between all codes
 *
 * \param [in]  global_comm       MPI global communicator
 * \param [in]  n_code            Number of codes on the current rank
 * \param [in]  code_names         Names of codes on the current rank (size = n_code)
 * \param [in]  is_coupled_rank   Is current rank used for coupling (size = n_code)
 * \param [in]  time_init         Time init (size = n_code)
 * \param [out] intra_comms        MPI intra communicators of each code
 *
 */

void 
CWP_Init
(
 const MPI_Comm           global_comm,
 const int                n_code,
 const char             **code_names,
 const CWP_Status_t      *is_coupled_rank,
 const double            *time_init,
 MPI_Comm                *intra_comms
);

/**
 *
 * \brief CWIPI completion
 *
 * This function finalize CWIPI
 *
 */

void 
CWP_Finalize
(
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
 * \brief Writing output to file.
 *
 * This function set the file for writing output.
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
 * \brief Number of codes known to CWIPI
 *
 * \return Number of codes
 *
 */

int
CWP_Codes_nb_get
(
);


/**
 * \brief list of codes known to CWIPI
 *
 * \return Names list of codes
 *
 */

const char **
CWP_Codes_list_get
(
void
);


/**
 * \brief Number of codes known to CWIPI
 *
 * \return Number of local codes
 *
 */

int
CWP_Loc_codes_nb_get
(
);


/**
 * \brief list of codes known to CWIPI
 *
 * \return Names list of local codes
 *
 */

const char **
CWP_Loc_codes_list_get
(
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
(
);

/*----------------------------------------------------------------------------*
 * General functions about coupling                                           *
 *----------------------------------------------------------------------------*/

/**
 * \brief Creating a coupling object.
 *
 * This function creates a coupling object and defines its properties.
 *
 * \param [in]  local_code_name     Local code name
 * \param [in]  cpl_id              Coupling identifier
 * \param [in]  coupled_code_name   Distant or local coupled code name
 * \param [in]  comm_type           Communication type
 * \param [in]  geom_algo           Geometric algorithm
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
 const CWP_Geom_t          geom_algo,
 const int                 n_part,
 const CWP_Displacement_t  displacement,   
 const CWP_Freq_t          recv_freq_type 
);


/**
 * \brief Initialize translation displacement.
 * 
 * This function defines the translation direction. The movement is updated 
 * before the end of the current time step by \ref CWP_Cpl_trans_update 
 * Two coupled codes have to define the same properties. The distant code is always
 * considered as a static interface 
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
 * \brief Define the next time step position.
 * 
 * This function compute the next time step position from a relative distance. If
 * it is a known position, Geometric algorithm is not reprocessed. Otherwise, 
 * geometric algorithm is launched from previous results. 
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
 * \brief Initialize translation displacement.
 * 
 * This function defines the rotation properties. The movement is updated 
 * before the end of the current time step by \ref CWP_Cpl_rotation_update. 
 * Two coupled codes have to define the same properties. The distant code is always
 * considered as a static interface 
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
 * \brief Define the next time step position.
 * 
 * This function compute the next time step position from a relative angle. If
 * it is a known position, Geometric algorithm is not reprocessed. Otherwise, 
 * geometric algorithm is launched from previous results. 
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
 * \brief Set storage properties                 
 *
 * This functions activates the storage of geometric results in case of rotation 
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
 * \brief Removing a coupling object
 *
 * This function delete a coupling abject
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
 * \brief Return the number of uncomputed targets
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
 const char *cpl_id
);

/**
 *
 * \brief Return uncomputed targets
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
 * \brief Return the number of computed targets
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
 * \brief Return computed targets
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
 * \brief Return distance from each target to the geometric interface                 
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 *
 * \return               Distance
 *
 */

const double *
CWP_Computed_tgts_dist_to_geom_get
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
 * the type of receiving frequency is \ref CWP_FREQ_RELATED_N_TIME_STEP
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
 * \brief Setting the next receiving time.
 *
 * This function set the next receiving time. It must be used when
 * the type of receiving frequency is \ref CWP_FREQ_ASYNCHRONOUS
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
 * \brief Setting the coupling time step.
 *
 * This function set the coupling time step. It must be used when
 * the type of receiving frequency is \ref CWP_FREQ_CPL_TIME_STEP
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  next_time_step        Coupling time step
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
 * Functions about geometry                                                   *
 *----------------------------------------------------------------------------*/

/**
 * \brief Computation geometry                                  
 *
 * This function compute geometry 
 *
 * \param [in]  local_code_name     Local code name
 * \param [in]  cpl_id              Coupling identifier
 * \param [out] n_uncomputed_tgt    Number of uncomputed target
 *
 */

void 
CWP_Geom_compute
(
 const char     *local_code_name,
 const char     *cpl_id,
 int            *n_uncomputed_tgt
);

/**
 * \brief Load geom results from it id
 *
 * This function set the location algorithm properties. It must be only used
 * when the type of geometric algorithm is \ref CWP_GEOM_LOCATION
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  fmt              Format with the syntax : "prop1, prop2, ..."
 * \param       ...              Values of each properties
 *
 */

void 
CWP_Geom_properties_set
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
 * \param [in]  format           Output format to visualize exchanged fields
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
(
 const char    *local_code_name,
 const char    *cpl_id,
 const int      freq,
 const char    *format,
 const char    *format_option
);

/*----------------------------------------------------------------------------*
 * Functions about User target points                                         *
 *----------------------------------------------------------------------------*/

/**
 * \brief Setting user target points
 *
 * This function must be called if the nature of receiving fields 
 * is \ref CWP_FIELD_VALUE_USER
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
 const int      n_pts,
 double         coord[]
);

/*----------------------------------------------------------------------------*
 * Functions about Mesh                                                    *
 *----------------------------------------------------------------------------*/

/**
 * \brief Setting vertices
 *
 * This function set partition vertices
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  n_pts            Number of points
 * \param [in]  coord            Coordinates (size = 3 * n_pts)          
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
 * \brief End setting of the mesh
 *
 * This function finalizes the mesh building after addition of block and coordinates setting.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 *
 */

void 
CWP_Mesh_interf_end_set
(
 const char        *local_code_name,
 const char        *cpl_id
);

/**
 * \brief Adding a connectivity block to the geometric support
 *
 * This function adds a connectivity block to the geometric support.
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
 * \param [in]  i_part           Current partition
 * \param [in]  block_type       Block type
 * \param [in]  n_elts           Number of elements
 * \param [in]  connec           Connectivity (size = n_vertex_elt * n_elts)          
 * \param [in]  global_num       Pointer to parent element number (or NULL)
 *
 */

void 
CWP_Mesh_interf_std_block_add
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const CWP_Block_t  block_type,
 const int          n_elts,
 int          connec[],
 CWP_g_num_t  global_num[]
);


/**
 * \brief Add a generic high order elements block
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  block_type       Block type
 * \param [in]  n_elts           Number of elements
 * \param [in]  order            Element order
 * \param [in]  connec           Connectivity (size = n_vertex_elt * n_elts)          
 * \param [in]  global_num       Pointer to parent element number (or NULL)
 *
 */


void 
CWP_Mesh_interf_h_order_block_add
(
 const char        *local_code_name,
 const char        *cpl_id,
 const int          i_part,
 const CWP_Block_t  block_type,
 const int          n_elts,
 const int          order, 
 int                connec[],
 CWP_g_num_t        global_num[]
);


/**
 * \brief Adding a polygon connectivity block to the mesh interface
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  n_elts           Number of elements
 * \param [in]  connec_idx       Connectivity index (connec_id[0] = 0 and 
 *                               size = n_elts + 1)          
 * \param [in]  connec           Connectivity (size = connec_idx[n_elts])          
 * \param [in]  parent_num       Pointer to parent element number (or NULL)
 *
 */

void 
CWP_Mesh_interf_f_poly_block_add
(
 const char             *local_code_name,
 const char             *cpl_id,
 const int               i_part,
 const int               n_elts,
 int                     connec_idx[],
 int                     connec[],
 CWP_g_num_t             parent_num[]
);

/**
 * \brief Adding a polyhedron connectivity block to the interface mesh
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  n_elts            Number of elements
 * \param [in]  cell_face_idx     Polyhedron to face index 
 *                                (cell_face_idx[0] = 0 and
 *                                size = n_elts + 1)
 * \param [in]  cell_face         Polyhedron to face connectivity 
 *                                (size = cell_face_idx[n_elts])
 * \param [in]  n_faces           Number of faces      
 * \param [in]  face_vtx_idx      Polyhedron face to vertex index 
 *                                (face_vtx_idx[0] = 0 and
 *                                 size = n_faces + 1
 * \param [in]  face_vtx          Polyhedron face to vertex connectivity
 *                                (size = face_vtx_idx[n_faces])
 * \param [in]  parent_num        Pointer to parent element number (or NULL)
 *
 */

void 
CWP_Mesh_interf_c_poly_block_add
(
 const char           *local_code_name,
 const char           *cpl_id,
 const int             i_part,
 const int             n_elts,
 int                   cell_face_idx[],
 int                   cell_face[],
 const int             n_faces,
 int                   face_vtx_idx[],
 int                   face_vtx[],
 CWP_g_num_t           parent_num[]
);

/**
 * \brief Interface mesh delation                                  
 * *
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
 * \brief Map a fvm nodal as mesh interface                                
 *
 * This function  map a fvm nodal as mesh interface
 *
 * \param [in] local_code_name   Local code name
 * \param [in] cpl_id            Coupling identifier
 * \param [in] i_part            Current partition
 * \param [in] fvmc_nodal        fvm nodal mes     
 *
 */

void 
CWP_Mesh_interf_shared_fvm_nodal
(
 const char   *local_code_name,
 const char   *cpl_id,
 const int     i_part,
 fvmc_nodal_t *fvmc_nodal
);


/**
 * \brief Define the volume interface mesh from a cell to face connectivity 
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  n_cells           Number of cells
 * \param [in]  cell_face_idx     Polyhedron to face index 
 *                                (src_poly_cell_face_idx[0] = 0 and
 *                                 size = n_elts + 1)
 * \param [in]  cell_face         Cell to face connectivity 
 *                                (size = cell_face_idx[n_elts])
 * \param [in]  n_faces           Number of faces      
 * \param [in]  face_vtx_idx      Polyhedron face to vertex index 
 *                                (face_vertex_idx[0] = 0 and
 *                                 size_idx = max(cell_face_connec) + 1)
 * \param [in]  face_vtx          Face to vertex connectivity
 *                                (size = face_vertex_idx[size_idx - 1])
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
 * \brief Define the surface interface mesh from a face to edge connectivity 
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  n_faces           Number of cells
 * \param [in]  face_edge_idx     Polygon to edge index 
 *                                (face_edge_idx[0] = 0 and
 *                                 size =  n_faces + 1)
 * \param [in]  face_edge         Face to edge connectivity 
 *                                (size = face_edge_idx[n_faces])
 * \param [in]  n_edges           Number of faces      
 * \param [in]  edge_vtx_idx      Polyhedron face to vertex index 
 *                                (edge_vtx_idx[0] = 0 and
 *                                 size_idx = max(edge__connec) + 1)
 * \param [in]  edge_vtx          Face to vertex connectivity
 *                                (size = edge_vtx_idx[n_edges])
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
 * \brief Creating a new field
 * 
 * \param [in] local_code_name Local code name
 * \param [in]  cpl_id         Coupling identifier
 * \param [in]  field_id       Field id
 * \param [in]  data_type      Data type          
 * \param [in]  storage        Storage type          
 * \param [in]  n_component    Number of componenent
 * \param [in]  nature         Nature
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
 const CWP_Field_value_t      nature,
 const CWP_Field_exch_t       exch_type,
 const CWP_Status_t           visu_status
);


/**
 *
 * \brief Set data mapping
 * 
 * \param [in] local_code_name   Local code name
 * \param [in] cpl_id            Coupling identifier
 * \param [in] field_id          Field identifier
 * \param [in] i_part            Current partition
 * \param [in] data              Storage array (Mapping)
 * 
 */

void
CWP_Field_mapping_set
(
 const char      *local_code_name,
 const char      *cpl_id,
 const char      *field_id,
 const int        i_part,
 double           data[]
);


/**
 *
 * \brief Set data mapping
 * 
 * TODO Define gradient storage
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
CWP_Field_gradient_mapping_set
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
 * \brief Get nunmber of field components
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
 * \brief Get field nature
 * 
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      Field nature
 * 
 */

CWP_Field_value_t
CWP_Field_location_get
(
 const char      *local_code_name,
 const char      *cpl_id,
 const char      *field_id
);


/**
 *
 * \brief Get field data type
 * 
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      Field data type
 * 
 */

 CWP_Type_t
 CWP_Field_type_get
 (
  const char                  *local_code_name,
  const char                  *cpl_id,
  const char                  *field_id
  );


/**
 *
 * \brief Get field data type
 * 
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      Field data type
 * 
 */

CWP_Type_t
CWP_Field_value_type_get
(
 const char      *local_code_name,
 const char      *cpl_id,
 const char      *field_id
);

/**
 *
 * \brief Get field storage type
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
 const char      *cpl_id,
 const char      *field_id
);

/**
 *
 * \brief Removing a field
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
 const char      *cpl_id,
 const char      *field_id
);

/*----------------------------------------------------------------------------*
 * Functions about exchange                                                   *
 *----------------------------------------------------------------------------*/


/**
 * \brief data exchange <b>(Not implemented yet)</b> 
 *
 * This function exchanges for each coupling depending on exchange frequency
 * 
 * \param [in] local_code_name      Local code name
 * \param [in] cpl_id     Coupling identifier
 *
 */

void
CWP_Exch
(
 const char *local_code_name,
 const char *cpl_id
 );

/**
 *
 * \brief Exchange data field with the coupled code with blocking 
 *        communications.
 *
 * This function exchanges interpolated fields between coupled codes. 
 * 
 * \warning  The size of tgt_field_id size is n_computed_tgt. 
 *           If \f$ n\_uncomputed\_tgt \ne n\_tgt\_pts \f$,
 *           user himself must set values for uncomputed target points.
 *
 * \param [in] local_code_name      Local code name
 * \param [in]  cpl_id              Coupling identifier
 * \param [in]  src_field_id        Source field id (0 -> no sending)
 * \param [in]  tgt_field_id        Target field id (0 -> no receiving)
 * \param [out] n_uncomputed_tgt    Number of uncomputed target for each partition (size n_part)
 *
 * \return                          Exchange status
 *
 */

CWP_Err_t 
CWP_Sendrecv
(
 const char   *local_code_name,
 const char   *cpl_id,
 const char   *src_field_id,
 const char   *tgt_field_id,
 int          *n_uncomputed_tgt
);

/**
 *
 * \brief Sending of data field to the coupled code with nonblocking 
 *        communications.
 *
 * This function sends interpolated field to the coupled code. 
 * 
 * \param [in] local_code_name  Local code name
 * \param [in]  cpl_id          Coupling identifier
 * \param [in]  src_field_id    Source field id
 *
 * \param [out] request         Request to call by \ref CWP_Wait_issend 
 *                              to wait the end of exchange
 *
 */

void 
CWP_Issend
(
 const char     *local_code_name,
 const char     *cpl_id,
 const char     *src_field_id,
 int            *request
);

/**
 *
 * \brief Receiving of Data field from the coupled code with nonblocking 
 *        communications.
 *
 * This function receives interpolated field from the coupled code 
 * 
 * \param [in] local_code_name  Local code name
 * \param [in]  cpl_id          Coupling identifier
 * \param [in]  tgt_field_id    Target field id
 *
 * \param [out] request         Request to call by \ref CWP_Wait_irecv  
 *                              to wait the end of exchange
 *
 */

void 
CWP_Irecv
(
 const char   *local_code_name,
 const char   *cpl_id,
 const char   *tgt_field_id,
 int          *request
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
 * \param [in] request          Request to wait the end of exchange
 *
 */

void 
CWP_Wait_issend
(
 const char  *local_code_name,
 const char  *cpl_id,
 int          request
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
 * \param [in] request          Request to wait the end of exchange
 *
 */

void 
CWP_Wait_irecv
(
 const char  *local_code_name,
 const char  *cpl_id,
 int          request
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
CWP_Interp_from_loc_set
(
 const char                 *local_code_name,
 const char                 *cpl_id,
 CWP_Interp_from_location_t fct
);

/**
 *
 * \brief Setting of a FORTRAN user interpolation from location.
 *
 * This function takes into account an user interpolation function written
 * in FORTRAN.
 * 
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] fct              Function
 *
 */

void 
CWP_Interp_from_loc_set_f
(
 const char *local_code_name,
 const char *cpl_id,
 void       *fct
);

/**
 *
 * \brief Setting of an user interpolation from intersection.
 *
 * This function takes into account an user interpolation function written with
 * void (*\ref CWP_Interp_from_intersec_t) interface.
 * 
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] fct              Function
 *
 */

void 
CWP_Interp_from_inter_set
(
 const char                *local_code_name,
 const char                *cpl_id,
 CWP_Interp_from_intersec_t fct
);

/**
 *
 * \brief Setting of a FORTRAN user interpolation from intersection.
 *
 * This function takes into account an user interpolation function written
 * in FORTRAN .
 * 
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] fct              Function
 *
 */

void 
CWP_Interp_from_inter_set_f
(
 const char *local_code_name,
 const char *cpl_id,
 void       *fct
);

/**
 *
 * \brief Setting of an user interpolation from closest points
 *
 * This function takes into account an user interpolation function written with
 *  void (*\ref CWP_Interp_from_intersec_t) interface.
 * 
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] fct              Function
 *
 */

void 
CWP_Interp_from_closest_set
(
 const char                     *local_code_name,
 const char                     *cpl_id,
 CWP_Interp_from_closest_pts_t   fct
);

/**
 *
 * \brief Setting of a FORTRAN user interpolation from closest points
 *
 * This function takes into account an user interpolation function written
 * in FORTRAN .
 * 
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] fct              Function
 *
 */

void 
CWP_Interp_from_closest_set_f
(
 const char *local_code_name,
 const char *cpl_id,
 void       *fct
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

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CWP_H__ */
