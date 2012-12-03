#ifndef __CWIPI_H__
#define __CWIPI_H__
/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011-2012  ONERA

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
 * \enum CWIPI_long_t
 * \brief Long int in cwipi 
 *
 * CWIPI_data_t describes the different ways to define coupling data 
 */

typedef long cwipi_long_int_t;

/*============================================================================
 * Enumeration definitions
 *============================================================================*/

/**
 * \enum CWIPI_data_t
 * \brief  How coupling data are defined 
 *
 * CWIPI_data_t describes the different ways to define coupling data 
 */

typedef enum {

  CWIPI_DATA_FROM_XML,      /*!< Data are defined in a XML data file */
  CWIPI_DATA_IN_CODE_SRC,   /*!< Data are Defined in code source */ 

} CWIPI_data_t;

/**
 * \enum CWIPI_communication_t
 * \brief Communication mode
 *
 * CWIPI_communication_t gives the different communication mode 
 */

typedef enum {

  CWIPI_COMM_PAR_WITH_PART,    /*!< Parallel communcation 
                                    on partitioned source */
  CWIPI_COMM_PAR_WITHOUT_PART, /*!< Parallel communcation 
                                    on unpartitioned source defined on 
                                    all process */
  CWIPI_COMM_SEQ,              /*!< Parallel communcation 
                                    on unpartitioned source defined on 
                                    master process */
  CWIPI_COMM_INTERNAL,         /*!< Internal communcation within a process */

} CWIPI_communication_t;

/**
 * \enum CWIPI_frequency_t
 * \brief  Echange frequency 
 *
 * CWIPI_frequency_t describes the different ways to define coupling data 
 */

typedef enum {

  CWIPI_FREQUENCY_NO,                /*!< Exchange unrelated to time */
  CWIPI_FREQUENCY_EACH_TIME_STEP,    /*!< Exchange at each time step */
  CWIPI_FREQUENCY_RELATED_TIME_STEP, /*!< Exchange frequency is linked to 
                                          time step  */ 
  CWIPI_FREQUENCY_ASYNCHRONOUS,      /*!< Exchange are asynchronous */ 

} CWIPI_frequency_t;

/**
 * \enum CWIPI_src_moving_t
 * \brief Active moving source
 *
 * CWIPI_src_moving_t active moving source (mesh or point cloud)  
 */

typedef enum {

  CWIPI_SRC_MOVING_ON,     /*!< Moving on */ 
  CWIPI_SRC_MOVING_OFF,    /*!< Moving off */ 

} CWIPI_source_moving_t;

/**
 * \enum CWIPI_field_nature_t
 * \brief Field nature
 *
 * CWIPI_field_nature_t gives dfieferent nature 
 */

typedef enum {

  CWIPI_FIELD_NATURE_CELL_CENTER,  /*!< Cell center field */
  CWIPI_FIELD_NATURE_CELL_VERTEX,  /*!< Cell vertex field */
  CWIPI_FIELD_NATURE_USER,         /*!< User defined field */ 

} CWIPI_field_nature_t ;

/**
 * \enum CWIPI_field_type_t
 * \brief Field type
 *
 * CWIPI_field_type_t gives types accepted by field
 */

typedef enum { 

  CWIPI_FIELD_TYPE_DOUBLE, /*!< Field type is double */

} CWIPI_field_type_t;

/**
 * \enum CWIPI_interpolation_t
 * \brief Interpolation type
 *
 * CWIPI_interpolation_t gives the different ways to interpolate
 */

typedef enum {

  CWIPI_INTERPOLATION_DEFAULT,  /*!< Default interpolation */
  CWIPI_INTERPOLATION_USER,     /*!< User interpolation */

} CWIPI_interpolation_t;

/**
 * \enum CWIPI_status_t
 * \brief Error codes
 *
 * CWIPI_status_t defines the different error codes 
 */

typedef enum {

  CWIPI_STATUS_OK,            /*!< Output without error */
  CWIPI_STATUS_ERROR,         /*!< output with error */

} CWIPI_status_t;

/**
 * \enum CWIPI_block_t
 * \brief Elements taken into account
 *
 * CWIPI_block_t defines elements taken into account  
 */

typedef enum {

  CWIPI_BLOCK_NODE,          /*!< Node */
  CWIPI_BLOCK_EDGE2,         /*!< Edge with two nodes */
  CWIPI_BLOCK_FACE_TRIA3,    /*!< Triangle with three nodes */
  CWIPI_BLOCK_FACE_QUAD4,    /*!< Quadrangle with three nodes */
  CWIPI_BLOCK_FACE_POLY,     /*!< Generic polygon */
  CWIPI_BLOCK_CELL_TETRA4,   /*!< Tetrahedron with four nodes */
  CWIPI_BLOCK_CELL_HEXA8,    /*!< Hexahedron with eight nodes */
  CWIPI_BLOCK_CELL_PRISM6,   /*!< Prism with six nodes */
  CWIPI_BLOCK_CELL_PYRAM5,   /*!< Pyramid with five nodes */
  CWIPI_BLOCK_CELL_POLY,     /*!< Generic polyhedron */

} CWIPI_block_t;

/**
 * \enum CWIPI_support_t
 * \brief Geomtric supports
 *
 * CWIPI_support_t gives different geometric supports on which source fields 
 * are defined  
 */

typedef enum {

  CWIPI_SUPPORT_MESH,         /*!< Mesh */
  CWIPI_SUPPORT_POINT_CLOUD,  /*!< Point cloud */

} CWIPI_support_t;

/**
 * \enum CWIPI_geometry_t
 * \brief Geomtric algorithms
 *
 * CWIPI_geometry_t gives different geometric algorithm on which interpolation 
 * method is based 
 */

typedef enum {

  CWIPI_GEOMETRY_CLOSEST_POINT, /*!< Closest points */
  CWIPI_GEOMETRY_INTERSECTION,  /*!< Meshes intersection */
  CWIPI_GEOMETRY_LOCATION,      /*!< Location into a mesh */

} CWIPI_geometry_t;

/**
 * \enum CWIPI_interface_t
 * \brief Coupling interfaces
 *
 * CWIPI_interface_t gives different coupling interfaces 
 */

typedef enum {

  CWIPI_INTERFACE_POINT,    /*!< Point interface */ 
  CWIPI_INTERFACE_LINEAR,   /*!< Linear interface */ 
  CWIPI_INTERFACE_SURFACE,  /*!< Surface interface */ 
  CWIPI_INTERFACE_VOLUME,   /*!< Volume interface */ 

} CWIPI_interface_t;

/*============================================================================
 * User interpolation type
 *============================================================================*/

/**
 * \typedef void (*CWIPI_interp_from_location_t)
 * \brief User interpolation function from location into a mesh.
 *
 * void (*CWIPI_interp_from_location_t) defines the user interpolation 
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
 * \param [in]  src_field_nature            source field nature
 * \param [in]  src_field                   source field
 *                                          (size depends on field type and stride)
 * \param [in]  src_field_nature            target field nature
 * \param [out] tgt_field                   target field
 *                                          (size = stride * n_tgt_pts)
 */

typedef void (*CWIPI_interp_from_location_t)
  (const int                   interface_type,
   const int                   n_src_vtcs,
   const int                   n_src_std_elts,
   const int                   n_src_poly,
   const int                   n_tgt_pts,
   const double                src_vtcs_coords[],
   const cwipi_long_t          src_parent_elts_num[],
   const cwipi_long_t          src_parent_vtcs_num[],
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
   const CWIPI_field_nature_t  src_field_nature,
   const void                 *src_field,
   const CWIPI_field_nature_t  tgt_field_nature,
   void                       *tgt_field
   );

/**
 * \typedef void (*CWIPI_interp_from_intersec_t)
 * \brief User interpolation function from intersection between meshes <b>(Not implemented yet)</b>
 *
 * void (*CWIPI_interp_from_intersec_t) defines the user interpolation 
 * interface to take into account an user interpolation from intersection 
 * between source and target meshes 
 *
 * \param [in]  interface_type              Interface type
 *
 */

typedef void (*CWIPI_interp_from_intersec_t)
  (const int interface_type
   );

/**
 * \typedef void (*CWIPI_interp_from_closest_pts_t)
 * \brief User interpolation function from closest points <b>(Not implemented yet)</b>
 *
 * void (*CWIPI_interp_from_closest_pts_t) defines the user interpolation 
 * interface to take into account an user interpolation from <i>n</i> closest
 * point 
 *
 * \param [in]  interface_type              Interface type
 *
 */

typedef void (*CWIPI_interp_from_closest_pts_t)
  (const int interface_type
   );

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes 
 *============================================================================*/

/*=============================================================================
 * Public function prototypes - with xml data
 *============================================================================*/

/*!
 * \brief Initialize CWIPI from xml data file
 *
 * This function create the MPI intra communicator for this code from
 * the MPI inter communicator that contains all code process. It is a
 * synchronization point between all codes.
 *
 * \param [in]  inter_comm   MPI inter communicator
 * \param [in]  data_file    xml_data_file
 * \param [out] intra_comm   MPI intra communicator
 *
 */

void 
CWIPI_init_from_xml
(const MPI_Comm           inter_comm,
 const char              *data_file,
 MPI_Comm                *intra_comm);

/*=============================================================================
 * Public function prototypes - without xml data
 *============================================================================*/

/*!
 * \brief Initialize CWIPI.
 *
 * This function create the MPI intra communicator for this code from
 * the MPI inter communicator that contains all code process. It is a
 * synchronization point between all codes
 *
 * \param [in]  inter_comm   MPI inter communicator
 * \param [in]  code_name    Name of this code
 * \param [in]  time_init    Time init
 * \param [out] intra_comm   MPI intra communicator
 *
 */

void 
CWIPI_init
(const MPI_Comm           inter_comm,
 const char              *code_name,
 const double            *time_init,
 MPI_Comm                *intra_comm);

/*!
 * \brief Writing output to fortran file.
 *
 * This function set the file fortran logical unit for writing output.
 *
 * \param [in]  iunit        File fortan logical unit
 *
 */

void 
PROCF (cwipi_output_to_fortran_unit_set, CWIPI_OUTPUT_TO_FORTRAN_UNIT_SET)
(int *iunit);

/*!
 * \brief Writing output to file.
 *
 * This function set the file for writing output.
 *
 * \param [in] output_file    Output file
 *
 */

void 
CWIPI_output_to_file_set
(FILE *output_file);

/*!
 * \brief Dump application properties.
 *
 * This function dump application properties.
 *
 */

void 
CWIPI_dump_appli_properties
(void);

/*!
 * \brief Creating a coupling object.
 *
 * This function dump application properties.
 *
 * \param [in]  coupling_name           <-- Coupling identifier
 * \param [in]  coupling_type           <-- Coupling type
 * \param [in]  coupled_application     <-- Coupled application name
 * \param [in]  entitiesDim             <-- Mesh entities dimension (1, 2 or 3)
 * \param [in]  tolerance               <-- Geometric tolerance to locate
 * \param [in]  mesh_type               <-- CWIPI_STATIC_MESH
 * \param [in]                              CWIPI_MOBILE_MESH (not implemented yet)
 * \param [in]  solver_type             <-- CWIPI_SOLVER_CELL_CENTER
 * \param [in]                              CWIPI_SOLVER_CELL_VERTEX
 * \param [in]  output_frequency        <-- Output frequency
 * \param [in]  output_format           <-- Output format to visualize exchanged fields
 *                                          on the coupled mesh. Choice between :
 *                                            - "EnSight Gold"
 *                                            - "MED_ficher"
 *                                            - "CGNS"
 * \param [in]  output_format_option    <-- Output options "opt1, opt2,
 *                             text             output text files
 *                             binary              output binary files (default)
 *                             big_endian          force binary files
 *                                                 to big-endian
 *                             discard_polygons    do not output polygons
 *                                                 or related values
 *                             discard_polyhedra   do not output polyhedra
 *                                                 or related values
 *                             divide_polygons     tesselate polygons
 *                                                 with triangles
 *                             divide_polyhedra    tesselate polyhedra
 *                                                 with tetrahedra and pyramids
 *                                                 (adding a vertex near
 *                                                 each polyhedron's center)
 *
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_create_coupling
(const char                 *id,
 const CWIPI_communication_t communication,
 const char                 *coupled_application,
 const CWIPI_geometry_t      geometry,
 const CWIPI_source_type_t   source,
 const double                coupling_time_step);

void 
CWIPI_time_step_update
(const char                 *id,
 const double                time_step);

void 
CWIPI_time_update
(const char                 *id);

void 
CWIPI_location_properties_def
(const char                 *id,
 const double                tolerance);

void 
CWIPI_closest_point_properties_def
(const char                 *id,
 const double                tolerance,
 const int                   n_closest_point);

void 
CWIPI_intersection_properties_def
(const char                 *id,
 const double                tolerance);

void 
CWIPI_visualization_def
(const char                 *id,
 const int                   output_frequency,
 const char                 *format,
 const char                 *format_option);

/*----------------------------------------------------------------------------
 *
 * Set points to locate. This function must be called if the points to locate
 * do not correspond to :
 *        - vertices for CELL_VERTEX nature
 *        - cell center for CELL_CENTER nature
 *
 * parameters:
 *   coupling_id        <-- coupling identifier
 *   n_points           <-- number of points to locate
 *   coordinates        <-- coordinates of points to locate (enterlaced)
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_user_target_points_def
(const char                 *coupling_id,
 const int                   n_points,
 double                      coordinate[]);

/*----------------------------------------------------------------------------
 *
 * Define the support mesh for a coupling. The connectivity is sorted if
 * necessary.
 *
 *
 * Order definition :
 *    1D : edges
 *    2D : triangles, quadrangles, polygons
 *    3D : tetrahedra, pyramids, prism, hexaedra
 *
 * Local connectivity for the following element type :
 *
 *  - edge :
 *
 *   1 x-------x 2
 *
 *  - triangle :
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - quadrangle :
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - tetrahedra :
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - pyramid :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - prism :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  -  hexaedra :
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * Order definition :
 *    1D : edges
 *    2D : triangles, quadrangles, polygons
 *    3D : tetrahedra, pyramids, prism, hexaedra
 *
 *
 * parameters:
 *   coupling_id        <-- coupling name
 *   n_vertex           <-- number of vertices
 *   n_elements         <-- number of elements
 *   coordinates        <-- vertex interlaced coordinates
 *   connectivity_index <-> element -> vertices index (O to n-1)
 *                          size: n_elements + 1
 *                          (out : ordered connectivity_index)
 *   connectivity       <-> element -> vertex connectivity (1 to n)
 *                          size: connectivity_index[n_elements]
 *                          (out : ordered connectivity)
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_define_mesh
(const char *coupling_id,
 const int n_vertex,
 const int n_element,
 double coordinates[],
 int connectivity_index[],
 int connectivity[]);

void 
CWIPI_shared_fvmc_nodal
(const char *coupling_name,
 void * fvmc_nodal);


/*----------------------------------------------------------------------------
 *
 * Add polyhedra to the mesh
 *
 * parameters:
 *   coupling_id                  <-- Coupling identifier
 *   n_elements                   <-- Polyhedra number to add
 *   face_index                   <-- Face index (0 to n-1)
 *                                    size : n_elements + 1
 *   cell_to_face_connectivity    <-- Polyhedra -> face (1 to n)
 *                                    size : face_index[n_elements]
 *   n_faces                      <-- Faces number
 *   face_connectivity_index      <-- Face connectivity index (0 to n-1)
 *                                    size : n_faces + 1
 *   face_connectivity            <-- Face connectivity (1 to n)
 *                                    size : face_connectivity_index[n_faces]
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_add_polyhedra
(const char *coupling_id,
 const int n_element,
 int face_index[],
 int cell_to_face_connectivity[],
 const int n_faces,
 int face_connectivity_index[],
 int face_connectivity[]);

/*----------------------------------------------------------------------------
 *
 * Distant points location
 * (synchronization point with the coupled application)
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void 
CWIPI_geometry_compute
(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Update location
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void 
CWIPI_source_delete
(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Exchanges from XML data
 *
 *----------------------------------------------------------------------------*/

void
CWIPI_exchange(void)

/*----------------------------------------------------------------------------
 *
 * Exchange data with the coupled application.
 * It is a synchronization point with the coupled application
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   exchange_name        <-- Exchange name
 *   exchange_type        <-- Exchange type (not implemented yet)
 *   stride               <-- Number of interlaced fields
 *   time_step            <-- Time step  (only for visualization)
 *   time_value           <-- Time value (only for visualization)
 *   sending_field_name   <-- Sending field name
 *   sending_field        <-- Sending field (NULL -> no sending)
 *   receiving_field_name <-- Receiving field name
 *   receiving_field      --> Receiving field
 *   n_not_located_points --> Number of not located points
 *
 * returns :
 *   CWIPI_exchange_status
 *
 *----------------------------------------------------------------------------*/

CWIPI_exchange_status_t 
CWIPI_sendrecv
(const char                *coupling_id,
 const char                *exchange_id,
 const int                  stride,
 const char                *sending_field_id,
 const char                *receiving_field_id,
 int                       *n_uncomputed_target);

/*----------------------------------------------------------------------------
 *
 * Send interpolated data to the coupled application. 
 * Non blocking comunication.
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   exchange_name        <-- Exchange name
 *   stride               <-- Number of interlaced field
 *   time_step            <-- Time step  (only for visualization)
 *   time_value           <-- Time value (only for visualization)
 *   sending_field_name   <-- Sending field name
 *   sending_field        <-- Sending field 
 *   request              --> Request
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_issend
(const char                *coupling_name,
 const char                *exchange_name,
 const int                  stride,
 const char                *sending_field_id,
 int                       *request);

/*----------------------------------------------------------------------------
 *
 * Receive interpolated data from the coupled application. 
 * Non blocking comunication. receiving_field is fully updated after 
 * CWIPI_wait_irecv calling
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   exchange_name        <-- Exchange name
 *   stride               <-- Number of interlaced field
 *   time_step            <-- Time step  (only for visualization)
 *   time_value           <-- Time value (only for visualization)
 *   receiving_field_name <-- Receiving field name
 *   receiving_field      <-- Receiving field 
 *   request              --> Request
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_irecv
(const char                *coupling_name,
 const char                *exchange_name,
 const int                  stride,
 const char                *receiving_field_id,
 int                       *request);

/*----------------------------------------------------------------------------
 *
 * Wait for CWIPI_issend. 
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   request              <-- Request
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_wait_issend
(const char  *coupling_name,
 int          request);

/*----------------------------------------------------------------------------
 *
 * Wait for CWIPI_irecv. 
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   request              <-- Request
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_wait_irecv
(const char  *coupling_name,
 int          request);

/*----------------------------------------------------------------------------
 *
 * Define a user interpolation from location
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   fct                  <-- Interpolation function
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_interp_from_location_set
(const char *coupling_id,
 CWIPI_interp_from_location_t fct);

/*----------------------------------------------------------------------------
 *
 * Define a FORTRAN user interpolation from mesh location
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   fct                  <-- Interpolation function
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_interp_from_location_set_f
(const char *coupling_id,
 void* fct);

/*----------------------------------------------------------------------------
 *
 * Define a user interpolation from mesh intersection
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   fct                  <-- Interpolation function
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_interp_from_intersect_set
(const char *coupling_id,
 CWIPI_interp_from_intersec_t fct);

/*----------------------------------------------------------------------------
 *
 * Define a FORTRAN user interpolation from mesh intersection
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   fct                  <-- Interpolation function
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_interp_from_intersect_set_f
(const char *coupling_id,
 void* fct);


/*----------------------------------------------------------------------------
 *
 * Define a user interpolation from closest points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   fct                  <-- Interpolation function
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_interp_from_closest_pts_set
(const char *coupling_id,
 CWIPI_inter_from_closest_pts_t fct);

/*----------------------------------------------------------------------------
 *
 * Define a FORTRAN user interpolation from closest points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   fct                  <-- Interpolation function
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_interp_from_closest_pts_set_f
(const char *coupling_id,
 void* fct);

/*----------------------------------------------------------------------------
 *
 * Delete a coupling
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_coupling_delete(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Finalize cwipi. This is a synchronization point between all applications
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_finalize(void);

/*----------------------------------------------------------------------------
 *
 * Get number of uncomputed target
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> Number of uncomputed targets
 *
 *----------------------------------------------------------------------------*/

int 
CWIPI_n_uncomputed_targets_get(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get uncomputed targets
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> Uncomputed targets
 *----------------------------------------------------------------------------*/

const int *
CWIPI_uncomputed_targets_get(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get number of not located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> Number of not located points
 *
 *----------------------------------------------------------------------------*/

int 
CWIPI_n_computed_targets_get(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get computed targets
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> Computed targets
 *----------------------------------------------------------------------------*/

const int *
CWIPI_computed_targets_get(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Add a integer control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_loc_int_ctrl_param_add(const char *name, int initial_value);

/*----------------------------------------------------------------------------
 *
 * Add a double control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_loc_dbl_ctrl_param_add(const char *name, double initial_value);

/*----------------------------------------------------------------------------
 *
 * Add a string control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_loc_str_ctrl_param_add(const char *name, const char *initial_value);

/*----------------------------------------------------------------------------
 *
 * Set a integer control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_loc_int_ctrl_param_set(const char *name, int value);

/*----------------------------------------------------------------------------
 *
 * Set a double control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_loc_dbl_ctrl_param_set(const char *name, double value);

/*----------------------------------------------------------------------------
 *
 * Set a string control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_loc_str_ctrl_param_set(const char *name, const char *value);

/*----------------------------------------------------------------------------
 *
 * Get a integer control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

int 
CWIPI_loc_int_ctrl_param_get(const char *name);

/*----------------------------------------------------------------------------
 *
 * Get a double control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

double 
CWIPI_loc_dbl_ctrl_param_get(const char *name);

/*----------------------------------------------------------------------------
 *
 * Get a string control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

const char* 
CWIPI_loc_str_ctrl_param_get(const char *name);

/*----------------------------------------------------------------------------
 *
 * Delete a current application int parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_loc_int_ctrl_param_del(const char *name);

/*----------------------------------------------------------------------------
 *
 * Delete a current application double parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_loc_dbl_ctrl_param_del(const char *name);

/*----------------------------------------------------------------------------
 *
 * Delete a current application string parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_loc_str_ctrl_param_del(const char *name);

/*----------------------------------------------------------------------------
 *
 * Get a integer control parameter of a other application
 *
 * parameters
 *    application_name       <-- application name
 *    name                   <-- parameter name
 *
 *----------------------------------------------------------------------------*/

int 
CWIPI_dis_int_ctrl_param_get
(const char *application_name,
 const char *name);

/*----------------------------------------------------------------------------
 *
 * Get a double control parameter of a other application
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

double 
CWIPI_dis_dbl_ctrl_param_get
(const char *application_name,
 const char *name);

/*----------------------------------------------------------------------------
 *
 * Get a string control parameter of a other application
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

const char* 
CWIPI_dis_str_ctrl_param_get
(const char *application_name,
 const char *name);

/*----------------------------------------------------------------------------
 *
 * Synchronize local control parameters with an other application.
 *  It is a synchronization point between two applications
 *
 * parameters
 *    application_name    <-- application name
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_ctrl_param_synch(const char *application_name);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CWIPI_H__ */
