#ifndef __CWIPI_H__
#define __CWIPI_H__
/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011  ONERA

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

//#include <fvmc_nodal.h> // 

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
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Data
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_DATA_FROM_XML,
  CWIPI_DATA_IN_CODE_SOURCE,

} CWIPI_data_t;

/*----------------------------------------------------------------------------
 * Communication
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_COMM_PAR_WITH_PART,
  CWIPI_COMM_PAR_WITHOUT_PART,
  CWIPI_COMM_SEQ,
  CWIPI_COMM_INTERNAL,

} CWIPI_communication_t;

/*----------------------------------------------------------------------------
 * Moving mesh 
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_MESH_MOVING_ON,
  CWIPI_MESH_MOVING_OFF,

} CWIPI_mesh_moving_t;

/*----------------------------------------------------------------------------
 * Solver type
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_FIELD_NATURE_CELL_CENTER,
  CWIPI_FIELD_NATURE_CELL_VERTEX,
  CWIPI_FIELD_NATURE_USER,

} CWIPI_field_nature_t;

/*----------------------------------------------------------------------------
 * Field type
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_FIELD_TYPE_FLOAT,
  CWIPI_FIELD_TYPE_DOUBLE,

} CWIPI_field_type_t;

/*----------------------------------------------------------------------------
 * Coupling interpolation type
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_INTERPOLATION_DEFAULT,
  CWIPI_INTERPOLATION_USER,

} CWIPI_interpolation_t;

/*----------------------------------------------------------------------------
 * Coupling exchange status
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_STATUS_OK,
  CWIPI_STATUS_ERROR,

} CWIPI_status_t;

/*----------------------------------------------------------------------------
 * Block type
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_BLOCK_NODE,
  CWIPI_BLOCK_EDGE2,
  CWIPI_BLOCK_FACE_TRIA3,
  CWIPI_BLOCK_FACE_TRIA6,
  CWIPI_BLOCK_FACE_QUAD4,
  CWIPI_BLOCK_FACE_POLY,
  CWIPI_BLOCK_CELL_TETRA4,
  CWIPI_BLOCK_CELL_HEXA8,
  CWIPI_BLOCK_CELL_PRISM6,
  CWIPI_BLOCK_CELL_PYRAM5,
  CWIPI_BLOCK_CELL_POLY,

} CWIPI_block_t;

/*----------------------------------------------------------------------------
 * Source
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_SOURCE_MESH,
  CWIPI_SOURCE_POINT_CLOUD,

} CWIPI_source_t;

/*----------------------------------------------------------------------------
 * Geometry
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_GEOMETRY_CLOSEST_POINT,
  CWIPI_GEOMETRY_INTERSECTION,
  CWIPI_GEOMETRY_LOCATION,

} CWIPI_geometry_t;

/*----------------------------------------------------------------------------
 * Entities
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_ENTITIES_POINT,
  CWIPI_ENTITIES_FACE,
  CWIPI_ENTITIES_CELL,

} CWIPI_entities_t;

/*----------------------------------------------------------------------------
 * Function pointer to define an user interpolation from mesh location 
 *
 * parameters:
 * ----------
 *
 * entities_dim                              <-- entities dimension of
 *                                               the local mesh (1, 2 or 3)
 * n_local_vertex                            <-- local mesh vertices number
 * n_local_element                           <-- local mesh elements number
 *                                               (without polyhedra)
 * n_local_polyhedra                         <-- local mesh elements number
 * n_distant_point                           <-- located distant point number
 * local_coordinates                         <-- local mesh vertex coordinates
 * local_parent_elt_num                      <-- pointer to parent element
 *                                               (or NULL if sorted elements)
 * local_connectivity_index                  <-- element -> vertices index
 *                                               (O to n-1)
 *                                               size:n_local_elements + 1
 * local_connectivity                        <-- element -> vertex connectivity
 *                                                       of the local mesh
 *                               size:local_connectivity_index[n_local_elements]
 * local_polyhedra_face_index                <-- polyhedra volume -> faces index
 * local_polyhedra_cell_to_face_connectivity <-- polyhedra -> face connectivity
 * local_polyhedra_face_connectivity_index   <-- polyhedra faces
 *                                               face -> vertices index
 * local_polyhedra_face_connectivity         <-- polyhedra
 *                                               face -> vertex connectivity
 * distant_points_coordinates                <-- distant point coordinates
 * distant_points_location                   <-- distant point location
 * distant_points_barycentric_coordinates_index
 *                                           <-- element -> barycentric coordinates
 *                                                (0 to n-1)
 *                                               size: n_distant_point + 1
 * distant_points_barycentric_coordinates    <-- distant point barycentric coordinates
 *                                             size: distant_points_barycentric_coordinates_index[n_distant_point]
 * stride                                    <-- interlaced field number
 * local_field                               <-- local field
 * distant_field                             --> distant field
 *
 *----------------------------------------------------------------------------*/

typedef void (*CWIPI_interp_from_location_t)
  (const int entities_dim,
   const int n_local_vertex,
   const int n_local_element,
   const int n_local_polhyedra,
   const int n_distant_point,
   const double local_coordinates[],
   const int local_connectivity_index[],
   const int local_connectivity[],
   const int local_polyhedra_face_index[],
   const int local_polyhedra_cell_to_face_connectivity[],
   const int local_polyhedra_face_connectivity_index[],
   const int local_polyhedra_face_connectivity[],
   const double distant_points_coordinates[],
   const int distant_points_location[],
   const float distant_points_distance[],
   const int distant_points_barycentric_coordinates_index[],
   const double distant_points_barycentric_coordinates[],
   const int stride,
   const CWIPI_solver_type_t  solver_type,
   const void *local_field,
   void *distant_field
   );

/*----------------------------------------------------------------------------
 * Function pointer to define an user interpolation from mesh intersection 
 *
 *                Not implemented yet
 *
 * parameters:
 * ----------
 *
 * entities_dim                              <-- entities dimension of
 *                                               the local mesh (1, 2 or 3)
 * n_source_vertex                            <-- local mesh vertices number
 * n_source_element                           <-- local mesh elements number
 *                                               (without polyhedra)
  *----------------------------------------------------------------------------*/

typedef void (*CWIPI_interp_from_intersec_t)
  (const int entities_dim,
   const int n_source_vertex,
   const int n_source_element
   );

/*----------------------------------------------------------------------------
 * Function pointer to define an user interpolation from closest points 
 *
 *                Not implemented yet
 *
 * parameters:
 * ----------
 *
 * entities_dim                              <-- entities dimension of
 *                                               the local mesh (1, 2 or 3)
 * n_source_vertex                            <-- local mesh vertices number
 * n_source_element                           <-- local mesh elements number
 *                                               (without polyhedra)
  *----------------------------------------------------------------------------*/

typedef void (*CWIPI_interp_from_closest_pts_t)
  (const int n_source_vertex,
   const int n_target_vertex);

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes 
 *============================================================================*/

/*!
 * \brief Initialize CWIPI.
 *
 * This function create the MPI intra communicator for this code from
 * the MPI inter communicator that contains all code process. It is a
 * synchronization point between all codes
 *
 * \param [in]  inter_comm   MPI inter communicator
 * \param [in]  code_name    name of this code
 * \param [in]  data_type    data type (read in xml data file or not)
 * \param [out] intra_comm   MPI intra communicator
 *
 */

void 
CWIPI_init
(const MPI_Comm           inter_comm,
 const char              *code_name,
 const CWIPI_data_t       data_type,
 const char              *data_file,
 MPI_Comm                *intra_comm);

/*=============================================================================
 * Public function prototypes - with xml data
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Initialize the cwipi library.
 * Redirect outputs in a file (Standard output with output_listing = NULL or
 * output_logical_unit = -1)
 * Create the current communicator application from 'common_comm'.
 *
 * parameters:
 *   common_comm       <-- Common MPI communicator
 *   data_from         <-- Where data are defined
 *   application_comm  --> Internal MPI communicator for the current
 *                         application
 *
 *         Synchronization point between all applications
 *----------------------------------------------------------------------------*/

void 
CWIPI_code_init
(const double            *time_init,
 const double            *time_step);
 MPI_Comm                *application_comm);

/*=============================================================================
 * Public function prototypes - without xml data
 *============================================================================*/


/*----------------------------------------------------------------------------
 *
 * Initialize the cwipi library.
 * Redirect outputs in a file (Standard output with output_listing = NULL or
 * output_logical_unit = -1)
 * Create the current communicator application from 'common_comm'.
 *
 * parameters:
 *   common_comm       <-- Common MPI communicator
 *   data_from         <-- Where data are defined
 *   application_comm  --> Internal MPI communicator for the current
 *                         application
 *
 *         Synchronization point between all applications
 *----------------------------------------------------------------------------*/

void 
CWIPI_code_init
(const MPI_Comm           common_comm,
 const char              *application_name,
 const CWIPI_data_t       data_from,
 const char              *data_file,
 const double            *time_init,
 const double            *time_step,
 MPI_Comm                *application_comm);

/*----------------------------------------------------------------------------
 *
 * Initialize the cwipi library.
 * Redirect outputs in a file (Standard output with output_listing = NULL or
 * output_logical_unit = -1)
 * Create the current communicator application from 'common_comm'.
 *
 * parameters:
 *   common_comm       <-- Common MPI communicator
 *   data_from         <-- Where data are defined
 *   application_comm  --> Internal MPI communicator for the current
 *                         application
 *
 *         Synchronization point between all applications
 *----------------------------------------------------------------------------*/

void 
CWIPI_code_init
(const MPI_Comm           common_comm,
 const char              *application_name,
 const CWIPI_data_t       data_from,
 const char              *data_file,
 const double            *time_init,
 const double            *time_step,
 MPI_Comm                *application_comm);

/*----------------------------------------------------------------------------
 *
 * Set up the file used for the output listing
 *
 * parameters:
 *   output_listing      <-- Output listing file (C function)
 *----------------------------------------------------------------------------*/

void 
PROCF (cwipi_set_output_logical_unit, CWIPI_SET_OUTPUT_LOGICAL_UNIT) 
(int *iunit);

void 
CWIPI_set_output_listing
(FILE *output_listing);

/*----------------------------------------------------------------------------
 *
 * Dump application properties
 *
 *----------------------------------------------------------------------------*/

void 
CWIPI_dump_appli_properties
(void);

/*----------------------------------------------------------------------------
 *
 * Create a coupling object
 *
 * parameters:
 *   coupling_name           <-- Coupling identifier
 *   coupling_type           <-- Coupling type
 *   coupled_application     <-- Coupled application name
 *   entitiesDim             <-- Mesh entities dimension (1, 2 or 3)
 *   tolerance               <-- Geometric tolerance to locate
 *   mesh_type               <-- CWIPI_STATIC_MESH
 *                               CWIPI_MOBILE_MESH (not implemented yet)
 *   solver_type             <-- CWIPI_SOLVER_CELL_CENTER
 *                               CWIPI_SOLVER_CELL_VERTEX
 *   output_frequency        <-- Output frequency
 *   output_format           <-- Output format to visualize exchanged fields
 *                               on the coupled mesh. Choice between :
 *                                 - "EnSight Gold"
 *                                 - "MED_ficher"
 *                                 - "CGNS"
 *   output_format_option    <-- Output options "opt1, opt2,
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
