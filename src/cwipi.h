#ifndef __CWIPI_H__
#define __CWIPI_H__

//Bug mpich2
#define MPICH_IGNORE_CXX_SEEK 1

#include <stdio.h>
#include <mpi.h>

#include <fvm_nodal.h>

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
 * MPI ranks used for the coupling
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
  CWIPI_COUPLING_PARALLEL_WITHOUT_PARTITIONING,
  CWIPI_COUPLING_SEQUENTIAL,

} cwipi_coupling_type_t;

/*----------------------------------------------------------------------------
 * MPI ranks used for the coupling
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_BASIC_INFO,
  CWIPI_DISTANT_MESH_INFO,

} cwipi_located_point_info_t;

/*----------------------------------------------------------------------------
 * Mesh type
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_STATIC_MESH,
  CWIPI_MOBILE_MESH,

} cwipi_mesh_type_t;

/*----------------------------------------------------------------------------
 * Solver type
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_SOLVER_CELL_CENTER,
  CWIPI_SOLVER_CELL_VERTEX,

} cwipi_solver_type_t;

/*----------------------------------------------------------------------------
 * Coupling type
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_FIELD_TYPE_FLOAT,
  CWIPI_FIELD_TYPE_DOUBLE,

} cwipi_field_type_t;

/*----------------------------------------------------------------------------
 * Coupling interpolation type
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_INTERPOLATION_DEFAULT,
  CWIPI_INTERPOLATION_USER,

} cwipi_interpolation_t;

/*----------------------------------------------------------------------------
 * Coupling exchange status
 *----------------------------------------------------------------------------*/

typedef enum {

  CWIPI_EXCHANGE_OK,
  CWIPI_EXCHANGE_BAD_RECEIVING,

} cwipi_exchange_status_t;

/*----------------------------------------------------------------------------
 * Function pointer to define an user interpolation method (callback)
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

typedef void (*cwipi_interpolation_fct_t)
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
   const int distant_points_barycentric_coordinates_index[],
   const double distant_points_barycentric_coordinates[],
   const int stride,
   const cwipi_solver_type_t  solver_type,
   const void *local_field,
   void *distant_field
   );

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
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
 *   application_name  <-- Current application name
 *   application_comm  --> Internal MPI communicator for the current
 *                         application
 *
 * It is a synchronization point between all applications
 *----------------------------------------------------------------------------*/

void cwipi_init
(const MPI_Comm                           common_comm,
 const char                               *application_name,
 MPI_Comm                                 *application_comm);

/*----------------------------------------------------------------------------
 *
 * Set up the file used for the output listing
 *
 * parameters:
 *   output_listing      <-- Output listing file (C function)
 *----------------------------------------------------------------------------*/

void PROCF (cwipi_set_output_listing_f, CWIPI_SET_OUTPUT_LISTING_F) (int iunit);

void cwipi_set_output_listing(FILE *output_listing);

/*----------------------------------------------------------------------------
 *
 * Add a integer control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void cwipi_add_local_int_control_parameter(const char *name, int initial_value);

/*----------------------------------------------------------------------------
 *
 * Add a double control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void cwipi_add_local_double_control_parameter(const char *name, double initial_value);

/*----------------------------------------------------------------------------
 *
 * Add a string control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void cwipi_add_local_string_control_parameter(const char *name, const char *initial_value);

/*----------------------------------------------------------------------------
 *
 * Set a integer control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void cwipi_set_local_int_control_parameter(const char *name, int value);

/*----------------------------------------------------------------------------
 *
 * Set a double control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void cwipi_set_local_double_control_parameter(const char *name, double value);

/*----------------------------------------------------------------------------
 *
 * Set a string control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void cwipi_set_local_string_control_parameter(const char *name, const char *value);

/*----------------------------------------------------------------------------
 *
 * Get a integer control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

int cwipi_get_local_int_control_parameter(const char *name);

/*----------------------------------------------------------------------------
 *
 * Get a double control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

double cwipi_get_local_double_control_parameter(const char *name);

/*----------------------------------------------------------------------------
 *
 * Get a string control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

const char* cwipi_get_local_string_control_parameter(const char *name);

/*----------------------------------------------------------------------------
 *
 * Delete a current application int parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void cwipi_delete_local_int_control_parameter(const char *name);

/*----------------------------------------------------------------------------
 *
 * Delete a current application double parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void cwipi_delete_local_double_control_parameter(const char *name);

/*----------------------------------------------------------------------------
 *
 * Delete a current application string parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void cwipi_delete_local_string_control_parameter(const char *name);

/*----------------------------------------------------------------------------
 *
 * Get a integer control parameter of a other application
 *
 * parameters
 *    application_name       <-- application name
 *    name                   <-- parameter name
 *
 *----------------------------------------------------------------------------*/

int cwipi_get_distant_int_control_parameter
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

double cwipi_get_distant_double_control_parameter
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

const char* cwipi_get_distant_string_control_parameter
(const char *application_name,
 const char *name);

/*----------------------------------------------------------------------------
 *
 * Synchronize local control parameters with an other application.
 *  It is a synchronization point with this second application
 *
 * parameters
 *    application_name    <-- application name
 *
 *----------------------------------------------------------------------------*/

void cwipi_synchronize_control_parameter(const char *application_name);

/*----------------------------------------------------------------------------
 *
 * Dump application properties
 *
 *----------------------------------------------------------------------------*/

void cwipi_dump_application_properties();

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

void cwipi_create_coupling
( const char  *coupling_name,
  const cwipi_coupling_type_t coupling_type,
  const char  *coupled_application,
  const int    entitiesDim,
  const double tolerance,
  const cwipi_mesh_type_t mesh_type,
  const cwipi_solver_type_t solver_type,
  const int    output_frequency,
  const char  *output_format,
  const char  *output_format_option);

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

void cwipi_set_points_to_locate
(const char  *coupling_id,
 const int    n_points,
 double       coordinate[]);

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

void cwipi_define_mesh(const char *coupling_id,
                       const int n_vertex,
                       const int n_element,
                       double coordinates[],
                       int connectivity_index[],
                       int connectivity[]);

void cwipi_shared_fvm_nodal(const char *coupling_name,
                            fvm_nodal_t* fvm_nodal);


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
 *   face_connectivity_index      <-- Face connectivity index (0 to n-1)
 *                                    size : n_faces + 1
 *   face_connectivity            <-- Face connectivity (1 to n)
 *                                    size : face_connectivity_index[n_faces]
 *
 *----------------------------------------------------------------------------*/

void cwipi_add_polyhedra(const char *coupling_id,
                             const int n_element,
                             int face_index[],
                             int cell_to_face_connectivity[],
                             int face_connectivity_index[],
                             int face_connectivity[]);

/*----------------------------------------------------------------------------
 *
 * Location completion.
 * It is a synchronization point with the coupled application
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void cwipi_locate (const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Set coupling info
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   info                 <-- Coupling info
 *----------------------------------------------------------------------------*/

void cwipi_set_info(const char *coupling_id, const cwipi_located_point_info_t info);

/*----------------------------------------------------------------------------
 *
 * Get located points location
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 * return
 *   located points location
 *----------------------------------------------------------------------------*/

const int *cwipi_get_distant_location (const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get barycentric coordinates index
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 * return
 *   barycentric coordinates index
 *----------------------------------------------------------------------------*/

const int *cwipi_get_distant_barycentric_coordinates_index (const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get barycentric coordinates index
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 * return
 *   barycentric coordinates
 *----------------------------------------------------------------------------*/

const double *cwipi_get_distant_barycentric_coordinates (const char *coupling_id);

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
 *   cwipi_exchange_status
 *
 *----------------------------------------------------------------------------*/

cwipi_exchange_status_t cwipi_exchange
(const char                          *coupling_id,
 const char                          *exchange_name,
 const int                            stride,
 const int                            time_step,
 const double                         time_value,
 const char                          *sending_field_name,
 const double                        *sending_field,
 char                                *receiving_field_name,
 double                              *receiving_field,
 int                                 *n_not_located_points);

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

void cwipi_issend
(const char                *coupling_name,
 const char                *exchange_name,
 const int                 tag,
 const int                 stride,
 const int                 time_step,
 const double              time_value,
 const char                *sending_field_name,
 const double              *sending_field,
 int                       *request);

/*----------------------------------------------------------------------------
 *
 * Receive interpolated data from the coupled application. 
 * Non blocking comunication. receiving_field is fully updated after 
 * cwipi_wait_irecv calling
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

void cwipi_irecv
(const char                *coupling_name,
 const char                *exchange_name,
 const int                 tag,
 const int                 stride,
 const int                 time_step,
 const double              time_value,
 char                *receiving_field_name,
 double                    *receiving_field,
 int                       *request);

/*----------------------------------------------------------------------------
 *
 * Wait for cwipi_issend. 
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   request              <-- Request
 *
 *----------------------------------------------------------------------------*/

void cwipi_wait_issend(const char  *coupling_name,
                       int          request);

/*----------------------------------------------------------------------------
 *
 * Wait for cwipi_irecv. 
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   request              <-- Request
 *
 *----------------------------------------------------------------------------*/

void cwipi_wait_irecv(const char  *coupling_name,
                      int          request);

/*----------------------------------------------------------------------------
 *
 * Define the interpolation function
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   fct                  <-- Interpolation function
 *
 *----------------------------------------------------------------------------*/

void cwipi_set_interpolation_function
(const char *coupling_id,
 cwipi_interpolation_fct_t fct);

/*----------------------------------------------------------------------------
 *
 * Delete a coupling
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 *----------------------------------------------------------------------------*/

void cwipi_delete_coupling(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Finalize cwipi. This is a synchronization point between all applications
 *
 *----------------------------------------------------------------------------*/

void cwipi_finalize();

/*----------------------------------------------------------------------------
 *
 * Get not located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> Not located points
 *----------------------------------------------------------------------------*/

const int * cwipi_get_not_located_points(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get number of located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> Number of located points
 *
 *----------------------------------------------------------------------------*/

int cwipi_get_n_located_points(const char *coupling_id);

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

int cwipi_get_n_not_located_points(const char *coupling_id);


/*----------------------------------------------------------------------------
 *
 * Get number of located distant point
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> Number of located distant points
 *
 *----------------------------------------------------------------------------*/

int cwipi_get_n_located_distant_points(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get distant elements that contain located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> Number of vertices
 *
 *----------------------------------------------------------------------------*/

const int *cwipi_get_element_containing(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get number of vertices of distant elements that contain located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> Number of vertices
 *
 *----------------------------------------------------------------------------*/

const int *cwipi_get_element_containing_n_vertex(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get vertices id of distant elements that contain located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> vertices id
 *
 *----------------------------------------------------------------------------*/

const int *cwipi_get_element_containing_vertex(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get vertices coordinates of distant elements that contain located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> Vertices coordinates
 *
 *----------------------------------------------------------------------------*/

const double *cwipi_get_element_containing_vertex_coords(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get barycentric coords in distant elements for located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> Barycentric coordinates
 *
 *----------------------------------------------------------------------------*/

const double *cwipi_get_element_containing_barycentric_coordinates(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * For each located point get the MPI rank of distant element
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 * return
 *                        --> MPI ranks
 *
 *----------------------------------------------------------------------------*/

const int *cwipi_get_element_containing_MPI_rank(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Exchange Fields on vertices of element containing each located point
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   sendingField         <-- Field defined on local mesh vertices
 *   receivingField       --> Field defined on vertices of distant
 *                            elements that contain each located point
 *   stride               <-- Number of field component
 *
 *----------------------------------------------------------------------------*/

void cwipi_exchange_cell_vertex_field_of_element_containing (const char *coupling_id,
                                                             double *sendingField,
                                                             double *receivingField,
                                                             const int stride);

/*----------------------------------------------------------------------------
 *
 * Exchange field on cells that contain each located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   sendingField         <-- Field defined on local mesh vertices
 *   receivingField       --> Field defined on vertices of distant
 *                            elements that contain each located point
 *   stride               <-- Number of field component
 *
 *----------------------------------------------------------------------------*/

void cwipi_exchange_cell_center_field_of_element_containing (const char *coupling_id,
                                                             double *sendingField,
                                                             double *receivingField,
                                                             const int stride);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CWIPI_H__ */
