#ifndef __COUPLINGS_H__
#define __COUPLINGS_H__

#include <stdio.h>
//Bug mpich2
#define MPICH_IGNORE_CXX_SEEK 1
#include <mpi.h>

/*=============================================================================
 * Macro definitions
 *============================================================================*/

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
 * Mesh type
 *----------------------------------------------------------------------------*/

typedef enum {
  
  COUPLINGS_STATIC_MESH,
  COUPLINGS_MOBILE_MESH, 
  
} couplings_mesh_type_t;

/*----------------------------------------------------------------------------
 * Solver type
 *----------------------------------------------------------------------------*/

typedef enum {
  
  COUPLINGS_SOLVER_CELL_CENTER, 
  COUPLINGS_SOLVER_CELL_VERTEX,
  
} couplings_solver_type_t;


/*----------------------------------------------------------------------------
 * Coupling dimension
 *----------------------------------------------------------------------------*/

typedef enum {
  
  COUPLINGS_FIELD_DIMENSION_SCALAR,
  COUPLINGS_FIELD_DIMENSION_INTERLACED_VECTOR,
  
} couplings_field_dimension_t;

/*----------------------------------------------------------------------------
 * Coupling type
 *----------------------------------------------------------------------------*/

typedef enum {
  
  COUPLINGS_FIELD_TYPE_FLOAT,
  COUPLINGS_FIELD_TYPE_DOUBLE,
  
} couplings_field_type_t;

/*----------------------------------------------------------------------------
 * Coupling interpolation type
 *----------------------------------------------------------------------------*/

typedef enum {
  
  COUPLINGS_INTERPOLATION_DEFAULT,
  COUPLINGS_INTERPOLATION_USER,
  
} couplings_interpolation_t;

/*----------------------------------------------------------------------------
 * Coupling non located point treatment
 *----------------------------------------------------------------------------*/

typedef enum {
  
  COUPLINGS_NOT_LOCATED_POINT_TREATMENT_STANDARD, 
  COUPLINGS_NOT_LOCATED_POINT_TREATNENT_USER,           
  
} couplings_not_located_point_treatment_t;

/*----------------------------------------------------------------------------
 * Coupling exchange status
 *----------------------------------------------------------------------------*/

typedef enum {
  
  COUPLINGS_EXCHANGE_OK,
  COUPLINGS_EXCHANGE_BAD_RECEIVING,           
  
} couplings_exchange_status_t;

/*----------------------------------------------------------------------------
 * Function pointer to define an user interpolation method (callback)
 *
 * parameters:
 *   dim                                    <-- space dimension (1, 2 or 3)
 *   n_local_vertex                         <-- number of vertices
 *   n_local_element                        <-- number of elements
 *   n_distant_point                        <-- number of located distant point
 *   local_coordinates                      <-- vertex coordinates
 *   local_parent_vertex_num                <-- pointer to parent vertex 
 *                                              numbers (or NULL)
 *   local_connectivity_index               <-- polygon face -> vertices index 
 *                                              (O to n-1)
 *                                              size: n_local_elements + 1
 *   local_connectivity                     <-- element -> vertex connectivity
 *                                       size: local_connectivity_index[n_local_elements] 
 *   distant_points_coordinates             <-- distant point coordinates
 *   distant_points_location                <-- distant point location
 *   distant_points_barycentric_coordinates <-- distant point barycentric_coordinates
 *   field_type                             <-- field type 
 *   field_dimension                        <-- field dimension
 *   field_nature                           <-- field dimension
 *   local_field                            <-- local field
 *   distant_field                          --> distant field
 *
 *----------------------------------------------------------------------------*/

typedef void (couplings_interpolation_fct_t) 
  (const int entities_dim,
   const int n_local_vertex,
   const int n_local_element,
   const int n_local_polhyedra,
   const int n_distant_point,
   const double local_coordinates[],
   const int *local_parent_vertex_num,
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
   const couplings_field_dimension_t data_dimension,
   const couplings_solver_type_t  solver_type,
   const void *local_field,
   void *distant_field
   );

/*----------------------------------------------------------------------------
 * Function pointer to define a user treatment for points that are not located
 * (callback)
 *
 * This function type updates the received field by taking into account 
 * the not located points
 * 
 * parameters:
 *   dim                                    <-- space dimension (1, 2 or 3)
 *   n_local_vertex                         <-- number of vertices
 *   n_local_element                        <-- number of elements
 *   n_distant_point                        <-- number of located distant point
 *   local_coordinates                      <-- vertex coordinates
 *   local_parent_vertex_num                <-- pointer to parent vertex 
 *                                              numbers (or NULL)
 *   local_connectivity_index               <-- polygon face -> vertices index 
 *                                              (O to n-1)
 *                                              size: n_local_elements + 1
 *   local_connectivity                     <-- element -> vertex connectivity
 *                                       size: local_connectivity_index[n_local_elements] 
 *   field_type                             <-- field type 
 *   field_dimension                        <-- field dimension
 *   n_points_to_locate                     <-- number of points to locate
 *   n_not_located_points                   <-- number of not located points
 *   not_located_point_list                 <-- not located points list
 *   point_to_locate_coordinates            <-- not located points coordinates
 *   received_field                         <-> received field
 *
 *----------------------------------------------------------------------------*/

typedef void (couplings_not_located_point_treatment_fct_t)
  (const int dim,
   const int n_local_vertex,
   const int n_local_element,
   const int n_distant_point,
   const double local_coordinates[],
   const int *local_parent_vertex_num,
   const int local_connectivity_index[],
   const int local_connectivity[],
   const couplings_field_type_t data_type,
   const couplings_field_dimension_t data_dimension,
   const int n_points_to_locate,
   const int n_located_points,
   const double not_located_point_list[],
   const double point_to_locate_coordinates[],
   void *received_field
); 

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Initialize the couplings library.
 * Redirect outputs in a file (Standard output with output_listing = NULL or 
 * output_logical_unit = -1)
 * Create the current communicator application from 'common_comm'.
 *
 * parameters:
 *   common_comm       <-- Common MPI communicator
 *   output_listing    <-- Output listing file 
 *   application_name  <-- Current application name 
 *   application_comm  --> Internal MPI communicator for the current 
 *                         application
 *
 * This is a synchronization between all applications 
 *----------------------------------------------------------------------------*/

void couplings_init
(const MPI_Comm common_comm,
 FILE           *output_listing, 
 const char     *application_name, 
 MPI_Comm       *application_comm);

/*----------------------------------------------------------------------------
 *
 * Add a integer control parameter
 * 
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void couplings_add_local_int_control_parameter(const char *name, int initial_value);

/*----------------------------------------------------------------------------
 *
 * Add a double control parameter
 * 
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void couplings_add_local_double_control_parameter(const char *name, double initial_value);

/*----------------------------------------------------------------------------
 *
 * Set a integer control parameter
 * 
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void couplings_set_local_int_control_parameter(const char *name, int value);

/*----------------------------------------------------------------------------
 *
 * Set a double control parameter
 * 
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void couplings_set_local_double_control_parameter(const char *name, double value);

/*----------------------------------------------------------------------------
 *
 * Get a integer control parameter of the current application
 * 
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

int couplings_get_local_int_control_parameter(const char *name);

/*----------------------------------------------------------------------------
 *
 * Get a double control parameter of the current application
 * 
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

double couplings_get_local_double_control_parameter(const char *name);

/*----------------------------------------------------------------------------
 *
 * Delete a current application int parameter
 * 
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void couplings_delete_local_int_control_parameter(const char *name);

/*----------------------------------------------------------------------------
 *
 * Delete a current application double parameter
 * 
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void couplings_delete_local_double_control_parameter(const char *name);

/*----------------------------------------------------------------------------
 *
 * Get a integer control parameter of a other application
 * 
 * parameters
 *    application_name       <-- application name
 *    name                   <-- parameter name
 *
 *----------------------------------------------------------------------------*/

int couplings_get_distant_int_control_parameter
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

double couplings_get_distant_double_control_parameter
(const char *application_name, 
 const char *name);

/*----------------------------------------------------------------------------
 *
 * Synchronise local control parameters with an other application.
 *  This is a synchornisation point with this second application
 * 
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void couplings_synchronise_control_parameter(const char *application_name);

/*----------------------------------------------------------------------------
 *
 * Dump application properties
 *
 *----------------------------------------------------------------------------*/

void couplings_dump_application_properties();

/*----------------------------------------------------------------------------
 *
 * Create a coupling object 
 *
 * parameters:
 *   coupling_name           <-- Coupling name
 *   dim                     <-- Mesh dim
 *   coupled_application     <-- Coupled application name
 *   solver_type             <-- Solver type
 *   output_format           <-- Output format to visualize exchanged fields 
 *                               on the coupled mesh. Choice between :
 *                                 - "EnSight Gold"
 *                                 - "MED_ficher"
 *                                 - "CGNS"
 *   output_format_option    <-- Outpout options
 *                             text                output text files
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

void couplings_create_coupling
( const char  *coupling_name,
  const char  *coupled_application,
  const int entitiesDim,
  const double tolerance,
  const couplings_mesh_type_t mesh_type,
  const couplings_solver_type_t solver_type, 
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
 *   coupling_id        <-- coupling identificator
 *   n_points           <-- number of points to locate
 *   coordinates        <-- coordinates of points to locate (enterlaced)
 *
 *----------------------------------------------------------------------------*/

void couplings_set_points_to_locate
(const char  *coupling_id,
 const int    n_points,
 double coordinate[]); 

/*----------------------------------------------------------------------------
 *
 * Define the support mesh for a coupling. The connectivity is ordered if 
 * necessary.
 *
 *  local connectivity for the following element type :
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
 * parameters:
 *   coupling_id        <-- coupling name
 *   ordered_element    <-- 1 : elements are ordered
 *                          0 : elements are not ordered 
 *                          (Warning : sorting is not yet implemented !)
 *                          
 *                          order definition :
 *                            1D : edges 
 *                            2D : triangles, quadrangles, polygons
 *                            3D : tetrahedra, pyramids, prism, hexaedra
 *   n_vertex           <-- number of vertex
 *   n_elements         <-- number of elements
 *   coordinates        <-- vertex enterlaced coordinates
 *   connectivity_index <-> element -> vertices index (O to n-1)
 *                          size: n_elements + 1 (out : ordered connectivity_index)
 *   connectivity       <-> element -> vertex connectivity
 *                          size: connectivity_index[n_elements] 
 *                          (out : ordered connectivity)
 *
 *----------------------------------------------------------------------------*/

void couplings_define_mesh(const char *coupling_id, 
                           const int n_vertex,
                           const int n_element,
                           const double coordinates[],
                           int connectivity_index[],
                           int connectivity[]);

void couplings_add_polyhedra(const char *coupling_id, 
                             const int n_element,
                             int face_index[],
                             int cell_to_face_connectivity[],
                             int face_connectivity_index[],
                             int face_connectivity[]);

/*----------------------------------------------------------------------------
 *
 * Exchange data with the coupled application. This is a synchronization point 
 * with the coupled application 
 *
 * parameters
 *   coupling_id          <-- Coupling identificator
 *   exchange_name        <-- Exchange name
 *   exchange_type        <-- Exchange type
 *   exchange_dimension   <-- Dimension of exchanged data : 
 *                            - COUPLINGS_DIMENSION_SCALAR
 *                            - COUPLINGS_DIMENSION_INTERLACED_VECTOR 
 *   time_step            <-- Time step  (only for visualization) 
 *   time_value           <-- Time value (only for visualization)
 *   sending_field_name   <-- Sending field name 
 *   sending_field        <-- Sending field (NULL -> no sending)
 *   receiving_field_name <-- Receiving field name 
 *   receiving_field      --> Receiving field
 *
 * returns :
 *   1 if data were received 
 *   0 else 
 *
 *----------------------------------------------------------------------------*/

couplings_exchange_status_t couplings_exchange
(const char                          *coupling_id,
 const char                          *exchange_name,
 const couplings_field_dimension_t    exchange_dimension, 
 const int                            time_step, 
 const double                         time_value,
 const char                          *sending_field_name,
 const double                        *sending_field, 
 char                                *receiving_field_name,
 double                              *receiving_field,
 int                                 *nNotLocatedPoints);

/*----------------------------------------------------------------------------
 *
 * Define the interpolation function 
 *
 * parameters
 *   coupling_id          <-- Coupling identificator
 *   fct                  <-- Interpolation function
 *
 *----------------------------------------------------------------------------*/

void couplings_set_interpolation_function
(const char *coupling_id,
 couplings_interpolation_fct_t * fct); 

/*----------------------------------------------------------------------------
 *
 * Define the non located point treatment function
 *
 * parameters
 *   coupling_id          <-- Coupling identificator
 *   fct                  <-- Non located point treatment function
 *
 *----------------------------------------------------------------------------*/

/*void couplings_set_not_located_point_treatment_function
(const char *coupling_id,
 couplings_not_located_point_treatment_fct_t *const fct); */

/*----------------------------------------------------------------------------
 *
 * Delete a coupling
 *
 * parameters
 *   coupling_id          <-- Coupling identificator
 *
 *----------------------------------------------------------------------------*/

void couplings_delete_coupling(const char *coupling_id); 

/*----------------------------------------------------------------------------
 *
 * Finalize couplings. This is a synchronization point between all applications 
 *
 *----------------------------------------------------------------------------*/

void couplings_finalize(); 

/*----------------------------------------------------------------------------
 *
 * Get not located points
 *
 * parameters
 *   coupling_id          <-- Coupling identificator
 *
 *----------------------------------------------------------------------------*/

const int * couplings_get_not_located_points(const char *coupling_id);

/*----------------------------------------------------------------------------
 *
 * Get not located points
 *
 * parameters
 *   coupling_id          <-- Coupling identificator
 *   
 * return
 *   locatedPoints        <-- Located points    
 *
 *----------------------------------------------------------------------------*/

const int * couplings_get_located_points(const char *coupling_id);


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __COUPLINGS_H__ */
