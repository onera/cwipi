#ifndef __COUPLINGS_CF_H__
#define __COUPLINGS_CF_H__

#include <stdio.h>
//Bug mpich2
#define MPICH_IGNORE_CXX_SEEK 1
#include <mpi.h>

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Macro for handling of different symbol names (underscored or not,
 * lowercase or uppercase) between C and Fortran, for link resolution.
 *----------------------------------------------------------------------------*/

#if !defined (__hpux)
#define PROCF(x, y) x##_
#else
#define PROCF(x, y) x
#endif

/*----------------------------------------------------------------------------
 * Macro used to handle automatic "Fortran string length" arguments
 * (not used by Code_Saturne calls, but set by many compilers).
 * Some compilers, like the Fujitsu VPP 5000 compiler, may not
 * support the variable length lists in mixed C/Fortran calls.
 *----------------------------------------------------------------------------*/

#if defined (__uxpv__)  /* Fujitsu VPP 5000 case */
#define ARGF_SUPP_CHAINE
#else
#define ARGF_SUPP_CHAINE , ...
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

void PROCF(couplings_init_cf, COUPLINGS_INIT_CF)
  (const int  *common_comm,
   const int  *output_logical_unit,
   const char *application_name,
   const int  *l_application_name,
   int        *application_comm
   ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Add a integer control parameter
 * 
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_add_local_int_control_parameter_cf, 
           COUPLINGS_ADD_LOCAL_INT_CONTROL_PARAMETER_CF)
  (const char *name, 
   const int  *l_name,
   int *initial_value
   ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Add a double control parameter
 * 
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_add_local_double_control_parameter_cf, 
           COUPLINGS_ADD_LOCAL_DOUBLE_CONTROL_PARAMETER_CF)
  (const char *name, 
   const int  *l_name,
   double *initial_value
   ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Set a integer control parameter
 * 
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_set_local_int_control_parameter_cf, 
           COUPLINGS_SET_LOCAL_INT_CONTROL_PARAMETER_CF)
  (const char *name, 
   const int  *l_name,
   int *initial_value
   ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Set a double control parameter
 * 
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_set_local_double_control_parameter_cf, 
           COUPLINGS_SET_LOCAL_DOUBLE_CONTROL_PARAMETER_CF)
  (const char *name, 
   const int  *l_name,
   double *initial_value
   ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Get a integer control parameter of the current application
 * 
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_get_local_int_control_parameter_cf, 
           COUPLINGS_GET_LOCAL_INT_CONTROL_PARAMETER_CF)
  (const char *name,
   const int  *l_name,
   int *value
   ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Get a double control parameter of the current application
 * 
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_get_local_double_control_parameter_cf,
           COUPLINGS_GET_LOCAL_DOUBLE_CONTROL_PARAMETER_CF)
  (const char *name,
   const int  *l_name,
   double *value
   ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Delete a current application parameter
 * 
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_delete_local_int_control_parameter_cf, 
           COUPLINGS_DELETE_LOCAL_INT_CONTROL_PARAMETER_CF)
  (const char *name,
   const int  *l_name
   ARGF_SUPP_CHAINE);


/*----------------------------------------------------------------------------
 *
 * Delete a current application parameter
 * 
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_delete_local_double_control_parameter_cf, 
           COUPLINGS_DELETE_LOCAL_DOUBLE_CONTROL_PARAMETER_CF)
  (const char *name,
   const int  *l_name
   ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Get a integer control parameter of a other application
 * 
 * parameters
 *    application_name       <-- application name
 *    name                   <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_get_distant_int_control_parameter_cf, 
           COUPLINGS_GET_DISTANT_INT_CONTROL_PARAMETER_CF)
  (const char *application_name,
   const int  *l_application_name,
   const char *name,
   const int  *l_name,
   int *value
   ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Get a double control parameter of a other application
 * 
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_get_distant_double_control_parameter_cf, 
           COUPLINGS_GET_DISTANT_DOUBLE_CONTROL_PARAMETER_CF)
  (const char *application_name,
   const int  *l_application_name,
   const char *name,
   const int  *l_name,
   double *value
   ARGF_SUPP_CHAINE);

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

void PROCF(couplings_synchronise_control_parameter_cf, 
           COUPLINGS_SYNCHRONISE_CONTROL_PARAMETER_CF)
  (const char *application_name,
   const int  *l_application_name
   ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Create a coupling object 
 *
 * parameters:
 *   coupled_application     <-- Coupled application name
 *   field_nature            <-- Nature of the current application fields
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
 * returns:
 *   The coupling id  
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_create_coupling_cf, 
           COUPLINGS_CREATE_COUPLING_CF)
  (const char *coupled_application,
   const int  *l_coupled_application,
   const int  *field_nature, 
   const char *output_format,
   const int  *l_output_format,
   const char *output_format_option,
   const int  *l_output_format_option,
   char  *coupling_id
   ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Set points to locate. This function must be called if the points to locate 
 * do not correspond to :
 *        - vertices for NATURE_NODE nature
 *        - cell center for NATURE_CELL_CENTER nature
 * 
 * parameters:
 *   coupling_id        <-- coupling identificator
 *   n_points           <-- number of points to locate
 *   coordinates        <-- coordinates of points to locate (enterlaced)
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_set_points_to_locate_cf, 
           COUPLINGS_SET_POINTS_TO_LOCATE_CF)
  (const char   *coupling_id,
   const int  *l_coupling_id,
   const int    *n_points,
   const double *coordinate
   ARGF_SUPP_CHAINE);
                                          
/*----------------------------------------------------------------------------
 *
 * Define the support mesh for a coupling. The connectivity is ordered if 
 * necessary. The connectivity order is :
 *        - 1D : edges
 *        - 2D : triangles, quadrangles, polygons
 *        - 3D : tetrahedra, pyramids, prism, hexaedra, polyhedra
 *
 * parameters:
 *   coupling_id        <-- coupling identificator
 *   dim                <-- space dimension (1, 2 or 3)
 *   n_vertex           <-- number of vertex
 *   n_elements         <-- number of elements
 *   coordinates        <-- vertex enterlaced coordinates
 *   parent_vertex_num  <-- pointer to parent vertex numbers (or NULL)
 *   connectivity_index <-> polygon face -> vertices index (O to n-1)
 *                          size: n_elements + 1
 *   connectivity       <-> element -> vertex connectivity
 *                          size: connectivity_index[n_elements] 
 * 
 *----------------------------------------------------------------------------*/

// Add couplings_update_mesh

void PROCF(couplings_define_mesh_cf, 
           COUPLINGS_DEFINE_MESH_CF)
  (const char *coupling_id,
   const int  *l_coupling_id,
   const int *dim,
   const int *n_vertex,
   const int *n_element,
   const double *coordinates,
   const int *shared_coordinates_status,
   int *parent_vertex_num,
   int *connectivity_index,
   int *connectivity
   ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 *  Locate points. This is a synchronization point with the coupled application 
 *
 *   coupling_id        <-- coupling identificator
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_locate_cf, 
           COUPLINGS_LOCATE_CF) 
(const char *coupling_id
 ARGF_SUPP_CHAINE); 

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

void PROCF(couplings_exchange_float_cf, 
           COUPLINGS_EXCHANGE_FLOAT_CF)
  (const char      *coupling_id,
   const int       *l_coupling_id,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *exchange_type, 
   const int       *exchange_dimension, 
   const int       *interpolation_type, 
   const int       *time_step, 
   const double    *time_value,
   const char      *sending_field_name,
   const int       *l_sending_field_name,
   const float     *sending_field, 
   const int       *l_receiving_field_name,
   char            *receiving_field_name,
   float           *receiving_field
   ARGF_SUPP_CHAINE);

void PROCF(couplings_exchange_double_cf, 
           COUPLINGS_EXCHANGE_DOUBLE_CF)
  (const char      *coupling_id,
   const int       *l_coupling_id,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *exchange_type, 
   const int       *exchange_dimension, 
   const int       *interpolation_type, 
   const int       *time_step, 
   const double    *time_value,
   const char      *sending_field_name,
   const int       *l_sending_field_name,
   const double    *sending_field, 
   const int       *l_receiving_field_name,
   char            *receiving_field_name,
   double          *receiving_field 
   ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Delete a coupling
 *
 * parameters
 *   coupling_id          <-- Coupling identificator
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_delete_coupling_cf, 
           COUPLINGS_DELETE_COUPLING_CF)
  (const char *coupling_id
   ARGF_SUPP_CHAINE);

/*----------------------------------------------------------------------------
 *
 * Finalize couplings. This is a synchronization point between all applications 
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_finalize_cf, 
           COUPLINGS_FINALIZE_CF) ();

/*----------------------------------------------------------------------------
 *
 * Dump application properties
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_dump_application_properties_f,
           COUPLINGS_DUMP_APPLICATION_PROPERTIES_F) ();

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __COUPLINGS_H__ */
