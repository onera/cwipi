/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_parall.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "couplings.h"
#include "applicationPropertiesDataBase.hxx"
#include "applicationPropertiesDataBase_i.hxx"
#include "couplingDataBase.hxx"
#include "couplingDataBase_i.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Output listing File (C printing)
 *----------------------------------------------------------------------------*/

static FILE* _couplings_output_listing;

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * bft_printf proxy setting for C interface 
 *
 *----------------------------------------------------------------------------*/

int _couplings_print_with_c
(
 const char     *const format,
       va_list         arg_ptr
)
{
  return vfprintf(_couplings_output_listing, format, arg_ptr);
}


/*============================================================================
 * Public function definitions
 *============================================================================*/


/*----------------------------------------------------------------------------
 *
 * Initialize the couplings library.
 * Redirect outputs in a file (Standard output with output_listing = NULL or 
 * output_logical_unit = -1)
 * Create the current communicator application from 'common_comm'.
 *
 * parameters:
 *   common_comm         <-- Common MPI communicator
 *   output_listing      <-- Output listing file (C function)
 *   output_logical_unit <-- Output listing logical unit (Fortran function)
 *   application_name    <-- Current application name 
 *   application_comm    --> Internal MPI communicator for the current 
 *                           application
 *
 * This is a synchronization point between all applications 
 *----------------------------------------------------------------------------*/

void couplings_init
(const MPI_Comm common_comm,
 FILE           *output_listing, 
 const char     *application_name, 
 MPI_Comm       *application_comm)

{
  _couplings_output_listing = output_listing; 
  bft_printf_proxy_set(_couplings_print_with_c);

  couplings::ApplicationPropertiesDataBase & properties = 
    couplings::ApplicationPropertiesDataBase::getInstance();
  // MPI_Comm *global_comm = new MPI_Comm[1];
  //*global_comm = common_comm;
  properties.init(application_name,
                  common_comm,
                  //*global_comm,
                  *application_comm);
}

/*----------------------------------------------------------------------------
 *
 * Add a integer control parameter
 * 
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void couplings_add_local_int_control_parameter(const char *name, int initial_value)
{
  couplings::ApplicationPropertiesDataBase & properties = 
    couplings::ApplicationPropertiesDataBase::getInstance();

  const std::string & nameStr = name; 
  properties.addLocalIntControlParameter(name, initial_value);
}

/*----------------------------------------------------------------------------
 *
 * Add a double control parameter
 * 
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void couplings_add_local_double_control_parameter
(const char *name, 
 double initial_value)
{
  couplings::ApplicationPropertiesDataBase & properties = 
    couplings::ApplicationPropertiesDataBase::getInstance();

  const std::string & nameStr = name; 
  properties.addLocalDoubleControlParameter(name, initial_value);
}

/*----------------------------------------------------------------------------
 *
 * Set a integer control parameter
 * 
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void couplings_set_local_int_control_parameter(const char *name, int value)
{
  couplings::ApplicationPropertiesDataBase & properties = 
    couplings::ApplicationPropertiesDataBase::getInstance();

  const std::string & nameStr = name; 
  properties.setLocalIntControlParameter(name, value);
}

/*----------------------------------------------------------------------------
 *
 * Set a double control parameter
 * 
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void couplings_set_local_double_control_parameter(const char *name, double value)
{
  couplings::ApplicationPropertiesDataBase & properties = 
    couplings::ApplicationPropertiesDataBase::getInstance();

  const std::string & nameStr = name; 
  properties.setLocalDoubleControlParameter(name, value);
}

/*----------------------------------------------------------------------------
 *
 * Get a integer control parameter of the current application
 * 
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

int couplings_get_local_int_control_parameter(const char *name)
{
  couplings::ApplicationPropertiesDataBase & properties = 
    couplings::ApplicationPropertiesDataBase::getInstance();

  const std::string & nameStr = name; 
  return properties.getLocalIntControlParameter(name);
}

/*----------------------------------------------------------------------------
 *
 * Get a double control parameter of the current application
 * 
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

double couplings_get_local_double_control_parameter(const char *name)
{
  couplings::ApplicationPropertiesDataBase & properties = 
    couplings::ApplicationPropertiesDataBase::getInstance();

  const std::string & nameStr = name; 
  return properties.getLocalDoubleControlParameter(name);
}

/*----------------------------------------------------------------------------
 *
 * Delete a current application Int parameter
 * 
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void couplings_delete_local_int_control_parameter(const char *name)
{
  couplings::ApplicationPropertiesDataBase & properties = 
    couplings::ApplicationPropertiesDataBase::getInstance();

  const std::string & nameStr = name; 
  properties.eraseLocalIntControlParameter(name);
}

/*----------------------------------------------------------------------------
 *
 * Delete a current application Int parameter
 * 
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void couplings_delete_local_double_control_parameter(const char *name)
{
  couplings::ApplicationPropertiesDataBase & properties = 
    couplings::ApplicationPropertiesDataBase::getInstance();

  const std::string & nameStr = name; 
  properties.eraseLocalDoubleControlParameter(name);
}

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
 const char *name)
{
  couplings::ApplicationPropertiesDataBase & properties = 
    couplings::ApplicationPropertiesDataBase::getInstance();

  const std::string &application_nameStr = application_name;
  const std::string &nameStr = name;
  return properties.getDistantIntControlParameter(application_nameStr, nameStr);
}

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
 const char *name)
{
  couplings::ApplicationPropertiesDataBase & properties = 
    couplings::ApplicationPropertiesDataBase::getInstance();

  const std::string &application_nameStr = application_name;
  const std::string &nameStr = name;
  return properties.getDistantDoubleControlParameter(application_nameStr, nameStr);
}

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

void couplings_synchronise_control_parameter(const char *application_name)
{
  couplings::ApplicationPropertiesDataBase & properties = 
    couplings::ApplicationPropertiesDataBase::getInstance();

  const std::string &application_nameStr = application_name;
  return properties.mergeParameters(application_nameStr);
}

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
 *                                 - "None"
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

void couplings_create_coupling
( const char  *coupling_name,
  const char  *coupled_application,
  const int entities_dim,
  const double tolerance,
  const couplings_mesh_type_t mesh_type,
  const couplings_solver_type_t solver_type, 
  const int    output_frequency,
  const char  *output_format,
  const char  *output_format_option)
{
  couplings::CouplingDataBase & couplingDataBase = 
    couplings::CouplingDataBase::getInstance();

  couplings::ApplicationPropertiesDataBase & properties = 
    couplings::ApplicationPropertiesDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;
  const std::string &coupled_application_str = coupled_application;

  couplingDataBase.createCoupling(coupling_name, 
                                  properties.getLocalApplicationProperties(),
                                  properties.getDistantApplicationProperties(coupled_application_str),
                                  entities_dim,
                                  tolerance,
                                  solver_type,
                                  output_frequency,
                                  output_format,
                                  output_format_option);
}

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

void couplings_set_points_to_locate
(const char  *coupling_name,
 const int    n_points,
 double coordinate[]) 
{
  couplings::CouplingDataBase & couplingDataBase = 
    couplings::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  couplings::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.setPointsToLocate(n_points, coordinate);
}
                                          
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

void couplings_define_mesh(const char *coupling_name, 
                           const int n_vertex,
                           const int n_element,
                           const double coordinates[],
                           int connectivity_index[],
                           int connectivity[])
{
  couplings::CouplingDataBase & couplingDataBase = 
    couplings::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  couplings::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.defineMesh(n_vertex, 
                      n_element, 
                      coordinates, 
                      connectivity_index, 
                      connectivity);
}

void couplings_add_polyhedra(const char *coupling_name, 
                             const int n_element,
                             int face_index[],
                             int cell_to_face_connectivity[],
                             int face_connectivity_index[],
                             int face_connectivity[])
{
  couplings::CouplingDataBase & couplingDataBase = 
    couplings::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  couplings::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.defineMeshAddPolyhedra(n_element,
                                  face_index,
                                  cell_to_face_connectivity,
                                  face_connectivity_index,
                                  face_connectivity);
}


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

int couplings_exchange
(const char                          *coupling_name,
 const char                          *exchange_name,
 const couplings_field_dimension_t    exchange_dimension, 
 /*const couplings_interpolation_t      interpolation_type,*/ 
 const int                            time_step, 
 const double                         time_value,
 const char                          *sending_field_name,
 const double                        *sending_field, 
 char                                *receiving_field_name,
 double                              *receiving_field)

{
  couplings::CouplingDataBase & couplingDataBase = 
    couplings::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  couplings::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  return coupling.exchange(exchange_name,
                           exchange_dimension, 
                           time_step, 
                           time_value,
                           sending_field_name,
                           sending_field, 
                           receiving_field_name,
                           receiving_field,
                           NULL);
}

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
(const char *coupling_name,
 couplings_interpolation_fct_t * fct)
{
  couplings::CouplingDataBase & couplingDataBase = 
    couplings::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  couplings::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.set_interpolation_function(fct);
}

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
 couplings_not_located_point_treatment_fct_t *const fct) 
{
  std::cout << "  couplings_set_not_located_point_treatment_function not yet implemented" << std::endl;
  }*/

/*----------------------------------------------------------------------------
 *
 * Delete a coupling
 *
 * parameters
 *   coupling_id          <-- Coupling identificator
 *
 *----------------------------------------------------------------------------*/

void couplings_delete_coupling(const char *coupling_name)
{
  couplings::CouplingDataBase & couplingDataBase = 
    couplings::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_name;

  couplingDataBase.deleteCoupling(coupling_name_str);
}

/*----------------------------------------------------------------------------
 *
 * Finalize couplings. This is a synchronization point between all applications 
 *
 *----------------------------------------------------------------------------*/

void couplings_finalize() 
{
  couplings::CouplingDataBase & couplingDataBase = 
    couplings::CouplingDataBase::getInstance();

  couplings::ApplicationPropertiesDataBase & properties = 
    couplings::ApplicationPropertiesDataBase::getInstance();

  const MPI_Comm globalComm = properties.getGlobalComm();


  bft_printf("Finalize couplings\n");
  couplingDataBase.kill();
  properties.kill();

  int flag = 0;
  MPI_Initialized(&flag);

  if (flag != 0) {
    bft_printf_flush();
    MPI_Barrier(globalComm);
    fvm_parall_set_mpi_comm(MPI_COMM_NULL);
    MPI_Finalize();
  }
  bft_printf("Finalize MPI\n");

}

/*----------------------------------------------------------------------------
 *
 * Dump application properties
 *
 *----------------------------------------------------------------------------*/

void couplings_dump_application_properties()
{
  couplings::ApplicationPropertiesDataBase & properties = 
    couplings::ApplicationPropertiesDataBase::getInstance();
  properties.dump();
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
