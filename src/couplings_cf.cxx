/*----------------------------------------------------------------------------
 * Standard C/C++ library headers
 *----------------------------------------------------------------------------*/

#include <cassert>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "couplings_cf.h"
#include "applicationPropertiesDataBase.hxx"
#include "couplings.h"

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

/*----------------------------------------------------------------------------
 * Fortran printing
 *----------------------------------------------------------------------------*/

  void PROCF (printfort, PRINTFORT) (char *buf_print_f, int *msgsize);

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

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Convert a fortran string to a C string
 *
 * parameters:
 *   application_name_f    <-- Fortran string
 *   l_application_name_f  <-- Fortran string length
 *
 * return:
 *   C string
 *
 *----------------------------------------------------------------------------*/

static char *_couplings_fortran_to_c_string(const char *application_name_f, 
                                  const int l_application_name_f)
{
  char *application_name_c = NULL;
  int imin = 0;
  int imax = 0;

  while (imin < l_application_name_f && 
         application_name_f[imin] == ' ') 
    imin++;

  while (imax < l_application_name_f && 
         application_name_f[l_application_name_f-imax-1] == ' ')
    imax++;

  imax = l_application_name_f-imax-1;

  assert(imax >= imin);

  if ((imax == l_application_name_f) || (imin == l_application_name_f)) {
    application_name_c = new char[1];
    application_name_c[0] = '\0';
  }
  else {
    int size = imax - imin + 2;
    application_name_c = new char[size];
    int index = 0;
    for (int k = imin; k <= imax; k++)
      application_name_c[index++] = application_name_f[k];
    application_name_c[index] = '\0';
  }

  return application_name_c;
}

/*----------------------------------------------------------------------------
 *
 * Set bft_printf proxy for Fortran interface
 *
 *----------------------------------------------------------------------------*/

int _couplings_print_with_fortran
(
 const char     *const format,
       va_list         arg_ptr
)
{
 int  line;
 int  msgsize;

 /* Tampon pour impressions depuis du code C : on imprime dans un chaîne
    de caractères, qui sera imprimée vers un fichier par du code Fortran.
    Une fois les impressions Fortran totalement remplacées par des impressions
    C, on pourra supprimer cette étape, mais elle est nécessaire pour l'instant
    afin de pouvoir utiliser les mêmes fichiers de sortie */

#undef BUF_PRINT_F_SIZE
#define BUF_PRINT_F_SIZE 16384

 static char buf_print_f[BUF_PRINT_F_SIZE];

 /* Impression dans le tampon */

#if defined  __STDC_VERSION__ 
  msgsize = vsprintf (buf_print_f, format, arg_ptr);
#else
  msgsize = vsnprintf (buf_print_f, BUF_PRINT_F_SIZE, format, arg_ptr);
#endif

  line = __LINE__ - 1;

  if (msgsize == -1 || msgsize > BUF_PRINT_F_SIZE - 1) {
    fprintf(stderr,
            "Fatal error: bft_printf() called on a message of size %d\n"
            "whereas the print buffer is of size %d.",
            msgsize, BUF_PRINT_F_SIZE);

    /* Try to force segmentation fault (to call signal handlers);
       as stack has most likely been corrupted, this is the most
       "similar" error that allows for portable handling. */
    {
      int *_force_err = NULL;
      *_force_err = 0;
    }
    exit(1);
  }

  /* Impression effective par le code Fortran */

  PROCF (printfort, PRINTFORT) (buf_print_f, &msgsize);

  return msgsize;
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
 * This is a synchronization between all applications 
 *----------------------------------------------------------------------------*/

void PROCF(couplings_init_cf, COUPLINGS_INIT_CF)
  (const int  *common_comm,
   const int  *output_logical_unit,
   const char *application_name_f,
   const int  *l_application_name,
   int        *application_comm
   ARGF_SUPP_CHAINE)
{
   char *application_name_c = _couplings_fortran_to_c_string(application_name_f, 
                                                             *l_application_name); 

   bft_printf_proxy_set(_couplings_print_with_fortran);

   couplings::ApplicationPropertiesDataBase & properties = 
     couplings::ApplicationPropertiesDataBase::getInstance();
   properties.init(application_name_c, 
                   static_cast <const MPI_Comm&> (*common_comm), 
                   static_cast <MPI_Comm&> (*application_comm));
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

void PROCF(couplings_add_local_int_control_parameter_cf, 
           COUPLINGS_ADD_LOCAL_INT_CONTROL_PARAMETER_CF)
  (const char *name, 
   const int  *l_name,
   int *initial_value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _couplings_fortran_to_c_string(name, *l_name);
  couplings_add_local_int_control_parameter(nameC, *initial_value);
  delete nameC;
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

void PROCF(couplings_add_local_double_control_parameter_cf, 
           COUPLINGS_ADD_LOCAL_DOUBLE_CONTROL_PARAMETER_CF)
  (const char *name, 
   const int  *l_name,
   double *initial_value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _couplings_fortran_to_c_string(name, *l_name);
  couplings_add_local_double_control_parameter(nameC, *initial_value);
  delete nameC;
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

void PROCF(couplings_set_local_int_control_parameter_cf, 
           COUPLINGS_SET_LOCAL_INT_CONTROL_PARAMETER_CF)
  (const char *name, 
   const int  *l_name,
   int *value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _couplings_fortran_to_c_string(name, *l_name);
  couplings_set_local_int_control_parameter(nameC, *value);
  delete nameC;
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

void PROCF(couplings_set_local_double_control_parameter_cf, 
           COUPLINGS_SET_LOCAL_DOUBLE_CONTROL_PARAMETER_CF)
  (const char *name, 
   const int  *l_name,
   double *value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _couplings_fortran_to_c_string(name, *l_name);
  couplings_set_local_double_control_parameter(nameC, *value);
  delete nameC;
}

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
   ARGF_SUPP_CHAINE)
{
  char* nameC = _couplings_fortran_to_c_string(name, *l_name);
  *value = couplings_get_local_int_control_parameter(nameC);
  delete nameC;
}

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
   ARGF_SUPP_CHAINE)
{
  char* nameC = _couplings_fortran_to_c_string(name, *l_name);
  *value = couplings_get_local_double_control_parameter(nameC);
  delete nameC;
}

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
   ARGF_SUPP_CHAINE)
{
  char* nameC = _couplings_fortran_to_c_string(name, *l_name);
  couplings_delete_local_int_control_parameter(nameC);
  delete nameC;
}

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
   ARGF_SUPP_CHAINE)
{
  char* nameC = _couplings_fortran_to_c_string(name, *l_name);
  couplings_delete_local_double_control_parameter(nameC);
  delete nameC;
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

void PROCF(couplings_get_distant_int_control_parameter_cf, 
           COUPLINGS_GET_DISTANT_INT_CONTROL_PARAMETER_CF)
  (const char *application_name,
   const int  *l_application_name,
   const char *name,
   const int  *l_name,
   int *value
   ARGF_SUPP_CHAINE)
{
  char *application_nameC = 
    _couplings_fortran_to_c_string(application_name, *l_application_name);
  char *nameC = _couplings_fortran_to_c_string(name, *l_name);

  *value = couplings_get_distant_int_control_parameter(application_nameC, 
                                                       nameC);

  delete nameC;
  delete application_nameC;
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

void PROCF(couplings_get_distant_double_control_parameter_cf, 
           COUPLINGS_GET_DISTANT_DOUBLE_CONTROL_PARAMETER_CF)
  (const char *application_name,
   const int  *l_application_name,
   const char *name,
   const int  *l_name,
   double *value
   ARGF_SUPP_CHAINE)
{
  char *application_nameC = 
    _couplings_fortran_to_c_string(application_name, *l_application_name);
  char *nameC = _couplings_fortran_to_c_string(name, *l_name);

  *value = couplings_get_distant_double_control_parameter(application_nameC, 
                                                          nameC);

  delete nameC;
  delete application_nameC;
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

void PROCF(couplings_synchronise_control_parameter_cf, 
           COUPLINGS_SYNCHRONISE_CONTROL_PARAMETER_CF)
  (const char *application_name,
   const int  *l_application_name
   ARGF_SUPP_CHAINE)
{
  char *application_nameC = 
    _couplings_fortran_to_c_string(application_name, *l_application_name);

  couplings_synchronise_control_parameter(application_nameC);

  delete application_nameC;
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
   ARGF_SUPP_CHAINE)
{
  std::cout << " couplings_create_coupling not yet implemented" << std::endl;
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

void PROCF(couplings_set_points_to_locate_cf, 
           COUPLINGS_SET_POINTS_TO_LOCATE_CF)
  (const char   *coupling_id,
   const int  *l_coupling_id,
   const int    *n_points,
   const double *coordinate
   ARGF_SUPP_CHAINE)
{
  std::cout << " couplings_set_points_to_locate_f not yet implemented" << std::endl;
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
   ARGF_SUPP_CHAINE)
{
  std::cout << " couplings_defne_mesh_f not yet implemented" << std::endl;
}

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
 ARGF_SUPP_CHAINE) 
{
  std::cout << " couplings_locate_f not yet implemented" << std::endl;
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
   ARGF_SUPP_CHAINE)
{
  std::cout << " couplings_exchange_float_f not yet implemented" << std::endl;
}

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
   ARGF_SUPP_CHAINE)
{
  std::cout << " couplings_exchange_double_f not yet implemented" << std::endl;
}

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
   ARGF_SUPP_CHAINE)
{
  std::cout << "  couplings_delete_coupling_f not yet implemented" << std::endl;
}

/*----------------------------------------------------------------------------
 *
 * Finalize couplings. This is a synchronization point between all applications 
 *
 *----------------------------------------------------------------------------*/


void PROCF(couplings_finalize_cf, 
           COUPLINGS_FINALIZE_CF) ()
{
  std::cout << "  couplings_finalize not yet implemented" << std::endl;
}

/*----------------------------------------------------------------------------
 *
 * Dump application properties
 *
 *----------------------------------------------------------------------------*/

void PROCF(couplings_dump_application_properties_f,
           COUPLINGS_DUMP_APPLICATION_PROPERTIES_F) ()
{
  couplings_dump_application_properties();
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
