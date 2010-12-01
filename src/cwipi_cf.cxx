/*----------------------------------------------------------------------------
 * Standard C/C++ library headers
 *----------------------------------------------------------------------------*/

#include <cassert>
#include <cstring>

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

#include "cwipi_cf.h"
#include "applicationPropertiesDataBase.hxx"
#include "couplingDataBase.hxx"
#include "couplingDataBase_i.hxx"
#include "cwipi.h"
#include "cwipi_config.h"
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

static char *_cwipi_fortran_to_c_string(const char *application_name_f,
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

int _cwipi_print_with_fortran
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
 * Initialize the cwipi library.
 * Redirect outputs in a file (Standard output with output_listing = NULL or
 * output_logical_unit = -1)
 * Create the current communicator application from 'common_comm'.
 *
 * parameters:
 *   common_comm         <-- Common MPI communicator
 *   application_name    <-- Current application name
 *   application_comm    --> Internal MPI communicator for the current
 *                           application
 *
 * This is a synchronization between all applications
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_init_cf, CWIPI_INIT_CF)
  (MPI_Fint  *common_fcomm,
   const char *application_name_f,
   const int  *l_application_name,
   MPI_Fint   *application_fcomm
   ARGF_SUPP_CHAINE)
{
  MPI_Comm common_comm = MPI_Comm_f2c(*common_fcomm);

  MPI_Comm application_comm = MPI_COMM_NULL;

  char *application_name_c = _cwipi_fortran_to_c_string(application_name_f,
                                                            *l_application_name);

  bft_printf("\ncwipi "CWIPI_VERSION" initializing\n");
  bft_printf("------------------------\n\n");

  cwipi::ApplicationPropertiesDataBase & properties =
    cwipi::ApplicationPropertiesDataBase::getInstance();

  application_comm = properties.init(application_name_c,
                                     common_comm);

  *application_fcomm = MPI_Comm_c2f(application_comm);


  delete[] application_name_c;
}

/*----------------------------------------------------------------------------
 *
 * Set up the file used for the output listing
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_set_output_listing_cf, CWIPI_SET_OUTPUT_LISTING_CF)()
{
  bft_printf_proxy_set(_cwipi_print_with_fortran);
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

void PROCF(cwipi_add_loc_int_ctrl_param_cf,
           CWIPI_ADD_LOC_INT_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   int *initial_value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  cwipi_add_local_int_control_parameter(nameC, *initial_value);
  delete[] nameC;
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

void PROCF(cwipi_add_loc_dbl_ctrl_param_cf,
           CWIPI_ADD_LOC_DBL_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   double *initial_value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  cwipi_add_local_double_control_parameter(nameC, *initial_value);
  delete[] nameC;
}

/*----------------------------------------------------------------------------
 *
 * Add a string control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    initial_value  <-- initial value
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_add_loc_str_ctrl_param_cf,
           CWIPI_ADD_LOC_STR_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   char *initial_value,
   int *l_value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  char* valueC = _cwipi_fortran_to_c_string(initial_value, *l_value);
  cwipi_add_local_string_control_parameter(nameC, valueC);
  delete[] nameC;
  delete[] valueC;
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

void PROCF(cwipi_set_loc_int_ctrl_param_cf,
           CWIPI_SET_LOC_INT_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   int *value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  cwipi_set_local_int_control_parameter(nameC, *value);
  delete[] nameC;
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

void PROCF(cwipi_set_loc_dbl_ctrl_param_cf,
           CWIPI_SET_LOC_DBL_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   double *value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  cwipi_set_local_double_control_parameter(nameC, *value);
  delete[] nameC;
}

/*----------------------------------------------------------------------------
 *
 * Set a string control parameter
 *
 * parameters
 *    name           <-- parameter name
 *    value          <-- value
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_set_loc_str_ctrl_param_cf,
           CWIPI_SET_LOC_STR_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   char *initial_value,
   int *l_value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  char* valueC = _cwipi_fortran_to_c_string(initial_value, *l_value);
  cwipi_set_local_string_control_parameter(nameC, valueC);
  delete[] nameC;
  delete[] valueC;
}

/*----------------------------------------------------------------------------
 *
 * Get a integer control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_loc_int_ctrl_param_cf,
           CWIPI_GET_LOC_INT_CRTL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   int *value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  *value = cwipi_get_local_int_control_parameter(nameC);
  delete[] nameC;
}

/*----------------------------------------------------------------------------
 *
 * Get a double control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_loc_dbl_ctrl_param_cf,
           CWIPI_GET_LOC_DBL_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   double *value
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  *value = cwipi_get_local_double_control_parameter(nameC);
  delete[] nameC;
}

/*----------------------------------------------------------------------------
 *
 * Get a string control parameter of the current application
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_loc_str_ctrl_param_cf,
           CWIPI_GET_LOC_STR_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name,
   char *value_str_f,
   int *l_value_str_f
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  const char* value_str_c = cwipi_get_local_string_control_parameter(nameC);
  const int l_value_str_c = strlen(value_str_c);

  int i = 0;
  while (i < l_value_str_c && i < *l_value_str_f) {
    value_str_f[i] = value_str_c[i];
    i+=1;
  }

  while (i < *l_value_str_f) {
    value_str_f[i] = ' ';
    i+=1;
  }
  delete[] nameC;
}

/*----------------------------------------------------------------------------
 *
 * Delete a current application parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_del_loc_int_ctrl_param_cf,
           CWIPI_DEL_LOC_INT_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  cwipi_delete_local_int_control_parameter(nameC);
  delete[] nameC;
}

/*----------------------------------------------------------------------------
 *
 * Delete a current application parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_del_loc_dbl_ctrl_param_cf,
           CWIPI_DEL_LOC_DBL_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  cwipi_delete_local_double_control_parameter(nameC);
  delete[] nameC;
}

/*----------------------------------------------------------------------------
 *
 * Delete a current application parameter
 *
 * parameters
 *    name           <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_del_loc_str_ctrl_param_cf,
           CWIPI_DEL_LOC_STR_CTRL_PARAM_CF)
  (const char *name,
   const int  *l_name
   ARGF_SUPP_CHAINE)
{
  char* nameC = _cwipi_fortran_to_c_string(name, *l_name);
  cwipi_delete_local_string_control_parameter(nameC);
  delete[] nameC;
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

void PROCF(cwipi_get_dis_int_ctrl_param_cf,
           CWIPI_GET_DIS_INT_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name,
   const char *name,
   const int  *l_name,
   int *value
   ARGF_SUPP_CHAINE)
{
  char *application_nameC =
    _cwipi_fortran_to_c_string(application_name, *l_application_name);
  char *nameC = _cwipi_fortran_to_c_string(name, *l_name);

  *value = cwipi_get_distant_int_control_parameter(application_nameC,
                                                       nameC);

  delete[] nameC;
  delete[] application_nameC;
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

void PROCF(cwipi_get_dis_dbl_ctrl_param_cf,
           CWIPI_GET_DIS_DBL_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name,
   const char *name,
   const int  *l_name,
   double *value
   ARGF_SUPP_CHAINE)
{
  char *application_nameC =
    _cwipi_fortran_to_c_string(application_name, *l_application_name);
  char *nameC = _cwipi_fortran_to_c_string(name, *l_name);

  *value = cwipi_get_distant_double_control_parameter(application_nameC,
                                                          nameC);

  delete[] nameC;
  delete[] application_nameC;
}
/*
 * ----------------------------------------------------------------------------
 *
 * Get a string control parameter of a other application
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_dis_str_ctrl_param_cf,
           CWIPI_GET_DIS_STR_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name,
   const char *name,
   const int  *l_name,
   char *value_str_f,
   int  *l_value_str_f
   ARGF_SUPP_CHAINE)

{
  char *application_nameC = _cwipi_fortran_to_c_string(application_name, *l_application_name);
  char *nameC = _cwipi_fortran_to_c_string(name, *l_name);

  const char *value_str_c = cwipi_get_distant_string_control_parameter(application_nameC, nameC);
  const int l_value_str_c = strlen(value_str_c);

  int i = 0;
  while (i < l_value_str_c && i < *l_value_str_f) {
    value_str_f[i] = value_str_c[i];
    i+=1;
  }

  while (i < *l_value_str_f) {
    value_str_f[i] = ' ';
    i+=1;
  }

  delete[] nameC;
  delete[] application_nameC;
}

/*----------------------------------------------------------------------------
 *
 * Synchronize local control parameters with an other application.
 *  This is a synchronization point with this second application
 *
 * parameters
 *    application_name    <-- application name
 *    name                <-- parameter name
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_synch_ctrl_param_cf,
           CWIPI_SYNCH_CTRL_PARAM_CF)
  (const char *application_name,
   const int  *l_application_name
   ARGF_SUPP_CHAINE)
{
  char *application_nameC =
    _cwipi_fortran_to_c_string(application_name, *l_application_name);

  cwipi_synchronize_control_parameter(application_nameC);

  delete[] application_nameC;
}

/*----------------------------------------------------------------------------
 *
 * Create a coupling object
 *
 * parameters:
 *   couplingName            <-- Coupling identifier
 *   couplingType            <-- Coupling type
 *   cplAppli                <-- Coupled application name
 *   entitiesDim             <-- Mesh entities dimension (1, 2 or 3)
 *   tolerance               <-- Geometric tolerance to locate
 *   meshT                   <-- CWIPI_STATIC_MESH
 *                               CWIPI_MOBILE_MESH (not implemented yet)
 *   solverT                 <-- CWIPI_SOLVER_CELL_CENTER
 *                               CWIPI_SOLVER_CELL_VERTEX
 *   outputFreq              <-- Output frequency
 *   outputFmt               <-- Output format to visualize exchanged fields
 *                               on the coupled mesh. Choice between :
 *                                 - "EnSight Gold"
 *                                 - "MED_ficher"
 *                                 - "CGNS"
 *   outputFmtOpt            <-- Output options
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
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_create_coupling_cf,
           CWIPI_CREATE_COUPLING_CF)
( const char *coupling_name,
  const int  *l_coupling_name,
  const int  *coupling_type,
  const char *coupled_application,
  const int  *l_coupled_application,
  const int  *entities_dim,
  const double *tolerance,
  const int *mesh_type,
  const int *solver_type,
  const int  * output_frequency,
  const char  *output_format,
  const int  *l_output_format,
  const char  *output_format_option,
  const int  *l_output_format_option
  ARGF_SUPP_CHAINE)

{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  char *coupled_applicationC =
    _cwipi_fortran_to_c_string(coupled_application, *l_coupled_application);

  char *output_formatC =
    _cwipi_fortran_to_c_string(output_format, *l_output_format);

  char *output_format_optionC =
    _cwipi_fortran_to_c_string(output_format_option, *l_output_format_option);

  cwipi_create_coupling(coupling_nameC,
                            (cwipi_coupling_type_t) *coupling_type,
                            coupled_applicationC,
                            *entities_dim,
                            *tolerance,
                            (cwipi_mesh_type_t) *mesh_type,
                            (cwipi_solver_type_t) *solver_type,
                            *output_frequency,
                            output_formatC,
                            output_format_optionC);

  delete[] coupling_nameC;
  delete[] coupled_applicationC;
  delete[] output_formatC;
  delete[] output_format_optionC;
}

/*----------------------------------------------------------------------------
 *
 * Set points to locate. This function must be called if the points to locate
 * do not correspond to :
 *        - vertices for NATURE_NODE nature
 *        - cell center for NATURE_CELL_CENTER nature
 *
 * parameters:
 *   coupling_id        <-- coupling identifier
 *   n_points           <-- number of points to locate
 *   coordinates        <-- coordinates of points to locate (enterlaced)
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_set_points_to_locate_cf,
           CWIPI_SET_POINTS_TO_LOCATE_CF)
  (const char   *coupling_name,
   const int  *l_coupling_name,
   const int    *n_points,
   double *coordinate
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi_set_points_to_locate(coupling_nameC,
                                 *n_points,
                                 coordinate);
  delete[] coupling_nameC;
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
 *   coupling_id        <-- coupling identifier
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

void PROCF(cwipi_define_mesh_cf,
           CWIPI_DEFINE_MESH_CF)
  (const char *coupling_name,
   const int  *l_coupling_name,
   const int *n_vertex,
   const int *n_element,
   double *coordinates,
   int *connectivity_index,
   int *connectivity
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi_define_mesh(coupling_nameC,
                        *n_vertex,
                        *n_element,
                        coordinates,
                        connectivity_index,
                        connectivity);
  delete[] coupling_nameC;
}

void PROCF(cwipi_add_polyhedra_cf,
           CWIPI_ADD_POLYHEDRA_CF)
  (const char *coupling_name,
   const int  *l_coupling_name,
   const int *n_element,
   int *face_index,
   int *cell_to_face_connectivity,
   int *face_connectivity_index,
   int *face_connectivity
   ARGF_SUPP_CHAINE)

{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi_add_polyhedra(coupling_name,
                          *n_element,
                          face_index,
                          cell_to_face_connectivity,
                          face_connectivity_index,
                          face_connectivity);
  delete[] coupling_nameC;
}

/*----------------------------------------------------------------------------
 *
 * Location completion.
 * This is a synchronization point with the coupled application
 *
 * parameters
 *   coupling_name        <-- Coupling identifier
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_locate_cf, CWIPI_LOCATE_CF) (const char *coupling_name,
                                                      const int  *l_coupling_name
                                                      ARGF_SUPP_CHAINE)
{

  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  coupling.locate();

  delete[] coupling_nameC;
}

/*----------------------------------------------------------------------------
 *
 * Get located points location
 *
 * parameters
 *   coupling_name        <-- Coupling identifier
 *   location             --> located points location
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_distant_location_cf,
           CWIPI_GET_DISTANT_LOCATION_CF) (const char *coupling_name,
                                       const int  *l_coupling_name,
                                       int *location
                                       ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  const unsigned int *locationC = (const unsigned int *) coupling.getDistantLocation();
  const int nDistantPoint = coupling.getNDistantPoint();
  for (int i = 0; i < nDistantPoint; i++)
    location[i] = locationC[i];

  delete[] coupling_nameC;
}

/*----------------------------------------------------------------------------
 *
 * Get number of located distant points
 *
 * parameters
 *   coupling_name           <-- Coupling identifier
 *   n_located_distant_Points --> Number of located distant points
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_n_located_dist_pts_cf,
           CWIPI_GET_N_LOCATED_DIST_PTS_CF) (const char *coupling_name,
                                                       const int  *l_coupling_name,
                                                       int *n_located_distant_Points
                                                       ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  *n_located_distant_Points = coupling.getNDistantPoint();

  delete[] coupling_nameC;
}


/*----------------------------------------------------------------------------
 *
 * Get located points barycentric coordinates index
 *
 * parameters
 *   coupling_name                <-- Coupling identifier
 *   barycentricCoordinatesIndex  --> located points barycentric
 *                                    coordinates index
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_dis_bary_coord_idx_cf,
           CWIPI_GET_DIS_BARY_COORD_IDX_CF) (const char *coupling_name,
                                       const int  *l_coupling_name,
                                       int *barycentricCoordinatesIndex
                                       ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  const int* barycentricCoordinatesIndexC = coupling.getDistantBarycentricCoordinatesIndex();
  const int nDistantPoint = coupling.getNDistantPoint();

  for (int i = 0; i < nDistantPoint + 1; i++)
    barycentricCoordinatesIndex[i] = barycentricCoordinatesIndexC[i];

  delete[] coupling_nameC;
}

/*----------------------------------------------------------------------------
 *
 * Get number of located points
 *
 * parameters
 *   coupling_name        <-- Coupling identifier
 *   n_located_Points     --> Number of located points
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_n_located_pts_cf,
           CWIPI_GET_N_LOCATED_PTS_CF) (const char *coupling_name,
                                               const int  *l_coupling_name,
                                               int *n_located_points
                                               ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  *n_located_points = coupling.getNLocatedPoint();

  delete[] coupling_nameC;
}

/*----------------------------------------------------------------------------
 *
 * Get number of not located points
 *
 * parameters
 *   coupling_name        <-- Coupling identifier
 *   n_not_located_Points --> Number of not located points
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_n_not_located_pts_cf,
           CWIPI_GET_N_NOT_LOCATED_PTS_CF) (const char *coupling_name,
                                                   const int  *l_coupling_name,
                                                   int *n_not_located_points
                                                   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  *n_not_located_points = coupling.getNNotlocatedPoint();

  delete[] coupling_nameC;
}

/*----------------------------------------------------------------------------
 *
 * Get located points barycentric coordinates
 *
 * parameters
 *   coupling_name                <-- Coupling identifier
 *   barycentricCoordinatesIndex  --> located points barycentric
 *                                    coordinates
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_dis_bary_coord_cf,
           CWIPI_GET_DIS_BARY_COORD_CF) (const char *coupling_name,
                                       const int  *l_coupling_name,
                                       double *barycentricCoordinates
                                       ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling =
    couplingDataBase.getCoupling(coupling_name_str);

  const int* barycentricCoordinatesIndexC = coupling.getDistantBarycentricCoordinatesIndex();
  const double* barycentricCoordinatesC = coupling.getDistantBarycentricCoordinates();
  const int nDistantPoint = coupling.getNDistantPoint();

  for (int i = 0; i < barycentricCoordinatesIndexC[nDistantPoint]; i++)
    barycentricCoordinates[i] = barycentricCoordinatesC[i];

  delete[] coupling_nameC;
}

/*----------------------------------------------------------------------------
 *
 * Exchange data with the coupled application. This is a synchronization point
 * with the coupled application
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   exchange_name        <-- Exchange name
 *   exchange_type        <-- Exchange type
 *   exchange_dimension   <-- Dimension of exchanged data :
 *                            - CWIPI_DIMENSION_SCALAR
 *                            - CWIPI_DIMENSION_INTERLACED_VECTOR
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

void PROCF(cwipi_exch_with_user_itp_cf,
           CWIPI_EXCH_WITH_USER_ITP_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *stride,
   const int       *n_step,
   const double    *time_value,
   const char      *sending_field_name,
   const int       *l_sending_field_name,
   const double    *sending_field,
   char            *receiving_field_name,
   const int       *l_receiving_field_name,
   double          *receiving_field,
   void            *ptFortranInterpolationFct,
   int             *n_not_located_points,
   int             *exchange_status
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  char *exchange_nameC =
    _cwipi_fortran_to_c_string(exchange_name, *l_exchange_name);

  char *sending_field_nameC =
    _cwipi_fortran_to_c_string(sending_field_name, *l_sending_field_name);

  char *receiving_field_nameC =
    _cwipi_fortran_to_c_string(receiving_field_name, *l_receiving_field_name);


  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  *exchange_status =  coupling.exchange(exchange_nameC,
                                        *stride,
                                        *n_step,
                                        *time_value,
                                        sending_field_nameC,
                                        sending_field,
                                        receiving_field_nameC,
                                        receiving_field,
                                        ptFortranInterpolationFct);

  *n_not_located_points = coupling.getNNotlocatedPoint();

  delete[] coupling_nameC;
  delete[] exchange_nameC;
  delete[] sending_field_nameC;
  delete[] receiving_field_nameC;
}


void PROCF(cwipi_exchange_cf,
           CWIPI_EXCHANGE_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *stride,
   const int       *n_step,
   const double    *time_value,
   const char      *sending_field_name,
   const int       *l_sending_field_name,
   const double    *sending_field,
   char            *receiving_field_name,
   const int       *l_receiving_field_name,
   double          *receiving_field,
   int             *n_not_located_points,
   int             *exchange_status
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  char *exchange_nameC =
    _cwipi_fortran_to_c_string(exchange_name, *l_exchange_name);

  char *sending_field_nameC =
    _cwipi_fortran_to_c_string(sending_field_name, *l_sending_field_name);

  char *receiving_field_nameC =
    _cwipi_fortran_to_c_string(receiving_field_name, *l_receiving_field_name);


  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  *exchange_status =  coupling.exchange(exchange_nameC,
                                        *stride,
                                        *n_step,
                                        *time_value,
                                        sending_field_nameC,
                                        sending_field,
                                        receiving_field_nameC,
                                        receiving_field,
                                        NULL);

  *n_not_located_points = coupling.getNNotlocatedPoint();

  delete[] coupling_nameC;
  delete[] exchange_nameC;
  delete[] sending_field_nameC;
  delete[] receiving_field_nameC;
}

void PROCF(cwipi_receive_cf,
           CWIPI_RECEIVE_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *stride,
   const int       *n_step,
   const double    *time_value,
   char            *receiving_field_name,
   const int       *l_receiving_field_name,
   double          *receiving_field,
   int             *n_not_located_points,
   int             *exchange_status
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  char *exchange_nameC =
    _cwipi_fortran_to_c_string(exchange_name, *l_exchange_name);

  char *receiving_field_nameC =
    _cwipi_fortran_to_c_string(receiving_field_name, *l_receiving_field_name);


  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  *exchange_status =  coupling.exchange(exchange_nameC,
                                        *stride,
                                        *n_step,
                                        *time_value,
                                        NULL,
                                        NULL,
                                        receiving_field_nameC,
                                        receiving_field,
                                        NULL);
  *n_not_located_points = coupling.getNNotlocatedPoint();

  delete[] coupling_nameC;
  delete[] exchange_nameC;
  delete[] receiving_field_nameC;
}

void PROCF(cwipi_send_with_user_itp_cf,
           CWIPI_SEND_WITH_USER_ITP_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *stride,
   const int       *n_step,
   const double    *time_value,
   const char      *sending_field_name,
   const int       *l_sending_field_name,
   const double    *sending_field,
   void            *ptFortranInterpolationFct,
   int             *exchange_status
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  char *exchange_nameC =
    _cwipi_fortran_to_c_string(exchange_name, *l_exchange_name);

  char *sending_field_nameC =
    _cwipi_fortran_to_c_string(sending_field_name, *l_sending_field_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  *exchange_status =  coupling.exchange(exchange_nameC,
                                        *stride,
                                        *n_step,
                                        *time_value,
                                        sending_field_nameC,
                                        sending_field,
                                        NULL,
                                        NULL,
                                        ptFortranInterpolationFct);
  delete[] coupling_nameC;
  delete[] exchange_nameC;
  delete[] sending_field_nameC;
}

void PROCF(cwipi_send_cf,
           CWIPI_SEND_CF)
  (const char      *coupling_name,
   const int       *l_coupling_name,
   const char      *exchange_name,
   const int       *l_exchange_name,
   const int       *stride,
   const int       *n_step,
   const double    *time_value,
   const char      *sending_field_name,
   const int       *l_sending_field_name,
   const double    *sending_field,
   int             *exchange_status
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  char *exchange_nameC =
    _cwipi_fortran_to_c_string(exchange_name, *l_exchange_name);

  char *sending_field_nameC =
    _cwipi_fortran_to_c_string(sending_field_name, *l_sending_field_name);

  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  *exchange_status =  coupling.exchange(exchange_nameC,
                                        *stride,
                                        *n_step,
                                        *time_value,
                                        sending_field_nameC,
                                        sending_field,
                                        NULL,
                                        NULL,
                                        NULL);
  delete[] coupling_nameC;
  delete[] exchange_nameC;
  delete[] sending_field_nameC;
}

/*----------------------------------------------------------------------------
 *
 * Delete a coupling
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_delete_coupling_cf,
           CWIPI_DELETE_COUPLING_CF)
  (const char *coupling_name,
   const int       *l_coupling_name
   ARGF_SUPP_CHAINE)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  cwipi_delete_coupling(coupling_nameC);
  delete[] coupling_nameC;
}

/*----------------------------------------------------------------------------
 *
 * Finalize cwipi. This is a synchronization point between all applications
 *
 *----------------------------------------------------------------------------*/


void PROCF(cwipi_finalize_cf,
           CWIPI_FINALIZE_CF) ()
{
  cwipi_finalize();
}

/*----------------------------------------------------------------------------
 *
 * Dump application properties
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_dump_appli_properties_cf,
           CWIPI_DUMP_APPLI_PROPERTIES_CF) ()
{
  cwipi_dump_application_properties();
}

/*----------------------------------------------------------------------------
 *
 * Get not located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   notLocatedPoints     --> Not located points
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_not_located_pts_cf,
           CWIPI_GET_NOT_LOCATED_PTS_CF)
  (const char *coupling_name,
   const int  *l_coupling_name,
   int *notLocatedPoints)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  const int n_not_located_points = coupling.getNNotlocatedPoint();
  const int *notLocatedPointsC = coupling.getNotlocatedPoint();

  for( int i = 0; i <  n_not_located_points; i++)
    notLocatedPoints[i] = notLocatedPointsC[i];

  delete[] coupling_nameC;

}

/*----------------------------------------------------------------------------
 *
 * Get located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   notLocatedPoints     --> Not located points
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_located_pts_cf,
           CWIPI_GET_LOCATED_PTS_CF)
  (const char *coupling_name,
   const int  *l_coupling_name,
   int *locatedPoints)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_name, *l_coupling_name);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  const int n_located_points = coupling.getNLocatedPoint();
  const int *locatedPointsC = coupling.getLocatedPoint();

  for( int i = 0; i < n_located_points; i++)
    locatedPoints[i] = locatedPointsC[i];

  delete[] coupling_nameC;

}

/*----------------------------------------------------------------------------
 *
 * Get distant elements that contain located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   elements             --> Element that contain located points
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_elt_cont_cf,
           CWIPI_GET_ELT_CONT_CF)
( const char *coupling_id,
  const int  *l_coupling_id,
  int *elements)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_id, *l_coupling_id);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  const int n_located_points = coupling.getNLocatedPoint();
  const int *elementsC = coupling.getElementContaining();

  for( int i = 0; i < n_located_points; i++)
    elements[i] = elementsC[i];

  delete[] coupling_nameC;
}

/*----------------------------------------------------------------------------
 *
 * Get number of vertices of distant elements that contain located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   n_vertices           --> Number of vertices of element that contain
 *                            located point
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_elt_cont_n_vtx_cf,
           CWIPI_GET_ELT_CONT_N_VTX_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 int *n_vertices)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_id, *l_coupling_id);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  const int n_located_points = coupling.getNLocatedPoint();
  const int *n_verticesC = coupling.getElementContainingNVertex();

  for( int i = 0; i < n_located_points + 1; i++)
    n_vertices[i] = n_verticesC[i];

  delete[] coupling_nameC;
}

/*----------------------------------------------------------------------------
 *
 * Get vertices id of distant elements that contain located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   vertices             --> Vertices id
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_elt_cont_vtx_cf,
           CWIPI_GET_ELT_CONT_VTX_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 int *vertices)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_id, *l_coupling_id);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  const int n_located_points = coupling.getNLocatedPoint();
  const int size = coupling.getElementContainingNVertex()[n_located_points];

  const int *verticesC = coupling.getElementContainingVertex();

  for( int i = 0; i < size; i++)
    vertices[i] = verticesC[i];

  delete[] coupling_nameC;
}

/*----------------------------------------------------------------------------
 *
 * Get vertices coordinates of distant elements that contain located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   coordinates          --> Vertices coordinates
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_elt_cont_vtx_coo_cf,
           CWIPI_GET_ELT_CONT_VTX_COO_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 double* coordinates)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_id, *l_coupling_id);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  const int n_located_points = coupling.getNLocatedPoint();
  const int size = coupling.getElementContainingNVertex()[n_located_points];

  const double *coordinatesC = coupling.getElementContainingVertexCoords();

  for( int i = 0; i < 3*size; i++)
    coordinates[i] = coordinatesC[i];

  delete[] coupling_nameC;
}

/*----------------------------------------------------------------------------
 *
 * Get barycentric coords in distant elements for located points
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   coordinates          --> Barycentric coordinates
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_elt_cont_bar_coo_cf,
           CWIPI_GET_ELT_CONT_BAR_COO_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 double *coordinates)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_id, *l_coupling_id);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  const int n_located_points = coupling.getNLocatedPoint();

  const double *coordinatesC = coupling.getElementContainingBarycentricCoordinates();

  for( int i = 0; i < n_located_points; i++)
    coordinates[i] = coordinatesC[i];

  delete[] coupling_nameC;
}

/*----------------------------------------------------------------------------
 *
 * For each located point get the MPI rank of distant element
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   MPIranks             --> MPI ranks that contains located point
 *
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_get_elt_cont_mpi_rank_cf,
           CWIPI_GET_ELT_CONT_MPI_RANK_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 int *MPIrank)
{
  cwipi::CouplingDataBase & couplingDataBase =
    cwipi::CouplingDataBase::getInstance();

  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_id, *l_coupling_id);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi::Coupling& coupling = couplingDataBase.getCoupling(coupling_name_str);

  const int n_located_points = coupling.getNLocatedPoint();
  const int *MPIrankC = coupling.getElementContainingMPIrank();

  for( int i = 0; i < n_located_points; i++)
    MPIrank[i] = MPIrankC[i];

  delete[] coupling_nameC;
}

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

void PROCF(cwipi_exch_cellvtxfd_eltcont_cf,
           CWIPI_EXCH_CELLVTXFD_ELTCONT_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 double *sendingField,
 double *receivingField,
 const int *stride)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_id, *l_coupling_id);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi_exchange_cell_vertex_field_of_element_containing(coupling_nameC,
                                                         sendingField,
                                                         receivingField,
                                                         *stride);

  delete[] coupling_nameC;
}

void PROCF(cwipi_send_cellvtxfd_eltcont_cf,
           CWIPI_SEND_CELLVTXFD_ELTCONT_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 double *sendingField,
 const int *stride)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_id, *l_coupling_id);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi_exchange_cell_vertex_field_of_element_containing(coupling_nameC,
                                                         sendingField,
                                                         NULL,
                                                         *stride);

  delete[] coupling_nameC;
}

void PROCF(cwipi_recv_cellvtxfd_eltcont_cf,
           CWIPI_RECV_CELLVTXFD_ELTCONT_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 double *receivingField,
 const int *stride)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_id, *l_coupling_id);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi_exchange_cell_vertex_field_of_element_containing(coupling_nameC,
                                                         NULL,
                                                         receivingField,
                                                         *stride);

  delete[] coupling_nameC;
}

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

void PROCF(cwipi_exch_cellcenfd_eltcont_cf,
           CWIPI_EXCH_CELLCENFD_ELTCONT_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 double *sendingField,
 double *receivingField,
 const int *stride)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_id, *l_coupling_id);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi_exchange_cell_center_field_of_element_containing(coupling_nameC,
                                                         sendingField,
                                                         receivingField,
                                                         *stride);
  delete[] coupling_nameC;
}

void PROCF(cwipi_send_cellcenfd_eltcont_cf,
           CWIPI_SEND_CELLCENFD_ELTCONT_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 double *sendingField,
 const int *stride)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_id, *l_coupling_id);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi_exchange_cell_center_field_of_element_containing(coupling_nameC,
                                                         sendingField,
                                                         NULL,
                                                         *stride);
  delete[] coupling_nameC;
}

void PROCF(cwipi_recv_cellcenfd_eltcont_cf,
           CWIPI_RECV_CELLCENFD_ELTCONT_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 double *receivingField,
 const int *stride)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_id, *l_coupling_id);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi_exchange_cell_center_field_of_element_containing(coupling_nameC,
                                                         NULL,
                                                         receivingField,
                                                         *stride);
  delete[] coupling_nameC;
}


/*----------------------------------------------------------------------------
 *
 * Set coupling info
 *
 * parameters
 *   coupling_id          <-- Coupling identifier
 *   info                 <-- Coupling info
 *----------------------------------------------------------------------------*/

void PROCF(cwipi_set_info_cf,
           CWIPI_SET_INFO_CF)
(const char *coupling_id,
 const int  *l_coupling_id,
 const int *info)
{
  char *coupling_nameC =
    _cwipi_fortran_to_c_string(coupling_id, *l_coupling_id);

  const std::string &coupling_name_str = coupling_nameC;

  cwipi_set_info(coupling_nameC, (cwipi_located_point_info_t) *info);
  delete[] coupling_nameC;

}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
