
/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <cstring>
#include <cstdlib>
#include <cassert>


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_logging.h"

#include "cwp.h"
#include "cwp_cf.h"
#include "cwipi_config.h"
#include "factory.hpp"
#include "codeProperties.hxx"
#include "codePropertiesDB.hxx"
#include "codePropertiesDB_i.hxx"
#include "couplingDB.hxx"
#include "couplingDB_i.hxx"
#include "coupling.hxx"
#include "coupling_i.hxx"
#include "commWithPart.hxx"
#include "commWithoutPart.hxx"
#include "commSeq.hxx"
#include "field.hxx"
#include "pdm.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "surfMeshGenerator.hxx"
#include "surfMeshGeneratorDB.hxx"
#include "spatialInterpClosestPoint.hxx"
#include "spatialInterpIntersection.hxx"
#include "spatialInterpLocationDistSurf.hxx"
#include "spatialInterpLocationMeshLocation.hxx"

#include "mesh.hxx"
#include "block.hxx"
#include "blockStd.hxx"
#include "blockFP.hxx"
#include "blockCP.hxx"
#include <algorithm>
#include <vector>

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


using namespace std;


/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

// TO DO: change to block_type insted of id
static int
_n_vtx_block_get
(
 int block_id
)
{
  switch (block_id)
  {
    case CWP_BLOCK_EDGE2:
      return 2;
      break;

    case CWP_BLOCK_FACE_TRIA3:
      return 3;
      break;

    case CWP_BLOCK_FACE_QUAD4:
      return 4;
      break;

    case CWP_BLOCK_CELL_TETRA4:
      return 4;
      break;

    case CWP_BLOCK_CELL_HEXA8:
      return 8;
      break;

    case CWP_BLOCK_CELL_PRISM6:
      return 6;
      break;

    case CWP_BLOCK_CELL_PYRAM5:
      return 5;
      break;

    default:
      PDM_error(__FILE__, __LINE__, 0, "Unkown block type %d\n", block_id);
  }
  return 0;
}

static cwipi::Coupling&
_cpl_get
(
 const char *local_code_name,
 const char *cpl_id
 )
{
  cwipi::CouplingDB & couplingDB =
    cwipi::CouplingDB::getInstance();

  cwipi::CodePropertiesDB & properties =
    cwipi::CodePropertiesDB::getInstance();

   const string &cpl_name_str = cpl_id;
   return couplingDB.couplingGet (properties.codePropertiesGet(string(local_code_name)),
                                 cpl_name_str);
}


static char *
_fortran_to_c_string (
  const char *application_name_f, 
  const int l_application_name_f
) 
{
  char *application_name_c;
  int imin = 0;
  int imax = 0;

  while (imin < l_application_name_f && application_name_f[imin] == ' ') {
    imin++;
  }
  while (imax < l_application_name_f && application_name_f[l_application_name_f - imax - 1] == ' ') {
    imax++;
  }

  imax = l_application_name_f - imax - 1;

  assert(imax >= imin);

  if ((imax == l_application_name_f) || (imin == l_application_name_f)) {
    application_name_c = new char[1];
    // application_name_c = (char *) malloc(sizeof(char) * 1);
    application_name_c[0] = '\0';
  }
  else {
    int size = imax - imin + 2;
    application_name_c = new char[size];
    // application_name_c = (char *) malloc(sizeof(char) * size);
    int index = 0;
    for (int k = imin ; k <= imax ; k++) {
      application_name_c[index++] = application_name_f[k];
    }
    application_name_c[index] = '\0';   
  } 

  return application_name_c;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 * \brief Initialize CWIPI.
 *
 * This function creates the MPI intra communicators of the codes from
 * the \p global_comm MPI communicator that contains all code ranks. This
 * function has to be called from all ranks contained in the \p global_comm.
 *
 * \param [in]  global_comm    MPI global communicator
 * \param [in]  n_code         Number of codes on the current rank
 * \param [in]  f_code_names   Names of codes on the current rank (size = \p n_code)
 * \param [in]  l_code_names   Length of code names on the current rank (size = \p n_code)
 * \param [in]  is_active_rank Is current rank have to be used by CWIPI (size = \p n_code)
 * \param [in]  time_init      Initial time (size = \p n_code)
 * \param [out] intra_comms    MPI intra communicators of each code (size = \p n_code)
 *
 */

void 
CWP_Init_cf (
  MPI_Fint f_global_comm, 
  const int n_code, 
  const char *f_code_names,
  const int *l_code_names, 
  const int *is_active_rank, 
  const double *time_init, 
  MPI_Fint *f_intra_comms
) 
{
  // Convert code names dealing with different size characters
  char **c_code_names = (char **) malloc(n_code * sizeof(char *));
  int idx = 0;
  for (int i = 0 ; i < n_code ; i++) {
    // c_code_names[i] = (char *) malloc((l_code_names[i] - 1) * sizeof(char));
    // memcpy(c_code_names[i], *f_code_names + i * l_code_names[i], l_code_names[i]); // TODO This is not perfectly right for the length of the f_code_name which could vary
    c_code_names[i] = _fortran_to_c_string(f_code_names + idx, l_code_names[i]);
    idx += l_code_names[i];
  }

  // Convert global MPI communicator
  MPI_Comm c_global_comm = MPI_Comm_f2c(f_global_comm);

  // Allocate local communicators in C
  MPI_Comm *c_intra_comms = (MPI_Comm *) malloc(n_code * sizeof(MPI_Comm));

  CWP_Init(c_global_comm, n_code, (const char **) c_code_names, (CWP_Status_t *) is_active_rank, time_init, c_intra_comms);

  // Convert local communicators to Fortran
  for (int i = 0 ; i < n_code ; i++) {
    f_intra_comms[i] = MPI_Comm_c2f(c_intra_comms[i]);

    delete [] c_code_names[i];
  }

  // delete [] c_code_names;
  free(c_code_names);
  // delete [] c_intra_comms;
  free(c_intra_comms);
}


/*----------------------------------------------------------------------------*
 * Functions about current code properties                                    *
 *----------------------------------------------------------------------------*/

/**
 * \brief Update code state.
 *
 * \param [in] local_code_name   Fortran local code name
 * \param [in] l_local_code_name Length of Fortran local code name
 * \param [in] state             State
 *
 */

void
CWP_State_update_cf
(
 const char* local_code_name,
 const int l_local_code_name,
 const CWP_State_t state
)
{
  char *c_local_code_name = _fortran_to_c_string(local_code_name, l_local_code_name);

  CWP_State_update (c_local_code_name,
                    state);

  delete [] c_local_code_name;
}


/**
 * \brief Update code time.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] l_local_code_name Length of Fortran local code name
 * \param [in]  current_time Current time
 *
 */

void
CWP_Time_update_cf
(
 const char* local_code_name,
 const int l_local_code_name,
 const double current_time
)
{
  char *c_local_code_name = _fortran_to_c_string(local_code_name, l_local_code_name);

  CWP_Time_update (c_local_code_name, current_time);
  
  delete [] c_local_code_name;
}


/**
 * \brief Define a user structure associated to a code
 *
 * This structure can be called into a callback
 *
 * \param [in] local_code_name  Local code name
 * \param [in] user_structure   User structure
 *
 */

void
CWP_User_structure_set_cf
(
 const char* local_code_name,
 const int   l_local_code_name,
       void* user_structure
)
{
  char *c_local_code_name = _fortran_to_c_string(local_code_name, l_local_code_name);

  CWP_User_structure_set (c_local_code_name, user_structure);

  delete [] c_local_code_name;
}


/**
 * \brief Return the user structure associated
 *
 * This structure can be called into a callback
 *
 * \param [in] local_code_name  Local code name
 *
 * \return  User structure
 *
 */

void *
CWP_User_structure_get_cf
(
 const char* local_code_name,
 const int   l_local_code_name
)
{
  char *c_local_code_name = _fortran_to_c_string(local_code_name, l_local_code_name);

  void *user_structure = CWP_User_structure_get (c_local_code_name);

  delete [] c_local_code_name;

  return user_structure;
}

/*----------------------------------------------------------------------------*
 * Functions about other code properties                               *
 *----------------------------------------------------------------------------*/


/**
 * \brief Return code state.
 *
 * \param [in] local_code_name   Fortran local code name
 * \param [in] l_local_code_name Length of Fortran local code name
 *
 * \return      Code state
 */

CWP_State_t
CWP_State_get_cf
(
 const char    *code_name,
 const int l_local_code_name
)
{
  char *c_local_code_name = _fortran_to_c_string(code_name, l_local_code_name);

  CWP_State_t res = CWP_State_get (c_local_code_name);

  delete [] c_local_code_name;

  return res;

}

/**
 * \brief Return code names and size of those names
 *
 * \param [in] code_list    Code names
 * \param [in] code_list_s  Size of those code name chars
 *
 */

void
CWP_Codes_list_get_cf
(
 const char ***code_list,
 int         **code_list_s,
 int          *n_codes
)
{
  *n_codes   = CWP_Codes_nb_get();
  *code_list = CWP_Codes_list_get();

  *code_list_s = (int *) malloc(sizeof(int) * (*n_codes));

  for (int i = 0; i < (*n_codes); i++) {
    (*code_list_s)[i] = strlen((*code_list)[i]);
  }
}

/**
 * \brief Return local code names and size of those names
 *
 * \param [in] loc_code_list    Code names
 * \param [in] loc_code_list_s  Size of those code name chars
 *
 */

void
CWP_Loc_codes_list_get_cf
(
 const char ***loc_code_list,
 int         **loc_code_list_s,
 int          *n_loc_codes
)
{
  *n_loc_codes = CWP_Loc_codes_nb_get();
  *loc_code_list   = CWP_Loc_codes_list_get();

  *loc_code_list_s = (int *) malloc(sizeof(int) * (*n_loc_codes));

  for (int i = 0; i < (*n_loc_codes); i++) {
    (*loc_code_list_s)[i] = strlen((*loc_code_list)[i]);
  }
}

/*----------------------------------------------------------------------------*
 * General functions about coupling                                           *
 *----------------------------------------------------------------------------*/

/**
 * \brief Create a coupling object and define its properties.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_coupled_code_name Distant or local coupled code name (Fortran)
 * \param [in]  l_coupled_code_name Length of Distant or local coupled code name 
 * \param [in]  comm_type           Communication type
 * \param [in]  spatial_interp      Spatial interpolation method
 * \param [in]  n_part              Number of interface partition
 * \param [in]  displacement        Mesh moving status
 * \param [in]  recv_freq_type      Type of receiving frequency
 *
 */


void 
CWP_Cpl_create_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id,
  const int l_cpl_id, 
  const char *f_coupled_code_name, 
  const int l_coupled_code_name, 
  const CWP_Interface_t entities_dim,
  const CWP_Comm_t comm_type,
  const CWP_Spatial_interp_t spatial_interp, 
  const int n_part, 
  const CWP_Dynamic_mesh_t displacement, 
  const CWP_Time_exch_t freq
) 
{
  char *c_local_code_name, *c_cpl_id, *c_coupled_code_name;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  c_coupled_code_name = _fortran_to_c_string(f_coupled_code_name, l_coupled_code_name);

  CWP_Cpl_create((const char *) c_local_code_name, 
                 (const char *) c_cpl_id, 
                 (const char *) c_coupled_code_name, 
                 entities_dim, comm_type, 
                 spatial_interp, 
                 n_part, 
                 displacement, 
                 freq);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_coupled_code_name;
}

/**
 *
 * \brief Delete a coupling object.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 *
 */

void 
CWP_Cpl_del_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id, 
  const int l_cpl_id
) 
{

  char *c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  char *c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);

  CWP_Cpl_del(c_local_code_name, c_cpl_id);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
}


/**
 *
 * \brief Return the number of uncomputed targets.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_field_id          Fortran Field identifier
 * \param [in]  l_field_id          Length of Fortran Field identifier
 * \param [in] i_part               Current partition
 *
 * \return                Number of uncomputed targets
 */

int 
CWP_N_uncomputed_tgts_get_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id, 
  const int l_cpl_id, 
  const char *f_field_id, 
  const int l_field_id, 
  int i_part
) 
{
 
  char *c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  char *c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  char *c_field_id = _fortran_to_c_string(f_field_id, l_field_id);

  int res = CWP_N_uncomputed_tgts_get(c_local_code_name, c_cpl_id, c_field_id, i_part);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_field_id;

  return res;
}


/**
 *
 * \brief Return uncomputed targets.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_field_id          Fortran Field identifier
 * \param [in]  l_field_id          Length of Fortran Field identifier
 * \param [in]  i_part              Current partition
 *
 * \return                Uncomputed targets
 */

const int *
CWP_Uncomputed_tgts_get_cf (
  const char *f_local_code_name,
  const int l_local_code_name,
  const char *f_cpl_id, 
  const int l_cpl_id, 
  const char *f_field_id, 
  const int l_field_id, 
  int i_part
) 
{

  char *c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  char *c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  char *c_field_id = _fortran_to_c_string(f_field_id, l_field_id);

  const int* res = CWP_Uncomputed_tgts_get(c_local_code_name, c_cpl_id, c_field_id, i_part);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_field_id;

  return res;
}

/**
 *
 * \brief Return the number of computed targets.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_field_id          Fortran Field identifier
 * \param [in]  l_field_id          Length of Fortran Field identifier
 * \param [in]   i_part             Current partition
 *
 * \return                Number of computed targets
 */

int 
CWP_N_computed_tgts_get_cf (
  const char *f_local_code_name, 
  const int l_local_code_name,
  const char *f_cpl_id, 
  const int l_cpl_id, 
  const char *f_field_id, 
  const int l_field_id, 
  int i_part
) 
{
  char *c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  char *c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  char *c_field_id = _fortran_to_c_string(f_field_id, l_field_id);

  int res = CWP_N_computed_tgts_get(c_local_code_name, c_cpl_id, c_field_id, i_part);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_field_id;

  return res;
}

/**
 *
 * \brief Return computed targets.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_field_id          Fortran Field identifier
 * \param [in]  l_field_id          Length of Fortran Field identifier
 * \param [in]  i_part              Current partition
 *
 * \return                Computed targets
 */

const int *
CWP_Computed_tgts_get_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id, 
  const int l_cpl_id, 
  const char *f_field_id, 
  const int l_field_id, 
  int i_part
)
{
  char *c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  char *c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  char *c_field_id = _fortran_to_c_string(f_field_id, l_field_id);

  const int *res = CWP_Computed_tgts_get(c_local_code_name, c_cpl_id, c_field_id, i_part);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_field_id;

  return res;
}

/**
 *
 * \brief Return the number of involved sources</b>
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_field_id          Fortran Field identifier
 * \param [in]  l_field_id          Length of Fortran Field identifier
 * \param [in]   i_part             Current partition
 *
 * \return                Number of involved sources
 */

int
CWP_N_involved_srcs_get_cf (
        const char *f_local_code_name,
        const int l_local_code_name,
        const char *f_cpl_id,
        const int l_cpl_id,
        const char *f_field_id,
        const int l_field_id,
        int i_part
                           )
{
  char *c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  char *c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  char *c_field_id = _fortran_to_c_string(f_field_id, l_field_id);

  int res = CWP_N_involved_srcs_get(c_local_code_name, c_cpl_id, c_field_id, i_part);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_field_id;

  return res;
}

/**
 *
 * \brief Return involved sources</b>
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_field_id          Fortran Field identifier
 * \param [in]  l_field_id          Length of Fortran Field identifier
 * \param [in]   i_part             Current partition
 *
 * \return                Involved sources
 */

const int *
CWP_Involved_srcs_get_cf (
        const char *f_local_code_name,
        const int l_local_code_name,
        const char *f_cpl_id,
        const int l_cpl_id,
        const char *f_field_id,
        const int l_field_id,
        int i_part
                         )
{
  char *c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  char *c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  char *c_field_id = _fortran_to_c_string(f_field_id, l_field_id);

  const int *res = CWP_Involved_srcs_get(c_local_code_name, c_cpl_id, c_field_id, i_part);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_field_id;

  return res;
}

/**
 * \brief Return distance from each target to the source interface. <b>(Not implemented yet)</b>
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 *
 * \return               Distance
 *
 */

const double *
CWP_Computed_tgts_dist_to_spatial_interp_get_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id,
  const int l_cpl_id
) 
{
  char *c_local_code_name, *c_cpl_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);

  const double *res = CWP_Computed_tgts_dist_to_spatial_interp_get(c_local_code_name, c_cpl_id);

  delete [] c_local_code_name;
  delete [] c_cpl_id;

  return res;
}


/**
 * \brief Compute spatial interpolation weights.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 *
 */

void 
CWP_Spatial_interp_weights_compute_cf (
  const char *f_local_code_name, 
  const int l_local_code_name,
  const char *f_cpl_id, 
  const int l_cpl_id
) 
{
  char *c_local_code_name, *c_cpl_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);

  CWP_Spatial_interp_weights_compute(c_local_code_name, c_cpl_id);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
}

/*----------------------------------------------------------------------------*
 * Functions about visualization                                              *
 *----------------------------------------------------------------------------*/

/**
 * \brief Enable visualization output.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  freq                Output frequency
 * \param [in]  format              Output format to visualize exchanged fieldsDouble
 *                                  on the coupled mesh. Choice between :
 *                                  - "EnSight Gold"
 * \param [in]  f_format_option     Fortran Output options "opt1, opt2, ..."
 *                                  - text : output text files
 *                                  - binary : output binary files (default)
 * \param [in]  l_format_option    Length of Fortran option 
 */

void CWP_Visu_set_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id, 
  const int l_cpl_id, 
  int freq, 
  CWP_Visu_format_t format, 
  const char *f_format_option, 
  const int l_format_option
) 
{
  char *c_local_code_name, *c_cpl_id, *c_format_option;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  c_format_option = _fortran_to_c_string(f_format_option, l_format_option);

  CWP_Visu_set(c_local_code_name, c_cpl_id, freq, format, c_format_option);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_format_option;
}

/*----------------------------------------------------------------------------*
 * Functions about User target points                                         *
 *----------------------------------------------------------------------------*/

/**
 * \brief Setting user target points.
 *
 * This function must be called if the degrees of freedom locations are
 * \ref CWP_DOF_LOCATION_USER
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  n_pts            Number of points
 * \param [in]  coord            Coordinates (size = 3 * n_pts)
 * \param [in]  g_num            global number or NUL (size = n_pts)
 *
 */

void CWP_User_tgt_pts_set_cf (
  const char *f_local_code_name,
  int l_local_code_name,
  const char *f_cpl_id,
  int l_cpl_id,
  int i_part,
  int n_pts,
  double coord[],
  CWP_g_num_t global_num[]
)
{
  char *c_local_code_name, *c_cpl_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);

  CWP_User_tgt_pts_set(c_local_code_name, c_cpl_id, i_part, n_pts, coord, global_num);
}


/*----------------------------------------------------------------------------*
 * Functions about Mesh                                                    *
 *----------------------------------------------------------------------------*/


/**
 * \brief Finalize interface mesh.
 *
 * This function computes the global numbers of mesh entities if they are
 * not provided.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 *
 */

void 
CWP_Mesh_interf_finalize_cf
(
  const char *f_local_code_name, 
  int l_local_code_name, 
  const char *f_cpl_id, 
  int l_cpl_id
) 
{
  char *c_local_code_name, *c_cpl_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);

  CWP_Mesh_interf_finalize(c_local_code_name, c_cpl_id);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
}


/**
 * \brief Set vertices.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  i_part              Current partition
 * \param [in]  n_pts               Number of points
 * \param [in]  coord               Coordinates (size = 3 * \p n_pts)
 * \param [in]  global_num          Pointer to parent element number (or NULL)
 *
 */

void 
CWP_Mesh_interf_vtx_set_cf(
  const char *f_local_code_name, 
  int l_local_code_name, 
  const char *f_cpl_id, 
  int l_cpl_id, int i_part, 
  int n_pts, double coord[], 
  CWP_g_num_t global_num[]
) 
{
  char *c_local_code_name, *c_cpl_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);

  CWP_Mesh_interf_vtx_set(c_local_code_name, c_cpl_id, i_part, n_pts, coord, global_num);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
}


/**
 * \brief Add a connectivity block to the interface mesh.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  block_type          Block type
 *
 * \return block identifier
 */

int 
CWP_Mesh_interf_block_add_cf (
  const char *f_local_code_name, 
  int l_local_code_name,
  const char *f_cpl_id, 
  int l_cpl_id, 
  CWP_Block_t block_type
) 
{
  char *c_local_code_name, *c_cpl_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);

  int id = CWP_Mesh_interf_block_add(c_local_code_name, c_cpl_id, block_type);

  delete [] c_local_code_name;
  delete [] c_cpl_id;

  return id;
}



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
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  i_part              Partition identifier
 * \param [in]  block_id            Block identifier
 * \param [in]  n_elts              Number of elements
 * \param [in]  connec              Connectivity (size = n_vertex_elt * n_elts)
 * \param [in]  global_num          Pointer to global element number (or NULL)
 *
 */

void 
CWP_Mesh_interf_block_std_set_cf (
  const char *f_local_code_name, 
  int l_local_code_name, 
  const char *f_cpl_id, 
  int l_cpl_id, 
  int i_part, 
  int block_id,
  int n_elts, 
  int connec[], 
  CWP_g_num_t global_num[]
) 
{
  char *c_local_code_name, *c_cpl_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);

  CWP_Mesh_interf_block_std_set(c_local_code_name, c_cpl_id, i_part, block_id, n_elts, connec, global_num);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
}


/**
 * \brief Get the properties of a standard block of the interface mesh.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Partition identifier
 * \param [in]  block_id         Block identifier
 * \param [out]  n_elts           Number of elements
 * \param [out]  connec           Connectivity (size = n_vertex_elt * n_elts)
 * \param [out]  global_num       Pointer to global element number (or NULL)
 */

void
CWP_Mesh_interf_block_std_get_cf
(
 const char         *f_local_code_name,
       int           l_local_code_name,
 const char         *f_cpl_id,
       int           l_cpl_id,
 const int           i_part,
 const int           block_id,
       int          *n_elts,
       int         **connec,
       CWP_g_num_t **global_num,
       int          *s_connec
)
{
  char *c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  char *c_cpl_id          = _fortran_to_c_string(f_cpl_id,          l_cpl_id);

  CWP_Mesh_interf_block_std_get(c_local_code_name,
                                c_cpl_id,
                                i_part,
                                block_id,
                                n_elts,
                                connec,
                                global_num);

  int n_vtx_block = _n_vtx_block_get(block_id);
  *s_connec = n_vtx_block * (*n_elts);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
}

/**
 * \brief Set the connectivity of a polygon block in a interface mesh partition.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  i_part              Current partition
 * \param [in]  block_id            Block identifier
 * \param [in]  n_elts              Number of elements
 * \param [in]  connec_idx          Connectivity index (\p connec_id[0] = 0 and
 *                                  size = \p n_elts + 1)
 * \param [in]  connec              Connectivity (size = \p connec_idx[\p n_elts])
 * \param [in]  global_num          Pointer to global element number (or NULL)
 *
 */

void 
CWP_Mesh_interf_f_poly_block_set_cf (
  const char *f_local_code_name, 
  int l_local_code_name, 
  const char *f_cpl_id, 
  int l_cpl_id, 
  int i_part, 
  int block_id, 
  int n_elts, 
  int connec_idx[], 
  int connec[], 
  CWP_g_num_t global_num[]
) 
{
  char *c_local_code_name, *c_cpl_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);

  CWP_Mesh_interf_f_poly_block_set(c_local_code_name, c_cpl_id, i_part, block_id, n_elts, connec_idx, connec, global_num);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
}


/**
 * \brief Get the properties of a polygon block of the interface mesh partition.
 *
 * \param [in]  local_code_name  Local code name
 * \param [in]  cpl_id           Coupling identifier
 * \param [in]  i_part           Current partition
 * \param [in]  block_id         Block identifier
 * \param [out]  n_elts           Number of elements
 * \param [out]  connec_idx       Connectivity index (\p connec_id[0] = 0 and
 *                               size = \p n_elts + 1)
 * \param [out]  connec           Connectivity (size = \p connec_idx[\p n_elts])
 * \param [out]  global_num       Pointer to global element number (or NULL)
 *
 */

void
CWP_Mesh_interf_f_poly_block_get_cf
(
 const char         *f_local_code_name,
       int           l_local_code_name,
 const char         *f_cpl_id,
       int           l_cpl_id,
 const int           i_part,
 const int           block_id,
       int          *n_elts,
       int         **connec_idx,
       int         **connec,
       CWP_g_num_t **global_num
)
{
  char *c_local_code_name, *c_cpl_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id          = _fortran_to_c_string(f_cpl_id,          l_cpl_id);

  CWP_Mesh_interf_f_poly_block_get(c_local_code_name, c_cpl_id, i_part, block_id, n_elts, connec_idx, connec, global_num);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
}

/**
 * \brief Adding a polyhedron connectivity block to the interface mesh.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  i_part              Current partition
 * \param [in]  block_id            Block identifier
 * \param [in]  n_elts              Number of elements
 * \param [in]  connec_cells_idx    Polyhedron to face index
 *                                  (\p src_poly_cell_face_idx[0] = 0 and
 *                                   size = \p n_elts + 1)
 * \param [in]  connec_cells        Polyhedron to face connectivity
 *                                  (size = \p cell_face_idx[\p n_elts])
 * \param [in]  n_faces             Number of faces
 * \param [in]  connec_faces_idx    Polyhedron face to vertex index
 *                                  (\p face_vertex_idx[0] = 0 and
 *                                   size = max(\p cell_face_connec) + 1)
 * \param [in]  connec_faces        Polyhedron face to vertex connectivity
 *                                  (size = \p face_vertex_idx[\p n_elts])
 * \param [in]  global_num          Pointer to global element number (or NULL)
 *
 */

void 
CWP_Mesh_interf_c_poly_block_set_cf (
  const char *f_local_code_name,
  int l_local_code_name,
  const char *f_cpl_id, 
  int l_cpl_id, 
  int i_part, 
  int block_id, 
  int n_elts, 
  int n_faces,
  int connec_faces_idx[], 
  int connec_faces[], 
  int connec_cells_idx[], 
  int connec_cells[], 
  CWP_g_num_t global_num[]
) 
{
  char *c_local_code_name, *c_cpl_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);

  CWP_Mesh_interf_c_poly_block_set(c_local_code_name, c_cpl_id, i_part, block_id, n_elts, n_faces, connec_faces_idx, connec_faces, connec_cells_idx, connec_cells, global_num);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
}


/**
 * \brief Get the properties of a polyhedron block of the interface mesh partition..
 *
 * \param [in]  local_code_name   Local code name
 * \param [in]  cpl_id            Coupling identifier
 * \param [in]  i_part            Current partition
 * \param [in]  block_id          Block identifier
 * \param [out]  n_elts            Number of elements
 * \param [out]  connec_cells_idx  Polyhedron to face index
 *                                (\p src_poly_cell_face_idx[0] = 0 and
 *                                 size = \p n_elts + 1)
 * \param [out]  connec_cells      Polyhedron to face connectivity
 *                                (size = \p cell_face_idx[\p n_elts])
 * \param [out]  n_faces           Number of faces
 * \param [out]  connec_faces_idx  Polyhedron face to vertex index
 *                                (\p face_vertex_idx[0] = 0 and
 *                                 size = max(\p cell_face_connec) + 1)
 * \param [out]  connec_faces      Polyhedron face to vertex connectivity
 *                                (size = \p face_vertex_idx[\p n_elts])
 * \param [out]  global_num        Pointer to global element number (or NULL)
 *
 */

void
CWP_Mesh_interf_c_poly_block_get_cf
(
 const char         *f_local_code_name,
       int           l_local_code_name,
 const char         *f_cpl_id,
       int           l_cpl_id,
 const int           i_part,
 const int           block_id,
       int          *n_elts,
       int          *n_faces,
       int         **connec_faces_idx,
       int         **connec_faces,
       int         **connec_cells_idx,
       int         **connec_cells,
       CWP_g_num_t **global_num
)
{
  char *c_local_code_name, *c_cpl_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id          = _fortran_to_c_string(f_cpl_id,          l_cpl_id);

  CWP_Mesh_interf_c_poly_block_get(c_local_code_name, c_cpl_id, i_part, block_id, n_elts, n_faces, connec_faces_idx, connec_faces, connec_cells_idx, connec_cells, global_num);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
}


/**
 * \brief Delete interface mesh.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 *
 */

void
CWP_Mesh_interf_del_cf (
        const char *f_local_code_name,
        const int l_local_code_name,
        const char *f_cpl_id,
        const int l_cpl_id
)
{

  char *c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  char *c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);

  CWP_Mesh_interf_del(c_local_code_name, c_cpl_id);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
}


/**
 * \brief Define the interface mesh from a cell to face connectivity.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  i_part              Current partition
 * \param [in]  n_cells             Number of cells
 * \param [in]  cell_face_idx       Polyhedron to face index
 *                                  (\p src_poly_cell_face_idx[0] = 0 and
 *                                   size = \p n_elts + 1)
 * \param [in]  cell_face           Cell to face connectivity
 *                                  (size = \p cell_face_idx[\p n_elts])
 * \param [in]  n_faces             Number of faces
 * \param [in]  face_vtx_idx        Polyhedron face to vertex index
 *                                  (\p face_vtx_idx[0] = 0 and
 *                                   size = \p n_faces + 1)
 * \param [in]  face_vtx            Face to vertex connectivity
 *                                  (size = \p face_vtx_idx[\p n_elts])
 * \param [in]  parent_num          Pointer to parent element number (or NULL)
 *
 */

void 
CWP_Mesh_interf_from_cellface_set_cf (
  const char *f_local_code_name, 
  int l_local_code_name, 
  const char *f_cpl_id, 
  int l_cpl_id, 
  int i_part, 
  int n_cells,
  int cell_face_idx[], 
  int cell_face[], 
  int n_faces, 
  int face_vtx_idx[], 
  int face_vtx[], 
  CWP_g_num_t parent_num[]
) 
{
  char *c_local_code_name, *c_cpl_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);

  CWP_Mesh_interf_from_cellface_set(c_local_code_name, c_cpl_id, i_part, n_cells, cell_face_idx, cell_face, n_faces, face_vtx_idx, face_vtx, parent_num);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
}


/**
 * \brief Define the surface interface mesh from a face to edge connectivity.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  i_part              Current partition
 * \param [in]  n_faces             Number of cells
 * \param [in]  face_edge_idx       Polygon to edge index
 *                                  (\p face_edge_idx[0] = 0 and
 *                                   size =  \p n_faces + 1)
 * \param [in]  face_edge           Face to edge connectivity
 *                                  (size = \p face_edge_idx[\p n_faces])
 * \param [in]  n_edges             Number of faces
 * \param [in]  edge_vtx_idx        Polyhedron face to vertex index
 *                                  (\p edge_vtx_idx[0] = 0 and
 *                                   size = \p n_edges + 1)
 * \param [in]  edge_vtx            Face to vertex connectivity
 *                                  (size = \p edge_vtx_idx[\p n_edges])
 * \param [in]  parent_num          Pointer to parent element number (or NULL)
 *
 */

void 
CWP_Mesh_interf_from_faceedge_set_cf (
  const char *f_local_code_name, 
  int l_local_code_name, 
  const char *f_cpl_id, 
  int l_cpl_id, 
  int i_part, 
  int n_faces,
  int face_edge_idx[], 
  int face_edge[], 
  int n_edges, 
  int edge_vtx_idx[], 
  int edge_vtx[], 
  CWP_g_num_t parent_num[]
) 
{
  char *c_local_code_name, *c_cpl_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);

  CWP_Mesh_interf_from_faceedge_set(c_local_code_name, 
                                    c_cpl_id, 
                                    i_part,
                                    n_faces, 
                                    face_edge_idx, 
                                    face_edge, 
                                    n_edges, 
                                    edge_vtx_idx, 
                                    edge_vtx, 
                                    parent_num);

  delete [] c_local_code_name;
  delete [] c_cpl_id;

}

/*----------------------------------------------------------------------------*
 * Functions about field                                                      *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Create a new field.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_field_id          Fortran Field identifier
 * \param [in]  l_field_id          Length of Fortran Field identifier
 * \param [in]  data_type           Data type
 * \param [in]  storage             Storage type
 * \param [in]  n_component         Number of component
 * \param [in]  target_location     Target location
 * \param [in]  exch_type           Exchange type
 * \param [in]  visu_status         Visualization status
 *
 */

void 
CWP_Field_create_cf
(
  const char *f_local_code_name, 
  int l_local_code_name, 
  const char *f_cpl_id, 
  int l_cpl_id, 
  const char *f_field_id, 
  int l_field_id, 
  CWP_Type_t data_type, 
  CWP_Field_storage_t storage, 
  int n_component, 
  CWP_Dof_location_t target_location, 
  CWP_Field_exch_t exch_type, 
  CWP_Status_t visu_status
) 
{
  char *c_local_code_name, *c_cpl_id, *c_field_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  c_field_id = _fortran_to_c_string(f_field_id, l_field_id);

  CWP_Field_create(c_local_code_name, 
                   c_cpl_id, 
                   c_field_id, 
                   data_type,
                   storage, 
                   n_component, 
                   target_location, 
                   exch_type, 
                   visu_status);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_field_id;

}


/**
 *
 * \brief Set field data.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  i_part              Current partition
 * \param [in]  data                Storage array (Mapping)
 *
 */

void 
CWP_Field_data_set_cf (
  const char *f_local_code_name, 
  int l_local_code_name, 
  const char *f_cpl_id, 
  int l_cpl_id, 
  const char *f_field_id, 
  int l_field_id, 
  int i_part, 
  int map_type,
  double data[]
) 
{
  char *c_local_code_name, *c_cpl_id, *c_field_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  c_field_id = _fortran_to_c_string(f_field_id, l_field_id);

  CWP_Field_data_set(c_local_code_name, c_cpl_id, c_field_id, i_part, (CWP_Field_map_t) map_type, data);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_field_id;

}


/**
 *
 * \brief Get number of field components.
 *  * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      number of field components
 *
 */

int
CWP_Field_n_component_get_cf
(
 const char      *f_local_code_name,
       int        l_local_code_name,
 const char      *f_cpl_id,
       int        l_cpl_id,
 const char      *f_field_id,
       int        l_field_id
)
{
  char *c_local_code_name, *c_cpl_id, *c_field_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id          = _fortran_to_c_string(f_cpl_id,          l_cpl_id);
  c_field_id        = _fortran_to_c_string(f_field_id,        l_field_id);

  int n_component = CWP_Field_n_component_get(c_local_code_name, c_cpl_id, c_field_id);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_field_id;

  return n_component;
}


/**
 *
 * \brief Get target degrees of freedom location.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] field_id         Field identifier
 *
 * \return                      Location of degrees of freedom
 *
 */

CWP_Dof_location_t
CWP_Field_target_dof_location_get_cf
(
 const char      *f_local_code_name,
       int        l_local_code_name,
 const char      *f_cpl_id,
       int        l_cpl_id,
 const char      *f_field_id,
       int        l_field_id
)
{
  char *c_local_code_name, *c_cpl_id, *c_field_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id          = _fortran_to_c_string(f_cpl_id,          l_cpl_id);
  c_field_id        = _fortran_to_c_string(f_field_id,        l_field_id);

  CWP_Dof_location_t dof_location = CWP_Field_target_dof_location_get(c_local_code_name, c_cpl_id, c_field_id);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_field_id;

  return dof_location;
}


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
CWP_Field_storage_get_cf
(
 const char      *f_local_code_name,
       int        l_local_code_name,
 const char      *f_cpl_id,
       int        l_cpl_id,
 const char      *f_field_id,
       int        l_field_id
)
{
  char *c_local_code_name, *c_cpl_id, *c_field_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id          = _fortran_to_c_string(f_cpl_id,          l_cpl_id);
  c_field_id        = _fortran_to_c_string(f_field_id,        l_field_id);

  CWP_Field_storage_t storage = CWP_Field_storage_get(c_local_code_name, c_cpl_id, c_field_id);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_field_id;

  return storage;
}


/**
 * \brief Delete a field.
 *
 * \param [in] local_code_name Local code name
 * \param [in]  cpl_id         Coupling identifier
 * \param [in]  field_id       Field identifier
 *
 */

void
CWP_Field_del_cf
(
 const char      *f_local_code_name,
       int        l_local_code_name,
 const char      *f_cpl_id,
       int        l_cpl_id,
 const char      *f_field_id,
       int        l_field_id
)
{
  char *c_local_code_name, *c_cpl_id, *c_field_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id          = _fortran_to_c_string(f_cpl_id,          l_cpl_id);
  c_field_id        = _fortran_to_c_string(f_field_id,        l_field_id);

  CWP_Field_del(c_local_code_name, c_cpl_id, c_field_id);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_field_id;
}

/*----------------------------------------------------------------------------*
 * Functions about exchange                                                   *
 *----------------------------------------------------------------------------*/

/**
 * \brief Send a spatially interpolated field to the coupled code with
 *        nonblocking communications.
 *
 * This function is independant of \ref CWP_Time_exch_t mode. The user has to
 * manually check the consistency of the exchanges.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_src_field_id      Fortran Source field id
 * \param [in]  l_src_field_id      Length of Source field id
 *
 *    
 */

void 
CWP_Field_issend_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id, 
  const int l_cpl_id, 
  const char *f_src_field_id, 
  const int l_src_field_id
) 
{
  char *c_local_code_name, *c_cpl_id, *c_src_field_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  c_src_field_id = _fortran_to_c_string(f_src_field_id, l_src_field_id);

  CWP_Field_issend(c_local_code_name, c_cpl_id, c_src_field_id);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_src_field_id;
}

/**
 *
 * \brief Receive a spatially interpolated field from the coupled code
 *        with nonblocking communications.
 *
 * This function is independant of \ref CWP_Time_exch_t mode. The user has to
 * manually check the consistency of the exchanges.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_tgt_field_id      Fortran Target field id
 * \param [in]  f_tgt_field_id      Length ofFortran Target field id
 *
 *
 */

void 
CWP_Field_irecv_cf(
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id, 
  const int l_cpl_id, 
  const char *f_tgt_field_id, 
  const int l_tgt_field_id
) 
{
  char *c_local_code_name, *c_cpl_id, *c_tgt_field_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  c_tgt_field_id = _fortran_to_c_string(f_tgt_field_id, l_tgt_field_id);

  CWP_Field_irecv(c_local_code_name, c_cpl_id, c_tgt_field_id);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_tgt_field_id;
}

/**
 *
 * \brief Wait the end of an exchange related to request from \ref CWP_Field_issend.
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_src_field_id      Fortran Source field id
 * \param [in]  l_src_field_id      Length of Source field id
 *
 */

void 
CWP_Field_wait_issend_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id,
  const int l_cpl_id, 
  const char *f_src_field_id, 
  const int l_src_field_id
) 
{
  char *c_local_code_name, *c_cpl_id, *c_src_field_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  c_src_field_id = _fortran_to_c_string(f_src_field_id, l_src_field_id);

  CWP_Field_wait_issend(c_local_code_name, c_cpl_id, c_src_field_id);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_src_field_id;
}

/**
 *
 * \brief Wait the end of an exchange related to request from \ref CWP_Field_irecv.
 *
 * This function waits the end of exchange related to request
 * from \ref CWP_Field_irecv
 *
 * \param [in]  f_local_code_name   Fortran local code name
 * \param [in]  l_local_code_name   Length of Fortran local code name
 * \param [in]  f_cpl_id            Fortran Coupling identifier
 * \param [in]  l_cpl_id            Length of Fortran coupling identifier
 * \param [in]  f_tgt_field_id      Fortran Target field id
 * \param [in]  f_tgt_field_id      Length ofFortran Target field id
 *
 */

void 
CWP_Field_wait_irecv_cf (
  const char *f_local_code_name, 
  const int l_local_code_name, 
  const char *f_cpl_id, 
  const int l_cpl_id, 
  const char *f_tgt_field_id, 
  const int l_tgt_field_id
) 
{
  char *c_local_code_name, *c_cpl_id, *c_tgt_field_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  c_tgt_field_id = _fortran_to_c_string(f_tgt_field_id, l_tgt_field_id);

  CWP_Field_wait_irecv(c_local_code_name, c_cpl_id, c_tgt_field_id);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_tgt_field_id;
}

/**
 *
 * \brief Setting of an user interpolation from location.
 *
 * This function takes into account an user interpolation function written with
 * void (*\ref CWP_Interp_from_location_t) interface.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] src_field_id     Source field id
 * \param [in] fct              Function
 *
 */

void
CWP_Interp_from_location_set_cf (
  const char *f_local_code_name,
  int l_local_code_name,
  const char *f_cpl_id,
  int l_cpl_id,
  const char *f_src_field_id,
  int l_src_field_id,
  void* ptInterpolationFct
)
{
  char *c_local_code_name, *c_cpl_id, *c_src_field_id;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  c_src_field_id = _fortran_to_c_string(f_src_field_id, l_src_field_id);

  CWP_Interp_from_location_set(c_local_code_name, c_cpl_id, c_src_field_id, (CWP_Interp_from_location_t) ptInterpolationFct);

  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_src_field_id;
}



void
CWP_Spatial_interp_property_set_cf
(
 const char     *f_local_code_name,
       int       l_local_code_name,
 const char     *f_cpl_id,
       int       l_cpl_id,
 const char     *f_property_name,
       int       l_property_name,
 const char     *f_property_type,
       int       l_property_type,
 const char     *f_property_value,
       int       l_property_value
)
{
  char *c_local_code_name, *c_cpl_id, *c_property_name, *c_property_type, *c_property_value;

  c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  c_cpl_id          = _fortran_to_c_string(f_cpl_id,          l_cpl_id);
  c_property_name   = _fortran_to_c_string(f_property_name,   l_property_name);
  c_property_type   = _fortran_to_c_string(f_property_type,   l_property_type);
  c_property_value  = _fortran_to_c_string(f_property_value,  l_property_value);

  CWP_Spatial_interp_property_set(c_local_code_name,
                                  c_cpl_id,
                                  c_property_name,
                                  c_property_type,
                                  c_property_value);
  delete [] c_local_code_name;
  delete [] c_cpl_id;
  delete [] c_property_name;
  delete [] c_property_type;
  delete [] c_property_value;
}


/*----------------------------------------------------------------------------*
 * Functions about control parameters                                         *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Add a new parameter and intialize it.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 * \param [in] initial_value    Initial value
 *
 */

void
CWP_Param_add_cf
(
 const char        *f_local_code_name,
 const int          l_local_code_name,
 const char        *f_param_name,
 const int          l_param_name,
 const CWP_Type_t   data_type,
 void              *initial_value
)
{
  char *c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  char *c_param_name      = _fortran_to_c_string(f_param_name,      l_param_name     );

  CWP_Param_add(c_local_code_name,
                c_param_name,
                data_type,
                initial_value);

  delete [] c_local_code_name;
  delete [] c_param_name;
}


/**
 *
 * \brief Set a parameter.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 * \param [in] value            Value
 *
 */

void
CWP_Param_set_cf
(
 const char        *f_local_code_name,
 const int          l_local_code_name,
 const char        *f_param_name,
 const int          l_param_name,
 const CWP_Type_t   data_type,
 void              *value
)
{
  char *c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  char *c_param_name      = _fortran_to_c_string(f_param_name,      l_param_name     );

  CWP_Param_set(c_local_code_name,
                c_param_name,
                data_type,
                value);

  delete [] c_local_code_name;
  delete [] c_param_name;
}

/**
 *
 * \brief Delete a parameter.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] param_name       Parameter name
 * \param [in] data_type        Parameter type
 *
 */

void
CWP_Param_del_cf
(
 const char       *f_local_code_name,
 const int         l_local_code_name,
 const char       *f_param_name,
 const int         l_param_name,
 const CWP_Type_t  data_type
)
{
  char *c_local_code_name = _fortran_to_c_string(f_local_code_name, l_local_code_name);
  char *c_param_name      = _fortran_to_c_string(f_param_name,      l_param_name     );

  CWP_Param_del(c_local_code_name,
                c_param_name,
                data_type);

  delete [] c_local_code_name;
  delete [] c_param_name;
}


/*----------------------------------------------------------------------------*
 * Functions about all code parameters                                        *
 *----------------------------------------------------------------------------*/

/**
 *
 * \brief Return the number of parameters for the code \p code_name.
 *
 * \param [in] code_name       Local or distant code name
 * \param [in] data_type       Parameter type,
 *
 * return  Number of parameters
 *
 */

int
CWP_Param_n_get_cf
(
 const char             *f_code_name,
 const int               l_code_name,
 const CWP_Type_t        data_type
)
{
  char *c_code_name = _fortran_to_c_string(f_code_name, l_code_name);

  int n_param = CWP_Param_n_get(c_code_name,
                                data_type);

  delete [] c_code_name;

  return n_param;
}


// /**
//  *
//  * \brief Return the list of parameters for the code \p code_name.
//  *
//  * \param [in]  code_name      Local or distant code name
//  * \param [in]  data_type      Parameter type,
//  * \param [out] nParam         Number of parameters
//  * \param [out] paramNames     Parameter names
//  *
//  *
//  */
//
// void
// CWP_Param_list_get
// (
//  const char             *code_name,
//  const CWP_Type_t        data_type,
//  int                    *nParam,
//  char                 ***paramNames
// )


/**
 *
 * \brief Is this \p code_name a parameter ?
 *
 * \param [in] code_name      Local or distant code name
 * \param [in] param_name     Parameter name
 * \param [in] data_type      Parameter type
 *
 * return  1 : true / 0 : false
 *
 */

int
CWP_Param_is_cf
(
 const char            *f_code_name,
 const int              l_code_name,
 const char            *f_param_name,
 const int              l_param_name,
 const CWP_Type_t       data_type
)
{
  char *c_code_name  = _fortran_to_c_string(f_code_name,  l_code_name);
  char *c_param_name = _fortran_to_c_string(f_param_name, l_param_name);

  int is_param = CWP_Param_is(c_code_name,
                              c_param_name,
                              data_type);

  delete [] c_code_name;
  delete [] c_param_name;

  return is_param;
}


/**
 *
 * \brief Return the parameter value of \p param_name on \p code_name.
 *
 * \param [in]  code_name  Local or distant code name
 * \param [in]  param_name Parameter name
 * \param [in]  data_type  Parameter type
 * \param [out] value      Parameter value
 *
 */

void
CWP_Param_get_cf
(
 const char       *f_code_name,
 const int         l_code_name,
 const char       *f_param_name,
 const int         l_param_name,
 const CWP_Type_t  data_type,
 void             *value
)
{
  char *c_code_name  = _fortran_to_c_string(f_code_name,  l_code_name);
  char *c_param_name = _fortran_to_c_string(f_param_name, l_param_name);

  CWP_Param_get(c_code_name,
                c_param_name,
                data_type,
                value);

  delete [] c_code_name;
  delete [] c_param_name;
}


// /**
//  *
//  * \brief Return the result of a reduce operation about a parameter.
//  *
//  * The parameter name has to be the same for all codes.
//  *
//  * \param [in]  op           Operation
//  * \param [in]  param_name   Parameter name
//  * \param [in]  data_type    Parameter type,
//  * \param [out] res          Result
//  * \param [in]  nCode        Number of codes
//  * \param       ...          Codes name
//  *
//  */
//
// void
// CWP_Param_reduce
// (
//  const CWP_Op_t    op,
//  const char       *param_name,
//  const CWP_Type_t  data_type,
//  void             *res,
//  const int         nCode,
//  ...
// );


/**
 *
 * \brief Lock access to local parameters from a distant code.
 *
 * \param [in]  code_name  Code to lock
 *
 */

void
CWP_Param_lock_cf
(
 const char *f_code_name,
 const int   l_code_name
 )
{
  char *c_code_name = _fortran_to_c_string(f_code_name, l_code_name);

  CWP_Param_lock(c_code_name);

  delete [] c_code_name;
}


/**
 *
 * \brief Unlock access to local parameters from a distant code.
 *
 * \param [in]  code_name  Code to unlock
 *
 */

void
CWP_Param_unlock_cf
(
 const char *f_code_name,
 const int   l_code_name
)
{
  char *c_code_name = _fortran_to_c_string(f_code_name, l_code_name);

  CWP_Param_unlock(c_code_name);

  delete [] c_code_name;
}

/**
 *
 * \brief Setting of an fortran user interpolation from location.
 *
 * This function takes into account a fortran user interpolation function written with
 * void (*\ref CWP_Interp_from_location_t) interface.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] src_field_id     Source field id
 * \param [in] fct              Function
 *
 */

void
CWP_Interp_from_location_fortran_set_cf
(
 const char *f_code_name,
 const int   l_code_name,
 const char *f_cpl_id,
 const int   l_cpl_id, 
 const char *f_src_field_id, 
 const int   l_src_field_id,
 void*       fct
)
{
  char *c_code_name    = _fortran_to_c_string(f_code_name, l_code_name);
  char *c_cpl_id       = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  char *c_src_field_id = _fortran_to_c_string(f_src_field_id, l_src_field_id);

  const string &field_name_str = c_src_field_id;

  cwipi::Coupling& cpl = _cpl_get(c_code_name, c_cpl_id);

  cpl.interpFortranFromLocSet(field_name_str, fct);

  delete [] c_code_name;
  delete [] c_cpl_id;
  delete [] c_src_field_id;

}

/**
 *
 * \brief Setting of an fortran user interpolation from location.
 *
 * This function takes into account a fortran user interpolation function written with
 * void (*\ref CWP_Interp_from_location_t) interface.
 *
 * \param [in] local_code_name  Local code name
 * \param [in] cpl_id           Coupling identifier
 * \param [in] src_field_id     Source field id
 * \param [in] fct              Function
 *
 */

void
CWP_Interp_from_location_fortran_unset_cf
(
 const char *f_code_name,
 const int   l_code_name,
 const char *f_cpl_id,
 const int   l_cpl_id, 
 const char *f_src_field_id, 
 const int   l_src_field_id
)
{
  char *c_code_name    = _fortran_to_c_string(f_code_name, l_code_name);
  char *c_cpl_id       = _fortran_to_c_string(f_cpl_id, l_cpl_id);
  char *c_src_field_id = _fortran_to_c_string(f_src_field_id, l_src_field_id);

  const string &field_name_str = c_src_field_id;

  cwipi::Coupling& cpl = _cpl_get(c_code_name, c_cpl_id);

  cpl.interpFortranFromLocUnSet(field_name_str);

  delete [] c_code_name;
  delete [] c_cpl_id;
  delete [] c_src_field_id;

}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
