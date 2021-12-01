
/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <cstring>
#include <cstdlib>
#include <cassert>


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/


#include "cwp.h"
#include "cwp_cf.h"

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
    application_name_c[0] = '\0';
  }
  else {
    int size = imax - imin + 2;
    application_name_c = new char[size];
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
  const char **f_code_names, 
  const int *l_code_names, 
  const int *is_active_rank, 
  const double *time_init, 
  MPI_Fint *f_intra_comms
) 
{
  // Convert code names dealing with different size characters
  char **c_code_names = (char **) malloc(n_code * sizeof(char *));
  for (int i = 0 ; i < n_code ; i++) {
    c_code_names[i] = (char *) malloc((l_code_names[i] - 1) * sizeof(char));
    memcpy(c_code_names[i], *f_code_names + i * l_code_names[i], l_code_names[i]); // TODO This is not perfectly right for the length of the f_code_name which could vary
  }

  // Convert global MPI communicator
  MPI_Comm c_global_comm = MPI_Comm_f2c(f_global_comm);

  // Allocate local communicators in C
  MPI_Comm *c_intra_comms = (MPI_Comm *) malloc(n_code * sizeof(MPI_Comm));

  CWP_Init(c_global_comm, n_code, (const char **) c_code_names, (CWP_Status_t *) is_active_rank, time_init, c_intra_comms);

  // Convert local communicators to Fortran
  for (int i = 0 ; i < n_code ; i++) {
    f_intra_comms[i] = MPI_Comm_c2f(c_intra_comms[i]);
  }

  delete [] c_code_names;
  delete [] c_intra_comms;
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
 * \brief Return the number of computed targets. <b>(Not implemented yet)</b>
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
 * \brief Return computed targets. <b>(Not implemented yet)</b>
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

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
