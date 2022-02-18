/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <sys/resource.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_timer.h"
#include "pdm_sort.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"

#include "pdm_point_cloud_gen.h"

#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/


/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Generate a uniformly random point cloud inside a cuboid.
 *
 * \param [in]   comm                   MPI Communicator id
 * \param [in]   gn_pts                 Global number of points in the cloud
 * \param [in]   x_min                  X-coordinate of the first cuboid corner
 * \param [in]   y_min                  Y-coordinate of the first cuboid corner
 * \param [in]   z_min                  Z-coordinate of the first cuboid corner
 * \param [in]   x_max                  X-coordinate of the opposite cuboid corner
 * \param [in]   y_max                  Y-coordinate of the opposite cuboid corner
 * \param [in]   z_max                  Z-coordinate of the opposite cuboid corner
 * \param [out]  ln_pts                 Local number of points in the cloud
 * \param [out]  coord                  XYZ-coordinates of the local points
 * \param [out]  g_num                  Global ids of the local points
 *
 */

void
PDM_point_cloud_gen_random
(
 PDM_MPI_Comm        comm,
 const PDM_g_num_t   gn_pts,
 const double        x_min,
 const double        y_min,
 const double        z_min,
 const double        x_max,
 const double        y_max,
 const double        z_max,
 int                *ln_pts,
 double            **coord,
 PDM_g_num_t       **g_num
 )
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_g_num_t *distrib_pts = PDM_compute_uniform_entity_distribution(comm, gn_pts);

  *ln_pts = (int) (distrib_pts[i_rank+1] - distrib_pts[i_rank]);


  /**
   * Coordinates
   */

  *coord = malloc (sizeof(double) * (*ln_pts) * 3);

  if (*coord == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Failed to allocate coords (size = %d * 3 * sizeof(double))\n", *ln_pts);
  }

  double origin[3] = {x_min, y_min, z_min};

  double length[3] = {
    x_max - x_min,
    y_max - y_min,
    z_max - z_min
  };

  double i_rand_max = 1. / ((double) RAND_MAX);

  for (int i = 0; i < *ln_pts; i++) {

    unsigned int seed = (unsigned int) (distrib_pts[i_rank] + i);
    srand(seed);

    for (int j = 0; j < 3; j++) {
      (*coord)[3*i + j] = origin[j] + length[j] * (double) rand() * i_rand_max;
    }
  }


  /**
   * Global numbers
   */
  if (1) {
    double _char_length = 1e-6 * PDM_MAX(length[0], PDM_MAX(length[1], length[2]));

    double *char_length = malloc(sizeof(double) * (*ln_pts));

    for (int i = 0; i < *ln_pts; i++) {
      char_length[i] = _char_length;
    }

    PDM_gen_gnum_t* gen_gnum = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, comm, PDM_OWNERSHIP_USER);

    PDM_gnum_set_from_coords (gen_gnum, 0, *ln_pts, *coord, char_length);

    PDM_gnum_compute (gen_gnum);

    *g_num = PDM_gnum_get (gen_gnum, 0);

    PDM_gnum_free (gen_gnum);
    free (char_length);
  }

  else {
    for (int i = 0; i < *ln_pts; i++) {
      (*g_num)[i] = distrib_pts[i_rank] + i;
    }
  }

  free (distrib_pts);
}

