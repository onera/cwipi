#ifndef __PDM_INTERPOLATE_FROM_MESH_LOCATION_PRIV_H__
#define __PDM_INTERPOLATE_FROM_MESH_LOCATION_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_location_priv.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef  __cplusplus
extern "C" {
#endif

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct _pdm_interpolate_from_mesh_location_t
 *
 * \brief  Data transfer from blocks to partitions
 *
 */

struct _pdm_interpolate_from_mesh_location_t {
  PDM_MPI_Comm   comm;                /*!< MSG communicator */
  int            n_rank;              /*!< Communicator size */
  int            i_rank;              /*!< Current rank in comm */

  int            n_part_src;          /*!< Number of partitions */
  int           *n_cell;              /*!< Number of elements for each partition */
  int           *n_face;              /*!< Number of elements for each partition */
  int           *n_vtx;               /*!< Number of elements for each partition */

  const int          **cell_face_idx;
  const int          **cell_face;
  const PDM_g_num_t  **cell_ln_to_gn;

  const int          **face_vtx_idx;
  const int          **face_vtx;
  const PDM_g_num_t  **face_ln_to_gn;

  const PDM_g_num_t  **vtx_ln_to_gn;
  const double       **coords;

  int            n_cloud_target;

  _points_in_element_t *points_in_elements; /*!< System CPU time */


} _pdm_interpolate_from_mesh_location_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

#ifdef  __cplusplus
}
#endif

#endif  /* __PDM_INTERPOLATE_FROM_MESH_LOCATION_PRIV_H__ */
