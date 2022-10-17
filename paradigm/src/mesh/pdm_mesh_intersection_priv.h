#ifndef __PDM_MESH_INTERSECTION_PRIV_H__
#define __PDM_MESH_INTERSECTION_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_dbbtree.h"
#include "pdm_surf_mesh.h"
#include "pdm_surf_part.h"
#include "pdm_surf_part_priv.h"
#include "pdm_timer.h"
#include "pdm_overlay.h"
#include "pdm_error.h"
#include "pdm_printf.h"
#include "pdm_part_mesh.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define NTIMER 23

/*=============================================================================
 * Static global variables
 *============================================================================*/

/**
 * \enum _ol_timer_step_t
 *
 */

typedef enum {

  BEGIN                                  = 0,
  REDISTRIBUTE_PTS_HILBERT               = 1,
  BUILD_COARSE_TREE_AND_EXTRACT          = 2,
  BUILD_BBOX_COARSE                      = 3,
  BBOX_COARSE_SOLICITATE                 = 4,
  EQUILIBRATE_WITH_SOLICITATON           = 5,
  EQUILIBRATE_WITH_SOLICITATON_TRANSFERT = 6,
  UPDATE_SOLICITATION_SEND               = 7,
  BUILD_LOCAL_TREE                       = 8,
  BUILD_SHARED_LOCAL_TREE                = 9,
  UPDATE_SOLICITATION_WAIT               = 10,
  LOCAL_SOLICITATE                       = 11,
  EQUILIBRATE_PB                         = 12,
  END                                    = 13

} _mesh_intersection_timer_step_t;


/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _PDM_ol_t
 * \brief  Overlay type
 *
 * _PDM_ol_t defines a overlaying structure
 *
 */
struct _pdm_mesh_intersection_t {

  double   project_coef;         /*!< Projection coefficient to define the overlay
                                      surface projection :
                                      If value == 0, the surface projection is MeshA
                                      If value == 1, the surface projection is MeshB
                                      If 0 < value < 1, the projection surface is an
                                      intermediate surface */
  double   vtx_car_length_tol;   /*!< Absolute tolerance used to define local geometric
                                      tolerance for vertex caracteristic lenght
                                      (tolerance > 0) */

  double   extents_tol;          /*!< Absolute tolerance used to define local geometric
                                      tolerance for vertex caracteristic lenght
                                      (tolerxance > 0) */

  double   same_plane_tol;       /*!< Absolute tolerance used to check if 2 surfaces
                                      are the same plane surface */


  PDM_MPI_Comm comm;             /*!< MPI communicator */

  PDM_mesh_intersection_kind_t intersect_kind;

  int               n_part_mesh_a;
  int               n_part_mesh_b;
  int               dim_mesh_a;
  int               dim_mesh_b;
  PDM_part_mesh_t  *mesh_a;       /*!< Mesh A */
  PDM_part_mesh_t  *mesh_b;       /*!< Mesh B */

  // _ol_mesh_t  *olMeshA;       /*!< Overlay Mesh A */
  // _ol_mesh_t  *olMeshB;       /*!< Overlay Mesh B */

  PDM_timer_t *timer;


  double times_elapsed[NTIMER]; /*!< Elapsed time */

  double times_cpu[NTIMER];     /*!< CPU time */

  double times_cpu_u[NTIMER];  /*!< User CPU time */

  double times_cpu_s[NTIMER];  /*!< System CPU time */

} _PDM_ol_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

// static PDM_Handles_t *olArray = NULL; /*!< Array to storage overlay identifiers */

/*=============================================================================
 * Static function definitions
 *============================================================================*/


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MESH_INTERSECTION_PRIV_H__ */
