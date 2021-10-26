#ifndef __PDM_TETRAHEDRON_H__
#define __PDM_TETRAHEDRON_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"


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
 * Type
 *============================================================================*/

/**
 * \enum PDM_tetrahedron_status_t
 * \brief Tetrahedron status type
 *
 */

typedef enum {

  PDM_TETRAHEDRON_INSIDE      = 0,  /*!< Inside  */
  PDM_TETRAHEDRON_OUTSIDE     = 1,  /*!< Outside */
  PDM_TETRAHEDRON_DEGENERATED = 2,  /*!< Degenerated */

} PDM_tetrahedron_status_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 * \brief Evaluates the position in a tetrahedron
 *
 * \param [in]  x                         Point coordinates to evaluate position
 * \param [in]  vtx_coord                 Tetrahedron vertices coordinates
 * \param [out] closest_point             Closest Point in Tetrahedron or NULL
 * \param [out] closest_point_dist2       Square of the distance
 * \param [out] closest_point_weights     Vertices weights or NULL
 *
 * \return      -1 if the tetrahedron is degenerate, 0 else
 *
 */

PDM_tetrahedron_status_t
PDM_tetrahedron_evaluate_position
(
 const double  x[3],
 const double  vtx_coord[12],
 double        closest_point[3],
 double       *closest_point_dist2,
 double        closest_point_weights[4]
 );


/**
 * \brief Computes tetrahedron barycenter
 *
 * \param [in]   pts     Tetrahedron vertices coordinates
 * \param [out]  bary    Barycenter
 *
 */

void
PDM_tetrahedron_compute_barycenter
(
 const double pts[12],
       double bary[3]
 );


void
PDM_tetrahedron_circumsphere
(
 const double  vtx_coord[12],
 double        center[3],
 double       *radius
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_TETRAHEDRON_H__ */
