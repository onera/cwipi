#ifndef __PDM_TRIANGLE_H__
#define __PDM_TRIANGLE_H__

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
 * \enum PDM_triangle_status_t
 * \brief Triangle status type
 *
 */

typedef enum {

  PDM_TRIANGLE_INSIDE      = 0,  /*!< Inside  */
  PDM_TRIANGLE_OUTSIDE     = 1,  /*!< Outside */
  PDM_TRIANGLE_DEGENERATED = 2,  /*!< Degenerated */

} PDM_triangle_status_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 * \brief Evaluates the position in a triangle
 *
 * \param [in]  x         Point coordinates to evaluate position
 * \param [in]  pts       Triangle vertices coordinates
 * \param [out] closest   Closest Point in Triangle or NULL
 * \param [out] min_dist2 Square of the distance
 * \param [out] weights   Vertices weights or NULL
 *
 * \return      \ref PDM_TRIANGLE_INSIDE or \ref PDM_TRIANGLE_OUTSIDE
 *              if the projected is in the triangle or not
 *
 */

/*  This function is derived from VTK                                      */
/*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen                */
/*  All rights reserved.                                                   */
/*  See Copyright.txt or http://www.kitware.com/Copyright.htm for details. */

PDM_triangle_status_t
PDM_triangle_evaluate_position
(
 const double  x[3],
 const double  pts[9],
       double *closestPoint,
       double *min_dist2,
       double *weights
 );


/**
 * \brief Computes triangle barycenter
 *
 * \param [in]   pts     Triangle vertices coordinates
 * \param [out]  bary    Barycenter
 *
 */

void
PDM_triangle_compute_barycenter
(
 const double pts[9],
       double bary[3]
);



/**
 * \brief Computes closest point on a triangle
 *
 * \param [in]   p              Point coordinates
 * \param [in]   v              Triangle vertices coordinates
 * \param [out]  closest_point  Closest point coordinates
 * \param [out]  min_dist2      Squared distance from point to triangle
 * \param [out]  weights        Barycentric coordinates of closest point
 *
 */

PDM_triangle_status_t
PDM_triangle_closest_point
(
 const double  p[3],
 const double  v[9],
 double       *closest_point,
 double       *min_dist2,
 double       *weights
 );


/**
 * \brief Computes the center and radius of a triangle's circumcircle
 *
 * \param [in]   vtx_coord  Triangle vertices coordinates
 * \param [out]  center     Circumcircle center
 * \param [out]  radius     Circumcircle radius
 *
 */

void
PDM_triangle_circumcircle
(
 const double  vtx_coord[9],
 double        center[3],
 double       *radius
 );


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_SURF_MESH_H__ */
