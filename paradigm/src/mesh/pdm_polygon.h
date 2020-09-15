#ifndef __PDM_POLYGON_H__
#define __PDM_POLYGON_H__

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
 * \enum PDM_polygon_status_t
 * \brief Polygon status type
 *
 */

typedef enum {

  PDM_POLYGON_INSIDE      = 0,  /*!< Inside  */
  PDM_POLYGON_OUTSIDE     = 1,  /*!< Outside */
  PDM_POLYGON_DEGENERATED = 2,  /*!< Degenerated */

} PDM_polygon_status_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 * \brief Get bounds
 *
 * \param [in]  numPts   Number of polygon vertices
 * \param [in]  pts      Polygon vertices coordinates
 *
 * \return      Bounds
 *
 */


double *
PDM_bounds_get
(
 const int     numPts,
 const double *pts
);


/**
 * \brief Evaluates the position in a polygon
 *
 * \param [in]  x        Point coordinates to evaluate position
 * \param [in]  numPts   Number of polygon vertices
 * \param [in]  pts      Polygon vertices coordinates
 * \param [out] closest  Closest Point in Polygon or NULL
 * \param [out] minDist2 Square of the distance
 *
 * \return      \ref PDM_POLYGON_INSIDE or \ref PDM_POLYGON_OUTSIDE
 *              if the projected is in the polygon or not
 *
 */

/*  This function is derived from VTK                                      */
/*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen                */
/*  All rights reserved.                                                   */
/*  See Copyright.txt or http://www.kitware.com/Copyright.htm for details. */

PDM_polygon_status_t
PDM_polygon_evaluate_position
(
 const double  x[3],
 const int     numPts,
 const double *pts,
 double        closestPoint[3],
 double        *minDist2
);


/**
 * \brief Computes polygon parametrization
 *
 * \param [in]  numPts  Number of polygon vertices
 * \param [in]  pts     Polygon vertices coordinates
 * \param [out] p0,     Origin vertex
 * \param [out] p10,    First edge vector
 * \param [out] l10,    First edge length
 * \param [out] p20,    First edge vector
 * \param [out] l20,    Second edge vector
 * \param [out] n       Normal
 *
 * \return      \ref PDM_TRUE except for a triangle
 *
 */

/*  This function is derived from VTK                                      */
/*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen                */
/*  All rights reserved.                                                   */
/*  See Copyright.txt or http://www.kitware.com/Copyright.htm for details. */

PDM_bool_t
PDM_polygon_parameterize
(
 const int     numPts,
 const double *pts,
 double       *p0,
 double       *p10,
 double       *l10,
 double       *p20,
 double       *l20,
 double       *n
);


/**
 * \brief Computes polygon parametrization
 *
 * \param [in]  x        Point coordinates to evaluate position
 * \param [in]  numPts  Number of polygon vertices
 * \param [in]  pts     Polygon vertices coordinates
 * \param [out] p0,     Origin vertex
 * \param [out] p10,    First edge vector
 * \param [out] l10,    First edge length
 * \param [out] p20,    First edge vector
 * \param [out] l20,    Second edge vector
 * \param [out] n       Normal
 *
 * \return      \ref Status inside, outside or degenerated
 *
 */

/*  This function is derived from VTK                                      */
/*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen                */
/*  All rights reserved.                                                   */
/*  See Copyright.txt or http://www.kitware.com/Copyright.htm for details. */

PDM_polygon_status_t
PDM_polygon_point_in
(
 const double  x[3],
 const int     numPts,
 const double *pts,
 double       *bounds,
 double       *n
);


/**
 * \brief Computes polygon barycenter
 *
 * \param [in]   numPts  Number of polygon vertices
 * \param [in]   pts     Polygon vertices coordinates
 * \param [out]  bary    Barycenter
 *
 *
 */

void
PDM_polygon_compute_barycenter
(
 const int numPts,
 const double *pts,
 double bary[3]
);
void
PDM_polygon_orthonormal_basis
(
 const int     n_vtx,
 const double *vtx_xyz,
 double        tangent_u[3],
 double        tangent_v[3],
 double        normal[3]
 );

void PDM_polygon_compute_uv_coordinates
(
 const int     n_pts,
 const double *xyz,
 const double *orig_xyz,
 const double *tangent_u,
 const double *tangent_v,
 double       *uv
 );


/**
 * \brief Test if a point is inside a 2d polygon using the Winding Number method
 *        (see http://geomalgorithms.com/a03-_inclusion.html)
 *
 * \param [in]  xy            Point (x,y)-coordinates
 * \param [in]  n_vtx         Number of polygon vertices
 * \param [in]  vtx_xy        Polygon vertices (x,y)-coordinates
 * \param [in]  char_length   Characteristic length (used to scale tolerance)
 * \param [in]  bounds        Bounds (xmin, xmax, ymin, ymax)
 *
 * \return      \ref Status inside, outside or degenerated
 *
 */


PDM_polygon_status_t PDM_polygon_point_in_2d
(
 const double  xy[2],
 const int     n_vtx,
 const double *vtx_xy,
 // const double  char_length,
 double        bounds[4]
 );

/**
 * \brief Test if a point is inside a 3d polygon using the Winding Number method
 *        (see http://geomalgorithms.com/a03-_inclusion.html)
 *
 * \param [in]  xyz           Point (x,y,z)-coordinates
 * \param [in]  n_vtx         Number of polygon vertices
 * \param [in]  vtx_xyz       Polygon vertices (x,y,z)-coordinates
 * \param [in]  char_length   Characteristic length (used to scale tolerance)
 * \param [in]  bounds        Bounds (xmin, xmax, ymin, ymax, zmin, zmax)
 *
 * \return      \ref Status inside, outside or degenerated
 *
 */

PDM_polygon_status_t PDM_polygon_point_in_3d
(
 const double  xyz[3],
 const int     n_vtx,
 const double *vtx_xyz,
 // const double  char_length,
 double        bounds[6],
 double        normal[3]
 );

int PDM_polygon_3d_to_2d
(
 const int    n_vtx,
 const double vtx_xyz[],
 double       vtx_uv[],
 const int    n_pts,
 const double pts_xyz[],
 double       pts_uv[],
 double       normal[3]
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_SURF_MESH_H__ */
