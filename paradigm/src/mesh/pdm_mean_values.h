#ifndef __PDM_MEAN_VALUES_H__
#define __PDM_MEAN_VALUES_H__

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

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Compute mean values of a list of points in a polygon
 *
 * \param [in]    n_pts            Number of points to locate
 * \param [in]    pts              xyz-Coordinates of points locate
 * \param [in]    n_vtx            Number of polygon vertices
 * \param [in]    poly_vtx         Polygon connectivity
 * \param [in]    mesh_vtx_coords  Coordinates of mesh vertices
 * \param [inout] mesh_vtx_coords  Mean value coordinates of points to locate
 *
 *
 */

void
PDM_mean_value_coordinates_polygon_2d
(
 const int    n_vtx,
 const double vtx_coord[],
 const int    n_pts,
 const double pts_coord[],
 double       mean_value_coord[]
 );

void
PDM_mean_value_coordinates_polygon_3d
(
 const int    n_vtx,
 const double vtx_coord[],
 const int    n_pts,
 const double pts_coord[],
 double       mean_value_coord[]
 );



void
PDM_mean_value_coordinates_polyhedron
(
 const int         n_vtx,
 const double      vtx_coord[],
 const PDM_l_num_t n_face,
 const PDM_l_num_t face_vtx_idx[],
 const PDM_l_num_t face_vtx[],
 const int         face_orientation[],
 const double      pt_coord[],
 double            mean_value_coord[]
 );


void
PDM_mean_values_polyhedron
(
 const int         n_vtx,
 const double      vtx_coord[],
 const PDM_l_num_t n_face,
 const PDM_l_num_t face_vtx_idx[],
 const PDM_l_num_t face_vtx[],
 const int         face_orientation[],
 const double      pt_coord[],
 double            weights[]
 );

void
PDM_mean_values_polygon
(
 const int         n_vtx,
 const PDM_l_num_t face_vtx[],
 const double      vtx_coord[],
 const double      pt_coord[],
 double            weights[]
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MEAN_VALUES_H__ */
