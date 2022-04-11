#ifndef __PDM_POINT_LOCATION_H__
#define __PDM_POINT_LOCATION_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mesh_nodal.h"

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
 * \brief Locate a set points inside a set of elements
 *
 * Elements are ordered by type (points, lines, triangles, quadrangles,
 * polygons, tetrahedra, pyramids, prisms, hexahedra, polyhedra).
 *
 * \param [in]   type_idx           Index for the element types (size = 11)
 * \param [in]   elt_vtx_idx        Index of the element-vertex connectivity
 * \param [in]   elt_vtx_coord      Coordinates of the elements' vertices
 * \param [in]   poly3d_face_idx    Index of the element-face connectivity (only for polyhedra)
 * \param [in]   face_vtx_idx       Index for the face-vertex connectivity
 * \param [in]   face_vtx           Face-vertex connectivity
 * \param [in]   face_orientation   Orientation of the faces
 * \param [in]   pts_idx            Index of points (size = n_elt + 1)
 * \param [in]   pts_coord          Coordinates of the points to locate
 * \param [in]   tolerance          Geometric tolerance
 * \param [out]  distance           Distance from points to elements (< 0 if inside, > 0 if outside)
 * \param [out]  projected_coord    Coordinates of the projection of the points on the elements
 * \param [out]  bar_coord_idx      Index for the mean-value coordinates of the projections
 * \param [out]  bar_coord          Mean-value coordinates of the projections
 *
 */

void
PDM_point_location_nodal
(
 const int           type_idx[],
 const int           elt_vtx_idx[],
 const double        elt_vtx_coord[],
 const PDM_l_num_t   poly3d_face_idx[],
 const PDM_l_num_t   face_vtx_idx[],
 const PDM_l_num_t   face_vtx[],
 const int           face_orientation[],
 const int           pts_idx[],
 const double        pts_coord[],
 const double        tolerance,
 double            **distance,
 double            **projected_coord,
 int               **bar_coord_idx,
 double            **bar_coord
 );


/**
 * \brief Compute hexahedron, pyramid, or prism parametric coordinates for a
 * given point.
 *
 * This function is adapted from the CGNS interpolation tool.
 *
 * \param [in]   elt_type       Type of element
 * \param [in]   point_coords   Point coordinates
 * \param [in]   vertex_coords  Pointer to element vertex coordinates
 * \param [in]   tolerance      Location tolerance factor
 * \param [out]  uvw            Parametric coordinates of point in element
 *
 */

PDM_bool_t
PDM_point_location_compute_uvw
(
 const PDM_Mesh_nodal_elt_t elt_type,
 const double               point_coords[3],
 const double               vertex_coords[],
 const double               tolerance,
 double                     uvw[3]
 );


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_POINT_LOCATION_H__ */
