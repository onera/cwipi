#ifndef __PDM_TRIANGULATE_H__
#define __PDM_TRIANGULATE_H__

/*============================================================================
 * Triangulation of a polygon
 *============================================================================*/
/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2005-2007  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*
 * Pointer to structure maintaining the state of the current triangulation;
 * the structure itself is private.
 */

typedef struct _PDM_triangulate_state_t PDM_triangulate_state_t;

/*
 * Describe how the resulting triangle connectivity is defined.
 */

typedef enum {

  PDM_TRIANGULATE_MESH_DEF,      /* Definition by mesh vertex numbers */
  PDM_TRIANGULATE_ELT_DEF        /* Definition by local (element) vertex
                                    position (1 to n) */

} PDM_triangulate_def_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a structure necessary to the polygon triangulation algorithm.
 *
 * parameters:
 *   n_vertices_max    <-- maximum expected number of vertices per polygon.
 *
 * returns:
 *   pointer to polygon triangulation state structure.
 *----------------------------------------------------------------------------*/

PDM_triangulate_state_t *
PDM_triangulate_state_create(const int  n_vertices_max);

/*----------------------------------------------------------------------------
 * Destroy a structure necessary to the polygon triangulation algorithm.
 *
 * parameters:
 *   this_state  <-> pointer to structure that should be destroyed.
 *
 * returns:
 *   NULL pointer.
 *----------------------------------------------------------------------------*/

PDM_triangulate_state_t *
PDM_triangulate_state_destroy(PDM_triangulate_state_t  *this_state);

/*----------------------------------------------------------------------------
 * Triangulate a polygonal face.
 *
 * For a polygon with n vertices, we should obtain a triangluation with
 * (n-2) triangles and (2n-3) edges. If the polygon_vertices argument
 * is NULL, 1, 2, ...,n local numbering is implied.
 *
 * parameters:
 *   dim               <-- spatial dimension (2 or 3).
 *   n_vertices        <-- number of vertices defining the polygon.
 *   coords            <-- coordinates of the triangulation's vertices.
 *   parent_vertex_num <-- optional indirection to vertex coordinates (1 to n).
 *   polygon_vertices  <-- polygon connectivity; size: n_vertices or empty.
 *   mode              <-- triangles connectivity by vertex number or
 *                         polygon vertex index (1 to n).
 *   triangle_vertices --> triangles connectivity;
 *                         size: (n_vertices - 2) * 3.
 *   state             <-> associated triangulation state structure.
 *
 * returns:
 *   number of resulting triangles.
 *----------------------------------------------------------------------------*/

int
PDM_triangulate_polygon(int                             dim,
                        int                             n_vertices,
                        const double                    coords[],
                        const PDM_l_num_t               parent_vertex_num[],
                        const PDM_l_num_t               polygon_vertices[],
                        PDM_triangulate_def_t           mode,
                        PDM_l_num_t                     triangle_vertices[],
                        PDM_triangulate_state_t  *const state);

/*----------------------------------------------------------------------------
 * Triangulate a quadrangle.
 *
 * A convex quadrangle is divided into two triangles along its shortest
 * diagonal. A non-convex quadrangle may only be divided along the diagonal
 * which lies inside the quadrangle.
 *
 * If the quadrangle_vertices argument is NULL, 1, 2, ...,n local numbering
 * is implied.
 *
 * parameters:
 *   dim                  <-- spatial dimension (2 or 3).
 *   coords               <-- coordinates of the triangulation's vertices.
 *   parent_vertex_num    <-- optional indirection to vertex coordinates
 *   quadrangle_vertices  <-- polygon connectivity; size: n_vertices or empty.
 *   triangle_vertices    --> triangles connectivity; size: 2 * 3.
 *
 * returns:
 *   number of resulting triangles.
 *----------------------------------------------------------------------------*/

int
PDM_triangulate_quadrangle(int                dim,
                           const double       coords[],
                           const PDM_l_num_t  parent_vertex_num[],
                           const PDM_l_num_t  quadrangle_vertices[],
                           PDM_l_num_t        triangle_vertices[]);


int
PDM_triangulate_pyramid (int               dim,
                         const double      coords[],
                         const PDM_l_num_t parent_vertex_num[],
                         const PDM_l_num_t pyramid_vertices[],
                         PDM_l_num_t       tetra_vertices[]);

int
PDM_triangulate_prism (int               dim,
                       const double      coords[],
                       const PDM_l_num_t parent_vertex_num[],
                       const PDM_l_num_t prism_vertices[],
                       PDM_l_num_t       tetrahedron_vertices[]);

int
PDM_triangulate_hexahedron (int               dim,
                            const double      coords[],
                            const PDM_l_num_t parent_vertex_num[],
                            const PDM_l_num_t hexa_vertices[],
                            PDM_l_num_t       tetrahedron_vertices[]);
/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_TRIANGULATE_H__ */
