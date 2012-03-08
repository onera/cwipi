/*============================================================================
 * Locate local points in a nodal representation associated with a mesh
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2007-2009  EDF

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

/*----------------------------------------------------------------------------*/

#include "fvmc_config.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bftc_error.h>
#include <bftc_mem.h>
#include <bftc_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_defs.h"
#include "fvmc_nodal.h"
#include "fvmc_nodal_priv.h"
#include "fvmc_triangulate.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_point_location.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/* In case <math.h> does not define HUGE_VAL, use a "safe" value */
#if !defined(HUGE_VAL)
#define HUGE_VAL 1.0e+30
#endif

#define MAXSTATICSUBDIVISION 1000

/* Geometric operation macros*/

enum {X, Y, Z};

#define _DOT_PRODUCT(vect1, vect2) \
  (vect1[X] * vect2[X] + vect1[Y] * vect2[Y] + vect1[Z] * vect2[Z])

#define _MODULE(vect) \
  sqrt(vect[X] * vect[X] + vect[Y] * vect[Y] + vect[Z] * vect[Z])

#define _CROSS_PRODUCT(prod_vect, vect1, vect2)  \
  (prod_vect[X] = vect1[Y] * vect2[Z] - vect2[Y] * vect1[Z], \
   prod_vect[Y] = vect2[X] * vect1[Z] - vect1[X] * vect2[Z], \
   prod_vect[Z] = vect1[X] * vect2[Y] - vect2[X] * vect1[Y])

#define _DOT_PRODUCT_2D(vect1, vect2) \
  (vect1[X] * vect2[X] + vect1[Y] * vect2[Y])

/*============================================================================
 * Type definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure defining a local octree (3d)
 *----------------------------------------------------------------------------*/

typedef struct {

  fvmc_lnum_t  octant_id[8];   /* Ids of sub-octants in octree array */
  fvmc_lnum_t  idx[9];         /* Start index of point list for each octant */
  fvmc_lnum_t  n_points;       /* Number of points in octree */

} _octant_t;

typedef struct {

  size_t       n_points;      /* Number of points in octree */
  size_t       n_nodes;       /* Current number of nodes in octree */
  size_t       n_nodes_max;   /* Maximum number of nodes in octree */
  size_t       n_static_subdivision;/* Number of static subdivisions 
                                       (all points in the same child node) */

  double       extents[6];    /* Associated extents */

  fvmc_lnum_t  *point_ids;     /* Id's of points sorted by octree
                                 (size: n_points + 1) */
  _octant_t   *nodes;         /* Array of octree nodes
                                 (size: n_nodes_max) */

} _octree_t;

/*----------------------------------------------------------------------------
 * Structure defining a local quadtree (2d)
 *----------------------------------------------------------------------------*/

typedef struct {

  fvmc_lnum_t  quadrant_id[4]; /* Id of sub-quadrants in quadtree array */
  fvmc_lnum_t  idx[5];         /* Start index of point list for each quadrant */
  fvmc_lnum_t  n_points;       /* Number of points in quadtree */

} _quadrant_t;

typedef struct {

  size_t        n_points;     /* Number of points in quadtree */
  size_t        n_nodes;      /* Current number of nodes in quadtree */
  size_t        n_nodes_max;  /* Maximum number of nodes in quadtree */
  size_t       n_static_subdivision; /* Number of static subdivisions 
                                       (all points in the same child node) */

  double        extents[4];   /* Associated extents */

  fvmc_lnum_t   *point_ids;    /* Id's of points sorted by quadtree
                                 (size: n_points + 1) */
  _quadrant_t  *nodes;        /* Array of quadtree nodes
                                 (size: n_nodes_max) */

} _quadtree_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

static double      _epsilon_denom = 1.e-30;       /* Minimum denominator */
static double      _epsilon_multi_point = 1.e-12; /* Minimum distance between 
                                                     2 points */

static fvmc_lnum_t  _octree_threshold = 4; /* Number of points in octree node
                                             under which the node is final */

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Test if two extents intersect
 *
 * parameters:
 *   dim             <-- spatial (coordinates) dimension
 *   extents_1       <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *   extents_2       <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *
 * returns:
 *   true if extents intersect, false otherwise
 *----------------------------------------------------------------------------*/

inline static _Bool
_intersect_extents(int           dim,
                   const double  extents_1[],
                   const double  extents_2[])
{
  int i;
  _Bool retval = true;

  for (i = 0; i < dim; i++) {
    if (   (extents_1[i] > extents_2[i + dim])
        || (extents_2[i] > extents_1[i + dim])) {
      retval = false;
      break;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Test if a point is within given extents
 *
 * parameters:
 *   dim             <-- spatial (coordinates) dimension
 *   coords          <-- coordinates: x, y, ...
 *                       size: dim
 *   extents         <-- extents: x_min, y_min, ..., x_max, y_max, ...
 *                       size: dim*2
 *
 * returns:
 *   true if point lies within extents, false otherwise
 *----------------------------------------------------------------------------*/

inline static _Bool
_within_extents(int                dim,
                const fvmc_coord_t  coords[],
                const double       extents[])
{
  int i;
  _Bool retval = true;

  for (i = 0; i < dim; i++) {
    if (   (coords[i] < extents[i])
        || (coords[i] > extents[i + dim])) {
      retval = false;
      break;
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Update element extents with a given vertex
 *
 * parameters:
 *   dim               <-- spatial (coordinates) dimension
 *   vertex_id         <-- vertex index (0 to n-1)
 *   parent_vertex_num <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   elt_extents       <-> extents associated with element:
 *                         x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *   elt_initialized   <-> are extents already initialized for this vertex
 *                         (for all element vertices except the first) ?
 *----------------------------------------------------------------------------*/

inline static void
_update_elt_extents(int                 dim,
                    fvmc_lnum_t          vertex_id,
                    const fvmc_lnum_t   *parent_vertex_num,
                    const fvmc_coord_t   vertex_coords[],
                    double              elt_extents[],
                    _Bool              *elt_initialized)
{
  fvmc_lnum_t  i, coord_idx;

  if (parent_vertex_num == NULL)
    coord_idx = vertex_id;
  else
    coord_idx = parent_vertex_num[vertex_id] - 1;

  if (*elt_initialized == false) {
    for (i = 0; i < dim; i++) {
      elt_extents[i]       = vertex_coords[(coord_idx * dim) + i];
      elt_extents[i + dim] = vertex_coords[(coord_idx * dim) + i];
    }
    *elt_initialized = true;
  }
  else {
    for (i = 0; i < dim; i++) {
      if (elt_extents[i]       > vertex_coords[(coord_idx * dim) + i])
        elt_extents[i]       = vertex_coords[(coord_idx * dim) + i];
      if (elt_extents[i + dim] < vertex_coords[(coord_idx * dim) + i])
        elt_extents[i + dim] = vertex_coords[(coord_idx * dim) + i];
    }

  }

}

/*----------------------------------------------------------------------------
 * Adjust element extents with search tolerance and update global extents
 *
 * parameters:
 *   dim         <-- spatial (coordinates) dimension
 *   elt_dim     <-- element dimension
 *   tolerance   <-- addition to local extents of each element:
 *                   extent = base_extent * (1 + tolerance)
 *   elt_extents <-> extents associated with element:
 *                   x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *----------------------------------------------------------------------------*/

inline static void
_elt_extents_finalize(int               dim,
                      int               elt_dim,
                      double            tolerance,
                      double  *restrict elt_extents)
{
  int i;
  double delta[3];

  for (i = 0; i < dim; i++)
    delta[i] = (elt_extents[i+dim] - elt_extents[i]) * tolerance;

  if (elt_dim < dim) {
    double delta_max = delta[0];  /* for 1d or 2d elements, ensure */
    for (i = 0; i < dim; i++) {   /* search extent "thickness" */
      if (delta[i] > delta_max)
        delta_max = delta[i];
    }
    for (i = 0; i < dim; i++)
      delta[i] = delta_max;
  }

  for (i = 0; i < dim; i++) {
    elt_extents[i]     = elt_extents[i]     - delta[i];
    elt_extents[i+dim] = elt_extents[i+dim] + delta[i];
  }

}

/*----------------------------------------------------------------------------
 * Compute extents of a point set
 *
 * parameters:
 *   dim          <-- space dimension of points to locate
 *   n_points     <-- number of points to locate
 *   point_index  <-- optional indirection array to point_coords
 *                    (1 to n_points numbering)
 *   point_coords <-- coordinates of points to locate
 *                    (dimension: dim * n_points)
 *   extents      --> extents associated with mesh:
 *                    x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *----------------------------------------------------------------------------*/

static void
_point_extents(const int            dim,
               const fvmc_lnum_t     n_points,
               const fvmc_lnum_t     point_index[],
               const fvmc_coord_t    point_coords[],
               double               extents[])
{
  int i;
  fvmc_lnum_t j, coord_idx;

  /* initialize extents in case mesh is empty or dim < 3 */
  for (i = 0; i < dim; i++) {
    extents[i]       =  HUGE_VAL;
    extents[i + dim] = -HUGE_VAL;
  }

  /* Compute extents */

  if (point_index != NULL) {

    for (j = 0; j < n_points; j++) {
      coord_idx = point_index[j] - 1;
      for (i = 0; i < dim; i++) {
        if (extents[i]       > point_coords[(coord_idx * dim) + i])
          extents[i]       = point_coords[(coord_idx * dim) + i];
        if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
          extents[i + dim] = point_coords[(coord_idx * dim) + i];
      }
      
    }
  }

  else {

    for (coord_idx = 0; coord_idx < n_points; coord_idx++) {
      for (i = 0; i < dim; i++) {
        if (extents[i]       > point_coords[(coord_idx * dim) + i])
          extents[i]       = point_coords[(coord_idx * dim) + i];
        if (extents[i + dim] < point_coords[(coord_idx * dim) + i])
          extents[i + dim] = point_coords[(coord_idx * dim) + i];
      }
    }
  }
}

/*----------------------------------------------------------------------------
 * Updates the location[] and distance[] arrays associated with a set
 * of 1d points for points that are in a given element extent, based only
 * on this extent only (temporary, unoptimzed location).
 *
 * parameters:
 *   elt_num         <-- number of element corresponding to extents
 *   extents         <-> extents associated with element:
 *                       x_min, x_max (size: 2)
 *   n_points        <-- number of points to locate
 *   point_coords    <-- point coordinates
 *   location        <-> number of element containing or closest to each
 *                       point (size: n_points)
 *   distance        <-> distance from point to element indicated by
 *                       location[]: < 0 if unlocated, 0 - 1 if inside,
 *                       > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_by_extents_1d(fvmc_lnum_t         elt_num,
                      const double       extents[],
                      fvmc_lnum_t         n_points,
                      const fvmc_coord_t  point_coords[],
                      fvmc_lnum_t         location[],
                      float              distance[])
{
  fvmc_lnum_t  i;

  /* For now, we base a minimal location test on the element extents */
  /* The behavior is quadradic, nothing is optimized yet */

  for (i = 0; i < n_points; i++) {

    double elt_coord_max = -1;
    double elt_coord = -1;

    double cur_coord = point_coords[i];

    elt_coord =   (cur_coord - 0.5*(extents[1] + extents[0]))
                / (            0.5*(extents[1] - extents[0]));

    elt_coord = FVMC_ABS(elt_coord);

    if (elt_coord > elt_coord_max)
      elt_coord_max = elt_coord;

    if (  (distance[i] < 0 && elt_coord_max < 1)
        || elt_coord_max < distance[i]) {

      location[i] = elt_num;
      distance[i] = elt_coord_max;

    }

  }

}

/*----------------------------------------------------------------------------
 * Updates the location[] and distance[] arrays associated with a set
 * of points for points that are in a given element extent, based on
 * a query of points with this extent.
 *
 * parameters:
 *   elt_num            <-- number of element corresponding to extents
 *   dim                <-- spatial dimension
 *   extents            <-> extents associated with element:
 *                          x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *   point_coords       <-- point coordinates (size > n_points_in_extent*dim)
 *   n_points_in_extent <-- number of points in extents
 *   points_in_extent   <-- ids of points in extents
 *   location           <-> number of element containing or closest to each
 *                          point (size: n_points)
 *   distance           <-> distance from point to element indicated by
 *                          location[]: < 0 if unlocated, 0 - 1 if inside,
 *                          > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_in_extents(const fvmc_lnum_t    elt_num,
                   const int           dim,
                   const double        extents[],
                   const fvmc_coord_t   point_coords[],
                   fvmc_lnum_t          n_points_in_extents,
                   const fvmc_lnum_t    points_in_extents[],
                   fvmc_lnum_t          location[],
                   float               distance[])
{
  fvmc_lnum_t  i, j, k;

  /* For now, we base a minimal location test on the element extents */
  /* The behavior is quadradic, nothing is optimized yet */

  for (i = 0; i < n_points_in_extents; i++) {

    double elt_coord_max = -1;
    double elt_coord = -1;

    j = points_in_extents[i];

    for (k = 0; k < dim; k++) {

      double cur_coord = point_coords[j*dim + k];

      elt_coord =   (cur_coord - 0.5*(extents[k+dim] + extents[k]))
                  / (            0.5*(extents[k+dim] - extents[k]));

      elt_coord = FVMC_ABS(elt_coord);

      if (elt_coord > elt_coord_max)
        elt_coord_max = elt_coord;

    }

    if (  (distance[j] < 0 && elt_coord_max < 1)
        || elt_coord_max < distance[j]) {

      location[j] = elt_num;
      distance[j] = elt_coord_max;

    }

  }

}

/*----------------------------------------------------------------------------
 * Compare points coordinates to check if it is a multi point
 *
 * parameters:
 *   dim                <-- Dimesion                   
 *   n_point            <-- Number of point to compare
 *   point_coords       <-- point coordinates
 *   point_idx          <-- point indexes
 *   count              <-- Number of point to compare
 *
 * return               
 *----------------------------------------------------------------------------*/

static int
_check_multi_point(const int           dim,
                   const fvmc_lnum_t    n_point,
                   const fvmc_coord_t   point_coords[],
                   const fvmc_lnum_t    point_idx[])
{
  fvmc_lnum_t i, j, k;
  double dist;

  for(i = 0; i < n_point; i++) {
    for(j = 0; j < n_point; j++) {
      dist = 0;
      for (k = 0; k < dim; k++) {
        dist += (point_coords[point_idx[j] *dim + k] - point_coords[point_idx[i] *dim + k]) *
                (point_coords[point_idx[j] *dim + k] - point_coords[point_idx[i] *dim + k]);
        if (dist > _epsilon_multi_point) 
          return 0;
      }
    }
  }
  return 1;
}

/*----------------------------------------------------------------------------
 * Build a local octree's leaves.
 *
 * parameters:
 *   extents            <-> extents associated with node:
 *                          x_min, y_min, z_min, x_max, y_max, z_max (size: 6)
 *   point_coords       <-- point coordinates
 *   point_ids_tmp      <-- temporary point indexes
 *   pos_tmp            <-- temporary point position in octree
 *   octree             <-> current octree structure
 *   point_range        <-> start and past-the end index in point_idx
 *                          for current node (size: 2)
 *----------------------------------------------------------------------------*/

static void
_build_octree_leaves(const double        extents[],
                     const fvmc_coord_t   point_coords[],
                     fvmc_lnum_t         *point_ids_tmp,
                     _octree_t          *octree,
                     fvmc_lnum_t          point_range[2])
{
  fvmc_lnum_t i, j, k, _n_nodes, _n_points, tmp_size;

  fvmc_lnum_t count[8], idx[9], octant_id[8];
  double mid[3], sub_extents[6];
  _octant_t  *_node;

  int octant_mask[3] = {4, 2, 1}; /* pow(2, 2), pow(2, 1), pow(2,0) */
  int check_multi_point;

  _n_nodes = octree->n_nodes;
  tmp_size = octree->n_nodes;


  /* Resize octree if necesary */

  if (octree->n_nodes >= octree->n_nodes_max) {
    if (octree->n_nodes == 0) {
      octree->n_nodes = 1;
      octree->n_nodes_max = 8;
    }
    octree->n_nodes_max *= 2;
    BFTC_REALLOC(octree->nodes, octree->n_nodes_max, _octant_t);
  }

  /* Number of points */

  _n_points = point_range[1] - point_range[0];

  /* Extents center */

  for (j = 0; j < 3; j++)
    mid[j]= (extents[j] + extents[j + 3]) * 0.5;

  for (j = 0; j < 8; j++) {
    count[j] = 0;
    octant_id[j] = -1;
  }

  /* Count points in each octant */

  for (i = point_range[0]; i < point_range[1]; i++) {

    for (j = 0, k = 0; j < 3; j++) {
      if (point_coords[octree->point_ids[i]*3 + j] > mid[j])
        k += octant_mask[j];
    }

    count[k] += 1;
  }

  /* Build index */

  idx[0] = 0;
  for (j = 0; j < 8; j++)
    idx[j+1] = idx[j] + count[j];

  for (j = 0; j < 8; j++)
    count[j] = 0;

  for (i = point_range[0], j = 0; i < point_range[1]; i++) {

    for (j = 0, k = 0; j < 3; j++) {
      if (point_coords[octree->point_ids[i]*3 + j] > mid[j])
        k += octant_mask[j];
    }

    point_ids_tmp[idx[k] + count[k]] = octree->point_ids[i];
    count[k] += 1;
  }
  
  /* Check if this subdivision is static
     and check coordinates to find multi point */

  for (j = 0; j < 8; j++) {
    if (count[j] == _n_points) {
      octree->n_static_subdivision += 1;    
      break;
    }
  }
   
  if (j == 8)
    octree->n_static_subdivision = 0;
  else {
    if (octree->n_static_subdivision >= MAXSTATICSUBDIVISION) { 
      check_multi_point = _check_multi_point(3,
                                             count[j],
                                             point_coords,
                                             point_ids_tmp + idx[j]);
      if (check_multi_point) {
        count[j] = 0; /* Stop the recursion */
      }
    }
  }

  for (i = point_range[0], j = 0; i < point_range[1]; i++, j++)
    octree->point_ids[i] = point_ids_tmp[j];

  for (i = 0; i < 9; i++)
    idx[i] = point_range[0] + idx[i];

  /* Build leaves recursively */

  for (i = 0; i < 8; i++) {

    if (count[i] > _octree_threshold) {

      tmp_size++;

      octant_id[i] = tmp_size;

      if (i < 4) {
        sub_extents[0] = extents[0];
        sub_extents[3] = mid[0];
      }
      else {
        sub_extents[0] = mid[0];
        sub_extents[3] = extents[3];
      }
      /* 1.0e-12 term in assert() used to allow for
         truncation error in for xmin = xmax case */
      assert(sub_extents[0] < sub_extents[3] + 1.0e-12);

      if (i%4 < 2) {
        sub_extents[1] = extents[1];
        sub_extents[4] = mid[1];
      }
      else {
        sub_extents[1] = mid[1];
        sub_extents[4] = extents[4];
      }
      assert(sub_extents[1] < sub_extents[4] + 1.0e-12);

      if (i%2 < 1) {
        sub_extents[2] = extents[2];
        sub_extents[5] = mid[2];
      }
      else {
        sub_extents[2] = mid[2];
        sub_extents[5] = extents[5];
      }
      assert(sub_extents[2] < sub_extents[5] + 1.0e-12);

      octree->n_nodes = tmp_size;
 
      _build_octree_leaves(sub_extents,
                           point_coords,
                           point_ids_tmp,
                           octree,
                           idx + i);

      tmp_size = octree->n_nodes;
    }
  }

  /* Finalize node */

  _node = octree->nodes + _n_nodes;

  for (i = 0; i < 9; i++)
    _node->idx[i] = idx[i];

  for (i = 0; i < 8; i++)
    _node->octant_id[i] = octant_id[i];

  _node->n_points = _n_points;
}

/*----------------------------------------------------------------------------
 * Build an octree structure to locate 3d points in mesh.
 *
 * parameters:
 *   n_points        <-- number of points to locate
 *   point_coords    <-- point coordinates
 *
 * returns:
 *   pointer to local octree structure
 *----------------------------------------------------------------------------*/

static _octree_t
_build_octree(fvmc_lnum_t         n_points,
              const fvmc_coord_t  point_coords[])
{
  size_t i;
  fvmc_lnum_t point_range[2];
  _octree_t _octree;

  int *point_ids_tmp = NULL;

  /* Initialization */

  point_range[0] = 0;
  point_range[1] = n_points;

  _octree.n_points = n_points;
  _octree.n_nodes = 0;
  _octree.n_nodes_max = 0;
  _octree.n_static_subdivision = 0;
  _octree.nodes = NULL;
  _octree.point_ids = NULL;

  if (n_points > 0) {

    _point_extents(3,
                   n_points,
                   NULL,
                   point_coords,
                   _octree.extents);

    BFTC_MALLOC(_octree.point_ids, _octree.n_points, fvmc_lnum_t);

    for (i = 0; i < _octree.n_points; i++)
      _octree.point_ids[i] = i;

    BFTC_MALLOC(point_ids_tmp, n_points, int);

    _build_octree_leaves(_octree.extents,
                         point_coords,
                         point_ids_tmp,
                         &_octree,
                         point_range);

    BFTC_FREE(point_ids_tmp);

  }

  return _octree;
}

/*----------------------------------------------------------------------------
 * Free an octree structure.
 *
 * parameters:
 *   octree <-> octree structure whose elements are to be freed
 *
 * returns:
 *   pointer to local octree structure
 *----------------------------------------------------------------------------*/

static void
_free_octree(_octree_t *octree)
{

  octree->n_points = 0;
  octree->n_nodes = 0;
  octree->n_nodes_max = 0;

  BFTC_FREE(octree->nodes);
  BFTC_FREE(octree->point_ids);
}

/*----------------------------------------------------------------------------
 * locate points in box defined by extents in an octant.
 *
 * parameters:
 *   search_extents  <-- extents associated with element:
 *   point_coords    <-- point coordinates
 *                       x_min, y_min, z_min, x_max, y_max, z_max (size: 6)
 *   octree          <-- point octree
 *   node_extents    <-- extents of octant
 *   node_id         <-- id of node in octree
 *   loc_point_ids   --> ids of located points (max size: octree->n_points)
 *   n_loc_points    --> number of points located
 *----------------------------------------------------------------------------*/

static void
_query_octree_node(const double        extents[],
                   const fvmc_coord_t   point_coords[],
                   const _octree_t    *octree,
                   const double        node_extents[],
                   int                 node_id,
                   fvmc_lnum_t         *loc_point_ids,
                   fvmc_lnum_t         *n_loc_points)
{
  _octant_t *node;
  int i, j, k;
  int dim = 3;
  double sub_extents[6], mid[3];

  node = octree->nodes + node_id;

  if (_intersect_extents(dim, node_extents, extents)) {

    for (j = 0; j < dim; j++)
      mid[j]= (node_extents[j] + node_extents[j + dim]) * 0.5;

    /* Loop on node leaves */

    for (i = 0; i < 8; i++) {

      /* Compute octant extents */

      if (i < 4) {
        sub_extents[0] = node_extents[0];
        sub_extents[3] = mid[0];
      }
      else {
        sub_extents[0] = mid[0];
        sub_extents[3] = node_extents[3];
      }

      if (i%4 < 2) {
        sub_extents[1] = node_extents[1];
        sub_extents[4] = mid[1];
      }
      else {
        sub_extents[1] = mid[1];
        sub_extents[4] = node_extents[4];
      }
      if (i%2 < 1) {
        sub_extents[2] = node_extents[2];
        sub_extents[5] = mid[2];
      }
      else {
        sub_extents[2] = mid[2];
        sub_extents[5] = node_extents[5];
      }

      /* Search recursively if octant is not final */

      if (node->octant_id[i] > -1)
        _query_octree_node(extents,
                           point_coords,
                           octree,
                           sub_extents,
                           node->octant_id[i],
                           loc_point_ids,
                           n_loc_points);

      else {

        /* Update list of points located */

        if (_intersect_extents(dim, sub_extents, extents)) {

          for (k = node->idx[i]; k < node->idx[i+1]; k++) {

            fvmc_lnum_t point_id = octree->point_ids[k];

            if (_within_extents(dim, point_coords + point_id*dim, extents)) {
              loc_point_ids[*n_loc_points] = point_id;
              (*n_loc_points)++;
            }

          }
        }

      }

    } /* End of loop on node leaves */

  }

}

/*----------------------------------------------------------------------------
 * locate points in box defined by extents using an octree.
 *
 * parameters:
 *   search_extents  <-- extents associated with element:
 *                       x_min, y_min, ..., x_max, y_max, ... (size: 2*dim)
 *   point_coords    <-- point coordinates
 *   octree          <-- point octree
 *   n_loc_points    --> number of points located
 *   loc_point_ids   --> ids of located points (max size: octree->n_points)
 *----------------------------------------------------------------------------*/

static void
_query_octree(const double        extents[],
              const fvmc_coord_t   point_coords[],
              const _octree_t    *octree,
              fvmc_lnum_t         *n_loc_points,
              fvmc_lnum_t          loc_point_ids[])
{
  *n_loc_points = 0;

  if (octree->n_points > 0)
    _query_octree_node(extents,
                       point_coords,
                       octree,
                       octree->extents,
                       0,
                       loc_point_ids,
                       n_loc_points);
}

/*----------------------------------------------------------------------------
 * Build a local quadtree's leaves.
 *
 * parameters:
 *   extents            <-> extents associated with node:
 *                          x_min, y_min, x_max, y_max, (size: 4)
 *   point_coords       <-- point coordinates
 *   point_ids_tmp      <-- temporary point indexes
 *   pos_tmp            <-- temporary point position in quadtree
 *   octree             <-> current quadtree structure
 *   point_range        <-> start and past-the end index in point_idx
 *                          for current node (size: 2)
 *----------------------------------------------------------------------------*/

static void
_build_quadtree_leaves(const double        extents[],
                       const fvmc_coord_t   point_coords[],
                       fvmc_lnum_t         *point_ids_tmp,
                       _quadtree_t        *quadtree,
                       fvmc_lnum_t          point_range[2])
{
  fvmc_lnum_t i, j, k, _n_nodes, _n_points, tmp_size;

  fvmc_lnum_t count[4], idx[5], quadrant_id[4];
  double mid[2], sub_extents[4];
  _quadrant_t  *_node;

  int quadrant_mask[2] = {2, 1}; /* pow(2, 1), pow(2,0) */
  int check_multi_point;

  _n_nodes = quadtree->n_nodes;
  tmp_size = quadtree->n_nodes;

  /* Resize quadtree if necesary */

  if (quadtree->n_nodes >= quadtree->n_nodes_max) {
    if (quadtree->n_nodes == 0) {
      quadtree->n_nodes = 1;
      quadtree->n_nodes_max = 4;
    }
    quadtree->n_nodes_max *= 2;
    BFTC_REALLOC(quadtree->nodes, quadtree->n_nodes_max, _quadrant_t);
  }

  /* Number of points */

  _n_points = point_range[1] - point_range[0];

  /* Extents center */

  for (j = 0; j < 2; j++)
    mid[j]= (extents[j] + extents[j + 2]) * 0.5;

  for (j = 0; j < 4; j++) {
    count [j] = 0;
    quadrant_id[j] = -1;
  }

  /* Count points in each quadrant */

  for (i = point_range[0]; i < point_range[1]; i++) {

    for (j = 0, k = 0; j < 2; j++) {
      if (point_coords[quadtree->point_ids[i]*2 + j] > mid[j])
        k += quadrant_mask[j];
    }

    count[k] += 1;
  }

  /* Build index */

  idx[0] = 0;
  for (j = 0; j < 4; j++)
    idx[j+1] = idx[j] + count[j];

  for (j = 0; j < 4; j++)
    count[j] = 0;

  for (i = point_range[0], j = 0; i < point_range[1]; i++) {

    for (j = 0, k = 0; j < 2; j++) {
      if (point_coords[quadtree->point_ids[i]*2 + j] > mid[j])
        k += quadrant_mask[j];
    }

    point_ids_tmp[idx[k] + count[k]] = quadtree->point_ids[i];
    count[k] += 1;
  }

  /* Check if this subdivision is static
     and check coordinates to find multi point */

  for (j = 0; j < 4; j++) {
    if (count[j] == _n_points) {
      quadtree->n_static_subdivision += 1;    
      break;
    }
  }
   
  if (j == 4)
    quadtree->n_static_subdivision = 0;
  else {
    if (quadtree->n_static_subdivision >= MAXSTATICSUBDIVISION) { 
      check_multi_point = _check_multi_point(2,
                                             count[j],
                                             point_coords,
                                             point_ids_tmp + idx[j]);
      if (check_multi_point) {
        count[j] = 0; /* Stop the recursion */
      }
    }
  }

  for (i = point_range[0], j = 0; i < point_range[1]; i++, j++)
    quadtree->point_ids[i] = point_ids_tmp[j];

  for (i = 0; i < 5; i++)
    idx[i] = point_range[0] + idx[i];

  /* Build leaves recursively */

  for (i = 0; i < 4; i++) {

    if (count[i] > _octree_threshold) {

      tmp_size++;

      quadrant_id[i] = tmp_size;

      if (i < 2) {
        sub_extents[0] = extents[0];
        sub_extents[2] = mid[0];
      }
      else {
        sub_extents[0] = mid[0];
        sub_extents[2] = extents[2];
      }
      assert(sub_extents[0] < sub_extents[2] + 1.0e-12);

      if (i%2 < 1) {
        sub_extents[1] = extents[1];
        sub_extents[3] = mid[1];
      }
      else {
        sub_extents[1] = mid[1];
        sub_extents[3] = extents[3];
      }
      assert(sub_extents[1] < sub_extents[3] + 1.0e-12);

      quadtree->n_nodes = tmp_size;

      _build_quadtree_leaves(sub_extents,
                             point_coords,
                             point_ids_tmp,
                             quadtree,
                             idx + i);

      tmp_size = quadtree->n_nodes;
    }

  }

  /* Finalize node */

  _node = quadtree->nodes + _n_nodes;

  for (i = 0; i < 5; i++)
    _node->idx[i] = idx[i];

  for (i = 0; i < 4; i++)
    _node->quadrant_id[i] = quadrant_id[i];

  _node->n_points = _n_points;
}

/*----------------------------------------------------------------------------
 * Build a quadtree structure to locate 2d points in mesh.
 *
 * parameters:
 *   n_points        <-- number of points to locate
 *   point_coords    <-- point coordinates
 *
 * returns:
 *   pointer to local quadtree structure
 *----------------------------------------------------------------------------*/

static _quadtree_t
_build_quadtree(fvmc_lnum_t         n_points,
                const fvmc_coord_t  point_coords[])
{
  size_t i;
  fvmc_lnum_t point_range[2];
  _quadtree_t _quadtree;

  int *point_ids_tmp = NULL;

  /* Initialization */

  point_range[0] = 0;
  point_range[1] = n_points;

  _quadtree.n_points = n_points;
  _quadtree.n_nodes = 0;
  _quadtree.n_nodes_max = 0;
  _quadtree.n_static_subdivision = 0;
  _quadtree.nodes = NULL;
  _quadtree.point_ids = NULL;

  if (n_points > 0) {

    _point_extents(2,
                   n_points,
                   NULL,
                   point_coords,
                   _quadtree.extents);

    BFTC_MALLOC(_quadtree.point_ids, _quadtree.n_points, fvmc_lnum_t);

    for (i = 0; i < _quadtree.n_points; i++)
      _quadtree.point_ids[i] = i;

    BFTC_MALLOC(point_ids_tmp, n_points, int);

    _build_quadtree_leaves(_quadtree.extents,
                           point_coords,
                           point_ids_tmp,
                           &_quadtree,
                           point_range);

    BFTC_FREE(point_ids_tmp);

  }

  return _quadtree;
}

/*----------------------------------------------------------------------------
 * Free a quadtree structure.
 *
 * parameters:
 *   quadtree <-> quadtree structure whose elements are to be freed
 *
 * returns:
 *   pointer to local quadtree structure
 *----------------------------------------------------------------------------*/

static void
_free_quadtree(_quadtree_t *quadtree)
{

  quadtree->n_points = 0;
  quadtree->n_nodes = 0;
  quadtree->n_nodes_max = 0;

  BFTC_FREE(quadtree->nodes);
  BFTC_FREE(quadtree->point_ids);
}

/*----------------------------------------------------------------------------
 * locate points in box defined by extents in a quadrant.
 *
 * parameters:
 *   search_extents  <-- extents associated with element:
 *   point_coords    <-- point coordinates
 *                       x_min, y_min, x_max, y_max (size: 4)
 *   quadtree        <-- point quadtree
 *   node_extents    <-- extents of quadrant
 *   node_id         <-- if of node in octree
 *   loc_point_ids   --> ids of located points (max size: quadtree->n_points)
 *   n_loc_points    --> number of points located
 *----------------------------------------------------------------------------*/

static void
_query_quadtree_node(const double        extents[],
                     const fvmc_coord_t   point_coords[],
                     const _quadtree_t  *quadtree,
                     const double        node_extents[],
                     int                 node_id,
                     fvmc_lnum_t         *loc_point_ids,
                     fvmc_lnum_t         *n_loc_points)
{
  _quadrant_t *node;
  int i, j, k;
  double sub_extents[4], mid[2];

  node = quadtree->nodes + node_id;

  if (_intersect_extents(2, node_extents, extents)) {

    for (j = 0; j < 2; j++)
      mid[j]= (node_extents[j] + node_extents[j + 2]) * 0.5;

    /* Loop on node leaves */

    for (i = 0; i < 4; i++) {

      /* Compute quadrant extents */

      if (i < 2) {
        sub_extents[0] = node_extents[0];
        sub_extents[2] = mid[0];
      }
      else {
        sub_extents[0] = mid[0];
        sub_extents[2] = node_extents[2];
      }

      if (i%2 < 1) {
        sub_extents[1] = node_extents[1];
        sub_extents[3] = mid[1];
      }
      else {
        sub_extents[1] = mid[1];
        sub_extents[3] = node_extents[3];
      }

      /* Search recursively if quadrant is not final */

      if (node->quadrant_id[i] > -1)
        _query_quadtree_node(extents,
                             point_coords,
                             quadtree,
                             sub_extents,
                             node->quadrant_id[i],
                             loc_point_ids,
                             n_loc_points);

      else {

        /* Update list of points located */

        if (_intersect_extents(2, sub_extents, extents)) {

          for (k = node->idx[i]; k < node->idx[i+1]; k++) {

            fvmc_lnum_t point_id = quadtree->point_ids[k];

            if (_within_extents(2, point_coords + point_id*2, extents)) {
              loc_point_ids[*n_loc_points] = point_id;
              (*n_loc_points)++;
            }

          }
        }

      }

    } /* End of loop on node leaves */

  }

}

/*----------------------------------------------------------------------------
 * locate points in box defined by extents using an octree.
 *
 * parameters:
 *   search_extents  <-- extents associated with element:
 *   point_coords    <-- point coordinates
 *                       x_min, y_min, x_max, y_max (size: 4)
 *   quadtree        <-- point quadtree
 *   n_loc_points    --> number of points located
 *   loc_point_ids   --> ids of located points (max size: octree->n_points)
 *----------------------------------------------------------------------------*/

static void
_query_quadtree(const double        extents[],
                const fvmc_coord_t   point_coords[],
                const _quadtree_t  *quadtree,
                fvmc_lnum_t         *n_loc_points,
                fvmc_lnum_t          loc_point_ids[])
{
  *n_loc_points = 0;

  if (quadtree->n_points > 0)
    _query_quadtree_node(extents,
                         point_coords,
                         quadtree,
                         quadtree->extents,
                         0,
                         loc_point_ids,
                         n_loc_points);
}

/*----------------------------------------------------------------------------
 * Locate points on a 3d edge, and update the location[] and distance[]
 * arrays associated with the point set.
 *
 * parameters:
 *   elt_num             <-- number of element corresponding to extents
 *   element_vertex_num  <-- element vertex numbers
 *   parent_vertex_num   <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords       <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extents <-- number of points in extents
 *   points_in_extents   <-- ids of points in extents
 *   tolerance           <-- associated tolerance (considered infinite if < 0)
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, absolute distance
 *                           to element if located (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_on_edge_3d(fvmc_lnum_t           elt_num,
                   const fvmc_lnum_t     element_vertex_num[],
                   const fvmc_lnum_t    *parent_vertex_num,
                   const fvmc_coord_t    vertex_coords[],
                   const fvmc_coord_t    point_coords[],
                   fvmc_lnum_t           n_points_in_extents,
                   const fvmc_lnum_t     points_in_extents[],
                   double               tolerance,
                   fvmc_lnum_t           location[],
                   float                distance[])
{
  fvmc_lnum_t  i, j, k, coord_idx_0, coord_idx_1;

  double u[3], v[3];
  double uv, len2, isop_0;
  double dist2, epsilon2, vertex_dist2;

  /* vertex index of the edge studied */

  if (parent_vertex_num == NULL) {
    coord_idx_0 = element_vertex_num[0] - 1;
    coord_idx_1 = element_vertex_num[1] - 1;
  }
  else {
    coord_idx_0 = parent_vertex_num[element_vertex_num[0] - 1] - 1;
    coord_idx_1 = parent_vertex_num[element_vertex_num[1] - 1] - 1;
  }

  /* Calculate edge vector and length */

  for (j = 0; j < 3; j++)
    u[j] =   vertex_coords[(coord_idx_1*3) + j]
           - vertex_coords[(coord_idx_0*3) + j];

  len2 = _DOT_PRODUCT(u, u);

  if (len2 < _epsilon_denom){
    bftc_printf("warning _locate_on_edge_3d : Reduce _epsilon_denom criteria : %12.5e < %12.5e\n", len2, _epsilon_denom);
    bftc_printf_flush();
    return;
  }

  else if (tolerance < 0.0)
    epsilon2 = HUGE_VAL;

  else
    epsilon2 = len2*tolerance*tolerance;

  /* Loop on points resulting from extent query */

  for (k = 0; k < n_points_in_extents; k++) {

    i =  points_in_extents[k];

    vertex_dist2 = distance[i]*distance[i];

    /* Calculate linear coordinates of projection of point on edge axis */

    for (j = 0; j < 3; j++)
      v[j] = point_coords[i*3 + j] - vertex_coords[(coord_idx_0*3) + j];

    uv = _DOT_PRODUCT(u, v);

    isop_0 = uv / len2;

    /* Set v to be the vector from the point to the closest point on
       the segment (if isop_0 < 0, v is already that vector) */

    if (isop_0 >= 1) {
      for (j = 0; j < 3; j++)
        v[j] = point_coords[i*3 + j] - vertex_coords[(coord_idx_1*3) + j];
    }
    else if (isop_0 > 0) {
      for (j = 0; j < 3; j++)
        v[j] -= isop_0*u[j];
    }

    /* Distance between point to locate and its projection */

    dist2 = _DOT_PRODUCT_2D(v, v);

    if (dist2 < epsilon2 && (dist2 < vertex_dist2 || distance[i] < 0.0)) {
      location[i] = elt_num;
      distance[i] = sqrt(dist2);
    }

  } /* End of loop on points resulting from extent query */

}

/*----------------------------------------------------------------------------
 * Locate points on a 2d edge, and update the location[] and distance[]
 * arrays associated with the point set.
 *
 * parameters:
 *   elt_num             <-- number of element corresponding to extents
 *   element_vertex_num  <-- element vertex numbers
 *   parent_vertex_num   <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords       <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extent  <-- number of points in extents
 *   points_in_extent    <-- ids of points in extents
 *   n_points_in_extent  <-- number of points in extents
 *   points_in_extent    <-- ids of points in extents
 *   tolerance           <-- associated tolerance (considered infinite if < 0)
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, absolute distance
 *                           to element if located (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_on_edge_2d(fvmc_lnum_t           elt_num,
                   const fvmc_lnum_t     element_vertex_num[],
                   const fvmc_lnum_t    *parent_vertex_num,
                   const fvmc_coord_t    vertex_coords[],
                   const fvmc_coord_t    point_coords[],
                   fvmc_lnum_t           n_points_in_extents,
                   const fvmc_lnum_t     points_in_extents[],
                   double               tolerance,
                   fvmc_lnum_t           location[],
                   float                distance[])
{
  fvmc_lnum_t  i, j, k, coord_idx_0, coord_idx_1;

  double u[2], v[2];
  double uv, len2, isop_0;
  double dist2, epsilon2, vertex_dist2;

  /* vertex index of the edge studied */

  if (parent_vertex_num == NULL) {
    coord_idx_0 = element_vertex_num[0] - 1;
    coord_idx_1 = element_vertex_num[1] - 1;
  }
  else {
    coord_idx_0 = parent_vertex_num[element_vertex_num[0] - 1] - 1;
    coord_idx_1 = parent_vertex_num[element_vertex_num[1] - 1] - 1;
  }

  /* Calculate edge vector and length */

  for (j = 0; j < 2; j++)
    u[j] =   vertex_coords[(coord_idx_1*2) + j]
           - vertex_coords[(coord_idx_0*2) + j];

  len2 = _DOT_PRODUCT_2D(u, u);

  if (len2 < _epsilon_denom){
    bftc_printf("warning _locate_on_edge_2d : Reduce _epsilon_denom criteria : %12.5e < %12.5e\n", len2, _epsilon_denom);
    bftc_printf_flush();
    return;
  }

  else if (tolerance < 0.0)
    epsilon2 = HUGE_VAL;

  else
    epsilon2 = len2*tolerance*tolerance;

  /* Loop on points resulting from extent query */

  for (k = 0; k < n_points_in_extents; k++) {

    i =  points_in_extents[k];

    vertex_dist2 = distance[i]*distance[i];

    /* Calculate linear coordinates of projection of point on edge axis */

    for (j = 0; j < 2; j++)
      v[j] = point_coords[i*2 + j] - vertex_coords[(coord_idx_0*2) + j];

    uv = u[0]*v[0] + u[1]*v[1];

    isop_0 = uv / len2;

    /* Set v to be the vector from the point to the closest point on
       the segment (if isop_0 < 0, v is already that vector) */

    if (isop_0 >= 1) {
      for (j = 0; j < 2; j++)
        v[j] = point_coords[i*2 + j] - vertex_coords[(coord_idx_1*2) + j];
    }
    else if (isop_0 > 0) {
      for (j = 0; j < 2; j++)
        v[j] -= isop_0*u[j];
    }

    /* Distance between point to locate and its projection */

    dist2 = _DOT_PRODUCT_2D(v, v);

    if (dist2 < epsilon2 && (dist2 < vertex_dist2 || distance[i] < 0.0)) {
      location[i] = elt_num;
      distance[i] = sqrt(dist2);
    }

  } /* End of loop on points resulting from extent query */

}

/*----------------------------------------------------------------------------
 * Locate points in a given set of 3d triangles, and updates the location[]
 * and distance[] arrays associated with the point set.
 *
 * This function is called for sets of triangles belonging to the subdivision
 * of a given 3d face. Barycentric coordinates are used to locate the
 * projection of points.
 *
 * parameters:
 *   elt_num             <-- number of element corresponding to extents
 *   n_triangles         <-- number of triangles
 *   triangle_vertices   <-- triangles connectivity; size: 2 * 3
 *   parent_vertex_num   <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords       <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extents <-- number of points in extents
 *   points_in_extents   <-- ids of points in extents
 *   tolerance           <-- associated tolerance (considered infinite if < 0)
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, absolute distance
 *                           to element if located (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_on_triangles_3d(fvmc_lnum_t           elt_num,
                        int                   n_triangles,
                        const fvmc_lnum_t     triangle_vertices[],
                        const fvmc_lnum_t    *parent_vertex_num,
                        const fvmc_coord_t    vertex_coords[],
                        const fvmc_coord_t    point_coords[],
                        fvmc_lnum_t           n_points_in_extents,
                        const fvmc_lnum_t     points_in_extents[],
                        const double         tolerance,
                        fvmc_lnum_t           location[],
                        float                distance[])
{
  fvmc_lnum_t  i, j, k, tria_id, coord_idx_0, coord_idx_1, coord_idx_2;

  double t[3], u[3], v[3], w[3], vect_tmp[3];
  double uu, vv, uv, ut, vt, ww, det, tmp_max;
  double epsilon2, dist2, vertex_dist2, isop_0, isop_1;

  double tolerance2 = tolerance*tolerance;

  /* Loop on element's sub-triangles */

  for (tria_id = 0; tria_id < n_triangles; tria_id++) {

    /* vertex index of the triangle studied */

    if (parent_vertex_num == NULL) {
      coord_idx_0 = triangle_vertices[tria_id*3]     - 1;
      coord_idx_1 = triangle_vertices[tria_id*3 + 1] - 1;
      coord_idx_2 = triangle_vertices[tria_id*3 + 2] - 1;
    }
    else {
      coord_idx_0 = parent_vertex_num[triangle_vertices[tria_id*3]    - 1] - 1;
      coord_idx_1 = parent_vertex_num[triangle_vertices[tria_id*3+ 1] - 1] - 1;
      coord_idx_2 = parent_vertex_num[triangle_vertices[tria_id*3+ 2] - 1] - 1;
    }

    /* Calculate triangle-constant values for barycentric coordinates */

    for (j = 0; j < 3; j++) {
      u[j] = - vertex_coords[(coord_idx_0*3) + j]
             + vertex_coords[(coord_idx_1*3) + j];
      v[j] = - vertex_coords[(coord_idx_0*3) + j]
             + vertex_coords[(coord_idx_2*3) + j];
      w[j] =   vertex_coords[(coord_idx_1*3) + j]
             - vertex_coords[(coord_idx_2*3) + j];
    }

    uu = _DOT_PRODUCT(u, u);
    vv = _DOT_PRODUCT(v, v);
    ww = _DOT_PRODUCT(w, w);
    uv = _DOT_PRODUCT(u, v);

    det = (uu*vv - uv*uv);

    if (det < _epsilon_denom){
      bftc_printf("warning _locate_on_triangles_2d : Reduce _epsilon_denom criteria : %12.5e < %12.5e\n", det, _epsilon_denom);
      bftc_printf_flush();
      continue;
    }

    /* epsilon2 is based on maximum edge length (squared) */

    tmp_max = FVMC_MAX(vv, ww);

    if (tolerance < 0.)
      epsilon2 = HUGE_VAL;
    else
      epsilon2 = FVMC_MAX(uu, tmp_max) * tolerance2;

    /* Loop on points resulting from extent query */

    for (k = 0; k < n_points_in_extents; k++) {

      i =  points_in_extents[k];

      vertex_dist2 = distance[i]*distance[i];

      /* Calculation of the barycenter coordinates for the projected node */

      for (j = 0; j < 3; j++)
        t[j] = - vertex_coords[(coord_idx_0*3) + j]
               + point_coords[i*3 + j];

      ut = _DOT_PRODUCT(u, t);
      vt = _DOT_PRODUCT(v, t);

      isop_0 = (ut*vv - vt*uv) / det;
      isop_1 = (uu*vt - uv*ut) / det;

      _CROSS_PRODUCT(vect_tmp, u, v);

      /* if the projected point is not on triangle, we project it
         on the nearest edge or node */

      if (isop_0 < 0.)
        isop_0 = 0.;

      if (isop_1 < 0.)
        isop_1 = 0.;

      if ((1.0 - isop_0 - isop_1) < 0.) {
        isop_0 = isop_0 / (isop_0 + isop_1);
        isop_1 = isop_1 / (isop_0 + isop_1);
      }

      /* re-use vect_tmp */

      for (j = 0; j < 3; j++)
        vect_tmp[j] =   vertex_coords[coord_idx_0*3 + j]
                      + u[j]*isop_0
                      + v[j]*isop_1
                      - point_coords[i*3 + j];

      /* distance between point to locate and its projection */

      dist2 = _DOT_PRODUCT(vect_tmp, vect_tmp);

      if (dist2 < epsilon2 && (dist2 < vertex_dist2 || distance[i] < 0.0)) {
        location[i] = elt_num;
        distance[i] = sqrt(dist2);
      }

    } /* End of loop on points resulting from extent query */

  } /* End of loop on element's sub-triangles */

}

/*----------------------------------------------------------------------------
 * Locate points in a given set of 2d triangles, and updates the location[]
 * and distance[] arrays associated with the point set.
 *
 * This function is called for sets of triangles belonging to the subdivision
 * of a given 2d face. Barycentric coordinates are used to locate the
 * projection of points.
 *
 * parameters:
 *   elt_num             <-- number of element corresponding to extents
 *   n_triangles         <-- number of triangles
 *   triangle_vertices   <-- triangles connectivity; size: 2 * 3
 *   parent_vertex_num   <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords       <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extent  <-- number of points in extents
 *   points_in_extent    <-- ids of points in extents
 *   n_points_in_extent  <-- number of points in extents
 *   points_in_extent    <-- ids of points in extents
 *   tolerance           <-- associated tolerance
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, 0 - 1 if inside,
 *                           > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_on_triangles_2d(fvmc_lnum_t           elt_num,
                        int                  n_triangles,
                        const fvmc_lnum_t     triangle_vertices[],
                        const fvmc_lnum_t    *parent_vertex_num,
                        const fvmc_coord_t    vertex_coords[],
                        const fvmc_coord_t    point_coords[],
                        fvmc_lnum_t           n_points_in_extents,
                        const fvmc_lnum_t     points_in_extents[],
                        double               tolerance,
                        fvmc_lnum_t           location[],
                        float                distance[])
{
  fvmc_lnum_t  i, j, k, tria_id, coord_idx_0, coord_idx_1, coord_idx_2;

  double t[2], u[2], v[2], shapef[3];
  double uu, vv, uv, ut, vt, det;
  double dist, max_dist, isop_0, isop_1;

  /* Loop on element's sub-triangles */

  for (tria_id = 0; tria_id < n_triangles; tria_id++) {

    /* vertex index of the triangle studied */

    if (parent_vertex_num == NULL) {
      coord_idx_0 = triangle_vertices[tria_id*3]     - 1;
      coord_idx_1 = triangle_vertices[tria_id*3 + 1] - 1;
      coord_idx_2 = triangle_vertices[tria_id*3 + 2] - 1;
    }
    else {
      coord_idx_0 = parent_vertex_num[triangle_vertices[tria_id*3]    - 1] - 1;
      coord_idx_1 = parent_vertex_num[triangle_vertices[tria_id*3+ 1] - 1] - 1;
      coord_idx_2 = parent_vertex_num[triangle_vertices[tria_id*3+ 2] - 1] - 1;
    }

    /* Calculate triangle-constant values for barycentric coordinates */

    for (j = 0; j < 2; j++) {
      u[j] = - vertex_coords[(coord_idx_0*2) + j]
             + vertex_coords[(coord_idx_1*2) + j];
      v[j] = - vertex_coords[(coord_idx_0*2) + j]
             + vertex_coords[(coord_idx_2*2) + j];
    }

    uu = _DOT_PRODUCT_2D(u, u);
    vv = _DOT_PRODUCT_2D(v, v);
    uv = _DOT_PRODUCT_2D(u, v);

    det = (uu*vv - uv*uv);

    if (det < _epsilon_denom) {
      bftc_printf("warning _locate_on_triangles_2d : Reduce _epsilon_denom criteria : %12.5e < %12.5e\n", det, _epsilon_denom);
      bftc_printf_flush();
      continue;
    }

    /* Loop on points resulting from extent query */

    for (k = 0; k < n_points_in_extents; k++) {

      i =  points_in_extents[k];

      /* Calculation of the barycenter coordinates for the projected node */

      for (j = 0; j < 2; j++)
        t[j] = - vertex_coords[(coord_idx_0*2) + j]
               + point_coords[i*2 + j];

      ut = u[0]*t[0] + u[1]*t[1];
      vt = v[0]*t[0] + v[1]*t[1];

      isop_0 = (ut*vv - vt*uv) / det;
      isop_1 = (uu*vt - uv*ut) / det;

      shapef[0] = 1. - isop_0 - isop_1;
      shapef[1] =      isop_0;
      shapef[2] =               isop_1;

      max_dist = -1.0;

      for (j = 0; j < 3; j++){

        dist = 2.*FVMC_ABS(shapef[j] - 0.5);

        if (max_dist < dist)
          max_dist = dist;
      }

      if (   (max_dist > -0.5 && max_dist < (1. + 2.*tolerance))
          && (max_dist < distance[i] || distance[i] < 0)) {
        location[i] = elt_num;
        distance[i] = max_dist;
      }

    } /* End of loop on points resulting from extent query */

  } /* End of loop on element's sub-triangles */

}

/*----------------------------------------------------------------------------
 * Locate points in a tetrahedron whose coordinates are pre-computed,
 * updating the location[] and distance[] arrays associated with a set
 * of points.
 *
 * parameters:
 *   elt_num             <-- element number
 *   tetra_coords[]      <-- tetrahedra vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extents <-- number of points in element extents
 *   points_in_extents   <-- ids of points in extents
 *   tolerance           <-- associated tolerance
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, 0 - 1 if inside,
 *                           > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_in_tetra(fvmc_lnum_t         elt_num,
                 fvmc_coord_t        tetra_coords[4][3],
                 const fvmc_coord_t  point_coords[],
                 fvmc_lnum_t         n_points_in_extents,
                 const fvmc_lnum_t   points_in_extents[],
                 double             tolerance,
                 fvmc_lnum_t         location[],
                 float              distance[])
{
  double vol6;
  double dist, max_dist;
  int i, j, k;

  double isop_0, isop_1, isop_2;
  double t00, t10, t20, t01, t02, t03, t11, t12, t13, t21, t22, t23;
  double v01[3], v02[3], v03[3], shapef[4];

  for (i = 0; i < 3; i++) {
    v01[i] = tetra_coords[1][i] - tetra_coords[0][i];
    v02[i] = tetra_coords[2][i] - tetra_coords[0][i];
    v03[i] = tetra_coords[3][i] - tetra_coords[0][i];
  }

  vol6 = fabs(  v01[0] * (v02[1]*v03[2] - v02[2]*v03[1])
              - v02[0] * (v01[1]*v03[2] - v01[2]*v03[1])
              + v03[0] * (v01[1]*v02[2] - v01[2]*v02[1]));

  if (vol6 < _epsilon_denom){
    bftc_printf("warning _locate_in_tetra : Reduce _epsilon_denom criteria : %12.5e < %12.5e\n", vol6, _epsilon_denom);
    bftc_printf_flush();
    return;
  }

  for (k = 0; k < n_points_in_extents; k++) {

    i = points_in_extents[k];

    t00  =   point_coords[i*3]     - tetra_coords[0][0];
    t10  =   point_coords[i*3 + 1] - tetra_coords[0][1];
    t20  =   point_coords[i*3 + 2] - tetra_coords[0][2];

    t01  = - tetra_coords[0][0] + tetra_coords[1][0];
    t02  = - tetra_coords[0][0] + tetra_coords[2][0];
    t03  = - tetra_coords[0][0] + tetra_coords[3][0];

    t11  = - tetra_coords[0][1] + tetra_coords[1][1];
    t12  = - tetra_coords[0][1] + tetra_coords[2][1];
    t13  = - tetra_coords[0][1] + tetra_coords[3][1];

    t21  = - tetra_coords[0][2] + tetra_coords[1][2];
    t22  = - tetra_coords[0][2] + tetra_coords[2][2];
    t23  = - tetra_coords[0][2] + tetra_coords[3][2];

    isop_0 = (  t00 * (t12*t23 - t13*t22)
              - t10 * (t02*t23 - t22*t03)
              + t20 * (t02*t13 - t12*t03)) / vol6;
    isop_1 = (- t00 * (t11*t23 - t13*t21)
              + t10 * (t01*t23 - t21*t03)
              - t20 * (t01*t13 - t03*t11)) / vol6;
    isop_2 = (  t00 * (t11*t22 - t21*t12)
              - t10 * (t01*t22 - t21*t02)
              + t20 * (t01*t12 - t11*t02)) / vol6;

    shapef[0] = 1. - isop_0 - isop_1 - isop_2;
    shapef[1] =      isop_0;
    shapef[2] =               isop_1;
    shapef[3] =                        isop_2;

    max_dist = -1.0;

    for (j = 0; j < 4; j++){

      dist = 2.*FVMC_ABS(shapef[j] - 0.5);

      if (max_dist < dist)
        max_dist = dist;
    }

    if (   (max_dist > -0.5 && max_dist < (1. + 2.*tolerance))
        && (max_dist < distance[i] || distance[i] < 0)) {
      location[i] = elt_num;
      distance[i] = max_dist;
    }

  }

}

/*---------------------------------------------------------------------------
 * Solve the equation "matrix.x = b" with Cramer's rule.
 *
 * parameters:
 *   m[3][3] <-- equation matrix
 *   b[3]    <-- b equation right hand side
 *   x[3]    <-> equation solution (unchanged if matrix is singular)
 *
 * returns:
 *   1 if matrix is singular, 0 otherwise
 *----------------------------------------------------------------------------*/

static int
_inverse_3x3(double  m[3][3],
             double  b[3],
             double  x[3])
{
  double det, det_inv, x0, x1, x2;

  det =   m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
        - m[1][0]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
        + m[2][0]*(m[0][1]*m[1][2] - m[1][1]*m[0][2]);

  if (FVMC_ABS(det) < _epsilon_denom)
    return 1;
  else
    det_inv = 1./det;

  /* Use local variables to ensure no aliasing */

  x0 = (  b[0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
        - b[1]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
        + b[2]*(m[0][1]*m[1][2] - m[1][1]*m[0][2])) * det_inv;

  x1 = (  m[0][0]*(b[1]*m[2][2] - b[2]*m[1][2])
        - m[1][0]*(b[0]*m[2][2] - b[2]*m[0][2])
        + m[2][0]*(b[0]*m[1][2] - b[1]*m[0][2])) * det_inv;

  x2 = (  m[0][0]*(m[1][1]*b[2] - m[2][1]*b[1])
        - m[1][0]*(m[0][1]*b[2] - m[2][1]*b[0])
        + m[2][0]*(m[0][1]*b[1] - m[1][1]*b[0])) * det_inv;

  /* Copy local variables to output */

  x[0] = x0; x[1] = x1; x[2] = x2;

  return 0;
}

/*----------------------------------------------------------------------------
 * Compute 3d shape functions and their derivatives given element
 * parametric coordinates.
 *
 * This function is adapted from the CGNS interpolation tool.
 *
 * parameters:
 *   elt_type    <-- type of element
 *   uvw[]       <-- parametric coordinates
 *   shapef[]    --> barycenter's coordinates
 *   deriv [][]  --> derivative of shape function
*----------------------------------------------------------------------------*/

static void
_compute_shapef_3d(fvmc_element_t  elt_type,
                   const double   uvw[3],
                   double         shapef[8],
                   double         deriv[8][3])

{
  switch (elt_type) {

  case FVMC_CELL_HEXA:

    shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[2] = uvw[0] * uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0]) * uvw[1] * (1.0 - uvw[2]);
    shapef[4] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * uvw[2];
    shapef[5] = uvw[0] * (1.0 - uvw[1]) * uvw[2];
    shapef[6] = uvw[0] * uvw[1] * uvw[2];
    shapef[7] = (1.0 - uvw[0]) * uvw[1] * uvw[2];

    if (deriv != NULL) {
      deriv[0][0] = -(1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[1][1] = -uvw[0] * (1.0 - uvw[2]);
      deriv[1][2] = -uvw[0] * (1.0 - uvw[1]);
      deriv[2][0] =  uvw[1] * (1.0 - uvw[2]);
      deriv[2][1] =  uvw[0] * (1.0 - uvw[2]);
      deriv[2][2] = -uvw[0] * uvw[1];
      deriv[3][0] = -uvw[1] * (1.0 - uvw[2]);
      deriv[3][1] =  (1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[3][2] = -(1.0 - uvw[0]) * uvw[1];
      deriv[4][0] = -(1.0 - uvw[1]) * uvw[2];
      deriv[4][1] = -(1.0 - uvw[0]) * uvw[2];
      deriv[4][2] =  (1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[5][0] =  (1.0 - uvw[1]) * uvw[2];
      deriv[5][1] = -uvw[0] * uvw[2];
      deriv[5][2] =  uvw[0] * (1.0 - uvw[1]);
      deriv[6][0] =  uvw[1] * uvw[2];
      deriv[6][1] =  uvw[0] * uvw[2];
      deriv[6][2] =  uvw[0] * uvw[1];
      deriv[7][0] = -uvw[1] * uvw[2];
      deriv[7][1] =  (1.0 - uvw[0]) * uvw[2];
      deriv[7][2] =  (1.0 - uvw[0]) * uvw[1];
    }

    break;

  case FVMC_CELL_PRISM:

    shapef[0] = (1.0 - uvw[0] - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[2]);
    shapef[2] = uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0] - uvw[1]) * uvw[2];
    shapef[4] = uvw[0] * uvw[2];
    shapef[5] = uvw[1] * uvw[2];

    if (deriv != NULL) {
      deriv[0][0] = -(1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0] - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[2]);
      deriv[1][1] =  0.0;
      deriv[1][2] = -uvw[0];
      deriv[2][0] =  0.0;
      deriv[2][1] =  (1.0 - uvw[2]);
      deriv[2][2] = -uvw[1];
      deriv[3][0] = -uvw[2];
      deriv[3][1] = -uvw[2];
      deriv[3][2] =  (1.0 - uvw[0] - uvw[1]);
      deriv[4][0] =  uvw[2];
      deriv[4][1] =  0.0;
      deriv[4][2] =  uvw[0];
      deriv[5][0] =  0.0;
      deriv[5][1] =  uvw[2];
      deriv[5][2] =  uvw[1];
    }

    break;

  case FVMC_CELL_PYRAM:

    shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[2] = uvw[0] * uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0]) * uvw[1] * (1.0 - uvw[2]);
    shapef[4] = uvw[2];

    if (deriv != NULL) {
      deriv[0][0] = -(1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[1][1] = -uvw[0] * (1.0 - uvw[2]);
      deriv[1][2] = -uvw[0] * (1.0 - uvw[1]);
      deriv[2][0] =  uvw[1] * (1.0 - uvw[2]);
      deriv[2][1] =  uvw[0] * (1.0 - uvw[2]);
      deriv[2][2] = -uvw[0] * uvw[1];
      deriv[3][0] = -uvw[1] * (1.0 - uvw[2]);
      deriv[3][1] =  (1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[3][2] = -(1.0 - uvw[0]) * uvw[1];
      deriv[4][0] =  0.0;
      deriv[4][1] =  0.0;
      deriv[4][2] =  1.0;
    }

    break;

  default:
    bftc_error(__FILE__, __LINE__, 0,
              _("_compute_shapef: unhandled element type %s\n"),
              fvmc_element_type_name[elt_type]);

  }

}

/*----------------------------------------------------------------------------
 * Compute hexahedron, pyramid, or prism parametric coordinates for a
 * given point.
 *
 * This function is adapted from the CGNS interpolation tool.
 *
 * parameters:
 *   elt_type            <-- type of element
 *   point_coords        <-- point coordinates
 *   vertex_coords[]     <-- pointer to element vertex coordinates
 *   tolerance           <-- location tolerance factor
 *   uvw[]               --> parametric coordinates of point in element
*----------------------------------------------------------------------------*/
static int
_compute_uvw(fvmc_element_t       elt_type,
             const fvmc_coord_t   point_coords[],
             double              vertex_coords[8][3],
             double              tolerance,
             double              uvw[3])

{
  int i, j, n_elt_vertices, iter;
  int max_iter = 20;
  double dist;
  double a[3][3], b[3], x[3], shapef[8], dw[8][3];

  n_elt_vertices = fvmc_nodal_n_vertices_element[elt_type];

  assert(   elt_type == FVMC_CELL_HEXA
         || elt_type == FVMC_CELL_PRISM
         || elt_type == FVMC_CELL_PYRAM);

  /* Use Newton-method to determine parametric coordinates and shape function */

  for (i = 0; i < 3; i++)
    uvw[i] = 0.5;

  for (iter = 0; iter < max_iter; iter++) {

    _compute_shapef_3d(elt_type, uvw, shapef, dw);

    b[0] = - point_coords[0];
    b[1] = - point_coords[1];
    b[2] = - point_coords[2];

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++)
        a[i][j] = 0.0;
    }

    for (i = 0; i < n_elt_vertices; i++) {

      b[0] += (shapef[i] * vertex_coords[i][0]);
      b[1] += (shapef[i] * vertex_coords[i][1]);
      b[2] += (shapef[i] * vertex_coords[i][2]);

      for (j = 0; j < 3; j++) {
        a[0][j] -= (dw[i][j] * vertex_coords[i][0]);
        a[1][j] -= (dw[i][j] * vertex_coords[i][1]);
        a[2][j] -= (dw[i][j] * vertex_coords[i][2]);
      }

    }

    if (_inverse_3x3(a, b, x))
      return 0;

    dist = 0.0;

    for (i = 0; i < 3; i++) {
      dist += x[i] * x[i];
      uvw[i] += x[i];
    }

    if (dist <= (tolerance * tolerance))
      return 1;

  }

  return 0;
}

/*----------------------------------------------------------------------------
 * Locate points in a given 3d cell (other than tetrahedra or polyhedra,
 * handlesd elsewhere), updating the location[] and distance[] arrays
 * associated with a set of points.
 *
 * parameters:
 *   elt_num             <-- element number
 *   elt_type            <-- type of element
 *   element_vertex_num  <-- element vertex numbers
 *   parent_vertex_num   <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords[]     <-- pointer to vertex coordinates
 *   point_coords        <-- point coordinates
 *   n_points_in_extents <-- number of points in element extents
 *   points_in_extents   <-- ids of points in extents
 *   tolerance           <-- associated tolerance
 *   location            <-> number of element containing or closest to each
 *                           point (size: n_points)
 *   distance            <-> distance from point to element indicated by
 *                           location[]: < 0 if unlocated, 0 - 1 if inside,
 *                           > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_locate_in_cell_3d(fvmc_lnum_t          elt_num,
                   fvmc_element_t       elt_type,
                   const fvmc_lnum_t    element_vertex_num[],
                   const fvmc_lnum_t   *parent_vertex_num,
                   const fvmc_coord_t   vertex_coords[],
                   const fvmc_coord_t   point_coords[],
                   fvmc_lnum_t          n_points_in_extents,
                   const fvmc_lnum_t    points_in_extents[],
                   double              tolerance,
                   fvmc_lnum_t          location[],
                   float               distance[])
{
  int i, j, k, n_vertices;
  fvmc_lnum_t coord_idx, vertex_id;

  double uvw[3], dist, shapef[8],max_dist;
  double  _vertex_coords[8][3];

  n_vertices = fvmc_nodal_n_vertices_element[elt_type];

  /* Initialize local element coordinates copy */

  for (vertex_id = 0; vertex_id < n_vertices; vertex_id++) {

    if (parent_vertex_num == NULL)
      coord_idx = element_vertex_num[vertex_id] -1;
    else
      coord_idx = parent_vertex_num[element_vertex_num[vertex_id] - 1] - 1;

    for (j = 0; j < 3; j++)
      _vertex_coords[vertex_id][j] = vertex_coords[(coord_idx * 3) + j];

  }

  /* Shape functions may be computed directly with tetrahedra */

  if (elt_type == FVMC_CELL_TETRA)

    _locate_in_tetra(elt_num,
                     _vertex_coords,
                     point_coords,
                     n_points_in_extents,
                     points_in_extents,
                     tolerance,
                     location,
                     distance);

  /* For cell shapes other than tetrahedra, find shape functions iteratively */

  else {

    for (k = 0; k < n_points_in_extents; k++) {

      i = points_in_extents[k];

      if (_compute_uvw(elt_type,
                       point_coords + 3*i,
                       _vertex_coords,
                       tolerance,
                       uvw)) {

        max_dist = -1.0;

        /* For hexahedra, no need to compute shape functions, as
           the 3 parametric coordinates are simpler to use */

        if (elt_type == FVMC_CELL_HEXA) {

          for (j = 0; j < 3; j++){

            dist = 2.*FVMC_ABS(uvw[j] - 0.5);

            if (max_dist < dist)
              max_dist = dist;
          }

        }

        /* For pyramids ands prisms, we need to compute shape functions */

        else {

          _compute_shapef_3d(elt_type, uvw, shapef, NULL);

          for (j = 0; j < n_vertices; j++){

            dist = 2.*FVMC_ABS(shapef[j] - 0.5);

          if (max_dist < dist)
            max_dist = dist;
          }

        }

        /* For all element types, update location and distance arrays */

        if (   (max_dist > -0.5 && max_dist < (1. + 2.*tolerance))
            && (max_dist < distance[i] || distance[i] < 0)) {
          location[i] = elt_num;
          distance[i] = max_dist;
        }

      }

    } /* End of loop on points in extents */

  }

}

/*----------------------------------------------------------------------------
 * Find elements in a given polyhedral section containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * Space and element dimensions are both equal to 3 here.
 *
 * parameters:
 *   this_section      <-- pointer to mesh section representation structure
 *   parent_vertex_num <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   tolerance         <-- associated tolerance
 *   base_element_num  <-- < 0 for location relative to parent element numbers,
 *                         number of elements in preceding sections of same
 *                         element dimension + 1 otherwise
 *   point_coords      <-- point coordinates
 *   octree            <-- point octree
 *   points_in_extents <-> array for query of ids of points in extents
 *                         (size: n_points, less usually needed)
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, 0 - 1 if inside,
 *                         > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_polyhedra_section_locate(const fvmc_nodal_section_t  *this_section,
                          const fvmc_lnum_t           *parent_vertex_num,
                          const fvmc_coord_t           vertex_coords[],
                          double                      tolerance,
                          fvmc_lnum_t                  base_element_num,
                          const fvmc_coord_t           point_coords[],
                          _octree_t                  *octree,
                          fvmc_lnum_t                  points_in_extents[],
                          fvmc_lnum_t                  location[],
                          float                       distance[])
{
  fvmc_lnum_t  i, j, k, n_vertices, face_id, vertex_id, elt_num;
  fvmc_coord_t  center[3];
  double elt_extents[6];

  double _tolerance = tolerance * 2; /* double tolerance, as polyhedra is
                                        split into tetrahedra, whose extents
                                        are 1/2 the polyhedron extents */

  fvmc_lnum_t n_vertices_max = 0;
  fvmc_lnum_t n_points_in_extents = 0;
  fvmc_lnum_t *triangle_vertices = NULL;
  fvmc_triangulate_state_t *state = NULL;

  assert(this_section->face_index != NULL);

  /* Return immediately if nothing to do for this rank */

  if (this_section->n_elements == 0)
    return;

  /* Counting loop on faces */

  for (i = 0; i < this_section->n_faces; i++) {
    n_vertices =   this_section->vertex_index[i + 1]
                 - this_section->vertex_index[i];
    if (n_vertices > n_vertices_max)
      n_vertices_max = n_vertices;
  }

  if (n_vertices_max < 3)
    return;

  BFTC_MALLOC(triangle_vertices, (n_vertices_max-2)*3, int);
  state = fvmc_triangulate_state_create(n_vertices_max);

  /* Loop on elements */

  for (i = 0; i < this_section->n_elements; i++) {

    _Bool elt_initialized = false;

    /* Compute extents */

    for (j = this_section->face_index[i];
         j < this_section->face_index[i + 1];
         j++) {
      face_id = FVMC_ABS(this_section->face_num[j]) - 1;
      for (k = this_section->vertex_index[face_id];
           k < this_section->vertex_index[face_id + 1];
           k++) {
        vertex_id = this_section->vertex_num[k] - 1;

        _update_elt_extents(3,
                            vertex_id,
                            parent_vertex_num,
                            vertex_coords,
                            elt_extents,
                            &elt_initialized);

      }
    }

    _elt_extents_finalize(3, 3, tolerance, elt_extents);

    if (base_element_num < 0) {
      if (this_section->parent_element_num != NULL)
        elt_num = this_section->parent_element_num[i];
      else
        elt_num = i + 1;
    }
    else
      elt_num = base_element_num + i;

    _query_octree(elt_extents,
                  point_coords,
                  octree,
                  &n_points_in_extents,
                  points_in_extents);

    if (n_points_in_extents < 1)
      continue;

    /* Compute psuedo-element center */

    for (j = 0; j < 3; j++)
      center[j] = (elt_extents[j] + elt_extents[j + 3]) * 0.5;

    /* Loop on element faces */

    for (j = this_section->face_index[i];
         j < this_section->face_index[i + 1];
         j++) {

      fvmc_lnum_t n_triangles;

      const fvmc_lnum_t *_vertex_num;

      face_id = FVMC_ABS(this_section->face_num[j]) - 1;

      n_vertices = (  this_section->vertex_index[face_id + 1]
                    - this_section->vertex_index[face_id]);

      _vertex_num = (  this_section->vertex_num
                     + this_section->vertex_index[face_id]);

      if (n_vertices == 4)

        n_triangles = fvmc_triangulate_quadrangle(3,
                                                 vertex_coords,
                                                 parent_vertex_num,
                                                 _vertex_num,
                                                 triangle_vertices);

      else if (n_vertices > 4)

        n_triangles = fvmc_triangulate_polygon(3,
                                              n_vertices,
                                              vertex_coords,
                                              parent_vertex_num,
                                              _vertex_num,
                                              FVMC_TRIANGULATE_MESH_DEF,
                                              triangle_vertices,
                                              state);

      else { /* n_vertices == 3 */

        n_triangles = 1;
        for (k = 0; k < 3; k++)
          triangle_vertices[k] = _vertex_num[k];

      }

      /* Loop on face triangles so as to loop on tetrahedra
         built by joining face triangles and psuedo-center */

      for (k = 0; k < n_triangles; k++) {

        fvmc_lnum_t l, coord_id[3];
        fvmc_coord_t tetra_coords[4][3];

        if (parent_vertex_num == NULL) {
          coord_id[0] = triangle_vertices[k*3    ] - 1;
          coord_id[1] = triangle_vertices[k*3 + 2] - 1;
          coord_id[2] = triangle_vertices[k*3 + 1] - 1;
        }
        else {
          coord_id[0] = parent_vertex_num[triangle_vertices[k*3    ] - 1] - 1;
          coord_id[1] = parent_vertex_num[triangle_vertices[k*3 + 2] - 1] - 1;
          coord_id[2] = parent_vertex_num[triangle_vertices[k*3 + 1] - 1] - 1;
        }

        for (l = 0; l < 3; l++) {
          tetra_coords[0][l] = vertex_coords[3*coord_id[0] + l];
          tetra_coords[1][l] = vertex_coords[3*coord_id[1] + l];
          tetra_coords[2][l] = vertex_coords[3*coord_id[2] + l];
          tetra_coords[3][l] = center[l];
        }

        _locate_in_tetra(elt_num,
                         tetra_coords,
                         point_coords,
                         n_points_in_extents,
                         points_in_extents,
                         _tolerance,
                         location,
                         distance);

      } /* End of loop on face triangles */

    } /* End of loop on element faces */

    _locate_in_extents(elt_num,
                       3,
                       elt_extents,
                       point_coords,
                       n_points_in_extents,
                       points_in_extents,
                       location,
                       distance);

  } /* End of loop on elements */

  BFTC_FREE(triangle_vertices);
  state = fvmc_triangulate_state_destroy(state);
}

/*----------------------------------------------------------------------------
 * Find elements in a given polygonal section containing 3d points: updates
 * the location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_section      <-- pointer to mesh section representation structure
 *   parent_vertex_num <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   tolerance         <-- associated tolerance
 *   base_element_num  <-- < 0 for location relative to parent element numbers,
 *                         number of elements in preceding sections of same
 *                         element dimension + 1 otherwise
 *   n_points          <-- number of points to locate
 *   point_coords      <-- point coordinates
 *   octree            <-- point octree
 *   points_in_extents <-- array for query of ids of points in extents
 *                         (size: n_points, less usually needed)
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, absolute distance
 *                         to element if located (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_polygons_section_locate_3d(const fvmc_nodal_section_t   *this_section,
                            const fvmc_lnum_t            *parent_vertex_num,
                            const fvmc_coord_t            vertex_coords[],
                            const double                 tolerance,
                            fvmc_lnum_t                   base_element_num,
                            const fvmc_coord_t            point_coords[],
                            _octree_t                   *octree,
                            fvmc_lnum_t                   points_in_extents[],
                            fvmc_lnum_t                   location[],
                            float                        distance[])
{
  fvmc_lnum_t  i, j, n_vertices, vertex_id, elt_num;
  int n_triangles;
  double elt_extents[6];

  int n_vertices_max = 0;
  fvmc_lnum_t n_points_in_extents = 0;

  fvmc_lnum_t *triangle_vertices = NULL;
  fvmc_triangulate_state_t *state = NULL;

  /* Return immediately if nothing to do for this rank */

  if (this_section->n_elements == 0)
    return;

  /* Counting loop on elements */

  for (i = 0; i < this_section->n_elements; i++) {

    n_vertices = (  this_section->vertex_index[i + 1]
                  - this_section->vertex_index[i]);

    if (n_vertices > n_vertices_max)
      n_vertices_max = n_vertices;

  }

  if (n_vertices_max < 3)
    return;

  BFTC_MALLOC(triangle_vertices, (n_vertices_max-2)*3, int);
  state = fvmc_triangulate_state_create(n_vertices_max);

  /* Main loop on elements */

  for (i = 0; i < this_section->n_elements; i++) {

    _Bool elt_initialized = false;

    for (j = this_section->vertex_index[i];
         j < this_section->vertex_index[i + 1];
         j++) {
      vertex_id = this_section->vertex_num[j] - 1;

      _update_elt_extents(3,
                          vertex_id,
                          parent_vertex_num,
                          vertex_coords,
                          elt_extents,
                          &elt_initialized);

    }

    _elt_extents_finalize(3, 2, tolerance, elt_extents);

    if (base_element_num < 0) {
      if (this_section->parent_element_num != NULL)
        elt_num = this_section->parent_element_num[i];
      else
        elt_num = i + 1;
    }
    else
      elt_num = base_element_num + i;

    _query_octree(elt_extents,
                  point_coords,
                  octree,
                  &n_points_in_extents,
                  points_in_extents);

    /* Triangulate polygon */

    n_vertices = (  this_section->vertex_index[i + 1]
                  - this_section->vertex_index[i]);
    vertex_id = this_section->vertex_index[i];

    n_triangles = fvmc_triangulate_polygon(3,
                                          n_vertices,
                                          vertex_coords,
                                          parent_vertex_num,
                                          (  this_section->vertex_num
                                           + vertex_id),
                                          FVMC_TRIANGULATE_MESH_DEF,
                                          triangle_vertices,
                                          state);

    /* Locate on triangulated polygon */

    _locate_on_triangles_3d(elt_num,
                            n_triangles,
                            triangle_vertices,
                            parent_vertex_num,
                            vertex_coords,
                            point_coords,
                            n_points_in_extents,
                            points_in_extents,
                            tolerance,
                            location,
                            distance);

  } /* End of loop on elements */

  BFTC_FREE(triangle_vertices);
  state = fvmc_triangulate_state_destroy(state);
}

/*----------------------------------------------------------------------------
 * Find elements in a given polygonal section closest to 3d points: updates
 * the location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_section      <-- pointer to mesh section representation structure
 *   parent_vertex_num <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   base_element_num  <-- < 0 for location relative to parent element numbers,
 *                         number of elements in preceding sections of same
 *                         element dimension + 1 otherwise
 *   point_coords      <-- point coordinates
 *   n_point_ids       <-- number of points to locate
 *   point_id          <-- ids of points to locate (size: n_points)
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, absolute distance
 *                         to element if located (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_polygons_section_closest_3d(const fvmc_nodal_section_t   *this_section,
                             const fvmc_lnum_t            *parent_vertex_num,
                             const fvmc_coord_t            vertex_coords[],
                             fvmc_lnum_t                   base_element_num,
                             const fvmc_coord_t            point_coords[],
                             fvmc_lnum_t                   n_point_ids,
                             const fvmc_lnum_t             point_id[],
                             fvmc_lnum_t                   location[],
                             float                        distance[])
{
  fvmc_lnum_t  i, n_vertices, vertex_id, elt_num;
  int n_triangles;

  int n_vertices_max = 0;

  fvmc_lnum_t *triangle_vertices = NULL;
  fvmc_triangulate_state_t *state = NULL;

  /* Return immediately if nothing to do for this rank */

  if (this_section->n_elements == 0)
    return;

  /* Counting loop on elements */

  for (i = 0; i < this_section->n_elements; i++) {

    n_vertices = (  this_section->vertex_index[i + 1]
                  - this_section->vertex_index[i]);

    if (n_vertices > n_vertices_max)
      n_vertices_max = n_vertices;

  }

  if (n_vertices_max < 3)
    return;

  BFTC_MALLOC(triangle_vertices, (n_vertices_max-2)*3, int);
  state = fvmc_triangulate_state_create(n_vertices_max);

  /* Main loop on elements */

  for (i = 0; i < this_section->n_elements; i++) {

    if (base_element_num < 0) {
      if (this_section->parent_element_num != NULL)
        elt_num = this_section->parent_element_num[i];
      else
        elt_num = i + 1;
    }
    else
      elt_num = base_element_num + i;

    /* Triangulate polygon */

    n_vertices = (  this_section->vertex_index[i + 1]
                  - this_section->vertex_index[i]);
    vertex_id = this_section->vertex_index[i];

    n_triangles = fvmc_triangulate_polygon(3,
                                          n_vertices,
                                          vertex_coords,
                                          parent_vertex_num,
                                          (  this_section->vertex_num
                                           + vertex_id),
                                          FVMC_TRIANGULATE_MESH_DEF,
                                          triangle_vertices,
                                          state);

    /* Locate on triangulated polygon */

    _locate_on_triangles_3d(elt_num,
                            n_triangles,
                            triangle_vertices,
                            parent_vertex_num,
                            vertex_coords,
                            point_coords,
                            n_point_ids,
                            point_id,
                            -1.,
                            location,
                            distance);

  } /* End of loop on elements */

  BFTC_FREE(triangle_vertices);
  state = fvmc_triangulate_state_destroy(state);
}

/*----------------------------------------------------------------------------
 * Find elements in a given section containing 3d points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_section      <-- pointer to mesh section representation structure
 *   parent_vertex_num <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   tolerance         <-- associated tolerance
 *   base_element_num  <-- < 0 for location relative to parent element numbers,
 *                         number of elements in preceding sections of same
 *                         element dimension + 1 otherwise
 *   point_coords      <-- point coordinates
 *   octree            <-- point octree
 *   points_in_extents <-- array for query of ids of points in extents
 *                         (size: octree->n_points, less usually needed)
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated; 0 - 1 if inside,
 *                         and > 1 if outside a volume element, or absolute
 *                         distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_nodal_section_locate_3d(const fvmc_nodal_section_t  *this_section,
                         const fvmc_lnum_t           *parent_vertex_num,
                         const fvmc_coord_t           vertex_coords[],
                         double                      tolerance,
                         fvmc_lnum_t                  base_element_num,
                         const fvmc_coord_t           point_coords[],
                         _octree_t                  *octree,
                         fvmc_lnum_t                  points_in_extents[],
                         fvmc_lnum_t                  location[],
                         float                       distance[])
{
  fvmc_lnum_t  i, j, vertex_id, elt_num, triangle_vertices[6];
  int n_triangles;
  double elt_extents[6];

  fvmc_lnum_t n_points_in_extents = 0;

  /* If section contains polyhedra */

  if (this_section->type == FVMC_CELL_POLY)

    _polyhedra_section_locate(this_section,
                              parent_vertex_num,
                              vertex_coords,
                              tolerance,
                              base_element_num,
                              point_coords,
                              octree,
                              points_in_extents,
                              location,
                              distance);

  /* If section contains polygons */

  else if (this_section->type == FVMC_FACE_POLY)

    _polygons_section_locate_3d(this_section,
                                parent_vertex_num,
                                vertex_coords,
                                tolerance,
                                base_element_num,
                                point_coords,
                                octree,
                                points_in_extents,
                                location,
                                distance);

  /* If section contains regular elements */

  else {

    for (i = 0; i < this_section->n_elements; i++) {

      _Bool elt_initialized = false;

      if (base_element_num < 0) {
        if (this_section->parent_element_num != NULL)
          elt_num = this_section->parent_element_num[i];
        else
          elt_num = i + 1;
      }
      else
        elt_num = base_element_num + i;

      for (j = 0; j < this_section->stride; j++) {

        vertex_id = this_section->vertex_num[i*this_section->stride + j] - 1;

        _update_elt_extents(3,
                            vertex_id,
                            parent_vertex_num,
                            vertex_coords,
                            elt_extents,
                            &elt_initialized);

      }

      _elt_extents_finalize(3,
                            this_section->entity_dim,
                            tolerance,
                            elt_extents);

      _query_octree(elt_extents,
                    point_coords,
                    octree,
                    &n_points_in_extents,
                    points_in_extents);

      if (this_section->entity_dim == 3)

        _locate_in_cell_3d(elt_num,
                           this_section->type,
                           this_section->vertex_num + i*this_section->stride,
                           parent_vertex_num,
                           vertex_coords,
                           point_coords,
                           n_points_in_extents,
                           points_in_extents,
                           tolerance,
                           location,
                           distance);

      else if (this_section->entity_dim == 2) {

        if (this_section->type == FVMC_FACE_QUAD)

          n_triangles = fvmc_triangulate_quadrangle(3,
                                                   vertex_coords,
                                                   parent_vertex_num,
                                                   (  this_section->vertex_num
                                                    + i*this_section->stride),
                                                   triangle_vertices);

        else {

          assert(this_section->type == FVMC_FACE_TRIA);

          n_triangles = 1;
          for (j = 0; j < 3; j++)
            triangle_vertices[j]
              = this_section->vertex_num[i*this_section->stride + j];


        }

        _locate_on_triangles_3d(elt_num,
                                n_triangles,
                                triangle_vertices,
                                parent_vertex_num,
                                vertex_coords,
                                point_coords,
                                n_points_in_extents,
                                points_in_extents,
                                tolerance,
                                location,
                                distance);
      }

      else if (this_section->entity_dim == 1) {

        assert(this_section->type == FVMC_EDGE);

        _locate_on_edge_3d(elt_num,
                           this_section->vertex_num + i*this_section->stride,
                           parent_vertex_num,
                           vertex_coords,
                           point_coords,
                           n_points_in_extents,
                           points_in_extents,
                           tolerance,
                           location,
                           distance);

      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Find elements in a given section closest to 3d points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_section      <-- pointer to mesh section representation structure
 *   parent_vertex_num <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   base_element_num  <-- < 0 for location relative to parent element numbers,
 *                         number of elements in preceding sections of same
 *                         element dimension + 1 otherwise
 *   point_coords      <-- point coordinates
 *   n_point_ids       <-- number of points to locate
 *   point_id          <-- ids of points to locate
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, or absolute
 *                         distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_nodal_section_closest_3d(const fvmc_nodal_section_t  *this_section,
                          const fvmc_lnum_t           *parent_vertex_num,
                          const fvmc_coord_t           vertex_coords[],
                          fvmc_lnum_t                  base_element_num,
                          const fvmc_coord_t           point_coords[],
                          fvmc_lnum_t                  n_point_ids,
                          const fvmc_lnum_t            point_ids[],
                          fvmc_lnum_t                  location[],
                          float                       distance[])
{
  fvmc_lnum_t  i, j, elt_num, triangle_vertices[6];
  int n_triangles;

  /* If section contains polygons */

  if (this_section->type == FVMC_FACE_POLY)

    _polygons_section_closest_3d(this_section,
                                 parent_vertex_num,
                                 vertex_coords,
                                 base_element_num,
                                 point_coords,
                                 n_point_ids,
                                 point_ids,
                                 location,
                                 distance);

  /* If section contains regular elements */

  else {

    for (i = 0; i < this_section->n_elements; i++) {

      if (base_element_num < 0) {
        if (this_section->parent_element_num != NULL)
          elt_num = this_section->parent_element_num[i];
        else
          elt_num = i + 1;
      }
      else
        elt_num = base_element_num + i;

      if (this_section->entity_dim == 2) {

        if (this_section->type == FVMC_FACE_QUAD)

          n_triangles = fvmc_triangulate_quadrangle(3,
                                                   vertex_coords,
                                                   parent_vertex_num,
                                                   (  this_section->vertex_num
                                                    + i*this_section->stride),
                                                   triangle_vertices);

        else {

          assert(this_section->type == FVMC_FACE_TRIA);

          n_triangles = 1;
          for (j = 0; j < 3; j++)
            triangle_vertices[j]
              = this_section->vertex_num[i*this_section->stride + j];


        }

        _locate_on_triangles_3d(elt_num,
                                n_triangles,
                                triangle_vertices,
                                parent_vertex_num,
                                vertex_coords,
                                point_coords,
                                n_point_ids,
                                point_ids,
                                (HUGE_VAL / 4.),
                                location,
                                distance);
      }

      else if (this_section->entity_dim == 1) {

        assert(this_section->type == FVMC_EDGE);

        _locate_on_edge_3d(elt_num,
                           this_section->vertex_num + i*this_section->stride,
                           parent_vertex_num,
                           vertex_coords,
                           point_coords,
                           n_point_ids,
                           point_ids,
                           -1.,
                           location,
                           distance);
      }
    }

  }
}

/*----------------------------------------------------------------------------
 * Find elements in a given section containing 2d points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_section      <-- pointer to mesh section representation structure
 *   parent_vertex_num <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   tolerance         <-- associated tolerance
 *   base_element_num  <-- < 0 for location relative to parent element numbers,
 *                         number of elements in preceding sections of same
 *                         element dimension + 1 otherwise
 *   point_coords      <-- point coordinates
 *   quadtree          <-- point quadtree
 *   points_in_extents <-- array for query of ids of points in extents
 *                         (size: quadtree->n_points, less usually needed)
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, 0 - 1 if inside,
 *                         > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_nodal_section_locate_2d(const fvmc_nodal_section_t  *this_section,
                         const fvmc_lnum_t           *parent_vertex_num,
                         const fvmc_coord_t           vertex_coords[],
                         double                      tolerance,
                         fvmc_lnum_t                  base_element_num,
                         const fvmc_coord_t           point_coords[],
                         _quadtree_t                *quadtree,
                         fvmc_lnum_t                  points_in_extents[],
                         fvmc_lnum_t                  location[],
                         float                       distance[])
{
  fvmc_lnum_t  i, j, vertex_id, elt_num, _triangle_vertices[6];
  int n_vertices;
  double elt_extents[4];

  int n_triangles = 0;
  int n_vertices_max = 0;
  fvmc_lnum_t n_points_in_extents = 0;
  fvmc_lnum_t *triangle_vertices = _triangle_vertices;
  fvmc_triangulate_state_t *state = NULL;

  /* Return immediately if nothing to do for this rank */

  if (this_section->n_elements == 0)
    return;

  /* Count maximum number of vertices */

  if (this_section->type == FVMC_FACE_POLY) {

    for (i = 0; i < this_section->n_elements; i++) {

      n_vertices = (  this_section->vertex_index[i + 1]
                    - this_section->vertex_index[i]);

      if (n_vertices > n_vertices_max)
        n_vertices_max = n_vertices;

    }

    if (n_vertices_max < 3)
      return;

    BFTC_MALLOC(triangle_vertices, (n_vertices_max-2)*3, int);
    state = fvmc_triangulate_state_create(n_vertices_max);

  }

  else if (this_section->type == FVMC_FACE_QUAD)
    n_vertices_max = 4;

  else if (this_section->type == FVMC_FACE_TRIA)
    n_vertices_max = 3;

  /* Main loop on elements */

  for (i = 0; i < this_section->n_elements; i++) {

    _Bool elt_initialized = false;

    if (this_section->type == FVMC_FACE_POLY) {

      for (j = this_section->vertex_index[i];
           j < this_section->vertex_index[i + 1];
           j++) {
        vertex_id = this_section->vertex_num[j] - 1;

        _update_elt_extents(2,
                            vertex_id,
                            parent_vertex_num,
                            vertex_coords,
                            elt_extents,
                            &elt_initialized);

      }

    }
    else {

      for (j = 0; j < this_section->stride; j++) {

        vertex_id = this_section->vertex_num[i*this_section->stride + j] - 1;

        _update_elt_extents(2,
                            vertex_id,
                            parent_vertex_num,
                            vertex_coords,
                            elt_extents,
                            &elt_initialized);

      }

    }

    _elt_extents_finalize(2,
                          this_section->entity_dim,
                          tolerance,
                          elt_extents);

    if (base_element_num < 0) {
      if (this_section->parent_element_num != NULL)
        elt_num = this_section->parent_element_num[i];
      else
        elt_num = i + 1;
    }
    else
      elt_num = base_element_num + i;

    _query_quadtree(elt_extents,
                    point_coords,
                    quadtree,
                    &n_points_in_extents,
                    points_in_extents);

    /* Divide all face types into triangles */

    if (this_section->type == FVMC_FACE_POLY) {

      /* Triangulate polygon */

      n_vertices = (  this_section->vertex_index[i + 1]
                    - this_section->vertex_index[i]);
      vertex_id = this_section->vertex_index[i];

      n_triangles = fvmc_triangulate_polygon(2,
                                            n_vertices,
                                            vertex_coords,
                                            parent_vertex_num,
                                            (  this_section->vertex_num
                                             + vertex_id),
                                            FVMC_TRIANGULATE_MESH_DEF,
                                            triangle_vertices,
                                            state);

    }
    else if (this_section->type == FVMC_FACE_QUAD) {

      /* Triangulate quadrangle */

      n_triangles = fvmc_triangulate_quadrangle(2,
                                               vertex_coords,
                                               parent_vertex_num,
                                               (  this_section->vertex_num
                                                + i*this_section->stride),
                                               triangle_vertices);

    }

    else if (this_section->type == FVMC_FACE_TRIA) {

      /* Already a triangle */

      n_triangles = 1;

      for (j = 0; j < 3; j++)
        triangle_vertices[j]
          = this_section->vertex_num[i*this_section->stride + j];

    }

    /* Locate on triangulated face */

    if (this_section->entity_dim == 2)

      _locate_on_triangles_2d(elt_num,
                              n_triangles,
                              triangle_vertices,
                              parent_vertex_num,
                              vertex_coords,
                              point_coords,
                              n_points_in_extents,
                              points_in_extents,
                              tolerance,
                              location,
                              distance);

    else if (this_section->entity_dim == 1) {

      assert(this_section->type == FVMC_EDGE);

      _locate_on_edge_2d(elt_num,
                         this_section->vertex_num + i*this_section->stride,
                         parent_vertex_num,
                         vertex_coords,
                         point_coords,
                         n_points_in_extents,
                         points_in_extents,
                         tolerance,
                         location,
                         distance);

    }

  } /* End of loop on elements */

  /* Free axiliary arrays and structures */

  if (triangle_vertices != _triangle_vertices)
    BFTC_FREE(triangle_vertices);

  if (state != NULL)
    state = fvmc_triangulate_state_destroy(state);
}

/*----------------------------------------------------------------------------
 * Find elements in a given section closest to 3d points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_section      <-- pointer to mesh section representation structure
 *   parent_vertex_num <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   base_element_num  <-- < 0 for location relative to parent element numbers,
 *                         number of elements in preceding sections of same
 *                         element dimension + 1 otherwise
 *   point_coords      <-- point coordinates
 *   n_point_ids       <-- number of points to locate
 *   point_id          <-- ids of points to locate
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, or absolute
 *                         distance to a line element (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_nodal_section_closest_2d(const fvmc_nodal_section_t  *this_section,
                          const fvmc_lnum_t           *parent_vertex_num,
                          const fvmc_coord_t           vertex_coords[],
                          fvmc_lnum_t                  base_element_num,
                          const fvmc_coord_t           point_coords[],
                          fvmc_lnum_t                  n_point_ids,
                          const fvmc_lnum_t            point_id[],
                          fvmc_lnum_t                  location[],
                          float                       distance[])
{
  fvmc_lnum_t  i, elt_num;

  /* Return immediately if nothing to do for this rank */

  if (   this_section->n_elements == 0
      || this_section->entity_dim != 1)
    return;

  assert(this_section->type == FVMC_EDGE);

  /* Main loop on elements */

  for (i = 0; i < this_section->n_elements; i++) {

    if (base_element_num < 0) {
      if (this_section->parent_element_num != NULL)
        elt_num = this_section->parent_element_num[i];
      else
        elt_num = i + 1;
    }
    else
      elt_num = base_element_num + i;

    /* Locate on edge */

    _locate_on_edge_2d(elt_num,
                       this_section->vertex_num + i*this_section->stride,
                       parent_vertex_num,
                       vertex_coords,
                       point_coords,
                       n_point_ids,
                       point_id,
                       -1.0,
                       location,
                       distance);

  } /* End of loop on elements */
}

/*----------------------------------------------------------------------------
 * Find elements in a given section containing 1d points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this section, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_section      <-- pointer to mesh section representation structure
 *   parent_vertex_num <-- pointer to parent vertex numbers (or NULL)
 *   vertex_coords     <-- pointer to vertex coordinates
 *   tolerance         <-- associated tolerance
 *   base_element_num  <-- < 0 for location relative to parent element numbers,
 *                         number of elements in preceding sections of same
 *                         element dimension + 1 otherwise
 *   n_points          <-- number of points to locate
 *   point_coords      <-- point coordinates
 *   points_in_extents <-- array for query of ids of points in extents
 *                         (size: n_points, less usually needed)
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, 0 - 1 if inside,
 *                         > 1 if outside (size: n_points)
 *----------------------------------------------------------------------------*/

static void
_nodal_section_locate_1d(const fvmc_nodal_section_t  *this_section,
                         const fvmc_lnum_t           *parent_vertex_num,
                         const fvmc_coord_t           vertex_coords[],
                         double                      tolerance,
                         fvmc_lnum_t                  base_element_num,
                         fvmc_lnum_t                  n_points,
                         const fvmc_coord_t           point_coords[],
                         fvmc_lnum_t                  location[],
                         float                       distance[])
{
  fvmc_lnum_t  i, j, vertex_id, elt_num;
  fvmc_coord_t edge_coords[2];
  double delta, elt_extents[2];

  for (i = 0; i < this_section->n_elements; i++) {

    if (base_element_num < 0) {
      if (this_section->parent_element_num != NULL)
        elt_num = this_section->parent_element_num[i];
      else
        elt_num = i + 1;
    }
    else
      elt_num = base_element_num + i;

    for (j = 0; j < 2; j++) {

      vertex_id = this_section->vertex_num[i*this_section->stride] - 1;

      if (parent_vertex_num == NULL)
        edge_coords[j] = vertex_coords[vertex_id];
      else
        edge_coords[j] = vertex_coords[parent_vertex_num[vertex_id] - 1];

    }

    if (edge_coords[0] < edge_coords[1]) {
      elt_extents[0] = edge_coords[0];
      elt_extents[1] = edge_coords[1];
    }
    else {
      elt_extents[0] = edge_coords[1];
      elt_extents[1] = edge_coords[0];
    }

    delta = (elt_extents[1] - elt_extents[0]) * tolerance;

    elt_extents[0] -= delta;
    elt_extents[1] += delta;

    _locate_by_extents_1d(elt_num,
                          elt_extents,
                          n_points,
                          point_coords,
                          location,
                          distance);

  }

}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Find elements in a given nodal mesh containing points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are in an element of this mesh, or closer to one
 * than to previously encountered elements.
 *
 * parameters:
 *   this_nodal        <-- pointer to nodal mesh representation structure
 *   tolerance         <-- associated tolerance
 *   locate_on_parents <-- location relative to parent element numbers if
 *                         true, id of element + 1 in concatenated sections
 *                         of same element dimension if false
 *   n_points          <-- number of points to locate
 *   point_coords      <-- point coordinates
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated; 0 - 1 if inside,
 *                         and > 1 if outside a volume element, or absolute
 *                         distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
fvmc_point_location_nodal(const fvmc_nodal_t  *this_nodal,
                         double              tolerance,
                         _Bool               locate_on_parents,
                         fvmc_lnum_t          n_points,
                         const fvmc_coord_t   point_coords[],
                         fvmc_lnum_t          location[],
                         float               distance[])
{
  int i;
  int max_entity_dim;
  fvmc_lnum_t   base_element_num;
  fvmc_lnum_t  *points_in_extents = NULL;

  if (this_nodal == NULL)
    return;

  if (locate_on_parents == true)
    base_element_num = -1;
  else
    base_element_num = 1;

  max_entity_dim = fvmc_nodal_get_max_entity_dim(this_nodal);

  /* Build point query list
     (max size: n_points, usually much less) */

  BFTC_MALLOC(points_in_extents, n_points, fvmc_lnum_t);

  /* Use octree for 3d point location */

  if (this_nodal->dim == 3) {

    _octree_t  octree = _build_octree(n_points, point_coords);

    /* Locate for all sections */

    for (i = 0; i < this_nodal->n_sections; i++) {

      const fvmc_nodal_section_t  *this_section = this_nodal->sections[i];

      if (this_section->entity_dim == max_entity_dim) {

        _nodal_section_locate_3d(this_section,
                                 this_nodal->parent_vertex_num,
                                 this_nodal->vertex_coords,
                                 tolerance,
                                 base_element_num,
                                 point_coords,
                                 &octree,
                                 points_in_extents,
                                 location,
                                 distance);

        if (base_element_num > -1)
          base_element_num += this_section->n_elements;

      }

    }

    _free_octree(&octree);
  }

  /* Use quadtree for 2d point location */

  else if (this_nodal->dim == 2) {

    _quadtree_t  quadtree = _build_quadtree(n_points, point_coords);

    /* Locate for all sections */

    for (i = 0; i < this_nodal->n_sections; i++) {

      const fvmc_nodal_section_t  *this_section = this_nodal->sections[i];

      if (this_section->entity_dim == max_entity_dim) {

        _nodal_section_locate_2d(this_section,
                                 this_nodal->parent_vertex_num,
                                 this_nodal->vertex_coords,
                                 tolerance,
                                 base_element_num,
                                 point_coords,
                                 &quadtree,
                                 points_in_extents,
                                 location,
                                 distance);

        if (base_element_num > -1)
          base_element_num += this_section->n_elements;

      }

    }

    _free_quadtree(&quadtree);
  }

  /* Use brute force for 1d point location */

  else if (this_nodal->dim == 1) {

    for (i = 0; i < this_nodal->n_sections; i++) {

      const fvmc_nodal_section_t  *this_section = this_nodal->sections[i];

      if (this_section->entity_dim == 1) {

        _nodal_section_locate_1d(this_section,
                                 this_nodal->parent_vertex_num,
                                 this_nodal->vertex_coords,
                                 tolerance,
                                 base_element_num,
                                 n_points,
                                 point_coords,
                                 location,
                                 distance);

        if (base_element_num > -1)
          base_element_num += this_section->n_elements;

      }

    }

  }

  BFTC_FREE(points_in_extents);
}

/*----------------------------------------------------------------------------
 * Find elements in a given nodal mesh closest to points: updates the
 * location[] and distance[] arrays associated with a set of points
 * for points that are closer to an element of this mesh than to previously
 * encountered elements.
 *
 * This function currently only handles elements of lower dimension than
 * the spatial dimension.
 *
 * parameters:
 *   this_nodal        <-- pointer to nodal mesh representation structure
 *   locate_on_parents <-- location relative to parent element numbers if
 *                         true, id of element + 1 in concatenated sections
 *                         of same element dimension if false
 *   n_points          <-- number of points to locate
 *   point_coords      <-- point coordinates
 *   location          <-> number of element containing or closest to each
 *                         point (size: n_points)
 *   distance          <-> distance from point to element indicated by
 *                         location[]: < 0 if unlocated, or absolute
 *                         distance to a surface element (size: n_points)
 *----------------------------------------------------------------------------*/

void
fvmc_point_location_closest_nodal(const fvmc_nodal_t  *this_nodal,
                                 _Bool               locate_on_parents,
                                 fvmc_lnum_t          n_points,
                                 const fvmc_coord_t   point_coords[],
                                 fvmc_lnum_t          location[],
                                 float               distance[])
{
  int i;
  int  max_entity_dim;
  fvmc_lnum_t  base_element_num;
  fvmc_lnum_t  *point_ids = NULL;

  if (this_nodal == NULL)
    return;

  if (locate_on_parents == true)
    base_element_num = -1;
  else
    base_element_num = 1;

  max_entity_dim = fvmc_nodal_get_max_entity_dim(this_nodal);

  if (max_entity_dim == this_nodal->dim)
    bftc_error(__FILE__, __LINE__, 0,
              _("Locating volume elements closest to points not handled yet"));

  if (this_nodal->dim > 1) {
    fvmc_lnum_t j;
    BFTC_MALLOC(point_ids, n_points, fvmc_lnum_t);
    for (j = 0; j < n_points; j++)
      point_ids[j] = j;
  }

  /* Use brute force for closest 3d point location */

  if (this_nodal->dim == 3) {

    /* Locate for all sections */

    for (i = 0; i < this_nodal->n_sections; i++) {

      const fvmc_nodal_section_t  *this_section = this_nodal->sections[i];

      if (this_section->entity_dim == max_entity_dim) {

        _nodal_section_closest_3d(this_section,
                                  this_nodal->parent_vertex_num,
                                  this_nodal->vertex_coords,
                                  base_element_num,
                                  point_coords,
                                  n_points,
                                  point_ids,
                                  location,
                                  distance);

        if (base_element_num > -1)
          base_element_num += this_section->n_elements;

      }

    }
  }

  /* Use brute force for closest 2d point location */

  else if (this_nodal->dim == 2) {

    /* Locate for all sections */

    for (i = 0; i < this_nodal->n_sections; i++) {

      const fvmc_nodal_section_t  *this_section = this_nodal->sections[i];

      if (this_section->entity_dim == max_entity_dim) {

        _nodal_section_closest_2d(this_section,
                                  this_nodal->parent_vertex_num,
                                  this_nodal->vertex_coords,
                                  base_element_num,
                                  point_coords,
                                  n_points,
                                  point_ids,
                                  location,
                                  distance);

        if (base_element_num > -1)
          base_element_num += this_section->n_elements;

      }

    }

  }

  if (point_ids != NULL)
    BFTC_FREE(point_ids);
}

/*----------------------------------------------------------------------------*/

#undef _DOT_PRODUCT
#undef _MODULE
#undef _CROSS_PRODUCT

#ifdef __cplusplus
}
#endif /* __cplusplus */
