#ifndef __PDM_BOX_TREE_H__
#define __PDM_BOX_TREE_H__

/*============================================================================
 * Search octrees and quadtrees of boxes.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_box.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _PDM_box_tree_data_t PDM_box_tree_data_t;
typedef struct _PDM_box_tree_t PDM_box_tree_t;

typedef enum {

  PDM_BOX_TREE_ASYNC_LEVEL,  /* Boxes are placed according to tree parameters,
                                and potentially at different levels */
  PDM_BOX_TREE_SYNC_LEVEL    /* All boxes are placed for all ranks at the
                                same level */

} PDM_box_tree_sync_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a PDM_box_tree_data_t structure and initialize it.
 *
 * returns:
 *   pointer to an empty PDM_box_tree_data_t structure.
 *----------------------------------------------------------------------------*/
PDM_box_tree_data_t *
PDM_box_tree_data_create(void);

/*----------------------------------------------------------------------------
 * Create a PDM_box_tree_t structure and initialize it.
 *
 * parameters:
 *  max_level     <-- max possible level
 *  threshold     <-- max number of  boxes linked to an octant if
 *                    max_level is not reached
 *  max_box_ratio <-- max n_linked_boxes / n_boxes ratio
 *
 * returns:
 *   pointer to an empty PDM_box_tree_t structure.
 *----------------------------------------------------------------------------*/

PDM_box_tree_t *
PDM_box_tree_create(int    max_level,
                    int    threshold,
                    float  max_box_ratio);

/*----------------------------------------------------------------------------
 * Destroy a PDM_box_tree_data_t structure.
 *
 * parameters:
 *   btd <-- pointer to pointer to PDM_box_tree_data_t structure to destroy
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_data_destroy(PDM_box_tree_data_t  *btd);

/*----------------------------------------------------------------------------
 * Destroy a PDM_box_tree_t structure.
 *
 * parameters:
 *   bt <-- pointer to pointer to PDM_box_tree_t structure to destroy
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_destroy(PDM_box_tree_t  **bt);

/*----------------------------------------------------------------------------
 * Get the deepest level allowed by the tree structure.
 *
 * parameters:
 *   bt <-- pointer to PDM_box_tree_t structure.
 *
 * returns:
 *   deepest allowed level of the tree
 *----------------------------------------------------------------------------*/

int
PDM_box_tree_get_max_level(const PDM_box_tree_t  *bt);

/*----------------------------------------------------------------------------
 * Assign a set of boxes to an empty PDM_box_tree_t structure.
 *
 * The box tree structure must have been created using to PDM_tree_create().
 *
 * The depth of the tree is adjusted so that a maximum of max_n_elts boxes
 * will be assigned to each leaf, unless this would require going beyond
 * the tree's maximum level.
 *
 * If max_level = -1, the highest level reachable is PDM_TREE_MAX_LEVEL but
 * there is no defined target level.
 *
 * parameters:
 *   bt         <-> pointer to PDM_box_tree_t structure.
 *   boxes      <-- pointer to the associated box set structure
 *   build_type <-- layout variant for building the tree structure
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_set_boxes(PDM_box_tree_t       *bt,
                       PDM_box_set_t  *boxes,
                       PDM_box_tree_sync_t   build_type);

void
PDM_box_tree_set_copied_rank_boxes(PDM_box_tree_t       *bt);


/*----------------------------------------------------------------------------
 * Compute an index based on Morton encoding to ensure a good distribution
 * of boxes among the participating ranks.
 *
 * parameters:
 *   bt         <-> pointer to PDM_box_tree_t structure.
 *   boxes      <-- pointer to the associated box set structure
 *
 * returns:
 *   pointer to newly created PDM_box_distrib_t structure.
 *----------------------------------------------------------------------------*/

PDM_box_distrib_t *
PDM_box_tree_get_distrib(PDM_box_tree_t        *bt,
                         const PDM_box_set_t   *boxes);

/*----------------------------------------------------------------------------
 * Build an indexed list on associated bt boxes to list boxes B intersections.
 *
 * The index and box_g_num arrays are allocated by this function,
 * and it is the caller's responsibility to free them.
 *
 * Upon return, box_index[i] points to the first position in box_g_num
 * relative to boxes intersecting box i of the boxes set, while
 * box_g_num contains the global numbers associated with those boxes.
 *
 * parameters:
 *   bt        <-- pointer to box tree structure to query
 *   boxesB    <-- pointer boxes that intersect associated tree boxes
 *   box_index --> pointer to the index array on bounding boxes
 *   box_g_num --> pointer to the list of intersecting bounding boxes
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_get_boxes_intersects(PDM_box_tree_t       *bt,
                                  const PDM_box_set_t  *boxesB,
                                  int                  *box_index[],
                                  int                  *box_l_num[]);

/*----------------------------------------------------------------------------
 * Build an indexed list on boxes to list intersections.
 *
 * The index and box_g_num arrays are allocated by this function,
 * and it is the caller's responsibility to free them.
 *
 * Upon return, box_index[i] points to the first position in box_g_num
 * relative to boxes intersecting box i of the boxes set, while
 * box_g_num contains the global numbers associated with those boxes.
 *
 * parameters:
 *   bt        <-- pointer to box tree structure to query
 *   box_index --> pointer to the index array on bounding boxes
 *   box_g_num --> pointer to the list of intersecting bounding boxes
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_get_intern_intersects(PDM_box_tree_t       *bt,
                                   int                  *box_index[],
                                   PDM_g_num_t           *box_g_num[]);

/*----------------------------------------------------------------------------
 * Get global box tree statistics.
 *
 * All fields returned are optional: if their argument is set to NULL,
 * the corresponding information will not be returned.
 *
 * For each field not set to NULL, 3 values are always returned:
 * the mean on all ranks (rounded to the closest integer), the minimum,
 * and the maximum value respectively.
 *
 * In serial mode, the mean, minimum, and maximum will be identical for most
 * fields, but all 3 values are returned nonetheless.
 *
 * Note that the theoretical memory use includes that of the associated
 * box set.
 *
 * parameters:
 *   bt                 <-- pointer to box tree structure
 *   depth              --> tree depth (max level used)
 *   n_leaves           --> number of leaves in the tree
 *   n_boxes            --> number of boxes in the tree
 *   n_threshold_leaves --> number of leaves where n_boxes > threshold
 *   n_leaf_boxes       --> number of boxes for a leaf
 *   mem_used           --> theoretical used memory
 *   mem_allocated      --> theoretical allocated memory
 *
 * returns:
 *   the spatial dimension associated with the box tree layout (3, 2, or 1)
 *----------------------------------------------------------------------------*/

int
PDM_box_tree_get_stats(const PDM_box_tree_t  *bt,
                       int                    depth[3],
                       int              n_leaves[3],
                       int              n_boxes[3],
                       int              n_threshold_leaves[3],
                       int              n_leaf_boxes[3],
                       size_t                 mem_used[3],
                       size_t                 mem_allocated[3]);

/*----------------------------------------------------------------------------
 * Display local statistics about a PDM_box_tree_t structure.
 *
 * parameters:
 *   bt <-- pointer to box tree structure
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_dump_statistics(const PDM_box_tree_t  *bt);

/*----------------------------------------------------------------------------
 * Dump an PDM_box_tree_t structure.
 *
 * parameters:
 *   bt <-- pointer to box tree structure
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_dump(PDM_box_tree_t  *bt);

/*----------------------------------------------------------------------------
 * Get minimum of maximum distance of boxes
 *
 * parameters:
 *   bt              <-- pointer to box tree structure
 *   n_pts           <-- Number of points
 *   pts             <-- Point coordinates
 *   box_id          --> Leaf box with the minimum of maximum distance
 *   box_max_dist[]  --> Maximum distance to box_id
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_min_dist_max_box
(
 PDM_box_tree_t  *bt,
 const int        n_pts,
 double          *pts,
 int             *box_id,
 double          *box_max_dist
 );


/*----------------------------------------------------------------------------
 * Get minimum of maximum distance of boxes
 *
 * parameters:
 *   bt                <-- pointer to box tree structure
 *   n_pts             <-- Number of points
 *   pts               <-- Point coordinates (size = 3 * n_pts)
 *   upper_bound_dist2 <-- Upper bound of the square of the distance (size = n_pts)
 *   i_boxes           --> Index of boxes (size = n_pts + 1)
 *   boxes             --> Boxes (size = i_boxes[n_pts])
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_closest_upper_bound_dist_boxes_get
(
 PDM_box_tree_t  *bt,
 const int        n_pts,
 double           pts[],
 double           upper_bound_dist2[],
 int             *i_boxes[],
 int             *boxes[]
 );

/*----------------------------------------------------------------------------
 * same as above but extended to support search in local box tree as well as in box trees from copied ranks
 *
 * if i_rank < 0 then search is performed in local box tree data,
 * otherwise search is performed in box tree data copied from proc with rank bt->copied_rank[i_rank]
 *
 *----------------------------------------------------------------------------*/

void
PDM_box_tree_closest_upper_bound_dist_boxes_get_v2
(
 PDM_box_tree_t  *bt,
 const int        i_rank,
 const int        n_pts,
 double           pts[],
 double           upper_bound_dist2[],
 int             *i_boxes[],
 int             *boxes[],
 const double    *d_opt
 );


/*----------------------------------------------------------------------------
 * Get minimum of maximum distance of boxes (extended version: rank tree data)
 *
 * parameters:
 *   bt                <-- pointer to box tree structure
 *   i_pts             <-- index of points (size = bt->n_copied_ranks+1)
 *   pts               <-- Point coordinates (size = 3 * n_pts) (with n_pts = i_pts[bt->n_copied_ranks])
 *   upper_bound_dist2 <-- Upper bound of the square of the distance (size = n_pts)
 *   i_boxes_rank      --> Index of boxes (size = bt->n_copied_ranks)
 *                            i_boxes_rank[r] is of size n_pts_rank[r] + 1 (with n_pts_rank[r] = i_pts[r+1] - i_pts[r])
 *   boxes_rank        --> Boxes (size = bt->n_copied_ranks)
 *                            boxes_rank[r] is of size i_boxes_rank[r][n_pts_rank[r]]
 *----------------------------------------------------------------------------*/
void
PDM_box_tree_closest_upper_bound_dist_boxes_get_from_copied_ranks
(
 PDM_box_tree_t  *bt,
 const int        i_pts[],
 double           pts[],
 double           upper_bound_dist2[],
 int            **i_boxes_rank[],
 int            **boxes_rank[]
 );




/*----------------------------------------------------------------------------
 * Send copies of box tree data from selected ranks to all other ranks for better load balancing
 *---------------------------------------------------------------------------*/
void
PDM_box_tree_copy_to_ranks
(
 PDM_box_tree_t *bt,
 int            *n_copied_ranks,
 int            *copied_ranks,
 int            *rank_copy_num
 );



void
PDM_box_tree_points_inside_boxes
(
 PDM_box_tree_t     *bt,
 const int           n_pts,
 const PDM_g_num_t   pts_g_num[],
 const double        pts_coord[],
 int               **pts_in_box_idx,
 PDM_g_num_t       **pts_in_box_g_num,
 double            **pts_in_box_coord
 );


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_BOX_TREE_H__ */
