/*
 * File:   pdm_octree.h
 * Author: equemera
 *
 * Created on November 8, 2017, 11:27 AM
 */

#ifndef PDM_OCTREE_SEQ_H
#define	PDM_OCTREE_SEQ_H

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#if 0
} /* Fake brace */
#endif
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_octree_seq_t PDM_octree_seq_t;

/**
 * \enum PDM_octree_seq_child_t
 * \brief Names of 8 children of a node
 *
 */

typedef enum {
  PDM_OCTREE_SEQ_NADIR,
  PDM_OCTREE_SEQ_ZENITH,
  PDM_OCTREE_SEQ_WEST,
  PDM_OCTREE_SEQ_EAST,
  PDM_OCTREE_SEQ_NORTH,
  PDM_OCTREE_SEQ_SOUTH,
} PDM_octree_seq_direction_t;

/**
 * \enum PDM_octree_seq_child_t
 * \brief Names of 8 children of a node
 *
 */

typedef enum {
  PDM_OCTREE_SEQ_NORTH_WEST_NADIR,
  PDM_OCTREE_SEQ_NORTH_WEST_ZENITH,
  PDM_OCTREE_SEQ_NORTH_EAST_NADIR,
  PDM_OCTREE_SEQ_NORTH_EAST_ZENITH,
  PDM_OCTREE_SEQ_SOUTH_WEST_NADIR,
  PDM_OCTREE_SEQ_SOUTH_WEST_ZENITH,
  PDM_OCTREE_SEQ_SOUTH_EAST_NADIR,
  PDM_OCTREE_SEQ_SOUTH_EAST_ZENITH,
} PDM_octree_seq_child_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create an octree structure
 *
 * \param [in]   n_point_cloud      Number of point cloud
 * \param [in]   depth_max          Maximum depth
 * \param [in]   points_in_leaf_max Maximum points in a leaf
 * \param [in]   tolerance          Relative geometric tolerance
 *
 * \return     Pointer to \ref PDM_octree_seq object
 */

PDM_octree_seq_t *
PDM_octree_seq_create
(
 const int    n_point_cloud,
 const int    depth_max,
 const int    points_in_leaf_max,
 const double tolerance
);


/**
 *
 * \brief Free an octree structure
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 *
 */

void
PDM_octree_seq_free
(
 PDM_octree_seq_t *octree
);


/**
 *
 * \brief Set a point cloud
 *
 * \param [in]   id                 Identifier
 * \param [in]   i_point_cloud      Number of point cloud
 * \param [in]   n_points           Maximum depth
 * \param [in]   coords             Point coordinates
 *
 */

void
PDM_octree_seq_point_cloud_set
(
 PDM_octree_seq_t *octree,
 const int         i_point_cloud,
 const int         n_points,
 const double     *coords
);


/**
 *
 * \brief Build octree
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 *
 */

void
PDM_octree_seq_build
(
 PDM_octree_seq_t *octree
);


/**
 *
 * \brief Get root node id
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 *
 * \return     Root node identifier (-1 if octree is not built)
 *
 */

int
PDM_octree_seq_root_node_id_get
(
 PDM_octree_seq_t *octree
);


/**
 *
 * \brief Get extents
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 *
 * \return     Extents
 *
 */

double *
PDM_octree_seq_extents_get
(
 PDM_octree_seq_t *octree
);


/**
 *
 * \brief Get ancestor node id
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 * \param [in]   node_id            Node identifier
 *
 * \return     Ancestor node identifier
 *
 */

int
PDM_octree_seq_ancestor_node_id_get
(
 PDM_octree_seq_t *octree,
 const int         node_id
);


/**
 *
 * \brief Get node extents
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 * \param [in]   node_id            Node identifier
 *
 * \return     Extents
 *
 */

const double *
PDM_octree_seq_node_extents_get
(
 PDM_octree_seq_t *octree,
 const int         node_id
);


/**
 *
 * \brief Get children of a node
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 * \param [in]   node_id            Node identifier
 * \param [in]   child              Children
 *
 * \return     Children node id
 *
 */

int
PDM_octree_seq_children_get
(
 PDM_octree_seq_t             *octree,
 const int                     node_id,
 const PDM_octree_seq_child_t  child
);


/**
 *
 * \brief Get Neighbor of node
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 * \param [in]   node_id            Node identifier
 * \param [in]   direction          Neighbor direction
 *
 * \return     Neighbor node id (-1 if no neighbor)
 *
 */

int
PDM_octree_seq_neighbor_get
(
 PDM_octree_seq_t                 *octree,
 const int                         node_id,
 const PDM_octree_seq_direction_t  direction
);

/**
 *
 * \brief Get the number of point inside a node
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 * \param [in]   node_id            Node identifier
 *
 * \return   Number of points
 *
 */

int
PDM_octree_seq_n_points_get
(
 PDM_octree_seq_t        *octree,
 const int                node_id
);


/**
 *
 * \brief Get indexes of points inside a node
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 * \param [in]   node_id            Node identifier
 * \param [out]  point_clouds_id    Point clouds number
 *                                  (size = Number of points inside the node)
 * \param [out]  point_indexes      Point indexes
 *                                  (size = Number of points inside the node)
 *
 */

void
PDM_octree_seq_points_get
(
 PDM_octree_seq_t        *octree,
 const int                node_id,
 int                    **point_clouds_id,
 int                    **point_indexes
);


/**
 *
 * \brief Is it a leaf
 *
 * \param [in]   octree             Pointer to \ref PDM_octree_seq object
 * \param [in]   node_id            Node identifier
 *
 * \return   1 or 0
 *
 */

int
PDM_octree_seq_leaf_is
(
 PDM_octree_seq_t        *octree,
 const int                node_id
);


/**
 *
 * \brief Look for closest points stored inside an octree
 *
 * \param [in]   octree                 Pointer to \ref PDM_octree_seq object
 * \param [in]   n_pts                  Number of points
 * \param [in]   pts                    Point Coordinates
 * \param [out]  closest_octree_pt_id   Closest point in octree index
 * \param [out]  closest_octree_pt_dist Closest point in octree distance
 *
 */

void
PDM_octree_seq_closest_point
(
PDM_octree_seq_t *octree,
const int         n_pts,
double           *pts,
int              *closest_octree_pt_id,
double           *closest_octree_pt_dist2
);


/**
 *
 * \brief Write octants in a VTK file
 *
 * \param [in]   octree                 Pointer to \ref PDM_octree_seq object
 * \param [in]   filename               Output file name
 *
 */

void PDM_octree_seq_write_octants
(
 PDM_octree_seq_t *octree,
 const char       *filename
 );

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_OCTREE_SEQ_H */

