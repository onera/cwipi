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
 * \param [in]   comm               MPI communicator
 *
 * \return     Identifier    
 */

int
PDM_octree_seq_create
(
 const int n_point_cloud,
 const int depth_max, 
 const int points_in_leaf_max,
 const double tolerance
);

//void
//PROCF (pdm_octree_seq_create, PDM_OCTREE_SEQ_CREATE)
//(
// const int *n_point_cloud,
// const int *depth_max, 
// const int *points_in_leaf_max,
// const double *tolerance, 
// const PDM_MPI_Fint *fcomm,
// const int *id
//);

/**
 *
 * \brief Free an octree structure   
 *
 * \param [in]   id                 Identifier 
 *  
 */

void
PDM_octree_seq_free
(
 const int          id
);

//void
//PROCF (pdm_octree_seq_free, PDM_OCTREE_SEQ_FREE)
//(
// const int          *id
//);


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
 const int          id,
 const int          i_point_cloud,
 const int          n_points,
 const double      *coords 
);

//void
//PROCF (pdm_octree_seq_point_cloud_set, PDM_OCTREE_SEQ_POINT_CLOUD_SET)
//(
// const int          *id
// const int          *i_point_cloud,
// const int          *n_points,
// const double       *coords 
//);


/**
 *
 * \brief Build octree  
 *
 * \param [in]   id                 Identifier 
 *
 */

void
PDM_octree_seq_build
(
 const int          id
);

//void
//PROCF (pdm_octree_seq_build, PDM_OCTREE_SEQ_BUILD)
//(
// const int          *id
//);

/**
 *
 * \brief Get root node id  
 *
 * \param [in]   id                 Identifier 
 *
 * \return     Root node identifier    
 * 
 */

int
PDM_octree_seq_root_node_id_get
(
 const int          id
);


/**
 *
 * \brief Get extents  
 *
 * \param [in]   id                 Identifier 
 *
 * \return     Extents    
 * 
 */

double *
PDM_octree_seq_extents_get
(
 const int          id
);

//void
//PROCF (pdm_octree_seq_root_node_id_get, PDM_OCTREE_SEQ_ROOT_NODE_ID_GET)
//(
// const int          *id,
// int                *root_node_id
//);


/**
 *
 * \brief Get ancestor node id  
 *
 * \param [in]   id                 Identifier 
 * \param [in]   node_id            Node identifier 
 *
 * \return     Ancestor node identifier    
 * 
 */

int
PDM_octree_seq_ancestor_node_id_get
(
 const int          id,
 const int          node_id
);

//void
//PROCF (pdm_octree_seq_ancestor_node_id_get, PDM_OCTREE_SEQ_ANCESTOR_NODE_ID_GET)
//(
// const int          *id,
// const int          *node_id, 
// int                *ancestor_node_id
//);


/**
 *
 * \brief Get node extents  
 *
 * \param [in]   id                 Identifier 
 * \param [in]   node_id            Node identifier 
 *
 * \return     Extents    
 * 
 */

const double *
PDM_octree_seq_node_extents_get
(
 const int          id,
 const int          node_id
);


/**
 *
 * \brief Get children of a node 
 *
 * \param [in]   id                 Identifier 
 * \param [in]   node_id            Node identifier 
 * \param [in]   child              Children 
 *
 * \return     Children node id    
 * 
 */

int
PDM_octree_seq_children_get
(
 const int                id,
 const int                node_id,
 const PDM_octree_seq_child_t child
);


/**
 *
 * \brief Get Neighbor of node 
 *
 * \param [in]   id                 Identifier 
 * \param [in]   node_id            Node identifier 
 * \param [in]   direction          Neighbor direction 
 *
 * \return     Neighbor node id    
 * 
 */

int
PDM_octree_seq_neighbor_get
(
 const int                    id,
 const int                    node_id,
 const PDM_octree_seq_direction_t direction
);

/**
 *
 * \brief Get the number of point inside a node 
 *
 * \param [in]   id                 Identifier 
 * \param [in]   node_id            Node identifier 
 *
 * \return   Number of points    
 * 
 */

int
PDM_octree_seq_n_points_get
(
 const int                id,
 const int                node_id
);


/**
 *
 * \brief Get indexes of points inside a node 
 *
 * \param [in]   id                 Identifier 
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
 const int                id,
 const int                node_id,
 int                    **point_clouds_id, 
 int                    **point_indexes 
);


/**
 *
 * \brief Is it a leaf 
 *
 * \param [in]   id                 Identifier 
 * \param [in]   node_id            Node identifier 
 *
 * \return   1 or 0    
 * 
 */

int
PDM_octree_seq_leaf_is
(
 const int                id,
 const int                node_id
);


/**
 *
 * Look for closest points stored inside an octree
 *
 * parameters:
 * \param [in]   id                     Identifier
 * \param [in]   n_pts                  Number of points
 * \param [in]   pts                    Point Coordinates
 * \param [out]  closest_octree_pt_id   Closest point in octree index (couple icloud, index)
 * \param [out]  closest_octree_pt_dist Closest point in octree distance
 *  
 */

void
PDM_octree_seq_closest_point
(
const int    id,
const int    n_pts,
double      *pts,
int         *closest_octree_pt_id,
double      *closest_octree_pt_dist2
);

#ifdef	__cplusplus
}
#endif

#endif	/* PDM_OCTREE_SEQ_H */

