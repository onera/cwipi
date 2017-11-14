/* 
 * File:   pdm_octree.h
 * Author: equemera
 *
 * Created on November 8, 2017, 11:27 AM
 */

#ifndef PDM_OCTREE_H
#define	PDM_OCTREE_H

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
 * \enum PDM_octree_child_t
 * \brief Names of 8 children of a node 
 *
 */

typedef enum {
  PDM_NORTH_WEST_NADIR,   
  PDM_NORTH_WEST_ZENITH,   
  PDM_NORTH_EAST_NADIR,   
  PDM_NORTH_EAST_ZENITH,   
  PDM_SOUTH_WEST_NADIR,   
  PDM_SOUTH_WEST_ZENITH,   
  PDM_SOUTH_EAST_NADIR,   
  PDM_SOUTH_EAST_ZENITH,   
} PDM_octree_child_t;

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
PDM_octree_create
(
 const int n_point_cloud,
 const int depth_max, 
 const int points_in_leaf_max,
 const double tolerance, 
 const PDM_MPI_Comm comm
);

//void
//PROCF (pdm_octree_create, PDM_OCTREE_CREATE)
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
PDM_octree_free
(
 const int          id
);

//void
//PROCF (pdm_octree_free, PDM_OCTREE_FREE)
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
PDM_octree_point_cloud_set
(
 const int          id,
 const int          i_point_cloud,
 const int          n_points,
 const double      *coords 
);

//void
//PROCF (pdm_octree_point_cloud_set, PDM_OCTREE_POINT_CLOUD_SET)
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
PDM_octree_build
(
 const int          id
);


/**
 *
 * \brief Processes extents  
 *
 * \param [in]   id                 Identifier 
 * \param [in]   i_proc             Process
 *
 */

const double *
PDM_octree_processes_extents_get
(
 const int          id,
 const int          i_proc
);

//void
//PROCF (pdm_octree_build, PDM_OCTREE_BUILD)
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
PDM_octree_root_node_id_get
(
 const int          id
);

//void
//PROCF (pdm_octree_root_node_id_get, PDM_OCTREE_ROOT_NODE_ID_GET)
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
PDM_octree_ancestor_node_id_get
(
 const int          id,
 const int          node_id
);

//void
//PROCF (pdm_octree_ancestor_node_id_get, PDM_OCTREE_ANCESTOR_NODE_ID_GET)
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
PDM_octree_node_extents_get
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
PDM_octree_children_get
(
 const int                id,
 const int                node_id,
 const PDM_octree_child_t child
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
PDM_octree_neighbor_get
(
 const int                id,
 const int                node_id,
 const PDM_octree_child_t child
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
PDM_octree_n_points_get
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

int
PDM_octree_points_get
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
PDM_octree_leaf_is
(
 const int                id,
 const int                node_id
);


#ifdef	__cplusplus
}
#endif

#endif	/* PDM_OCTREE_H */

