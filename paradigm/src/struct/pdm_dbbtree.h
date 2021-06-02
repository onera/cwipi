#ifndef __PDM_DBBTREE_H__
#define __PDM_DBBTREE_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm.h"
#include "pdm_box.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/**
 * \struct PDM_dbbtree_t
 * \brief  Distributed boundary box tree
 *
 *  PDM_dbbtree_t defines a distributed boundary box tree
 *
 */

typedef struct _PDM_dbbtree_t PDM_dbbtree_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 * \brief Return an intialized \ref PDM_dbbtree_t structure
 *
 * This function returns an initialized \ref PDM_dbbtree_t structure
 *
 * \param [in]  comm             Associated communicator
 * \param [in]  dim              boxes dimension
 * \param [in]  global_extents   Globals of elements to storage into the tree
 *                               (automatic computation if NULL)
 *
 * \return      A new initialized \ref PDM_dbbtree_t structure
 *
 */

PDM_dbbtree_t *
PDM_dbbtree_create
(
 PDM_MPI_Comm  comm,
 int           dim,
 double       *global_extents
 );



/**
 * \brief Free a \ref PDM_dbbtree_t structure
 *
 * \param [in]  dbbt   Pointer to a distributed bounding box tree
 *
 * \return      NULL
 *
 */

PDM_dbbtree_t *
PDM_dbbtree_free
(
 PDM_dbbtree_t     *dbbt
 );


/**
 * \brief Assign a set of boxes to an empty \ref PDM_dbbtree_t structure.
 *
 * This function assigns a set of boxes to an empty \ref PDM_dbbtree_t structure.
 *
 * \param [in]  dbbt     Pointer to a distributed bounding box tree
 * \param [in]  n_part    Number of partitions
 * \param [in]  nElts    Number of elements of each partition
 * \param [in]  extents  Extents of each element of each partition
 * \param [in]  gNum     Global number of each element of each partition
 *
 * \return associated \ref PDM_box_set_t structure distributed according to
 * the tree location
 *
 */

PDM_box_set_t  *
PDM_dbbtree_boxes_set
(
 PDM_dbbtree_t     *dbbt,
 const int          n_part,
 const int         *nElts,
 const double     **extents,
 const PDM_g_num_t **gNum
 );


/**
 * \brief Assign boxes to intersect to the tree.
 *
 * This function  assigns boxes to intersect to the tree.
 *
 * \param [in]  dbbt     Pointer to a distributed bounding box tree
 * \param [in]  n_part    Number of partitions
 * \param [in]  nElts    Number of elements of each partition
 * \param [in]  extents  Extents of each element of each partition
 * \param [in]  gNum     Global number of each element of each partition
 * \param [out] box_index Pointer to the index array on associated tree bounding boxeq
 * \param [out] box_g_num Pointer to the list of intersecting bounding boxes
 *
 * \return associated \ref PDM_box_set_t structure distributed according
 * to the tree intersection
 *
 */


PDM_box_set_t  *
PDM_dbbtree_intersect_boxes_set
(
 PDM_dbbtree_t    *dbbt,
 const int         n_part,
 const int        *nElts,
 const double     **extents,
 const PDM_g_num_t **gNum,
 int              *box_index[],
 int              *box_l_num[]
 );


/**
 *
 * Get the boxes closer than the upper bound distance
 *
 *   \param [in] bt               Pointer to box tree structure
 *   \param [in] n_pts            Number of points
 *   \param [in] pts              Point coordinates (size = 3 * n_pts)
 *   \param [in] pts_g_num        Point global numbers
 *   \param [in] upper_bound_dist2 Upper bound of the squer of the distance (size = n_pts)
 *   \param [out] box_index       Index of boxes (size = n_pts + 1)
 *   \param [out] box_g_num       Global num of boxes (size = i_boxes[n_pts])
 *
 */

void
PDM_dbbtree_closest_upper_bound_dist_boxes_get
(
 PDM_dbbtree_t    *dbbt,
 const int        n_pts,
 double           pts[],
 PDM_g_num_t      pts_g_num[],
 double           upper_bound_dist2[],
 int             *box_index[],
 PDM_g_num_t     *box_g_num[]
 );


/**
 *
 * Get the location of a point cloud
 *
 *   \param [in] bt                    Pointer to box tree structure
 *   \param [in] n_pts                 Number of points
 *   \param [in] pts_coord             Point coordinates (size = 3 * n_pts)
 *   \param [in] pts_g_num             Point global numbers (size = n_pts)
 *   \param [in] n_boxes               Number of boxes
 *   \param [in] box_g_num             Global num of boxes (size = n_boxes)
 *   \param [inout] pts_in_box_idx     Index of points in boxes (size = n_boxes+1, allocated inside function)
 *   \param [inout] pts_in_box_g_num   Global num of points in boxes (size = pts_in_box_idx[n_boxes], allocated inside function)
 *   \param [inout] pts_in_box_coord   Coordinates of points in boxes (size = 3*pts_in_box_idx[n_boxes], allocated inside function)
 *
 */

void PDM_dbbtree_points_inside_boxes
(
 PDM_dbbtree_t      *dbbt,
 const int           n_pts,
 PDM_g_num_t         pts_g_num[],
 double              pts_coord[],
 const int           n_boxes,
 const PDM_g_num_t   box_g_num[],
 int               **pts_in_box_idx,
 PDM_g_num_t       **pts_in_box_g_num,
 double            **pts_in_box_coord
 );


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_DBBTREE_H__ */
