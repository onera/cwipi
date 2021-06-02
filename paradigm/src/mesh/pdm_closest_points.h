#ifndef PDM_CLOSEST_POINTS_H
#define PDM_CLOSEST_POINTS_H

/*
  This file is part of the ParaDiGM library.

  Copyright (C) 2019       ONERA

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
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct _pdm_closest_point_t PDM_closest_point_t;

/*============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Create a structure to look for the closest points of a point cloud
 * (target cloud) in an other point cloud (source cloud)
 *
 * \param [in]   comm           MPI communicator
 * \param [in]   n_closest      Number of closest source points to find for each
 *                              target point
 *
 * \return     Identifier
 *
 */

PDM_closest_point_t*
PDM_closest_points_create
(
 const PDM_MPI_Comm    comm,
 const int             n_closest,
 const PDM_ownership_t owner
);

PDM_closest_point_t*
PDM_closest_points_create_cf
(
 const PDM_MPI_Fint     comm,
 const int              n_closest,
 const PDM_ownership_t  owner
);


/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   id                Identifier
 * \param [in]   n_part_cloud_src  Number of partitions of the source cloud
 * \param [in]   n_part_cloud_tgt  Number of partitions od the target cloud
 *
 */

void
PDM_closest_points_n_part_cloud_set
(
       PDM_closest_point_t* cls,
 const int                  n_part_cloud_src,
 const int                  n_part_cloud_tgt
);


/**
 *
 * \brief Set the target point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number
 *
 */

void
PDM_closest_points_tgt_cloud_set
(
       PDM_closest_point_t *cls,
 const int                  i_part,
 const int                  n_points,
       double              *coords,
       PDM_g_num_t         *gnum
);


/**
 *
 * \brief Set the source point cloud
 *
 * \param [in]   id              Identifier
 * \param [in]   i_part          Index of partition
 * \param [in]   n_points        Number of points
 * \param [in]   coords          Point coordinates
 * \param [in]   gnum            Point global number
 *
 */

void
PDM_closest_points_src_cloud_set
(
       PDM_closest_point_t *cls,
 const int                  i_part,
 const int                  n_points,
       double              *coords,
       PDM_g_num_t         *gnum
);

/**
 *
 * \brief Look for closest points
 *
 * \param [in]   id  Identifier
 *
 */

void
PDM_closest_points_compute
(
 PDM_closest_point_t *cls
);


/**
 *
 * \brief Get mesh distance
 *
 * \param [in]   id                    Identifier
 * \param [in]   i_part_tgt            Index of partition of the cloud
 * \param [out]  closest_src_g_num     Global number of the closest element (size = n_closest * n_tgt_points)
 * \param [out]  closest_src_distance  Distance (size = n_closest * n_tgt_points)
 *
 */

void
PDM_closest_points_get
(
       PDM_closest_point_t  *cls,
 const int                   i_part_tgt,
       PDM_g_num_t         **closest_src_gnum,
       double              **closest_src_distance
);


/**
 *
 * \brief Free a distance mesh structure
 *
 * \param [in]  id       Identifier
 * \param [in]  partial  if partial is equal to 0, all data are removed.
 *                       Otherwise, results are kept.
 *
 */

void
PDM_closest_points_free
(
 PDM_closest_point_t  *cls
);


/**
 *
 * \brief  Dump elapsed an CPU time
 *
 * \param [in]  id       Identifier
 *
 */

void
PDM_closest_points_dump_times
(
 PDM_closest_point_t  *cls
);

/**
 *
 * \brief Get mesh distance
 *
 * \param [in]   id                 Identifier
 * \param [in]   i_part_src         Index of partition of the cloud
 * \param [out]  tgt_in_src_idx     For each src point the number of target localised  (size = n_src_points )
 * \param [out]  tgt_in_src         For each src point the globla number of target point located (size = tgt_in_src_idx[n_src_points] )
 *
 */

void
PDM_closest_points_tgt_in_src_get
(
       PDM_closest_point_t  *cls,
 const int                   i_part_src,
       int                 **tgt_in_src_idx,
       PDM_g_num_t         **tgt_in_src
);


void
PDM_transform_to_parent_gnum
(
       PDM_g_num_t  *results,
 const int           n_results,
 const PDM_g_num_t  *ln_to_gn,
 const PDM_g_num_t  *parent_ln_to_gn,
 const int           n_elmt,
       PDM_MPI_Comm  comm
);


/**
 *
 * \brief  transfert _closest_pts var as it seems this static var is not readable
 *          when we switch to the nvcc compiler
 *
 */

PDM_closest_point_t*
PDM_closest_points_closest_transfert
(
  PDM_closest_point_t  *cls
);

#ifdef	__cplusplus
}
#endif

#endif // PDM_CLOSEST_POINTS_H
