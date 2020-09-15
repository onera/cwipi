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

int
PDM_closest_points_create
(
 const PDM_MPI_Comm comm,
 const int          n_closest
);

void
PDM_closest_points_create_cf 
(
 const PDM_MPI_Fint comm,
 const int          n_closest,
 int *id
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
 const int  id,
 const int  n_part_cloud_src,
 const int  n_part_cloud_tgt
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
 const int          id,
 const int          i_part,
 const int          n_points,
       double      *coords,
       PDM_g_num_t *gnum
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
 const int          id,
 const int          i_part,
 const int          n_points,
       double      *coords,
       PDM_g_num_t *gnum
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
 const int id
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
 const int        id,
 const int        i_part_tgt,
 PDM_g_num_t    **closest_src_gnum,
       double   **closest_src_distance
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
 const int id,
 const int partial
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
 const int id
 );

#ifdef	__cplusplus
}
#endif

#endif // PDM_CLOSEST_POINTS_H
