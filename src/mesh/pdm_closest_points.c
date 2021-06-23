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
 * System headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_timer.h"
#include "pdm_closest_points.h"
#include "pdm_closest_points_priv.h"
#include "pdm_para_octree.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_array.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"
#include "pdm_logging.h"

/*----------------------------------------------------------------------------*/

#ifdef	__cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif


/*============================================================================
 * Macro definitions
 *============================================================================*/


/*============================================================================
 * Type definitions
 *============================================================================*/

#define NTIMER 2

/**
 * \enum _timer_step_t
 *
 */

typedef enum {

  BEGIN    = 0,
  END      = 1,

} _timer_step_t;

/*============================================================================
 * Global variable
 *============================================================================*/

//static int idebug = 0;

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Reverse result
 *
 * \param [in]   id                    Identifier
 *
 */
static void
_closest_points_reverse_results
(
 PDM_closest_point_t  *cls
)
{

  assert (cls->tgt_cloud->closest_src_gnum != NULL);
  assert (cls->tgt_cloud->closest_src_dist != NULL);

  int* n_points = (int * ) malloc( cls->tgt_cloud->n_part * sizeof(int));
  PDM_g_num_t **tgt_g_num   = (PDM_g_num_t ** ) malloc( cls->tgt_cloud->n_part * sizeof(PDM_g_num_t *));
  int         **tgt_g_num_n = (int         ** ) malloc( cls->tgt_cloud->n_part * sizeof(int         *));
  for (int i_part = 0; i_part < cls->tgt_cloud->n_part; i_part++) {
    n_points[i_part] = cls->tgt_cloud->n_points[i_part] * cls->n_closest;
    tgt_g_num  [i_part] = (PDM_g_num_t * ) malloc( n_points[i_part] * sizeof(PDM_g_num_t));
    tgt_g_num_n[i_part] = (int         * ) malloc( n_points[i_part] * sizeof(int        ));

    // PDM_log_trace_array_long(cls->tgt_cloud->closest_src_gnum[i_part], cls->tgt_cloud->n_points[i_part], "cls->tgt_cloud->closest_src_gnum:: " );

    for(int i = 0; i < cls->tgt_cloud->n_points[i_part]; ++i) {
      for(int ii = 0; ii < cls->n_closest; ++ii) {
        int idx = i * cls->n_closest + ii;
        tgt_g_num  [i_part][idx] = cls->tgt_cloud->gnum[i_part][i];
        tgt_g_num_n[i_part][idx] = 1;
      }
    }
  }

  /*
   * Compute the total number of target point to setup properly the part_to_block partial
   */
  PDM_g_num_t n_g_src = 0;
  for (int i_part = 0; i_part < cls->src_cloud->n_part; i_part++) {
    for(int i = 0; i < cls->src_cloud->n_points[i_part]; ++i) {
      n_g_src = PDM_MAX(n_g_src, cls->src_cloud->gnum[i_part][i]);
    }
  }
  PDM_g_num_t _n_g_src = 0;
  PDM_MPI_Allreduce (&n_g_src, &_n_g_src, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, cls->comm);
  n_g_src = _n_g_src;

  /*
   *  First part to block to map in global numbering of src all target associate
   */
  PDM_part_to_block_t *ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_MERGE,
                                                      1.,
                                                      cls->tgt_cloud->closest_src_gnum,
                                                      NULL,
                                                      n_points,
                                                      cls->tgt_cloud->n_part,
                                                      cls->comm);

  /*
   * For each target (associated to a src gnum prepare an exhange of the current target g_num)
   */
  int *block_tgt_in_src_n = NULL;
  PDM_g_num_t *block_tgt_in_src_g_num = NULL;
  int blk_size = PDM_part_to_block_exch (ptb,
                                         sizeof(PDM_g_num_t),
                                         PDM_STRIDE_VAR,
                                        1,
                                        tgt_g_num_n,
                              (void **) tgt_g_num,
                                        &block_tgt_in_src_n,
                              (void **) &block_tgt_in_src_g_num);

  if(0 == 1) {
    int block_n_elt = PDM_part_to_block_n_elt_block_get (ptb);
    PDM_log_trace_array_int(block_tgt_in_src_n     , block_n_elt, "block_tgt_in_src_n:: " );
    PDM_log_trace_array_long(block_tgt_in_src_g_num, blk_size   , "block_tgt_in_src_g_num:: " );
  }

  for (int i_part = 0; i_part < cls->tgt_cloud->n_part; i_part++) {
    free(tgt_g_num  [i_part]);
    free(tgt_g_num_n[i_part]);
  }
  free(tgt_g_num  );
  free(tgt_g_num_n);
  free(n_points);

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(cls->comm, &i_rank);
  PDM_MPI_Comm_size(cls->comm, &n_rank);

  PDM_g_num_t *_block_distrib_idx = PDM_part_to_block_adapt_partial_block_to_block(ptb,
                                                                                   &block_tgt_in_src_n, /* Realloc inside */
                                                                                   n_g_src);

  PDM_part_to_block_free(ptb);

  PDM_block_to_part_t *btp = PDM_block_to_part_create(_block_distrib_idx,
                               (const PDM_g_num_t **) cls->src_cloud->gnum,
                                                      cls->src_cloud->n_points,
                                                      cls->src_cloud->n_part,
                                                      cls->comm);


  int** tgt_in_src_n;
  PDM_block_to_part_exch2(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR,
                          block_tgt_in_src_n,
                          block_tgt_in_src_g_num,
                         &tgt_in_src_n,
              (void ***) &cls->src_cloud->tgt_in_src);


  cls->src_cloud->tgt_in_src_idx = (int **) malloc( cls->src_cloud->n_part * sizeof(int *));
  for (int i_part = 0; i_part < cls->src_cloud->n_part; i_part++) {
    // PDM_log_trace_array_int(tgt_in_src_n[i_part]     , cls->src_cloud->n_points[i_part], "cls->src_cloud->n_points[i_part]:: " );
    cls->src_cloud->tgt_in_src_idx[i_part] = PDM_array_new_idx_from_sizes_int(tgt_in_src_n[i_part], cls->src_cloud->n_points[i_part]);
    free(tgt_in_src_n[i_part]);
  }
  free(tgt_in_src_n);

  PDM_block_to_part_free(btp);
  free(_block_distrib_idx);
  free(block_tgt_in_src_n);
  free(block_tgt_in_src_g_num);
}


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
)
{
  PDM_closest_point_t *closest = (PDM_closest_point_t *) malloc(sizeof(PDM_closest_point_t));

  closest->comm                         = comm;
  closest->owner                        = owner;
  closest->results_is_getted            = PDM_FALSE;
  closest->tgt_in_src_results_is_getted = PDM_FALSE;

  closest->n_closest = n_closest;
  closest->src_cloud = NULL;
  closest->tgt_cloud = NULL;

  closest->timer = PDM_timer_create ();

  for (int i = 0; i < NTIMER; i++) {
    closest->times_elapsed[i] = 0.;
    closest->times_cpu[i] = 0.;
    closest->times_cpu_u[i] = 0.;
    closest->times_cpu_s[i] = 0.;
  }

  return closest;
}

PDM_closest_point_t*
PDM_closest_points_create_cf
(
 const PDM_MPI_Fint     comm,
 const int              n_closest,
 const PDM_ownership_t  owner
)
{
  const PDM_MPI_Comm _comm        = PDM_MPI_Comm_f2c(comm);

  return PDM_closest_points_create(_comm, n_closest, owner);
}


/**
 *
 * \brief Set the number of partitions of a point cloud
 *
 * \param [in]   id                Identifier
 * \param [in]   n_part_cloud_src  Number of partitions of the source cloud
 * \param [in]   n_part_cloud_tgt  Number of partitions of the target cloud
 *
 */

void
PDM_closest_points_n_part_cloud_set
(
       PDM_closest_point_t* cls,
 const int                  n_part_cloud_src,
 const int                  n_part_cloud_tgt
)
{
  assert(cls->src_cloud == NULL);
  assert(cls->tgt_cloud == NULL);

  cls->src_cloud = malloc (sizeof(_src_point_cloud_t));
  cls->tgt_cloud = malloc (sizeof(_tgt_point_cloud_t));

  cls->src_cloud->n_part   = n_part_cloud_src;
  cls->src_cloud->coords   = malloc (sizeof(double      *) * n_part_cloud_src);
  cls->src_cloud->gnum     = malloc (sizeof(PDM_g_num_t *) * n_part_cloud_src);
  cls->src_cloud->n_points = malloc (sizeof(int          ) * n_part_cloud_src);

  cls->src_cloud->tgt_in_src_idx = NULL;
  cls->src_cloud->tgt_in_src     = NULL;

  cls->tgt_cloud->n_part           = n_part_cloud_tgt;
  cls->tgt_cloud->coords           = malloc (sizeof(double      *) * n_part_cloud_tgt);
  cls->tgt_cloud->gnum             = malloc (sizeof(PDM_g_num_t *) * n_part_cloud_tgt);
  cls->tgt_cloud->n_points         = malloc (sizeof(int          ) * n_part_cloud_tgt);
  cls->tgt_cloud->closest_src_gnum = NULL;
  cls->tgt_cloud->closest_src_dist = NULL;
}


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
)
{
  assert(cls->tgt_cloud != NULL);
  cls->tgt_cloud->n_points[i_part] = n_points;
  cls->tgt_cloud->coords[i_part] = coords;
  cls->tgt_cloud->gnum[i_part] = gnum;
}


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
)
{
  assert(cls->src_cloud != NULL);
  cls->src_cloud->n_points[i_part] = n_points;
  cls->src_cloud->coords  [i_part] = coords;
  cls->src_cloud->gnum    [i_part] = gnum;
}

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
)
{

  cls->times_elapsed[BEGIN] = PDM_timer_elapsed(cls->timer);
  cls->times_cpu[BEGIN]     = PDM_timer_cpu(cls->timer);
  cls->times_cpu_u[BEGIN]   = PDM_timer_cpu_user(cls->timer);
  cls->times_cpu_s[BEGIN]   = PDM_timer_cpu_sys(cls->timer);

  PDM_timer_resume(cls->timer);


  int i_rank;
  PDM_MPI_Comm_rank (cls->comm, &i_rank);

  const int depth_max = 31;
  const int points_in_leaf_max = cls->n_closest;
  const int build_leaf_neighbours = 1;


  /* Create empty parallel octree structure */
  int octree_id = PDM_para_octree_create (cls->src_cloud->n_part,
                                          depth_max,
                                          points_in_leaf_max,
                                          build_leaf_neighbours,
                                          cls->comm);


  /* Set source point clouds */
  for (int i_part = 0; i_part < cls->src_cloud->n_part; i_part++) {
    PDM_para_octree_point_cloud_set (octree_id,
                                     i_part,
                                     cls->src_cloud->n_points[i_part],
                                     cls->src_cloud->coords[i_part],
                                     cls->src_cloud->gnum[i_part]);
  }

  /* Build parallel octree */
  PDM_para_octree_build (octree_id, NULL);
  //PDM_para_octree_dump (octree_id);
  PDM_para_octree_dump_times (octree_id);
  //<--


  /* Concatenate partitions */
  int n_tgt = 0;
  for (int i_part = 0; i_part < cls->tgt_cloud->n_part; i_part++)
    n_tgt += cls->tgt_cloud->n_points[i_part];

  double      *tgt_coord = malloc (sizeof(double)      * n_tgt * 3);
  PDM_g_num_t *tgt_g_num = malloc (sizeof(PDM_g_num_t) * n_tgt);
  PDM_g_num_t *closest_src_gnum = malloc (sizeof(PDM_g_num_t) * n_tgt * cls->n_closest);
  double      *closest_src_dist = malloc (sizeof(double)      * n_tgt * cls->n_closest);

  n_tgt = 0;
  for (int i_part = 0; i_part < cls->tgt_cloud->n_part; i_part++) {
    for (int i = 0; i < cls->tgt_cloud->n_points[i_part]; i++) {
      for (int j = 0; j < 3; j++)
        tgt_coord[3*(n_tgt + i) + j] = cls->tgt_cloud->coords[i_part][3*i + j];
      tgt_g_num[n_tgt + i] = cls->tgt_cloud->gnum[i_part][i];
    }
    n_tgt += cls->tgt_cloud->n_points[i_part];
  }

  /* Search closest source points from target points */
  if (cls->n_closest == 1) {
    PDM_para_octree_single_closest_point (octree_id,
                                          n_tgt,
                                          tgt_coord,
                                          tgt_g_num,
                                          closest_src_gnum,
                                          closest_src_dist);
  } else {
    PDM_para_octree_closest_points (octree_id,
                                    cls->n_closest,
                                    n_tgt,
                                    tgt_coord,
                                    tgt_g_num,
                                    closest_src_gnum,
                                    closest_src_dist);
  }


  // PDM_log_trace_array_long(tgt_g_num, n_tgt, "tgt_g_num:: " );
  // PDM_log_trace_array_double(tgt_coord, 3 * n_tgt, "tgt_coord:: " );
  // PDM_log_trace_array_long(closest_src_gnum, n_tgt * cls->n_closest, "closest_src_gnum:: " );
  // PDM_log_trace_array_double(closest_src_dist, n_tgt * cls->n_closest, "closest_src_dist:: " );

  /* Restore partitions */
  free (tgt_coord);
  free (tgt_g_num);
  n_tgt = 0;

  cls->tgt_cloud->closest_src_gnum = malloc (sizeof(PDM_g_num_t *) * cls->tgt_cloud->n_part);
  cls->tgt_cloud->closest_src_dist = malloc (sizeof(double *)      * cls->tgt_cloud->n_part);

  for (int i_part = 0; i_part < cls->tgt_cloud->n_part; i_part++) {
    int s_closest_src = cls->n_closest * cls->tgt_cloud->n_points[i_part];

    cls->tgt_cloud->closest_src_gnum[i_part] = malloc (sizeof(PDM_g_num_t) * s_closest_src);
    cls->tgt_cloud->closest_src_dist[i_part] = malloc (sizeof(double)      * s_closest_src);

    for (int i = 0; i < cls->tgt_cloud->n_points[i_part]; i++) {
      for (int j = 0; j < cls->n_closest; j++) {
        cls->tgt_cloud->closest_src_gnum[i_part][cls->n_closest*i+j] =
          closest_src_gnum[n_tgt + cls->n_closest*i + j];

        cls->tgt_cloud->closest_src_dist[i_part][cls->n_closest*i+j] =
          closest_src_dist[n_tgt + cls->n_closest*i + j];
      }
    }
    n_tgt += cls->n_closest * cls->tgt_cloud->n_points[i_part];
  }
  free (closest_src_gnum);
  free (closest_src_dist);



  //-->GPU
  /* Free parallel octree */
  PDM_para_octree_free (octree_id);
  //<--

  _closest_points_reverse_results(cls);


  PDM_timer_hang_on(cls->timer);

  cls->times_elapsed[END] = PDM_timer_elapsed(cls->timer);
  cls->times_cpu[END]     = PDM_timer_cpu(cls->timer);
  cls->times_cpu_u[END]   = PDM_timer_cpu_user(cls->timer);
  cls->times_cpu_s[END]   = PDM_timer_cpu_sys(cls->timer);

  PDM_timer_resume(cls->timer);
}


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
)
{

  assert (cls->tgt_cloud->closest_src_gnum != NULL);
  assert (cls->tgt_cloud->closest_src_dist != NULL);

  *closest_src_gnum     = cls->tgt_cloud->closest_src_gnum[i_part_tgt];
  *closest_src_distance = cls->tgt_cloud->closest_src_dist[i_part_tgt];

  cls->results_is_getted = PDM_TRUE;
}


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
)
{

  assert (cls->src_cloud->tgt_in_src_idx != NULL);
  assert (cls->src_cloud->tgt_in_src != NULL);

  *tgt_in_src_idx = cls->src_cloud->tgt_in_src_idx[i_part_src];
  *tgt_in_src     = cls->src_cloud->tgt_in_src    [i_part_src];

  // int size = cls->src_cloud->n_points[i_part_src];
  // PDM_log_trace_array_long(cls->src_cloud->tgt_in_src_idx[i_part_src], size, "get -> tgt_in_src_idx :: " );
  // PDM_log_trace_array_long(cls->src_cloud->tgt_in_src    [i_part_src], (*tgt_in_src_idx)[size], "get -> tgt_in_src :: " );

  cls->tgt_in_src_results_is_getted = PDM_TRUE;
}



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
)
{

  if(( cls->owner == PDM_OWNERSHIP_KEEP ) ||
     ( cls->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !cls->results_is_getted)){
    if (cls->tgt_cloud->closest_src_gnum != NULL) {
      for (int j = 0; j < cls->tgt_cloud->n_part ; j++) {
        if (cls->tgt_cloud->closest_src_gnum[j] != NULL) {
          free (cls->tgt_cloud->closest_src_gnum[j]);
        }
      }
    }

    if (cls->tgt_cloud->closest_src_dist != NULL) {
      for (int j = 0; j < cls->tgt_cloud->n_part ; j++) {
        if (cls->tgt_cloud->closest_src_dist[j] != NULL) {
          free (cls->tgt_cloud->closest_src_dist[j]);
        }
      }
    }
  }

  free (cls->tgt_cloud->closest_src_gnum);
  free (cls->tgt_cloud->closest_src_dist);

  if(( cls->owner == PDM_OWNERSHIP_KEEP ) ||
     ( cls->owner == PDM_OWNERSHIP_UNGET_RESULT_IS_FREE && !cls->tgt_in_src_results_is_getted)||
     ( cls->owner == PDM_OWNERSHIP_USER                 && !cls->tgt_in_src_results_is_getted)){ // Dernière condition pour le python essentiellement ou si un utilisateur n'a pas besoin de ce résultats

    if (cls->src_cloud->tgt_in_src_idx != NULL) {
      for (int j = 0; j < cls->src_cloud->n_part ; j++) {
        if (cls->src_cloud->tgt_in_src_idx[j] != NULL) {
          free (cls->src_cloud->tgt_in_src_idx[j]);
          free (cls->src_cloud->tgt_in_src[j]);
        }
      }
    }
  }

  free (cls->src_cloud->tgt_in_src_idx);
  free (cls->src_cloud->tgt_in_src);

  if (cls->tgt_cloud->gnum != NULL) {
    free (cls->tgt_cloud->gnum);
  }
  if (cls->tgt_cloud->coords != NULL) {
    free (cls->tgt_cloud->coords);
  }
  if (cls->tgt_cloud->n_points != NULL) {
    free (cls->tgt_cloud->n_points);
  }
  if (cls->tgt_cloud != NULL) {
    free (cls->tgt_cloud);
  }


  if (cls->src_cloud->gnum != NULL) {
    free (cls->src_cloud->gnum);
  }
  if (cls->src_cloud->coords != NULL) {
    free (cls->src_cloud->coords);
  }
  if (cls->src_cloud->n_points != NULL) {
    free (cls->src_cloud->n_points);
  }
  if (cls->src_cloud != NULL) {
    free (cls->src_cloud);
  }

  PDM_timer_free(cls->timer);

  free (cls);

}


/**
 *
 * \brief  Dump elapsed and CPU time
 *
 * \param [in]  id       Identifier
 *
 */

void
PDM_closest_points_dump_times
(
PDM_closest_point_t  *cls
)
{
  double t1 = cls->times_elapsed[END] - cls->times_elapsed[BEGIN];
  double t2 = cls->times_cpu[END] - cls->times_cpu[BEGIN];

  double t1max;
  PDM_MPI_Allreduce (&t1, &t1max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, cls->comm);

  double t2max;
  PDM_MPI_Allreduce (&t2, &t2max, 1, PDM_MPI_DOUBLE, PDM_MPI_MAX, cls->comm);

  int rank;
  PDM_MPI_Comm_rank (cls->comm, &rank);

  if (rank == 0) {

    PDM_printf( "closest_points timer : all (elapsed and cpu) : %12.5es %12.5es\n",
                t1max, t2max);
  }
}


void
PDM_transform_to_parent_gnum
(
       PDM_g_num_t  *results,
 const int           n_results,
 const PDM_g_num_t  *ln_to_gn,
 const PDM_g_num_t  *parent_ln_to_gn,
 const int           n_elmt,
       PDM_MPI_Comm  comm
)
{
  PDM_part_to_block_t* ptb = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                      PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                      1.,
                                ( PDM_g_num_t **)    &ln_to_gn,
                                                      NULL,
                                          ( int *)   &n_elmt,
                                                      1,
                                                      comm);

  int         *block_stride = NULL;
  PDM_g_num_t *block_parent = NULL;
  int s_block_data = PDM_part_to_block_exch (ptb,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          1,
                          NULL,
               (void **) &parent_ln_to_gn,
                         &block_stride,
               (void **) &block_parent);

  PDM_g_num_t *block_distrib_idx = PDM_part_to_block_distrib_index_get (ptb);

  if(0 == 1){
    PDM_log_trace_array_long(block_parent, s_block_data, "block_parent :: " );
  }

  PDM_block_to_part_t *btp = PDM_block_to_part_create(block_distrib_idx,
                               (const PDM_g_num_t **) &results,
                                                      &n_results,
                                                      1,
                                                      comm);

  int stride_one = 1;
  PDM_block_to_part_exch(btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_CST,
                          &stride_one,
                          block_parent,
                          NULL,
              (void **) &results);


  PDM_part_to_block_free(ptb);
  PDM_block_to_part_free(btp);
  free(block_parent);
}

/**
 *
 * \brief  transfert _closest_pts var as it seems this static var is not readable
 *          when we switch to the nvcc compiler
 *
 */

PDM_closest_point_t *
PDM_closest_points_closest_transfert
(
  PDM_closest_point_t  *cls
)
{
  return cls;
}

#ifdef	__cplusplus
}
#endif
#undef NTIMER
