/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/
#include <sys/resource.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_doctree_priv.h"
#include "pdm_doctree.h"
#include "pdm_timer.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_sort.h"
#include "pdm_box.h"
#include "pdm_box_tree.h"
#include "pdm_box_priv.h"
#include "pdm_array.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_binary_search.h"
#include "pdm_octree_seq.h"
#include "pdm_kdtree_seq.h"
#include "pdm_vtk.h"


#ifdef __cplusplus
extern "C"
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Static function definitions
 *============================================================================*/

static
void
_redistribute_pts_geom
(
 PDM_doctree_t        *doct,
 PDM_part_to_block_t **ptb_out,
 double              **dpts_coords_out
)
{

  /*
   * Redistribute all pts and impose hilbert ordering
   */
  int **weight = malloc(doct->n_part_cloud * sizeof(int *));
  for(int i_part = 0; i_part < doct->n_part_cloud; ++i_part) {
    weight[i_part] = malloc(doct->n_point_cloud[i_part] * sizeof(int));
    for(int i = 0; i < doct->n_point_cloud[i_part]; ++i) {
      weight[i_part][i] = 1;
    }
  }


  PDM_part_to_block_t* ptb = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_CLEANUP, // A voir avec merge mais attention au init_location
                                                           1.,
                                                           PDM_PART_GEOM_HILBERT,
                                                           doct->pts_coords,
                                                           doct->pts_g_num,
                                                           weight,
                                                           doct->n_point_cloud,
                                                           doct->n_part_cloud,
                                                           doct->comm);

  for(int i_part = 0; i_part < doct->n_part_cloud; ++i_part) {
    free(weight[i_part]);
  }
  free(weight);

  /* Il faut le faire en MERGE mais enlever les doublons de coords */
  double *blk_pts_coord = NULL;
  PDM_part_to_block_exch(ptb,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
               (void **) doct->pts_coords,
                         NULL,
               (void **) &blk_pts_coord);


  /* Transport init_location - Attention au merge du ptb à faire */
  int have_init_location = 1;
  for(int i = 0; i < doct->n_part_cloud; ++i) {
    if(doct->pts_init_location[i] == NULL) {
      have_init_location = 0;
    }
  }

  if(have_init_location == 1) {
    int *blk_init_location_pts_n = NULL;
    int *blk_init_location_pts   = NULL;

    int **stride_one = malloc(doct->n_part_cloud * sizeof(int *));
    for(int i_part = 0; i_part < doct->n_part_cloud; ++i_part) {
      stride_one[i_part] = malloc(doct->n_point_cloud[i_part] * sizeof(int));
      for(int i = 0; i < doct->n_point_cloud[i_part]; ++i) {
        stride_one[i_part][i] = 1;
      }
    }

    PDM_part_to_block_exch(ptb,
                           3 * sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           1,
                           stride_one,
                 (void **) doct->pts_init_location,
                           &blk_init_location_pts_n,
                 (void **) &blk_init_location_pts);

    for(int i_part = 0; i_part < doct->n_part_cloud; ++i_part) {
      free(stride_one[i_part]);
    }
    free(stride_one);


    free(blk_init_location_pts_n);
    free(blk_init_location_pts);

  }


  // int          n_parent    = PDM_part_to_block_n_elt_block_get  (ptb);
  // PDM_g_num_t* parent_gnum = PDM_part_to_block_block_gnum_get   (ptb);
  // PDM_g_num_t* distrib_pts = PDM_part_to_block_distrib_index_get(ptb);

  *ptb_out         = ptb;
  *dpts_coords_out = blk_pts_coord;

}

// void
// _exchange_global_tree
// (
//  PDM_doctree_t        *doct,
//  int                   n_coarse_box,
//  double               *coarse_box_extents,
//  int                  *box_n_pts
// )
// {

// }


/*=============================================================================
 * Public function definitions
 *============================================================================*/

PDM_doctree_t*
PDM_doctree_create
(
 PDM_MPI_Comm              comm,
 int                       dim,
 int                       n_part_cloud,
 double                   *global_extents,
 PDM_doctree_local_tree_t  local_tree_kind
)
{
  PDM_doctree_t* doct = (PDM_doctree_t *) malloc(sizeof(PDM_doctree_t));

  doct->comm = comm;
  doct->dim  = dim;

  doct->local_tree_kind = local_tree_kind;

  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE) {
    doct->coarse_depth_max          = 5;
    doct->coarse_points_in_leaf_max = 60;

    doct->local_depth_max          = 31;
    doct->local_points_in_leaf_max = 30;
    doct->local_tolerance          = 1e-6;

  } else if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE){
    doct->coarse_depth_max          = 6;
    doct->coarse_points_in_leaf_max = 60;

    doct->local_depth_max          = 31;
    doct->local_points_in_leaf_max = 30;
    doct->local_tolerance          = 1e-6;
  } else {
    abort();
  }

  doct->coarse_octree   = NULL;
  doct->local_octree    = NULL;
  doct->shmem_octree    = NULL;

  doct->coarse_kdtree   = NULL;
  doct->local_kdtree    = NULL;
  doct->shmem_kdtree    = NULL;

  doct->comm_dist_graph = PDM_MPI_COMM_NULL;
  doct->n_degree_in     = 0;
  doct->neighbor_in     = NULL;

  doct->comm_shared     = PDM_MPI_COMM_NULL;

  PDM_UNUSED(global_extents);

  doct->n_part_cloud      = n_part_cloud;
  doct->n_point_cloud     = malloc(n_part_cloud * sizeof(int          ));
  doct->pts_g_num         = malloc(n_part_cloud * sizeof(PDM_g_num_t *));
  doct->pts_coords        = malloc(n_part_cloud * sizeof(double      *));
  doct->pts_init_location = malloc(n_part_cloud * sizeof(int         *));

  for(int i = 0; i < n_part_cloud; ++i) {
    doct->n_point_cloud    [i] = 0;
    doct->pts_g_num        [i] = NULL;
    doct->pts_coords       [i] = NULL;
    doct->pts_init_location[i] = NULL;
  }

  return doct;
}


void
PDM_doctree_build
(
 PDM_doctree_t     *doct
)
{
  int n_rank, i_rank;
  PDM_MPI_Comm_rank (doct->comm, &i_rank);
  PDM_MPI_Comm_size (doct->comm, &n_rank);

  /*
   * Prepare graphe comm for hybrid MPI-MPI
   */
  PDM_MPI_setup_hybrid_dist_comm_graph(doct->comm,
                                       &doct->comm_shared,
                                       &doct->comm_dist_graph,
                                       &doct->n_degree_in,
                                       &doct->neighbor_in);

  int n_rank_in_shm, i_rank_in_shm;
  PDM_MPI_Comm_rank (doct->comm_shared, &i_rank_in_shm);
  PDM_MPI_Comm_size (doct->comm_shared, &n_rank_in_shm);

  /*
   * Redistribute all pts
   */
  double              *blk_pts_coord = NULL;
  PDM_part_to_block_t *ptb           = NULL;
  _redistribute_pts_geom(doct, &ptb, &blk_pts_coord);

  PDM_g_num_t* distrib_pts = PDM_part_to_block_distrib_index_get(ptb);

  /*
   * Step 2 : Build local coarse tree
   */
  int dn_pts = distrib_pts[i_rank+1] - distrib_pts[i_rank];
  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE) {

    assert(doct->coarse_octree == NULL);
    doct->coarse_octree = PDM_octree_seq_create(1, // n_point_cloud
                                          doct->coarse_depth_max,
                                          doct->coarse_points_in_leaf_max,
                                          doct->local_tolerance);

    PDM_octree_seq_point_cloud_set(doct->coarse_octree,
                                   0,
                                   dn_pts,
                                   blk_pts_coord);

    PDM_octree_seq_build(doct->coarse_octree);


  } else if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE){

    assert(doct->coarse_kdtree == NULL);
    doct->coarse_kdtree = PDM_kdtree_seq_create(1, // n_point_cloud
                                                doct->coarse_depth_max,
                                                doct->coarse_points_in_leaf_max,
                                                doct->local_tolerance);

    PDM_kdtree_seq_point_cloud_set(doct->coarse_kdtree,
                                   0,
                                   dn_pts,
                                   blk_pts_coord);

    PDM_kdtree_seq_build(doct->coarse_kdtree);
  }

  /*
   * Extract extents on all local_tree
   */
  int n_coarse_box = 0;
  double *coarse_box_extents = NULL;
  int    *coarse_box_n_pts   = NULL; // Number of point in boxes
  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE) {
    int n_depth_per_proc = 2;
    PDM_octree_seq_extract_extent(doct->coarse_octree,
                                  0,
                                  n_depth_per_proc,
                                  &n_coarse_box,
                                  &coarse_box_extents,
                                  &coarse_box_n_pts);
  } else if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE){
    int n_depth_per_proc = 6; // 2^(depth)
    PDM_kdtree_seq_extract_extent(doct->coarse_kdtree,
                                  0,
                                  n_depth_per_proc,
                                  &n_coarse_box,
                                  &coarse_box_extents,
                                  &coarse_box_n_pts);
  }

  /*
   * Equilibrate among nodes/numa - To reduce memory footprint we set up data in shared memory
   */
  int *lrecv_count = malloc(doct->n_degree_in * sizeof(int));
  PDM_MPI_Neighbor_allgather(&n_coarse_box , 1, PDM_MPI_INT,
                             lrecv_count   , 1, PDM_MPI_INT, doct->comm_dist_graph);

  PDM_mpi_win_shared_t* wshared_local_nodes_n   = PDM_mpi_win_shared_create(n_rank  , sizeof(int), doct->comm_shared);
  PDM_mpi_win_shared_t* wshared_local_nodes_idx = PDM_mpi_win_shared_create(n_rank+1, sizeof(int), doct->comm_shared);
  int *shared_local_nodes_n   = PDM_mpi_win_shared_get(wshared_local_nodes_n);
  int *shared_local_nodes_idx = PDM_mpi_win_shared_get(wshared_local_nodes_idx);
  PDM_mpi_win_shared_lock_all (0, wshared_local_nodes_n  );
  PDM_mpi_win_shared_lock_all (0, wshared_local_nodes_idx);

  for(int i = 0; i < doct->n_degree_in; ++i) {
    shared_local_nodes_n[doct->neighbor_in[i]] = lrecv_count[i];
  }
  PDM_MPI_Barrier(doct->comm_shared);

  if(i_rank_in_shm == 0) {
    shared_local_nodes_idx[0] = 0;
    for(int i = 0; i < n_rank; ++i) {
      shared_local_nodes_idx[i+1] = shared_local_nodes_idx[i] + shared_local_nodes_n[i];
    }
  }
  PDM_MPI_Barrier(doct->comm_shared);

  if(1 == 1) {
    PDM_log_trace_array_int(shared_local_nodes_n  , n_rank , "shared_local_nodes_n   ::");
    PDM_log_trace_array_int(shared_local_nodes_idx, n_rank+1 , "shared_local_nodes_idx ::");
    PDM_log_trace_array_int(doct->neighbor_in, doct->n_degree_in , "doct->neighbor_in ::");
  }

  // Hook local recv_shift
  int *recv_shift = malloc(doct->n_degree_in * sizeof(int));
  for(int i = 0; i < doct->n_degree_in; ++i) {
    recv_shift[i] = shared_local_nodes_idx[doct->neighbor_in[i]];
  }

  /*
   * Exchange extents
   */
  PDM_mpi_win_shared_t* wshared_coarse_box_n_pts   = PDM_mpi_win_shared_create(    shared_local_nodes_idx[n_rank], sizeof(int)   , doct->comm_shared);
  PDM_mpi_win_shared_t* wshared_coarse_box_extents = PDM_mpi_win_shared_create(6 * shared_local_nodes_idx[n_rank], sizeof(double), doct->comm_shared);
  int    *shared_coarse_box_n_pts   = PDM_mpi_win_shared_get(wshared_coarse_box_n_pts  );
  double *shared_coarse_box_extents = PDM_mpi_win_shared_get(wshared_coarse_box_extents);
  PDM_mpi_win_shared_lock_all (0, wshared_coarse_box_n_pts  );
  PDM_mpi_win_shared_lock_all (0, wshared_coarse_box_extents);

  PDM_MPI_Neighbor_allgatherv(coarse_box_n_pts       , n_coarse_box, PDM_MPI_INT,
                              shared_coarse_box_n_pts, lrecv_count  , recv_shift, PDM_MPI_INT, doct->comm_dist_graph);
  PDM_MPI_Barrier(doct->comm_shared);


  /* Update */
  for(int i = 0; i < doct->n_degree_in; ++i) {
    recv_shift [i]  = shared_local_nodes_idx[doct->neighbor_in[i]]*6;
    lrecv_count[i] *= 6;
  }

  PDM_MPI_Neighbor_allgatherv(coarse_box_extents, 6 * n_coarse_box, PDM_MPI_DOUBLE,
                              shared_coarse_box_extents, lrecv_count        , recv_shift, PDM_MPI_DOUBLE, doct->comm_dist_graph);
  PDM_MPI_Barrier(doct->comm_shared);

  /*
   * Solicitate
   */
  int n_shared_boxes = shared_local_nodes_idx[n_rank];
  PDM_g_num_t* distrib_shared_boxes = PDM_compute_uniform_entity_distribution(doct->comm_shared, n_shared_boxes);

  PDM_box_set_t  *box_set   = NULL;
  PDM_box_tree_t *bt_shared = NULL;
  int   max_boxes_leaf_shared = 10; // Max number of boxes in a leaf for coarse shared BBTree
  int   max_tree_depth_shared = 6; // Max tree depth for coarse shared BBTree
  float max_box_ratio_shared  = 5; // Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes)

  PDM_MPI_Comm comm_alone;
  PDM_MPI_Comm_split(doct->comm, i_rank, 0, &(comm_alone));
  const int n_info_location = 3;
  int *init_location_proc = PDM_array_zeros_int (n_info_location * n_shared_boxes);

  PDM_mpi_win_shared_t* wshared_coarse_boxes_gnum = PDM_mpi_win_shared_create(    n_shared_boxes, sizeof(PDM_g_num_t), doct->comm_shared);
  PDM_mpi_win_shared_t* wshared_coarse_box_center        = PDM_mpi_win_shared_create(3 * n_shared_boxes, sizeof(double     ), doct->comm_shared);
  PDM_g_num_t    *shared_coarse_boxes_gnum   = PDM_mpi_win_shared_get(wshared_coarse_boxes_gnum  );
  double         *shared_box_center          = PDM_mpi_win_shared_get(wshared_coarse_box_center  );

  PDM_mpi_win_shared_lock_all (0, wshared_coarse_boxes_gnum  );
  PDM_mpi_win_shared_lock_all (0, wshared_coarse_box_center  );
  for (int i = distrib_shared_boxes[i_rank_in_shm]; i < distrib_shared_boxes[i_rank_in_shm+1]; i++) {
    shared_coarse_boxes_gnum[i] = i + 1;

    shared_box_center[3*i  ] = 0.5 * (shared_coarse_box_extents[6*i  ] + shared_coarse_box_extents[6*i+3]);
    shared_box_center[3*i+1] = 0.5 * (shared_coarse_box_extents[6*i+1] + shared_coarse_box_extents[6*i+4]);
    shared_box_center[3*i+2] = 0.5 * (shared_coarse_box_extents[6*i+2] + shared_coarse_box_extents[6*i+5]);

  }
  PDM_MPI_Barrier(doct->comm_shared);


  // double *shared_coarse_box_extents = PDM_mpi_win_shared_get(wshared_coarse_box_extents);
  // PDM_mpi_win_shared_lock_all (0, wshared_coarse_box_extents);
  int dn_box_shared = distrib_shared_boxes[i_rank_in_shm+1] - distrib_shared_boxes[i_rank_in_shm];
  box_set = PDM_box_set_create(3,
                               0,  // No normalization to preserve initial extents
                               0,  // No projection to preserve initial extents
                               n_shared_boxes,
                               shared_coarse_boxes_gnum,
                               shared_coarse_box_extents,
                               1,
                               &n_shared_boxes,
                               init_location_proc,
                               comm_alone);

  bt_shared = PDM_box_tree_create (max_tree_depth_shared,
                                   max_boxes_leaf_shared,
                                   max_box_ratio_shared);

  PDM_box_tree_set_boxes (bt_shared,
                          box_set,
                          PDM_BOX_TREE_ASYNC_LEVEL);

  free(init_location_proc);

  int* coarse_tree_box_to_box_idx = NULL;
  int* coarse_tree_box_to_box     = NULL;
  if(doct->solicitation_kind == PDM_TREE_SOLICITATION_BOXES_POINTS) {
    assert(doct->n_part == 1);
    PDM_box_tree_intersect_boxes_boxes2(bt_shared,
                                        -1,
                                        doct->n_entity[0],
                                        doct->entity_coords[0],
                                        &coarse_tree_box_to_box_idx,
                                        &coarse_tree_box_to_box);

    if(1 == 1) {
      PDM_log_trace_connectivity_int(coarse_tree_box_to_box_idx,
                                     coarse_tree_box_to_box,
                                     n_shared_boxes,
                                     "coarse_tree_box_to_box : ");
    }

  } else {
    abort();
  }

  /*
   * Pour chaque shared box on connait le poids de la solitation
   *    -> part_to_block sur les shared
   *  Optim = can shared coarse_boxes_gnum -> distribution par noeuds
   */
  int    *weight = malloc(n_shared_boxes * sizeof(int));
  for(int i = 0; i < n_shared_boxes; ++i) {
    weight[i] = coarse_tree_box_to_box_idx[i+1] - coarse_tree_box_to_box_idx[i];
  }

  /*
   * Equilibrate boxes leaf with solicitation
   */
  PDM_part_to_block_t* ptb_equi_box = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                    PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                                    1.,
                                                                    PDM_PART_GEOM_HILBERT,
                                                                    &shared_box_center,
                                                                    &shared_coarse_boxes_gnum,
                                                                    &weight,
                                                                    &n_shared_boxes,
                                                                    1,
                                                                    doct->comm);

  int          dn_equi_tree     = PDM_part_to_block_n_elt_block_get  (ptb_equi_box);
  PDM_g_num_t* parent_tree_gnum = PDM_part_to_block_block_gnum_get   (ptb_equi_box);
  PDM_g_num_t* distrib_tree     = PDM_part_to_block_distrib_index_get(ptb_equi_box);

  if(1 == 1) {
    PDM_log_trace_array_long(parent_tree_gnum, dn_equi_tree , "parent_tree_gnum :: ");
    PDM_log_trace_array_long(distrib_tree    , n_rank+1, "distrib_tree : ");
  }

  // Exchange n_pts
  // int *equi_box_n_pts = NULL;
  // PDM_part_to_block_exch(ptb,
  //                        sizeof(int),
  //                        PDM_STRIDE_CST_INTERLACED,
  //                        1,
  //                        NULL,
  //              (void **) box_n_pts,
  //                        NULL,
  //              (void **) &equi_box_n_pts);
  // if(1 == 1) {
  //   PDM_log_trace_array_int(equi_box_n_pts, dn_equi_tree, "equi_box_n_pts :");
  // }



  free(weight);

  PDM_mpi_win_shared_unlock_all(wshared_coarse_boxes_gnum);
  PDM_mpi_win_shared_unlock_all(wshared_coarse_box_center);
  PDM_mpi_win_shared_free (wshared_coarse_boxes_gnum);
  PDM_mpi_win_shared_free (wshared_coarse_box_center);

  /*
   * Update permutation of pre-solicitation
   */
  PDM_MPI_Neighbor_allgather(&dn_equi_tree, 1, PDM_MPI_INT,
                             lrecv_count  , 1, PDM_MPI_INT, doct->comm_dist_graph);

  for(int i = 0; i < doct->n_degree_in; ++i) {
    shared_local_nodes_n[doct->neighbor_in[i]] = lrecv_count[i];
  }
  PDM_MPI_Barrier(doct->comm_shared);

  if(i_rank_in_shm == 0) {
    shared_local_nodes_idx[0] = 0;
    for(int i = 0; i < n_rank; ++i) {
      shared_local_nodes_idx[i+1] = shared_local_nodes_idx[i] + shared_local_nodes_n[i];
    }
  }
  PDM_MPI_Barrier(doct->comm_shared);

  // Hook local recv_shift
  for(int i = 0; i < doct->n_degree_in; ++i) {
    recv_shift[i] = shared_local_nodes_idx[doct->neighbor_in[i]];
  }

  /*
   * Exchange dparent_gnum
   */
  PDM_mpi_win_shared_t* wshared_parent_tree_gnum   = PDM_mpi_win_shared_create(    shared_local_nodes_idx[n_rank], sizeof(int)   , doct->comm_shared);
  int    *shared_parent_tree_gnum   = PDM_mpi_win_shared_get(wshared_parent_tree_gnum  );
  PDM_mpi_win_shared_lock_all (0, wshared_parent_tree_gnum  );

  PDM_MPI_Neighbor_allgatherv(parent_tree_gnum       , dn_equi_tree, PDM__PDM_MPI_G_NUM,
                              shared_parent_tree_gnum, lrecv_count , recv_shift, PDM__PDM_MPI_G_NUM, doct->comm_dist_graph);
  PDM_MPI_Barrier(doct->comm_shared);

  if(0 == 1) {
    PDM_log_trace_array_long(shared_parent_tree_gnum, n_shared_boxes, "shared_parent_tree_gnum :");
  }


  /*
   * Setup partitioning
   */
  PDM_g_num_t* impli_distrib_tree = PDM_compute_entity_distribution(doct->comm, n_coarse_box);

  if(1 == 1) {
    PDM_log_trace_array_long(impli_distrib_tree    , n_rank+1, "impli_distrib_tree : ");
  }
  PDM_block_to_part_t* btp = PDM_block_to_part_create(impli_distrib_tree,
                               (const PDM_g_num_t **) &parent_tree_gnum,
                                                      &dn_equi_tree,
                                                      1,
                                                      doct->comm);




  /*
   * Prepare buffer
   */
  int    *extract_box_id    = NULL;
  int n_pts_tot = 0;
  for(int i = 0; i < n_coarse_box; ++i ) {
    int node_id = extract_box_id[i];
    if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE) {
      n_pts_tot += PDM_octree_seq_n_points_get(doct->coarse_octree, node_id);
    } else if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
      abort();
      // n_pts_tot += PDM_kdtree_seq_n_points_get(doct->coarse_kdtree, node_id);
    }
  }

  int dn_blk = impli_distrib_tree[i_rank+1] - impli_distrib_tree[i_rank];
  assert(dn_blk == n_coarse_box);

  double *reorder_blk_coord_send = malloc(3 * n_pts_tot * sizeof(double));

  int idx_write = 0;
  for(int i = 0; i < n_coarse_box; ++i ) {
    int node_id = extract_box_id[i];

    int  n_pts = 0;
    int *point_clouds_id = NULL;
    int *point_indexes   = NULL;

    if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE) {
      n_pts = PDM_octree_seq_n_points_get(doct->coarse_octree, node_id);
      PDM_octree_seq_points_get(doct->coarse_octree, node_id, &point_clouds_id, &point_indexes);
    } else if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE) {
      abort();
      // n_pts = PDM_kdtree_seq_points_get(doct->coarse_kdtree, node_id, &point_clouds_id, &point_indexes);
    }

    for(int i_pt = 0; i_pt < n_pts; ++i_pt) {
      int orig_pt = point_indexes[i_pt];
      reorder_blk_coord_send[3*idx_write  ] = blk_pts_coord[3*orig_pt  ];
      reorder_blk_coord_send[3*idx_write+1] = blk_pts_coord[3*orig_pt+1];
      reorder_blk_coord_send[3*idx_write+2] = blk_pts_coord[3*orig_pt+2];
    }
  }

  /*
   * Coordinates
   */
  double **tmp_equi_pts_coords = NULL;
  int    **tmp_equi_n_pts      = NULL;
  PDM_block_to_part_exch(btp,
                         3 * sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         coarse_box_n_pts,
                         reorder_blk_coord_send,
                         &tmp_equi_n_pts,
              (void ***) &tmp_equi_pts_coords);
  double *equi_pts_coords = tmp_equi_pts_coords[0];
  int    *equi_n_pts      = tmp_equi_n_pts[0];
  free(tmp_equi_pts_coords);
  free(tmp_equi_n_pts);


  int equi_n_pts_tot = 0;
  for(int i = 0; i < dn_equi_tree; ++i) {
    equi_n_pts_tot += equi_n_pts[i];
  }

  /*
   * g_num
   */
  PDM_g_num_t* equi_pts_gnum = NULL;

  /*
   * Init location -> Attention si gnum de points dupliqué -> variable
   */
  PDM_g_num_t* equi_pts_init_location = NULL;



  free(reorder_blk_coord_send);
  PDM_block_to_part_free(btp);
  PDM_part_to_block_free(ptb_equi_box);

  free(distrib_shared_boxes);
  free(impli_distrib_tree);
  free(coarse_tree_box_to_box_idx);
  free(coarse_tree_box_to_box);

  PDM_MPI_Comm_free(&comm_alone);

  PDM_box_set_destroy (&box_set);
  PDM_box_tree_destroy (&bt_shared);

  PDM_mpi_win_shared_unlock_all(wshared_coarse_box_n_pts);
  PDM_mpi_win_shared_unlock_all(wshared_coarse_box_extents);

  PDM_mpi_win_shared_unlock_all (wshared_local_nodes_n  );
  PDM_mpi_win_shared_unlock_all (wshared_local_nodes_idx);

  PDM_mpi_win_shared_free (wshared_local_nodes_n  );
  PDM_mpi_win_shared_free (wshared_local_nodes_idx);

  PDM_mpi_win_shared_free(wshared_coarse_box_n_pts);
  PDM_mpi_win_shared_free(wshared_coarse_box_extents);


  PDM_mpi_win_shared_unlock_all (wshared_parent_tree_gnum  );
  PDM_mpi_win_shared_free (wshared_parent_tree_gnum  );

  free(recv_shift);
  free(lrecv_count);

  free(coarse_box_n_pts);
  free(coarse_box_extents);

  free(blk_pts_coord);

  PDM_part_to_block_free(ptb);

  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE) {
    if(doct->coarse_octree != NULL) {
      PDM_octree_seq_free(doct->coarse_octree);
    }
  } else if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE){
    if(doct->coarse_kdtree != NULL) {
      PDM_kdtree_seq_free(doct->coarse_kdtree);
    }
  }

  // Preparation of send count and box_rank/box_rank_idx
  // for(int i = 0; i < n_rank; ++i) {
  //   for(int j = shared_all_rank_idx[i]; j < shared_all_rank_idx[i+1]; ++j) {
  //     send_count[i] += shared_to_box_idx[j+1] - shared_to_box_idx[j];
  //   }
  // }

  /*
   *  Update shared_to_box
   *  re-Solicitation is already done - Transfer update shared to box !!!
   *  Transfer boxes in shared maner + asynchronous with dist_comm_graph
   *  Pour avoir l'info if faut echanger le parent gnum -> Allgatherv shared to minimize memory footprint
   *  Transfer of box_extents + gnum + init_location
   */



  /*
   * All pts are redistribute to equilibrate solicitation
   *   We can now rebuild a finer tree to finalize solicitation
   */

  /*
   * Step 3 : Build local tree
   */
  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE) {

    assert(doct->coarse_octree == NULL);
    doct->local_octree = PDM_octree_seq_create(1, // n_point_cloud
                                          doct->local_depth_max,
                                          doct->local_points_in_leaf_max,
                                          doct->local_tolerance);

    PDM_octree_seq_point_cloud_set(doct->local_octree,
                                   0,
                                   equi_n_pts_tot,
                                   equi_pts_coords);

    PDM_octree_seq_build(doct->local_octree);


  } else if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE){

    assert(doct->local_kdtree == NULL);
    doct->local_kdtree = PDM_kdtree_seq_create(1, // n_point_cloud
                                               doct->local_depth_max,
                                               doct->local_points_in_leaf_max,
                                               doct->local_tolerance);

    PDM_kdtree_seq_point_cloud_set(doct->local_kdtree,
                                   0,
                                   equi_n_pts_tot,
                                   equi_pts_coords);

    PDM_kdtree_seq_build(doct->local_kdtree);
  }





  free(equi_pts_coords);
  free(equi_n_pts);

   /*
    * Le nouveau tri donne le lien old_to_new_rank (car on permutera les bbox a peu de choses près)
    * Une fois qu'on connait le tri, on peut faire l'échange en asynchrone pdt que l'arbre se construit ?
    *
    */

  /*
   *  On pourrait faire du hilbert sur des niveaux de feuilles avec gnum = node_id implicitement odered
   *  Avec du sampling
   *  On encode en hilbert le millieu des feuilles -> meme en kd-tree ca marchera
   *  On fait la pré-soliciation pour avoir des poids par boîtes
   *  On fait part_to_block sur des child_id implictement hilbert puis on echange le contenu des noeuds
   *   La stride = le nombre de points -> Ca fait des echanges mais osef !!!
   *   L'algo ressemble ENORMEMENT à PDM_part_assemble_partitions
   *   La pré-solicitation nous permet également d'avoir le lien grossier
   *   Si on le preserve on n'a plus a interoger l'octree !!!
   *   Pdt le transfert de la pré-solicitation --> On construit l'arbre fin
   *   On peut également utiliser l'info du g_child_id dans lequel on est solicité pour preconditionné la recherche local (on gagnera 3/4 niveaux)
   *   A affiner avec le double niveau node / numa
   *   Reprendre le Allgatherv du para_octree sur le comm circulaire
   */

  /*
   * Si l'octree ne change pas l'ordre des points -> On fait une allocation shared
   */

  /*
   *  L'octree seq fait un tri indirect sur les pts, c'est maybe possible de le faire shared inplace !
   *  Sinon on fait un chapeau pour le faire shared !!!!!
   */


  PDM_MPI_Comm_free(&doct->comm_dist_graph);
  PDM_MPI_Comm_free(&doct->comm_shared);
  free(doct->neighbor_in);

  /*
   * A revoir quand on fera l'équilibrage
   */
  // doct->local_octree = coarse_octree;
  // doct->local_kdtree = coarse_kdtree;
}

void
PDM_doctree_point_set
(
 PDM_doctree_t     *doct,
 const int          i_part_cloud,
 const int          n_points,
 const int         *pts_init_location,
 const PDM_g_num_t *pts_g_num,
 const double      *pts_coords
)
{
  assert(i_part_cloud < doct->n_part_cloud);

  doct->n_point_cloud    [i_part_cloud] = n_points;
  doct->pts_g_num        [i_part_cloud] = (PDM_g_num_t *) pts_g_num;
  doct->pts_coords       [i_part_cloud] = (double      *) pts_coords;
  doct->pts_init_location[i_part_cloud] = (int         *) pts_init_location;
}

void
PDM_doctree_solicitation_set
(
 PDM_doctree_t             *doct,
 PDM_tree_solicitation_t    solicitation_kind,
 int                        n_part,
 int                       *n_entity,
 int                      **init_location_entity,
 PDM_g_num_t              **entity_gnum,
 double                   **entity_coords
)
{
  doct->solicitation_kind    = solicitation_kind;
  doct->n_part               = n_part;
  doct->n_entity             = n_entity;
  doct->init_location_entity = init_location_entity;
  doct->entity_gnum          = entity_gnum;
  doct->entity_coords        = entity_coords;
}

// void
// PDM_doctree_results_get
// (
//  PDM_doctree_t      *doct,
//  int                *init_location_entity,
//  int                *dentity1_entity2_idx,
//  PDM_g_num_t        *dentity1_entity2,
//  int                *
//  double             *entity_coords
// )
// {

// }


void
PDM_doctree_free
(
  PDM_doctree_t   *doct
)
{
  free(doct->n_point_cloud    );
  free(doct->pts_g_num        );
  free(doct->pts_coords       );
  free(doct->pts_init_location);

  if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_OCTREE) {
    if(doct->local_octree != NULL) {
      PDM_octree_seq_free(doct->local_octree);
    }
  } else if(doct->local_tree_kind == PDM_DOCTREE_LOCAL_TREE_KDTREE){
    if(doct->local_kdtree != NULL) {
      PDM_kdtree_seq_free(doct->local_kdtree);
    }
  }

  free(doct);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
