#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#include "pdm.h"
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_octree.h"
#include "pdm_logging.h"
#include "pdm_distrib.h"
#include "pdm_vtk.h"
#include "pdm_doctree.h"
#include "pdm_box_gen.h"
#include "pdm_sphere_surf_gen.h"
#include "pdm_part_to_block.h"
#include "pdm_point_tree_seq.h"
#include "pdm_box.h"
#include "pdm_box_tree.h"
#include "pdm_block_to_part.h"
#include "pdm_box_priv.h"

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */

static void
_usage
(
int exit_code
)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -n      <level>  Number of points (default : 10).\n\n"
     "  -radius <level>  Radius of domain (default : 10).\n\n"
     "  -local           Number of points is local (default : global).\n\n"
     "  -rand            Random definition of point coordinates (default : false).\n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nPts   Number of points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *nPts,
 double        *radius,
 int           *local,
 int           *rand
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _nPts = atol(argv[i]);
        *nPts = (PDM_g_num_t) _nPts;
      }
    }

    else if (strcmp(argv[i], "-radius") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *radius = atof(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-local") == 0) {
      *local = 1;
    }

    else if (strcmp(argv[i], "-rand") == 0) {
      *rand = 1;
    }

    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}

static
void
_adaptative_tree
(
  int           n_pts,
  double       *pts_coord,
  PDM_g_num_t  *pts_gnum,
  int           n_box,
  double       *box_extents,
  PDM_g_num_t  *box_gnum,
  PDM_MPI_Comm  comm
)
{
  int i_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (comm, &n_rank);
  /*
   * Hilbert of points AND boxes
   */
  int    *weight_pts = malloc(    n_pts * sizeof(int   ));
  int    *weight_box = malloc(    n_box * sizeof(int   ));
  double *box_center = malloc(3 * n_box * sizeof(double));
  for(int i = 0; i < n_pts; ++i) {
    weight_pts[i] = 1;
  }
  for(int i = 0; i < n_box; ++i) {
    weight_box[i] = 1;
    box_center[3*i  ] = 0.5 * (box_extents[6*i  ] + box_extents[6*i+3]);
    box_center[3*i+1] = 0.5 * (box_extents[6*i+1] + box_extents[6*i+4]);
    box_center[3*i+2] = 0.5 * (box_extents[6*i+2] + box_extents[6*i+5]);
  }

  PDM_part_to_block_t* ptb_pts = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                               PDM_PART_TO_BLOCK_POST_CLEANUP, // A voir avec merge mais attention au init_location
                                                               1.,
                                                               PDM_PART_GEOM_HILBERT,
                                                               &pts_coord,
                                                               &pts_gnum,
                                                               &weight_pts,
                                                               &n_pts,
                                                               1,
                                                               comm);
  free(weight_pts);

  PDM_part_to_block_t* ptb_box = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                               PDM_PART_TO_BLOCK_POST_CLEANUP, // A voir avec merge mais attention au init_location
                                                               1.,
                                                               PDM_PART_GEOM_HILBERT,
                                                               &box_center,
                                                               &box_gnum,
                                                               &weight_box,
                                                               &n_box,
                                                               1,
                                                               comm);
  free(weight_box);
  free(box_center);


  double *blk_pts_coord = NULL;
  PDM_part_to_block_exch(ptb_pts,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
               (void **) &pts_coord,
                         NULL,
               (void **) &blk_pts_coord);


  double *blk_box_extents = NULL;
  PDM_part_to_block_exch(ptb_box,
                         6 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
               (void **) &box_extents,
                         NULL,
               (void **) &blk_box_extents);

  int dn_pts = PDM_part_to_block_n_elt_block_get(ptb_pts);
  int dn_box = PDM_part_to_block_n_elt_block_get(ptb_box);

  /* Global extents */
  double g_global_extents[6];
  double local_extents[6];

  for (int i = 0; i < 3; i++) {
    local_extents[i]     =  DBL_MAX;
    local_extents[i + 3] = -DBL_MAX;
  }

  for (int i = 0; i < dn_box; i++) {
    for (int j = 0; j < 3; j++) {
      local_extents[j]     = PDM_MIN(local_extents[j      ], box_extents[i*3*2 + j      ]);
      local_extents[j+3] = PDM_MAX(local_extents[j + 3], box_extents[i*3*2 + j + 3]);
    }
  }

  for (int i = 0; i < dn_pts; i++) {
    for (int  j = 0; j < 3; j++) {
      if (blk_pts_coord[i*3 + j] < local_extents[j]) {
        local_extents[j] = blk_pts_coord[i*3 + j];
      }
      if (blk_pts_coord[i*3 + j] > local_extents[j + 3]) {
        local_extents[j + 3] = blk_pts_coord[i*3 + j];
      }
    }
  }

  PDM_MPI_Allreduce(local_extents  , g_global_extents    , 3, PDM_MPI_DOUBLE, PDM_MPI_MIN, comm);
  PDM_MPI_Allreduce(local_extents+3, g_global_extents + 3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, comm);

  double s[3];
  double d[3];
  for (int j = 0; j < 3; j++) {
    s[j] = g_global_extents[j];
    d[j] = g_global_extents[j+3] - g_global_extents[j];
  }

  if(1 == 1) {
    char filename[999];
    sprintf(filename, "g_global_extents.vtk");
    PDM_g_num_t one = 1;
    PDM_vtk_write_boxes(filename,
                        1,
                        g_global_extents,
                        &one);
  }


  // Plus besoin du gnum pour l'instant ....
  PDM_part_to_block_free(ptb_pts);
  PDM_part_to_block_free(ptb_box);

  PDM_MPI_Comm comm_alone;
  PDM_MPI_Comm_split(comm, i_rank, 0, &(comm_alone));

  PDM_MPI_Datatype mpi_extent_type;
  PDM_MPI_Type_create_contiguous(6, PDM_MPI_DOUBLE, &mpi_extent_type);
  PDM_MPI_Type_commit(&mpi_extent_type);

  /*
   * Iterative algorithm
   */
  int n_iter = 3;
  for(int i_iter = 0; i_iter < n_iter; ++i_iter) {

    // _adaptative_tree_intersect(dn_pts,
    //                            blk_pts_coord,
    //                            dn_box,
    //                            blk_box_extents,
    //                            comm);

    /*
     * Create octree of point
     */
    PDM_point_tree_seq_t* coarse_tree_pts = PDM_point_tree_seq_create(PDM_DOCTREE_LOCAL_TREE_OCTREE,
                                                                      1, // depth_max
                                                                      1,
                                                                      1e-8);
    PDM_point_tree_seq_point_cloud_set(coarse_tree_pts, dn_pts, blk_pts_coord);
    PDM_point_tree_seq_build(coarse_tree_pts);
    if(1 == 1) {
      char filename[999];
      sprintf(filename, "out_coarse_tree_%i.vtk", i_rank);
      PDM_point_tree_seq_write_nodes(coarse_tree_pts, filename);
    }

    /*
     * Creation of box tree
     */
    PDM_g_num_t *blk_box_gnum      = malloc(dn_box * sizeof(PDM_g_num_t));
    int         *init_location_box = malloc(3 * dn_box * sizeof(int));
    for(int i = 0; i < dn_box; ++i) {
      blk_box_gnum[i] = i + 1;
      init_location_box[3*i  ] = 0;
      init_location_box[3*i+1] = 0;
      init_location_box[3*i+2] = i;
    }

    PDM_box_set_t  *box_set = PDM_box_set_create(3,
                                                 0,  // No normalization to preserve initial extents
                                                 0,  // No projection to preserve initial extents
                                                 dn_box,
                                                 blk_box_gnum,
                                                 blk_box_extents,
                                                 1,
                                                 &dn_box,
                                                 init_location_box,
                                                 comm_alone);
    memcpy (box_set->d, d, sizeof(double) * 3);
    memcpy (box_set->s, s, sizeof(double) * 3);

    int   max_boxes_leaf_shared = 10; // Max number of boxes in a leaf for coarse shared BBTree
    int   max_tree_depth_shared = 1;  // Max tree depth for coarse shared BBTree
    float max_box_ratio_shared  = 5;  // Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes)
    PDM_box_tree_t* bt_shared = PDM_box_tree_create (max_tree_depth_shared,
                                                     max_boxes_leaf_shared,
                                                     max_box_ratio_shared);

    PDM_box_tree_set_boxes (bt_shared,
                            box_set,
                            PDM_BOX_TREE_ASYNC_LEVEL);

    free(blk_box_gnum);
    free(init_location_box);

    if(1 == 1) {
      char filename[999];
      sprintf(filename, "bt_shared_%i.vtk", i_rank);
      PDM_box_tree_write_vtk(filename, bt_shared, -1, 0);
    }


    /*
     * On a tout les arbres locaux, maintenant on doit interoger les arbres en parallèle
     */


    /*
     * Extract extents on all local_tree
     */
    int     n_coarse_pts_box       = 0;
    int    *coarse_pts_box_id      = NULL;
    double *coarse_pts_box_extents = NULL;
    int    *coarse_pts_box_n_pts   = NULL; // Number of point in boxes

    PDM_point_tree_seq_extract_nodes(coarse_tree_pts,
                                     0,
                                     1, // Depth
                                     &n_coarse_pts_box,
                                     &coarse_pts_box_id,
                                     &coarse_pts_box_extents,
                                     &coarse_pts_box_n_pts);

    int *n_g_coarse_pts_box = malloc(n_rank * sizeof(int));
    PDM_MPI_Allgather (&n_coarse_pts_box , 1, PDM_MPI_INT,
                       n_g_coarse_pts_box, 1, PDM_MPI_INT, comm);

    int *g_extract_boxes_idx = (int *) malloc (sizeof(int) * (n_rank+1));
    g_extract_boxes_idx[0] = 0;
    for(int i = 0; i < n_rank; ++i) {
      g_extract_boxes_idx[i+1] = g_extract_boxes_idx[i] + n_g_coarse_pts_box[i];
    }
    double *g_coarse_pts_box_extents = malloc(6 * g_extract_boxes_idx[n_rank] * sizeof(double));

    PDM_MPI_Allgatherv(coarse_pts_box_extents  , n_coarse_pts_box, mpi_extent_type,
                       g_coarse_pts_box_extents, n_g_coarse_pts_box,
                       g_extract_boxes_idx,
                       mpi_extent_type, comm);

    int *g_coarse_pts_box_id = malloc( g_extract_boxes_idx[n_rank] * sizeof(int));
    PDM_MPI_Allgatherv(coarse_pts_box_id  , n_coarse_pts_box, PDM_MPI_INT,
                       g_coarse_pts_box_id, n_g_coarse_pts_box,
                       g_extract_boxes_idx,
                       PDM_MPI_INT, comm);



    /*
     * Build tree
     */
    PDM_g_num_t *coarse_pts_box_gnum         = malloc(    g_extract_boxes_idx[n_rank] * sizeof(PDM_g_num_t));
    int         *init_location_coase_pts_box = malloc(3 * g_extract_boxes_idx[n_rank] * sizeof(int        ));
    for(int i = 0; i < g_extract_boxes_idx[n_rank]; ++i) {
      // coarse_pts_box_gnum[i] = g_coarse_pts_box_id[i] + 1;
      coarse_pts_box_gnum[i] = g_coarse_pts_box_id[i]; // On suppose que root = 0 donc g_id lineraire a partir de 1
      init_location_coase_pts_box[3*i  ] = 0;
      init_location_coase_pts_box[3*i+1] = 0;
      init_location_coase_pts_box[3*i+2] = i;
    }
    free(g_coarse_pts_box_id);

    PDM_box_set_t  *coarse_pts_box_set = PDM_box_set_create(3,
                                                            0,  // No normalization to preserve initial extents
                                                            0,  // No projection to preserve initial extents
                                                            g_extract_boxes_idx[n_rank],
                                                            coarse_pts_box_gnum,
                                                            g_coarse_pts_box_extents,
                                                            1,
                                                            &g_extract_boxes_idx[n_rank],
                                                            init_location_coase_pts_box,
                                                            comm_alone);
    memcpy (coarse_pts_box_set->d, d, sizeof(double) * 3);
    memcpy (coarse_pts_box_set->s, s, sizeof(double) * 3);

    int   max_boxes_leaf_coarse = 1;   // Max number of boxes in a leaf for coarse coarse BBTree
    int   max_tree_depth_coarse = 31;  // Max tree depth for coarse coarse BBTree
    float max_box_ratio_coarse  = 5;   // Max ratio for local BBTree (nConnectedBoxe < ratio * nBoxes)
    PDM_box_tree_t* coarse_pts_bt_shared = PDM_box_tree_create (max_tree_depth_coarse,
                                                                max_boxes_leaf_coarse,
                                                                max_box_ratio_coarse);

    PDM_box_tree_set_boxes (coarse_pts_bt_shared,
                            coarse_pts_box_set,
                            PDM_BOX_TREE_ASYNC_LEVEL);

    free(coarse_pts_box_gnum);
    free(init_location_coase_pts_box);

    if(1 == 1) {
      char filename[999];
      sprintf(filename, "coarse_pts_box_set_%i.vtk", i_rank);
      PDM_box_tree_write_vtk(filename, coarse_pts_bt_shared, -1, 0);

      sprintf(filename, "coarse_pts_box_%i.vtk", i_rank);
      PDM_vtk_write_boxes(filename,
                          g_extract_boxes_idx[n_rank],
                          g_coarse_pts_box_extents,
                          NULL);
    }


    /*
     * Extract coarse box for box_tree
     */
    int n_coarse_box_box = 0;
    double *coarse_box_extents = NULL;
    int  n_extract_child = 0;
    int *extract_child_id = NULL;
    PDM_box_tree_extract_extents(bt_shared,
                                 0,
                                 1,
                                 &n_coarse_box_box,
                                 &coarse_box_extents,
                                 &n_extract_child,
                                 &extract_child_id);
    log_trace("n_coarse_box_box = %i \n", n_coarse_box_box);
    log_trace("n_extract_child  = %i \n", n_extract_child);
    PDM_log_trace_array_int(extract_child_id, n_extract_child, "extract_child_id ::");
    free(extract_child_id);

    if(1 == 1) {
      char filename[999];
      sprintf(filename, "coarse_box_box_%i.vtk", i_rank);
      PDM_vtk_write_boxes(filename,
                          n_coarse_box_box,
                          coarse_box_extents,
                          NULL);
    }

    /*
     * Maintenant on veut interoger l'arbre global des points avec la solicitation des boites
     *  On peut avoir la connectivité : coarse_pts_box_to_coarse_box
     *  On devra dans tout les cas échangé des choses
     */
    int *coarse_box_to_coarse_box_pts_idx = NULL;
    int *coarse_box_to_coarse_box_pts = NULL;
    PDM_box_tree_intersect_boxes_boxes(coarse_pts_bt_shared,
                                       -1,
                                       n_coarse_box_box,
                                       coarse_box_extents,
                                       &coarse_box_to_coarse_box_pts_idx,
                                       &coarse_box_to_coarse_box_pts);

    free(coarse_box_extents);

    PDM_log_trace_connectivity_int(coarse_box_to_coarse_box_pts_idx, coarse_box_to_coarse_box_pts, n_coarse_box_box, "coarse_box_to_coarse_box_pts ::");

    const int         *box_origin   = PDM_box_set_origin_get(coarse_pts_box_set);
    const PDM_g_num_t *box_pts_gnum = PDM_box_set_get_g_num (coarse_pts_box_set);
    PDM_log_trace_array_int(box_origin, 3 * g_extract_boxes_idx[n_rank], "box_origin") ;

    /*
     * Equilibrate coarse_box_pts -> map on box_gnum_id (of pts )
     *  Puis chaque connexion correspond à une feuille de l'arbre de boites,
     *  Qu'on envoie au futur boites de pts
     *   Donc il faut créer une partition de boites à partir du box_tree -> On veut tout les boxes_id par feuilles en gros
     *    On est malin -> Si aucune boites intersecte l'arbre de pts on met stride = 0
     *  Il faut également update le coarse_box_to_coarse_box_pts (avec le nouveau numero de feuilles et l'envoyé)
     *
     *
     *  A la reception des boites on connait donc deja leur lien avec l'octree de pts grossier
     */
    int n_connect_box_coarse_box_pts = coarse_box_to_coarse_box_pts_idx[n_coarse_box_box];
    PDM_g_num_t  *coarse_box_to_coarse_box_pts_gnum = malloc(n_connect_box_coarse_box_pts * sizeof(PDM_g_num_t));
    double       *weight                            = malloc(n_connect_box_coarse_box_pts * sizeof(double     ));

    for(int i = 0; i < n_connect_box_coarse_box_pts; ++i) {
      coarse_box_to_coarse_box_pts_gnum[i] = box_pts_gnum[coarse_box_to_coarse_box_pts[i]]; // + 1;
      weight[i] = 1.;
    }


    PDM_part_to_block_t* ptb_equi_pts_box = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                                     PDM_PART_TO_BLOCK_POST_MERGE,
                                                                     1.,
                                                                     &coarse_box_to_coarse_box_pts_gnum,
                                                                     &weight,
                                                                     &n_connect_box_coarse_box_pts,
                                                                     1,
                                                                     comm);



    free(coarse_box_to_coarse_box_pts_gnum);
    free(weight);

    int n_next_coarse_box_pts = PDM_part_to_block_n_elt_block_get(ptb_equi_pts_box);
    PDM_g_num_t* next_gnum_box_pts = PDM_part_to_block_block_gnum_get(ptb_equi_pts_box);

    PDM_log_trace_array_long(next_gnum_box_pts, n_next_coarse_box_pts, "next_gnum_box_pts :");

    /*
     * Il faut maintenant recupérer les pts des boites de pts !
     *   block_to_part avec gnum = block_to_part_get_gnum()
     */
    PDM_g_num_t* distrib_coarse_box_pts = PDM_compute_entity_distribution(comm, n_coarse_pts_box);
    PDM_block_to_part_t* btp = PDM_block_to_part_create(distrib_coarse_box_pts,
                                (const PDM_g_num_t **)  &next_gnum_box_pts,
                                                        &n_next_coarse_box_pts,
                                                        1,
                                                        comm);

    /*
     * Extraction depuis l'octree et envoi
     */
    double* tree_coords = NULL;
    PDM_point_tree_seq_sorted_points_get(coarse_tree_pts, &tree_coords);

    /*
     * Copie pour l'instant
     */
    int point_range[2];
    int *n_pts_in_leaf = malloc(n_coarse_pts_box * sizeof(int));

    // PDM_sort_int(coarse_pts_box_id, NULL, n_coarse_pts_box);
    PDM_log_trace_array_int(coarse_pts_box_id, n_coarse_pts_box, "coarse_pts_box_id ::");
    assert(n_coarse_pts_box == 8);
    for(int i = 0; i < n_coarse_pts_box; ++i) {
      n_pts_in_leaf[coarse_pts_box_id[i]-1] = PDM_point_tree_seq_point_range_get(coarse_tree_pts, coarse_pts_box_id[i], point_range);
      log_trace(" coarse_pts_box_id[%i] = %i | range = %i / %i \n", i, coarse_pts_box_id[i], point_range[0], point_range[1]);
    }

    // double *reorder_tree_coord = malloc(3 * dn_pts * sizeof(double));


    if(1 == 1) {
      char filename[999];
      sprintf(filename, "tree_coords_%i.vtk", i_rank);
      PDM_vtk_write_point_cloud(filename,
                                dn_pts,
                                tree_coords,
                                NULL,
                                NULL);
    }

    int    **tmp_next_pts_coords_n = NULL;
    double **tmp_next_pts_coords   = NULL;
    PDM_block_to_part_exch(btp,
                           3 * sizeof(double),
                           PDM_STRIDE_VAR_INTERLACED,
                           n_pts_in_leaf,
                           tree_coords,
                           &tmp_next_pts_coords_n,
                (void ***) &tmp_next_pts_coords);
    int    *next_pts_coords_n = tmp_next_pts_coords_n[0];
    double *next_pts_coords   = tmp_next_pts_coords  [0];
    free(tmp_next_pts_coords_n);
    free(tmp_next_pts_coords  );

    dn_pts = 0;
    for(int i = 0; i < n_next_coarse_box_pts; ++i) {
      dn_pts += next_pts_coords_n[i];
    }

    free(next_pts_coords_n);
    free(n_pts_in_leaf);

    free(blk_pts_coord);
    blk_pts_coord = next_pts_coords;

    if(1 == 1) {

      char filename[999];
      sprintf(filename, "blk_pts_coord_%i.vtk", i_rank);
      PDM_vtk_write_point_cloud(filename,
                                dn_pts,
                                blk_pts_coord,
                                NULL,
                                NULL);
    }



    // exit(1);

    PDM_point_tree_seq_free(coarse_tree_pts);

    /*
     * Envoie du lien implicite entre les boites et les pts --> Sous communicateur par feuilles ?
     *  --> Construire un neighbor --> plus smart :) !!!
     */
    free(distrib_coarse_box_pts);
    PDM_block_to_part_free(btp);

    PDM_part_to_block_free(ptb_equi_pts_box);


    PDM_box_set_destroy (&coarse_pts_box_set);
    PDM_box_tree_destroy(&coarse_pts_bt_shared);

    free(g_extract_boxes_idx);
    free(coarse_pts_box_id     );
    free(coarse_pts_box_extents);
    free(coarse_pts_box_n_pts  );
    free(n_g_coarse_pts_box);
    free(g_coarse_pts_box_extents);
    free(coarse_box_to_coarse_box_pts_idx);
    free(coarse_box_to_coarse_box_pts);
    PDM_box_set_destroy (&box_set);
    PDM_box_tree_destroy(&bt_shared);

  }

  PDM_MPI_Comm_free(&comm_alone);
  PDM_MPI_Type_free(&mpi_extent_type);

  free(blk_pts_coord);
  free(blk_box_extents);

}


/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */

int
main
(
int argc,
char *argv[]
)
{

  PDM_MPI_Init (&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int n_rank;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &n_rank);

  PDM_g_num_t nPts   = 10;
  double radius = 0.5;
  int local = 0;
  int rand = 0;

  _read_args(argc,
             argv,
             &nPts,
             &radius,
             &local,
             &rand);

  /* Initialize random */

  if (rand) {
    srand(time(NULL));
  }
  else {
    srand(i_rank);
  }

  /* Random point cloud */
  /* Generate src and tgt point clouds */
  int          n_src     = 0;
  double      *src_coord = NULL;

  const double x_center = 0.;
  const double y_center = 0.15;
  const double z_center = 0.85;

  int         *dback_face_vtx_idx = NULL;
  PDM_g_num_t *dback_face_vtx     = NULL;
  PDM_g_num_t *back_distrib_vtx   = NULL;
  PDM_g_num_t *back_distrib_face  = NULL;
  PDM_sphere_surf_icosphere_gen(comm,
                                nPts,
                                x_center,
                                y_center,
                                z_center,
                                radius,
                                &src_coord,
                                &dback_face_vtx_idx,
                                &dback_face_vtx,
                                &back_distrib_vtx,
                                &back_distrib_face);

  n_src = back_distrib_vtx[i_rank+1] - back_distrib_vtx[i_rank];
  PDM_g_num_t *src_g_num = malloc(n_src * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_src; ++i) {
    src_g_num[i] = back_distrib_vtx[i_rank] + i + 1;
  }

  free(dback_face_vtx_idx);
  free(dback_face_vtx    );
  free(back_distrib_vtx  );
  free(back_distrib_face );

  if(1 == 1) {
    char filename[999];
    sprintf(filename, "sphere_cloud_%i.vtk", i_rank);
    PDM_vtk_write_point_cloud(filename,
                              n_src,
                              src_coord,
                              src_g_num,
                              NULL);
  }




  PDM_doctree_local_tree_t local_tree_kind = PDM_DOCTREE_LOCAL_TREE_OCTREE;
  // PDM_doctree_local_tree_t local_tree_kind = PDM_DOCTREE_LOCAL_TREE_KDTREE;
  // PDM_doctree_t *doct = PDM_doctree_create(comm,
  //                                          3,
  //                                          1,
  //                                          NULL, // global_extents
  //                                          local_tree_kind);

  int *init_location_pts = malloc(3 * n_src * sizeof(int));
  for(int i = 0; i < n_src; ++i) {
    init_location_pts[3*i  ] = i_rank;
    init_location_pts[3*i+1] = 0; // i_part
    init_location_pts[3*i+2] = i;
  }

  // PDM_doctree_point_set(doct,
  //                       0,
  //                       n_src,
  //                       init_location_pts,
  //                       src_g_num,
  //                       src_coord);
  int n_box   = 0;
  int n_vtx_x = 12;
  int n_vtx_y = 12;
  int n_vtx_z = 12;
  double      *box_extents = NULL;
  PDM_g_num_t *box_gnum    = NULL;
  PDM_box_gen_cartesian(comm,
                        n_vtx_x,
                        n_vtx_y,
                        n_vtx_z,
                        -0., -0., -0.,
                        2., 2., 2.,
                        &n_box,
                        &box_extents,
                        &box_gnum);

  if(1 == 1) {
    char filename[999];
    sprintf(filename, "boxes_%i.vtk", i_rank);
    PDM_vtk_write_boxes(filename,
                        n_box,
                        box_extents,
                        box_gnum);
  }

  int *init_location_box = malloc(3 * n_box * sizeof(int));
  for(int i = 0; i < n_box; ++i) {
    init_location_box[3*i  ] = i_rank;
    init_location_box[3*i+1] = 0; // i_part
    init_location_box[3*i+2] = i;
  }

  // PDM_doctree_solicitation_set(doct,
  //                              PDM_TREE_SOLICITATION_BOXES_POINTS,
  //                              1,
  //                              &n_box,
  //                              &init_location_box,
  //                              &box_gnum,
  //                              &box_extents);

  // PDM_doctree_build(doct);

  _adaptative_tree(n_src,
                   src_coord,
                   src_g_num,
                   n_box,
                   box_extents,
                   box_gnum,
                   comm);

  // if(1 == 1) {

  //   int         *box_pts_idx = NULL;
  //   PDM_g_num_t *box_pts     = NULL;
  //   double      *pts_coord   = NULL;
  //   PDM_doctree_results_in_orig_frame_get(doct,
  //                                         n_box,
  //                                         box_gnum,
  //                                         &box_pts_idx,
  //                                         &box_pts,
  //                                         &pts_coord);

  //   free(box_pts_idx);
  //   free(box_pts    );
  //   free(pts_coord  );
  // }


  // PDM_doctree_free(doct);


  free(box_gnum);
  free(box_extents);
  free(init_location_box);
  free(init_location_pts);

  free (src_coord);
  free (src_g_num);

  if (i_rank == 0) {
    PDM_printf ("-- End\n");
    fflush(stdout);
  }


  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;

}
