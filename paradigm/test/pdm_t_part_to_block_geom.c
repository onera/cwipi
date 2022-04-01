#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_gnum.h"
#include "pdm_distrib.h"
#include "pdm_point_cloud_gen.h"
#include "pdm_closest_points.h"
#include "pdm_part_to_block.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_version.h"
#include "pdm_hilbert.h"

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
     "  -s       <level> Number of Source points (default : 10).\n\n"
     "  -h               This message.\n\n");

  exit (exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]    argc   Number of arguments
 * \param [in]    argv   Arguments
 * \param [inout] nClosest   Number of closest points
 * \param [inout] n_src   Number of Source points
 * \param [inout] nTgt   Number of Target points
 * \param [inout] ls     Low scalability
 * \param [inout] length Length of domains
 *
 */

static void
_read_args
(
 int            argc,
 char         **argv,
 PDM_g_num_t   *n_src
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);
    else if (strcmp(argv[i], "-s") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _n_src = atol(argv[i]);
        *n_src = (PDM_g_num_t) _n_src;
      }
    }
    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
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

  int i_rank, n_rank;
  PDM_MPI_Comm_rank (comm, &i_rank);
  PDM_MPI_Comm_size (comm, &n_rank);

  PDM_g_num_t gn_src = 10;
  double      radius = 1;

  _read_args(argc,
             argv,
             &gn_src);

  if (i_rank == 0) {
    PDM_printf ("%Parametres : \n");
    PDM_printf ("  - n_rank           : %d\n", n_rank);
    PDM_printf ("  - n_src            : "PDM_FMT_G_NUM"\n", gn_src);
  }


  /* Generate src and tgt point clouds */
  int          n_src     = 0;
  double      *src_coord = NULL;
  PDM_g_num_t *src_g_num = NULL;
  PDM_point_cloud_gen_random (comm,
                              gn_src,
                              -radius, -radius, -radius,
                              radius, radius, radius,
                              &n_src,
                              &src_coord,
                              &src_g_num);

  int *weight =  malloc( n_src * sizeof(int));
  for(int i = 0; i < n_src; ++i) {
    weight[i] = 1;
  }

  PDM_part_to_block_t* ptb = PDM_part_to_block_geom_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                           PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                           1.,
                                                           PDM_PART_GEOM_HILBERT,
                                                           &src_coord,
                                                           &src_g_num,
                                                           &weight,
                                                           &n_src,
                                                           1,
                                                           comm);

  double *blk_src_coord = NULL;
  PDM_part_to_block_exch(ptb,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
               (void **) &src_coord,
                         NULL,
               (void **) &blk_src_coord);

  PDM_g_num_t *blk_check_gnum = NULL;
  PDM_part_to_block_exch(ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
               (void **) &src_g_num,
                         NULL,
               (void **) &blk_check_gnum);


  int          n_parent    = PDM_part_to_block_n_elt_block_get  (ptb);
  PDM_g_num_t* parent_gnum = PDM_part_to_block_block_gnum_get   (ptb);
  PDM_g_num_t* distrib_pts = PDM_part_to_block_distrib_index_get(ptb);

  if(0 == 1) {
    PDM_log_trace_array_long(distrib_pts, n_rank+1, "distrib_pts :: ");
    PDM_log_trace_array_long(parent_gnum, n_parent, "parent_gnum :: ");
    PDM_log_trace_array_long(blk_check_gnum, n_parent, "blk_check_gnum :: ");
  }

  if(1 == 1) {
    char filename[999];
    sprintf(filename, "origin_pvtx_coord_%2.2d.vtk", i_rank);

    PDM_vtk_write_point_cloud(filename,
                              n_src,
                              src_coord,
                              src_g_num,
                              NULL);


    PDM_g_num_t *debug_gnum = malloc( n_parent * sizeof(PDM_g_num_t));

    for(int i = 0; i < n_parent; ++i) {
      debug_gnum[i] = distrib_pts[i_rank] + i + 1;
    }

    sprintf(filename, "redistrib_pvtx_coord_%2.2d.vtk", i_rank);
    // PDM_vtk_write_point_cloud(filename,
    //                           n_parent,
    //                           blk_src_coord,
    //                           parent_gnum,
    //                           NULL);

    // To debug hilbert ordering
    PDM_vtk_write_point_cloud(filename,
                              n_parent,
                              blk_src_coord,
                              debug_gnum,
                              NULL);

    free(debug_gnum);
  }

  /*
   * Check ordering
   */
  int dim = 3;
  double extents[2*dim]; /** DIM x 2**/

  /** Get EXTENTS **/
  PDM_hilbert_get_coord_extents_par(dim, n_parent, blk_src_coord, extents, comm);

  PDM_hilbert_code_t* hilbert_codes = (PDM_hilbert_code_t * ) malloc(n_parent * sizeof(PDM_hilbert_code_t));
  PDM_hilbert_encode_coords(dim, PDM_HILBERT_CS, extents, n_parent, blk_src_coord, hilbert_codes);

  PDM_hilbert_code_t first = 0.;
  if(n_parent > 0) {
    first = hilbert_codes[0];
  }
  for(int i = 1; i < n_parent; ++i) {
    assert(first <= hilbert_codes[i]);
    first = hilbert_codes[i];
  }

  // PDM_log_trace_array_double(hilbert_codes, n_parent, "hilbert_codes ::");


  PDM_part_to_block_free(ptb);

  /* Free */
  free (hilbert_codes);
  free (blk_src_coord);
  free (blk_check_gnum);
  free (weight);
  free (src_coord);
  free (src_g_num);

  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);
  if (i_rank == 0) {
    PDM_printf ("-- End\n");
  }


  PDM_MPI_Finalize ();

  return 0;
}
