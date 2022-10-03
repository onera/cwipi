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
#include "pdm_version.h"

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
     "  -c       <level> Number of closest points (default : 10).\n\n"
     "  -s       <level> Number of Source points (default : 10).\n\n"
     "  -t       <level> Number of Target points (default : 10).\n\n"
     "  -radius <level>  Radius of domain (default : 10).\n\n"
     "  -local           Number of points is local (default : global).\n\n"
     "  -rand            Random definition of point coordinates (default : false).\n\n"
     "  -clumps          Source points distributed in clumps around target points (default : false).\n\n"
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
 * \param [inout] nSrc   Number of Source points
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
 int           *nClosest,
 PDM_g_num_t   *nSrc,
 PDM_g_num_t   *nTgt,
 double        *radius
 )
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-c") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *nClosest = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-s") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _nSrc = atol(argv[i]);
        *nSrc = (PDM_g_num_t) _nSrc;
      }
    }

    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        long _nTgt = atol(argv[i]);
        *nTgt = (PDM_g_num_t) _nTgt;
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

  if (i_rank == 0) {
    char *version = PDM_version_get();

    printf("Version de ParaDiGM : %s\n", version);
    free(version);
  }

  int         n_closest_points = 10;
  PDM_g_num_t gn_src           = 10;
  PDM_g_num_t gn_tgt           = 10;
  double      radius           = 10.;

  _read_args(argc,
             argv,
             &n_closest_points,
             &gn_src,
             &gn_tgt,
             &radius);

  if (i_rank == 0) {
    PDM_printf ("%Parametres : \n");
    PDM_printf ("  - n_rank           : %d\n", n_rank);
    PDM_printf ("  - n_closest_points : %d\n", n_closest_points);
    PDM_printf ("  - n_src            : "PDM_FMT_G_NUM"\n", gn_src);
    PDM_printf ("  - n_tgt            : "PDM_FMT_G_NUM"\n", gn_tgt);
    PDM_printf ("  - radius           : %f\n", radius);
  }


  /* Generate src and tgt point clouds */
  int          n_src;
  double      *src_coord;
  PDM_g_num_t *src_g_num;
  PDM_point_cloud_gen_random (comm,
                              0, // seed
                              0, // geometric_g_num
                              gn_src,
                              -radius, -radius, -radius,
                              radius, radius, radius,
                              &n_src,
                              &src_coord,
                              &src_g_num);


  int          n_tgt;
  double      *tgt_coord;
  PDM_g_num_t *tgt_g_num;
  PDM_point_cloud_gen_random (comm,
                              123456789, // seed
                              0,         // geometric_g_num
                              gn_tgt,
                              -radius, -radius, -radius,
                              radius, radius, radius,
                              &n_tgt,
                              &tgt_coord,
                              &tgt_g_num);


  PDM_closest_point_t* clsp = PDM_closest_points_create (PDM_MPI_COMM_WORLD,
                                                         n_closest_points,
                                                         PDM_OWNERSHIP_KEEP);

  PDM_closest_points_n_part_cloud_set (clsp,
                                       1,
                                       1);

  PDM_closest_points_src_cloud_set (clsp,
                                    0,
                                    n_src,
                                    src_coord,
                                    src_g_num);

  PDM_closest_points_tgt_cloud_set (clsp,
                                    0,
                                    n_tgt,
                                    tgt_coord,
                                    tgt_g_num);


  PDM_closest_points_compute (clsp);


  PDM_closest_points_dump_times (clsp);

  PDM_g_num_t *closest_src_gnum = NULL;
  double      *closest_src_dist = NULL;

  PDM_closest_points_get (clsp,
                          0,
                          &closest_src_gnum,
                          &closest_src_dist);


  if (0 == 1) {
    printf("\n\n============================\n\n");

    for (int i = 0; i < n_tgt; i++) {
      printf("Target point #%d ("PDM_FMT_G_NUM") [%f, %f, %f]\n", i, tgt_g_num[i],
             tgt_coord[3*i], tgt_coord[3*i+1], tgt_coord[3*i+2]);
      for (int j = 0; j < n_closest_points; j++)
        printf("\t%d:\t"PDM_FMT_G_NUM"\t%f\n",
               j+1,
               closest_src_gnum[n_closest_points*i + j],
               closest_src_dist[n_closest_points*i + j]);
      printf("\n\n");
    }


    printf("============================\n\n");
  }

  PDM_closest_points_free (clsp);





  /* Free */

  free (src_coord);
  free (src_g_num);
  free (tgt_coord);
  free (tgt_g_num);

  if (i_rank == 0) {

    PDM_printf ("-- End\n");

  }

  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;
}
