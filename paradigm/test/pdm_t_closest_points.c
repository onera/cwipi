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
 double        *radius,
 int           *local,
 int           *rand,
 int           *clumps
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

    else if (strcmp(argv[i], "-local") == 0) {
      *local = 1;
    }

    else if (strcmp(argv[i], "-rand") == 0) {
      *rand = 1;
    }

    else if (strcmp(argv[i], "-clumps") == 0) {
      *clumps = 1;
    }


    else {
      _usage(EXIT_FAILURE);
    }
    i++;
  }
}

/**
 *
 * \brief  Random value
 *
 *  \return a random double in [-1, 1]
 */

static double
_random01
(
 void
)
{
  int sign;
  int rsigna = rand();
  int rsignb = rand();
  sign = (rsigna - rsignb) / PDM_ABS (rsigna - rsignb);
  double resultat = sign*((double)rand())/((double)RAND_MAX);
  return resultat;
}


static void
_gen_clouds_random2
(
 const int      nSrc,
 const int      nTgt_l,
 const double   radius,
 const int      numProcs,
 const int      i_rank,
 double       **src_coord,
 double       **tgt_coord,
 int           *nSrc_l
 )
{
  *nSrc_l = (int) (nSrc/numProcs);
  if (i_rank < nSrc%numProcs) {
    (*nSrc_l) += 1;
  }

  *src_coord = malloc (sizeof(double) * 3 * (*nSrc_l));
  double *_src_coord = *src_coord;
  double x;
  int idx = 0;
  for (int i = 0; i < numProcs*(*nSrc_l); i++) {
    for (int j = 0; j < 3; j++) {
      x = _random01() * radius;
      if (i%numProcs == i_rank) {
        _src_coord[idx++] = x;
      }
    }
  }


  *tgt_coord = malloc (sizeof(double) * 3 * nTgt_l);
  double *_tgt_coord = *tgt_coord;
  for (int i = 0; i < nTgt_l; i++) {
    for (int j = 0; j < 3; j++) {
      _tgt_coord[3*i+j] = _random01() * radius;
    }
  }
}


/* FIXME: _gen_clouds_random a conserver ? */
   /* static void */
/* _gen_clouds_random */
/* ( */
/*  const int      nSrc_l, */
/*  const int      nTgt_l, */
/*  const double   radius, */
/*  double       **src_coord, */
/*  double       **tgt_coord */
/*  ) */
/* { */

/*   *src_coord = malloc (sizeof(double) * 3 * nSrc_l); */
/*   double *_src_coord = *src_coord; */
/*   for (int i = 0; i < nSrc_l; i++) { */
/*     for (int j = 0; j < 3; j++) { */
/*       _src_coord[3*i+j] = _random01() * radius; */
/*     } */
/*   } */

/*   *tgt_coord = malloc (sizeof(double) * 3 * nTgt_l); */
/*   double *_tgt_coord = *tgt_coord; */
/*   for (int i = 0; i < nTgt_l; i++) { */
/*     for (int j = 0; j < 3; j++) { */
/*       _tgt_coord[3*i+j] = _random01() * radius; */
/*     } */
/*   } */
/* } */











static void
_gen_clouds_clumps
(
 const int      n_closest_points,
 const int      nTgt,
 const double   radius,
 const int      numProcs,
 const int      i_rank,
 const double   clump_scale,
 double       **src_coord,
 double       **tgt_coord,
 int           *nTgt_l,
 int           *nSrc_l
 )
{
  *nTgt_l = (int) nTgt/numProcs;
  *nSrc_l = n_closest_points * (*nTgt_l);

  const double clump_radius = clump_scale * radius / pow(nTgt, 1./3);

  *tgt_coord = malloc (sizeof(double) * 3 * (*nTgt_l));
  *src_coord = malloc (sizeof(double) * 3 * (*nSrc_l));

  double *_tgt_coord = *tgt_coord;
  double *_src_coord = *src_coord;
  double x;

  int ii = 0;
  int idx = 0;
  for (int i = 0; i < (*nTgt_l)*numProcs; i++) {
    for (int j = 0; j < 3; j++) {
      x = _random01() * radius;
      if (i%numProcs == i_rank) {
        _tgt_coord[3*ii+j] = x;
      }
    }

    for (int k = 0; k < n_closest_points; k++) {
      for (int j = 0; j < 3; j++) {
        x = _random01() * clump_radius;

        if (i%numProcs == i_rank)
          _src_coord[idx++] = _tgt_coord[3*ii+j] + x;
      }
    }

    if (i%numProcs == i_rank)
      ii++;
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

  int i_rank;
  PDM_MPI_Comm_rank (PDM_MPI_COMM_WORLD, &i_rank);

  int numProcs;
  PDM_MPI_Comm_size (PDM_MPI_COMM_WORLD, &numProcs);


  char *version = PDM_version_get();

  printf("Version de ParaDiGM : %s\n", version);
  free(version);

  int n_closest_points = 10;
  PDM_g_num_t nSrc = 10;
  PDM_g_num_t nTgt = 10;
  double radius = 10.;
  int local = 0;
  int rand = 0;
  int clumps = 0;

  _read_args(argc,
             argv,
             &n_closest_points,
             &nSrc,
             &nTgt,
             &radius,
             &local,
             &rand,
             &clumps);

  /* Initialize random */

  if (rand) {
    srand(time(NULL));
  }
  else {
    srand(i_rank);
  }


  /* Define the numbers of Source/Target points */

  int _n_src_l = 0, _n_tgt_l = 0;
  if (local) {
    _n_src_l = (int) nSrc;
    _n_tgt_l = (int) nTgt;
  }
  else {
    _n_src_l = (int) (nSrc/numProcs);
    _n_tgt_l = (int) (nTgt/numProcs);
    if (i_rank < nSrc%numProcs) {
      _n_src_l += 1;
    }
    if (i_rank < nTgt%numProcs) {
      _n_tgt_l += 1;
    }
  }

  double *src_coords = NULL;
  double *tgt_coords = NULL;


  if (clumps) {
    const double clump_scale = 0.3;

    _gen_clouds_clumps (n_closest_points,
                        nTgt,
                        radius,
                        numProcs,
                        i_rank,
                        clump_scale,
                        &src_coords,
                        &tgt_coords,
                        &_n_tgt_l,
                        &_n_src_l);
  } else {
    /*_gen_clouds_random (_n_src_l,
      _n_tgt_l,
      radius,
      &src_coords,
      &tgt_coords);*/
    _gen_clouds_random2 (nSrc,
                         _n_tgt_l,
                         radius,
                         numProcs,
                         i_rank,
                         &src_coords,
                         &tgt_coords,
                         &_n_src_l);
  }


  //n_closest_points = PDM_MIN (n_closest_points, nSrc);

  /* Source points definition  */
  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD, PDM_OWNERSHIP_USER);

  double *src_char_length = malloc(sizeof(double) * _n_src_l);

  for (int i = 0; i < _n_src_l; i++) {
    src_char_length[i] = radius * 1.e-6;
  }

  PDM_gnum_set_from_coords (gen_gnum, 0, _n_src_l, src_coords, src_char_length);

  PDM_gnum_compute (gen_gnum);

  PDM_g_num_t *src_gnum = PDM_gnum_get(gen_gnum, 0);

  PDM_gnum_free (gen_gnum);



  /* Target points definition */
  gen_gnum = PDM_gnum_create (3, 1, PDM_FALSE, 1e-3, PDM_MPI_COMM_WORLD, PDM_OWNERSHIP_USER);

  double *tgt_char_length = malloc(sizeof(double) * _n_tgt_l);

  for (int i = 0; i < _n_tgt_l; i++) {
    tgt_char_length[i] = radius * 1.e-6;
  }

  PDM_gnum_set_from_coords (gen_gnum, 0, _n_tgt_l, tgt_coords, tgt_char_length);

  PDM_gnum_compute (gen_gnum);

  PDM_g_num_t *tgt_gnum = PDM_gnum_get(gen_gnum, 0);

  PDM_gnum_free (gen_gnum);



  PDM_closest_point_t* clsp = PDM_closest_points_create (PDM_MPI_COMM_WORLD,
                                                         n_closest_points,
                                                         PDM_OWNERSHIP_KEEP);

  PDM_closest_points_n_part_cloud_set (clsp,
                                       1,
                                       1);

  PDM_closest_points_src_cloud_set (clsp,
                                    0,
                                    _n_src_l,
                                    src_coords,
                                    src_gnum);

  PDM_closest_points_tgt_cloud_set (clsp,
                                    0,
                                    _n_tgt_l,
                                    tgt_coords,
                                    tgt_gnum);


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

    for (int i = 0; i < _n_tgt_l; i++) {
      printf("Target point #%d ("PDM_FMT_G_NUM") [%f, %f, %f]\n", i, tgt_gnum[i],
             tgt_coords[3*i], tgt_coords[3*i+1], tgt_coords[3*i+2]);
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

  free (src_coords);
  free (src_char_length);
  free (src_gnum);
  free (tgt_coords);
  free (tgt_char_length);
  free (tgt_gnum);

  if (i_rank == 0) {

    PDM_printf ("\nfin Test\n");

  }

  PDM_MPI_Barrier (PDM_MPI_COMM_WORLD);

  PDM_MPI_Finalize ();

  return 0;
}
