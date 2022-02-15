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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_part.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_distrib.h"

#include "pdm_dcube_nodal_gen.h"
#include "pdm_domain_interface.h"

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Usage
 *
 */

static void
_usage(int exit_code)
{
  PDM_printf
  ("\n"
   "  Usage: \n\n"
   "  -n      <level>  Number of vertices on the cube side.\n\n"
   "  -l      <level>  Cube length.\n\n"
   "  -n_part <level>  Number of partitions par process.\n\n"
   "  -post            Ensight outputs (only if n_part == 1). \n\n"
   "  -parmetis        Call ParMETIS.\n\n"
   "  -pt-scocth       Call PT-Scotch.\n\n"
   "  -h               This message.\n\n");

  exit(exit_code);
}


/**
 *
 * \brief  Read arguments from the command line
 *
 * \param [in]      argc     Number of arguments
 * \param [in]      argv     Arguments
 * \param [inout]   n_vtx_seg  Number of vertices on the cube side
 * \param [inout]   length   Cube length
 * \param [inout]   n_part   Number of partitions par process
 * \param [inout]   post     Ensight outputs status
 * \param [inout]   method   Partitioner (1 ParMETIS, 2 Pt-Scotch)
 *
 */

static void
_read_args(int            argc,
           char         **argv,
           PDM_g_num_t   *nx,
           PDM_g_num_t   *ny,
           PDM_g_num_t   *nz,
           int           *n_dom_i,
           int           *n_dom_j,
           int           *n_dom_k,
           int           *periodic_i,
           int           *periodic_j,
           int           *periodic_k,
           int           *t_elt,
           double        *length,
           int           *n_part,
           int           *post,
           int           *method)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
        *ny = (PDM_g_num_t) _n;
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nx = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-ny") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *ny = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-nz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        long _n = atol(argv[i]);
        *nz = (PDM_g_num_t) _n;
      }
    }
    else if (strcmp(argv[i], "-ni") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_dom_i = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-nj") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_dom_j = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-nk") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_dom_k = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-pi") == 0) {
      *periodic_i = 1;
    }
    else if (strcmp(argv[i], "-pj") == 0) {
      *periodic_j = 1;
    }
    else if (strcmp(argv[i], "-pk") == 0) {
      *periodic_k = 1;
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *t_elt = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-l") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *length = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-post") == 0) {
      *post = 1;
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *method = 2;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *method = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}



/**
 *
 * \brief  Main
 *
 */

int main
(
 int   argc,
 char *argv[]
 )
{
  /*
   *  Set default values
   */
  PDM_g_num_t          nx         = 10;
  PDM_g_num_t          ny         = 10;
  PDM_g_num_t          nz         = 10;
  int                  n_dom_i    = 1;
  int                  n_dom_j    = 1;
  int                  n_dom_k    = 1;
  int                  periodic_i = 0;
  int                  periodic_j = 0;
  int                  periodic_k = 0;
  double               length     = 1.;
  int                  n_part     = 1;
  int                  post       = 0;
  PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_TRIA3;
  // 2 -> tria
  // 3 -> quad
  // 5 -> tetra
  // 6 -> pyramid
  // 7 -> prism
  // 8 -> hexa

#ifdef PDM_HAVE_PARMETIS
  PDM_part_split_t method  = PDM_PART_SPLIT_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_part_split_t method  = PDM_PART_SPLIT_PTSCOTCH;
#endif
#endif

  /*
   *  Read args
   */

  _read_args(argc,
             argv,
             &nx,
             &ny,
             &nz,
             &n_dom_i,
             &n_dom_j,
             &n_dom_k,
             &periodic_i,
             &periodic_j,
             &periodic_k,
             (int *) &t_elt,
             &length,
             &n_part,
             &post,
             (int *) &method);

  /*
   *  Init
   */
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);



  /*
   *  Initialize structs
   */
  int n_interface =
  n_dom_j*n_dom_k*(n_dom_i - 1 + periodic_i) +
  n_dom_k*n_dom_i*(n_dom_j - 1 + periodic_j) +
  n_dom_i*n_dom_j*(n_dom_k - 1 + periodic_k);

  int n_domain = n_dom_i * n_dom_j * n_dom_k;

  PDM_dcube_nodal_t **dcube = (PDM_dcube_nodal_t **) malloc(sizeof(PDM_dcube_nodal_t *) * n_domain);
  PDM_dmesh_nodal_t **dmn   = (PDM_dmesh_nodal_t **) malloc(sizeof(PDM_dmesh_nodal_t *) * n_domain);

  const int order = 1;

  int i_domain = 0;

  for (int k = 0; k < n_dom_k; k++) {
    for (int j = 0; j < n_dom_j; j++) {
      for (int i = 0; i < n_dom_i; i++) {
        dcube[i_domain] = PDM_dcube_nodal_gen_create(comm,
                                                     nx,
                                                     ny,
                                                     nz,
                                                     length,
                                                     length*i,
                                                     length*j,
                                                     length*k,
                                                     t_elt,
                                                     order,
                                                     PDM_OWNERSHIP_KEEP);
        PDM_dcube_nodal_gen_build (dcube[i_domain]);
        dmn[i_domain] = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube[i_domain]);
        PDM_dmesh_nodal_generate_distribution(dmn[i_domain]);

        i_domain++;
      }
    }
  }

  PDM_domain_interface_t *dom_intrf = PDM_domain_interface_create(n_interface,
                                                                  n_domain,
                                                                  PDM_DOMAIN_INTERFACE_MULT_NO,
                                                                  PDM_OWNERSHIP_KEEP,
                                                                  comm);

  /*
   *  Interfaces
   */
  int          *interface_dn  = (int          *) malloc(sizeof(int          ) * n_interface);
  PDM_g_num_t **interface_ids = (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_interface);
  int         **interface_dom = (int         **) malloc(sizeof(int         *) * n_interface);

  int i_interface = 0;

  /* i-direction */
  PDM_g_num_t *distrib_i = PDM_compute_uniform_entity_distribution(comm,
                                                                   ny * nz);

  for (int i = 0; i < n_dom_i - 1 + periodic_i; i++) {
    for (int k = 0; k < n_dom_k; k++) {
      for (int j = 0; j < n_dom_j; j++) {

        int i_domain1 = i + n_dom_i*(j + n_dom_j*k);// +1?
        int i_domain2 = (i+1)%n_dom_i + n_dom_i*(j + n_dom_j*k);// +1?

        interface_dn[i_interface] = (int) (distrib_i[i_rank+1] - distrib_i[i_rank]);

        interface_dom[i_interface] = (int *) malloc(sizeof(int) * 2);
        interface_dom[i_interface][0] = i_domain1;
        interface_dom[i_interface][1] = i_domain2;

        interface_ids[i_interface] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 2 * interface_dn[i_interface]);

        for (int idx = 0; idx < interface_dn[i_interface]; idx++) {
          PDM_g_num_t g = distrib_i[i_rank] + idx;

          PDM_g_num_t jj = g % ny;
          PDM_g_num_t kk = g / ny;

          interface_ids[i_interface][2*idx  ] = 1 + nx + nx*(jj + ny*kk);
          interface_ids[i_interface][2*idx+1] = 1 +      nx*(jj + ny*kk);
        }

        i_interface++;
      }
    }
  }
  free(distrib_i);



  /* j-direction */
  PDM_g_num_t *distrib_j = PDM_compute_uniform_entity_distribution(comm,
                                                                   ny * nz);

  for (int j = 0; j < n_dom_j - 1 + periodic_j; j++) {
    for (int i = 0; i < n_dom_i; i++) {
      for (int k = 0; k < n_dom_k; k++) {

        int i_domain1 = i + n_dom_i*(j + n_dom_j*k);// +1?
        int i_domain2 = i + n_dom_i*((j+1)%n_dom_j + n_dom_j*k);// +1?

        interface_dn[i_interface] = (int) (distrib_j[i_rank+1] - distrib_j[i_rank]);

        interface_dom[i_interface] = (int *) malloc(sizeof(int) * 2);
        interface_dom[i_interface][0] = i_domain1;
        interface_dom[i_interface][1] = i_domain2;

        interface_ids[i_interface] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 2 * interface_dn[i_interface]);

        for (int idx = 0; idx < interface_dn[i_interface]; idx++) {
          PDM_g_num_t g = distrib_j[i_rank] + idx;

          PDM_g_num_t kk = g % nz;
          PDM_g_num_t ii = g / nz;

          interface_ids[i_interface][2*idx  ] = 1 + ii + nx*(ny + ny*kk);
          interface_ids[i_interface][2*idx+1] = 1 + ii + nx*(     ny*kk);
        }

        i_interface++;
      }
    }
  }
  free(distrib_j);



  /* k-direction */
  PDM_g_num_t *distrib_k = PDM_compute_uniform_entity_distribution(comm,
                                                                   ny * nz);

  for (int k = 0; k < n_dom_k - 1 + periodic_k; k++) {
    for (int j = 0; j < n_dom_j; j++) {
      for (int i = 0; i < n_dom_i; i++) {

        int i_domain1 = i + n_dom_i*(j + n_dom_j*k);// +1?
        int i_domain2 = i + n_dom_i*(j + n_dom_j*(k+1)%n_dom_k);// +1?

        interface_dn[i_interface] = (int) (distrib_k[i_rank+1] - distrib_k[i_rank]);

        interface_dom[i_interface] = (int *) malloc(sizeof(int) * 2);
        interface_dom[i_interface][0] = i_domain1;
        interface_dom[i_interface][1] = i_domain2;

        interface_ids[i_interface] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 2 * interface_dn[i_interface]);

        for (int idx = 0; idx < interface_dn[i_interface]; idx++) {
          PDM_g_num_t g = distrib_k[i_rank] + idx;

          PDM_g_num_t ii = g % nx;
          PDM_g_num_t jj = g / nx;

          interface_ids[i_interface][2*idx  ] = 1 + ii + nx*(jj + ny*nz);
          interface_ids[i_interface][2*idx+1] = 1 + ii + nx*(jj        );
        }

        i_interface++;
      }
    }
  }
  free(distrib_k);

  // printf("i_interface = %d / %d\n", i_interface, n_interface);



  PDM_domain_interface_set (dom_intrf,
                            PDM_BOUND_TYPE_VTX,
                            interface_dn,
                            interface_ids,
                            interface_dom);


  /*
   *  Free memory
   */
  for (int i = 0; i < n_domain; i++) {
    PDM_dcube_nodal_gen_free(dcube[i]);
    // PDM_dcube_nodal_gen_free(dmn[i]);
  }

  PDM_domain_interface_free(dom_intrf);

  for (int i = 0; i < n_interface; i++) {
    free(interface_ids[i]);
    free(interface_dom[i]);
  }
  free(interface_dn);
  free(interface_ids);
  free(interface_dom);
  free(dcube);
  free(dmn);

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;
}
