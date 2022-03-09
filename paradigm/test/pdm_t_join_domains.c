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
#include "pdm_multipart.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_mesh_nodal.h"

#include "pdm_dcube_nodal_gen.h"
#include "pdm_domain_interface.h"
#include "pdm_dmesh_nodal_to_dmesh.h"

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
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -pi -t 8
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
  PDM_split_dual_t method  = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t method  = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
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
  int dim = PDM_Mesh_nodal_elt_dim_get(t_elt);

  if (dim == 2) {
    n_dom_k    = 1;
    nz         = 1;
    periodic_k = 0;
  }

  // int n_interface =
  // n_dom_j*n_dom_k*(n_dom_i - 1 + periodic_i) +
  // n_dom_k*n_dom_i*(n_dom_j - 1 + periodic_j) +
  // n_dom_i*n_dom_j*(n_dom_k - 1 + periodic_k);

  int n_domain = n_dom_i * n_dom_j * n_dom_k;

  // PDM_dcube_nodal_t **dcube = (PDM_dcube_nodal_t **) malloc(sizeof(PDM_dcube_nodal_t *) * n_domain);
  PDM_dmesh_nodal_t **dmn   = (PDM_dmesh_nodal_t **) malloc(sizeof(PDM_dmesh_nodal_t *) * n_domain);
  PDM_dcube_nodal_t **dcube = NULL;
  PDM_domain_interface_t *dom_intrf = NULL;

  const int order = 1;

  PDM_dcube_nodal_cart_topo(comm,
                            n_dom_i,
                            n_dom_j,
                            n_dom_k,
                            periodic_i,
                            periodic_j,
                            periodic_k,
                            nx,
                            ny,
                            nz,
                            length,
                            0.,
                            0.,
                            0.,
                            t_elt,
                            order,
                            &dcube,
                            &dom_intrf,
                            PDM_OWNERSHIP_KEEP);
  // /* Check */
  // if (n_rank == 1) {

  //   double vect[4][3] = {
  //     {0.,             0.,             0.},
  //     {length*n_dom_i, 0.,             0.},
  //     {0,              length*n_dom_j, 0.},
  //     {0.,             0.,             length*n_dom_k}
  //   };

  //   double **all_coord = (double **) malloc(sizeof(double *) * n_domain);
  //   for (int i = 0; i < n_domain; i++) {
  //     all_coord[i] = PDM_DMesh_nodal_vtx_get(dmn[i]);
  //   }

  //   for (i_interface = 0; i_interface < n_interface; i_interface++) {

  //     int i_domain1 = interface_dom[i_interface][0];// - 1?
  //     int i_domain2 = interface_dom[i_interface][1];// - 1?

  //     for (int i = 0; i < interface_dn[i_interface]; i++) {

  //       PDM_g_num_t ivtx_1 = interface_ids[i_interface][2*i]   - 1;
  //       PDM_g_num_t ivtx_2 = interface_ids[i_interface][2*i+1] - 1;

  //       // if(ivtx_1 < 0 || ivtx_1 >= PDM_DMesh_nodal_n_vtx_get(dmn[i_domain1])) {
  //       //   printf("error : interface %d, domain 1:%d, pair %d, vtx "PDM_FMT_G_NUM"/%d\n",
  //       //          i_interface, i_domain1, i, ivtx_1+1,
  //       //          PDM_DMesh_nodal_n_vtx_get(dmn[i_domain1]));
  //       // }
  //       // if(ivtx_2 < 0 || ivtx_2 >= PDM_DMesh_nodal_n_vtx_get(dmn[i_domain2])) {
  //       //   printf("error : interface %d, domain 2:%d, pair %d, vtx "PDM_FMT_G_NUM"/%d\n",
  //       //          i_interface, i_domain2, i, ivtx_2+1,
  //       //          PDM_DMesh_nodal_n_vtx_get(dmn[i_domain2]));
  //       // }

  //       double *vtx_1 = all_coord[i_domain1] + 3*ivtx_1;
  //       double *vtx_2 = all_coord[i_domain2] + 3*ivtx_2;

  //       double dist2 = 0;
  //       for (int j = 0; j < 3; j++) {
  //         double delta = vtx_1[j] - (vtx_2[j] + vect[i_period[i_interface]][j]);
  //         dist2 += delta * delta;
  //       }

  //       // assert(dist2 < 1e-9);
  //       if (dist2 > 1e-9) {
  //         printf("!! interface %d, domains %d %d, vtx ("PDM_FMT_G_NUM": %f %f %f) ("PDM_FMT_G_NUM": %f %f %f), dist = %f\n",
  //                i_interface, i_domain1, i_domain2,
  //                ivtx_1+1,
  //                vtx_1[0], vtx_1[1], vtx_1[2],
  //                ivtx_2+1,
  //                vtx_2[0], vtx_2[1], vtx_2[2],
  //                sqrt(dist2));
  //       }
  //     }
  //   }
  //   free (all_coord);
  // }

  /*
   * Partitionnement
   */
  int* n_part_by_domain = (int *) malloc(n_domain * sizeof(int));
  for(int i = 0; i < n_domain; ++i) {
    n_part_by_domain[i] = n_part;
  }
  PDM_multipart_t *mpart_id = PDM_multipart_create(n_domain,
                                                   n_part_by_domain,
                                                   PDM_FALSE,
                                                   method,
                                                   PDM_PART_SIZE_HOMOGENEOUS,
                                                   NULL,
                                                   comm,
                                                   PDM_OWNERSHIP_KEEP);

  for (int i = 0; i < n_domain; i++) {
    dmn[i] = PDM_dcube_nodal_gen_dmesh_nodal_get(dcube[i]);
    PDM_dmesh_nodal_generate_distribution(dmn[i]);
    PDM_multipart_register_dmesh_nodal(mpart_id, i, dmn[i]);
  }

  PDM_multipart_run_ppart(mpart_id);



  /*
   * Create dmesh_nodal_to_dmesh to setup face and edge
   */
  // PDM_dmesh_nodal_to_dmesh_t* dmn_to_dm = PDM_dmesh_nodal_to_dmesh_create(n_domain, comm, PDM_OWNERSHIP_KEEP);
  // for (int i = 0; i < n_domain; i++) {
  //   PDM_dmesh_nodal_to_dmesh_add_dmesh_nodal(dmn_to_dm, i, dmn[i]);
  // }

  // PDM_dmesh_nodal_to_dmesh_compute(dmn_to_dm,
  //                                  PDM_DMESH_NODAL_TO_DMESH_TRANSFORM_TO_FACE,
  //                                  PDM_DMESH_NODAL_TO_DMESH_TRANSLATE_GROUP_TO_FACE);

  // PDM_dmesh_t **dm = malloc(n_domain * sizeof(PDM_dmesh_t *));

  // int          *dn_vtx        = malloc( n_domain * sizeof(int          ));
  // int          *dn_face       = malloc( n_domain * sizeof(int          ));
  // int         **dface_vtx_idx = malloc( n_domain * sizeof(int         *));
  // PDM_g_num_t **dface_vtx     = malloc( n_domain * sizeof(PDM_g_num_t *));

  // for (int i = 0; i < n_domain; i++) {
  //   PDM_dmesh_nodal_to_dmesh_get_dmesh(dmn_to_dm, 0, &dm[i]);

  //   int _dn_cell, _dn_face, _dn_edge, _dn_vtx, _n_bnd, _n_join;
  //   PDM_dmesh_dims_get(dm[i], &_dn_cell, &_dn_face, &_dn_edge, &_dn_vtx, &_n_bnd, &_n_join);

  //   dn_vtx [i] = _dn_vtx;
  //   dn_face[i] = PDM_dmesh_connectivity_get(dm[i],
  //                                           PDM_CONNECTIVITY_TYPE_FACE_VTX,
  //                                           &dface_vtx[i],
  //                                           &dface_vtx_idx[i],
  //                                           PDM_OWNERSHIP_KEEP);
  //   assert(dface_vtx_idx[i] != NULL);
  //   assert(dface_vtx[i] != NULL);
  // }

  // /*
  //  * Transform interface by vtx by interface by faces
  //  */
  // PDM_domain_interface_translate_vtx2face(dom_intrf,
  //                                         dn_vtx,
  //                                         dn_face,
  //                                         dface_vtx_idx,
  //                                         dface_vtx);


  /*
   *  Prepare pointer by domain and by part
   */
  int           *pn_n_part     = (int           *) malloc( n_domain * sizeof(int          *));
  int          **pn_vtx        = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t ***pvtx_ln_to_gn = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    pn_n_part    [i_dom] = n_part;
    pn_vtx       [i_dom] = (int          *) malloc( n_part * sizeof(int          ));
    pvtx_ln_to_gn[i_dom] = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    for(int i_part = 0; i_part < n_part; ++i_part) {
      pn_vtx[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart_id,
                                                                 i_dom,
                                                                 i_part,
                                                                 PDM_MESH_ENTITY_VERTEX,
                                                                 &pvtx_ln_to_gn[i_dom][i_part],
                                                                 PDM_OWNERSHIP_KEEP);
    }
  }


  PDM_part_domain_interface_t* pdi = PDM_domain_interface_to_part_domain_interface(dom_intrf,
                                                                                   pn_n_part,
                                                                                   NULL,
                                                                                   NULL,
                                                                                   pn_vtx,
                                                                                   NULL,
                                                                                   NULL,
                                                                                   pvtx_ln_to_gn);



  /*
   *  Pour le debug : extration des faces avec le extract part + vtk
   */

  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    free(pn_vtx       [i_dom]);
    free(pvtx_ln_to_gn[i_dom]);
  }
  free(pn_vtx       );
  free(pvtx_ln_to_gn);
  free(pn_n_part);

  /*
   *  Free memory
   */
  // free(i_period);

  for (int i = 0; i < n_domain; i++) {
    PDM_dcube_nodal_gen_free(dcube[i]);
    // PDM_dcube_nodal_gen_free(dmn[i]);
  }
  PDM_multipart_free(mpart_id);

  PDM_part_domain_interface_free(pdi);
  PDM_UNUSED(pdi);

  // free(dm);
  // free(dn_vtx);
  // free(dn_face);
  // free(dface_vtx_idx);
  // free(dface_vtx);

  // PDM_dmesh_nodal_to_dmesh_free(dmn_to_dm);
  PDM_domain_interface_free(dom_intrf);
  free(n_part_by_domain);

  // for (int i = 0; i < n_interface; i++) {
  //   free(interface_ids[i]);
  //   free(interface_dom[i]);
  // }
  // free(interface_dn);
  // free(interface_ids);
  // free(interface_dom);
  free(dcube);
  free(dmn);

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;
}
