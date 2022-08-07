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
#include "pdm_config.h"
#include "pdm_mpi.h"
#include "pdm_multipart.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_distrib.h"
#include "pdm_array.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_extension.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_vtk.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_logging.h"
#include "pdm_order.h"
#include "pdm_binary_search.h"

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

static void
_compute_face_vtx
(
 int   n_face,
 int  *face_edge_idx,
 int  *face_edge,
 int  *edge_vtx,
 int **face_vtx
 )
{
  int dbg = 0;

  *face_vtx = malloc (sizeof(int) * face_edge_idx[n_face]);

  for(int i_face = 0; i_face <  face_edge_idx[n_face]; ++i_face) {
    (*face_vtx)[i_face] = 10000;
  }

  int n_edge = 0;
  for (int i = 0; i < face_edge_idx[n_face]; i++) {
    n_edge = PDM_MAX(n_edge, PDM_ABS(face_edge[i]));
  }


  int *edge_tag = PDM_array_zeros_int(n_edge);

  for (int iface = 0; iface < n_face; iface++) {
    int *_face_vtx  = *face_vtx  + face_edge_idx[iface];
    int *_face_edge =  face_edge + face_edge_idx[iface];

    if (dbg) {
      log_trace("\nFace %d\n", iface);
      for (int idx_edge = face_edge_idx[iface]; idx_edge < face_edge_idx[iface+1]; idx_edge++) {
        int iedge = PDM_ABS(face_edge[idx_edge]) - 1;
        log_trace("  edge %d: %d %d\n",
                  face_edge[idx_edge],
                  edge_vtx[2*iedge], edge_vtx[2*iedge+1]);
      }
    }

    int _n_edge = face_edge_idx[iface+1] - face_edge_idx[iface];
    // first edge
    int iedge = PDM_ABS(_face_edge[0]) - 1;
    edge_tag[iedge] = 1;
    _face_vtx[0] = PDM_ABS(edge_vtx[2*iedge  ]);
    _face_vtx[1] = PDM_ABS(edge_vtx[2*iedge+1]);

    for (int i = 2; i < _n_edge; i++) {

      for (int j = 1; j < _n_edge; j++) {
        iedge = PDM_ABS(_face_edge[j]) - 1;

        if (edge_tag[iedge]) {
          continue;
        }

        if (edge_vtx[2*iedge] == _face_vtx[i-1]) {
          _face_vtx[i] = PDM_ABS(edge_vtx[2*iedge+1]);
          edge_tag[iedge] = 1;
          break;
        }
        else if (edge_vtx[2*iedge+1] == _face_vtx[i-1]) {
          _face_vtx[i] = PDM_ABS(edge_vtx[2*iedge]);
          edge_tag[iedge] = 1;
          break;
        }
      }
    }

    if (dbg) {
      log_trace("  face_vtx = ");
      for (int ivtx = 0; ivtx < face_edge_idx[iface+1] - face_edge_idx[iface]; ivtx++) {
        log_trace("%d ", _face_vtx[ivtx]);
      }
      log_trace("\n");
    }

    // reset tags
    for (int i = 0; i < _n_edge; i++) {
      iedge = PDM_ABS(_face_edge[i]) - 1;
      edge_tag[iedge] = 0;
    }
  }
  free(edge_tag);
}



/**
 *
 * \brief  Main
 *
 */
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -pi -t 8
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -t 8 -ni 2 -pi 1
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -t 8 -pi
// mpirun -np 1 ./test/pdm_t_join_domains -nx 3 -ny 3 -nz 2 -t 8 -pi -pj
// Repro LS89 : mpirun -np 2 ./test/pdm_t_join_domains -nx 5 -ny 4 -nz 2 -pt-scotch -pj
int main
(
 int   argc,
 char *argv[]
)
{
  // PDM_MPI_Init(&argc, &argv);

  // int size          = 4;
  // int array[2*4] = {4, 0, 2, 1, 2, 0, 1, 1};
  // int* lorder = malloc(size * sizeof(int));
  // PDM_order_lnum_s(array, 2, lorder, size);

  // PDM_order_array(size, 2 * sizeof(int), lorder, array);

  // PDM_log_trace_array_int(lorder,   size, "lorder : ");
  // PDM_log_trace_array_int(array, 2 * size, "Sorted array : ");


  // // Search
  // int elt_to_find1[] = {1, 1};
  // int elt_to_find2[] = {2, 0};
  // int elt_to_find3[] = {2, 1};
  // int elt_to_find4[] = {4, 0};
  // int elt_to_find5[] = {4, 1};
  // int pos1 = PDM_order_binary_search_long(elt_to_find1, array, 2, size);
  // int pos2 = PDM_order_binary_search_long(elt_to_find2, array, 2, size);
  // int pos3 = PDM_order_binary_search_long(elt_to_find3, array, 2, size);
  // int pos4 = PDM_order_binary_search_long(elt_to_find4, array, 2, size);
  // int pos5 = PDM_order_binary_search_long(elt_to_find5, array, 2, size);

  // printf("pos1 = %i \n", pos1);
  // printf("pos2 = %i \n", pos2);
  // printf("pos3 = %i \n", pos3);
  // printf("pos4 = %i \n", pos4);
  // printf("pos5 = %i \n", pos5);

  // free(lorder);

  // PDM_MPI_Finalize();
  // return 0;
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
  // PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_TRIA3;
  PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_QUAD4;
  // PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_HEXA8;
  // PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_TETRA4;
  // PDM_Mesh_nodal_elt_t t_elt      = PDM_MESH_NODAL_PRISM6;
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

  // const int renum_properties_cell[6] = {1024, 0, 1, 64, 3, 1};
  // const int renum_properties_cell[6] = {12, 0, 1, 64, 3, 1};
  // const int renum_properties_cell[6] = {256, 0, 1, 64, 1, 1};
  // const int renum_properties_cell[6] = {16, 0, 1, 64, 1, 1};
  // PDM_multipart_set_reordering_options(mpart_id, -1, "PDM_PART_RENUM_CELL_HPC",
  //                                           (void * ) renum_properties_cell,
  //                                                    "PDM_PART_RENUM_FACE_NONE");
  PDM_multipart_set_reordering_options(mpart_id,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");
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
  // int          *dn_edge       = malloc( n_domain * sizeof(int          ));
  // int         **dface_vtx_idx = malloc( n_domain * sizeof(int         *));
  // PDM_g_num_t **dface_vtx     = malloc( n_domain * sizeof(PDM_g_num_t *));
  // PDM_g_num_t **dedge_vtx     = malloc( n_domain * sizeof(PDM_g_num_t *));
  // int         **dedge_vtx_idx = malloc( n_domain * sizeof(int *        ));
 
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

  //   dedge_vtx    [i] = NULL;
  //   dedge_vtx_idx[i] = NULL;
  //   dn_edge[i] = PDM_dmesh_connectivity_get(dm[i],
  //                                           PDM_CONNECTIVITY_TYPE_EDGE_VTX,
  //                                           &dedge_vtx[i],
  //                                           &dedge_vtx_idx[i],
  //                                           PDM_OWNERSHIP_KEEP);
  //   assert(dedge_vtx_idx[i] != NULL);
  //   assert(dedge_vtx[i] != NULL);

  //   // PDM_log_trace_array_int(dedge_vtx_idx[i], dn_edge[i], "dedge_vtx_idx[i] :: ");
  //   // PDM_log_trace_array_long(dedge_vtx[i], dedge_vtx_idx[i][dn_edge[i]], "dedge_vtx[i] :: ");

  // }

  /*
   * Transform interface by vtx by interface by faces
   */
  // PDM_domain_interface_translate_vtx2face(dom_intrf,
  //                                         dn_vtx,
  //                                         dn_face,
  //                                         dface_vtx_idx,
  //                                         dface_vtx);

  // PDM_domain_interface_translate_vtx2edge(dom_intrf,
  //                                         dn_vtx,
  //                                         dn_edge,
  //                                         dedge_vtx_idx,
  //                                         dedge_vtx);

  /*
   *  Prepare pointer by domain and by part
   */
  int           *pn_n_part      = (int           *) malloc( n_domain * sizeof(int          *));
  int          **pn_face        = (int          **) malloc( n_domain * sizeof(int          *));
  int          **pn_edge        = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t ***pface_ln_to_gn = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  PDM_g_num_t ***pedge_ln_to_gn = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  int          **pn_vtx         = (int          **) malloc( n_domain * sizeof(int          *));
  PDM_g_num_t ***pvtx_ln_to_gn  = (PDM_g_num_t ***) malloc( n_domain * sizeof(PDM_g_num_t **));
  int         ***pface_vtx      = (int         ***) malloc( n_domain * sizeof(int         **));
  int         ***pedge_vtx_idx  = (int         ***) malloc( n_domain * sizeof(int         **));
  int         ***pedge_vtx      = (int         ***) malloc( n_domain * sizeof(int         **));
  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    pn_n_part    [i_dom] = n_part;
    pn_face       [i_dom] = (int          *) malloc( n_part * sizeof(int          ));
    pn_edge       [i_dom] = (int          *) malloc( n_part * sizeof(int          ));
    pface_ln_to_gn[i_dom] = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    pedge_ln_to_gn[i_dom] = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    pn_vtx        [i_dom] = (int          *) malloc( n_part * sizeof(int          ));
    pvtx_ln_to_gn [i_dom] = (PDM_g_num_t **) malloc( n_part * sizeof(PDM_g_num_t *));
    pface_vtx     [i_dom] = (int         **) malloc( n_part * sizeof(int         *));
    pedge_vtx_idx [i_dom] = (int         **) malloc( n_part * sizeof(int         *));
    pedge_vtx     [i_dom] = (int         **) malloc( n_part * sizeof(int         *));

    for (int i_part = 0; i_part < pn_n_part[i_dom]; i_part++) {
      pn_face[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart_id,
                                                               i_dom,
                                                               i_part,
                                                               PDM_MESH_ENTITY_FACE,
                                                               &pface_ln_to_gn[i_dom][i_part],
                                                               PDM_OWNERSHIP_KEEP);
      pn_edge[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart_id,
                                                               i_dom,
                                                               i_part,
                                                               PDM_MESH_ENTITY_EDGE,
                                                               &pedge_ln_to_gn[i_dom][i_part],
                                                               PDM_OWNERSHIP_KEEP);
      pn_vtx[i_dom][i_part] = PDM_multipart_part_ln_to_gn_get(mpart_id,
                                                              i_dom,
                                                              i_part,
                                                              PDM_MESH_ENTITY_VERTEX,
                                                              &pvtx_ln_to_gn[i_dom][i_part],
                                                              PDM_OWNERSHIP_KEEP);

      PDM_multipart_part_connectivity_get(mpart_id,
                                          i_dom,
                                          i_part,
                                          PDM_CONNECTIVITY_TYPE_EDGE_VTX,
                                          &pedge_vtx     [i_dom][i_part],
                                          &pedge_vtx_idx [i_dom][i_part],
                                          PDM_OWNERSHIP_KEEP);
      assert(pedge_vtx_idx [i_dom][i_part] == NULL);
      pedge_vtx_idx [i_dom][i_part] = malloc((pn_edge[i_dom][i_part]+1) * sizeof(int));

      for(int i = 0; i < pn_edge[i_dom][i_part]+1; ++i) {
        pedge_vtx_idx [i_dom][i_part][i] = 2*i;
      }

    }
  }

  /*
   * Deduction en partition du graphe entre domaine
   */
  PDM_part_domain_interface_t* pdi = PDM_domain_interface_to_part_domain_interface(dom_intrf,
                                                                                   pn_n_part,
                                                                                   pn_face,
                                                                                   pn_edge,
                                                                                   pn_vtx,
                                                                                   pface_ln_to_gn,
                                                                                   pedge_ln_to_gn,
                                                                                   pvtx_ln_to_gn);

  PDM_part_domain_interface_translate(pdi,
                                      PDM_BOUND_TYPE_VTX,
                                      PDM_BOUND_TYPE_EDGE,
                                      pn_n_part,
                                      pn_vtx,
                                      pn_edge,
                                      pedge_ln_to_gn,
                                      pedge_vtx_idx,
                                      pedge_vtx);



  for (int i = 0; i < n_domain; i++) {
    PDM_dcube_nodal_gen_free(dcube[i]);
    // PDM_dcube_nodal_gen_free(dmn[i]);
  }
  PDM_multipart_free(mpart_id);

  if(pdi != NULL) {
    PDM_part_domain_interface_free(pdi);
  }
  PDM_UNUSED(pdi);

  // free(dm);
  // free(dn_vtx);
  // free(dn_face);
  // free(dn_edge);
  // free(dface_vtx_idx);
  // free(dface_vtx);
  // free(dedge_vtx_idx);
  // free(dedge_vtx);

  // PDM_dmesh_nodal_to_dmesh_free(dmn_to_dm);
  PDM_domain_interface_free(dom_intrf);
  free(n_part_by_domain);
  free(dcube);
  // free(dmn);

  for (int i_dom = 0; i_dom < n_domain; i_dom++) {
    for (int i_part = 0; i_part < pn_n_part[i_dom]; i_part++) {
      free(pedge_vtx_idx [i_dom][i_part]);
    }
    free(pn_face       [i_dom]);
    free(pn_edge       [i_dom]);
    free(pface_ln_to_gn[i_dom]);
    free(pedge_ln_to_gn[i_dom]);
    free(pn_vtx        [i_dom]);
    free(pvtx_ln_to_gn [i_dom]);
    free(pface_vtx     [i_dom]);
    free(pedge_vtx_idx [i_dom]);
    free(pedge_vtx     [i_dom]);
  }

  free(pn_n_part     );
  free(pn_face       );
  free(pn_edge       );
  free(pface_ln_to_gn);
  free(pedge_ln_to_gn);
  free(pn_vtx        );
  free(pvtx_ln_to_gn );
  free(pface_vtx     );
  free(pedge_vtx_idx );
  free(pedge_vtx     );

  if (i_rank == 0) {
    printf("-- End\n");
    fflush(stdout);
  }

  PDM_MPI_Finalize();

  return 0;
}
